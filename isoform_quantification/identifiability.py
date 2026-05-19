"""
Identifiability analysis for miniQuant multi-platform isoform quantification.

For each gene with T isoforms, computes two sets of identifiability metrics:

1. Restricted Jacobian k-value (sigma_max / sigma_min_pos):
     LR:    J_L = A^(L) B
     SR:    J_S = A^(S) J_phi B
     Hybrid: J_H = vstack([J_L, J_S])

2. Restricted Fisher information (multinomial model):
     F_L = J_L^T W_L J_L,   W_L = N_L * diag(1 / p^(L)),  p^(L) = A^(L) theta
     F_S = J_S^T W_S J_S,   W_S = N_S * diag(1 / p^(S)),  p^(S) = A^(S) phi(theta)
     F_H = F_L + F_S
   Metrics: condition number (lambda_max / lambda_min_pos) and lambda_min.

where B is the T x (T-1) tangent basis for {v : 1^T v = 0},
phi(theta) = D theta / (1^T D theta),  D = diag(isoform SR effective lengths),
and theta is taken from EM quantification results (Isoform_abundance.out).

If no EM results are found, only the original k-values (SVD of A) are output.

Outputs: identifiability.tsv  (one row per gene)
"""

import math
import os
import csv
import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed


# ---------------------------------------------------------------------------
# Step 1: Tangent space basis
# ---------------------------------------------------------------------------

def compute_tangent_basis(T):
    """
    Return a T x (T-1) matrix B whose columns form an orthonormal basis
    for the tangent space {v in R^T : 1^T v = 0}.
    Uses the Helmert contrast matrix (closed-form, unique):
        b_k = 1/sqrt(k(k+1)) * (1,...,1, -k, 0,...,0),  k = 1,...,T-1
    The choice of basis does not affect any computed metric (eigenvalues,
    condition number, SE/CI) because all quantities are invariant under
    orthogonal rotation within the tangent space.
    """
    B = np.zeros((T, T - 1))
    for k in range(1, T):
        B[:k, k-1] = 1.0 / np.sqrt(k * (k + 1))
        B[k,  k-1] = -k  / np.sqrt(k * (k + 1))
    return B  # shape (T, T-1), orthonormal columns


# ---------------------------------------------------------------------------
# Step 2: SR isoform effective lengths
# ---------------------------------------------------------------------------

def compute_isoform_eff_len_sr(isoform_lengths, sr_read_len):
    """
    Compute SR isoform effective length for each isoform t:
        l_tilde_t = max(l_t - r + 1, 1)
    where l_t is the transcript length and r is the SR read length.
    """
    return np.maximum(isoform_lengths - sr_read_len + 1, 1.0)


# ---------------------------------------------------------------------------
# Step 4: Jacobian of phi
# ---------------------------------------------------------------------------

def compute_Jphi(theta_hat, l_tilde):
    """
    T x T Jacobian of phi(theta) = D theta / (1^T D theta),  D = diag(l_tilde).

    J_phi_{ij} = D_i/s * delta_{ij} - (D_i theta_i)(D_j) / s^2
    where s = 1^T D theta_hat = sum_t l_tilde_t * theta_hat_t.
    """
    s = float(l_tilde @ theta_hat)
    if s == 0.0:
        s = 1.0
    Jphi = np.diag(l_tilde / s) - np.outer(l_tilde * theta_hat, l_tilde) / (s * s)
    return Jphi


# ---------------------------------------------------------------------------
# Step 5: Identifiability metrics from restricted Jacobian / Fisher
# ---------------------------------------------------------------------------

def _metrics_from_Jtan(J_tan):
    """
    Compute k-value = sigma_max / sigma_min_positive and the two singular values
    from the restricted Jacobian J_tan.
    Returns (k_value, sigma_max, sigma_min_pos).
    """
    sv = np.linalg.svd(J_tan, compute_uv=False)
    rank = np.linalg.matrix_rank(J_tan)
    if rank == 0 or sv[0] == 0.0:
        return math.nan, math.nan, math.nan
    sigma_max = float(sv[0])
    sigma_min_pos = float(sv[rank - 1])
    k_value = sigma_max / sigma_min_pos if sigma_min_pos > 0 else math.inf
    return k_value, sigma_max, sigma_min_pos


_FISHER_EIG_TOL = 1e-5   # eigenvalue threshold for positive-definiteness

def _metrics_from_Fisher(F_tan):
    """
    Compute Fisher condition number and identifiability flag from the restricted
    Fisher matrix F_tan (symmetric, (T-1) x (T-1)).

    Returns (fisher_lambda_cond, lambda_min, fisher_identifiable, fisher_det):
      fisher_lambda_cond : log10(lambda_max / lambda_min)  (nan if lambda_min == 0)
      lambda_min          : smallest eigenvalue thresholded at _FISHER_EIG_TOL;
                            eigenvalues < _FISHER_EIG_TOL are set to 0.
                            0 for non-identifiable matrices, actual value otherwise.
      fisher_identifiable : True if F_tan is positive definite (all eigenvalues
                            >= _FISHER_EIG_TOL), else False.
      fisher_det          : log10(|det(F_tan)| + 1) — D-optimality metric in log10 scale;
                            0 when singular (det=0); nan on numerical sign error.
    """
    eigvals = np.linalg.eigvalsh(F_tan)          # ascending order
    eigvals_clipped = np.maximum(eigvals, 0.0)    # PSD guarantee; clip fp negatives
    lambda_max = float(eigvals_clipped[-1])
    sign, logabsdet = np.linalg.slogdet(F_tan)
    if sign < 0:
        fisher_det = math.nan   # sign < 0: fp overflow artifact for PSD matrix
    else:
        # log10(|det| + 1); logaddexp(logabsdet, 0) = ln(exp(logabsdet)+1) stably
        # sign=0 → logabsdet=-inf → logaddexp(-inf,0)=0 → fisher_det=0
        fisher_det = float(np.logaddexp(logabsdet, 0.0) / math.log(10))
    # eigenvalues < _FISHER_EIG_TOL are treated as 0; lambda_min = 0 for non-identifiable
    thresholded = np.where(eigvals_clipped < _FISHER_EIG_TOL, 0.0, eigvals_clipped)
    lambda_min = float(thresholded[0])
    fisher_lambda_cond = math.log10(lambda_max / lambda_min) if lambda_min > 0 else math.nan
    is_identifiable = bool(lambda_min > 0)   # positive definite: lambda_min > 0 iff all eigenvalues >= _FISHER_EIG_TOL
    return fisher_lambda_cond, lambda_min, is_identifiable, fisher_det


# ---------------------------------------------------------------------------
# SE and CI from restricted Fisher
# ---------------------------------------------------------------------------

_Z_95 = 1.959964   # z_{0.975} for 95% CI

def _compute_se_ci(F_tan, B, theta_hat, gamma=0.05):
    """
    Compute per-isoform SE and truncated Wald CI via the Moore-Penrose pseudoinverse.
    When F_tan is positive definite (invertible), the pseudoinverse equals the regular
    inverse and all isoforms receive finite SE; the pseudoinverse formulation is used
    uniformly for numerical stability and to handle the singular case without branching.

    Computes the pseudoinverse covariance:
        Sigma_theta^dagger = B @ F_tan^dagger @ B.T
    where F_tan = B^T F_X B and F_tan^dagger uses eigendecomposition with _FISHER_EIG_TOL:
      - Eigenvalue > _FISHER_EIG_TOL : 1/lambda  (identifiable direction)
      - Eigenvalue <= _FISHER_EIG_TOL : 0         (null direction)

    An isoform t is marked null (SE = nan) when its tangent-space representation
    b_t = B^T e_t has any projection onto the null space of F_tan (threshold = 0).
    When F_tan is positive definite (fisher_identifiable=True), no null space exists and
    all isoforms have finite SE regardless of the magnitude of diag_var.

    CI_t = [max(0, theta_t - z*SE_t), min(1, theta_t + z*SE_t)]  (nan/nan if SE = nan)

    Returns (se, ci_lo, ci_hi) each of shape (T,).
    """
    T = B.shape[0]
    eigvals, eigvecs = np.linalg.eigh(F_tan)          # ascending, symmetric
    eigvals = np.maximum(eigvals, 0.0)                 # F_tan is PSD by construction; clip fp negatives
    inv_eigvals = np.where(eigvals > _FISHER_EIG_TOL,
                           1.0 / np.maximum(eigvals, _FISHER_EIG_TOL),
                           0.0)   # Moore-Penrose pseudoinverse: 0 for null directions
    # Sigma_theta^dagger = B (B^T F_X B)^dagger B^T
    F_tan_pinv = eigvecs @ np.diag(inv_eigvals) @ eigvecs.T
    Sigma_theta_pinv = B @ F_tan_pinv @ B.T
    diag_var = np.diag(Sigma_theta_pinv)

    # Null isoforms: b_t = B^T e_t lies predominantly in the null space of F_tan.
    # When F_tan is positive definite, its null space is empty and is_null is all-False.
    null_mask = eigvals <= _FISHER_EIG_TOL
    if np.any(null_mask):
        null_vecs = eigvecs[:, null_mask]          # (T-1) x n_null
        BT = B.T                                   # (T-1) x T
        null_proj = np.sum((null_vecs.T @ BT) ** 2, axis=0)   # shape (T,)
        b_norm_sq = np.sum(BT ** 2, axis=0)
        is_null = null_proj / np.maximum(b_norm_sq, 1e-30) > 0
    else:
        is_null = np.zeros(T, dtype=bool)

    se = np.where(is_null, np.nan, np.sqrt(np.maximum(diag_var, 0.0)))

    z = _Z_95 if gamma == 0.05 else float(
        -np.log(gamma / 2) ** 0.5 * np.sqrt(2))   # fallback approx
    safe_se = np.where(is_null, 0.0, se)           # avoid nan propagation in arithmetic
    ci_lo = np.where(is_null, np.nan, np.maximum(0.0, theta_hat - z * safe_se))
    ci_hi = np.where(is_null, np.nan, np.minimum(1.0, theta_hat + z * safe_se))
    return se, ci_lo, ci_hi


# ---------------------------------------------------------------------------
# Per-gene computation
# ---------------------------------------------------------------------------

def compute_gene_identifiability(sr_mds, lr_mds, sr_theoretical_mds=None,
                                  theta_hat_em=None, n_lrs=None, n_srs=None):
    """
    Compute identifiability metrics for a single gene, supporting multiple samples.

    Parameters
    ----------
    sr_mds             : list of dict  SR matrix_dicts, one per sample (can be empty)
    lr_mds             : list of dict  LR matrix_dicts, one per sample (can be empty)
    sr_theoretical_mds : list of dict or None  theoretical SR matrices (same length as
                         sr_mds); None entries fall back to the corresponding sr_md
    theta_hat_em       : np.ndarray or None  isoform relative abundances (sums to 1)
    n_lrs              : list of int  reads per LR sample (scales F_L; default 1 each)
    n_srs              : list of int  reads per SR sample (scales F_S; default 1 each)

    Returns
    -------
    dict with subset of keys 'LR', 'SR', 'Hybrid' depending on available platforms.
    Jacobians from all samples are stacked for k-values; Fisher matrices are summed
    with per-sample N weights.
    """
    if sr_theoretical_mds is None:
        sr_theoretical_mds = [None] * len(sr_mds)
    sr_ref_mds = [th if th is not None else act
                  for th, act in zip(sr_theoretical_mds, sr_mds)]
    if n_lrs is None:
        n_lrs = [1] * len(lr_mds)
    if n_srs is None:
        n_srs = [1] * len(sr_ref_mds)

    has_lr = len(lr_mds) > 0
    has_sr = len(sr_ref_mds) > 0
    if not has_lr and not has_sr:
        return None

    T = (lr_mds[0]['isoform_region_matrix'].shape[1] if has_lr
         else sr_ref_mds[0]['isoform_region_matrix'].shape[1])
    if T <= 1:
        return None

    # ---- k_orig（取第一个样本的 SVD condition number） ----
    result = {}
    if has_lr:
        cond_LR = lr_mds[0].get('condition_number', (math.nan,) * 4)
        result['LR'] = {'k_orig': float(cond_LR[2]) if cond_LR[2] is not None else math.nan}
    if has_sr:
        cond_SR = sr_ref_mds[0].get('condition_number', (math.nan,) * 4)
        result['SR'] = {'k_orig': float(cond_SR[2]) if cond_SR[2] is not None else math.nan}

    if theta_hat_em is None:
        return result

    B = compute_tangent_basis(T)
    _EPS = 1e-12

    # ---- LR：各样本 Jacobian 堆叠，Fisher 加权求和 ----
    J_L_list = []
    F_L = np.zeros((T - 1, T - 1))
    F_L_list = []
    for lr_md, n_lr in zip(lr_mds, n_lrs):
        A_L = lr_md['isoform_region_matrix']
        J_Lk = A_L @ B
        J_L_list.append(J_Lk)
        p_Lk = np.maximum(A_L @ theta_hat_em, _EPS)
        F_Lk = max(n_lr, 1) * (J_Lk.T @ ((1.0 / p_Lk)[:, None] * J_Lk))
        F_L += F_Lk
        F_L_list.append(F_Lk)

    if has_lr:
        # 多样本时额外输出每个样本的单独指标
        if len(lr_mds) > 1:
            for k, (lr_md, J_Lk, F_Lk) in enumerate(zip(lr_mds, J_L_list, F_L_list)):
                cond_k = lr_md.get('condition_number', (math.nan,) * 4)
                k_orig_k = float(cond_k[2]) if cond_k[2] is not None else math.nan
                kv_k, sm_k, smi_k = _metrics_from_Jtan(J_Lk)
                flc_k, flmin_k, fi_k, fd_k = _metrics_from_Fisher(F_Lk)
                result[f'LR_{k+1}'] = {
                    'k_orig': k_orig_k,
                    'k_value': kv_k, 'sigma_max': sm_k, 'sigma_min': smi_k,
                    'fisher_lambda_cond': flc_k, 'fisher_lambda_min': flmin_k,
                    'fisher_identifiable': fi_k, 'fisher_det': fd_k,
                    'se_ci': _compute_se_ci(F_Lk, B, theta_hat_em),
                }
        # 合并 LR 结果
        J_L = np.vstack(J_L_list)
        k_LR, smax_LR, smin_LR = _metrics_from_Jtan(J_L)
        flc_LR, flmin_LR, fi_LR, fd_LR = _metrics_from_Fisher(F_L)
        result['LR'].update({
            'k_value': k_LR, 'sigma_max': smax_LR, 'sigma_min': smin_LR,
            'fisher_lambda_cond': flc_LR, 'fisher_lambda_min': flmin_LR,
            'fisher_identifiable': fi_LR, 'fisher_det': fd_LR,
            'se_ci': _compute_se_ci(F_L, B, theta_hat_em),
        })

    # ---- SR：各样本 Jacobian 堆叠，Fisher 加权求和 ----
    J_S_list = []
    F_S = np.zeros((T - 1, T - 1))
    F_S_list = []
    for sr_ref_md, n_sr in zip(sr_ref_mds, n_srs):
        A_S = sr_ref_md['isoform_region_matrix']
        l_tilde = compute_isoform_eff_len_sr(
            sr_ref_md['isoform_lengths'], sr_ref_md['sr_read_len'])
        Jphi = compute_Jphi(theta_hat_em, l_tilde)
        J_Sk = A_S @ Jphi @ B
        J_S_list.append(J_Sk)
        s = float(l_tilde @ theta_hat_em)
        psi = (l_tilde * theta_hat_em) / s if s > 0 else np.ones(T) / T
        p_Sk = np.maximum(A_S @ psi, _EPS)
        F_Sk = max(n_sr, 1) * (J_Sk.T @ ((1.0 / p_Sk)[:, None] * J_Sk))
        F_S += F_Sk
        F_S_list.append(F_Sk)

    if has_sr:
        # 多样本时额外输出每个样本的单独指标
        if len(sr_ref_mds) > 1:
            for k, (sr_md, J_Sk, F_Sk) in enumerate(zip(sr_ref_mds, J_S_list, F_S_list)):
                cond_k = sr_md.get('condition_number', (math.nan,) * 4)
                k_orig_k = float(cond_k[2]) if cond_k[2] is not None else math.nan
                kv_k, sm_k, smi_k = _metrics_from_Jtan(J_Sk)
                flc_k, fsmin_k, fi_k, fd_k = _metrics_from_Fisher(F_Sk)
                result[f'SR_{k+1}'] = {
                    'k_orig': k_orig_k,
                    'k_value': kv_k, 'sigma_max': sm_k, 'sigma_min': smi_k,
                    'fisher_lambda_cond': flc_k, 'fisher_lambda_min': fsmin_k,
                    'fisher_identifiable': fi_k, 'fisher_det': fd_k,
                    'se_ci': _compute_se_ci(F_Sk, B, theta_hat_em),
                }
        # 合并 SR 结果
        J_S = np.vstack(J_S_list)
        k_SR, smax_SR, smin_SR = _metrics_from_Jtan(J_S)
        flc_SR, fsmin_SR, fi_SR, fd_SR = _metrics_from_Fisher(F_S)
        result['SR'].update({
            'k_value': k_SR, 'sigma_max': smax_SR, 'sigma_min': smin_SR,
            'fisher_lambda_cond': flc_SR, 'fisher_lambda_min': fsmin_SR,
            'fisher_identifiable': fi_SR, 'fisher_det': fd_SR,
            'se_ci': _compute_se_ci(F_S, B, theta_hat_em),
        })

    # ---- Hybrid（只要总样本数 > 1 就计算：含多LR、多SR、或LR+SR）----
    n_total_samples = (len(lr_mds) if has_lr else 0) + (len(sr_mds) if has_sr else 0)
    if n_total_samples > 1:
        all_J = J_L_list + J_S_list
        J_H = np.vstack(all_J)
        k_H, smax_H, smin_H = _metrics_from_Jtan(J_H)
        F_H = F_L + F_S  # F_L/F_S 均从零初始化，无该平台时保持零矩阵
        flc_H, fhmin_H, fi_H, fd_H = _metrics_from_Fisher(F_H)
        result['Hybrid'] = {
            'k_value': k_H, 'sigma_max': smax_H, 'sigma_min': smin_H,
            'fisher_lambda_cond': flc_H, 'fisher_lambda_min': fhmin_H,
            'fisher_identifiable': fi_H, 'fisher_det': fd_H,
            'se_ci': _compute_se_ci(F_H, B, theta_hat_em),
        }

    return result


# ---------------------------------------------------------------------------
# EM result loader
# ---------------------------------------------------------------------------

def _load_em_tpm(output_path):
    """
    Read EM quantification results and return {gene_name: {isoform_name: tpm}}.
    Tries the following files in order (community mode first):
        1. Isoform_abundance.out  -- community mode (columns: Isoform, Gene, TPM, ...)
        2. expression_isoform.out -- regular mode   (columns: Isoform, Gene, Chr, ..., TPM, Alpha)
    Returns None if neither file exists.
    """
    for fname in ('Isoform_abundance.out', 'expression_isoform.out'):
        expr_file = os.path.join(output_path, fname)
        if os.path.isfile(expr_file):
            gene_isoform_tpm = {}
            with open(expr_file, 'r') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    gene = row.get('Gene', '')
                    iso  = row.get('Isoform', '')
                    if not gene or not iso:
                        continue
                    try:
                        tpm = float(row['TPM'])
                    except (ValueError, KeyError):
                        tpm = 0.0
                    gene_isoform_tpm.setdefault(gene, {})[iso] = tpm
            print(f'[INFO] EM TPM loaded from {fname}', flush=True)
            return gene_isoform_tpm
    return None


def _load_gene_read_counts(output_path):
    """
    Load per-gene, per-sample read counts from Isoform_abundance.out.

    For sample k: num_expected_SRs_k = N_{S,k} * theta_t.
    Summing over all isoforms in a gene gives N_{S,k} for that sample.

    Returns {gene_name: (n_lr_list, n_sr_list)} where each list is ordered by
    sample index (matching lr_matrix_input / sr_matrix_input order).
    Returns None if file not found or has no count columns.
    """
    fpath = os.path.join(output_path, 'Isoform_abundance.out')
    if not os.path.isfile(fpath):
        return None

    with open(fpath, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        fieldnames = reader.fieldnames or []
        sr_cols = sorted(c for c in fieldnames if c.startswith('num_expected_SRs_'))
        lr_cols = sorted(c for c in fieldnames if c.startswith('num_expected_LRs_'))
        if not sr_cols and not lr_cols:
            return None
        n_sr_samples = len(sr_cols)
        n_lr_samples = len(lr_cols)
        gene_n_sr = {}   # gene -> [n_s1, n_s2, ...]
        gene_n_lr = {}   # gene -> [n_l1, n_l2, ...]
        for row in reader:
            gene = row.get('Gene', '')
            if not gene:
                continue
            if gene not in gene_n_sr:
                gene_n_sr[gene] = [0.0] * n_sr_samples
                gene_n_lr[gene] = [0.0] * n_lr_samples
            for i, c in enumerate(sr_cols):
                gene_n_sr[gene][i] += float(row.get(c, 0) or 0)
            for i, c in enumerate(lr_cols):
                gene_n_lr[gene][i] += float(row.get(c, 0) or 0)

    all_genes = set(gene_n_sr) | set(gene_n_lr)
    if not all_genes:
        return None
    result = {
        g: (
            [int(round(v)) for v in gene_n_lr.get(g, [])],
            [int(round(v)) for v in gene_n_sr.get(g, [])],
        )
        for g in all_genes
    }
    print(f'[INFO] Gene read counts loaded from Isoform_abundance.out '
          f'({len(result)} genes, {n_lr_samples} LR sample(s), {n_sr_samples} SR sample(s)).',
          flush=True)
    return result


def _build_theta_hat_em(gene_name, isoform_names_indics, gene_isoform_tpm):
    """
    Build theta_hat_em (shape T,) from EM TPM, aligned to isoform_names_indics order.
    Returns None if gene not found or all TPMs are zero.
    """
    iso_tpm = gene_isoform_tpm.get(gene_name)
    if iso_tpm is None:
        return None
    T = len(isoform_names_indics)
    theta = np.zeros(T)
    for iso_name, idx in isoform_names_indics.items():
        theta[idx] = iso_tpm.get(iso_name, 0.0)
    s = theta.sum()
    if s <= 0:
        return None
    return theta / s


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def compute_and_output_identifiability(output_path, sr_matrix_input, lr_matrix_input, sr_theoretical_input=None, threads=1):
    """
    Compute identifiability metrics for every gene and write identifiability.tsv.

    Parameters
    ----------
    output_path      : str  directory to write identifiability.tsv
    sr_matrix_input  : dict or list of dicts  (chr -> gene -> matrix_dict) for SR
    lr_matrix_input  : dict or list of dicts  (chr -> gene -> matrix_dict) for LR
    sr_theoretical_input : dict or None  (chr -> gene -> theoretical matrix_dict for SR)
    
    This function is read-only with respect to the input dicts and does NOT
    affect quantification results.
    """
    print('[INFO] Computing identifiability metrics...', flush=True)

    sr_list = sr_matrix_input if isinstance(sr_matrix_input, list) else [sr_matrix_input]
    lr_list = lr_matrix_input if isinstance(lr_matrix_input, list) else [lr_matrix_input]

    # Use first-sample matrices; treat None as empty dict
    sr_all = sr_list[0] if sr_list[0] is not None else {}
    lr_all = lr_list[0] if lr_list[0] is not None else {}
    sr_theoretical_all = sr_theoretical_input

    # Load EM quantification results if available
    gene_isoform_tpm = _load_em_tpm(output_path)
    has_em = gene_isoform_tpm is not None
    if has_em:
        print('[INFO] expression_isoform.out detected; using EM theta for restricted Jacobian.', flush=True)
    else:
        print('[INFO] expression_isoform.out not found; outputting k_orig only.', flush=True)

    # Load per-gene read counts for Fisher matrix scaling (N_L, N_S)
    gene_read_counts = _load_gene_read_counts(output_path) if has_em else None

    def _fmt(v):
        return 'nan' if math.isnan(v) or math.isinf(v) else f'{v:.6g}'

    def _fse(v):
        return 'null' if (not np.isfinite(v)) else f'{v:.6g}'

    # Collect union of (chr, gene) pairs across all samples of both platforms
    all_gene_keys = set()
    for s_all in sr_list:
        if s_all is None:
            continue
        for chr_name in s_all:
            for gene_name in s_all[chr_name]:
                all_gene_keys.add((chr_name, gene_name))
    for l_all in lr_list:
        if l_all is None:
            continue
        for chr_name in l_all:
            for gene_name in l_all[chr_name]:
                all_gene_keys.add((chr_name, gene_name))

    def _process_one_gene(chr_name, gene_name):
        """Per-gene worker; returns (row, se_ci_rows_gene, is_zero_expr) or None."""
        # 收集每个 LR 样本的 matrix_dict
        lr_mds_gene = [l_all[chr_name][gene_name]
                       for l_all in lr_list
                       if l_all is not None
                       and chr_name in l_all
                       and gene_name in l_all[chr_name]]

        # 收集每个 SR 样本的 matrix_dict 及对应 theoretical md
        sr_mds_gene = []
        sr_theoretical_mds_gene = []
        for s_all in sr_list:
            if s_all is None or chr_name not in s_all or gene_name not in s_all[chr_name]:
                continue
            sr_mds_gene.append(s_all[chr_name][gene_name])
            th = (sr_theoretical_all[chr_name][gene_name]
                  if sr_theoretical_all is not None
                  and chr_name in sr_theoretical_all
                  and gene_name in sr_theoretical_all[chr_name]
                  else None)
            sr_theoretical_mds_gene.append(th)

        ref_md = (lr_mds_gene[0] if lr_mds_gene
                  else sr_mds_gene[0] if sr_mds_gene else None)
        if ref_md is None:
            return None

        theta_hat_em = None
        is_zero_expr = False
        if has_em:
            iso_indics = ref_md.get('isoform_names_indics', {})
            theta_hat_em = _build_theta_hat_em(gene_name, iso_indics, gene_isoform_tpm)
            if theta_hat_em is None and gene_name in gene_isoform_tpm:
                is_zero_expr = True

        counts = gene_read_counts.get(gene_name, ([], [])) if gene_read_counts else ([], [])
        n_lrs_gene = [max(n, 1) for n in counts[0]] if counts[0] else [1] * len(lr_mds_gene)
        n_srs_gene = [max(n, 1) for n in counts[1]] if counts[1] else [1] * len(sr_mds_gene)
        n_lrs_gene = (n_lrs_gene + [1] * len(lr_mds_gene))[:len(lr_mds_gene)]
        n_srs_gene = (n_srs_gene + [1] * len(sr_mds_gene))[:len(sr_mds_gene)]

        try:
            metrics = compute_gene_identifiability(
                sr_mds_gene, lr_mds_gene, sr_theoretical_mds_gene,
                theta_hat_em, n_lrs=n_lrs_gene, n_srs=n_srs_gene)
        except Exception as exc:
            print(f'[WARN] identifiability skipped for {gene_name}: {exc}', flush=True)
            return None

        if metrics is None:
            return None

        ref_for_T = sr_mds_gene[0] if sr_mds_gene else lr_mds_gene[0]
        T = ref_for_T['isoform_region_matrix'].shape[1]

        # 只输出 metrics 中实际存在的平台对应的列
        has_LR     = 'LR'     in metrics
        has_SR     = 'SR'     in metrics
        has_Hybrid = 'Hybrid' in metrics
        # 分样本键：LR_1, LR_2, ... 和 SR_1, SR_2, ...
        lr_per_sample_keys = sorted([k for k in metrics if k.startswith('LR_') and k[3:].isdigit()])
        sr_per_sample_keys = sorted([k for k in metrics if k.startswith('SR_') and k[3:].isdigit()])
        # 只要有分样本键就不再单独输出该平台合并列，合并结果统一用 Hybrid_* 表示
        lr_platform_redundant = bool(lr_per_sample_keys)
        sr_platform_redundant = bool(sr_per_sample_keys)

        row = {'gene': gene_name, 'chr': chr_name, 'n_isoforms': T}
        # 分样本 k_orig 先输出（不依赖 EM）
        for tag_k in lr_per_sample_keys:
            row[f'{tag_k}_k_orig'] = _fmt(metrics[tag_k].get('k_orig', math.nan))
        for tag_k in sr_per_sample_keys:
            row[f'{tag_k}_k_orig'] = _fmt(metrics[tag_k].get('k_orig', math.nan))
        # 平台合并 k_orig（仅在非冗余时输出）
        if has_LR and not lr_platform_redundant:
            row['LR_k_orig'] = _fmt(metrics['LR']['k_orig'])
        if has_SR and not sr_platform_redundant:
            row['SR_k_orig'] = _fmt(metrics['SR']['k_orig'])

        if has_em and theta_hat_em is not None:
            # 分样本 Fisher 指标先输出
            for tag_k in lr_per_sample_keys:
                m_k = metrics[tag_k]
                row.update({
                    f'{tag_k}_rJacobi_sigma_min_pos':      _fmt(m_k.get('sigma_min', math.nan)),
                    f'{tag_k}_rJacobi_k_value':            _fmt(m_k.get('k_value', math.nan)),
                    f'{tag_k}_rfisher_lambda_min':         _fmt(m_k.get('fisher_lambda_min', math.nan)),
                    f'{tag_k}_rfisher_lambda_cond':        _fmt(m_k.get('fisher_lambda_cond', math.nan)),
                    f'{tag_k}_rfisher_det':                _fmt(m_k.get('fisher_det', math.nan)),
                    f'{tag_k}_rfisher_identifiable':       m_k.get('fisher_identifiable', 'nan'),
                })
            for tag_k in sr_per_sample_keys:
                m_k = metrics[tag_k]
                row.update({
                    f'{tag_k}_rJacobi_sigma_min_pos':      _fmt(m_k.get('sigma_min', math.nan)),
                    f'{tag_k}_rJacobi_k_value':            _fmt(m_k.get('k_value', math.nan)),
                    f'{tag_k}_rfisher_lambda_min':         _fmt(m_k.get('fisher_lambda_min', math.nan)),
                    f'{tag_k}_rfisher_lambda_cond':        _fmt(m_k.get('fisher_lambda_cond', math.nan)),
                    f'{tag_k}_rfisher_det':                _fmt(m_k.get('fisher_det', math.nan)),
                    f'{tag_k}_rfisher_identifiable':       m_k.get('fisher_identifiable', 'nan'),
                })
            # 平台合并指标（仅在非冗余时输出）
            if has_LR and not lr_platform_redundant:
                lr_m = metrics['LR']
                row.update({
                    'LR_rJacobi_sigma_min_pos':  _fmt(lr_m.get('sigma_min', math.nan)),
                    'LR_rJacobi_k_value':        _fmt(lr_m.get('k_value', math.nan)),
                    'LR_rfisher_lambda_min':     _fmt(lr_m.get('fisher_lambda_min', math.nan)),
                    'LR_rfisher_lambda_cond':    _fmt(lr_m.get('fisher_lambda_cond', math.nan)),
                    'LR_rfisher_det':            _fmt(lr_m.get('fisher_det', math.nan)),
                    'LR_rfisher_identifiable':   lr_m.get('fisher_identifiable', 'nan'),
                })
            if has_SR and not sr_platform_redundant:
                sr_m = metrics['SR']
                row.update({
                    'SR_rJacobi_sigma_min_pos':  _fmt(sr_m.get('sigma_min', math.nan)),
                    'SR_rJacobi_k_value':        _fmt(sr_m.get('k_value', math.nan)),
                    'SR_rfisher_lambda_min':     _fmt(sr_m.get('fisher_lambda_min', math.nan)),
                    'SR_rfisher_lambda_cond':    _fmt(sr_m.get('fisher_lambda_cond', math.nan)),
                    'SR_rfisher_det':            _fmt(sr_m.get('fisher_det', math.nan)),
                    'SR_rfisher_identifiable':   sr_m.get('fisher_identifiable', 'nan'),
                })
            # Hybrid 放最后
            if has_Hybrid:
                h = metrics['Hybrid']
                row.update({
                    'Hybrid_rJacobi_sigma_min_pos':  _fmt(h.get('sigma_min', math.nan)),
                    'Hybrid_rJacobi_k_value':        _fmt(h.get('k_value', math.nan)),
                    'Hybrid_rfisher_lambda_min':     _fmt(h.get('fisher_lambda_min', math.nan)),
                    'Hybrid_rfisher_lambda_cond':    _fmt(h.get('fisher_lambda_cond', math.nan)),
                    'Hybrid_rfisher_det':            _fmt(h.get('fisher_det', math.nan)),
                    'Hybrid_rfisher_identifiable':   h.get('fisher_identifiable', 'nan'),
                })

        # ---- 收集 per-isoform SE / CI 行 ----
        se_ci_rows_gene = []
        if has_em:
            iso_indics = ref_md.get('isoform_names_indics', {})
            idx_to_name = {v: k for k, v in iso_indics.items()}
            if theta_hat_em is not None:
                active_tags = (
                    [(tag_k, metrics[tag_k].get('se_ci')) for tag_k in lr_per_sample_keys] +
                    [(tag_k, metrics[tag_k].get('se_ci')) for tag_k in sr_per_sample_keys] +
                    ([] if lr_platform_redundant else [('LR', metrics['LR'].get('se_ci'))] if has_LR else []) +
                    ([] if sr_platform_redundant else [('SR', metrics['SR'].get('se_ci'))] if has_SR else []) +
                    ([('Hybrid', metrics['Hybrid'].get('se_ci'))] if has_Hybrid else [])
                )
                for t in range(T):
                    iso_name = idx_to_name.get(t, f'isoform_{t}')
                    iso_row = {
                        'gene':      gene_name,
                        'chr':       chr_name,
                        'isoform':   iso_name,
                        'theta_hat': f'{theta_hat_em[t]:.6g}',
                    }
                    for tag, sc in active_tags:
                        if sc is not None:
                            se, ci_lo, ci_hi = sc
                            iso_row[f'{tag}_SE']    = _fse(se[t])
                            iso_row[f'{tag}_CI_lo'] = _fse(ci_lo[t])
                            iso_row[f'{tag}_CI_hi'] = _fse(ci_hi[t])
                        else:
                            iso_row[f'{tag}_SE']    = 'null'
                            iso_row[f'{tag}_CI_lo'] = 'null'
                            iso_row[f'{tag}_CI_hi'] = 'null'
                    se_ci_rows_gene.append(iso_row)
            else:
                for t in range(T):
                    iso_name = idx_to_name.get(t, f'isoform_{t}')
                    se_ci_rows_gene.append({
                        'gene':      gene_name,
                        'chr':       chr_name,
                        'isoform':   iso_name,
                        'theta_hat': '0' if is_zero_expr else 'null',
                    })

        return row, se_ci_rows_gene, is_zero_expr

    # ---- 并行计算每个基因 ----
    rows = []
    se_ci_rows = []
    zero_expression_gene_keys = set()
    gene_keys_sorted = sorted(all_gene_keys)
    results_dict = {}
    with ThreadPoolExecutor(max_workers=max(1, threads)) as executor:
        future_to_key = {
            executor.submit(_process_one_gene, chr_name, gene_name): (chr_name, gene_name)
            for chr_name, gene_name in gene_keys_sorted
        }
        for future in as_completed(future_to_key):
            key = future_to_key[future]
            result = future.result()
            if result is not None:
                results_dict[key] = result

    # 按原始排序顺序收集结果，保证输出行顺序确定
    for key in gene_keys_sorted:
        if key not in results_dict:
            continue
        row, se_ci_rows_gene, is_zero_expr = results_dict[key]
        if is_zero_expr:
            zero_expression_gene_keys.add(key)
        rows.append(row)
        se_ci_rows.extend(se_ci_rows_gene)

    out_file = os.path.join(output_path, 'identifiability.tsv')
    if rows:
        # Collect all fieldnames from all rows (union, preserving first-seen order)
        fieldnames = list(rows[0].keys())
        for r in rows[1:]:
            for k in r:
                if k not in fieldnames:
                    fieldnames.append(k)
        with open(out_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t', extrasaction='ignore')
            writer.writeheader()
            for r in rows:
                is_zero = (r['chr'], r['gene']) in zero_expression_gene_keys
                writer.writerow({k: r.get(k, 'zero' if is_zero else 'nan') for k in fieldnames})
        print(f'[INFO] Identifiability metrics written to {out_file} ({len(rows)} genes)',
              flush=True)
    else:
        print('[INFO] No genes with sufficient data for identifiability analysis.', flush=True)

    # ---- 写 isoform_SE_CI.tsv ----
    se_ci_file = os.path.join(output_path, 'isoform_SE_CI.tsv')
    if se_ci_rows:
        # 列名从实际行数据中动态收集（保留首次出现顺序）
        se_ci_fields = list(se_ci_rows[0].keys())
        for r in se_ci_rows[1:]:
            for k in r:
                if k not in se_ci_fields:
                    se_ci_fields.append(k)
        with open(se_ci_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=se_ci_fields, delimiter='\t',
                                    extrasaction='ignore')
            writer.writeheader()
            for r in se_ci_rows:
                is_zero = (r['chr'], r['gene']) in zero_expression_gene_keys
                fill = 'zero' if is_zero else 'null'
                writer.writerow({k: r.get(k, fill) for k in se_ci_fields})
        print(f'[INFO] Isoform SE/CI written to {se_ci_file} ({len(se_ci_rows)} isoforms)',
              flush=True)
    else:
        print('[INFO] No isoform SE/CI computed (no invertible Fisher matrix found).', flush=True)

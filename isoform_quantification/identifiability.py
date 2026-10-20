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

Outputs: identifiability.tsv, isoform_SE_CI.tsv  (one row per gene / per isoform)

Optional --filterTPM: additionally recomputes both files after dropping, per gene,
isoforms with theta_hat < 1e-5, writing identifiability_filterTPM.tsv and
isoform_SE_CI_filterTPM.tsv alongside the unfiltered files (which are unaffected).
"""

import math
import os
import csv
import hashlib
import json
import sys
import tempfile
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
    if len(sv) == 0 or sv[0] == 0.0:
        return math.nan, math.nan, math.nan
    tol = sv[0] * max(J_tan.shape) * np.finfo(sv.dtype).eps
    rank = int(np.sum(sv > tol))
    if rank == 0:
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
      fisher_det          : log10(|det(F_tan)|^(1/T-1)) — D-optimality metric (log10 of
                            T-th root of absolute determinant,
                            T = number of isoforms = F_tan.shape[0]+1);
                            nan when singular (det=0) or on numerical sign error.
    """
    eigvals = np.linalg.eigvalsh(F_tan)          # ascending order
    eigvals_clipped = np.maximum(eigvals, 0.0)    # PSD guarantee; clip fp negatives
    lambda_max = float(eigvals_clipped[-1])
    sign, logabsdet = np.linalg.slogdet(F_tan)
    T = F_tan.shape[0] + 1                        # number of isoforms
    if sign <= 0:
        fisher_det = math.nan   # det <= 0: singular or fp artifact
    else:
        fisher_det = float(logabsdet / ((T-1) * math.log(10)))
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
        np.sqrt(-2.0 * np.log(gamma / 2.0)))       # fallback approx
    safe_se = np.where(is_null, 0.0, se)           # avoid nan propagation in arithmetic
    ci_lo = np.where(is_null, np.nan, np.maximum(0.0, theta_hat - z * safe_se))
    ci_hi = np.where(is_null, np.nan, np.minimum(1.0, theta_hat + z * safe_se))
    return se, ci_lo, ci_hi


def _fisher_metrics_and_se_ci(F_tan, B, theta_hat):
    """
    Combined replacement for _metrics_from_Fisher + _compute_se_ci.
    Performs a single eigh decomposition instead of three separate calls
    (eigvalsh, slogdet, eigh), halving the linear-algebra cost per (gene, platform).

    Returns
    -------
    fisher_lambda_cond, lambda_min, is_identifiable, fisher_det : same as _metrics_from_Fisher
    se_ci : (se, ci_lo, ci_hi) tuple,  same as _compute_se_ci
    """
    T_iso = B.shape[0]                                      # number of isoforms

    # Single decomposition — replaces eigvalsh + slogdet + eigh
    eigvals, eigvecs = np.linalg.eigh(F_tan)                # ascending order, symmetric
    eigvals_clipped = np.maximum(eigvals, 0.0)              # PSD by construction; clip fp negatives

    lambda_max = float(eigvals_clipped[-1])

    # log det from eigenvalues (replaces slogdet; singular when any raw eigenvalue <= 0)
    if np.any(eigvals <= 0.0):
        fisher_det = math.nan
    else:
        fisher_det = float(np.sum(np.log(eigvals)) / ((T_iso-1) * math.log(10)))

    thresholded = np.where(eigvals_clipped < _FISHER_EIG_TOL, 0.0, eigvals_clipped)
    lambda_min = float(thresholded[0])
    fisher_lambda_cond = math.log10(lambda_max / lambda_min) if lambda_min > 0 else math.nan
    is_identifiable = bool(lambda_min > 0)

    # SE/CI — reuses eigvals_clipped and eigvecs from the decomposition above
    inv_eigvals = np.where(eigvals_clipped > _FISHER_EIG_TOL,
                           1.0 / np.maximum(eigvals_clipped, _FISHER_EIG_TOL),
                           0.0)
    F_tan_pinv = eigvecs @ np.diag(inv_eigvals) @ eigvecs.T
    Sigma_theta_pinv = B @ F_tan_pinv @ B.T
    diag_var = np.diag(Sigma_theta_pinv)

    null_mask = eigvals_clipped <= _FISHER_EIG_TOL
    if np.any(null_mask):
        null_vecs = eigvecs[:, null_mask]
        BT = B.T
        null_proj = np.sum((null_vecs.T @ BT) ** 2, axis=0)
        b_norm_sq = np.sum(BT ** 2, axis=0)
        is_null = null_proj / np.maximum(b_norm_sq, 1e-30) > 0
    else:
        is_null = np.zeros(T_iso, dtype=bool)

    se = np.where(is_null, np.nan, np.sqrt(np.maximum(diag_var, 0.0)))
    safe_se = np.where(is_null, 0.0, se)
    ci_lo = np.where(is_null, np.nan, np.maximum(0.0, theta_hat - _Z_95 * safe_se))
    ci_hi = np.where(is_null, np.nan, np.minimum(1.0, theta_hat + _Z_95 * safe_se))

    return fisher_lambda_cond, lambda_min, is_identifiable, fisher_det, (se, ci_lo, ci_hi)


# ---------------------------------------------------------------------------
# Per-gene computation
# ---------------------------------------------------------------------------

def compute_gene_identifiability(sr_mds, lr_mds, sr_theoretical_mds=None,
                                  theta_hat_em=None, n_lrs=None, n_srs=None,
                                  _diagnostics=None):
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
    _diagnostics       : dict or None  optional collector for lossless diagnostic arrays;
                         does not change the returned metrics or calculations

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

    result = {}
    if has_lr:
        result['LR'] = {}
    if has_sr:
        result['SR'] = {}

    if theta_hat_em is None:
        return result

    B = compute_tangent_basis(T)
    _EPS = 1e-12
    if _diagnostics is not None:
        ref_md = lr_mds[0] if has_lr else sr_mds[0]
        _diagnostics.update({
            'B': B, 'theta_hat': theta_hat_em,
            'isoform_names': _diagnostic_isoform_names(ref_md),
            'fisher_eig_tol': _FISHER_EIG_TOL, 'probability_floor': _EPS,
            'z_95': _Z_95,
        })

    # ---- LR：各样本计算 Fisher 矩阵，加权求和 ----
    F_L = np.zeros((T - 1, T - 1))
    F_L_list = []
    for lr_md, n_lr in zip(lr_mds, n_lrs):
        A_L = lr_md['isoform_region_matrix']
        J_Lk = A_L @ B
        p_Lk = np.maximum(A_L @ theta_hat_em, _EPS)
        F_Lk = max(n_lr, 1) * (J_Lk.T @ ((1.0 / p_Lk)[:, None] * J_Lk))
        F_L += F_Lk
        F_L_list.append(F_Lk)
        if _diagnostics is not None:
            tag = f'LR_{len(F_L_list)}'
            _diagnostics.update({
                f'{tag}_A': A_L, f'{tag}_J_tan': J_Lk,
                f'{tag}_p_used': p_Lk, f'{tag}_N_used': max(n_lr, 1),
                f'{tag}_F_tan': F_Lk,
                f'{tag}_isoform_names': _diagnostic_isoform_names(lr_md),
            })

    if has_lr:
        # 多样本时额外输出每个样本的单独指标
        if len(lr_mds) > 1:
            for k, F_Lk in enumerate(F_L_list):
                flc_k, flmin_k, fi_k, fd_k, se_ci_k = _fisher_metrics_and_se_ci(F_Lk, B, theta_hat_em)
                result[f'LR_{k+1}'] = {
                    'fisher_lambda_cond': flc_k, 'fisher_lambda_min': flmin_k,
                    'fisher_identifiable': fi_k, 'fisher_det': fd_k,
                    'se_ci': se_ci_k,
                }
        # 合并 LR 结果
        flc_LR, flmin_LR, fi_LR, fd_LR, se_ci_LR = _fisher_metrics_and_se_ci(F_L, B, theta_hat_em)
        result['LR'].update({
            'fisher_lambda_cond': flc_LR, 'fisher_lambda_min': flmin_LR,
            'fisher_identifiable': fi_LR, 'fisher_det': fd_LR,
            'se_ci': se_ci_LR,
        })

    # ---- SR：各样本计算 Fisher 矩阵，加权求和 ----
    F_S = np.zeros((T - 1, T - 1))
    F_S_list = []
    for sr_ref_md, n_sr in zip(sr_ref_mds, n_srs):
        A_S = sr_ref_md['isoform_region_matrix']
        l_tilde = compute_isoform_eff_len_sr(
            sr_ref_md['isoform_lengths'], sr_ref_md['sr_read_len'])
        Jphi = compute_Jphi(theta_hat_em, l_tilde)
        J_Sk = A_S @ Jphi @ B
        s = float(l_tilde @ theta_hat_em)
        psi = (l_tilde * theta_hat_em) / s if s > 0 else np.ones(T) / T
        p_Sk = np.maximum(A_S @ psi, _EPS)
        F_Sk = max(n_sr, 1) * (J_Sk.T @ ((1.0 / p_Sk)[:, None] * J_Sk))
        F_S += F_Sk
        F_S_list.append(F_Sk)
        if _diagnostics is not None:
            tag = f'SR_{len(F_S_list)}'
            _diagnostics.update({
                f'{tag}_A': A_S, f'{tag}_J_tan': J_Sk,
                f'{tag}_p_used': p_Sk, f'{tag}_N_used': max(n_sr, 1),
                f'{tag}_F_tan': F_Sk, f'{tag}_phi': psi,
                f'{tag}_J_phi': Jphi, f'{tag}_effective_lengths': l_tilde,
                f'{tag}_isoform_lengths': sr_ref_md['isoform_lengths'],
                f'{tag}_sr_read_len': sr_ref_md['sr_read_len'],
                f'{tag}_isoform_names': _diagnostic_isoform_names(sr_ref_md),
                f'{tag}_uses_theoretical': sr_theoretical_mds[len(F_S_list)-1] is not None,
            })

    if has_sr:
        # 多样本时额外输出每个样本的单独指标
        if len(sr_ref_mds) > 1:
            for k, F_Sk in enumerate(F_S_list):
                flc_k, fsmin_k, fi_k, fd_k, se_ci_k = _fisher_metrics_and_se_ci(F_Sk, B, theta_hat_em)
                result[f'SR_{k+1}'] = {
                    'fisher_lambda_cond': flc_k, 'fisher_lambda_min': fsmin_k,
                    'fisher_identifiable': fi_k, 'fisher_det': fd_k,
                    'se_ci': se_ci_k,
                }
        # 合并 SR 结果
        flc_SR, fsmin_SR, fi_SR, fd_SR, se_ci_SR = _fisher_metrics_and_se_ci(F_S, B, theta_hat_em)
        result['SR'].update({
            'fisher_lambda_cond': flc_SR, 'fisher_lambda_min': fsmin_SR,
            'fisher_identifiable': fi_SR, 'fisher_det': fd_SR,
            'se_ci': se_ci_SR,
        })

    # ---- Hybrid（只要总样本数 > 1 就计算：含多LR、多SR、或LR+SR）----
    n_total_samples = (len(lr_mds) if has_lr else 0) + (len(sr_mds) if has_sr else 0)
    if n_total_samples > 1:
        F_H = F_L + F_S  # F_L/F_S 均从零初始化，无该平台时保持零矩阵
        if _diagnostics is not None:
            _diagnostics['Hybrid_F_tan'] = F_H
        flc_H, fhmin_H, fi_H, fd_H, se_ci_H = _fisher_metrics_and_se_ci(F_H, B, theta_hat_em)
        result['Hybrid'] = {
            'fisher_lambda_cond': flc_H, 'fisher_lambda_min': fhmin_H,
            'fisher_identifiable': fi_H, 'fisher_det': fd_H,
            'se_ci': se_ci_H,
        }

    if _diagnostics is not None:
        if has_lr:
            _diagnostics['LR_F_tan'] = F_L
        if has_sr:
            _diagnostics['SR_F_tan'] = F_S
    return result


def _diagnostic_isoform_names(md):
    """Return names in actual column-index order, not dictionary insertion order."""
    idx_to_name = {idx: name for name, idx in md.get('isoform_names_indics', {}).items()}
    return np.asarray([idx_to_name.get(t, f'isoform_{t}')
                       for t in range(md['isoform_region_matrix'].shape[1])], dtype=str)


def _write_identifiability_diagnostic(directory, chr_name, gene_name, arrays,
                                     scale_by_N, filtered, filter_tpm_threshold,
                                     lr_source_indices, sr_source_indices):
    """Write arrays without rounding or object pickles; errors propagate to the caller."""
    if not arrays or 'theta_hat' not in arrays:
        return
    payload = dict(arrays)
    payload.update({
        'schema_version': 1, 'gene': gene_name, 'chr': chr_name,
        'scale_by_N': scale_by_N, 'filtered': filtered,
        'filter_tpm_threshold': filter_tpm_threshold,
        'LR_source_sample_indices': np.asarray(lr_source_indices, dtype=np.int64),
        'SR_source_sample_indices': np.asarray(sr_source_indices, dtype=np.int64),
    })
    # Gene identifiers may contain path separators; the payload retains their exact text.
    key = json.dumps([chr_name, gene_name], ensure_ascii=False).encode('utf-8')
    stem = hashlib.sha256(key).hexdigest()
    suffix = '_filterTPM' if filtered else ''
    path = os.path.join(directory, f'{stem}{suffix}.npz')
    with open(path, 'xb') as stream:
        np.savez_compressed(stream, **payload)


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
# --filterTPM: optional low-abundance isoform filtering before identifiability
# ---------------------------------------------------------------------------

_FILTER_TPM_THRESHOLD = 1e-5


def _tpm_keep_mask(theta_hat, threshold=_FILTER_TPM_THRESHOLD):
    """
    Boolean mask selecting isoforms with theta_hat >= threshold.
    Returns None if fewer than 2 isoforms would remain (identifiability is undefined
    for T < 2, so filtering is skipped for that gene rather than dropping it).
    """
    keep_mask = theta_hat >= threshold
    if keep_mask.sum() < 2:
        return None
    return keep_mask


def _filter_md(md, keep_mask):
    """
    Restrict matrix_dict `md` to the isoforms selected by keep_mask (aligned to
    isoform_names_indics values, which are 0..T-1 contiguous column indices).

    theta_hat is NOT renormalised to sum to 1 after filtering: since threshold is on
    the order of 1e-5, the renormalisation error is far below _FISHER_EIG_TOL (1e-5)
    and other numerical tolerances already used in this module, so it is dropped as a
    deliberate, documented approximation rather than silently ignored.
    """
    if md is None:
        return None
    old_items = sorted(md['isoform_names_indics'].items(), key=lambda kv: kv[1])
    kept_names = [name for name, idx in old_items if keep_mask[idx]]
    new_md = {
        'isoform_region_matrix': md['isoform_region_matrix'][:, keep_mask],
        'isoform_names_indics': {name: i for i, name in enumerate(kept_names)},
    }
    if 'isoform_lengths' in md:
        new_md['isoform_lengths'] = np.asarray(md['isoform_lengths'])[keep_mask]
    if 'sr_read_len' in md:
        new_md['sr_read_len'] = md['sr_read_len']
    return new_md


def _fmt(v):
    return 'nan' if math.isnan(v) or math.isinf(v) else f'{v:.6g}'


def _fse(v):
    return 'null' if (not np.isfinite(v)) else f'{v:.6g}'


def _build_metrics_row(gene_name, chr_name, metrics, extra_fields=None):
    """
    Build one identifiability.tsv row from a `metrics` dict (as returned by
    compute_gene_identifiability), shared by the unfiltered and --filterTPM paths.

    Returns (row, has_LR, has_SR, has_Hybrid, lr_per_sample_keys, sr_per_sample_keys,
    lr_redundant, sr_redundant) — the trailing fields are also needed by
    _build_se_ci_rows to know which platform tags are active.
    """
    has_LR, has_SR, has_Hybrid = 'LR' in metrics, 'SR' in metrics, 'Hybrid' in metrics
    lr_per_sample_keys = sorted(k for k in metrics if k.startswith('LR_') and k[3:].isdigit())
    sr_per_sample_keys = sorted(k for k in metrics if k.startswith('SR_') and k[3:].isdigit())
    lr_redundant = bool(lr_per_sample_keys)
    sr_redundant = bool(sr_per_sample_keys)

    row = {'gene': gene_name, 'chr': chr_name}
    if extra_fields:
        row.update(extra_fields)

    def _add(tag, m):
        row.update({
            f'{tag}_rfisher_lambda_min':   _fmt(m.get('fisher_lambda_min', math.nan)),
            f'{tag}_rfisher_lambda_cond':  _fmt(m.get('fisher_lambda_cond', math.nan)),
            f'{tag}_rfisher_det':          _fmt(m.get('fisher_det', math.nan)),
            f'{tag}_rfisher_identifiable': m.get('fisher_identifiable', 'nan'),
        })

    for tag_k in lr_per_sample_keys:
        _add(tag_k, metrics[tag_k])
    for tag_k in sr_per_sample_keys:
        _add(tag_k, metrics[tag_k])
    if has_LR and not lr_redundant:
        _add('LR', metrics['LR'])
    if has_SR and not sr_redundant:
        _add('SR', metrics['SR'])
    if has_Hybrid:
        _add('Hybrid', metrics['Hybrid'])

    return row, has_LR, has_SR, has_Hybrid, lr_per_sample_keys, sr_per_sample_keys, lr_redundant, sr_redundant


def _build_se_ci_rows(gene_name, chr_name, T, theta_hat, idx_to_name, metrics,
                       has_LR, has_SR, has_Hybrid,
                       lr_per_sample_keys, sr_per_sample_keys, lr_redundant, sr_redundant):
    """Build per-isoform SE/CI rows, shared by the unfiltered and --filterTPM paths."""
    active_tags = (
        [(tag_k, metrics[tag_k].get('se_ci')) for tag_k in lr_per_sample_keys] +
        [(tag_k, metrics[tag_k].get('se_ci')) for tag_k in sr_per_sample_keys] +
        ([] if lr_redundant else [('LR', metrics['LR'].get('se_ci'))] if has_LR else []) +
        ([] if sr_redundant else [('SR', metrics['SR'].get('se_ci'))] if has_SR else []) +
        ([('Hybrid', metrics['Hybrid'].get('se_ci'))] if has_Hybrid else [])
    )
    rows = []
    for t in range(T):
        iso_name = idx_to_name.get(t, f'isoform_{t}')
        iso_row = {'gene': gene_name, 'chr': chr_name, 'isoform': iso_name,
                   'theta_hat': f'{theta_hat[t]:.6g}'}
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
        rows.append(iso_row)
    return rows


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def compute_and_output_identifiability(output_path, sr_matrix_input, lr_matrix_input,
                                        sr_theoretical_input=None, threads=1,
                                        _preloaded_tpm=None, _preloaded_counts=None,
                                        scale_by_N=True, filter_tpm=False,
                                        filter_tpm_threshold=_FILTER_TPM_THRESHOLD,
                                        output_identifiability_diagnostics=False):
    """
    Compute identifiability metrics for every gene and write identifiability.tsv.

    Parameters
    ----------
    output_path          : str   directory to write identifiability.tsv
    sr_matrix_input      : dict or list of dicts  (chr -> gene -> matrix_dict) for SR
    lr_matrix_input      : dict or list of dicts  (chr -> gene -> matrix_dict) for LR;
                           pass {} to skip LR entirely.
    sr_theoretical_input : dict or None  (chr -> gene -> theoretical matrix_dict for SR)
    _preloaded_tpm       : dict or None  {gene_name: {isoform_name: abundance}} —
                           when provided, bypasses reading expression_isoform.out /
                           Isoform_abundance.out from disk (used by --identi_data mode).
    _preloaded_counts    : dict or None  {gene_name: ([n_lr_list], [n_sr_list])} —
                           when provided, bypasses reading gene read counts from disk.
    scale_by_N           : bool  if True (default), Fisher matrices are scaled by per-gene
                           read counts N; if False, N=1 for all samples and results are
                           written to identifiability_no_N.tsv / isoform_SE_CI_no_N.tsv.
    filter_tpm            : bool  if True, ADDITIONALLY compute identifiability after
                           dropping isoforms with theta_hat < filter_tpm_threshold
                           (per gene, from the same EM theta used for the unfiltered
                           metrics) and write it to identifiability_filterTPM.tsv /
                           isoform_SE_CI_filterTPM.tsv (with the _no_N suffix applied
                           the same way as the unfiltered files). The unfiltered output
                           is always written unchanged; this only adds extra files.
    filter_tpm_threshold : float  relative-abundance threshold for filter_tpm
                           (default 1e-5). theta_hat is not renormalised after
                           filtering — see _filter_md.
    output_identifiability_diagnostics : bool  default False; write lossless per-gene
                           NPZ files into a new run directory below
                           output_path/identifiability_diagnostics (see IDENTIFIABILITY_DIAGNOSTICS.md).

    This function is read-only with respect to the input dicts and does NOT
    affect quantification results.
    """
    suffix = '' if scale_by_N else '_no_N'
    print(f'[INFO] Computing identifiability metrics (scale_by_N={scale_by_N})...', flush=True)

    sr_list = sr_matrix_input if isinstance(sr_matrix_input, list) else [sr_matrix_input]
    lr_list = lr_matrix_input if isinstance(lr_matrix_input, list) else [lr_matrix_input]

    # Use first-sample matrices; treat None as empty dict
    sr_all = sr_list[0] if sr_list[0] is not None else {}
    lr_all = lr_list[0] if lr_list[0] is not None else {}
    sr_theoretical_all = sr_theoretical_input

    # Load EM quantification results if available
    if _preloaded_tpm is not None:
        gene_isoform_tpm = _preloaded_tpm
        has_em = True
        print('[INFO] Using pre-loaded TPM data from feature_data.', flush=True)
    else:
        gene_isoform_tpm = _load_em_tpm(output_path)
        has_em = gene_isoform_tpm is not None
        if has_em:
            print('[INFO] expression_isoform.out detected; using EM theta for restricted Fisher metrics.', flush=True)
        else:
            print('[INFO] expression_isoform.out not found; Fisher metrics will not be computed.', flush=True)

    # Load per-gene read counts for Fisher matrix scaling (N_L, N_S)
    if not scale_by_N:
        gene_read_counts = None   # N=1 for all samples; skip loading
    elif _preloaded_counts is not None:
        gene_read_counts = _preloaded_counts
    else:
        gene_read_counts = _load_gene_read_counts(output_path) if has_em else None

    diagnostic_dir = None
    if output_identifiability_diagnostics:
        diagnostic_root = os.path.join(output_path, 'identifiability_diagnostics')
        os.makedirs(diagnostic_root, exist_ok=True)
        diagnostic_dir = tempfile.mkdtemp(
            prefix='with_N_' if scale_by_N else 'no_N_', dir=diagnostic_root)
        with open(__file__, 'rb') as source:
            source_hash = hashlib.sha256(source.read()).hexdigest()
        diagnostic_metadata = {
            'schema_version': 1, 'status': 'started',
            'python_version': sys.version, 'numpy_version': np.__version__,
            'argv': sys.argv, 'identifiability_source_sha256': source_hash,
            'scale_by_N': scale_by_N, 'filter_tpm': filter_tpm,
            'filter_tpm_threshold': filter_tpm_threshold,
            'fisher_eig_tol': _FISHER_EIG_TOL,
            'has_abundance_data': has_em,
            'counts_available': gene_read_counts is not None,
        }
        with open(os.path.join(diagnostic_dir, 'metadata.json'), 'x', encoding='utf-8') as stream:
            json.dump(diagnostic_metadata, stream, ensure_ascii=False, indent=2)
        print(f'[INFO] High-precision diagnostics: {diagnostic_dir}', flush=True)

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
        """
        Per-gene worker; returns None, or
        (row, se_ci_rows_gene, filtered_row, filtered_se_ci_rows_gene) where the
        filtered_* entries are None/[] unless filter_tpm is enabled and at least
        one isoform was actually dropped for this gene.
        """
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
        not_in_em = False
        if has_em:
            iso_indics = ref_md.get('isoform_names_indics', {})
            theta_hat_em = _build_theta_hat_em(gene_name, iso_indics, gene_isoform_tpm)
            if theta_hat_em is None:
                if gene_name in gene_isoform_tpm:
                    is_zero_expr = True   # gene present in EM output, all isoform TPM == 0
                else:
                    not_in_em = True      # gene absent from EM output entirely

        counts = gene_read_counts.get(gene_name, ([], [])) if gene_read_counts else ([], [])
        n_lrs_gene = [max(n, 1) for n in counts[0]] if counts[0] else [1] * len(lr_mds_gene)
        n_srs_gene = [max(n, 1) for n in counts[1]] if counts[1] else [1] * len(sr_mds_gene)
        n_lrs_gene = (n_lrs_gene + [1] * len(lr_mds_gene))[:len(lr_mds_gene)]
        n_srs_gene = (n_srs_gene + [1] * len(sr_mds_gene))[:len(sr_mds_gene)]

        diagnostics = {} if diagnostic_dir is not None else None
        try:
            metrics = compute_gene_identifiability(
                sr_mds_gene, lr_mds_gene, sr_theoretical_mds_gene,
                theta_hat_em, n_lrs=n_lrs_gene, n_srs=n_srs_gene,
                _diagnostics=diagnostics)
        except Exception as exc:
            print(f'[WARN] identifiability skipped for {gene_name}: {exc}', flush=True)
            return None

        if metrics is None:
            return None

        if diagnostic_dir is not None:
            lr_source_indices = [k for k, d in enumerate(lr_list, 1)
                                 if d is not None and gene_name in d.get(chr_name, {})]
            sr_source_indices = [k for k, d in enumerate(sr_list, 1)
                                 if d is not None and gene_name in d.get(chr_name, {})]
            _write_identifiability_diagnostic(
                diagnostic_dir, chr_name, gene_name, diagnostics, scale_by_N,
                False, filter_tpm_threshold, lr_source_indices, sr_source_indices)

        ref_for_T = sr_mds_gene[0] if sr_mds_gene else lr_mds_gene[0]
        T = ref_for_T['isoform_region_matrix'].shape[1]

        # 只输出 metrics 中实际存在的平台对应的列
        row = {'gene': gene_name, 'chr': chr_name, 'n_isoforms': T}
        has_LR = has_SR = has_Hybrid = False
        lr_per_sample_keys, sr_per_sample_keys = [], []
        lr_platform_redundant = sr_platform_redundant = False

        if has_em and theta_hat_em is not None:
            (row, has_LR, has_SR, has_Hybrid,
             lr_per_sample_keys, sr_per_sample_keys,
             lr_platform_redundant, sr_platform_redundant) = _build_metrics_row(
                gene_name, chr_name, metrics, extra_fields={'n_isoforms': T})
        elif is_zero_expr or not_in_em:
            has_LR, has_SR = 'LR' in metrics, 'SR' in metrics

            def _mark_zero_expr(tag):
                row[f'{tag}_rfisher_lambda_min']   = 'zero'
                row[f'{tag}_rfisher_lambda_cond']  = 'zero'
                row[f'{tag}_rfisher_det']          = 'zero'
                row[f'{tag}_rfisher_identifiable'] = 'zero'

            if has_LR:
                _mark_zero_expr('LR')
            if has_SR:
                _mark_zero_expr('SR')
            if len(lr_mds_gene) + len(sr_mds_gene) > 1:
                _mark_zero_expr('Hybrid')

        # ---- 收集 per-isoform SE / CI 行 ----
        se_ci_rows_gene = []
        if has_em:
            iso_indics = ref_md.get('isoform_names_indics', {})
            idx_to_name = {v: k for k, v in iso_indics.items()}
            if theta_hat_em is not None:
                se_ci_rows_gene = _build_se_ci_rows(
                    gene_name, chr_name, T, theta_hat_em, idx_to_name, metrics,
                    has_LR, has_SR, has_Hybrid,
                    lr_per_sample_keys, sr_per_sample_keys,
                    lr_platform_redundant, sr_platform_redundant)
            else:
                for t in range(T):
                    iso_name = idx_to_name.get(t, f'isoform_{t}')
                    se_ci_rows_gene.append({
                        'gene':      gene_name,
                        'chr':       chr_name,
                        'isoform':   iso_name,
                        'theta_hat': 'null',
                    })

        # ---- --filterTPM：额外用去除低丰度转录本(theta_hat < filter_tpm_threshold)
        #      后的矩阵重算一遍可识别性，原始（未过滤）结果不受影响 ----
        filtered_row = None
        filtered_se_ci_rows_gene = []
        if filter_tpm and has_em and theta_hat_em is not None:
            keep_mask = _tpm_keep_mask(theta_hat_em, filter_tpm_threshold)
            if keep_mask is not None and not keep_mask.all():
                iso_indics = ref_md.get('isoform_names_indics', {})
                idx_to_name = {v: k for k, v in iso_indics.items()}
                removed_names = [idx_to_name.get(t, f'isoform_{t}')
                                  for t in range(T) if not keep_mask[t]]

                lr_mds_f = [_filter_md(md, keep_mask) for md in lr_mds_gene]
                sr_mds_f = [_filter_md(md, keep_mask) for md in sr_mds_gene]
                sr_theoretical_mds_f = [_filter_md(md, keep_mask) for md in sr_theoretical_mds_gene]
                theta_hat_f = theta_hat_em[keep_mask]
                T_f = int(keep_mask.sum())

                diagnostics_f = {} if diagnostic_dir is not None else None
                try:
                    metrics_f = compute_gene_identifiability(
                        sr_mds_f, lr_mds_f, sr_theoretical_mds_f,
                        theta_hat_f, n_lrs=n_lrs_gene, n_srs=n_srs_gene,
                        _diagnostics=diagnostics_f)
                except Exception as exc:
                    print(f'[WARN] filterTPM identifiability skipped for {gene_name}: {exc}', flush=True)
                    metrics_f = None

                if metrics_f is not None:
                    if diagnostic_dir is not None:
                        diagnostics_f['kept_original_indices'] = np.flatnonzero(keep_mask)
                        _write_identifiability_diagnostic(
                            diagnostic_dir, chr_name, gene_name, diagnostics_f,
                            scale_by_N, True, filter_tpm_threshold,
                            lr_source_indices, sr_source_indices)
                    extra = {
                        'n_isoforms': T_f,
                        'n_isoforms_removed': T - T_f,
                        'removed_isoforms': ';'.join(removed_names),
                    }
                    (filtered_row, has_LR_f, has_SR_f, has_Hybrid_f,
                     lr_keys_f, sr_keys_f, lr_redund_f, sr_redund_f) = _build_metrics_row(
                        gene_name, chr_name, metrics_f, extra_fields=extra)

                    ref_md_f = lr_mds_f[0] if lr_mds_f else sr_mds_f[0]
                    idx_to_name_f = {v: k for k, v in ref_md_f['isoform_names_indics'].items()}
                    filtered_se_ci_rows_gene = _build_se_ci_rows(
                        gene_name, chr_name, T_f, theta_hat_f, idx_to_name_f, metrics_f,
                        has_LR_f, has_SR_f, has_Hybrid_f,
                        lr_keys_f, sr_keys_f, lr_redund_f, sr_redund_f)

        return row, se_ci_rows_gene, filtered_row, filtered_se_ci_rows_gene

    def _write_tsv(rows_, out_file, default_val, label):
        if not rows_:
            print(f'[INFO] No {label} to write.', flush=True)
            return
        fieldnames = list(rows_[0].keys())
        for r in rows_[1:]:
            for k in r:
                if k not in fieldnames:
                    fieldnames.append(k)
        with open(out_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t', extrasaction='ignore')
            writer.writeheader()
            for r in rows_:
                writer.writerow({k: r.get(k, default_val) for k in fieldnames})
        print(f'[INFO] {label} written to {out_file} ({len(rows_)} rows)', flush=True)

    # ---- 并行计算每个基因 ----
    rows = []
    se_ci_rows = []
    filtered_rows = []
    filtered_se_ci_rows = []
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
        row, se_ci_rows_gene, filtered_row, filtered_se_ci_rows_gene = results_dict[key]
        rows.append(row)
        se_ci_rows.extend(se_ci_rows_gene)
        if filtered_row is not None:
            filtered_rows.append(filtered_row)
            filtered_se_ci_rows.extend(filtered_se_ci_rows_gene)

    if not rows:
        print('[INFO] No genes with sufficient data for identifiability analysis.', flush=True)
    _write_tsv(rows, os.path.join(output_path, f'identifiability{suffix}.tsv'),
               'nan', 'Identifiability metrics')
    _write_tsv(se_ci_rows, os.path.join(output_path, f'isoform_SE_CI{suffix}.tsv'),
               'null', 'Isoform SE/CI')

    # ---- --filterTPM：额外写出去除低丰度转录本后的可识别性结果 ----
    # 仅对实际发生了过滤的基因输出行；未触发过滤的基因（theta_hat 全部
    # >= 阈值，或过滤后剩余异构体不足 2 个）不出现在这两个文件里。
    if filter_tpm:
        _write_tsv(filtered_rows, os.path.join(output_path, f'identifiability_filterTPM{suffix}.tsv'),
                   'nan', f'Identifiability metrics (filterTPM, threshold={filter_tpm_threshold:g})')
        _write_tsv(filtered_se_ci_rows, os.path.join(output_path, f'isoform_SE_CI_filterTPM{suffix}.tsv'),
                   'null', f'Isoform SE/CI (filterTPM, threshold={filter_tpm_threshold:g})')

    if diagnostic_dir is not None:
        diagnostic_metadata['status'] = 'complete'
        diagnostic_metadata['genes_in_metrics_output'] = len(rows)
        diagnostic_metadata['npz_files'] = sum(
            name.endswith('.npz') for name in os.listdir(diagnostic_dir))
        with open(os.path.join(diagnostic_dir, 'metadata.json'), 'w', encoding='utf-8') as stream:
            json.dump(diagnostic_metadata, stream, ensure_ascii=False, indent=2)


# ---------------------------------------------------------------------------
# --identi_data mode: load feature_data.tsv
# ---------------------------------------------------------------------------

def load_feature_data(feature_data_path):
    """
    Parse feature_data.tsv (tab-separated) for the --identi_data interface.

    Expected columns (header line is auto-detected and skipped):
        1. transcript_name   – transcript identifier (must match GTF transcript_id)
        2. gene_name         – gene identifier (must match GTF gene_id)
        3. gene_count        – number of reads in the gene region (same value for
                               all rows of the same gene; used as N for Fisher scaling)
        4. relative_abundance – relative abundance of each transcript within its gene
                               (non-negative; values across transcripts of the same gene
                               should sum to 1, but any non-negative scale is accepted
                               since the code re-normalises to sum → 1 per gene)
        5. read_length       – integer read length used to construct matrix A
                               (same for every row; used as SR read length r)

    Returns
    -------
    gene_info : dict
        { gene_name:
            { 'transcripts' : [transcript_name, ...],         # ordered by file appearance
              'count'        : int,                            # gene-level read count
              'abundances'   : {transcript_name: float},      # relative abundance per transcript
              'read_length'  : int                            # SR read length
            }
        }
    global_read_length : int
        The single read-length value shared by all rows (150 if file is empty).
    """
    gene_info = {}
    global_read_length = None

    with open(feature_data_path, 'r') as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 5:
                continue
            # Auto-skip a header row: col-5 must be numeric
            try:
                float(parts[4])
            except ValueError:
                if lineno == 1:
                    continue          # skip header
                raise ValueError(
                    f'[ERROR] feature_data.tsv line {lineno}: '
                    f'column 5 is not numeric ("{parts[4]}")')

            transcript_name = parts[0].strip()
            gene_name       = parts[1].strip()
            try:
                gene_count  = float(parts[2])
            except ValueError:
                gene_count  = 0.0
            try:
                abund       = float(parts[3])
            except ValueError:
                abund       = 0.0
            try:
                read_len    = int(float(parts[4]))
            except ValueError:
                read_len    = 150

            if global_read_length is None:
                global_read_length = read_len

            if gene_name not in gene_info:
                gene_info[gene_name] = {
                    'transcripts': [],
                    'count':       int(round(gene_count)),
                    'abundances':  {},
                    'read_length': read_len,
                }
            gene_info[gene_name]['transcripts'].append(transcript_name)
            gene_info[gene_name]['abundances'][transcript_name] = abund

    if global_read_length is None:
        global_read_length = 150

    print(f'[INFO] feature_data loaded: {len(gene_info)} genes, '
          f'read_length={global_read_length}', flush=True)
    return gene_info, global_read_length

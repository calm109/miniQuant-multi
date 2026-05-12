# miniQuant — Multi-Platform Isoform Quantification

miniQuant_multi is a multi-platform RNA isoform quantification tool that integrates sequencing data generated from different platforms or multiple datasets from the same platform, including long-read (LR) and short-read (SR) sequencing. It uses a mixed Bayesian network framework, with parameter estimation performed via the expectation-maximization (EM) algorithm. miniQuant-multi also provides an identifiability analysis module, TrEESR, which computes identifiability indicators to quantify how well isoforms can be distinguished from one another.

## Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Data Preparation](#data-preparation)
- [Quick Start](#quick-start)
- [Subcommand: `quantify`](#subcommand-quantify)
- [Subcommand: `cal_K_value`](#subcommand-cal_k_value)
- [Running Modes](#running-modes)
- [Output Files](#output-files)
- [Example Workflow](#example-workflow)

---

## Features

- **Multi-platform support**: accepts Oxford Nanopore, PacBio (long reads) and Illumina (short reads) SAM files simultaneously
- **Multi-sample support**: multiple LR and/or SR SAM files can be provided; a community-based parallel EM algorithm resolves each gene community independently
- **Identifiability analysis (TrEESR)**: computes k-values and condition numbers that quantify how distinguishable isoforms are, both from pure annotation structure and from actual read data
- **Adaptive alpha**: a pretrained model (`cDNA-ONT`, `dRNA-ONT`, `cDNA-PacBio`) automatically predicts the optimal LR/SR balance parameter α for each gene community
- **Sample weighting**: supports equal weights (default), user-specified weights, or automatic quality-based weights derived from unique mapping rates

---

## Installation

### Option 1: Conda environment (recommended)

```bash
conda env create -f environment.yml
conda activate miniQuant
```

### Option 2: pip

```bash
pip install -r requirements.txt
```

### Option 3: Singularity container

```bash
singularity pull miniQuant.sif docker://your_registry/miniQuant:latest
singularity exec miniQuant.sif python /path/to/isoform_quantification/main.py quantify ...
```

---

## Data Preparation

### Long reads — genome alignment with minimap2

```bash
# cDNA / ONT direct RNA
minimap2 -ax splice -uf --secondary=no \
    reference_genome.fa long_reads.fastq.gz > LR.sam

# PacBio Iso-Seq
minimap2 -ax splice:hq --secondary=no \
    reference_genome.fa long_reads.fastq.gz > LR.sam
```

### Short reads — genome alignment with STAR or bowtie2

```bash
# STAR (genome alignment, recommended)
STAR --runMode genomeGenerate --genomeDir star_index \
     --genomeFastaFiles reference_genome.fa \
     --sjdbGTFfile annotation.gtf

STAR --genomeDir star_index \
     --readFilesIn reads_R1.fastq.gz reads_R2.fastq.gz \
     --readFilesCommand zcat \
     --outSAMtype BAM SortedByCoordinate \
     --outFileNamePrefix SR_

# bowtie2 (transcriptome alignment, also supported)
bowtie2-build transcriptome.fa bowtie2_index
bowtie2 -x bowtie2_index -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -S SR.sam
```

> **Note**: The SAM/BAM file must be coordinate-sorted or read-name sorted. Both genome-aligned and transcriptome-aligned SAM files are supported; miniQuant auto-detects the alignment mode.

---

## Quick Start

```bash
cd isoform_quantification

# Hybrid mode (LR + SR, single sample each)
python main.py quantify \
    -gtf annotation.gtf \
    -lrsam LR.sam \
    -srsam SR.sam \
    -o output/ \
    --EM_choice hybrid \
    --pretrained_model_path cDNA-ONT

# LR-only
python main.py quantify \
    -gtf annotation.gtf \
    -lrsam LR.sam \
    -o output/

# SR-only
python main.py quantify \
    -gtf annotation.gtf \
    -srsam SR.sam \
    -o output/ \
    --EM_choice SR

# Identifiability analysis only
python main.py cal_K_value \
    -gtf annotation.gtf \
    -lrsam LR.sam \
    -srsam SR.sam \
    -o output/
```

---

## Subcommand: `quantify`

```
python main.py quantify -gtf <GTF> -o <OUTPUT> [options]
```

### Required arguments

| Argument | Description |
|---|---|
| `-gtf` / `--gtf_annotation_path` | Path to GTF annotation file |
| `-o` / `--output_path` | Output directory |

### Input data (at least one required)

| Argument | Description |
|---|---|
| `-lrsam` / `--long_read_sam_path` | Path(s) to long-read SAM file(s). Multiple files: `-lrsam LR1.sam LR2.sam` |
| `-srsam` / `--short_read_sam_path` | Path(s) to short-read SAM file(s). Multiple files: `-srsam SR1.sam SR2.sam` |

### Key optional arguments

| Argument | Default | Description |
|---|---|---|
| `--EM_choice` | `LR` | Quantification mode: `LR` (long-read only), `SR` (short-read only), `hybrid` (LR + SR). Automatically set to `hybrid` when both are provided. |
| `--pretrained_model_path` | `cDNA-ONT` | Pretrained model for adaptive α: `cDNA-ONT`, `dRNA-ONT`, `cDNA-PacBio`, or a custom path |
| `--alpha` | `adaptive` | LR/SR balance (0=SR only, 1=LR only). `adaptive` uses the pretrained model per gene community |
| `-t` / `--threads` | `1` | Number of threads |
| `--filtering` | `False` | Filter very short long reads (`True`/`False`) |
| `--EM_SR_num_iters` | `200` | Number of SR EM iterations |
| `--isoform_start_end_site_tolerance` | `20` | Tolerance (bp) for matching LR start/end sites to isoform boundaries |
| `--junction_site_tolerance` | `5` | Tolerance (bp) for matching splice junction sites |
| `--LR_cond_prob_calc` | `form_2` | LR conditional probability model: `form_1` or `form_2` |
| `--add_full_length_region` | `all` | Whether to keep zero-read-count regions in LR matrix A: `all`, `nonfullrank`, `none` |
| `--normalize_lr_A` | `True` | Column-normalize LR design matrix A |
| `--normalize_sr_A` | `True` | Column-normalize SR design matrix A |
| `--sr_region_selection` | `read_length` | SR region filtering: `read_length` (filter by read length coverage), `real_data` (filter by actual reads) |
| `--kde_lr` | `False` | Train a KDE model from LR data to weight the LR A matrix |
| `--keep_kde_lr` | `False` | Keep the trained KDE model file after the run |
| `--output_matrix_info` | `False` | Output per-gene A matrices and b vectors to `matrix_info/` |

### Multi-sample weighting arguments

| Argument | Default | Description |
|---|---|---|
| `--lr_weights` | (equal) | Per-sample weights for LR files (same order as `-lrsam`), e.g. `--lr_weights 0.3 0.3` |
| `--sr_weights` | (equal) | Per-sample weights for SR files (same order as `-srsam`) |
| `--use_quality_weights` | `False` | Automatically weight samples by unique mapping rate |
| `--normalize_q` | `False` | Normalize each sample's q by its sum before weighting in multi-sample E-step |

> When neither `--lr_weights`/`--sr_weights` nor `--use_quality_weights` is set, all samples receive equal weight.  
> When `--lr_weights`/`--sr_weights` are given, they take precedence over `--use_quality_weights` for the respective group.

---

## Subcommand: `cal_K_value`

Computes identifiability metrics (k-values, condition numbers) for all genes without running quantification.

```
python main.py cal_K_value -gtf <GTF> -o <OUTPUT> [options]
```

### Required arguments

| Argument | Description |
|---|---|
| `-gtf` / `--gtf_annotation_path` | Path to GTF annotation file |
| `-o` / `--output_path` | Output directory |

### Key optional arguments

| Argument | Default | Description |
|---|---|---|
| `-lrsam` / `--long_read_sam_path` | — | Long-read SAM file(s) (used to build data-driven LR A matrix) |
| `-srsam` / `--short_read_sam_path` | — | Short-read SAM file(s) (used to build data-driven SR A matrix) |
| `-t` / `--threads` | `1` | Number of threads |
| `--sr_region_selection` | `read_length` | SR region selection mode: `read_length` (filter by actual read length), `real_data` (filter by observed reads) |
| `--lr_region_selection` | `read_length` | LR region selection mode: `read_length` (filter by global median LR read length), `real_data` (use all observed regions) |
| `--add_full_length_region` | `nonfullrank` | Whether to add zero-count full-length regions to LR A matrix: `all`, `nonfullrank`, `none` |
| `--normalize_lr_A` | `True` | Column-normalize LR A matrix |
| `--normalize_sr_A` | `True` | Column-normalize SR A matrix |
| `--same_struc_isoform_handling` | `merge` | How to handle isoforms with identical exon structure: `merge`, `keep` |
| `--keep_sr_exon_region` | `nonfullrank` | Whether to keep zero-count exon regions in SR real_data mode: `nonfullrank`, `all`, `none` |
| `--output_matrix_info` | `False` | Output per-gene A matrices to `matrix_info/` |
| `--kde_lr` | `True` | Train KDE model from LR data for LR A matrix weighting |
| `--singular_values_tol` | `0` | Tolerance for treating small singular values as zero |

---

## Running Modes

### Routing logic

The `quantify` subcommand automatically selects the computation path:

| Condition | Mode | Algorithm |
|---|---|---|
| Only SR provided (any number) | SR / multi-SR | `EM_hybrid_multi` with empty LR placeholder |
| Only LR provided, single file | LR-only | `EM_hybrid` (LIQA-modified) |
| Only LR provided, multiple files | multi-LR | `EM_hybrid_multi` |
| Both LR and SR provided | Hybrid | `EM_hybrid_multi` |

**`is_multi` flag**: routing to `EM_hybrid_multi` (community-based parallel EM) occurs when:
- More than one SR file is provided, OR
- More than one LR file is provided, OR
- Any SR file is provided (even a single SR file uses the multi-sample path)

### Community-based EM

For multi-sample mode, reads and isoforms form a bipartite graph; connected components (gene communities) are solved independently in parallel. This allows efficient parallel computation and reduces interference between unrelated genes.

### Alpha (LR/SR balance)

- `alpha = 1.0`: LR-only quantification
- `alpha = 0.0`: SR-only quantification
- `0 < alpha < 1`: hybrid mode, weighting the contribution of LR vs SR
- `adaptive` (default): a pretrained XGBoost model predicts optimal α per gene community based on data characteristics

---

## Output Files

### `quantify` output

| File | Description |
|---|---|
| `Isoform_abundance.out` | Main result: isoform expression (TPM and read count) for each sample |
| `matrix_info/` | Per-gene A matrix and b vector files (when `--output_matrix_info True`) |

**`Isoform_abundance.out` columns**:

| Column | Description |
|---|---|
| `isoform_id` | Isoform transcript ID |
| `gene_id` | Gene ID |
| `TPM` | Transcripts per million |
| `read_count` | Estimated read count |
| `length` | Transcript length (bp) |

### `cal_K_value` output

| File | Description |
|---|---|
| `kvalues.out` | Per-gene k-values: k_orig (pure annotation structure), k_SR, k_LR (data-driven), k_orig_tangent, k_tangent (restricted Jacobian-based) |
| `kvalues_gene.out` | Gene-level summary (max k across isoforms) |
| `kvalues_isoform.out` | Isoform-level k-value breakdown |
| `identifiability.tsv` | Full identifiability metrics table with condition numbers and matrix ranks |

**k-value interpretation**:

| k-value | Meaning |
|---|---|
| `k = 1` | Perfectly identifiable (orthogonal isoform signatures) |
| `1 < k < 10` | Well identifiable |
| `k > 100` | Poorly identifiable; isoform expression estimates unreliable |
| `k = inf` / `NaN` | Not identifiable (rank-deficient matrix) |

> Genes with empty region sets (e.g., single-exon genes) appear in `kvalues.out` but are excluded from `identifiability.tsv`.

---

## Example Workflow

### Step 1: Quantify, then compute identifiability

Use `run_miniQuant_quantify_then_kvalue.sh` as a template (adapts SR/LR SAM lists and output directory):

```bash
# Edit the script to set your paths
vim run_miniQuant_quantify_then_kvalue.sh

# Submit to cluster (SLURM) or run locally
bash run_miniQuant_quantify_then_kvalue.sh
```

The script:
1. Runs `quantify` (EM-based isoform abundance estimation)
2. Runs `cal_K_value` using the default `read_length` region selection for identifiability analysis

### Step 2: Identifiability analysis only

```bash
python isoform_quantification/main.py cal_K_value \
    -gtf annotation.gtf \
    -lrsam LR1.sam LR2.sam \
    -srsam SR1.sam SR2.sam \
    -o output/ \
    -t 8 \
    --output_matrix_info True
```

### Multi-sample hybrid quantification

```bash
python isoform_quantification/main.py quantify \
    -gtf annotation.gtf \
    -lrsam LR_rep1.sam LR_rep2.sam \
    -srsam SR_rep1.sam SR_rep2.sam \
    -o output/ \
    -t 16 \
    --EM_choice hybrid \
    --pretrained_model_path cDNA-ONT \
    --alpha adaptive \
    --use_quality_weights \
    --kde_lr \
    --output_matrix_info True
```

---

## Dependencies

| Package | Version |
|---|---|
| numpy | ≥ 1.17 |
| scipy | ≥ 1.5 |
| pandas | ≥ 1.1 |
| pysam | ≥ 0.21 |
| scikit-learn | ≥ 1.0 |
| xgboost | ≥ 2.0 |
| intervaltree | ≥ 3.1 |
| Cython | ≥ 0.29 |

See `requirements.txt` for exact pinned versions.

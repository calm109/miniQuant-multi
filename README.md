# miniQuant-multi — Multi-Platform Isoform Quantification

miniQuant-multi is a multi-platform RNA isoform quantification tool that integrates sequencing data generated from different platforms or multiple datasets from the same platform, including long-read (LR) and short-read (SR) sequencing. It uses a mixed Bayesian network framework, with parameter estimation performed via the expectation-maximization (EM) algorithm. miniQuant-multi also provides an identifiability analysis function that computes identifiability indicators to quantify how well isoforms can be distinguished from one another.

## Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Data Preparation](#data-preparation)
- [Quick Start](#quick-start)
- [Subcommand: `quantify`](#subcommand-quantify)
- [Subcommand: `cal_K_value`](#subcommand-cal_k_value)
- [Running Modes](#running-modes)
- [Output Files](#output-files)

---

## Features

- **Flexible input support**: accepts Illumina, Oxford Nanopore, and PacBio SAM files, allowing single-platform, multi-platform, or multi-sample analyses.
- **Joint isoform quantification**: integrates multiple LR and/or SR datasets using a mixed Bayesian network with EM-based parameter estimation.
- **Community-based parallel EM**: partitions genes into independent communities and resolves each community separately for efficient inference.
- **Flexible sample weighting**: supports equal (default), user-defined, or quality-based sample weights estimated from unique mapping rates.
- **Isoform identifiability analysis**: uses `TrEESR.py` and `identifiability.py` to compute identifiability indicators that measure isoform distinguishability from annotation structure and read data.

---

## Installation

### Create conda environment (recommended)

```bash
conda env create -f environment.yml
conda activate miniQuant-multi
```
---

## Data Preparation

* Install `minimap2` (v2.24) and `bowtie2` (v2.4.1).

### Long-read alignment data mapped to the reference genome in SAM format
* Use `minimap2` to map long reads data.
```bash
# cDNA / ONT direct RNA
minimap2 -t 3 -ax splice -uf --secondary=no \
    reference_genome.fa long_reads.fastq.gz > LR.sam

# PacBio Iso-Seq
minimap2 -t 3 -ax splice:hq --secondary=no \
    reference_genome.fa long_reads.fastq.gz > LR.sam
```

### Short-read alignment data mapped to the reference transcriptome in SAM format.  
* Use `Bowtie2` to map short reads data.
```bash
# bowtie2 (transcriptome alignment)
bowtie2-build -f transcriptome.fa bowtie2_index

bowtie2 -q --phred33 --sensitive --dpad 0 --gbar 99999999 --mp 1,1 --np 1 --score-min L,0,-0.1 -I 1 -X 1000 --no-mixed --no-discordant -p 10 -k 200 \
-x bowtie2_index -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz > SR.sam
```

---

## Quick Start

```bash
# Quantify analysis
MINIQUANT_DIR="/path/miniQuant/isoform_quantification"
GTF="/path/annotation.gtf"
OUTPUT_DIR="/path/output"
THREADS=8
SR_SAMS=("/path/to/SR1.sam" "/path/to/SR2.sam")
LR_SAMS=("/path/to/LR1.sam" "/path/to/LR2.sam")

[[ ${#LR_SAMS[@]} -gt 0 ]] && LR_ARGS="-lrsam ${LR_SAMS[*]}" || LR_ARGS=""
[[ ${#SR_SAMS[@]} -gt 0 ]] && SR_ARGS="-srsam ${SR_SAMS[*]}" || SR_ARGS=""
LR_ARGS=""
if [[ ${#LR_SAMS[@]} -gt 0 ]]; then
    LR_ARGS="-lrsam ${LR_SAMS[*]}"
fi
SR_ARGS=""
if [[ ${#SR_SAMS[@]} -gt 0 ]]; then
    SR_ARGS="-srsam ${SR_SAMS[*]}"
fi
mkdir -p "$OUTPUT_DIR"

python "$MINIQUANT_DIR/main.py" quantify \
    -gtf "$GTF" \
    -o "$OUTPUT_DIR" \
    -t "$THREADS" \
    $LR_ARGS \
    $SR_ARGS

# Identifiability analysis
python "$MINIQUANT_DIR/main.py" cal_K_value \
    -gtf "$GTF" \
    -o "$OUTPUT_DIR" \
    -t "$THREADS" \
    $SR_ARGS \
    $LR_ARGS
```

---

## Subcommand: `quantify`

```
python main.py quantify -gtf <GTF> -o <OUTPUT> [options]
```

### Required arguments

| Argument | Description |
|---|---|
| -gtf / --gtf_annotation_path | Path to GTF annotation file |
| -o / --output_path | Output directory |

### Input data (at least one required)

| Argument | Description |
|---|---|
| `-lrsam` / `--long_read_sam_path` | Path(s) to long-read SAM file(s). Multiple files: `-lrsam LR1.sam LR2.sam` |
| `-srsam` / `--short_read_sam_path` | Path(s) to short-read SAM file(s). Multiple files: `-srsam SR1.sam SR2.sam` |

### Key optional arguments

| Argument | Default | Description |
|---|---|---|
| `-t` / `--threads` | `1` | Number of threads |
| `--EM_SR_num_iters` | `200` | Number of SR EM iterations |
| `--isoform_start_end_site_tolerance` | `20` | Tolerance (bp) for matching LR start/end sites to isoform boundaries |
| `--junction_site_tolerance` | `5` | Tolerance (bp) for matching splice junction sites |
| `--lr_weights` | equal | Per-sample weights for LR files (same order as `-lrsam`), e.g. `--lr_weights 0.3 0.3`. By default all samples (LR + SR) share equal weights that sum to 1. |
| `--sr_weights` | equal | Per-sample weights for SR files (same order as `-srsam`), e.g. `--sr_weights 0.4`. By default all samples (LR + SR) share equal weights that sum to 1. |
---

## Subcommand: `cal_K_value`

Computes identifiability metrics (least eigenvalue, condition numbers, standard deviation, and confidence interval) for all genes after running quantification.

```
python main.py cal_K_value -gtf <GTF> -o <OUTPUT> [options]
```

### Required arguments

| Argument | Description |
|---|---|
| `-gtf` / `--gtf_annotation_path` | Path to GTF annotation file |
| `-o` / `--output_path` | Output directory |

### Input data (at least one required)

| Argument | Description |
|---|---|
| `-lrsam` / `--long_read_sam_path` | Long-read SAM file(s) (used to build data-driven LR A matrix) |
| `-srsam` / `--short_read_sam_path` | Short-read SAM file(s) (used to build data-driven SR A matrix) |

### Key optional arguments

| Argument | Default | Description |
|---|---|---|
| `-t` / `--threads` | `1` | Number of threads |
| `--sr_region_selection` | `read_length` | SR region selection mode: `read_length` (filter by actual read length), `real_data` (filter by observed reads) |
| `--lr_region_selection` | `read_length` | LR region selection mode: `read_length` (filter by global median LR read length), `real_data` (use all observed regions) |
| `--add_full_length_region` | `nonfullrank` | Whether to add zero-count full-length regions to LR A matrix: `all`, `nonfullrank`, `none` |
| `--keep_sr_exon_region` | `nonfullrank` | Whether to keep zero-count exon regions in SR real_data mode: `nonfullrank`, `all`, `none` |
| `--singular_values_tol` | `0` | Tolerance for treating small singular values as zero |
| `--add_full_length_region` | `all` | Whether to keep zero-read-count regions in LR matrix A: `all`, `nonfullrank`, `none` |
| `--normalize_lr_A` | `True` | Column-normalize LR design matrix A |
| `--normalize_sr_A` | `True` | Column-normalize SR design matrix A |
| `--output_matrix_info` | `False` | Output per-gene A matrices and b vectors to `matrix_info/` |

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


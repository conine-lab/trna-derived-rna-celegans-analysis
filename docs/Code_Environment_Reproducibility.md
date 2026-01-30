# Code & Environment Reproducibility

This repository contains analysis and plotting code supporting the manuscript figures.  
The instructions below describe how to recreate the exact software environment used to run the scripts and verify that everything works as expected.

---

## What is this file?

This file is written in **Markdown (`.md`)**, a plain-text format that can be opened **anywhere**:

- Any text editor (TextEdit, Notepad, VS Code)
- GitHub / GitLab (renders automatically)
- Zenodo (renders automatically)
- Positron / RStudio / Jupyter
- Word processors (via copy–paste)

No special software is required.

---

## Tested environment

- **Python:** 3.9.2  
- **OS:** macOS (also tested on Linux-compatible environments)

---

## Required packages (pinned versions)

The scripts were developed and tested using the following versions:

```text
pandas==2.3.2
numpy==2.0.2
matplotlib==3.9.4
seaborn==0.13.2
biopython==1.85
plotnine==0.13.6
mizani==0.11.4
statsmodels==0.14.6
patsy==1.0.2
scipy==1.13.1
```

These versions are listed in `requirements.txt`.

---

## External tools

Some scripts rely on external command-line tools that are **not installed via pip**:

- **samtools** (≥ 1.10) — BAM filtering and indexing  
- **deeptools** (`bamCoverage`) — per-base coverage computation  

These tools must be installed and available on the system `PATH`.

---

## Create a clean virtual environment

Using the Python 3.9 interpreter:

```bash
/Library/Frameworks/Python.framework/Versions/3.9/bin/python3 -m venv trf_zenodo_env
source trf_zenodo_env/bin/activate
```

Upgrade packaging tools:

```bash
python -m pip install --upgrade pip setuptools wheel
```

Install dependencies:

```bash
pip install -r requirements.txt
```

---

## Verify the environment

Confirm package versions:

```bash
python -c "import importlib.metadata as md; pkgs=['pandas','numpy','matplotlib','seaborn','biopython','plotnine','mizani','statsmodels','patsy','scipy']; print({p: md.version(p) for p in pkgs})"
```

Test imports:

```bash
python -c "import pandas, numpy, matplotlib, seaborn, Bio, plotnine, scipy, statsmodels, mizani, patsy; print('all imports ok')"
```

---

## Script run order

Recommended execution order:

1. `make_samplesheet_template.py`  
2. `fiveprime_threeprime_bar.py`  
3. `read_length.py`  
4. `dynamic_diverging_bar.py`  
5. `filter_bams_and_coverage.sh`  

---

## Scripts

### `make_samplesheet_template.py`

Generates a **sample sheet template** from BAM filenames.

**Inputs**

- BAM files specified via `--bam-glob`

**Outputs**

A CSV template with columns:

- `bam` – path to BAM file  
- `genotype` – genotype or condition (may be inferred; otherwise blank)  
- `rep` – replicate identifier (may be inferred; otherwise blank)  

Example:

```bash
python make_samplesheet_template.py --bam-glob "*.bam" --out samplesheet_template.csv
```

If genotype or replicate values cannot be inferred from filenames, the corresponding
fields are left blank and should be filled manually before downstream analysis.

---

### `fiveprime_threeprime_bar.py`

Splits tRNA alignments into 5′ and 3′ categories, summarizes counts using a sample sheet,
and optionally generates stacked percentage bar plots.

**Inputs**

- Sample sheet CSV with columns: `bam`, `genotype`, `rep`
- BAM files referenced in the sample sheet

**Outputs**

- Summary CSV of 5′ and 3′ counts per genotype and replicate  
- SVG plot comparing 5′ and 3′ fractions (unless `--no-plot` is specified)

Example (help):

```bash
python 5prime_3prime_bar.py --samples-csv samples_template.csv
```

---

### `read_length.py`

Processes FASTQ files to compute read-length and 5′ nucleotide statistics.

**Inputs**

- FASTQ files named `PREFIX_1.fastq`, `PREFIX_2.fastq`, …  
- Length range specified via `--length-range MIN MAX`

**Outputs**

- Per-replicate TSV summary files  
- Averaged TSV summary file 
- Vectorized plot showing average read length and 5′ nucleotide composition  

Example:

```bash
python read_length.py -p SAMPLE_PREFIX -l 18 40 -r 4
```

---

### `diverging_bar.py`

Generates stacked **diverging horizontal bar plots** comparing two samples across
length bins, stratified by starting nucleotide.

**Inputs**

- Left input TSV/CSV (`--left`)  
- Right input TSV/CSV (`--right`)  

Both input files must contain columns specifying length, starting nucleotide, and
a numeric value (e.g. average TPM).

**Output**

- PDF or SVG figure specified via `--out`

Example:

```bash
python dynamic_diverging_bar.py \
  --left tRNA@ce_maleAVG.tsv \
  --right tRNA@ce_spermAVG.tsv \
  --left-label "Male (-)" \
  --right-label "Sperm (+)" \
  --length-min 18 \
  --length-max 40 \
  --ytick-every 1 \
  --out tRNA@male_spermAVG1040.pdf \
  --title "Small RNA Size Distribution"

```

---

### `filter_bams_and_coverage.sh`

Filters BAM files by reference name, indexes the filtered BAMs, and computes
per-base coverage.

**Inputs**

- Directory containing BAM files  
- Reference name or pattern used for filtering  

**Outputs**

- Filtered BAM files (`*_filtered.bam`)  
- BAM index files (`*_filtered.bai`)  
- Per-base coverage files (`*_summary_unnorm.tsv`)  

**Coverage output format**

Per-base coverage files are written in **bedGraph format** and therefore **do not
include a header line**.

Columns are, in order:

1. `chrom` – reference name  
2. `start` – 0-based start coordinate  
3. `end` – end coordinate  
4. `coverage` – per-base coverage (unnormalized)

This format is compatible with standard genome tools (e.g. IGV, UCSC utilities).
When loading these files into analysis software (e.g. pandas), column names
should be assigned explicitly.

Example:

```bash
./filter_bams_and_coverage.sh tRNA/sense tRNA/sense/filtered/GluCTC GluCTC
```

---

## Notes
- Final figure layout and annotation were performed using standard Prism software.
- When relevant, normalization was performed manually prior to plotting, as described in the manuscript.
- The scripts reproduce the underlying data processing and plotting logic.
- Newer versions of bamCoverage have included functionality to expand intervals that have zero coverage.

---

## Citation

If you use this code, please cite the associated manuscript.

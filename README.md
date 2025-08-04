# STRiker

A Python tool for detecting short tandem repeats (STRs), analyzing motif patterns, and visualizing repeat number distributions in long read sequencing data (e.g., BAM files from Oxford Nanopore sequencing).

---

## 📦 Features

- Detects STR motifs from reference sequence (e.g., GRCh38) and aligned BAM files 
- Generates:
  - PDF reports with:
    - Visualization of STR motifs pattern for each region
    - Kernel Density Estimation (KDE) plots for read lengths for each region
    - Histograms of repeat numbers for each 
  - Read length KDE plots
  - Repeat number histograms

---

## 🚀 Quick Start

### 1. Clone the repository

```bash
git clone https://github.com/BaeLab/STRiker
cd STRiker
```

### 2. Set up environment (conda recommended)


```bash
conda env create  --name STRiker -f environment.yaml
conda activate STRiker
```


---

## 📂 Input Format

### 1. `input.bam` (required)
- BAM file containing aligned long reads. BAM files should be indexed and sorted.
- Aligned and indexed BAM file (`.bai` required)

### 2. `loci.csv` (required)
```csv
gene,chr,start,end,known_motif,pathogenic_expansion_number
HTT,chr4,3074876,3074941,CAG,35
ATXN8,chr13,70139383,70139428,CAG/TAG,73
```
Loci with multiple known motifs can be separated by `/` in the `known_motif` column.


---

## 🛠 Change parameters
You can change the motif-finding parameters by modifying the `__init__.py` file in the `config` directory. \
The parameters include:
- `min_repeat_length`: Minimum length of the repeat motif to be considered




## 🧬 Usage

<!-- ```bash
python striker.py \
    --bam input.bam \
    --regions STR_regions.json \
    --output results/
``` -->
```bash
python STRiker.py <loci.csv> <reference.fasta> <bam_file>
```



<!-- Arguments:
- `--bam`: Input BAM file path
- `--regions`: JSON file with STR region metadata
- `--output`: Folder where PDF and summary files will be saved
- `--read-threshold`: (Optional) Minimum number of valid reads per gene to be included in plot -->

---

## 🖼 Output

Striker generates two directories in the input folder:
- `gene_panel_output/`: Contains PDF reports and summary files for each gene
- `motif_results/`: Contains consolidated motif counts and summary statistics in `xlsx` format
---




## 📄 License

This project is licensed under CC BY-NC-SA 4.0 - Non-commercial use only.

### Commercial Use
For commercial licensing inquiries, please contact bbakgosu@snu.ac.kr

## Usage Rights
- ✅ Personal use
- ✅ Educational use  
- ✅ Research use
- ❌ Commercial use (requires separate license)
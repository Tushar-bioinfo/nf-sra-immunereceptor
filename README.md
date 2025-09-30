# SRA to BAM: Immune Receptor Extraction Pipeline        

This pipeline automates the conversion of SRA accessions into sorted, indexed BAM files focused on immune receptor regions (IGH, IGK, IGL, TRA, TRB, TRG). It is built using **Nextflow DSL2** for modularity, supports SLURM-based HPC or local execution, and ensures all intermediate files are cleaned up automatically.

---

## Project Structure

```
sra_pipeline/
├── main.nf
├── nextflow.config
├── modules/
│   ├── prefetch.nf
│   ├── fasterq_dump.nf
│   ├── star_align.nf
│   ├── sort_bam.nf
│   ├── slice_receptors.nf
│   └── merge_index.nf
├── assets/
│   ├── manifest.csv              # Contains Run accessions (SRA IDs)
│   ├── key.ngc                   # dbGaP authorization key
│   └── ref/                      # STAR genome index directory
└── results/
    └── final_bams/              # Final output BAM + index
```

---

## Requirements

- **Nextflow** v24+
- SLURM or local setup
- Tools in PATH:
  - `sratoolkit` (prefetch, fasterq-dump)
  - `STAR`
  - `sambamba`
  - `samtools`

---

## Input Files

### 1. `assets/manifest.csv`

A CSV file with a header and a `Run` column containing SRA accession IDs:

```csv
Run
SRR12345678
SRR87654321
```

### 2. `assets/key.ngc`

Your NCBI dbGaP `.ngc` access key.

Download via NCBI dbGaP portal → place in `assets/key.ngc`.

### 3. `assets/ref/`

Path to STAR genome index (e.g., GRCh38) directory.

Generate via:

```bash
STAR --runThreadN 8 --runMode genomeGenerate      --genomeDir assets/ref      --genomeFastaFiles genome.fa      --sjdbGTFfile annotation.gtf      --sjdbOverhang 100
```

---

## Run the Pipeline

Environment Setup

This pipeline supports three options for managing dependencies and environments:

---

### Option 1: Conda (Default)
No container needed — just use your existing Conda environment.

#### Step 1: Create environment

```bash
conda env create -f environment.yml
conda activate immune-receptor-pipeline
```

Make sure all tools like `fasterq-dump`, `STAR`, `sambamba`, and `samtools` are available after activation.

#### Step 2: Run pipeline with Conda

```bash
nextflow run main.nf -profile local
```

> Conda will be used as the environment manager. All tools must be pre-installed or handled via your `environment.yml`.

---

### Option 2: Docker (for Local Use)

If you're running locally (not HPC), you can use Docker to encapsulate all dependencies.

#### Step 1: Build Docker image

```bash
docker build -t immune-pipeline .
```

#### Step 2: Run pipeline with Docker

```bash
nextflow run main.nf -profile docker
```

> This uses the `Dockerfile` and mounts the current working directory. No need to install tools manually.

---

### Option 3: Singularity (Recommended for HPC)

If you are running on HPC (e.g., with SLURM), Singularity ensures full compatibility.

#### Step 1: Enable Singularity in `nextflow.config`

This is already configured:

```groovy
profiles {
  slurm {
    process.executor = 'slurm'
    singularity.enabled = true
    singularity.autoMounts = true
  }
}
```

#### Step 2: Build or pull container

If using Docker on your local machine, convert to a Singularity image:

```bash
singularity build immune_pipeline.sif docker-daemon://immune-pipeline:latest
```

> Or use `--with-singularity <image>` if already available.

#### Step 3: Run on HPC

```bash
nextflow run main.nf -profile slurm
```

---
## Output

You’ll get sorted, indexed BAM files for each accession, focused on immune receptors:

```
results/final_bams/
├── SRR12345678.bam
├── SRR12345678.bam.bai
├── SRR87654321.bam
└── SRR87654321.bam.bai
```

These BAM files can be processed further using the **VDJ mining pipeline developed by the Blanck Lab**, which extracts **amino acid CDR3 sequences** and outputs results as `.csv` files. [Link](https://github.com/arpansahoo/vdj)

---

## Cleanup Logic

All intermediate files (FASTQ, temporary STAR output, unsorted BAMs, receptor-specific slices) are automatically deleted after final BAM generation. Each step passes only what's needed to the next, optimizing space and speed.

---

## Process Summary

| Step            | Tool            | Output                            |
|-----------------|-----------------|-----------------------------------|
| Prefetch        | `prefetch`      | `.sra` file                       |
| FASTQ Extraction| `fasterq-dump`  | paired FASTQ                      |
| Alignment       | `STAR`          | unsorted BAM                      |
| Sorting         | `sambamba`      | sorted + indexed BAM              |
| Region Slicing  | `samtools`      | receptor-specific BAMs            |
| Merging         | `sambamba merge`| final BAM + index                 |

---



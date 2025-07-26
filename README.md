# SRA to BAM: Immune Receptor Extraction Pipeline (Nextflow)

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

### SLURM (HPC)

```bash
nextflow run main.nf -profile slurm
```

### Local Machine

```bash
nextflow run main.nf -profile local
```

### Custom Parameters

```bash
nextflow run main.nf   --manifest assets/my_manifest.csv   --ngc_key assets/my_key.ngc   --ref_dir assets/custom_ref/
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



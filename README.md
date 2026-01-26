# VirusScan

**Version:** 1.2  
**Author:** Song Cao  
**Contact:** scao@wustl.edu  
**Original Release:** April 25, 2016  
**Last Updated:** 2026 (LSF chained job dependency support)

---

## Citation

If you use VirusScan, please cite:

> **Song Cao**, Michael C. Wendl, Matthew A. Wyczalkowski, Kristine Wylie, Kai Ye, Reyka Jayasinghe, Mingchao Xie, Song Wu, Beifang Niu, Robert Grubb III, Kimberly J. Johnson, Hiram Gay, Ken Chen, Janet S. Rader, John F. DiPersio, Feng Chen, and Li Ding.  
> **Divergent viral presentation among human tumors and adjacent normal tissues.**  
> *Scientific Reports*, 2016, **6:28294**.

---

## Overview

**VirusScan** is a fully automated and modular pipeline for detecting known viral sequences from NGS data.  
It is designed to run on **LSF (IBM Spectrum LSF / Platform LSF)** clusters using Docker.

All required tools are bundled inside the Docker image used by the pipeline.

---

## Key Features

- End-to-end viral detection from BAM files  
- RepeatMasker-based low-complexity filtering  
- Host genome BLAST filtering  
- Viral BLASTN detection and assignment  
- Uses NCBI taxonomy dump files directly  
- Per-sample chained LSF job dependencies  
- Supports full runs and partial restarts  

---

## What’s New

### Per-sample Job Dependency Chaining

Pipeline steps for each sample are submitted sequentially using LSF job dependencies:

```
-w "done(<previous_job>)"
```

Each sample runs independently, but steps within a sample never start before the previous step finishes.

---

## Requirements

### Runtime Environment

All software dependencies are included in the Docker image.

No local installation is required for:
- RepeatMasker
- BLAST / BLAST+
- Perl modules (BioPerl, etc.)
- VirusScan helper scripts

### Reference Data

Users must provide or configure paths to:
- Viral reference databases
- NCBI nt database (if used)
- NCBI taxonomy dump files (`nodes.dmp`, `names.dmp`, `merged.dmp`)
- Host reference genome (if applicable)

---

## Installation

```bash
git clone https://github.com/ding-lab/VirusScan.git
cd VirusScan
```

---

## Input Directory Structure

```
run_folder/
├── sample1/
│   └── sample1.bam
├── sample2/
│   └── sample2.bam
```

**Important:**  
The BAM file name must match the sample directory name.

---

## Usage

```bash
perl VirusScan.pl <run_folder> <step_number>
```

- `run_folder`: directory containing sample subdirectories  
- `step_number`: pipeline step to execute  

### Run the entire pipeline

```bash
perl VirusScan.pl <run_folder> 0
```

---

## Pipeline Steps

| Step | Description |
|------|------------|
| 0 | Run all steps (1–14) |
| 1 | Extract unmapped non-human reads |
| 2 | Split files for RepeatMasker |
| 3 | Submit RepeatMasker job array |
| 4 | Sequence quality control |
| 5 | Split files for human genome BLAST |
| 6 | Submit human genome BLAST |
| 7 | Parse human genome BLAST results |
| 8 | Pool and split files for BLASTN |
| 9 | Submit BLASTN job array |
| 10 | Parse BLASTN results |
| 11 | Generate BLASTN summary |
| 12 | Assignment report per sample |
| 13 | Assignment summary per sample |
| 14 | Generate final run-level report |

---

## Chained Execution Shortcuts

| Step | Runs |
|------|------|
| 22 | Steps 2–14 |
| 23 | Steps 3–14 |
| 24 | Steps 4–14 |
| 25 | Steps 5–14 |
| 26 | Steps 6–14 |
| 27 | Steps 7–14 |
| 28 | Steps 8–14 |
| 29 | Steps 9–14 |
| 30 | Steps 10–14 |
| 31 | Steps 11–14 |
| 32 | Steps 12–14 |
| 33 | Steps 13–14 |

---

## Running on LSF with Docker

VirusScan submits all jobs using Docker-enabled LSF queues.

Example submission pattern:

```bash
bsub -q <queue> -n <threads> \
  -R "select[mem>30000] rusage[mem=30000]" -M 30000000 \
  -a 'docker(<docker_image>)' \
  -o <out.log> -e <err.log> \
  bash <job_script.sh>
```

---

## Contact

**Song Cao**  
Email: scao@wustl.edu

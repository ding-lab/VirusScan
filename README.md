# VirusScan

**Version:** 1.2  
**Author:** Song Cao  
**Contact:** scao@wustl.edu  
**Original Release:** April 25, 2016  
**Last Updated:** November 2026 (`--sre` re-run option, run report dependency, BLASTN parsing fix)

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
- Run report waits for all samples to finish  
- Supports full runs, partial restarts, and re-parsing of BLASTN results  

---

## What’s New

### Per-sample Job Dependency Chaining

Pipeline steps for each sample are submitted sequentially using LSF job dependencies:

```
-w "done(<previous_job>)"
```

Each sample runs independently, but steps within a sample never start before the previous step finishes.

### Run Report Waits for All Samples

In chained runs (step 0 and steps 22–33), the run-level report (step 14) starts only after the last job (step 13) of **every** sample is done:

```
-w "done(<sample1_job>) && done(<sample2_job>) && ..."
```

- If any sample job fails, the report job stays pending. Kill the pending jobs and restart from the failed step.
- If any job cannot be submitted, the report job is not submitted and the affected samples are listed.

### Re-run Option for Step 10 (`--sre`)

By default, step 10 skips files that already have a complete parsed result. Use `--sre 1` to parse them again.

### BLASTN Parsing Fix

- Merged nt records (one sequence shared by several GIs) are now assigned a taxonomy ID. Previously the lookup failed for these records, and viral reads were reported as `Ambiguous`.
- A hit without a taxonomy ID is no longer counted as a hit to another species.

Results parsed before this fix should be re-parsed with `--step 30 --sre 1`.

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

Run the script from the VirusScan directory:

```bash
cd VirusScan
perl VirusScan.pl --rdir <run_folder> --log <log_folder> --q <queue> \
  --users <compute_username> --groupname <job_group> --step <step_number> [--sre 0|1]
```

| Option | Description |
|--------|-------------|
| `--rdir` | Full path of the directory containing sample subdirectories |
| `--log` | Full path of the directory for job scripts, LSF logs, and the taxonomy index (`taxo`) |
| `--q` | LSF queue name |
| `--users` | Compute username, used in the LSF job group `/<users>/<groupname>` |
| `--groupname` | LSF job group name |
| `--step` | Pipeline step to execute (see below) |
| `--sre` | Re-run status for step 10: `0` skips complete parsed results (default), `1` removes existing parsed results and parses again |

All options except `--sre` are required.

### Run the entire pipeline

```bash
perl VirusScan.pl --rdir /path/to/run_folder --log /path/to/run_folder.log --q general \
  --users <compute_username> --groupname <job_group> --step 0
```

### Re-parse BLASTN results and update all reports

```bash
perl VirusScan.pl --rdir /path/to/run_folder --log /path/to/run_folder.log --q general \
  --users <compute_username> --groupname <job_group> --step 30 --sre 1
```

With `--sre 1`, the old `.blastn.parsed` and `.BNfiltered.fa` files are deleted when each job starts.

### Continue after some step 10 jobs failed

```bash
perl VirusScan.pl --rdir /path/to/run_folder --log /path/to/run_folder.log --q general \
  --users <compute_username> --groupname <job_group> --step 30
```

Only missing or incomplete parsed results are redone.

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
| 34 | Step 14 only |

Steps 14 and 34 do not wait for other jobs. Run them only after all samples are finished.
The final report is written to `<run_folder>/Analysis_Report_<run_folder_name>.tsv`.

---

## Taxonomy Index

Step 10 uses a taxonomy index in `<log_folder>/taxo`. It is built from `nodes.dmp` and `names.dmp` the first time step 10 runs, and reused afterwards.

- Do not delete or move `taxo` before re-running step 10. Many jobs start at the same time, and building the index from several jobs can corrupt it.
- A corrupt index gives wrong or empty names in the lineage, e.g. an unrelated species inside a viral lineage, or `;;`.
- To check the index, confirm that the number of lines in `id2names` is the largest taxonomy ID in `names.dmp` plus 1, and that it has no NUL bytes:

```bash
wc -l < <log_folder>/taxo/id2names
cut -f1 names.dmp | sort -n | tail -1
tr -cd '\000' < <log_folder>/taxo/id2names | wc -c
```

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

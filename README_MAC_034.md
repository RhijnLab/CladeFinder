# CladeFinder v0.3.4 — macOS Docker-isolated runner

This release addresses macOS Docker Desktop bind-mount `ENOENT` errors by **not mounting any host directory into Docker**. The `cladefinder-mac` executable uses `docker cp` to transfer a copy of the reference database into a stopped/started container, runs CladeFinder against Docker-local files, checks outputs, and copies them back to the host. The original `database.json` and other source files are not overwritten by preparation.

## Requirements

- macOS, Docker Desktop running, Python 3 (macOS host), enough Docker Desktop disk space for the database and FASTQ files.
- Existing `my-clade-data/clade_db/database.json`, `reference.nex`, and reference database from v0.3.
- This package does **not** include your 4.5 GB VCF or reference database.

## Build

```sh
cd /path/to/extracted/cladefinder
dot_clean .
find . -name '._*' -type f -delete
docker build -t cladefinder:0.3.4 .
docker run --rm cladefinder:0.3.4 --help
```

## Prepare

From any directory:

```sh
/path/to/cladefinder/cladefinder-mac prepare /path/to/my-clade-data 4
```

The runner copies `clade_db` to Docker, regenerates the 513-reference SNP alignment and Nexus topology in Docker, optimises branch lengths with RAxML-NG, then publishes verified files to the host. It publishes `placement_manifest.json` last. Back up `clade_db` first for extra protection. No need to rebuild the VCF database.

## Validate EPA-ng (5 tip-withheld isolates)

```sh
/path/to/cladefinder/cladefinder-mac validate-placement /path/to/my-clade-data 5 4
```

Output: `my-clade-data/placement-validation/`.

## Classify paired-end FASTQ

```sh
/path/to/cladefinder/cladefinder-mac classify /path/to/my-clade-data \
  --r1 /path/to/sample_R1.fastq.gz --r2 /path/to/sample_R2.fastq.gz \
  --sample sample --threads 4 --phylogeny auto
```

For single-end reads omit `--r2`. Output: `my-clade-data/results/sample/`.

## Scientific and operational limitations

- EPA-ng uses `GTR+G` while RAxML-NG fits `GTR+G+ASC_LEWIS`; the likelihood models are not fully consistent. Phylogenetic support values must be considered exploratory, **not calibrated confidence probabilities**.
- Internal leave-one-out on a SNP panel chosen using all isolates is not an independent evaluation. Do not use for clinical decisions without external validation.
- SNP distance and phylogenetic placement are separate outputs; disagreement must be reviewed.
- `docker cp` transfers files, which requires free disk space and can take time for FASTQ inputs.
- Host-side publication uses atomic per-file replacement, but not a transaction spanning all files. Manifest is published last.
- This workflow has not been end-to-end verified on the user's Mac. Do not claim otherwise.

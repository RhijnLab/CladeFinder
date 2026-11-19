# CladeFinder: phylogenetic and SNP-distance cluster classification

Research software for assigning *Aspergillus fumigatus* isolates to one of five pre-existing clusters. Supports FASTQ input, a reference VCF, metadata, a Nexus reference topology, SNP nearest-neighbour classification and experimental EPA-ng phylogenetic placement.

## Important status

**This is meant for research purposes only, and not as a validated clinical diagnostic.** Full 513-isolate internal tip-withheld validation achieved 485/513 correct (94.5%), 512/513 assigned, macro one-vs-rest AUROC 0.960. The SNP panel was selected from all 513 reference strains before cross-validation

## Installation (macOS)

Install Docker Desktop, then in this repository:

```sh
dot_clean .
find . -name '._*' -type f -delete
docker build -t cladefinder:0.3.7 .
```

Put your `clade_db` directory in a data directory and use the included `cladefinder-mac` wrapper. The wrapper copies files into Docker rather than bind-mounting the macOS filesystem.

```sh
./cladefinder-mac prepare /absolute/path/to/my-clade-data 4
./cladefinder-mac validate-placement /absolute/path/to/my-clade-data 5 4
./cladefinder-mac validate-placement /absolute/path/to/my-clade-data 513 4
python3 roc_report.py /absolute/path/to/my-clade-data/placement-validation/placement_leave_one_out.tsv --output roc-results
```

Paired-end FASTQ:

```sh
./cladefinder-mac classify /absolute/path/to/my-clade-data --r1 /absolute/path/to/sample_R1.fastq.gz --r2 /absolute/path/to/sample_R2.fastq.gz --sample sample --threads 4 --phylogeny auto
```

Single-end: omit `--r2`.

The wrapper defaults to the Docker image tag `cladefinder:0.3.7`. `cladefinder --help` lists database-building commands; the large reference VCF stays on your machine.

## Repository files

- `src/cladefinder/`: classification and placement pipeline
- `cladefinder-mac`: isolated Docker workflow for macOS
- `roc_report.py`: five one-vs-rest ROC curves from placement validation
- `tests/`: unit tests
- `environment.yml`, `Dockerfile`, `pyproject.toml`: reproducible installation
- `docs/validation/`: supplied internal validation summary and figures

## License

MIT. See `LICENSE`.

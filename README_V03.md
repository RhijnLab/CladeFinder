# CladeFinder v0.3 — two independent classification methods

**Status: research prototype, not validated for clinical or diagnostic use.**

The existing `build`, `validate` and `classify` commands remain. `classify`
now outputs `classification.json` with (1) `predicted_cluster` from SNP
nearest-neighbour distances and (2) `phylogenetic_placement.predicted_cluster`
from EPA-ng placement onto the supplied Nexus **topology**. The methods can
**disagree**; neither is silently substituted for the other.

## Docker (macOS)

```bash
# Inside the unpacked cladefinder source directory:
dot_clean .
find . -name '._*' -type f -delete
docker build -t cladefinder:0.3 .
docker run --rm cladefinder:0.3 --help
```

On Apple Silicon Docker Desktop may use linux/arm64; if a dependency fails to
resolve, try `docker build --platform linux/amd64 -t cladefinder:0.3 .`
(the emulated build is slower). Docker Desktop should have sufficient RAM and
free disk space; RAxML-NG on 513 × 20,000 sites can be resource-intensive.

### Reuse existing v0.2 database

In your `my-clade-data` directory:

```bash
docker run --rm -v "$PWD:/data" cladefinder:0.3 tree-report \
  --database /data/clade_db --output /data/tree_report

docker run --rm -v "$PWD:/data" cladefinder:0.3 prepare-placement \
  --database /data/clade_db --threads 4
```

The second command reconstructs `reference_snps.fasta` from the saved v0.2
reference SNP strings, strips BEAST time-scale branch lengths from the Nexus
tree, and fits branch lengths on that topology using RAxML-NG. It then saves
`placement_tree.nwk`. It **does not** need the 4.5 GB VCF again.

Classify paired-end FASTQs:

```bash
docker run --rm -v "$PWD:/data" cladefinder:0.3 classify \
  --r1 /data/fastq/sample_R1.fastq.gz \
  --r2 /data/fastq/sample_R2.fastq.gz \
  --database /data/clade_db --output /data/results/sample \
  --phylogeny required --threads 4
```

For single-end FASTQ omit `--r2`. `--phylogeny auto` (default) returns an
explicit unavailable status if placement cannot run; `required` fails instead.

### Validate both methods

SNP-distance leave-one-out (fast):

```bash
docker run --rm -v "$PWD:/data" cladefinder:0.3 validate \
  --database /data/clade_db --output /data/validation-v03
```

True withheld-tip EPA-ng validation (slow: tree fitting per fold):

```bash
docker run --rm -v "$PWD:/data" cladefinder:0.3 validate-placement \
  --database /data/clade_db --output /data/placement-validation \
  --limit 20 --threads 4
```

`--limit 20` is an initial smoke test, not a statistically representative
estimate. Full 513-fold validation requires `--limit 513` and substantial CPU.
This procedure withholds each tip from the reference tree and alignment before
placing it, avoiding the invalid practice of placing a tip on a tree that still
contains itself. The selected SNP sites themselves were ascertained using all
isolates; this is **not** independent validation. For a proper held-out
benchmark, use separate isolates and FASTQs not used in site selection.

## Output files

- `tree_report/tree_clusters.tsv`: cluster monophyly, MRCA contamination,
  number of separate pure subtrees.
- `clade_db/reference_snps.fasta`: reference SNP alignment.
- `clade_db/placement_tree.nwk`: fixed Nexus topology with fitted branches.
- `results/sample/classification.json`: **both** SNP and tree predictions.
- `results/sample/placement/epa_result.jplace`: raw EPA-ng placements.
- `placement-validation/placement_leave_one_out.tsv`: withheld-tip results.

For the original installation and VCF build instructions, see `README.md`.

# CladeFinder v0.3.3 — macOS Docker Desktop file-handling fix

This patch replaces `prepare-placement` with a scratch-first implementation. It generates the alignment and Nexus topology, runs RAxML-NG, and verifies the tree and model **inside the container's own temporary filesystem**, then copies completed artifacts to the bind-mounted database. This avoids the observed intermittent `ENOENT` on newly created temporary files inside the macOS Docker bind mount.

## Build

From the extracted `cladefinder` source directory:

```bash
dot_clean .
find . -name '._*' -type f -delete
docker build -t cladefinder:0.3.3 .
```

## Prepare

From the folder containing your existing `clade_db` directory:

```bash
cp -R clade_db clade_db_backup_before_033
docker run --rm -v "$PWD:/data" cladefinder:0.3.3 prepare-placement --database /data/clade_db --threads 4
```

## Verify

```bash
ls -lh clade_db/{database.json,reference_snps.fasta,nexus_topology.nwk,placement_tree.nwk,placement_model.txt,placement_manifest.json}
docker run --rm -v "$PWD:/data" cladefinder:0.3.3 tree-report --database /data/clade_db --output /data/tree_report
```

## Validate first (experimental)

```bash
docker run --rm -v "$PWD:/data" cladefinder:0.3.3 validate-placement --database /data/clade_db --output /data/placement-validation --limit 5 --threads 4
```

**Scientific limitation:** RAxML fits the SNP-only alignment using `GTR+G+ASC_LEWIS`, but EPA-ng is invoked with `GTR+G`; its likelihood weights are exploratory and are **not** calibrated probabilities. This release fixes the preparation file-handling bug; it does not establish phylogenetic classification accuracy. Do not use the output for clinical decisions.

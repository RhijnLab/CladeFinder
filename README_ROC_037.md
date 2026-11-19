# CladeFinder v0.3.7: ROC-enabled validation

This release adds five cluster score columns (`score_cluster_1` ... `score_cluster_5`) and `unresolved_weight` to `placement_leave_one_out.tsv`. Scores are summed EPA-ng placement likelihood weight ratios over cluster-resolving edges, **not calibrated probabilities**. Zero means no weight was mapped to that cluster, not proof of exclusion. It also creates the EPA-ng output directory before invoking EPA-ng.

## macOS commands

In the extracted source directory:

```bash
dot_clean .
find . -name '._*' -type f -delete
docker build -t cladefinder:0.3.7 .
```

In your data directory (`/Volumes/LaCie/my-clade-data`):

```bash
/path/to/extracted/cladefinder/cladefinder-mac validate-placement "$PWD" 5 4
head -n 2 placement-validation/placement_leave_one_out.tsv
```

Confirm `score_cluster_1` ... `score_cluster_5` appear in the header. Then run all 513 folds:

```bash
/path/to/extracted/cladefinder/cladefinder-mac validate-placement "$PWD" 513 4
python3 /path/to/extracted/cladefinder/roc_report.py placement-validation/placement_leave_one_out.tsv --output roc-results
```

The ROC tool uses only Python's standard library and produces `roc-results/roc_curves.svg` and `roc-results/roc_summary.json`. **Do not interpret the initial five-fold validation as ROC evidence:** all five examples are Cluster 3, so one-vs-rest AUROC cannot be calculated for every cluster.

The full run may take hours, and this runner does not provide robust resume/checkpoint functionality. The SNP sites were selected using all reference isolates, so the leave-one-out evaluation is not completely independent. RAxML uses ASC_LEWIS while EPA-ng uses GTR+G; these models differ. Independent FASTQ validation is still needed.

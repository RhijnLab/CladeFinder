# CladeFinder v0.3.5 — leave-one-out invariant-site fix

This release fixes the RAxML-NG error `ASC_LEWIS ... contains N invariant sites` during `validate-placement`.

For each withheld isolate, the validator now selects variable A/C/G/T sites **using only the remaining reference isolates**, and applies the same indices to the query. The output TSV includes `original_sites`, `retained_sites`, and `removed_invariant_or_ambiguous_sites` for every fold. This avoids leaking the withheld query into site selection.

## macOS installation

Extract the ZIP, open Terminal in the extracted `cladefinder` directory, then run:

```sh
dot_clean .
find . -name '._*' -type f -delete
docker build -t cladefinder:0.3.5 .
```

Change to the directory containing your existing `clade_db`, e.g. `/Volumes/LaCie/my-clade-data`, and run the included macOS wrapper using its actual path:

```sh
/path/to/cladefinder/cladefinder-mac validate-placement "$PWD" 5 4
```

The script defaults to image `cladefinder:0.3.5`. It creates a temporary container, copies the existing database into it, runs the analysis without a host bind mount, and copies results to `placement-validation/` only if validation completes successfully. It does not rebuild the VCF database.

The existing `prepare` command remains available:

```sh
/path/to/cladefinder/cladefinder-mac prepare "$PWD" 4
```

## Checks

```sh
cat placement-validation/placement_validation.json
head -n 8 placement-validation/placement_leave_one_out.tsv
```

## Scope and scientific limitations

- RAxML-NG uses `GTR+G+ASC_LEWIS`, while EPA-ng currently uses `GTR+G`. Their likelihood models are not ascertainment-equivalent; placement likelihood weights are exploratory, not calibrated cluster probabilities.
- Site selection for the original 20,000-SNP database used the full dataset. Fold-specific filtering removes newly invariant sites but does not undo that original panel ascertainment. This is internal leave-one-out testing, not an independent FASTQ validation.
- EPA-ng may expose other runtime or model issues. A passing Python test suite does not establish that the complete Docker placement workflow has passed.
- Clade 1 is non-monophyletic in the supplied Nexus topology; edge-based cluster assignments should be interpreted cautiously.

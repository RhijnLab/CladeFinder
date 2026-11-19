#!/bin/sh
set -eu
# Run from the extracted cladefinder repository after creating an empty GitHub repository.
: "${GITHUB_REPO_URL:?Set GITHUB_REPO_URL to your GitHub repository HTTPS or SSH URL}"
git init
git add .
git commit -m 'Initial CladeFinder research release with internal validation'
git branch -M main
git remote add origin "$GITHUB_REPO_URL"
git push -u origin main

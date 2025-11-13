#!/bin/bash
# Simple script to sync the published docs to the docs/ folder for local viewing

echo "Syncing results/publish_docs/ to docs/..."
rsync -av --delete results/publish_docs/ docs/

echo ""
echo "✓ Docs updated!"
echo ""
echo "Summary plots location:"
echo "  - docs/htmls/phenotypes_overlaid.html"
echo "  - docs/index.html (see 'Integrated summary plots' section)"

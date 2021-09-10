#!/usr/bin/env bash

TEMP=tests/golden/test_output.ann
REFERENCE=tests/golden/output_GCF_000280675.1_ASM28067v1_genomic.ann


python src/main.py \
  --gff3 tests/golden/GCF_000280675.1_ASM28067v1_genomic.gff.gz \
  --fasta tests/golden/GCF_000280675.1_ASM28067v1_genomic.fna.gz \
  --transl_table=4 \
  --metadata ./tests/golden/metadata.toml \
  --output "$TEMP" \
  2> /dev/null

if [[ "$(diff "$TEMP" "$REFERENCE")" ]]; then
    echo "[FAIL] Failed in golden test"
else
    echo "[SUCCESS]"
fi

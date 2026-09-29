#!/usr/bin/env bash

set -euo pipefail

max_events="${1:--1}"
sample_file="config/phase4_samples.txt"

if ! command -v root >/dev/null 2>&1; then
  echo "ERROR: ROOT is not available in PATH." >&2
  exit 1
fi

while IFS='|' read -r label record url; do
  if [[ -z "${label}" || "${label}" == \#* ]]; then
    continue
  fi

  echo "Computing normalization for ${label} (record ${record})"
  root -l -b -q \
    "macros/computeNormalization.C(\"${url}\",\"${label}\",${max_events})"
done < "${sample_file}"


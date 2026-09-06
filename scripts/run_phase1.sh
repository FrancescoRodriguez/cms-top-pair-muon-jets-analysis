#!/usr/bin/env bash

set -euo pipefail

max_events="${1:-1000}"
sample_file="config/phase1_samples.txt"

if ! command -v root >/dev/null 2>&1; then
  echo "ERROR: ROOT is not available in PATH." >&2
  exit 1
fi

while IFS='|' read -r label sample_type record url; do
  if [[ -z "${label}" || "${label}" == \#* ]]; then
    continue
  fi

  echo "Inspecting ${label} (${sample_type}, record ${record})"
  root -l -b -q \
    "macros/inspectNanoAOD.C(\"${url}\",\"${label}\",${max_events})"
done < "${sample_file}"


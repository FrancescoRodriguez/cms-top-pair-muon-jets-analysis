#!/usr/bin/env bash

set -euo pipefail

max_events="${1:--1}"

if [[ ! -s config/filelists/data_Run2016G.txt ]]; then
  echo "Input file lists not found. Run ./scripts/fetch_file_lists.sh first." >&2
  exit 1
fi

root -l -b -q "macros/measureInclusive.C(${max_events})"


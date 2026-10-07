#!/usr/bin/env bash

set -u

checkpoint_dir="output/checkpoints"
mc_file_limit="${1:--1}"
mkdir -p "$checkpoint_dir"
failed=0

process_list() {
  local label="$1"
  local kind="$2"
  local list="$3"
  local limit="${4:--1}"
  local index=0

  while IFS= read -r url; do
    [[ -z "$url" || "${url:0:1}" == "#" ]] && continue
    index=$((index + 1))
    if [[ "$limit" -gt 0 && "$index" -gt "$limit" ]]; then
      break
    fi
    local number
    number=$(printf '%04d' "$index")
    local output="${checkpoint_dir}/${label}_${number}.txt"
    if [[ -s "$output" ]] && grep -q '^status=ok$' "$output"; then
      echo "[skip] ${label} file ${index}: already completed"
      continue
    fi

    local success=0
    for attempt in 1 2 3; do
      echo "[run] ${label} file ${index}, attempt ${attempt}"
      root -l -b -q \
        "macros/processOneFile.C(\"${url}\",\"${kind}\",\"${output}\")"
      if [[ -s "$output" ]] && grep -q '^status=ok$' "$output"; then
        success=1
        break
      fi
    done
    if [[ "$success" -eq 0 ]]; then
      echo "[failed] ${label} file ${index}; continue with the other files" >&2
      failed=1
    fi
  done < "$list"
}

process_list dataG data config/filelists/data_Run2016G.txt -1
process_list dataH data config/filelists/data_Run2016H.txt -1
process_list tt tt config/filelists/ttbar_semileptonic.txt "$mc_file_limit"
process_list w w config/filelists/wjets.txt "$mc_file_limit"
process_list dy dy config/filelists/dyjets.txt "$mc_file_limit"

if [[ "$failed" -ne 0 ]]; then
  echo "Some files failed. Run this script again: completed files will be skipped." >&2
  exit 1
fi

if [[ "$mc_file_limit" -gt 0 ]]; then
  root -l -b -q \
    "macros/combineCheckpoints.C(\"${checkpoint_dir}\",${mc_file_limit},${mc_file_limit},${mc_file_limit})"
else
  root -l -b -q 'macros/combineCheckpoints.C()'
fi

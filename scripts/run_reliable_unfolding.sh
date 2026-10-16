#!/usr/bin/env bash

set -u

mc_files="${1:-10}"
checkpoint_dir="output/unfolding_checkpoints"
mkdir -p "$checkpoint_dir"
failed=0

process_list() {
  local label="$1"
  local kind="$2"
  local list="$3"
  local limit="$4"
  local index=0

  while IFS= read -r url; do
    [[ -z "$url" || "${url:0:1}" == "#" ]] && continue
    index=$((index + 1))
    if [[ "$limit" -gt 0 && "$index" -gt "$limit" ]]; then
      break
    fi
    local number
    number=$(printf '%04d' "$index")
    local output="${checkpoint_dir}/${label}_${number}.root"
    if [[ -s "$output" ]]; then
      echo "[skip] ${label} file ${index}: checkpoint present"
      continue
    fi

    local success=0
    for attempt in 1 2 3; do
      echo "[run] ${label} file ${index}, attempt ${attempt}"
      root -l -b -q \
        "macros/processUnfoldingFile.C(\"${url}\",\"${kind}\",\"${output}\")"
      if [[ -s "$output" ]]; then
        success=1
        break
      fi
    done
    if [[ "$success" -eq 0 ]]; then
      echo "[failed] ${label} file ${index}" >&2
      failed=1
    fi
  done < "$list"
}

process_list dataG data config/filelists/data_Run2016G.txt -1
process_list dataH data config/filelists/data_Run2016H.txt -1
process_list tt tt config/filelists/ttbar_semileptonic.txt "$mc_files"
process_list w w config/filelists/wjets.txt "$mc_files"
process_list dy dy config/filelists/dyjets.txt "$mc_files"

if [[ "$failed" -ne 0 ]]; then
  echo "Some files failed. Run the same command again; completed files are skipped." >&2
  exit 1
fi

root -l -b -q \
  "macros/combineUnfolding.C(\"${checkpoint_dir}\",${mc_files},${mc_files},${mc_files},4)"

#!/usr/bin/env bash

set -euo pipefail

output_directory="config/filelists"
mkdir -p "${output_directory}"

fetch_record_files() {
  local label="$1"
  local record="$2"
  local temporary_file
  temporary_file="$(mktemp)"

  curl -sS -L "https://opendata.cern.ch/api/records/${record}" \
    | jq -r '.metadata._file_indices[]?.files[]?.uri,
             .metadata.files[]?.uri' \
    | awk 'NF && !seen[$0]++' > "${temporary_file}"

  if [[ ! -s "${temporary_file}" ]]; then
    echo "ERROR: no files found for CERN record ${record}." >&2
    exit 1
  fi
  mv "${temporary_file}" "${output_directory}/${label}.txt"
  echo "${label}: $(wc -l < "${output_directory}/${label}.txt") files"
}

fetch_record_files "data_Run2016G" 30530
fetch_record_files "data_Run2016H" 30563
fetch_record_files "ttbar_semileptonic" 67993
fetch_record_files "wjets" 69747
fetch_record_files "dyjets" 35671

curl -sS -L \
  "https://opendata.cern.ch/record/14220/files/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt" \
  -o "${output_directory}/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"

echo "Downloaded the CMS 2016 certified-luminosity JSON."


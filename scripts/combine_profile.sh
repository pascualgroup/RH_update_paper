#!/bin/bash

set -euo pipefail
shopt -s nullglob

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 RUN_ID" >&2
    exit 1
fi

RUN_ID=$1

combine_files() {
    local -n files_ref=$1
    local output=$2
    local tool=$3
    
    if [[ ${#files_ref[@]} -eq 0 ]]; then
        echo "No files matched for ${output}" >&2
        exit 1
    fi
    
    "$tool" "${files_ref[@]}" > "${output}"
}

result_fit_files=( *"results_${RUN_ID}"*fit.csv )
result_refined_files=( *"results_${RUN_ID}"*refined.csv )
result_refined_second_files=( *"results_${RUN_ID}"*refined_second.csv )

manifest_fit_files=( *"manifest_${RUN_ID}"*fit.json )
manifest_refined_files=( *"manifest_${RUN_ID}"*refined.json )
manifest_refined_second_files=( *"manifest_${RUN_ID}"*refined_second.json )

combine_files result_fit_files "combined_results_${RUN_ID}_fit.csv" csvstack
combine_files result_refined_files "combined_results_${RUN_ID}_refined.csv" csvstack
combine_files result_refined_second_files "combined_results_${RUN_ID}_refined_second.csv" csvstack

combine_files manifest_fit_files "combined_manifest_${RUN_ID}_fit.json" cat
combine_files manifest_refined_files "combined_manifest_${RUN_ID}_refined.json" cat
combine_files manifest_refined_second_files "combined_manifest_${RUN_ID}_refined_second.json" cat

rm -- "${result_fit_files[@]}" \
"${result_refined_files[@]}" \
"${result_refined_second_files[@]}" \
"${manifest_fit_files[@]}" \
"${manifest_refined_files[@]}" \
"${manifest_refined_second_files[@]}"

echo "Combined and removed source files for RUN_ID=${RUN_ID}"

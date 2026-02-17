#!/usr/bin/env bash
# Submit an array run (scripts/run_fitting.sh) then submit the aggregator
# (scripts/combine_array_outputs.sh) with a dependency so it runs after the
# array completes. Both jobs receive the same RUN_ID.
#
# Usage:
#   scripts/submit_run_and_combine.sh [--run-id RUNID] [--array 1-100] [--no-export]
#
# If --run-id is not provided a timestamp + short random suffix will be created.
set -euo pipefail

BASEDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_SCRIPT="$BASEDIR/run_fitting.sh"
COMBINE_SCRIPT="$BASEDIR/combine_array_outputs.sh"

print_usage() {
  cat <<'USAGE'
Usage: submit_run_and_combine.sh [--run-id RUNID] [--array RANGE]

Options:
  --run-id RUNID   Provide a RUN_ID to use for both jobs. If omitted a value
                    like 20251002T120000Z_ab12 will be generated.
  --array RANGE    Optional array range to pass to sbatch (e.g. 1-100).
  --help           Show this help and exit.

Example:
  ./scripts/submit_run_and_combine.sh --run-id 20251002T120000Z_42 --array 1-50
USAGE
}

RUN_ID=""
ARRAY_RANGE=""

while [[ ${#} -gt 0 ]]; do
    case "$1" in
        --run-id)
        RUN_ID="$2"; shift 2;;
        --array)
        ARRAY_RANGE="$2"; shift 2;;
        --help|-h)
        print_usage; exit 0;;
        *)
        echo "Unknown arg: $1" >&2; print_usage; exit 2;;
    esac
done

if [ -z "$RUN_ID" ]; then
    # generate timestamp + short random hex (openssl if available, else date+pid)
    TS=$(date -u +%Y%m%dT%H%M%SZ)
    if command -v openssl >/dev/null 2>&1; then
        SUF=$(openssl rand -hex 2)
    else
        SUF=$((RANDOM % 10000))
    fi
    RUN_ID="${TS}_${SUF}"
fi

echo "Submitting run_fitting array with RUN_ID=${RUN_ID}"
echo "Submitted in $(date --rfc-3339=seconds)"
# Build sbatch args incrementally; avoid commas inside array elements
SBATCH_COMMON=(--export "ALL,RUN_ID=${RUN_ID}")
if [ -n "$ARRAY_RANGE" ]; then
    SBATCH_COMMON+=(--array "${ARRAY_RANGE}")
fi

# Submit the array job; capture the job id (parsable)
ARRAY_JOBID=$(sbatch --parsable "${SBATCH_COMMON[@]}" "$RUN_SCRIPT")
if [ -z "$ARRAY_JOBID" ]; then
    echo "Failed to submit array job" >&2
    exit 3
fi
echo "Submitted array job id: ${ARRAY_JOBID}"

echo "Submitting combine job dependent on array job completion (afterany:${ARRAY_JOBID})"
COMBINE_JOBID=$(sbatch --parsable --dependency=afterany:${ARRAY_JOBID} --export=ALL,RUN_ID=${RUN_ID} "$COMBINE_SCRIPT")
if [ -z "$COMBINE_JOBID" ]; then
    echo "Failed to submit combine job" >&2
    exit 4
fi

cat <<-MSG
Submitted:
  Array job id:    ${ARRAY_JOBID}
  Combine job id:  ${COMBINE_JOBID}
  RUN_ID:          ${RUN_ID}

To monitor:
  squeue -j ${ARRAY_JOBID}
  squeue -j ${COMBINE_JOBID}
MSG

exit 0

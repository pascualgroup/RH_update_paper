#!/bin/bash
#SBATCH --job-name=VC_combine
#SBATCH --output=out/combine_%A.out
#SBATCH --error=error/combine_%A.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4GB
#SBATCH --cpus-per-task=1
#SBATCH --account=torch_pr_134_general

set -euo pipefail

########################################
# UPDATED: shared scratch setup + sync with $HOME
########################################
SCRATCH_BASE="${SCRATCH_BASE:-/scratch}"
PROJECT_NAME="${PROJECT_NAME:-$(basename "$SLURM_SUBMIT_DIR")}"
SCRATCH_DIR="${SCRATCH_BASE}/${USER}/${PROJECT_NAME}"
SOURCE_DIR="${SOURCE_DIR:-$HOME/${PROJECT_NAME}}"

mkdir -p "$SCRATCH_DIR"

if [ ! -d "$SOURCE_DIR" ]; then
    echo "[$(date)] ERROR: SOURCE_DIR not found: $SOURCE_DIR" >&2
    exit 1
fi

echo "[$(date)] Using SCRATCH_DIR=$SCRATCH_DIR"
echo "[$(date)] Syncing repo (incl. results) from $SOURCE_DIR to scratch"

# Here we *do* copy results/manifests, since the combiner needs them
rsync -a --delete \
--exclude='.git' \
--exclude='out' \
--exclude='error' \
"$SOURCE_DIR/" "$SCRATCH_DIR/"

cd "$SCRATCH_DIR"
mkdir -p results/aggregated

########################################
# EXISTING: module + RUN_ID parsing
########################################
module unload r || true
module load r/4.5.1 || true

RUN_ID_ENV="${RUN_ID:-}"
ARG_RUN_ID="${1:-}"
ARRAY_ARG="${2:-}"

if [ -n "$ARG_RUN_ID" ]; then
    RUN_ID="$ARG_RUN_ID"
    elif [ -n "$RUN_ID_ENV" ]; then
    RUN_ID="$RUN_ID_ENV"
else
    echo "ERROR: RUN_ID must be provided via environment or first argument" >&2
    exit 1
fi

# optional array id may be passed as second arg (e.g. A5 or 5)
if [ -n "$ARRAY_ARG" ]; then
    if [[ "$ARRAY_ARG" =~ ^A[0-9]+$ ]]; then
        ARRAY_ID="$ARRAY_ARG"
    else
        ARRAY_ID="A${ARRAY_ARG}"
    fi
else
    ARRAY_ID=""
fi

echo "Combining results for RUN_ID=$RUN_ID ARRAY=$ARRAY_ID"

########################################
# EXISTING: R invocation (unchanged)
########################################
RSCRIPT="scripts/combine_array_outputs.R"
RARGS=("--run_id=${RUN_ID}")
if [ -n "$ARRAY_ID" ]; then
    RARGS+=("--array=${ARRAY_ID}")
fi
RARGS+=("--results-dir=results" "--manifests-dir=manifests" "--out-dir=results/aggregated")

echo "Running: Rscript $RSCRIPT ${RARGS[*]}"
Rscript "$RSCRIPT" "${RARGS[@]}"

########################################
# UPDATED: sync aggregated results back to $HOME repo
########################################
echo "[$(date)] Syncing aggregated results back to $SOURCE_DIR"

rsync -a results/aggregated/ "$SOURCE_DIR/results/aggregated/"

########################################
# UPDATED: delete per-array outputs for this RUN_ID
########################################
echo "[$(date)] Cleaning per-array outputs for RUN_ID=$RUN_ID in $SOURCE_DIR"

# Remove per-array result files (but keep aggregated/)
# adjust patterns as needed to match your naming scheme
find "$SOURCE_DIR/results" -maxdepth 1 -type f \
-name "*${RUN_ID}_*" -print -delete

# Remove manifests for this RUN_ID
if [ -d "$SOURCE_DIR/manifests" ]; then
    find "$SOURCE_DIR/manifests" -maxdepth 1 -type f \
    -name "*${RUN_ID}_*" -print -delete
fi

echo "Combine job finished"

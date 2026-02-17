#!/bin/bash
#SBATCH --job-name=spline_model
#SBATCH --output=out/out_%A_%a.out
#SBATCH --error=error/error_%A_%a.err
#SBATCH --time=4:00:00
#SBATCH --array=1-100
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --mem-per-cpu=1GB
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lsf9196@nyu.edu
#SBATCH --account=torch_pr_134_general

########################################
# UPDATED: shared scratch setup + sync with $HOME
########################################
# Base scratch location (change if your cluster uses another path)
SCRATCH_BASE="${SCRATCH_BASE:-/scratch}"

# Use a single shared scratch repo (not per-array task)
PROJECT_NAME="${PROJECT_NAME:-$(basename "$SLURM_SUBMIT_DIR")}"
SCRATCH_DIR="${SCRATCH_BASE}/${USER}/${PROJECT_NAME}"
SOURCE_DIR="${SOURCE_DIR:-$HOME/${PROJECT_NAME}}"

mkdir -p "$SCRATCH_DIR"

if [ ! -d "$SOURCE_DIR" ]; then
    echo "[$(date)] ERROR: SOURCE_DIR not found: $SOURCE_DIR" >&2
    exit 1
fi

echo "[$(date)] Using SCRATCH_DIR=$SCRATCH_DIR"
echo "[$(date)] Syncing repo from $SOURCE_DIR to scratch"

# One-time sync guard to avoid 60 parallel rsyncs
SYNC_MARKER="${SCRATCH_DIR}/.sync_done"
LOCK_DIR="${SCRATCH_DIR}/.sync_lock"

if [ ! -f "$SYNC_MARKER" ]; then
    while ! mkdir "$LOCK_DIR" 2>/dev/null; do
        sleep 2
    done
    if [ ! -f "$SYNC_MARKER" ]; then
        rsync -a --delete \
        --exclude='.git' \
        --exclude='out' \
        --exclude='error' \
        "$SOURCE_DIR/" "$SCRATCH_DIR/"
        touch "$SYNC_MARKER"
    fi
    rmdir "$LOCK_DIR"
fi

# Work from the shared scratch repo root
cd "$SCRATCH_DIR"

# Ensure output dirs exist in scratch
mkdir -p out error results manifests

########################################
# EXISTING: module + env setup
########################################
module unload r

module load r/4.5.1

if [ -z "${RUN_ID:-}" ]; then
    RUN_ID="$(date -u +%Y%m%dT%H%M%SZ)_${SLURM_ARRAY_JOB_ID}"
    export RUN_ID
fi
END_YEAR=2015
export END_YEAR

########################################
# CHANGED: run R inside scratch/out
########################################
R CMD BATCH \
"--vanilla --args Surat inf_exponent \
Max_Temp 5 9 dataset_Surat_2023.csv" \
scripts/run_fitting.R > out/out_${SLURM_ARRAY_TASK_ID}.Rout

########################################
# UPDATED: sync results back to $HOME repo
########################################
echo "[$(date)] Syncing outputs back to source dir: $SOURCE_DIR"

# Per-task Rout
rsync -a out/ "$SOURCE_DIR/out/"

# If your R code writes results/manifests per task, pull them back too
rsync -a results/   "$SOURCE_DIR/results/"   || true
rsync -a manifests/ "$SOURCE_DIR/manifests/" || true

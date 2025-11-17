#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mem=2gb
#SBATCH --time=00:10:00
#SBATCH --signal=B:USR1@60
#SBATCH --job-name=c_dmtcp_pi
#SBATCH --output=c_pi_out.txt
#SBATCH --open-mode=append

# Log function for timestamped output
log () {
    echo "$(date) jobid=$SLURM_JOB_ID " "$@"
}

log "Job starting"

# Load required modules
module load spack
spack load dmtcp
module load gcc/11.2.0

# Handle job ending signal - checkpoint and requeue
handle_job_end () {
    log "Job is ending. Checkpointing..."
    dmtcp_command --bcheckpoint
    log "Checkpoint complete."
    log "Spawning new job..."
    sbatch c_monte_carlo_pi.sh
    log "Exiting..."
    exit 0
}

# Register signal handler
trap handle_job_end SIGUSR1

# Set checkpoint directory
export DMTCP_CHECKPOINT_DIR=/scratch/$USER/c-pi-checkpoint

# Check if restoring from checkpoint or starting fresh
if [ -d "$DMTCP_CHECKPOINT_DIR" ]; then
    log "Checkpoint directory exists. Restoring from checkpoint..."
    dmtcp_restart --interval 120 $DMTCP_CHECKPOINT_DIR/*.dmtcp &
    wait
else
    log "No checkpoint directory. Creating new and starting computation..."
    mkdir "$DMTCP_CHECKPOINT_DIR"
    
    # Compile the C program if not already compiled
    if [ ! -f "./monte_carlo_pi" ]; then
        log "Compiling C program..."
        gcc -o monte_carlo_pi monte_carlo_pi.c -lm
    fi
    
    log "Launching computation with DMTCP..."
    dmtcp_launch --interval 120 ./monte_carlo_pi &
    wait
fi

log "Job ran to completion!"

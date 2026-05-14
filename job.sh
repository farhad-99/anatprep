#!/bin/bash
#SBATCH --job-name=anatprep
#SBATCH --output=batch_logs/%x_%j.out
#SBATCH --time=36:00:00
#SBATCH --error=batch_logs/%x_%j.err
#SBATCH --mem=30G
#SBATCH --cpus-per-task=48

# Load required modules
module load python fsl


mkdir -p batch_logs

# Copy current working directory (codebase) to /tmp
cd /nfs/khan/trainees/ffallahz/anatprep/

# Run your training or main script (edit as needed)
echo "Running experiment..."
pixi run anatprep /nfs/trident3/mri/prado/ki3/bids/ /nfs/khan/trainees/ffallahz/ki3/derivatives/anatprep/ participant --filter_mri reconstruction:match=avgecho --template_mri DSURQE --cores all





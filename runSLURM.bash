#!/bin/bash
#SBATCH --job-name=fastq2EZbakR
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=35:00:00

module purge

eval "$(mamba shell hook --shell bash)"
mamba activate fastq2EZbakR
if [ $? -ne 0 ]; then
    echo "Failed to activate conda environment. Exiting."
    exit 1
fi
echo "Mamba environment activated successfully."

conda config --set channel_priority strict

snakemake --cores all --use-conda --conda-frontend conda --conda-prefix $HOME/local_conda_envs --rerun-incomplete --rerun-triggers mtime --verbose --keep-going
if [ $? -ne 0 ]; then
    echo "Snakemake workflow failed. Exiting."
    exit 1
fi

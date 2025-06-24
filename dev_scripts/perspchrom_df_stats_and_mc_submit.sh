#!/bin/bash

#SBATCH --job-name=chromevsim   # This is the name of the parent job that controls the subjobs in
#SBATCH --cpus-per-task=1
#SBATCH --mem=7GB
#SBATCH --time=10-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

module load snakemake/7.32.4-3.12.1
module load ete3

SNAKEFILE=/scratch/molevo/dts/odp/dev_scripts/perspchrom_df_stats_and_mc.snakefile

snakemake --snakefile ${SNAKEFILE} \
    --jobname=SM.{name}.{jobid}.sh \
    --cluster "sbatch \
    --mem={resources.mem_mb} \
    --cpus-per-task={threads} \
    --time={resources.runtime} \
    --output=log/slurm_{rule}-%A.out \
    --err=log/slurm_{rule}-%A.err" \
    --jobs 10 --rerun-incomplete

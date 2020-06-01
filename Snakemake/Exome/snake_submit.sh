#!/bin/bash

#SBATCH --job-name=snake_rnaseq_pipeline
#SBATCH --mem=1G
#SBATCH -t 120:00:00

module load python3

snakemake -s Snakefile \
--latency-wait 100 -j 1 \
--cluster 'sbatch --job-name={params.jobname} --error="{params.jobname}-%j.out" -t {params.runtime} --mem={params.mem} -c {threads}'
#--cluster-status ./find_job_status.py

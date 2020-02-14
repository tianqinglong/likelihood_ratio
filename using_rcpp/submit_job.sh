#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=25:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 16 processor core(s) per node
#SBATCH --job-name="stable genius"
#SBATCH --mail-user=qltian@protonmail.com   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module --ignore-cache load gcc
module --ignore-cache load gsl
module --ignore-cache load r

Rscript run_on_cluster.R $1 $2 $3 $4

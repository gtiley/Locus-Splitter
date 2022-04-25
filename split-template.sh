#!/bin/bash
#SBATCH --job-name=__RUNID__
#SBATCH --output=__LOGFILE__
#SBATCH --mail-user=YOUR_EMAIL
#SBATCH --mail-type=FAIL
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=2000M
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=YOUR_PARTITION_OR_DELETE
module load blast/2.10.0

#!/bin/bash
#SBATCH -p cpuq
#SBATCH --job-name=FastQC
#SBATCH --cpus-per-task=8
#SBATCH --mem=5G
#SBATCH --time=12:00:00
#SBATCH --output=qc_fastqs.%j.out
#SBATCH --error=qc_fastqs.%j.err

#my home institution does not use slurm, it uses grid engine, so I copied the above directly from
#the scripts by Lauren Mills

################################################################


#activate the relevant conda environment
source /grid/genomicscourse/home/shared/conda_2025/miniconda3/bin/activate
conda activate multiqc

#go to the right directory
cd ~/CSHL_Chromatin_Workshop_2025/data/subset/

#FastQC: ~3 -4 min total. 
#this will run for every file that ends in .fastq

fastqc *.fastq

#MultiQC: ~ 1 min
#this will run only once but use every .zip file
#and combine them into a report

multiqc *.zip

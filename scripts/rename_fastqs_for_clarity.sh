#!/bin/bash
#SBATCH -p cpuq
#SBATCH --job-name=rename_fastqs
#SBATCH --cpus-per-task=8
#SBATCH --mem=5G
#SBATCH --time=12:00:00
#SBATCH --output=fastq_rename.%j.out
#SBATCH --error=fastq_rename.%j.err

#my home institution does not use slurm, it uses grid engine, so I copied the above directly from
#the scripts by Lauren Mills

################################################################

#This will rename the fastq files for clarity
#The library names are going to be on the right column but without the .fastq file extension, which are also the contents
#of the sample.txt file
#it also creates and moves them into a data folder to keep things neat.

cd ~/CSHL_Chromatin_Workshop_2025
datadir=$(pwd)
cd ${datadir}

mkdir -p data/

mv SRR5063143.fastq data/SRR5063143_naive_H3K27ac.fastq
mv SRR5063144.fastq data/SRR5063144_naive_H3K27ac.fastq
mv SRR5063149.fastq data/SRR5063149_naive_H3K4me3.fastq
mv SRR5063150.fastq data/SRR5063150_naive_H3K4me3.fastq
mv SRR5063153.fastq data/SRR5063153_naive_input.fastq
mv SRR5063154.fastq data/SRR5063154_naive_input.fastq

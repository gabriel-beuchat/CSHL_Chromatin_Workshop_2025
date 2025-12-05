#!/bin/bash
#SBATCH -p cpuq
#SBATCH --job-name=make_index
#SBATCH --cpus-per-task=8
#SBATCH --mem=5G
#SBATCH --time=12:00:00
#SBATCH --output=index_maker.%j.out
#SBATCH --error=index_maker.%j.err

#my home institution does not use slurm, it uses grid engine, so I copied the above directly from
#the scripts by Lauren Mills

################################################################

source /grid/genomicscourse/home/shared/conda_2025/miniconda3/bin/activate
conda activate chromap
cd ~/CSHL_Chromatin_Workshop_2025/genome/
#this will build the index from the genome, and call it "index"
chromap -i -r hg38_chr22.fa -o index
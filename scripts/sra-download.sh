#!/bin/bash
#SBATCH -p cpuq
#SBATCH --job-name=donload_sras
#SBATCH --cpus-per-task=8
#SBATCH --mem=5G
#SBATCH --time=12:00:00
#SBATCH --output=down_sra.%j.out
#SBATCH --error=down_sra.%j.err

#my home institution does not use slurm, it uses grid engine, so I copied the above directly from
#the scripts by Lauren Mills

################################################################

cd ~/CSHL_Chromatin_Workshop_2025
datadir=$(pwd)
cd ${datadir}

source /grid/genomicscourse/home/shared/conda_2025/miniconda3/bin/activate
conda activate sra_tools
prefetch --option-file sra.txt

#this downloads the list of accessions in the sra.txt document from SRA, which is where a lot of
#data that gets published in papers gets deposited.
#the prefetch command will download a compressed nonsnensical version of the file which must then
#be extracted with the fasterq-dump function. To get at that function, check the split.sh script.

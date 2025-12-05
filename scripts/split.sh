#!/bin/bash
#SBATCH -p cpuq
#SBATCH --job-name=split_sras
#SBATCH --cpus-per-task=8
#SBATCH --mem=5G
#SBATCH --time=12:00:00
#SBATCH --output=divide.%j.out
#SBATCH --error=divide.%j.err

#my home institution does not use slurm, it uses grid engine, so I copied the above directly from
#the scripts by Lauren Mills

################################################################


source /grid/genomicscourse/home/shared/conda_2025/miniconda3/bin/activate
conda activate sra_tools

cd ~/CSHL_Chromatin_Workshop_2025
datadir=$(pwd)
cd ${datadir}

for SAMPLE_ID in `cat sra.txt`; do

cd ${SAMPLE_ID}

fasterq-dump --split-files ${SAMPLE_ID}.sra
mv *.fastq ../.
cd ../
rm -r ${SAMPLE_ID}
done 

#fasterq-dump will extract all the fastq files from the downloaded accessions.
#after fasterq-dump extract the files it will move them one directory up and then delete the directory
#that has the .sra file in it to save space.




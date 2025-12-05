#!/bin/bash
#SBATCH -p cpuq
#SBATCH --job-name=get_genome_sizes
#SBATCH --cpus-per-task=8
#SBATCH --mem=5G
#SBATCH --time=12:00:00
#SBATCH --output=get_sizes.%j.out
#SBATCH --error=get_sizes.%j.err

#my home institution does not use slurm, it uses grid engine, so I copied the above directly from
#the scripts by Lauren Mills

################################################################

#to convert to bam you will need to know the sizes of each chromosome.

cd ~/CSHL_Chromatin_Workshop_2025
datadir=$(pwd)
cd ${datadir}/genome
source /grid/genomicscourse/home/shared/conda_2025/miniconda3/bin/activate
conda activate basic_tools

#samtools faidx will make an index for the genome's fasta file.
samtools faidx hg38_chr22.fa
#cut is a unix utility that takes just the first and second fields in each line (in this case)
    #and outputs them to the sizes.genome file, which will be formatted:
        #chr1   12345
        #chr2   12345
        #chr3   12345
cut -f1,2 hg38_chr22.fa.fai > sizes.genome
#this just prints it to make sure it worked well.
cat sizes.genome
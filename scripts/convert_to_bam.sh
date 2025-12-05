#!/bin/bash
#SBATCH -p cpuq
#SBATCH --job-name=to_bams
#SBATCH --cpus-per-task=8
#SBATCH --mem=5G
#SBATCH --time=12:00:00
#SBATCH --output=to_bams.%j.out
#SBATCH --error=to_bams.%j.err

#my home institution does not use slurm, it uses grid engine, so I copied the above directly from
#the scripts by Lauren Mills

################################################################

#THESE ARE INCOMPLETE PATHS, YOU NEED TO ADD TO THEIR BEGINNING
datadir="~/CSHL_Chromatin_Workshop_2025"
dir=${datadir}/data/subset
genome=${datadir}/genome/hg38_chr22.fa #genome
index=${datadir}/genome/index #index
list=${datadir}/data/subset/sample.txt
size=${datadir}/genome/sizes.genome #size of chrm
##

cd $dir

source /grid/genomicscourse/home/shared/conda_2025/miniconda3/bin/activate
conda activate basic_tools

#list was previously set in the variables at the beginning as "sample.txt" and is a file that contains library names one per line.
# the backticks tells unix to run whatever is inside them and replace the spot in the line with the output.
# it then assings those to SAMPLE_ID, once for each time in the loop.
#for each of those it uses bedtools to convert them to bam.
for SAMPLE_ID in `cat $list`; do
#convert alignments to BAM
        bedtools bedtobam  -i ${SAMPLE_ID}_chromap.bed -g $size > ${SAMPLE_ID}_chromap.bam #input files are the same from chromap
done

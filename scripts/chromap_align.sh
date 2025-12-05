#!/bin/bash
#SBATCH -p cpuq
#SBATCH --job-name=align
#SBATCH --cpus-per-task=8
#SBATCH --mem=5G
#SBATCH --time=12:00:00
#SBATCH --output=align.%j.out
#SBATCH --error=align.%j.err

#my home institution does not use slurm, it uses grid engine, so I copied the above directly from
#the scripts by Lauren Mills

################################################################

#THESE ARE NOT ABSOLUTE PATHS YOU NEED TO ADD TO THEM
datadir="~/CSHL_Chromatin_Workshop_2025"
dir="${datadir}/data/subset"
genome=${datadir}/genome/hg38_chr22.fa #genome
index=${datadir}/genome/index #index
list=${datadir}/data/subset/sample.txt

source /grid/genomicscourse/home/shared/conda_2025/miniconda3/bin/activate
conda activate chromap

cd $dir
for SAMPLE_ID in `cat $list`; do
# map using chromap, output is bed file
chromap --preset chip -x $index -r $genome -q 20 --min-read-length 10   -1  ${SAMPLE_ID}.fastq  -o  ${SAMPLE_ID}_chromap.bed
done
conda activate basic_tools
#dos2unix basically just converts all the End Of Line (EOL) marker from DOS/Windows format to unix, so from \r\n to \n
#the reason this is here is just because when I tried to run "convert to bam" later, there was an error that was fixed by this.
dos2unix *.bed

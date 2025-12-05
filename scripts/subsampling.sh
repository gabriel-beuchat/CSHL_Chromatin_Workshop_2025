#!/bin/bash
#SBATCH -p cpuq
#SBATCH --job-name=subsampling
#SBATCH --cpus-per-task=8
#SBATCH --mem=5G
#SBATCH --time=12:00:00
#SBATCH --output=shrink.%j.out
#SBATCH --error=shrink.%j.err

#my home institution does not use slurm, it uses grid engine, so I copied the above directly from
#the scripts by Lauren Mills

################################################################
##select random fastq - subsample
##############################################################

cd ~/CSHL_Chromatin_Workshop_2025
datadir=$(pwd)
cd ${datadir}

##
out="subset"
##
mkdir -p data/$out

source /grid/genomicscourse/home/shared/conda_2025/miniconda3/bin/activate
cd data
conda activate sra_tools

for i in *.fastq; do
echo $i
seqtk sample -s100 $i 3000000 > $out/${i}
done

#this will "randomly" take 3 million reads from each file that matches *.fastq and deposit it
#in the $out directory with the same name.

#a note: the custom names that are more informative present in the later scripts were added by hand using the mv command.
#a script that will also do this is saved under "rename_fastqs_for_clarity.sh"
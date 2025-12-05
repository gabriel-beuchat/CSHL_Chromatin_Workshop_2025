#!/bin/bash
#SBATCH -p cpuq
#SBATCH --job-name=peak_graphs
#SBATCH --cpus-per-task=8
#SBATCH --mem=5G
#SBATCH --time=12:00:00
#SBATCH --output=peak_graphs.%j.out
#SBATCH --error=peak_graphs.%j.err

#my home institution does not use slurm, it uses grid engine, so I copied the above directly from
#the scripts by Lauren Mills

################################################################


datadir="~/CSHL_Chromatin_Workshop_2025"
dir=${datadir}/data/subset

cd ${dir}

mkdir -p consensus_matrixes/
mkdir -p deeptools_graphs/

source /grid/genomicscourse/home/shared/conda_2025/miniconda3/bin/activate
conda activate deepTools

input_combinations=(
"SRR5063143_naive_H3K27ac_treat.bam,SRR5063153_naive_input_treat.bam"
"SRR5063144_naive_H3K27ac_treat.bam,SRR5063154_naive_input_treat.bam"
"SRR5063149_naive_H3K4me3_treat.bam,SRR5063153_naive_input_treat.bam"
"SRR5063150_naive_H3K4me3_treat.bam,SRR5063154_naive_input_treat.bam"
)

for files in ${input_combinations[@]}; do

IFS=',' read -r file1 file2 <<< $files
#for more information on the above line check the call_peaks.sh script.
bamCompare -b1 ${file1} -b2 ${file2} -o ${file1//\.bam/_differential}.bw
#bamCompare is a deeptools command that makes a bigwig by essentially subtracting one bam
    #from the other.
peakfile=${file1//_treat.bam/_peaks\.Peak}
bw_file=$(echo ${file1//\.bam/_differential}.bw)
#these just set vairables to point to relevant files for ease and readability

computeMatrix scale-regions -p ${NSLOTS} -S ${bw_file} -R ${peakfile} -b 3000 -a 3000 \
        -o consensus_matrixes/${peakfile//_peaks\.Peak/\.matrix}
#this takes one or more bigwig files and one or more bed files, and essentially stakcs and scales the bed file regions horizontally
    #so they all appear the same size
    #so they are all the same width, and calculates the depth at a bin sliding across the bed region
    #the -b 3000 command tells computeMatrix to expand the bed regions upstream by 3000 bp,
    #the -a does the same thing but downstream of the bed region.
    #-o determines the output as usual.
    #one note: this script is VERY slow if not paralellized, so make sure to include a -p ${NSLOTS} argument.

plotHeatmap -m consensus_matrixes/${peakfile//_peaks\.Peak/\.matrix} -o deeptools_graphs/${peakfile//_peaks\.Peak/\.pdf} \
        --dpi 300 --startLabel "Peak Start" --endLabel "Peak End" -x "Distance" --heatmapWidth 12 --regionsLabel "Peaks"

#this guy makes the heatmap from the consensus matrix that was just made in the previous line. --dpi is dots per inch, it
    #essentially controls the quality. the others are all labels and such for customizing the figure.
    #the suffix on the -o will determine the format for output.

done

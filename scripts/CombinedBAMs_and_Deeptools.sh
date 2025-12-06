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
cd ~/CSHL_Chromatin_Workshop_2025
datadir=$(pwd)
dir="${datadir}/data/subset"
genome=${datadir}/genome/hg38_chr22.fa #genome
index=${datadir}/genome/index #index
list=${datadir}/data/subset/sample.txt

source /grid/genomicscourse/home/shared/conda_2025/miniconda3/bin/activate

#merge the bam files together to make totals:
conda activate basic_tools
cd ${dir}
#do the merging
samtools merge -f Merged_H3K27ac.bam SRR5063143_naive_H3K27ac_treat.bam SRR5063144_naive_H3K27ac_treat.bam
samtools merge -f Merged_H3K4me3.bam SRR5063149_naive_H3K4me3_treat.bam SRR5063150_naive_H3K4me3_treat.bam
samtools merge -f Merged_Input.bam SRR5063154_naive_input_treat.bam SRR5063153_naive_input_treat.bam
#re sort the files (always sort bams cause they usually need it)
samtools sort Merged_H3K27ac.bam -@ ${SLURM_CPUS_ON_NODE} -o tempfile
mv tempfile Merged_H3K27ac.bam
samtools sort Merged_H3K4me3.bam -@ ${SLURM_CPUS_ON_NODE} -o tempfile
mv tempfile Merged_H3K4me3.bam
samtools sort Merged_Input.bam -@ ${SLURM_CPUS_ON_NODE} -o tempfile
mv tempfile Merged_Input.bam


#making the indexes (always make indexes you usually need them)
samtools index Merged_H3K27ac.bam
samtools index Merged_H3K4me3.bam
samtools index Merged_Input.bam

#calling consensus peaks
conda activate macs3

file1="Merged_H3K27ac.bam"
file2="Merged_Input.bam"

macs3 callpeak  -t  $file1 -c $file2 \
        -f BAM  -g 50818468 --nomodel --shift -100 --extsize 200 \
        -n ${file1%.bam} --broad --keep-dup all \
        --outdir . 2> ${file1%.bam}_broad_macs3.log

file1="Merged_H3K4me3.bam"
file2="Merged_Input.bam"

macs3 callpeak  -t  $file1 -c $file2 \
        -f BAM  -g 50818468  --shift -100 --extsize 200 \
        -n ${file1%.bam} --keep-dup all \
        --outdir . 2> ${file1%.bam}_narrow_macs3.log

#MACS outputs the peak files in .narrowPeak or .broadPeak format.
#I would usually recommend clearly keeping the labels on these so you don't lose track,
#but for downstream convenience we will rename them today.
for peakfile in `ls *.broadPeak`
do
mv ${peakfile} ${peakfile//\.broadPeak/\.Peak}
done

for peakfile in `ls *.narrowPeak`
do
mv ${peakfile} ${peakfile//\.narrowPeak/\.Peak}
done

#display numbers of peaks by counting the file lines.
wc -l *.Peak

#make bigwigs:
conda activate deepTools
for merged_file in $(ls Merged*.bam); do \
bamCoverage -p max -b ${merged_file}  --normalizeUsing RPKM  -v  -o ${merged_file//\.bam/._norm.bw} ##you can use this on the genome browser
  #normalizing is usually a good idea so you can compare the bigwigs to each other, RPKM is the standard method people use but there are others
    #with strengths and weaknesses, google will help you here.
done

#make comparison heatmaps
#first do a bamcompare for each test
H3K4me3="Merged_H3K4me3.bam"
H3K27ac="Merged_H3K27ac.bam"
Input="Merged_Input.bam"

bamCompare -b1 ${H3K27ac} -b2 ${Input} -o ${H3K27ac//\.bam/_differential}.bw
bamCompare -b1 ${H3K4me3} -b2 ${Input} -o ${H3K4me3//\.bam/_differential}.bw

#combine peaks WITHOUT SORTING THEM
#usually you sort everything but if you do not sort them you will be able to see differences between your treatments
cut -f1,2,3 Merged*.Peak > Concat_Merged_Peaks.bed


#then compute a combined matrix for both
computeMatrix reference-point -p ${SLURM_CPUS_ON_NODE} --referencePoint "center" -S Merged*_differential.bw -R Concat_Merged_Peaks.bed -b 3000 -a 3000 \
	-o consensus_matrixes/Concat_Merged_Peaks.matrix --sortRegions "keep"

#then plot it
plotHeatmap -m consensus_matrixes/Concat_Merged_Peaks.matrix -o deepTools_graphs/Concat_Merged_Heatmap.pdf \
        --dpi 300 --startLabel "Peak Start" --endLabel "Peak End" -x "Distance" --heatmapWidth 12 --regionsLabel "Peaks" --sortRegions "keep"

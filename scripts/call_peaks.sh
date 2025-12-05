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

source /grid/genomicscourse/home/shared/conda_2025/miniconda3/bin/activate
cd ~/CSHL_Chromatin_Workshop_2025/data/subset/
conda activate macs3

input_combinations_broad=(
"SRR5063143_naive_H3K27ac_treat.bam,SRR5063153_naive_input_treat.bam"
"SRR5063144_naive_H3K27ac_treat.bam,SRR5063154_naive_input_treat.bam"
)

input_combinations_narrow=(
"SRR5063149_naive_H3K4me3_treat.bam,SRR5063153_naive_input_treat.bam"
"SRR5063150_naive_H3K4me3_treat.bam,SRR5063154_naive_input_treat.bam"
)


for files in ${input_combinations_broad[@]}; do

#the below line is handy for reading in tables, for example. In this case the "table" is generated
#in the previous lines, but you could imagine them written as a comma-separated file like a .csv
    #the <<< is basically telling unix to run whatever is to the left of it on whatever is to the right of it.
    #IFS stands for Internal Field Separator, and is just a way to tell unix the delimiter it should use to split
    #the input. read -r is used to take input from the user and assign it to a variable (or more).
    #because you are running it before the <<< you are telling read to treat $files as the user input and assign it to
    #variables $file1 and $file2, and to split it based on the IFS of ','

IFS=',' read -r file1 file2 <<< $files
##macs2-broad
macs3 callpeak  -t  $file1 -c $file2 \
        -f BAM  -g 50818468 --nomodel --shift -100 --extsize 200 \
        -n ${file1%_treat.bam} --broad \
        --outdir . 2> ${file1%_treat.bam}_broad_macs3.log
done

#macs3 will call peaks based on your bam(s). for more detail you can look at the macs manuals.
    #in this specific case we are telling it to use the human genome size with -g hs, we are telling it
    #not to make a model of the binding, because there are some assays that result in reads offset from
    #the actual site of binding for whatever you are testing for, but chip seq on histone modifications
    #is not one of them. (chip seq on transcription factors is, though)
    #instead we are going to manually shift and extend the read to make up for the fact that we are using single end
    #data instead of paired end data, so that we can have a longer insert sizes around the place we are sequencing.
    #--outdir . is just a way unix uses to say "the current directory" or in the case of moving files "the current file name."
    #the 2> *.log is telling the command to redirect the error logs to that file.

    #In this specific case, we need to provide the -c control bam, but not all assays require this. All chip seq does but
    #for example ATACseq does not. If you simply do no include a -c flag and file it will not use one.

for files in ${input_combinations_narrow[@]}; do

IFS=',' read -r file1 file2 <<< $files
##macs2-narrow
macs3 callpeak  -t  $file1 -c $file2 \
        -f BAM  -g 50818468  --shift -100 --extsize 200 \
        -n ${file1%_treat.bam}  \
        --outdir . 2> ${file1%_treat.bam}_narrow_macs3.log

done

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


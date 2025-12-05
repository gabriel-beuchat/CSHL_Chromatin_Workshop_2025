#!/bin/bash
#SBATCH -p cpuq
#SBATCH --job-name=FRiP
#SBATCH --cpus-per-task=8
#SBATCH --mem=5G
#SBATCH --time=12:00:00
#SBATCH --output=frip_job.%j.out
#SBATCH --error=frip_job.%j.err

#my home institution does not use slurm, it uses grid engine, so I copied the above directly from
#the scripts by Lauren Mills

################################################################

source /grid/genomicscourse/home/shared/conda_2025/miniconda3/bin/activate
conda activate basic_tools

#this line should return an error the first time you run it, that's okay it's just to clear any 
#existing report files since the >> appends to the file but does not replace it.
rm FRiP_Scores_Report.txt

for peakfile in `ls *.Peak`; do
#the command in the tickmarks this time lists every file that ends with .Peak

readfile=${peakfile//_peaks.Peak/_sorted_chromap.bam}
#this line calls the variable $peakfile you just made, finds _peaks.Peak within it,
    #and replaces it with _chromap_sorted.bam,
    #which are the reads that made up the peak file.

reads=$(samtools view -c ${readfile})
#samtools view -c just counts the number of reads in the sam/bam file
reads_peaks=$(bedtools intersect -u -a ${readfile} -b ${peakfile} -ubam | samtools view -c)
#the bedtools intersect tells you the reads in -a that intersect with the peaks in -b,
    #the -u and -ubam flag means it should be printed as an unsorted bam,
    #then the samtools view -c counts THOSE to get the reads in peaks
frip_score=$(echo "scale = 6; ${reads_peaks} / ${reads}" | bc)
#to get the frip score we just divide one by the other. The only complication is that
    #bash does not deal with floating point numbers (numbers that have decimals on them)
    #so a simple $reads / $reads_peaks won't work. Instead we need to pass the expression through
    #a pipe to bc, which can. we do this by building our expression and printing it with echo.
    #scale = 6 is saying we should calculate and print to 6 decimal places.

echo -e "${peakfile//_peaks.Peak/}\t${frip_score}" >> FRiP_Scores_Report.txt
#this line just prints the name of the library (so peakfile but without the _peaks.Peak ending)
    #and the frip score separated by a tab.
done
#this line just prints out the report.
cat FRiP_Scores_Report.txt
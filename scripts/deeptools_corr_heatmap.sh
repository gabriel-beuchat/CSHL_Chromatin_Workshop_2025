#!/bin/bash
#SBATCH -p cpuq
#SBATCH --job-name=correlations
#SBATCH --cpus-per-task=8
#SBATCH --mem=5G
#SBATCH --time=12:00:00
#SBATCH --output=correlations.%j.out
#SBATCH --error=correlations.%j.err

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

multiBigwigSummary bins -b *norm.bw -o bw_corr.npz -p ${NSLOTS}
#the above command takes all the bigwigs listed after -b separated by a space (automatically done here
    #with the *norm.bw), splits them into bins, figures out the pairwise correlations for each bin for each bigwig
    #averages the correlation among all the bin pairwise comparisons for each bigwig (so you get one correlation
    #value for each bigwig to bigwig comparison) and stores that in the output matrix.
    #-p ${NSLOTS} is for the number of parallel threads available to speed this calculation up.

plotCorrelation -in bw_corr.npz -c spearman -p heatmap --plotNumbers -o deeptools_graphs/correlation_heatmap.pdf
#plotCorrelation -in bw_corr.npz -c spearman -p scatterplot -o deeptools_graphs/correlation_scatterplot.pdf

#This just plots a heatmap and clusters them hierarchicaly.
#the commented out scatterplot also works and essentially plots the same thing, but I find it less useful.

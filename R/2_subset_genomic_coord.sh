#!/usr/bin/bash

## bash wrapper to subset genomic coordinates from copyKat output based on given  CNV threshold ##

## Author: Surbhi Sona ##

#SBATCH --mail-type=END
#SBATCH --mail-user=sonas@ccf.org
#SBATCH --job-name=loss
#SBATCH -n 20 
#SBATCH --array=1-17
#SBATCH --mem-per-cpu=3900

#SBATCH -o subset_gain_per_clade%J.out
#SBATCH -e subset_gain_per_clade%J.err

module load R/4.1.0

clade=$SLURM_ARRAY_TASK_ID

### n=8 ###

# Gain

#Rscript /home/sonas/beegfs/results/scRNA/bladder_copyKat_heatmap/misc/copykat_subset_genomic_coord_per_clade.R /home/sonas/beegfs/results/scRNA/bladder_copyKat_heatmap/matrix/recurrent_pairs_matrix.rds /home/sonas/beegfs/results/scRNA/bladder_copyKat_heatmap/genes/gain_loss_subset/primary_recurrent_per_clade2/ primary_recurrent gain /home/sonas/beegfs/results/scRNA/bladder_copyKat_heatmap/hclust/recurrent_pairs_Hclust.rds 10 ${clade} 

# Loss

#Rscript /home/sonas/beegfs/results/scRNA/bladder_copyKat_heatmap/misc/copykat_subset_genomic_coord_per_clade.R /home/sonas/beegfs/results/scRNA/bladder_copyKat_heatmap/matrix/recurrent_pairs_matrix.rds /home/sonas/beegfs/results/scRNA/bladder_copyKat_heatmap/genes/gain_loss_subset/primary_recurrent_per_clade2/ primary_recurrent loss /home/sonas/beegfs/results/scRNA/bladder_copyKat_heatmap/hclust/recurrent_pairs_Hclust.rds 10 ${clade}

### Sample 17 ###

# Gain

#Rscript /home/sonas/beegfs/results/scRNA/bladder_copyKat_heatmap/misc/copykat_subset_genomic_coord_per_clade.R /home/sonas/beegfs/results/scRNA/bladder_copyKat_heatmap/matrix/samples17_CNA_clean_matrix.rds /home/sonas/beegfs/results/scRNA/bladder_copyKat_heatmap/genes/gain_loss_subset/samples17_per_clade/ samples17 gain /home/sonas/beegfs/results/scRNA/bladder_copyKat_heatmap/hclust/samples17_CNA_clean_Hclust.rds 17 ${clade}

# Loss

Rscript /home/sonas/beegfs/results/scRNA/bladder_copyKat_heatmap/misc/copykat_subset_genomic_coord_per_clade.R /home/sonas/beegfs/results/scRNA/bladder_copyKat_heatmap/matrix/samples17_CNA_clean_matrix.rds /home/sonas/beegfs/results/scRNA/bladder_copyKat_heatmap/genes/gain_loss_subset/samples17_per_clade/ samples17 loss /home/sonas/beegfs/results/scRNA/bladder_copyKat_heatmap/hclust/samples17_CNA_clean_Hclust.rds 17 ${clade}




#### Archive #####

# Save cutoffs in a file and run them in parallel

export cutoff_file='/home/sonas/beegfs/results/scRNA/bladder_copyKat_heatmap/script/cutoff_g.txt'

#cutoff_g=$(cat ${cutoff_file} | nl -w1 -s ' ' | grep "^$SLURM_ARRAY_TASK_ID " | cut -f2 -d ' ')



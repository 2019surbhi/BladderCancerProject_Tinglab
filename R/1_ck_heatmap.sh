#!/usr/bin/bash

#SBATCH --mail-type=END
#SBATCH --mail-user=sonas@ccf.org
#SBATCH --job-name=final_heatmap
#SBATCH -n 40 
#SBATCH --mem-per-cpu=3900

#SBATCH -o save_mat_recurrent_pairs.out
#SBATCH -e save_mat_recurrent_pairs.err

module load R/4.1.0 

#Rscript final_copyKat_heatmap.R 6,7,35,10,36,43,44,45 recurrent_pairs 6 saved

#Rscript final_copyKat_heatmap.R 6,7,35,10,36,43,44,45 recurrent_pairs 8 saved 

Rscript final_copyKat_heatmap.R 6,7,35,10,36,43,44,45 recurrent_pairs 10 saved

#Rscript save_mat.R 6,7,35,10,36,43,44,45 recurrent_pairs 8 notsaved

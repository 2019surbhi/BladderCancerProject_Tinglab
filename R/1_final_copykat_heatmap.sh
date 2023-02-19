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

#Rscript final_copyKat_heatmap.R 6,7,35,10,36,43,44,45 recurrent_pairs 10 saved

#Rscript save_mat.R 6,7,35,10,36,43,44,45 recurrent_pairs 8 notsaved

#Rscript final_copyKat_heatmap.R 1,6,7,8,9,10,11,13,17,19,21,23,24,25,29,31,33,34,35,36,39,40,42,43,44,45,46 samples27 10 saved


sub=0.85 #(90% fails because of R vector limitation)

mode='subsampled'

out=/home/sonas/beegfs/results/scRNA/bladder_copyKat_heatmap/final_heatmap/${sub}/
mkdir ${out}

Rscript final_copyKat_heatmap.R 1,6,7,8,9,10,11,13,17,19,21,23,24,25,29,31,33,34,35,36,39,40,42,43,44,45,46 samples27 10 ${mode} ${sub}

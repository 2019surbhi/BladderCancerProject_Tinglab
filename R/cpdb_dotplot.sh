dir='/Volumes/tingalab/Surbhi/PROJECTS_tinglab_drive/BLADDER/cellphonedb/group_comparison/cpdb_outputs100/'
mean='means.txt'
pval='pvalues.txt'
input='/Users/sonas/Desktop/latest_results/cellphonedb/dotplots/dotplot_input/'
out='/Users/sonas/Desktop/latest_results/cellphonedb/dotplots/dotplots_output/'

cell_prefix='CD4_T_Central_Memory2'
cols=CD4_T_Central_Memory2_imm_cols.txt
rows=imm-imm_all_sig_int_rows.txt


# Naive
pre='naive'
fname=${pre}_${cell_prefix}_dotplot.pdf

cellphonedb plot dot_plot --means-path ${dir}/Naive/${mean} --pvalues-path ${dir}/Naive/${pval} --output-path ${out} --output-name ${fname} --rows ${input}${rows} --columns ${input}${cols} --verbose

# Recurrent
pre='recurrent'
fname=${pre}_${cell_prefix}_dotplot.pdf



cellphonedb plot dot_plot --means-path ${dir}/Recurrent/${mean} --pvalues-path ${dir}/Recurrent/${pval} --output-path ${out} --output-name ${fname} --rows ${input}${rows} --columns ${input}${cols} --verbose


### CD8 effector and Uro ###

cell_prefix='CD8_T_effector_uro_cols'
cols=CD8_T_effector_uro_cols.txt
rows=CD8_T_effector-uro_top20_int_rows.txt



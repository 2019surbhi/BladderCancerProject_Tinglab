out_dir='/Users/sonas/Documents/PROJECTS/BLADDER/RESULTS/cellphonedb/'
prefix='Naive_w'
counts=${out_dir}${prefix}_counts.txt
meta=${out_dir}${prefix}_metadata.txt
cellphonedb method statistical_analysis ${meta} ${counts} --counts-data gene_name --project-name ${prefix} --output-path ${out_dir} --verbose --pvalue 0.001

# Run cellphonedb

input_dir='/Users/sonas/Documents/PROJECTS/BLADDER/RESULTS/cellphonedb/inputs/'
out_dir='/Users/sonas/Documents/PROJECTS/BLADDER/RESULTS/cellphonedb/group_comparison/'

prefix='Naive_w'
counts=${input_dir}${prefix}/
meta=${input_dir}${prefix}/${prefix}_metadata.txt

cellphonedb method statistical_analysis ${meta} ${counts} --counts-data gene_name --project-name ${prefix} --output-path ${out_dir} --verbose --pvalue 0.001


# Generate plots

cellphonedb plot heatmap_plot ${meta} --pvalues-path ${out_dir}${prefix}/pvalues.txt  --output-path ${out_dir}${prefix}/ --count-name  ${prefix}_cellphonedb_heatmap.pdf --verbose --pvalue 0.001


# Troubleshooting: Makesafe upgrade resulted in error "ImportError: cannot import name 'soft_unicode' from 'markupsafe' (/opt/anaconda3/envs/cpdb/lib/python3.7/site-packages/markupsafe/__init__.py)" because of missing dependencies, so downgradaed thatto older version using following command:

pip install markupsafe==2.0.1

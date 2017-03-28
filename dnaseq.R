project_path = '~/projects/## PROJECT ##'
data_path = file.path(project_path, 'neolution-prep')

setwd(data_path)
source(file.path(data_path, 'ngs-tools', 'seqdata_helpers.R'))


registerDoMC(2)


# Variant calling ---------------------------------------------------------
bam_files = list.files(path = file.path(data_path, '1c_dnaseq_data/bam'),
											 pattern = '\\.bam')

# perform variant calling





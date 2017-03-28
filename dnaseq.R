project_path = '~/projects/## PROJECT ##'
data_path = file.path(project_path, 'neolution-prep')

setwd(data_path)
source(file.path(data_path, 'ngs-tools', 'seqdata_helpers.R'))


registerDoMC(2)

bam_files = list.files(path = file.path(data_path, 'bam'),
											 pattern = '\\.bam')

# perform variant calling
callGermlineAndSomaticVariants(normal_bam = ,
															 tumor_bam = ,
															 ref_genome = )

data_path = '~/projects/## PROJECT ##'

setwd(data_path)
source('seqdata_helpers.R')

registerDoMC(2)

bam_files = list.files(path = file.path(data_path, 'bam'),
											 pattern = '\\.bam')

# perform variant calling
callGermlineAndSomaticVariants(normal_bam = ,
															 tumor_bam = ,
															 ref_genome = )

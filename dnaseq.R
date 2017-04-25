project_path = '~/projects/## PROJECT ##'
data_path = file.path(project_path, 'neolution-prep')

setwd(data_path)
source(file.path(data_path, 'ngs-tools', 'seqdata_helpers.R'))


registerDoMC(2)


# Variant calling ---------------------------------------------------------
# base quality recalibration
bam_files = list.files(path = file.path(data_path, '1c_dnaseq_data/bam'),
											 pattern = '\\.bam$',
											 full.names = TRUE) %>% grep('-bq', ., invert = T, value = T)

invisible(sapply(bam_files,
								 function(file) performBaseQualityRecalibrationUsingGatk(bam = file,
								 																												bed_regions = )))

# perform variant calling
bam_files = list.files(path = file.path(data_path, '1c_dnaseq_data/bam'),
											 pattern = '\\-bq.bam$',
											 full.names = TRUE)

invisible(mapply(function(normal, tumor) callSomaticVariantsUsingGatkMutect2(normal_bam = ,
								 																									 tumor_bam = ,
								 																									 bed_regions = ),
								 bam_files[seq(1, length(bam_files), 2)],
								 bam_files[seq(2, length(bam_files), 2)])) ## check that this works out for you paired samples!

slopCoordinatesUsingBedtools(vcf = )

callGermlineVariantsUsingGatkHaplotypeCaller(normal_bam = ,
																						 n_threads = 1,
																						 bed_regions = )



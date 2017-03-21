data_path = '~/projects/## PROJECT ##'

setwd(data_path)
source('seqdata_helpers.R')
# perform quality checking
foreach(i = seq(1, length(fastq_files), 1)) %dopar% {
	performFastQP(file = fastq_files[i])
}

registerDoMC(2)

# perform alignments
fastq_files = list.files(path = file.path(data_path, 'fastq'),
												 full.names = TRUE)

foreach(i = 1:length(fastq_files)) %dopar% {
	performSTARAlignment(filename_one = fastq_files[i],
											 filename_two = NULL,
											 output_path = file.path(data_path, 'bam'),
											 quant_mode = 'salmon',
											 execute = T)
}


# perform quantifications
bam_files = list.files(path = file.path(data_path, 'bam'),
											 pattern = 'Aligned\\.toTranscriptome\\.out\\.bam',
											 recursive = TRUE,
											 full.names = TRUE)

foreach(i = 1:length(bam_files)) %dopar% {
	performSalmonQuantification(filename = bam_files[i],
															output_path = file.path(data_path, 'processed_salmon'),
															execute = T)
}


# merge ENSG identifiers to ENST in salmon output
quant_files = list.files(path = file.path(data_path, 'processed_salmon'),
												 recursive = T,
												 full.names = T,
												 pattern = '\\.sf')

quant_data = mergeEnsgInfo(quant_file = quant_files,
													 gtf_path = tool_options$salmon$gtf_annotation,
													 enst_ensg_table_path = '~/resources/ensembl_81/ensg-enst_table/Homo_sapiens.GRCh38.81.gtf_ensg-enst-table.tsv',
													 aggregate_by_ensg = TRUE)

invisible(sapply(1:length(quant_data),
								 function(idx) write.table(x = quant_data[[idx]],
								 													file = file.path(data_path, 'processed_salmon', paste(names(quant_data)[idx], 'salmon-quant-by-ensg.tsv', sep = '_')),
								 													sep = '\t',
								 													row.names = F)
))

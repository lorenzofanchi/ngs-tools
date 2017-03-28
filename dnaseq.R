project_path = '~/projects/## PROJECT ##'
data_path = file.path(project_path, 'neolution-prep')

setwd(data_path)
source(file.path(data_path, 'ngs-tools', 'seqdata_helpers.R'))


registerDoMC(2)


# Variant calling ---------------------------------------------------------
bam_files = list.files(path = file.path(data_path, '1c_dnaseq_data/bam'),
											 pattern = '\\.bam')

# perform variant calling




# Variant annotation ------------------------------------------------------
vcf_files = list.files(path = file.path(data_path, '1a_variants/vcf'),
                       pattern = '\\.vcf')

# sort VCFs, if necessary
createFastaSequenceDictionary(fasta = tool_options$salmon$fasta_dna)

invisible(sapply(vcf_files, function(file) sortVcfUsingVcfSorter(vcf = file,
                                                                 seq_dict = tool_options$general$fasta_dna_dict,
                                                                 execute = TRUE)))

# annotate variants with dbSNP ids
vcf_files = list.files(path = file.path(data_path, '1a_variants/vcf'),
                       pattern = '-sorted\\.vcf')

invisible(sapply(vcf_files, function(file) performGatkVariantAnnotation(vcf = file)))

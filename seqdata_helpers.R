## functions for processing rnaseq data
suppressPackageStartupMessages(library(doMC))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(rtracklayer))

# general parameters
tool_paths = list(general = list(gffread = 'gffread',
																 samtools = '~/libs/samtools-1.3/bin/samtools',
																 picard = '~/libs/picard-tools-2.9.0/picard.jar',
																 trimmomatic = '~/libs/Trimmomatic-0.36/trimmomatic-0.36.jar',
																 vcfsorter = '~/libs/vcfsorter/vcfsorter.pl'),

									align = list(bowtie2 = '~/libs/bowtie2-2.2.8/bowtie2',
															 mixcr = '~/libs/mixcr-2.1.1/mixcr',
															 tophat2 = '~/libs/tophat-2.1.1.Linux_x86_64/tophat2',
															 star = '~/libs/STAR-2.5.3a/bin/Linux_x86_64/STAR'),

									quality_check = list(fastqc = '~/libs/FastQC/fastqc',
																			 rseqc = '~/libs/RSeQC-2.6.3/scripts'),

									quantify = list(cufflinks = "~/libs/cufflinks-2.2.1-patched/cufflinks",
																	salmon = '~/libs/Salmon-0.8.2_linux_x86_64/bin/salmon'),

									variant_calling = list(gatk = '~/libs/GATK-3.7/GenomeAnalysisTK.jar',
									                       varscan2 = '~/libs/varscan2-2.4.3/VarScan.v2.4.3.jar')
									)

tool_options = list(general = list(parallel_threads = 18),

										cufflinks = list(gtf_annotation = '~/resources',
																		 gtf_mask_annotation = '~/resources'),

										rseqc = list(bed_reference = '~/resources/hg19_Ensembl.bed.gz'),

										salmon = list(fasta_dna = '~/resources/ensembl_81/fasta_dna/Homo_sapiens.GRCh38.81.dna.primary_assembly.fa',
																	gtf_annotation = '~/resources/ensembl_81/gtf/Homo_sapiens.GRCh38.81_no-contigs.gtf',
																	fasta_transcripts = '~/resources/ensembl_81/fasta_transcripts/Homo_sapiens.GRCh38.81.transcripts_no-contigs.fa'),

										star = list(index = '~/resources/ensembl_81/star_index',
																read_files_command = 'zcat',
																out_sam_type = 'BAM SortedByCoordinate', # options are: 'None', 'SAM', 'BAM' | For latter two, with additional: 'Unsorted', 'SortedByCoordinate'
																quant_mode = 'TranscriptomeSAM'),

										tophat2 = list(library_type = 'fr-unstranded',
																	 sam_strand_field = 'intronMotif', # for unstranded Illumina RNAseq use "intronMotif", for stranded use "None"
																	 bowtie2_index = '~/resources',
																	 bowtie2_transcriptome_index = '~/resources'),

										trimmomatic = list(phred_encoding = '-phred33',
																			 adaptor_sequences = '~/libs/Trimmomatic-0.36/adapters/Illumina_TruSeq_Adapters.fa'))


# sort bam files by name using samtools
sortBamByNameUsingSamtools = function(file, execute = TRUE) {
  command = paste(tool_paths$general$samtools,
                  'sort',
                  '-n',
                  file,
                  '-o', file.path(dirname(file), paste0(gsub('.bam', '', basename(file)), '-sorted.bam')))

  if (execute) {
    system(command = command,
           wait = TRUE)
  } else {
    message(command)
  }
}

# convert bam files to fastq using picard
convertBamToFastqUsingPicard = function(file, execute = TRUE) {
	command_samtools = paste(tool_paths$general$samtools,
													 'view',
													 '-h',
													 '-f',
													 '0x2',
													 file)

	command_picard = paste('java',
												 '-Xmx6g',
												 '-jar', tool_paths$general$picard,
												 'SamToFastq',
												 # paste0('INPUT=', file),
												 'INPUT=/dev/stdin',
												 paste0('FASTQ=', file.path(dirname(file), paste0(gsub('.bam', '', basename(file)), '_1.fastq'))),
												 paste0('SECOND_END_FASTQ=', file.path(dirname(file), paste0(gsub('.bam', '', basename(file)), '_2.fastq'))),
												 paste0('UNPAIRED_FASTQ=', file.path(dirname(file), paste0(gsub('.bam', '', basename(file)), '_unpaired.fastq')))
	)

	command = paste(command_samtools, command_picard, sep = ' | ')

  if (execute) {
    system(command = command,
           wait = TRUE)
  } else {
    message(command)
  }
}

# convert bam files to fastq using bedtools
convertBamToFastqUsingBedtools = function(file, execute = TRUE) {
	command = paste('bedtools',
									'bamtofastq',
									'-i', file,
									'-fq', file.path(dirname(file), paste0(gsub('.bam', '', basename(file)), '_1.fastq')),
									'-fq2', file.path(dirname(file), paste0(gsub('.bam', '', basename(file)), '_2.fastq'))
	)

	if (execute) {
		system(command = command,
					 wait = TRUE)
	} else {
		message(command)
	}
}


# Quality control ---------------------------------------------------------

# perform QC using FastQC
performFastQC = function(file, execute = TRUE) {
	dir.create(file.path(dirname(file),
											 gsub('\\..+$','', basename(file))), showWarnings = F)

	command = paste(tool_paths$quality_check$fastqc,
									'-t', tool_options$general$parallel_threads,
									paste0('--outdir=', file.path(dirname(file),
																								gsub('\\..+$','', basename(file)))),
									file)

	if (grepl('\\.gz$', file)) {
		command = paste('zcat', file, '|', tool_paths$quality_check$fastqc, paste0("--outdir=", file.path(dirname(file),
																																																			gsub('\\..+$','', basename(file)))), 'stdin')
	}

	if (execute) {
		system(command = command,
					 wait = TRUE)
	} else {
		message(command)
	}
}

# perform QC using FastQP
performFastQP = function(file, execute = TRUE) {
	command = paste('cd', dirname(file),';',
									tool_paths$quality_check$fastqp,
									'-o', gsub('\\..+$', '', basename(file)),
									file)

	if (execute) {
		system(command = command,
					 wait = TRUE)
	} else {
		message(command)
	}
}

# trim adapters & low Q seq using Trimmomatic
performTrimmingUsingTrimmomatic = function(filename_one, filename_two, output_path, execute = TRUE) {
  command = paste('java',
                  '-jar', tool_paths$general$trimmomatic,
                  'PE',
                  '-threads', tool_options$general$parallel_threads,
                  tool_options$trimmomatic$phred_encoding,
                  '-trimlog', file.path(output_path, paste0(basename(output_path), ".log")),
                  paste(filename_one, filename_two),
                  '-baseout', file.path(output_path, paste0(gsub('_1.fastq|_R1_001.fastq', '', basename(filename_one)), '_filtered.fq')),
                  'LEADING:3',
                  'TRAILING:3',
                  paste0('ILLUMINACLIP:', tool_options$trimmomatic$adaptor_sequences, ':2:40:15'),
                  'SLIDINGWINDOW:4:15',
                  'MINLEN:16'
                  )

  if (execute) {
    system(command = paste("nice -n 19", command),
           intern = TRUE,
           wait = TRUE)
  } else {
    message(paste("nice -n 19", command))
  }
}

# infer strandedness using RSeQC
inferStrandednessUsingRSeQC = function(bam_file, execute = TRUE) {
	command = paste('python',
									file.path(tool_paths$quality_check$rseqc, 'infer_experiment.py'),
									'-i', bam_file,
									'-r', tool_options$rseqc$bed_reference)

	if (execute) {
		system(command = paste("nice -n 19", command),
					 intern = TRUE,
					 wait = TRUE)
	} else {
		message(paste("nice -n 19", command))
	}
}


# Read alignment ----------------------------------------------------------

# align using tophat2
performTophat2Alignment = function(filename_one, filename_two, output_path, execute = TRUE) {
  command = paste(tool_paths$align$tophat2,
                  '-p', tool_options$general$parallel_threads,
                  '-m', '1',
                  '--no-novel-juncs',
                  '--no-novel-indels',
                  '--no-coverage-search',
                  '--library-type', tool_options$tophat2$library_type,
                  '-o', dirname(filename_one),
                  paste0('--transcriptome-index=', tool_options$tophat2$bowtie2_transcriptome_index),
  								tool_options$tophat2$bowtie2_index,
                  paste(filename_one, filename_two, sep = ','))

  if (execute) {
    system(command = paste("nice -n 19", command),
           intern = TRUE,
           wait = TRUE)
  } else {
    message(paste("nice -n 19", command))
  }
}

# align using star
performSTARAlignment = function(filename_one, filename_two = '', output_path, quant_mode = 'salmon', ram_limit = NULL, additional_args = NULL, execute = TRUE) {
	dir.create(file.path(output_path, gsub('_L[0-9]{3}.+|_merged.+|\\.[^.]+$', '', basename(filename_one))),
						 showWarnings = F)

  command = paste(tool_paths$align$star,
                  '--runThreadN', tool_options$general$parallel_threads,
                  '--genomeDir', tool_options$star$index,
  								if (grepl(pattern = '\\.gz$', x = filename_one)) {paste('--readFilesCommand', tool_options$star$read_files_command)},
  								'--readFilesIn', paste(filename_one, filename_two, sep = ' '),
  								'--outFileNamePrefix', file.path(output_path, gsub('_L[0-9]{3}.+|merged.+|\\.[^.]+$', '', basename(filename_one)), gsub('L[0-9]{3}.+|merged.+|\\.[^.]+$', '', basename(filename_one))),
  								'--outSAMtype', tool_options$star$out_sam_type,
  								switch(EXPR = quant_mode,
  											 'none' = '',
  											 'salmon' = paste('--quantMode', tool_options$star$quant_mode)),
  								if (!is.null(ram_limit)) {paste('--limitBAMsortRAM', format(ram_limit, scientific = F))},
  								if (!is.null(additional_args)) {paste(additional_args)}
  								)

  if (execute) {
    system(command = paste('nice -n 19', command),
           intern = FALSE,
           wait = TRUE)
  } else {
    message(paste("nice -n 19", command, '\n'))
  }
}

performMiXcrAlignment = function(filename_one, filename_two = NULL, species = 'hsa', mode = 'default', source = 'transcriptomic', execute = TRUE) {
	command = paste(tool_paths$align$mixcr,
									'align',
									'--threads', tool_options$general$parallel_threads,
									'--parameters', mode,
									'--species', species,
									if (mode == 'rna-seq') {'-OallowPartialAlignments=true'},
									if (source == 'genomic') {'-OvParameters.geneFeatureToAlign=VGeneWithP'},
									'--report', paste0(gsub('[LS][0-9]{2,3}.+', '', filename_one), 'log.txt'),
									paste(filename_one, filename_two),
									paste0(gsub('[LS][0-9]{2,3}.+', '', filename_one), 'alignments.vdjca')
	)

	if (execute) {
		system(command = paste("nice -n 19", command),
					 intern = FALSE,
					 wait = TRUE)
	} else {
		message(paste("nice -n 19", command, '\n'))
	}
}

performMiXcrContigAssembly = function(alignments, execute = TRUE) {
	command = paste(tool_paths$align$mixcr,
									'assemblePartial',
									alignments,
									paste0(gsub('\\.[^.]+$', '', alignments), '_rescued.vdjca')
	)

	if (execute) {
		system(command = paste("nice -n 19", command),
					 intern = FALSE,
					 wait = TRUE)
	} else {
		message(paste("nice -n 19", command, '\n'))
	}
}

performMiXcrAlignmentExtension = function(rescued_alignments, execute = TRUE) {
	command = paste(tool_paths$align$mixcr,
									'extendAlignments',
									rescued_alignments,
									paste0(gsub('\\.[^.]+$', '', rescued_alignments), '_extended.vdjca')
	)

	if (execute) {
		system(command = paste("nice -n 19", command),
					 intern = FALSE,
					 wait = TRUE)
	} else {
		message(paste("nice -n 19", command, '\n'))
	}
}

performMiXcrClonotypeAssembly = function(alignments, execute = TRUE) {
	command = paste(tool_paths$align$mixcr,
									'assemble',
									alignments,
									paste0(gsub('\\.[^.]+$', '', alignments), '_clones.clns')
	)

	if (execute) {
		system(command = paste("nice -n 19", command),
					 intern = FALSE,
					 wait = TRUE)
	} else {
		message(paste("nice -n 19", command, '\n'))
	}
}

performMiXcrCloneExport = function(clones, chain = NULL, execute = TRUE) {
	command = paste(tool_paths$align$mixcr,
									'exportClones',
									if (is.character(chain)) {paste('-c', chain)},
									clones,
									if (is.character(chain)) {
										paste0(gsub('\\.[^.]+$', '', clones), '_', chain, '.txt')
									} else {
										paste0(gsub('\\.[^.]+$', '', clones), '.txt')
									}
	)

	if (execute) {
		system(command = paste("nice -n 19", command),
					 intern = FALSE,
					 wait = TRUE)
	} else {
		message(paste("nice -n 19", command, '\n'))
	}
}

# Variant calling ---------------------------------------------------------

callGermlineAndSomaticVariants = function(normal_bam, tumor_bam, ref_genome) {
	if (any(!sapply(c(normal_bam, tumor_bam, ref_genome), file.exists))) {
		stop('Please check bam file and/or ref genome paths')
	}

	normal_pileup = paste(tool_paths$general$samtools,
												'mpileup',
												'-q', 1,
												'-f', ref_genome,
												normal_bam)
	tumor_pileup = paste(tool_paths$general$samtools,
												'mpileup',
												'-q', 1,
												'-f', ref_genome,
												tumor_bam)

	command = paste('java -jar', tool_paths$variant_calling$varscan2, 'somatic',
									'<(', normal_pileup,')', '<(', tumor_pileup, ')', basename(tumor_bam))

	if (execute) {
		system(command = paste("nice -n 19", command),
					 intern = FALSE,
					 wait = TRUE)
	} else {
		message(paste("nice -n 19", command, '\n'))
	}
}

# RNA quantification ------------------------------------------------------

# quantify using cufflinks
performCufflinksQuantification = function(filename, output_path, execute = TRUE) {
  command = paste(tool_paths$quantify$cufflinks,
                  "-p", tool_options$general$parallel_threads,
                  "-q",
                  "--library-type", tool_options$tophat2$library_type,
                  "-G", tool_options$cufflinks$gtf_annotation,
                  "-M", tool_options$cufflinks$gtf_mask_annotation,
                  "-o", output_path,
                  filename)

  if (execute) {
    system(command = paste("nice -n 19", command),
           intern = FALSE,
           wait = TRUE)
  } else {
    message(paste("nice -n 19", command, '\n'))
  }
}

# quantify using salmon
performSalmonQuantification = function(filename, output_path, execute = TRUE) {
	dir.create(file.path(output_path, gsub('_{0,1}Aligned.+', '', basename(filename))),
						 showWarnings = F)

	command = paste(tool_paths$quantify$salmon,
									'quant',
									'-p', tool_options$general$parallel_threads,
									'-l A',
									'-t', tool_options$salmon$fasta_transcripts,
									'-a', filename,
									'-o', file.path(output_path, gsub('_{0,1}Aligned.+', '', basename(filename))))

	if (execute) {
		system(command = paste("nice -n 19", command),
					 intern = FALSE,
					 wait = TRUE)
	} else {
		message(paste("nice -n 19", command, '\n'))
	}
}

generateTranscriptomeFasta = function(execute = F) {
	command = paste(tool_paths$general$gffread,
									'-w', tool_options$salmon$fasta_transcripts,
									'-g', tool_options$salmon$fasta_dna, tool_options$salmon$gtf_annotation)

	if (execute) {
		system(command = command,
					 intern = F,
					 wait = T)
	} else {
		message(paste("nice -n 19", command, '\n'))
	}
}

generateStarIndex = function(execute = F) {
	command = paste(tool_paths$align$star,
									'--runThreadN 32',
									'--runMode genomeGenerate',
									'--genomeDir', tool_options$star$index,
									'--genomeFastaFiles', tool_options$salmon$fasta_dna,
									'--sjdbGTFfile', tool_options$salmon$gtf_annotation)

	if (execute) {
		system(command = command,
					 intern = F,
					 wait = T)
	} else {
		message(paste("nice -n 19", command, '\n'))
	}
}

# merge ensg info to salmon output
mergeEnsgInfo = function(quant_file, enst_ensg_table_path = '', gtf_path = NULL, aggregate_by_ensg = TRUE) {
	if (file.exists(enst_ensg_table_path)) {
		message('Loading conversion table from: ', enst_ensg_table_path)
		ensg_enst_table = fread(enst_ensg_table_path)
	} else if (is.null(gtf_path)) {
		stop('ENST-ENSG path does not exist. Please supply correct path or path to GTF annotation file to generate ENST-ENSG table')
	} else if (file.exists(gtf_path) & enst_ensg_table_path != '') {
		message('Generating new ENST-ENSG conversion table from: ', gtf_path)

		gtf = readGFF(gtf_path)

		ensg_enst_table = unique(x = subset(x = gtf,
																				subset = !is.na(gtf$transcript_id),
																				select = c('gene_id', 'transcript_id')),
														 by = c('gene_id', 'transcript_id'))

		dir.create(dirname(enst_ensg_table_path), showWarnings = F)
		write.table(x = ensg_enst_table, file = enst_ensg_table_path, sep = '\t', row.names = F)
	} else {
		stop('Please supply enst_ensg_table_path & gtf_path to generate ENST-ENSG conversion table')
	}

	quant_data = setNames(object = lapply(quant_file, fread, col.names = c('transcript_id', 'transcript_length_bp', 'effective_length', 'tpm', 'read_number')),
												nm = sapply(quant_file, function(path) unlist(strsplit(x = path, split = '/'))[length(unlist(strsplit(x = path, split = '/'))) - 1], USE.NAMES = F))

	quant_data = lapply(quant_data,
											function(dt) {
												data = merge(x = dt,
																		 y = ensg_enst_table,
																		 by = 'transcript_id',
																		 all.x = TRUE)
												setcolorder(x = data,
																		neworder = c('gene_id', 'transcript_id', 'transcript_length_bp', 'effective_length', 'tpm', 'read_number'))
												return(data)
											})

	if (aggregate_by_ensg) {
		quant_data = lapply(quant_data,
												function(dt) {
													aggregate(tpm ~ gene_id, dt, sum)
												})
	}
	return(quant_data)
}


# Annotation --------------------------------------------------------------

sortVcfUsingPicard = function(vcf, seq_dict = NULL, execute = FALSE) {
  command = paste('java -jar', tool_paths$general$picard,
                  'SortVcf',
                  paste0('INPUT=', vcf),
                  if (!is.null(seq_dict)) {paste0('SEQUENCE_DICTIONARY=', seq_dict)},
                  paste0('OUTPUT=', file.path(dirname(vcf), paste0(gsub('[.][^.]+$', '', basename(vcf)), '-sorted.vcf')))
  )
  
  if (execute) {
    system(command = paste("nice -n 19", command),
           intern = FALSE,
           wait = TRUE)
  } else {
    message(paste("nohup nice -n 19", command, '&\n'))
  }
}

sortVcfUsingVcfSorter = function(vcf, seq_dict, execute = TRUE) {
  command = paste('perl ', tool_paths$general$vcfsorter,
                  seq_dict,
                  vcf,
                  '>',
                  file.path(dirname(vcf), paste0(gsub('[.][^.]+$', '', basename(vcf)), '-sorted.vcf'))
  )
  
  if (execute) {
    system(command = paste("nice -n 19", command),
           intern = FALSE,
           wait = TRUE)
  } else {
    message(paste("nice -n 19", command, '\n'))
  }
}

createFastaSequenceDictionary = function(fasta, execute = FALSE) {
  command = paste('java -jar', tool_paths$general$picard,
                  'CreateSequenceDictionary',
                  paste0('R=', fasta),
                  paste0('O=', file.path(dirname(fasta), paste0(gsub('[.][^.]+$', '', basename(fasta)), '.dict')))
  )
  
  if (execute) {
    system(command = paste("nice -n 19", command),
           intern = FALSE,
           wait = TRUE)
  } else {
    message(paste("nohup nice -n 19", command, '&\n'))
  }
}

performGatkVariantAnnotation = function(vcf, snp_db, execute = FALSE) {
  command = paste('java -jar', tool_paths$variant_calling$gatk,
                  '--analysis_type', 'VariantAnnotator',
                  '--reference_sequence', tool_options$salmon$fasta_dna,
                  '--variant', vcf,
                  '--dbsnp', snp_db,
                  '--out', file.path(dirname(vcf), paste0(gsub('[.][^.]+$', '', basename(vcf)), '-dbsnp.vcf'))
                  )
  
  if (execute) {
    system(command = paste("nice -n 19", command),
           intern = FALSE,
           wait = TRUE)
  } else {
    message(paste("nice -n 19", command, '\n'))
  }
}

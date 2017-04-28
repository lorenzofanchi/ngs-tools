## functions for processing rnaseq data
# install/load required packages
if (!require('pacman')) install.packages('pacman')
library(pacman)

required_packages = c('data.table', 'doMC', 'dplyr', 'dtplyr', 'naturalsort', 'parallel', 'pbapply', 'rtracklayer', 'stringr')

if (!p_isinstalled('rtracklayer')) {
  source('https://bioconductor.org/biocLite.R')
  biocLite('rtracklayer', suppressUpdates = T)
}

pacman::p_load(char = required_packages)

# general parameters
tool_paths = list(general = list(bedtools = '~/libs/bedtools-2.26.0/bin/bedtools',
                                 gffread = 'gffread',
                                 samtools = '~/libs/samtools-1.4/bin/samtools',
                                 picard = '~/libs/picard-tools-2.9.0/picard.jar',
                                 snpsift = '~/libs/snpEff/SnpSift.jar',
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

tool_options = list(general = list(parallel_threads = 18,

                                   fasta_dna = '~/resources/ensembl_88/fasta_dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa',
                                   fasta_dna_dict = '~/resources/ensembl_88/fasta_dna/Homo_sapiens.GRCh38.dna.primary_assembly.dict',
                                   gtf_annotation = '~/resources/ensembl_88/gtf/Homo_sapiens.GRCh38.88.gtf',
                                   fasta_transcripts = '~/resources/ensembl_88/fasta_transcripts/Homo_sapiens.GRCh38.88.transcripts.fa',

                                   snp_db = '~/resources/dbsnp_grch38/All_20161122.vcf.gz',
                                   cosmic_db = '~/resources/cosmic_80/CosmicCodingMuts.vcf.gz'),

                    cufflinks = list(gtf_annotation = '~/resources',
                                     gtf_mask_annotation = '~/resources'),

                    rseqc = list(bed_reference = '~/resources/hg19_Ensembl.bed.gz'),

                    star = list(index = '~/resources/ensembl_88/star_index',
                                read_files_command = 'zcat',
                                out_sam_type = 'BAM SortedByCoordinate', # options are: 'None', 'SAM', 'BAM' | For latter two, with additional: 'Unsorted', 'SortedByCoordinate'
                                quant_mode = 'TranscriptomeSAM'),

                    tophat2 = list(library_type = 'fr-unstranded',
                                   sam_strand_field = 'intronMotif', # for unstranded Illumina RNAseq use "intronMotif", for stranded use "None"
                                   bowtie2_index = '~/resources',
                                   bowtie2_transcriptome_index = '~/resources'),

                    trimmomatic = list(phred_encoding = '-phred33',
                                       adaptor_sequences = '~/libs/Trimmomatic-0.36/adapters/Illumina_TruSeq_Adapters.fa'))



# Convenience functions ---------------------------------------------------

# wrapper for executing commands (or not)
commandWrapper = function(command, nice = 19, intern = FALSE, wait = TRUE, nohup_out = 'nohup.out', execute) {
  if (is.numeric(nice)) {command = paste('nice -n', nice, command)}

  if (!execute) {command = paste('nohup', command, '>', nohup_out, '&\n')}

  if (execute) {
    system(command = command,
           intern = intern,
           wait = wait)
  } else {
    message(command)
  }
}

# Format conversion -------------------------------------------------------

# sort bam files by name using samtools
sortBamByNameUsingSamtools = function(file, execute = FALSE) {
  command = paste(tool_paths$general$samtools,
                  'sort',
                  '-n',
                  file,
                  '-o', file.path(dirname(file), paste0(gsub('.bam', '', basename(file)), '-sorted.bam')))

  commandWrapper(command = command, execute = execute)
}

# convert bam files to fastq using picard
convertBamToFastqUsingPicard = function(file, execute = FALSE) {
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

  commandWrapper(command = command, execute = execute)
}

# convert bam files to fastq using bedtools
convertBamToFastqUsingBedtools = function(file, execute = FALSE) {
  command = paste('bedtools',
                  'bamtofastq',
                  '-i', file,
                  '-fq', file.path(dirname(file), paste0(gsub('.bam', '', basename(file)), '_1.fastq')),
                  '-fq2', file.path(dirname(file), paste0(gsub('.bam', '', basename(file)), '_2.fastq'))
  )

  commandWrapper(command = command, execute = execute)
}


# Quality control ---------------------------------------------------------

# perform QC using FastQC
performFastQC = function(file, execute = FALSE) {
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

  commandWrapper(command = command, execute = execute)
}

# perform QC using FastQP
performFastQP = function(file, execute = FALSE) {
  command = paste('cd', dirname(file),';',
                  tool_paths$quality_check$fastqp,
                  '-o', gsub('\\..+$', '', basename(file)),
                  file)

  commandWrapper(command = command, execute = execute)
}

# trim adapters & low Q seq using Trimmomatic
performTrimmingUsingTrimmomatic = function(filename_one, filename_two, output_path, execute = FALSE) {
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

  commandWrapper(command = command, execute = execute)
}

# infer strandedness using RSeQC
inferStrandednessUsingRSeQC = function(bam_file, execute = FALSE) {
  command = paste('python',
                  file.path(tool_paths$quality_check$rseqc, 'infer_experiment.py'),
                  '-i', bam_file,
                  '-r', tool_options$rseqc$bed_reference)

  commandWrapper(command = command, execute = execute)
}


# Read alignment ----------------------------------------------------------

# align using tophat2
performTophat2Alignment = function(filename_one, filename_two, output_path, execute = FALSE) {
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

  commandWrapper(command = command, execute = execute)
}

# align using star
performSTARAlignment = function(filename_one, filename_two = '', output_path, quant_mode = 'salmon', ram_limit = NULL, additional_args = NULL, execute = FALSE) {
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

  commandWrapper(command = command, execute = execute)
}

performMiXcrAlignment = function(filename_one, filename_two = NULL, species = 'hsa', mode = 'default', source = 'transcriptomic', execute = FALSE) {
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

  commandWrapper(command = command, execute = execute)
}

performMiXcrContigAssembly = function(alignments, execute = FALSE) {
  command = paste(tool_paths$align$mixcr,
                  'assemblePartial',
                  alignments,
                  paste0(gsub('\\.[^.]+$', '', alignments), '_rescued.vdjca')
  )

  commandWrapper(command = command, execute = execute)
}

performMiXcrAlignmentExtension = function(rescued_alignments, execute = FALSE) {
  command = paste(tool_paths$align$mixcr,
                  'extendAlignments',
                  rescued_alignments,
                  paste0(gsub('\\.[^.]+$', '', rescued_alignments), '_extended.vdjca')
  )

  commandWrapper(command = command, execute = execute)
}

performMiXcrClonotypeAssembly = function(alignments, execute = FALSE) {
  command = paste(tool_paths$align$mixcr,
                  'assemble',
                  alignments,
                  paste0(gsub('\\.[^.]+$', '', alignments), '_clones.clns')
  )

  commandWrapper(command = command, execute = execute)
}

performMiXcrCloneExport = function(clones, chain = NULL, execute = FALSE) {
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

  commandWrapper(command = command, execute = execute)
}

# Variant calling ---------------------------------------------------------

performBaseQualityRecalibrationUsingGatk = function(bam,
                                                    ref_genome = tool_options$general$fasta_dna,
                                                    db_snp = tool_options$general$snp_db,
                                                    n_threads = 1,
                                                    execute = FALSE) {
  # perform base quality recalibration
  command_bq = paste('java -Xmx4g -jar', tool_paths$variant_calling$gatk,
                     '-nt', n_threads,
                     '-T', 'BaseRecalibrator',
                     '-R', ref_genome,
                     '-I', bam,
                     '-knownSites', db_snp,
                     '-o', gsub('\\.bam$', '_recal_data.table', bam))

  # Prints all reads that have a mapping quality above zero
  command_pr = paste('java -Xmx4g -jar', tool_paths$variant_calling$gatk,
                     '-nt', n_threads,
                     '-T', 'PrintReads',
                     '-R', ref_genome,
                     '-I', bam,
                     '-o', gsub('\\.bam$', '-bq.bam', bam),
                     '--read_filter', 'MappingQualityZero',
                     '--BQSR', gsub('\\.bam', '_recal_data.table', bam))

  command = paste(command_bq, command_pr, sep = ' && ')

  commandWrapper(command = command, execute = execute)
}

callGermlineVariantsUsingGatkHaplotypeCaller = function(normal_bam,
                                                        output_vcf,
                                                        ref_genome = tool_options$general$fasta_dna,
                                                        bed_regions = NULL,
                                                        bed_slop = 150,
                                                        db_snp = tool_options$general$snp_db,
                                                        n_threads = 1,
                                                        execute = FALSE) {
  if (any(!sapply(c(normal_bam, ref_genome, db_snp), file.exists))) {
    stop('Please check bam file and/or ref_genome/dbsnp/cosmic paths')
  }

  dir.create(file.path(dirname(output_vcf)),
             showWarnings = F)

  command = paste('java -Xmx8g -jar', tool_paths$variant_calling$gatk,
                  '-nt', n_threads,
                  '-T', 'HaplotypeCaller',
                  '-R', ref_genome,
                  if (!is.null(bed_regions)) { paste('-L', bed_regions) },
                  if (!is.null(bed_regions)) { paste('--interval_padding', bed_slop) },
                  '-I', normal_bam,
                  '-o', output_vcf,
                  '--dbsnp', db_snp,
                  '--maxNumHaplotypesInPopulation', 96,
                  '--dontUseSoftClippedBases',
                  '--annotateNDA')

  commandWrapper(command = command, wait = F, execute = execute)
}

callSomaticVariantsUsingGatkMutect2 = function(normal_bam,
                                               tumor_bam,
                                               output_vcf,
                                               ref_genome = tool_options$general$fasta_dna,
                                               bed_regions = NULL,
                                               bed_slop = 150,
                                               db_snp = tool_options$general$snp_db,
                                               cosmic_db = tool_options$general$cosmic_db,
                                               n_threads = 1,
                                               execute = FALSE) {
  if (any(!sapply(c(normal_bam, tumor_bam, ref_genome, db_snp, cosmic_db), file.exists))) {
    stop('Please check bam file and/or ref_genome/dbsnp/cosmic paths')
  }

  dir.create(file.path(dirname(output_vcf)),
             showWarnings = F)

  command = paste('java -Xmx8g -jar', tool_paths$variant_calling$gatk,
                  '-nt', n_threads,
                  '-T', 'MuTect2',
                  '-R', ref_genome,
                  if (!is.null(bed_regions)) { paste('-L', bed_regions) },
                  if (!is.null(bed_regions)) { paste('--interval_padding', bed_slop) },
                  '-I:tumor', tumor_bam,
                  '-I:normal', normal_bam,
                  '-U', 'ALLOW_SEQ_DICT_INCOMPATIBILITY',
                  '-o', output_vcf,
                  '--dbsnp', db_snp,
                  '--cosmic', cosmic_db,
                  '--tumor_lod', 8,
                  '--max_alt_alleles_in_normal_count', 10,
                  '--max_alt_alleles_in_normal_qscore_sum', 400,
                  '--maxNumHaplotypesInPopulation', 96,
                  '--dontUseSoftClippedBases',
                  '--annotateNDA')

  commandWrapper(command = command, wait = F, execute = execute)
}

mergeVcfs = function(vcfs = NULL, somatic_vcf = NULL, germline_vcf = NULL, qual_cutoff = 100) {
  # check which files to use
  files = c(vcfs, somatic_vcf, germline_vcf)
  if (!all(sapply(files, file.exists)) | all(is.null(c(vcfs, somatic_vcf, germline_vcf)))) { stop('Please check input file paths; some/all are missing') }

  if (length(somatic_vcf) > 1 & length(germline_vcf) > 1) {
    if (length(somatic_vcf) != length(germline_vcf)) {
      stop('Equal number of germline and somatic variants must be provided')
    } else {
      warning('Joining germline and somatic variants in order provided')
    }
  }

  # read VCFs, separate headers and variant calls
  vcfs = setNames(object = lapply(files,
                                  function(file) {
                                    all_data = readLines(file)
                                    header_data = all_data[grepl(pattern = '^#', x =  all_data)]
                                    variant_data = fread(input = paste0(all_data[!grepl(pattern = '^#', x =  all_data)], collapse = '\n'), sep = '\t', na.strings = c('NA', 'N.A.', '.', ''))

                                    setnames(x = variant_data, unlist(strsplit(x = header_data[length(header_data)], split = '\t')))

                                    if (length(names(variant_data)) == 10) {
                                      # replace germline sample name with 'NORMAL' to allow merging germline and tumor VCFs
                                      message('Replacing "', names(variant_data)[10], '" with "NORMAL" in VCF "', basename(file), '" header')
                                      setnames(x = variant_data, old = names(variant_data)[10], new = 'NORMAL')
                                      variant_data[, TUMOR := NA]

                                      # add gs_id to germline variants
                                      variant_data[, ID := paste0('gs', 1:.N, ';', ID)]
                                    }

                                    return(list(headers = header_data, variants = variant_data))
                                  }),
                  nm = basename(files))

  if (!is.null(vcfs)) {
    groups = data.table(file = names(vcfs),
                        cluster = cutree(hclust(as.dist(adist(names(vcfs)))),
                                         h = 0.1 * mean(nchar(names(vcfs)))))
  } else {
    groups = data.table(file = names(vcfs),
                        cluster = rep(1:length(somatic_vcf), 2))
  }

  # merge unfiltered vcfs
  merged_vcfs = setNames(object = lapply(groups[, unique(cluster)],
                                         function(clstr) {
                                           cluster_files = groups[cluster == clstr, file]
                                           list(headers = unique(c(grep('^##fileformat', vcfs[[cluster_files[1]]]$headers, value = T),
                                                                   unlist(sapply(cluster_files,
                                                                                 function(filename) grep('^##FILTER', vcfs[[filename]]$headers, value = T), USE.NAMES = F)),
                                                                   unlist(sapply(cluster_files,
                                                                                 function(filename) grep('^##FORMAT', vcfs[[filename]]$headers, value = T), USE.NAMES = F)),
                                                                   unlist(sapply(cluster_files,
                                                                                 function(filename) grep('^##(fileformat|FILTER|FORMAT|INFO|SAMPLE|contig|reference)|^#CHROM', vcfs[[filename]]$headers, value = T, invert = T, ignore.case = T), USE.NAMES = F)),
                                                                   unlist(sapply(cluster_files,
                                                                                 function(filename) grep('^##INFO', vcfs[[filename]]$headers, value = T), USE.NAMES = F)),
                                                                   unlist(sapply(cluster_files,
                                                                                 function(filename) grep('^##SAMPLE', vcfs[[filename]]$headers, value = T), USE.NAMES = F)),
                                                                   unlist(sapply(cluster_files,
                                                                                 function(filename) grep('^##contig', vcfs[[filename]]$headers, value = T), USE.NAMES = F)),
                                                                   unlist(sapply(cluster_files,
                                                                                 function(filename) grep('^##reference', vcfs[[filename]]$headers, value = T), USE.NAMES = F)),
                                                                   vcfs[[cluster_files[1]]]$headers[length(vcfs[[cluster_files[1]]]$headers)])),
                                                variants = rbindlist(lapply(cluster_files,
                                                                            function(filename) vcfs[[filename]]$variants), use.names = TRUE))
                                         }),
                         nm = gsub('_\\d+.vcf', '', groups[!duplicated(cluster), file]))

  dir.create(path = dirname(files[1]),
             showWarnings = F)

  # write unfiltered vcfs
  invisible(sapply(1:length(merged_vcfs),
                   function(i) {
                     writeLines(text = merged_vcfs[[i]]$headers,
                                con = if (!is.null(vcfs)) { file.path(dirname(files[1]), paste0(names(merged_vcfs)[i], '.vcf')) }
                                else { file.path(dirname(files[1]), paste0(names(merged_vcfs)[i], '-complete-unfiltered.vcf')) },
                                sep = '\n')
                     write.table(x = merged_vcfs[[i]]$variants,
                                 file = if (!is.null(vcfs)) { file.path(dirname(files[1]), paste0(names(merged_vcfs)[i], '.vcf')) }
                                 else { file.path(dirname(files[1]), paste0(names(merged_vcfs)[i], '-complete-unfiltered.vcf')) },
                                 append = T,
                                 quote = F,
                                 sep = '\t',
                                 na = '.',
                                 row.names = F,
                                 col.names = F)
                   }))

  if (is.null(vcfs)) {
    # write filtered vcfs
    invisible(sapply(1:length(merged_vcfs),
                     function(i) {
                       writeLines(text = merged_vcfs[[i]]$headers,
                                  con = file.path(dirname(files[1]), paste0(names(merged_vcfs)[i], '-complete.vcf')),
                                  sep = '\n')
                       write.table(x = merged_vcfs[[i]]$variants[(is.na(FILTER) | FILTER == 'PASS') # filter somatic variants for 'FILTER == PASS', take germline variants along with 'FILTER == NA'
                                                                 & (is.na(QUAL) | QUAL >= qual_cutoff)], # filter germline variants 'QUAL >= qual_cutoff', take somatic variants along with 'QUAL == NA'
                                   file = file.path(dirname(files[1]), paste0(names(merged_vcfs)[i], '-complete.vcf')),
                                   append = T,
                                   quote = F,
                                   sep = '\t',
                                   na = '.',
                                   row.names = F,
                                   col.names = F)
                     }))
  }
}

slopCoordinatesUsingBedtools = function(vcf, n_bases = 200, ref_genome = tool_options$general$fasta_dna, execute = TRUE) {
  command_snpsift = paste('java -jar',
                          tool_paths$general$snpsift,
                          'filter',
                          '\"FILTER = \'PASS\'\"',
                          vcf)

  command_awk = paste('cut',
                      '-f1,2',
                      '| grep -v \"#\"',
                      '| awk \'{OFS=\"\t\"; print $1,$2,$2}\'')

  command_slop = paste(tool_paths$general$bedtools,
                       'slop',
                       '-b', n_bases,
                       '-i', 'stdin',
                       '-g', paste0(ref_genome, '.fai'))

  command_merge = paste('sort -k1,1 -k2,2n |',
                        tool_paths$general$bedtools,
                        'merge',
                        '>', gsub('\\.vcf$', '.bed', vcf))

  command = paste(command_snpsift, command_awk, command_slop, command_merge, sep = ' | ')

  commandWrapper(command = command, execute = execute)
}

splitBed = function(bed, n_split = 10, execute = TRUE) {
  bed_data = fread(bed, header = F, col.names = c('chromosome', 'start', 'end'))

  split_index = round(seq(from = 1,
                          to =  length(bed_data[, unique(chromosome)]),
                          length.out = n_split + 1))

  split_index[length(split_index)] = split_index[length(split_index)] + 1

  bed_split = lapply(seq(1:(length(split_index) - 1)),
                     function(i) {
                       bed_data[grepl(paste0('^',bed_data[, unique(chromosome)][split_index[i]:(split_index[i + 1] - 1)], collapse = '|', '$'), chromosome), ]
                     })

  if (!all(unlist(sapply(bed_split, function(dt) dt[, unique(chromosome)])) %in% bed_data[, unique(chromosome)])) {
    stop('Splitting error, chromosomes missing from split bed files')
  }

  invisible(sapply(seq(1, length(bed_split)),
                   function(i) {
                     write.table(x = bed_split[[i]],
                                 file = file.path(dirname(bed), gsub('\\.bed', paste0('_', i, '.bed'), basename(bed))),
                                 append = F,
                                 quote = F,
                                 sep = '\t',
                                 row.names = F,
                                 col.names = F)
                   }))
}

# RNA quantification ------------------------------------------------------

# quantify using cufflinks
performCufflinksQuantification = function(filename, output_path, execute = FALSE) {
  command = paste(tool_paths$quantify$cufflinks,
                  "-p", tool_options$general$parallel_threads,
                  "-q",
                  "--library-type", tool_options$tophat2$library_type,
                  "-G", tool_options$cufflinks$gtf_annotation,
                  "-M", tool_options$cufflinks$gtf_mask_annotation,
                  "-o", output_path,
                  filename)

  commandWrapper(command = command, execute = execute)
}

# quantify using salmon
performSalmonQuantification = function(filename, output_path, execute = FALSE) {
  dir.create(file.path(output_path, gsub('_{0,1}Aligned.+', '', basename(filename))),
             showWarnings = F)

  command = paste(tool_paths$quantify$salmon,
                  'quant',
                  '-p', tool_options$general$parallel_threads,
                  '-l A',
                  '-t', tool_options$general$fasta_transcripts,
                  '-a', filename,
                  '-o', file.path(output_path, gsub('_{0,1}Aligned.+', '', basename(filename))))

  commandWrapper(command = command, execute = execute)
}

generateTranscriptomeFasta = function(execute = FALSE) {
  command = paste(tool_paths$general$gffread,
                  '-w', tool_options$general$fasta_transcripts,
                  '-g', tool_options$general$fasta_dna, tool_options$general$gtf_annotation)

  commandWrapper(command = command, execute = execute)
}

generateStarIndex = function(execute = FALSE) {
  command = paste(tool_paths$align$star,
                  '--runThreadN 32',
                  '--runMode genomeGenerate',
                  '--genomeDir', tool_options$star$index,
                  '--genomeFastaFiles', tool_options$general$fasta_dna,
                  '--sjdbGTFfile', tool_options$general$gtf_annotation)

  commandWrapper(command = command, execute = execute)
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

  quant_data = setNames(object = mclapply(quant_file, fread, col.names = c('transcript_id', 'transcript_length_bp', 'effective_length', 'tpm', 'read_number'), mc.cores = 10),
                        nm = sapply(quant_file, function(path) unlist(strsplit(x = path, split = '/'))[length(unlist(strsplit(x = path, split = '/'))) - 1], USE.NAMES = F))

  quant_data = pblapply(quant_data,
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
    quant_data = pblapply(quant_data,
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

  commandWrapper(command = command, execute = execute)
}

sortVcfUsingVcfSorter = function(vcf, seq_dict, execute = TRUE) {
  command = paste('perl ', tool_paths$general$vcfsorter,
                  seq_dict,
                  vcf,
                  '>',
                  file.path(dirname(vcf), paste0(gsub('[.][^.]+$', '', basename(vcf)), '-sorted.vcf'))
  )

  commandWrapper(command = command, execute = execute)
}

createFastaSequenceDictionary = function(fasta, execute = FALSE) {
  command = paste('java -jar', tool_paths$general$picard,
                  'CreateSequenceDictionary',
                  paste0('R=', fasta),
                  paste0('O=', file.path(dirname(fasta), paste0(gsub('[.][^.]+$', '', basename(fasta)), '.dict')))
  )

  commandWrapper(command = command, execute = execute)
}

performGatkVariantAnnotation = function(vcf, annotation_db = 'dbsnp', execute = FALSE) {
  command = paste('java -jar', tool_paths$variant_calling$gatk,
                  '--analysis_type', 'VariantAnnotator',
                  '--reference_sequence', tool_options$general$fasta_dna,
                  '--variant', vcf,
                  '--dbsnp', if (annotation_db == 'dbsnp') {tool_options$general$snp_db} else if (annotation_db == 'cosmic') {tool_options$general$cosmic_db},
                  '--out', file.path(dirname(vcf), paste0(gsub('[.][^.]+$', '', basename(vcf)),
                                                          if (annotation_db == 'dbsnp') {'-dbsnp.vcf'} else if (annotation_db == 'cosmic') {'-cosmic.vcf'}))
  )

  commandWrapper(command = command, execute = execute)
}

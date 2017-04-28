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
                 function(file) performBaseQualityRecalibrationUsingGatk(bam = file)))

# split bed file to perform parallel MuTect2 runs
splitBed('')

# inventorize bed & bam files
bed_files = list.files(path = file.path('~/resources/exome_bed'),
                       pattern = '_\\d+\\.bed$',
                       full.names = TRUE) %>% naturalsort(.)

bam_files = list.files(path = file.path(data_path, '1c_dnaseq_data/bam'),
                       pattern = '-bq\\.bam$',
                       full.names = TRUE) %>% naturalsort(.)

# make table with all permutations
mutect_runtable = data.table(bam_normal = bam_files[seq(1, length(bam_files), 2)], # !!!@@@ check that these are your NORMAL bam files @@@!!!
                             bam_tumor = bam_files[seq(2, length(bam_files), 2)]) # !!!@@@ check that these are your TUMOR bam files @@@!!!
mutect_runtable = cbind(mutect_runtable[rep(1:.N, length(bed_files))] %>% arrange(bam_tumor), bed_files)
mutect_runtable[, vcf_tumor := file.path(dirname(gsub('\\/bam\\/', '\\/vcf\\/', bam_tumor)),
                                         paste0(gsub('\\.bam', '', basename(bam_tumor)), '_', rep(1:(.N/length(unique(bam_tumor))), length(unique(bam_tumor))), '.vcf'))]

# perform somatic variant calling
invisible(mapply(function(bam_normal, bam_tumor, bed, vcf) {
  callSomaticVariantsUsingGatkMutect2(normal_bam = bam_normal,
                                      tumor_bam = bam_tumor,
                                      output_vcf = vcf,
                                      bed_regions = bed)
},
mutect_runtable[, bam_normal],
mutect_runtable[, bam_tumor],
mutect_runtable[, bed],
mutect_runtable[, vcf_tumor]))

# merge somatic vcfs
vcf_files = list.files(path = file.path(data_path, '1c_dnaseq_data/vcf'),
                       pattern = '_\\d+.vcf$',
                       full.names = TRUE) %>% naturalsort(.)

mergeVcfs(vcfs = vcf_files)

# slop coordinates of somatic variants, to determine search area for germline variants
vcf_files = list.files(path = file.path(data_path, '1c_dnaseq_data/vcf'),
                       pattern = '-bq\\.vcf$',
                       full.names = TRUE)

sapply(vcf_files, function(file) slopCoordinatesUsingBedtools(vcf = file))

# perform germline variant calling
bed_slop_files = list.files(path = file.path(data_path, '1c_dnaseq_data/vcf'),
                            pattern = '\\.bed$',
                            full.names = TRUE)

invisible(mapply(function(bam_normal, bed, vcf) {
  callGermlineVariantsUsingGatkHaplotypeCaller(normal_bam = bam_normal,
                                               output_vcf = vcf,
                                               bed_regions = bed)
},
mutect_runtable[, unique(bam_normal)],
bed_slop_files, # !!!@@@ make sure the order of these files, corresponds with the NORMAL bam order @@@!!!
file.path(dirname(gsub('\\/bam\\/', '\\/vcf\\/', mutect_runtable[, unique(bam_normal)])),
          basename(gsub('\\.bam', '.vcf', mutect_runtable[, unique(bam_normal)])))))

# merge somatic and germline variants
mergeVcfs(somatic_vcf = ,
          germline_vcf = )

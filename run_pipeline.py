import subprocess
import os
import sys
import pysam
import matplotlib.pyplot as plt
import pysamstats

from funcs import index_ref_genome, run_bowtie2_alignment, sam_to_bam, merge_bam_files, index_sam, filter_chromosomes, \
    mask_low_coverage, create_bed_for_low_quality_reads, call_variants, filter_bcf, calculate_snp_stats, mark_area_around_indel,\
    make_consensus, zip_files, convert_fasta_to_psmcfa, run_psmc, plot_psmc, vcf2smc, sort_bam, add_read_groups, mark_duplicates, merge_bcf, change_sample_names

# Define the paths and file names
working_dir = "/home/data1/ohad/panthera/snow_leopard_alignment/sample_12491/preprocessing"
sample_12491 = working_dir + "/variants_12491_filtered_DP12-90QUAL40.bcf"
sample_9611 = working_dir + "/variants_9611_filtered_DP12-90QUAL40.bcf"
ref = working_dir + "/reference_snow_leopard.fasta"
regions_txt = "/home/data1/ohad/panthera/ready_files/order/regions.txt"
psmc_location = "/home/data1/ohad/panthera/alignment/psmc/psmc"

samples = ["12491", "9611"]

##### index reference genome
# index_ref_genome(ref_genome)

##### align sequence to reference genome
# index_suffix = "ref_genome_index"
# run_bowtie2_alignment(pattern='*_1.fq.gz.filtered.gz', index_suffix=index_suffix,threads=50)

#### convert sam to bam
# sam_to_bam(working_dir, threads=50)

#### merge bam files
# merge_bam_files(working_dir, 50, 'merged', 'merged.bam')

##### index merged bam file
# index_sam('9611_merged_OnlyChr.bam',threads=100)

#### filter chromosomes
# filter_chromosomes('9611_merged.bam', chromosomes)

#### mark low coverage reads in bed file
# sample = 12491
# working_dir = f'/home/data1/ohad/panthera/snow_leopard_alignment/sample_{sample}/preprocessing'
# bam_file = working_dir + f'/{sample}_merged_OnlyChr.bam'
# mask_low_coverage(bam_file, min_depth=12, max_depth=80)

#### mark low quality reads
# for sample in samples:
#     working_dir = f'/home/data1/ohad/panthera/snow_leopard_alignment/sample_{sample}/preprocessing'
#     bam_file = working_dir + f'/{sample}_merged_OnlyChr.bam'
#     create_bed_for_low_quality_reads(bam_file, working_dir,rms_threshold=25)


##### call variants
# call_variants('12491_merged_OnlyChr.bam', ref_genome, 'variants_12491.bcf', threads=80)

##### filter bcf
# bcf_file = working_dir + '/variants_12491.bcf'
# filter_bcf(bcf_file, min_qual=40, min_depth=12,max_depth=90)

#### mark areas around indels
for sample in samples:
    working_dir = f'/home/data1/ohad/panthera/snow_leopard_alignment/sample_{sample}/preprocessing'
    bcf_file = working_dir + f'/variants_{sample}.bcf'
    output_file = 'mark_indels_area.bed'
    mark_area_around_indel(bcf_file,working_dir, output_file)




# # Calculate initial SNP stats
# for sample in samples:
#     input_bcf = f"/home/data1/panthera/ready_files/order/{sample}_no_filter.bcf"
#     print(f"Initial SNP statistics {sample}:")
#     calculate_snp_stats(input_bcf)
#



#######filter bcf files for quality and depth and chromosome X and MT removed
# for sample in samples:
#     input_bcf = f"/home/data1/panthera/ready_files/order/{sample}_no_filter.bcf"
#     output_filtered_bcf = f"/home/data1/panthera/ready_files/order/{sample}_filtered.bcf"
#     regions_txt = f"/home/data1/panthera/ready_files/order/regions.txt"
#     output_sex_mt_removed_bcf = f"/home/data1/panthera/ready_files/order/{sample}_filtered_autosom.bcf"
#
#     print(f"Processing {sample}...")
#
#     # Filter VCF
#     # filter_vcf(input_bcf, output_filtered_bcf)
#
#     # Remove sex and mitochondrial chromosomes
#     remove_sex_mt_chromosomes(output_filtered_bcf, output_sex_mt_removed_bcf, regions_txt)
#
#     print(f"Finished processing {sample}.\n")
#
######### Index the BCF files
# index_bcf(sample_9611)
# index_bcf(sample_12491)


############### PSMC preprocessing and analysis
############### change directory to the psmc directory

# mask_file = working_dir + "/mask_regions.bed"
# for sample in samples:
#     input_bcf = working_dir + f"/variants_{sample}_filtered_DP12-90QUAL40.bcf"
#     make_consensus(input_bcf, ref, mask_file)
#     print(f"Finished processing {sample}.\n")

# how soft masking changes the results?
# how het looks in the fasta?



#################### zip files
# for sample in samples:
#     input = f"/home/data1/panthera/ready_files/order/{sample}_filtered_autosom_consensus.fa"
#     zip_files(input)
#     print(f"Finished zipping {sample}.\n")

########## convert fasta to psmcfa as the required format for psmc
# for sample in samples:
#     input_fa = f"/home/data1/panthera/ready_files/order/{sample}_filtered_autosom_consensus.fa.gz"
#     convert_fasta_to_psmcfa(input_fa)
#     print(f"Finished converting {sample}.\n")

######### run psmc
# for sample in samples:
#     input_psmcfa = f"/home/data1/ohad/panthera/ready_files/order/psmc/{sample}_filtered_autosom_consensus.psmcfa"
#     run_psmc(input_psmcfa,p='1+1+1+1+1+1+1+1+1+1+20+24+10', t=15, r=5, N=25)
#     print(f"Finished running PSMC for {sample}.\n")

########## plot psmc
# for sample in samples:
#     input_psmc = f"/home/data1/ohad/panthera/ready_files/order/psmc/{sample}_filtered_autosom_consensus_10-1.psmc"
#     plot_psmc(input_psmc)
#     print(f"Finished plotting PSMC for {sample}.\n")


####### MSMC
# working_dir = "/home/data1/ohad/panthera/ready_files/order/msmc"
# input_vcf = os.path.join(working_dir, 'merged_filtered_autosom_name.vcf.gz')
# output_smc = os.path.join(working_dir, 'merged_filtered_autosom_A2.smc.gz')
#
# vcf2smc(working_dir=working_dir,
#         input_vcf="/home/data1/ohad/panthera/ready_files/order/msmc/merged_filtered_autosom_name.vcf.gz",
#         output_smc="/home/data1/ohad/panthera/ready_files/order/msmc/merged_filtered_autosom_A2.smc.gz",
#         samples= 'sample_9611 sample_12491',
#         missing_cutoff=5000,
#         chr='A2',
#         pop_samples='sample_9611,sample_12491')














##########miscellaneous
##########sort bam
# for sample in samples:
#     input_bam = "/home/data1/panthera/ready_files/order/merged_{sample}.bam
#     sort_bam(input_bam)
#     print(f"Finished sorting BAM for {sample}.\n")

########add read groups
# for sample in samples:
#     input_bam = f"/home/data1/panthera/ready_files/order/merged_{sample}.bam"
#     add_read_groups(input_bam)
#     print(f"Finished adding read groups for {sample}.\n")

##########mark duplicates
# for sample in samples:
#     input_bam = f"/home/data1/panthera/ready_files/order/merged_{sample}_rg.bam"
#     mark_duplicates(input_bam)
#     print(f"Finished marking duplicates for {sample}.\n")

#########merge bcf files
# sample_12491 = f"/home/data1/ohad/panthera/ready_files/order/msmc/9611_filtered_autosom.bcf"
# sample_9611 = f"/home/data1/ohad/panthera/ready_files/order/msmc/12491_filtered_autosom.bcf"
# output = f"/home/data1/ohad/panthera/ready_files/order/msmc/merged_filtered_autosom.vcf.gz"
# merge_bcf(sample_12491, sample_9611, output, output_type='z')

#########change sample names
# sample_names = f"home/data1/ohad/panthera/ready_files/order/msmc/sample_names.txt"
# input_bcf = f"/home/data1/ohad/panthera/ready_files/order/msmc/merged_filtered_autosom.vcf.gz"
# output_bcf = f"/home/data1/ohad/panthera/ready_files/order/msmc/merged_filtered_autosom_renamed.vcf.gz"
# change_sample_names(input_bcf=input_bcf, output_bcf= output_bcf, sample_names=sample_names)



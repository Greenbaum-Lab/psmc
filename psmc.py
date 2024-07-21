import subprocess
import os
this is a test
# Define the paths and file names
sample_12491 = "/home/data1/ohad/panthera/ready_files/order/12491_no_filter.bcf"
sample_9611 = "/home/data1/ohad/panthera/ready_files/order/9611_no_filter.bcf"
output_filtered_bcf = "/home/data1/ohad/panthera/ready_files/order/12491_filtered.bcf"
regions_txt = "/home/data1/ohad/panthera/ready_files/order/regions.txt"
output_sex_mt_removed_bcf = "/home/data1/ohad/panthera/ready_files/order/12491_sex_mt_removed.bcf"
ref = "/home/data1/ohad/panthera/ready_files/order/felis_catus_ref.fa"
psmc_location = "/home/data1/ohad/panthera/alignment/psmc/psmc"

def run_command(command):
    """Run a shell command and capture the output."""
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        raise Exception(f"Error running command: {command}\n{stderr.decode()}")
    return stdout.decode()

def calculate_snp_stats(bcf_file):
    """Calculate the number of total SNPs and the percentage of heterozygous SNPs."""
    # Count total SNPs
    total_snps_cmd = f"bcftools view -v snps -H {bcf_file} | wc -l"
    total_snps = int(run_command(total_snps_cmd).strip())

    # Count heterozygous SNPs
    het_snps_cmd = f"bcftools view -v snps -g het -H {bcf_file} | wc -l"
    het_snps = int(run_command(het_snps_cmd).strip())

    # Calculate heterozygous percentage
    heterozygous_percentage = (het_snps / total_snps) * 100 if total_snps > 0 else 0

    print(f"Total SNPs: {total_snps}")
    print(f"Heterozygous SNPs: {het_snps} ({heterozygous_percentage:.2f}%)")

def filter_vcf(input_bcf, output_bcf,q=50,dp_min=10,dp_max=50,keep='snps'):
    """Filter VCF based on quality and depth."""
    filter_cmd = f"bcftools filter -i 'QUAL > {q} & INFO/DP > {dp_min} & INFO/DP < {dp_max}' {input_bcf} | bcftools view -v {keep} -o {output_bcf} -O b"
    run_command(filter_cmd)
    print(f"Filtered VCF saved to {output_bcf}")

def create_regions_file(input_bcf, regions_txt):
    """Create regions.txt excluding X and MT chromosomes"""
    create_regions_cmd = f"bcftools view -h {input_bcf} | grep '^##contig' | grep -v 'ID=X' | grep -v 'ID=MT' | awk -F'[=,]' '{{print $3}}' > {regions_txt}"
    run_command(create_regions_cmd)
    print(f"Regions file created at {regions_txt}")

def keep_autosom(input_bcf, output_bcf, regions_txt):
    """Remove sex (X) and mitochondrial (MT) chromosomes."""
    chr = "^X, MT"
    # Filter the VCF to remove X and MT chromosomes
    filter_regions_cmd = f"bcftools view -t {chr} -o {output_bcf} -O b {input_bcf}"
    run_command(filter_regions_cmd)
    print(f"Sex and mitochondrial chromosomes removed. Output saved to {output_bcf}")

def index_bcf(bcf_file):
    """Index the BCF file."""
    index_cmd = f"bcftools index {bcf_file}"
    run_command(index_cmd)
    print(f"BCF file indexed: {bcf_file}")


def make_consensus(input_bcf, ref):
    """Create a consensus FASTA sequence from a BCF file."""
    output_fasta = input_bcf.replace(".bcf", "_consensus.fa")
    consensus_cmd = f"bcftools consensus -f {ref} -I {input_bcf} -o {output_fasta}"
    run_command(consensus_cmd)
    print(f"Consensus FASTA sequence saved to {output_fasta}")

def zip_files(input):
    """Zip files"""
    zip_cmd = f"gzip {input}"
    run_command(zip_cmd)
    print(f"Files zipped")


def convert_fasta_to_psmcfa(input_fa):
    """Convert FASTA to PSMC format."""
    output_psmcfa = input_fa.replace(".fa", ".psmcfa")
    convert_cmd = f"/home/data1/panthera/alignment/psmc/utils/fq2psmcfa {input_fa} > {output_psmcfa}"
    run_command(convert_cmd)
    print(f"PSMC format saved to {output_psmcfa}")


# =4+25*2+4+6
def run_psmc(input_psmcfa, p, t=15, r=5, N=25):
    """Run PSMC. output should have 10 recombinatio nevents per time interval"""
    # Get the base name of the input file
    base_name = os.path.basename(input_psmcfa)
    # Replace the '.psmcfa' extension with '.psmc' in the base name
    output_base_name = base_name.replace('.psmcfa', '_10-1.psmc')
    # Join the current directory with the 'psmc' directory and the output base name to get the output file path
    output_psmc = os.path.join(os.getcwd(), 'psmc', output_base_name)
    psmc_cmd = f"/home/data1/ohad/panthera/alignment/psmc/psmc -N{N} -t{t} -r{r} -p '{p}' -o {output_psmc} {input_psmcfa}"
    run_command(psmc_cmd)
    print(f"PSMC output saved to {output_psmc}")


def plot_psmc(input_psmc,g=5,u=1.1e-9):
    """Plot the PSMC results."""
    output_plot = input_psmc.replace('.psmc', '_plot')
    plot_cmd = f"/home/data1/ohad/panthera/alignment/psmc/utils/psmc_plot.pl " \
               f"-u {u} -g {g} {output_plot} {input_psmc}"
    run_command(plot_cmd)
    print("PSMC plot generated.")

def sort_bam(input_bam):
    """Sort the BAM file."""
    output_bam = input_bam.replace(".bam", "_sorted.bam")
    sort_cmd = f"samtools sort -o {output_bam} {input_bam}"
    run_command(sort_cmd)
    print(f"BAM file sorted: {output_bam}")


def mark_duplicates(input_bam,remove='true'):
    """Mark duplicates in a sorted BAM file."""
    output_bam = input_bam.replace(".bam", "_rmdup.bam")
    output_metrics = input_bam.replace(".bam", "_metrics.txt")
    mark_cmd = f"picard MarkDuplicates I={input_bam} O={output_bam}" \
               f" M={output_metrics} REMOVE_DUPLICATES={remove} CREATE_INDEX=true"
    run_command(mark_cmd)
    print(f"Duplicates marked in BAM file: {output_bam}")


def add_read_groups(input_bam):
    """Add or replace read groups in a BAM file."""
    output_bam = input_bam.replace(".bam", "_rg.bam")
    add_rg_cmd = f"picard AddOrReplaceReadGroups I={input_bam} O={output_bam}" \
                 f" RGID=ID_1 RGLB=LB_1  RGPL=illumina RGPU=unit1 RGSM=SM_1"
    run_command(add_rg_cmd)
    print(f"Read groups added to BAM file: {output_bam}")


def merge_bcf(sample1, sample2, output_name, output_type='b' or 'z' or 'v' or 'u'):
    """Merge two BCF files."""
    merge_cmd = f"bcftools merge --merge 'none' --missing-to-ref -O {output_type} -o {output_name} {sample1} {sample2}"
    run_command(merge_cmd)
    print(f"BCF files merged: {output_name}")


def change_sample_names(input_bcf, output_bcf, sample_names):
    """Change sample names in a vcf/bcf file."""
    change_names_cmd = f"bcftools reheader -s {sample_names} -o {output_bcf} {input_bcf}"
    run_command(change_names_cmd)
    print(f"Sample names changed. Output saved to {output_bcf}")




# smc++ vcf2smc -d sample_9611 sample_12491 --missing-cutoff=100000 merged_filtered_autosom_name.vcf.gz merged_filtered_autosom_A2.smc.gz A2 :sample_9611,sample_12491

def vcf2smc(working_dir, input_vcf, output_smc, samples, missing_cutoff,chr, pop_samples):
    """Convert a VCF file to SMC++ format.
    Args:"""
    # Change directory and print the current working directory
    change_directory = f"cd {working_dir} && pwd"
    current_dir = run_command(change_directory)
    print(current_dir)

    vcf2smc_cmd = f"smc++ vcf2smc -d {samples} --missing-cutoff={missing_cutoff} {input_vcf} {output_smc} {chr} pop:{pop_samples}"
    run_command(vcf2smc_cmd)
    print(f"SMC++ format saved to {output_smc}")




samples = ["12491", "9611"]



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
########## Index the BCF files
# sample_12491 = "/home/data1/panthera/ready_files/order/12491_filtered_autosombcf"
# sample_9611 = "/home/data1/panthera/ready_files/order/9611_filtered_autosom.bcf"
# index_bcf(sample_9611)
# index_bcf(sample_12491)


###############make consensus sequence
# for sample in samples:
#     input_bcf = f"/home/data1/panthera/ready_files/order/{sample}_filtered_autosom.bcf"
#     make_consensus(input_bcf, ref)
#     print(f"Finished processing {sample}.\n")

####################zip files
# for sample in samples:
#     input = f"/home/data1/panthera/ready_files/order/{sample}_filtered_autosom_consensus.fa"
#     zip_files(input)
#     print(f"Finished zipping {sample}.\n")

##########convert fasta to psmcfa as the required format for psmc
# for sample in samples:
#     input_fa = f"/home/data1/panthera/ready_files/order/{sample}_filtered_autosom_consensus.fa.gz"
#     convert_fasta_to_psmcfa(input_fa)
#     print(f"Finished converting {sample}.\n")

#########run psmc
# for sample in samples:
#     input_psmcfa = f"/home/data1/ohad/panthera/ready_files/order/psmc/{sample}_filtered_autosom_consensus.psmcfa"
#     run_psmc(input_psmcfa,p='1+1+1+1+1+1+1+1+1+1+20+24+10', t=15, r=5, N=25)
#     print(f"Finished running PSMC for {sample}.\n")

##########plot psmc
# for sample in samples:
#     input_psmc = f"/home/data1/ohad/panthera/ready_files/order/psmc/{sample}_filtered_autosom_consensus_10-1.psmc"
#     plot_psmc(input_psmc)
#     print(f"Finished plotting PSMC for {sample}.\n")

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

# msmc
# need to phase data (seperate into maternal and paternal haplotypes)


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


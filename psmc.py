import subprocess
import os

# Define the paths and file names
sample_12491 = "/home/data1/panthera/ready_files/order/12491_no_filter.bcf"
sample_9611 = "/home/data1/panthera/ready_files/order/9611_no_filter.bcf"
output_filtered_bcf = "/home/data1/panthera/ready_files/order/12491_filtered.bcf"
regions_txt = "/home/data1/panthera/ready_files/order/regions.txt"
output_sex_mt_removed_bcf = "/home/data1/panthera/ready_files/order/12491_sex_mt_removed.bcf"
ref = "/home/data1/panthera/ready_files/order/felis_catus_ref.fa"


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

def remove_sex_mt_chromosomes(input_bcf, output_bcf, regions_txt):
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



def run_psmc(input_psmcfa):
    """Run PSMC."""
    # Get the base name of the input file
    base_name = os.path.basename(input_psmcfa)
    # Replace the '.psmcfa' extension with '.psmc' in the base name
    output_base_name = base_name.replace('.psmcfa', '.psmc')
    # Join the current directory with the 'psmc' directory and the output base name to get the output file path
    output_psmc = os.path.join(os.getcwd(), 'psmc', output_base_name)

    psmc_cmd = f"/home/data1/panthera/alignment/psmc/psmc -N25 -t15 -r5 -p '4+25*2+4+6' -o {output_psmc} {input_psmcfa}"
    run_command(psmc_cmd)
    print(f"PSMC output saved to {output_psmc}")


# # Calculate initial SNP stats
# sample_12491 = "/home/data1/panthera/ready_files/order/12491_no_filter.bcf"
# sample_9611 = "/home/data1/panthera/ready_files/order/9611_no_filter.bcf"
# print("Initial SNP statistics 12491:")
# calculate_snp_stats(sample_12491)
#
# # Calculate initial SNP stats
# print("Initial SNP statistics 9611:")
# calculate_snp_stats(sample_9611)

#######filter bcf files for quality and depth and chromosome X and MT removed
samples = ["12491", "9611"]
# for sample in samples:
#     input_bcf = f"/home/data1/panthera/ready_files/order/{sample}_no_filter.bcf"
#     output_filtered_bcf = f"/home/data1/panthera/ready_files/order/{sample}_filtered.bcf"
#     regions_txt = f"/home/data1/panthera/ready_files/order/regions.txt"
#     output_sex_mt_removed_bcf = f"/home/data1/panthera/ready_files/order/{sample}_filtered_sex_mt_removed.bcf"
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
# sample_12491 = "/home/data1/panthera/ready_files/order/12491_filtered_sex_mt_removed.bcf"
# sample_9611 = "/home/data1/panthera/ready_files/order/9611_filtered_sex_mt_removed.bcf"
# index_bcf(sample_9611)
# index_bcf(sample_12491)


###############make consensus sequence
# for sample in samples:
#     input_bcf = f"/home/data1/panthera/ready_files/order/{sample}_filtered_sex_mt_removed.bcf"
#     make_consensus(input_bcf, ref)
#     print(f"Finished processing {sample}.\n")

####################zip files
# for sample in samples:
#     input = f"/home/data1/panthera/ready_files/order/{sample}_filtered_sex_mt_removed_consensus.fa"
#     zip_files(input)
#     print(f"Finished zipping {sample}.\n")

##########convert fasta to psmcfa
# for sample in samples:
#     input_fa = f"/home/data1/panthera/ready_files/order/{sample}_filtered_sex_mt_removed_consensus.fa.gz"
#     convert_fasta_to_psmcfa(input_fa)
#     print(f"Finished converting {sample}.\n")

##########run psmc
for sample in samples:
    input_psmcfa = f"/home/data1/panthera/ready_files/order/{sample}_filtered_sex_mt_removed_consensus.psmcfa"
    run_psmc(input_psmcfa)
    print(f"Finished running PSMC for {sample}.\n")

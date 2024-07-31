import subprocess
import os
import sys
import pysam
import matplotlib.pyplot as plt
import pysamstats

working_dir = '/home/data1/ohad/panthera/snow_leopard_alignment/sample_12491/test_chr'
ref_genome = f"{working_dir}/reference_snow_leopard.fasta"
bcf_file = working_dir + '/chr_default.bcf'
# chromosoms_list = working_dir + '/chromosomes.txt'
# with open(chromosoms_list, 'r') as f:
#     chromosomes = f.read().splitlines()


def run_command(command,directory=working_dir):
    """Run a shell command, capture the output, and print it."""
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True,cwd=directory)
    stdout, stderr = process.communicate()

    # Decode and print stdout and stderr
    stdout_decoded = stdout.decode()
    stderr_decoded = stderr.decode()
    print(stdout_decoded)
    if stderr_decoded:
        print(stderr_decoded, file=sys.stderr)

    if process.returncode != 0:
        raise Exception(f"Error running command: {command}\n{stderr_decoded}")
    return stdout_decoded


def index_ref_genome(ref_genome):
    """Index the reference genome for alignment."""
    output_path = os.path.dirname(ref_genome)
    align_cmd = f"bowtie2-build -f {ref_genome} {output_path}/ref_genome_index --threads 30"
    run_command(align_cmd)
    print("Reference genome indexed.")


def run_bowtie2_alignment(pattern: str, index_suffix: str,threads=50):
    """Align paired-end reads using Bowtie2.
    Args:
        pattern: A glob pattern to match the FASTQ files. e.g. "*_1.fq.gz.filtered.gz"
    """

    # Use the `ls` command to list files matching the pattern and split the output to get a list of file paths
    fastq_files_1 = run_command(f'ls {pattern}')
    fastq_files_1 = fastq_files_1.strip().split('\n')

    # Adjust the pattern to find the second set of files and repeat the process
    pattern2 = pattern.replace("_1", "_2")
    fastq_files_2 = run_command(f'ls {pattern2}')
    fastq_files_2 = fastq_files_2.strip().split('\n')

    # Iterate over each pair of files
    for fq1, fq2 in zip(fastq_files_1, fastq_files_2):
        output = os.path.join(working_dir, f"{os.path.basename(fq1).replace('_1.fq.gz.filtered.gz', '_aligned.sam')}")
        print(f"Aligning {fq1} and {fq2} to {output}")

        align_cmd = f"bowtie2 -x {index_suffix} -1 {fq1} -2 {fq2} -S {output} --threads {threads}"
        run_command(align_cmd)


def index_sam(sam_file, threads=50):
    """Index the BCF file."""
    index_cmd = f"samtools index {sam_file} -@ {threads}"
    run_command(index_cmd)
    print(f"SAM file indexed: {sam_file}")


def sam_to_bam(working_dir, threads):
    """Convert SAM files to BAM files. Then sort by coordinate and index the BAM files."""
    sam_files = run_command(f'ls *.sam')
    sam_files = sam_files.strip().split('\n')

    for sam_file in sam_files:
        # Define the output BAM file name (replace .sam with .sorted.bam)
        bam_file = f"{sam_file.replace('.sam', '_sorted.bam')}"
        print(f"Converting and sorting {sam_file} to {bam_file}...")
        # Convert SAM to BAM and sort
        convert = f"samtools view -bS {sam_file} -@ {threads}| samtools sort -o {bam_file} -@ {threads}"
        run_command(convert)

        index_cmd = f"samtools index {bam_file} -@ {threads}"
        run_command(index_cmd)
        print(f"BAM file indexed: {bam_file}")



def merge_bam_files(working_dir, threads,output, bam_pattern):
    """Merge BAM files into a single file."""
    bam_files = run_command(f'ls *{bam_pattern}')
    bam_files = bam_files.strip().replace('\n', ' ')
    print(bam_files)
    # Define the output BAM file name
    merged_bam = f"{output}.bam"
    merge_cmd = f"samtools merge -@ {threads} {merged_bam} *{bam_files}"
    run_command(merge_cmd)
    print(f"BAM files merged into {merged_bam}")



def mask_low_coverage(bam_file, min_depth=12, max_depth=80):
    """Mask regions with low coverage in the BAM file."""
    # Define the output BED file name
    output_bed = "mark_low_coverage.bed"
    mask_cmd = f"bedtools genomecov -ibam {bam_file} -bga | awk '$4 < {min_depth} || $4 > {max_depth}' > {output_bed}"
    run_command(mask_cmd)
    print(f"Low coverage regions masked: {output_bed}")

def create_bed_for_low_quality_reads(bam_file, working_dir, rms_threshold: 25):
    bam_file = pysam.AlignmentFile(bam_file)
    with open(f"{working_dir}/mark_RMS.bed", 'a') as mask_rms:
        # create bed file tab delimited for low quality reads
        for record in pysamstats.stat_mapq(bam_file):
            if record['rms_mapq'] < rms_threshold:
                mask_rms.write(f"{record['chrom']}\t{record['pos']}\t{record['pos']+1}\n")


def filter_chromosomes(input_bam, chromosomes):
    """Filter the input BAM file to include only the chromosomes in the list.
    Args:chrosmos: list of chromosomes to keep in txt format"""
    # Create individual BAM files for each chromosome
    output_bam = f"{input_bam.replace('.bam', '_OnlyChr.bam')}"
    chromosomes = ' '.join(chromosomes)
    print(f"Filtering chromosomes in {input_bam} to {output_bam}...")
    filter_cmd = f"samtools view -b {input_bam} {chromosomes} -o {output_bam} -@ 50"
    run_command(filter_cmd)


def call_variants(bam_file, ref_genome, output_bcf, threads=50):
    """Call variants using bcftools mpileup and call."""
    # Define the output BCF file name
    mpileup_cmd = (f"bcftools mpileup -Ou --threads {threads} -a 'FORMAT/DP,FORMAT/AD,INFO/AD' "
                   f"-f {ref_genome} {bam_file} | bcftools call --threads {threads} -mv -Ob -o {output_bcf}")
    run_command(mpileup_cmd)
    print(f"Variants called: {output_bcf}")



def mark_area_around_indel(bcf_file, working_dir, output_file):
    """Mark 10bp around each indel in the BCF file.
    need to left align (normalize) the indels first"""

    vcf = pysam.VariantFile(bcf_file, "rb")

    with open(f'{working_dir}/{output_file}', 'w') as indel_pos:
        for record in vcf:
            if record.info.get('INDEL', False):

                indel_pos.write(f"{record.chrom}\t{record.pos-10}\t{record.pos-1}\n")
                indel_pos.write(f"{record.chrom}\t{record.pos+1}\t{record.pos+10}\n")
    vcf.close()




def filter_bcf(input_bcf, min_qual=40, min_depth=12,max_depth=90):
    """Filter the input BCF file based on quality and depth."""
    output_bcf_path = f"{input_bcf.replace('.bcf', f'_filtered_DP{min_depth}-{max_depth}QUAL{min_qual}.bcf')}"
    input_bcf = pysam.VariantFile(input_bcf,threads=50)
    output_bcf = pysam.VariantFile(output_bcf_path, 'w', header=input_bcf.header)
    print(f"Filtering {input_bcf} to {output_bcf_path}...")
    for record in input_bcf:
        if record.qual >= min_qual and min_depth <= record.info['DP'] <= max_depth:
            output_bcf.write(record)
    output_bcf.close()
    print(f"BCF file filtered: {output_bcf_path}")


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


def make_consensus(input_bcf, ref, mask_file):
    """Create a consensus FASTA sequence from a BCF file.
    I=IUPAC code for all genotypes"""
    output_fasta = input_bcf.replace(".bcf", "_consensus.fa")
    consensus_cmd = f"bcftools consensus -m {mask_file} -f {ref} -I {input_bcf} -o {output_fasta}"
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
    convert_cmd = f"/home/data1/ohad/panthera/alignment/psmc/utils/fq2psmcfa {input_fa} > {output_psmcfa}"
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



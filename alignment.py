import subprocess
import os
import sys

ref_genome = "/home/data1/ohad/panthera/leopard_alignment/reference_snow_leopard.fasta"
working_dir = '/home/data1/ohad/panthera/snow_leopard_alignment/sample_9611'

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



def filter_chromosomes(input_bam, chromosomes):
    """Filter the input BAM file to include only the chromosomes in the list."""
    # Create individual BAM files for each chromosome
    output_bam = f"{input_bam.replace('.bam', '_OnlyChr.bam')}"
    print(f"Filtering chromosomes in {input_bam} to {output_bam}...")
    for chromosome in chromosomes:
        temp = f"temp_{chromosome}.bam"
        filter_cmd = f"samtools view -b {input_bam} {chromosome} -o {temp} -@ 50"
        run_command(filter_cmd)
        print(f"Chromosome {chromosome} was bammed.")
    # Merge the individual BAM files into the final output BAM file
    merge_cmd = f"samtools merge {output_bam} out.Hic_asm_*.bam -@ 50"
    run_command(merge_cmd)
    print("Chromosomes filtered.")



def filter_chromosomes(input_bam, chromosomes):
    """Filter the input BAM file to include only the chromosomes in the list.
    Args:chrosmos: list of chromosomes to keep in txt format"""
    # Create individual BAM files for each chromosome
    output_bam = f"{input_bam.replace('.bam', '_OnlyChr.bam')}"
    chromosomes = ' '.join(chromosomes)
    print(f"Filtering chromosomes in {input_bam} to {output_bam}...")
    filter_cmd = f"samtools view -b {input_bam} {chromosomes} -o {output_bam} -@ 50"
    run_command(filter_cmd)


##### index reference genome
# index_ref_genome(ref_genome)

##### align sequence to reference genome
# index_suffix = "ref_genome_index"
# run_bowtie2_alignment(pattern='*_1.fq.gz.filtered.gz', index_suffix=index_suffix,threads=50)

#### convert sam to bam
# sam_to_bam(working_dir, threads=50)
#
# #### merge bam files
# merge_bam_files(working_dir, 50, 'merged', 'aligned_sorted.bam')
#
##### index merged bam file
# index_sam('9611_merged_OnlyChr.bam',threads=100)

#### filter chromosomes
# chromosoms_list = working_dir + '/chromosomes.txt'
# with open(chromosoms_list, 'r') as f:
#     chromosomes = f.read().splitlines()
# filter_chromosomes('9611_merged.bam', chromosomes)

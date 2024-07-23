import subprocess
import os
import sys

ref_genome = "/home/data1/ohad/panthera/leopard_alignment/reference_snow_leopard.fasta"
working_dir = '/home/data1/ohad/panthera/leopard_alignment'

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


working_dir = '/home/data1/ohad/panthera/leopard_alignment'
def run_bowtie2_alignment(pattern: str, ref_genome: str, index_suffix: str):
    """Align paired-end reads using Bowtie2.
    Args:
        pattern: A glob pattern to match the FASTQ files. e.g. "*_1.fq.gz.filtered.gz"
    """

    # Use the `ls` command to list files matching the pattern and split the output to get a list of file paths
    fastq_files_1_output = run_command(f'ls {pattern}')
    fastq_files_1 = fastq_files_1_output.strip().split('\n')

    # Adjust the pattern to find the second set of files and repeat the process
    pattern2 = pattern.replace("_1", "_2")
    fastq_files_2_output = run_command(f'ls {pattern2}')
    fastq_files_2 = fastq_files_2_output.strip().split('\n')

    # Iterate over each pair of files
    for fq1, fq2 in zip(fastq_files_1, fastq_files_2):
        output = os.path.join(working_dir, f"{os.path.basename(fq1).replace('_1.fq.gz.filtered.gz', '_aligned.sam')}")
        print(output)
        print(f"Aligning {fq1} and {fq2}")

        align_cmd = f"bowtie2  --threads=50 -x={index_suffix} -1={fq1} -2={fq2} -S={output}"
        run_command(align_cmd)


def index_sam(sam_file):
    """Index the BCF file."""
    index_cmd = f"samtools index {sam_file}"
    run_command(index_cmd)
    print(f"SAM file indexed: {sam_file}")

def sam_to_bam(working_dir):
    """Convert SAM files to BAM files. Then sort by coordinate and index the BAM files."""

    sam_files = os.listdir(working_dir)
    print(sam_files)

    for sam_file in sam_files:
        # Define the output BAM file name (replace .sam with .sorted.bam)
        bam_file = f"{sam_file.replace('.sam', '_sorted.bam')}"
        print(f"Converting and sorting {sam_file} to {bam_file}...")
        # Convert SAM to BAM and sort
        convert = f"samtools view -bS {sam_file} | samtools sort -o {bam_file}"
        run_command(convert)

        index_cmd = f"samtools index {bam_file}"
        run_command(index_cmd)
        print(f"BAM file indexed: {bam_file}")



# index_ref_genome(ref_genome)

index_suffix = "ref_genome_index"
# run_bowtie2_alignment(pattern='*_1.fq.gz.filtered.gz', ref_genome=ref_genome,
#                       index_suffix=index_suffix)

##### convert sam to bam
sam_to_bam(working_dir)

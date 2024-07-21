import subprocess
import os


def run_command(command):
    """Run a shell command, capture the output, and print it."""
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
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


ref_genome = "/home/data1/ohad/panthera/leopard_alignment/reference_snow_leopard.fasta"


def index_ref_genome(ref_genome):
    """Index the reference genome for alignment."""
    output_path = os.path.dirname(ref_genome)
    align_cmd = f"bowtie2-build {ref_genome} {output_path}/ref_genome_index"
    run_command(align_cmd)
    print("Reference genome indexed.")


index_ref_genome(ref_genome)


import pysam
import matplotlib.pyplot as plt

working_dir = '/home/data1/ohad/panthera/snow_leopard_alignment/sample_9611'
ref_genome = f"{working_dir}/reference_snow_leopard.fasta"
bcf_path = working_dir + '/variants_9611.bcf'
chromosoms_list = working_dir + '/chromosomes.txt'
with open(chromosoms_list, 'r') as f:
    chromosomes = f.read().splitlines()


def plot_coverage(bcf_path):

    bcf = pysam.VariantFile(bcf_path)
    data = []

    # Extract coverage depth from each record
    for record in bcf:
        if 'DP' in record.info:
            data.append(record.info['DP'])


    # Plotting the histogram
    plt.figure(figsize=(10, 6))
    plt.hist(data, bins=100)
    plt.xlim(0, 100)
    plt.ylabel('Frequency')
    plt.show()

def plot_qual(bcf_path):

    bcf = pysam.VariantFile(bcf_path)
    data = []

    # Extract qual from each record
    for record in bcf:
        data.append(record.qual)

    # Plotting the histogram
    plt.figure(figsize=(10, 6))
    plt.hist(data, bins=100)
    plt.ylabel('Frequency')
    plt.show()


def plot_coverage_by_position(bcf_path,chromosomes):
    """ i want to plot the coverage by position of a bam file to see the quality in
    the telomeres.     # Extract coverage depth from each record for each chromosome
"""
    bcf = pysam.VariantFile(bcf_path)

    for chr in chromosomes:
        coverage = []
        pos = []

        for record in bcf.fetch(chr):
            coverage.append(record.info['DP'])
            pos.append(record.pos)
        # Plotting the histogram
        plt.figure(figsize=(10, 6))
        plt.stem(pos, coverage)
        plt.title(f'Coverage by position for {chr}')
        plt.ylabel('Coverage')
        plt.xlabel('Position')
        # max_pos = max(pos)  # Get the maximum position value
        # min_range = max(0, max_pos - 1000000)  # Subtract 1,000,000 but not less than 0
        # plt.xlim(min_range, max_pos)  # Set the limits to the last 1,000,000 values
        plt.show()
        plt.clf()


##### plot coverage
plot_coverage_by_position(bcf_path,chromosomes)
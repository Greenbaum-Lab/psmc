import pysam
import matplotlib.pyplot as plt

# Path to your VCF file
working_dir = '/home/data1/ohad/panthera/snow_leopard_alignment/sample_12491'
sample_vcf = f"{working_dir}/chr_default.vcf"
# Open VCF file
vcf = pysam.VariantFile(working_dir)

# List to hold coverage depths
depths = []

# Extract coverage depth from each record
for record in vcf:
    if 'DP' in record.info:
        depths.append(record.info['DP'])

# Plotting the histogram
plt.figure(figsize=(10, 6))
plt.hist(depths, bins=50, color='blue', alpha=0.7)
plt.title('Histogram of Coverage Depth')
plt.xlabel('Depth of Coverage')
plt.ylabel('Frequency')
plt.grid(True)
plt.show()

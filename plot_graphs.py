import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.colorbar import ColorbarBase
import os

# Function to parse VCF file and extract SNP data
def parse_vcf(vcf_file):
    snps = []
    with open(vcf_file, 'r') as file:
        for line in file:
            if line.startswith('#'):
                if line.startswith('#CHROM'):
                    header = line.strip().split('\t')
                    sample_names = header[9:]
                    sample_names = [os.path.basename(sample)[0:7] for sample in sample_names]  # Extract filename only
                else:
                    continue
            fields = line.strip().split('\t')
            chrom, pos, id, ref, alt, qual, fil, info, fmt, *genotypes = fields
            genotypes = [gt.split(':')[0] for gt in genotypes]  # Extract genotype data

            snp = {
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "genotypes": genotypes
            }
            snps.append(snp)

    return snps, sample_names


# Function to prepare data for heatmap
def prepare_data(snps, samples):

    data = np.zeros((len(snps), len(samples)))
    
    for i, snp in enumerate(snps):
        for j, genotype in enumerate(snp["genotypes"]):
            if genotype == '0/0':
                data[i, j] = 0  # Homozygous reference
            elif genotype == '0/1' or genotype == '1/0':
                data[i, j] = 1  # Heterozygous
            elif genotype == '1/1':
                data[i, j] = 2  # Homozygous alternate
            else:
                data[i, j] = -1  # Missing data
    return data

# Function to generate heatmap
def generate_heatmap(data, samples):
    
    colors = ["white", "blue", "yellow", "red"]
    cmap = ListedColormap(colors)
    plt.figure(figsize=(10, 8))
    plt.imshow(data, cmap=cmap, aspect='auto', interpolation='nearest')
    # Add colorbar with discrete values
    legend_labels = ['Missing data', 'Homozygous Ref', 'Heterozygous', 'Homozygous Alt']
    legend_handles = [plt.Rectangle((0, 0), 1, 1, color=color) for color in colors]
    plt.legend(legend_handles, legend_labels, loc='upper right')
    plt.title('SNP Heatmap')
    plt.xlabel('Samples')
    plt.xticks(np.arange(len(samples)), samples, rotation=90)
    plt.ylabel('SNP Position')
    # Add yticks with SNP positions instead of indices once every 10 SNPs
    plt.yticks(np.arange(0, data.shape[0], 30), [snp['pos'] for snp in snps[::30]])
    plt.tight_layout()
    plt.savefig('heatmap.png', dpi=300, bbox_inches='tight')  # Save the figure with higher resolution and tight bounding box
    #plt.show()


if __name__ == "__main__":
    # Example usage
    vcf_file = 'genes.vcf'
    snps, samples = parse_vcf(vcf_file)
    data = prepare_data(snps, samples)
    generate_heatmap(data, samples)




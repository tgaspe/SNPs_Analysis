import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.colorbar import ColorbarBase

# Function to parse VCF file and extract SNP data
def parse_vcf(vcf_file):
    snps = []
    with open(vcf_file, 'r') as file:
        for line in file:
            if line.startswith('#'):
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
    return snps


# Function to prepare data for heatmap
def prepare_data(snps):

    samples = [ "sample" + str(i)  for i in range(len(snps[0]["genotypes"]))] 
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
    
    return data, samples

# Function to generate heatmap
def generate_heatmap(data, samples):

    colors = ["white", "blue", "yellow", "red"]
    cmap = ListedColormap(colors)
    plt.figure(figsize=(10, 8))
    plt.imshow(data, cmap=cmap, aspect='auto', interpolation='nearest')
    
    # Add colorbar with discrete values
    colorbar = plt.colorbar(orientation='vertical', ticks=[-1, 0, 1, 2], format='%1i', spacing='uniform')
    colorbar.set_ticklabels(['Missing data', 'Homozygous Ref', 'Heterozygous', 'Homozygous Alt'])
    colorbar.set_label('Genotype Type')

    plt.title('SNP Heatmap')
    plt.xlabel('Samples')
    plt.ylabel('SNP Position')
    plt.xticks(np.arange(len(samples)), samples, rotation=90)
    plt.tight_layout()
    plt.savefig('heatmap.png')


if __name__ == "__main__":
    # Example usage
    vcf_file = 'snp_call_snakemake/genes.vcf'
    snps = parse_vcf(vcf_file)
    data, samples = prepare_data(snps)
    #print(data)
    #print(samples)
    generate_heatmap(data, samples)




#import pandas as pd
import numpy as np


"""
TODO: 
- parse VCF file
- create DATA matrix
- generate Heatmap

"""

def parse_vcf(filename):
    # Create a SNP list
    snps_list = []
    f = open(filename)
    for line in f:
        if "#" in line:
            continue
        else:
            # Split line into fields
            data = line.split("\t")
            snp = { # Create dictionary
                "position": data[1],
                "ref_allele": data[3],
                "alt_allele": data[4],
                "quality": data[5],
                "info": data[7],
                "format": data[8],
                "genotype": data[9:],
            }
            snps_list.append(snp)
    f.close()

    return snps_list
    


# Step 3: Create a Data Matrix
def create_data_matrix(filtered_snps, samples):
    # Initialize an empty DataFrame to store the data matrix
    data_matrix = pd.DataFrame(index=samples)

    # Iterate over filtered SNPs and populate the data matrix
    for snp in filtered_snps:
        # Extract SNP information
        snp_name = snp['name']
        snp_genotypes = snp['genotypes']

        # Create a column in the data matrix for the SNP
        data_matrix[snp_name] = [snp_genotypes.get(sample, None) for sample in samples]

    return data_matrix


if __name__ == "__main__":
    # Example usage:
    snps_list = parse_vcf("genes.vcf")

    #data_matrix = create_data_matrix(filtered_snps, samples)
    print(snps_list[0])
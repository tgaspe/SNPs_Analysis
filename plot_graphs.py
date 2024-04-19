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
            else:
                fields = line.strip().split('\t')
                chrom, pos, id, ref, alt, qual, fil, info, fmt, *genotypes = fields
                genotypes = [gt.split(':')[0] for gt in genotypes]  # Extract genotype data
                
                _, func, impact, gene, *_ = info.split("|")  # Extract fields from the INFO column and ignore the rest of the list items

                snp = {
                    "pos": int(pos),
                    "ref": ref,
                    "alt": alt,
                    "qual": qual,
                    "func": func, # functional class is the 2nd field in the INFO column
                    "impact": impact, # impact class is the 3rd field in the INFO column
                    "gene": gene,  # gene of interest is the 4th field in the INFO column
                    "genotypes": genotypes
                }
                snps.append(snp)

    return snps, sample_names

# Function to count the number of SNPs in each functional class
def count_functional_classes(snps):
    functional_classes = {
        'synonymous_variant': 0,
        'missense_variant': 0,
        'stop_gained': 0,
        'frameshift_variant': 0,
        'inframe_insertion': 0,
        'inframe_deletion': 0,
        'splice_acceptor_variant': 0,
        'splice_donor_variant': 0,
        'start_lost': 0,
        'stop_lost': 0,
        'stop_retained_variant': 0,
        'start_retained_variant': 0
    }
    for snp in snps:
        if snp["func"] in functional_classes:
            functional_classes[snp["func"]] += 1
        else:
            functional_classes[snp["func"]] = 1
    return filter_dictionary(functional_classes)

# Helper method: Filter dict keys with zero values
def filter_dictionary(dictionary):
    filtered_dict = {k: v for k, v in dictionary.items() if v > 0}
    return filtered_dict

# Function to count the number of SNPs in each impact class
def count_impact_classes(snps):
    impact_classes = {
        'MODIFIER': 0,
        'LOW': 0,
        'MODERATE': 0,
        'HIGH': 0
    }
    for snp in snps:
        if snp["impact"] in impact_classes:
            impact_classes[snp["impact"]] += 1
        else:
            impact_classes[snp["impact"]] = 1
    return impact_classes

# Bar chart for impact classes
def generate_bar_chart(impact_classes, title, xlabel, ylabel, filename):
    plt.figure(figsize=(8, 6))
    plt.bar(impact_classes.keys(), impact_classes.values())
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(filename, dpi=300)
    #plt.show()

# Function count the number of SNPs per sample and gene (APP, SOD1, DYRK1A)
# def heatmap_matrix(snps, gene_list, samples):
#     data = np.zeros((3, len(samples)))
#     for snp in snps:
#         #for gene in gene_list:
#         if snp["gene"] == "APP":
#             for j, genotype in enumerate(snp["genotypes"]):
#                 # Homozygous alternative
#                 if genotype == '1/1':
#                     data[0, j] += 1  
#                 # Heterozygous
#                 elif genotype == '0/1' or genotype == '1/0':
#                     data[0, j] += 1  
#                 else: # Missing data or homozygous reference
#                     continue 
#         elif snp["gene"] == "SOD1":
#             for j, genotype in enumerate(snp["genotypes"]):
#                 # Homozygous alternative
#                 if genotype == '1/1':
#                     data[1, j] += 1  
#                 # Heterozygous
#                 elif genotype == '0/1' or genotype == '1/0':
#                     data[1, j] += 1  
#                 else:
#                     continue
#         elif snp["gene"] == "DYRK1A":
#             for j, genotype in enumerate(snp["genotypes"]):
#                 # Homozygous alternative
#                 if genotype == '1/1':
#                     data[2, j] += 1  
#                 # Heterozygous
#                 elif genotype == '0/1' or genotype == '1/0':
#                     data[2, j] += 1  
#                 else:
#                     continue
#     return data

# Function count the number of SNPs per sample and gene (APP, SOD1, DYRK1A)
def heatmap_matrix(snps, gene_list, samples_list):
    # Generating a matrix filed with zeros. rows = genes, columns = samples
    data = np.zeros((len(gene_list), len(samples_list)))
    for snp in snps:
        for gene in gene_list:
            if snp["gene"] == gene:
                for j, genotype in enumerate(snp["genotypes"]):
                    # Homozygous alternative
                    if genotype == '1/1':
                        data[0, j] += 1  
                    # Heterozygous
                    elif genotype == '0/1' or genotype == '1/0':
                        data[0, j] += 1  
                    else: # Missing data or homozygous reference
                        continue 
    return data

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
                data[i, j] = -1 # Missing data
    return data

# General function to generate heatmap
# def generate_heatmap(dictionary, title, xlabel, ylabel, filename):
#     plt.figure(figsize=(10, 8))
#     plt.imshow(dictionary, cmap="coolwarm", aspect='auto', interpolation='nearest')
#     plt.title(title)
#     plt.xlabel(xlabel)
#     plt.xticks(np.arange(len(dictionary.keys())), dictionary.keys(), rotation=90)
#     plt.ylabel(ylabel)
#     for i in range(dictionary.keys()):
#         for j in range(data.shape[1]):
#             plt.text(j, i, data[i, j], ha='center', va='center', color='black')
#     plt.tight_layout()
#     plt.savefig('heatmap1.png', dpi=300, bbox_inches='tight')

# Function to generate heatmap
def generate_heatmap(data, samples, filename='heatmap.png'):
    
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
    plt.savefig(filename, dpi=300, bbox_inches='tight')  # Save the figure with higher resolution and tight bounding box
    #plt.show()

def generate_heatmap1(data, samples, filename='heatmap.png'):
    plt.figure(figsize=(10, 8))
    plt.imshow(data, cmap="coolwarm", aspect='auto', interpolation='nearest')
    # Add colorbar with discrete values
    
    plt.title('SNP Heatmap')
    # Add data values in the cells of the heatmap
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            plt.text(j, i, data[i, j], ha='center', va='center', color='black')
    plt.xlabel('Samples')
    plt.xticks(np.arange(len(samples)), samples, rotation=90)
    plt.ylabel('Genes')
    plt.yticks(np.arange(len(data)), range(len(data)))
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    #plt.show()

# Bar chart 
def generate_bar_chart(dictionary, title, xlabel, ylabel, filename, rotation=90):
    plt.figure(figsize=(8, 6))
    plt.bar(dictionary.keys(), dictionary.values())
    plt.xticks(rotation=rotation)  # Rotate x-axis labels by 90 degrees
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(filename, dpi=300)
    #plt.show()

# Pie chart
def generate_pie_chart(dictionary, title, filename):
    plt.figure(figsize=(8, 8))
    plt.pie(dictionary.values(), labels=dictionary.keys(), autopct='%1.1f%%', startangle=140)
    plt.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    plt.title(title)
    plt.savefig(filename, dpi=300)
    #plt.show()



if __name__ == "__main__":
    
    # Load data
    vcf_file = 'snp_call_snakemake/genes.vcf'
    snps, sample_names = parse_vcf(vcf_file)

    # Create dictionaries for functional and impact classes
    func_dict = count_functional_classes(snps)
    impact_dict = count_impact_classes(snps)
    
    # Create a heatmap matrix for the genes APP, SOD1, DYRK1A
    matrix = heatmap_matrix(snps, ['APP', 'SOD1', 'DYRK1A'], sample_names)
    data = prepare_data(snps, sample_names)
    
    # Generate heatmap
    generate_heatmap(data, sample_names, "heatmap1.png")
    generate_heatmap1(matrix, sample_names, "heatmap2.png")
    
    # Generate pie charts
    generate_pie_chart(impact_dict, 'Impact Classes Distribution', 'pie_impact_classes.png')
    generate_pie_chart(func_dict, 'Functional Classes Distribution', 'pie_functional_classes.png')
    
    # Generate bar charts
    generate_bar_chart(impact_dict, 'Impact Classes Distribution', 'Impact Class', 'Number of SNPs', 'impact_classes.png')
    generate_bar_chart(func_dict, 'Functional Classes Distribution', 'Functional Class', 'Number of SNPs', 'functional_classes.png')

    
    # data = prepare_data(snps, samples)
    # data1 = prepare_data1(snps, samples)

    # generate_pie_chart(count_impact_classes(snps), 'Impact Classes Distribution', 'pie_impact_classes.png')
    # generate_pie_chart(filtered_func_dict, 'Functional Classes Distribution', 'pie_functional_classes.png')
    # generate_bar_chart(count_impact_classes(snps), 'Impact Classes Distribution', 'Impact Class', 'Number of SNPs', 'impact_classes.png')
    # generate_bar_chart(count_functional_classes(snps), 'Functional Classes Distribution', 'Functional Class', 'Number of SNPs', 'functional_classes.png')
    # generate_bar_chart(filtered_func_dict, 'Functional Classes Distribution', 'Functional Class', 'Number of SNPs', '1111functional_classes.png')
    # print(count_functional_classes(snps))

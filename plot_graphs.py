import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.colorbar import ColorbarBase
import os
import sys

# Function to parse VCF file and extract SNP data
def parse_vcf(vcf_file):
    snps = []
    sample_names = []
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


# Function count the number of SNPs per sample and gene (APP, SOD1, DYRK1A)
def heatmap_matrix(snps, gene_list, samples_list):
    # Generating a matrix filed with zeros. rows = genes, columns = samples
    data = np.zeros((len(gene_list), len(samples_list)))
    for snp in snps:
        for i, gene in enumerate(gene_list):
            if snp["gene"] == gene:
                for j, genotype in enumerate(snp["genotypes"]):
                    # Homozygous alternative
                    if genotype == '1/1':
                        data[i, j] += 1  
                    # Heterozygous
                    elif genotype == '0/1' or genotype == '1/0':
                        data[i, j] += 1  
                    else: # Missing data or homozygous reference
                        continue 
    return data


# Heatmap
def generate_heatmap1(data, samples, filename='heatmap.png'):
    plt.rcParams['font.family'] = 'monospace'
    plt.figure(figsize=(12, 9))
    plt.imshow(data, cmap="coolwarm", aspect='auto', interpolation='nearest')
      # Change the font family to Arial
    plt.title('SNP Heatmap', fontweight='bold', fontsize=24, pad=20)
    # Add data values in the cells of the heatmap
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            plt.text(j, i, data[i, j], ha='center', va='center', color='black')
    plt.xlabel('Samples', fontweight='bold', fontsize=14, labelpad=15)
    plt.xticks(np.arange(len(samples)), samples, rotation=90)
    plt.ylabel('Genes', fontweight='bold', fontsize=14, labelpad=15)
    plt.yticks(np.arange(len(data)), ['APP', 'SOD1', 'DYRK1A'])
    # add a colorbar
    cbar = plt.colorbar()
    cbar.set_label('Number of SNPs', fontweight='bold', fontsize=14, labelpad=15)

    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    #plt.show()


# Bar chart 
def generate_bar_chart(dictionary, title, xlabel, ylabel, filename, rotation=90):
    # Change the font family
    plt.rcParams['font.family'] = 'monospace'
    # Truncate long labels
    truncated_labels = [key[:10] for key in dictionary.keys()] 
    plt.figure(figsize=(8, 6))
    plt.bar(truncated_labels, dictionary.values(), color='green')
    plt.xticks(rotation=rotation)  # Rotate x-axis labels by 90 degrees
    plt.title(title, fontweight='bold', fontsize=20)
    plt.xlabel(xlabel, fontweight='bold', fontsize=12, labelpad=5)
    plt.ylabel(ylabel, fontweight='bold', fontsize=12, labelpad=5)
    plt.savefig(filename, dpi=300)
    #plt.show()


# Pie chart
def generate_pie_chart(dictionary, title, filename):
    # Change the font family
    plt.rcParams['font.family'] = 'monospace'
    plt.figure(figsize=(8, 8))
    patches, texts, autotexts = plt.pie(dictionary.values(), autopct=my_autopct, startangle=140)
    plt.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    plt.legend(patches, dictionary.keys(), loc='lower right', fontsize='small')
    plt.title(title, fontweight='bold', fontsize=20)
    for autotext in autotexts:
        autotext.set_horizontalalignment('center')
        autotext.set_verticalalignment('center')
    plt.tight_layout()  # Adjust layout to prevent overlapping
    plt.savefig(filename, dpi=300)
    #plt.show()


# Custom autopct function to display percentages only if they are greater than 1.5%
def my_autopct(pct):
    return '{:.2f}%'.format(pct) if pct >= 1.5 else ''


if __name__ == "__main__":
    
    # Usage Message
    if len(sys.argv) != 3:
        print("Usage: python3 plot_graphs.py <argument> <vcf_file>")
        exit(1)
    
    # Parse VCF file
    vcf_file = sys.argv[2]
    snps, sample_names = parse_vcf(vcf_file)

    if (sys.argv[1] == 'impact'):
        # Generate pie and bar chart for impact classes
        generate_pie_chart(count_impact_classes(snps), 'SNPs Impact Distribution', 'images/impact_pie_chart.png')
        generate_bar_chart(count_impact_classes(snps), 'SNPs Impact Distribution', 'Impact Class', 'Number of SNPs', 'images/impact_bar_chart.png', 0)
    
    elif (sys.argv[1] == 'functional'):
        # Generate pie and bar chart for functional classes
        snps_classes = count_functional_classes(snps)
        # Summing up percentages below the threshold
        total = sum(snps_classes.values())
        small_total = sum(val for key, val in snps_classes.items() if val/total < 0.015)
        # Creating a modified dictionary with 'Others'
        modified_dictionary = {key: val for key, val in snps_classes.items() if val/total >= 0.015}
        modified_dictionary['Others'] = small_total
        # Generate pie and bar chart
        generate_pie_chart(modified_dictionary, 'SNPs Classes Distribution', 'images/snp_classes_pie_chart.png')
        generate_bar_chart(modified_dictionary, 'SNPs Classes Distribution', 'Functional Class', 'Number of SNPs', 'images/snp_classes_bar_chart.png', 0)
    
    elif (sys.argv[1] == 'heatmap'):
        # Create a heatmap matrix for the genes APP, SOD1, DYRK1A
        matrix = heatmap_matrix(snps, ['APP', 'SOD1', 'DYRK1A'], sample_names)
        generate_heatmap1(matrix, sample_names, "heatmap.png")
    
    else:
        print("Invalid argument. Please enter 'impact', 'functional' or 'heatmap' as an argument")
        exit(1)
    
    print("Done!")
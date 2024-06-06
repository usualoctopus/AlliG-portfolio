# This script is a scalable, dynamic workflow for preparing gene lists from edgeR output for use in Rank-rank hypergeometric overlap (RRHO) analysis.

### Multiple files can be provided at the same time and labeled dynamically (in my case cell line and hormones are different for each gene list) using a metadata csv file as input. Each gene list is separated and sorted into positive and negative logFC in order to assign directionality to each log10(p-value): positive values correspond to negative associations and vice versa. The gene list is then recombined and saved to an output directory as a two column file: Genes and log10(p-value) with assigned directionality.  For ease of use downstream, all lists are then combined into one .csv file, sorted by Gene name.

#### Data input for RRHO can be transformed in several ways, but log10(p-value) with directionality assignment from FC is a common form and originially presented here: Seema B. Plaisier, Richard Taschereau, Justin A. Wong, Thomas G. Graeber, Rank–rank hypergeometric overlap: identification of statistically significant overlap between gene-expression signatures, Nucleic Acids Research, Volume 38, Issue 17, 1 September 2010, Page e169, https://doi.org/10.1093/nar/gkq636)


```python
import pandas as pd
import numpy as np
import os


# Create a function to read in a single gene list, add -log10(p-value) column, 
# separated and sort lists into neg and pos logFC, assign directionality to logPValue, recombine to a single gene list
def process_gene_list(file_path, cell_line, hormone):
    # Read in the data
    gene_data = pd.read_csv(file_path)

    # Rename the first column to "Gene" and add a column with -log10(p-value)
    gene_data = gene_data.rename(columns={gene_data.columns[0]: "Gene"})
    gene_data['logPValue'] = -gene_data['PValue'].apply(lambda p: np.log10(p))

    # Extract only genes negatively correlated with the hormone
    gene_data_neg = gene_data[gene_data['logFC'] < 0].copy()
    gene_data_neg['logPValue'] = -gene_data_neg['logPValue']
    gene_data_neg = gene_data_neg.sort_values(by='logPValue', ascending=False)

    # Extract only genes positively correlated with the hormone
    gene_data_pos = gene_data[gene_data['logFC'] > 0].copy()
    gene_data_pos = gene_data_pos.sort_values(by='logPValue', ascending=False)

    # Combine the sorted positive and negative lists
    gene_data_sorted = pd.concat([gene_data_pos, gene_data_neg])

    # Select the Gene and logPValue columns and rename the logPValue column
    logPValue_colname = f"{cell_line}.{hormone}.logPValue"
    gene_data_sorted_selected = gene_data_sorted[['Gene', 'logPValue']]
    gene_data_sorted_selected = gene_data_sorted_selected.rename(columns={'logPValue': logPValue_colname})

    return gene_data_sorted_selected


# Create a function to process multiple gene lists using metadata and the function above, save results
def process_multiple_gene_lists(metadata_file, output_dir):
    # Read the metadata file
    metadata = pd.read_csv(metadata_file)

    # Ensure the output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    processed_files = []

    for i, row in metadata.iterrows():
        file_path = row['file_path']
        cell_line = row['cell_line']
        hormone = row['hormone']

        # Process the gene list
        processed_data = process_gene_list(file_path, cell_line, hormone)

        # Extract the base file name without extension
        base_file_name = os.path.splitext(os.path.basename(file_path))[0]

        # Construct the new output file name
        output_file_name = f"{cell_line}_{hormone}_DGE_sorted_log10pRank.csv"
        output_file = os.path.join(output_dir, output_file_name)

        # Save the processed data
        processed_data.to_csv(output_file, index=False)

        # Add the processed data to the list
        processed_files.append(processed_data)

    return processed_files


# Function to combine multiple datasets and sort by Gene name
def combine_gene_lists(processed_files):
    # Start with the first data frame
    combined_data = processed_files[0]

    # Iteratively perform inner joins with the remaining data frames
    for i in range(1, len(processed_files)):
        combined_data = combined_data.merge(processed_files[i], on='Gene')

    # Remove duplicate rows
    combined_data = combined_data.drop_duplicates()

    # Sort the combined data by Gene name
    combined_data = combined_data.sort_values(by='Gene')

    return combined_data


# Example usage
metadata_file = "~/Documents/Documents_AlliMacBookPro/Neuronal_cell_lines/Analyses/RRHO/Cell_line_metadata_RRHO.csv"
output_dir = "~/Documents/Documents_AlliMacBookPro/Neuronal_cell_lines/Analyses/RRHO/DL_UL_SH_NPCcontrol/Processed_genelists"

# Process multiple gene lists
processed_files = process_multiple_gene_lists(metadata_file, output_dir)

# Combine the processed gene lists and sort by Gene name
combined_data = combine_gene_lists(processed_files)

# Save the combined and sorted data
combined_data.to_csv(os.path.join(output_dir, "DLULSHNSCcontrol_AllHormone_combined_gene_lists.csv"), index=False)

```

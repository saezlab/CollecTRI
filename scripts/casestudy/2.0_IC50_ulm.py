import pandas as pd
import decoupler as dc

CollecTRI = pd.read_csv("CollecTRI.csv")
Dorothea_ABC = pd.read_csv("dorothea_ABC.csv")

# Load GDSC2 data to filter for cell lines of interest
GDSC2_data = pd.read_csv("GDSC2_fitted_dose_response.csv", delimiter="\t")

# Read gene counts data from CSV file
gene_counts = pd.read_csv("rnaseq_fpkm_20220624_flt.csv")

# Remove rows with missing values
gene_counts = gene_counts.dropna()

# Convert columns 2 to ncol(gene_counts) to numeric
gene_counts.iloc[:, 2:] = gene_counts.iloc[:, 2:].apply(pd.to_numeric)


# Get unique cell line IDs from GDSC2 data
cell_lines = GDSC2_data['SANGER_MODEL_ID'].unique()
common_cell_lines = set(cell_lines) & set(gene_counts.columns)

# Filter gene counts for selected cell lines
GDSC2_expr = gene_counts[['X'] + list(common_cell_lines)].copy()
GDSC2_expr = GDSC2_expr.drop_duplicates(subset='X', keep='first')
GDSC2_expr = GDSC2_expr.set_index('X')
GDSC2_expr = GDSC2_expr.T

acts_CollecTRI, pvals_CollecTRI = dc.run_ulm(GDSC2_expr, CollecTRI, min_n=5)

acts_CollecTRI = acts_CollecTRI.reset_index()
pvals_CollecTRI = pvals_CollecTRI.reset_index()

# Melt the 'acts' DataFrame to long format
acts_CollecTRI_long = pd.melt(acts_CollecTRI, id_vars=['index'], var_name='TF', value_name='Activity')
acts_CollecTRI_long.columns = ['Sample', 'TF', 'Activity']

# Melt the 'pvals' DataFrame to long format
pvals_CollecTRI_long = pd.melt(pvals_CollecTRI, id_vars=['index'], var_name='TF', value_name='Pval')
pvals_CollecTRI_long.columns = ['Sample', 'TF', 'Pval']

# Merge the two long DataFrames based on 'Sample' and 'TF' columns
IC50_ulm_CollecTRI = pd.merge(acts_CollecTRI_long, pvals_CollecTRI_long, on=['Sample', 'TF'])

IC50_ulm_CollecTRI.to_csv('IC50_ulm_CollecTRI.csv', index=False)

acts_Dorothea, pvals_Dorothea = dc.run_ulm(GDSC2_expr, Dorothea_ABC, min_n=5)

acts_Dorothea = acts_Dorothea.reset_index()
pvals_Dorothea = pvals_Dorothea.reset_index()

# Melt the 'acts' DataFrame to long format
acts_Dorothea_long = pd.melt(acts_Dorothea, id_vars=['index'], var_name='TF', value_name='Activity')
acts_Dorothea_long.columns = ['Sample', 'TF', 'Activity']

# Melt the 'pvals' DataFrame to long format
pvals_Dorothea_long = pd.melt(pvals_Dorothea, id_vars=['index'], var_name='TF', value_name='Pval')
pvals_Dorothea_long.columns = ['Sample', 'TF', 'Pval']

# Merge the two long DataFrames based on 'Sample' and 'TF' columns
IC50_ulm_Dorothea = pd.merge(acts_Dorothea_long, pvals_Dorothea_long, on=['Sample', 'TF'])

IC50_ulm_Dorothea.to_csv('IC50_ulm_Dorothea.csv', index=False)
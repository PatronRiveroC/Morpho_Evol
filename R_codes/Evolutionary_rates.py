
### Title: Rates of trait evolution ####
### Author: Patrón-Rivero, C. ####
### Date: 25/02/2026 ###
### Project: "Integration and modular patterns shaping the morphological evolution in Neotropical hognose pitvipers Porthidium" ###

# ------------------------------------------------------------------------------------------------ #

# Libraries #

# ------------------------------------------------------------------------------------------------ #

import os
import numpy as np
import pandas as pd
import dendropy
import matplotlib.pyplot as plt
from tqdm import tqdm

# ------------------------------------------------------------------------------------------------ #

# Settings #

# ------------------------------------------------------------------------------------------------ #
 
BASE_PATH = r"D:\1_morpho_evol"
CSV_FILE = "adu_traits_mean_pca_spp.csv"
TREE_FILE = "P_tree.newick"

TOTAL_GENERATIONS = 20_000_000 
THINNING = 2000                  
BURN_IN_PERCENT = 0.25           

# ------------------------------------------------------------------------------------------------ #

# Auxiliary Functions #

# ------------------------------------------------------------------------------------------------ #
 
def calculate_ess_manual(series):
    series = np.array(series)
    n = len(series)
    if n < 2: return 0
    
    mean = np.mean(series)
    var = np.var(series)
    
    if var == 0: return 0
    
    centered = series - mean
    corr = np.correlate(centered, centered, mode='full')
    corr = corr[n-1:] 
    acf = corr / (var * n)
    
    sum_rho = 0
    for rho in acf[1:]: 
        if rho < 0:
            break
        sum_rho += rho
        
    ess = n / (1 + 2 * sum_rho)
    return ess

def log_likelihood_bm(traits_data, sigma_sq, vcv_inv, vcv_det, n_taxa):
    if sigma_sq <= 0: return -np.inf
    log_det = n_taxa * np.log(sigma_sq) + np.log(vcv_det)
    centered_data = traits_data - np.mean(traits_data) 
    term_quad = (centered_data.T @ (vcv_inv / sigma_sq) @ centered_data)
    ll = -0.5 * (log_det + term_quad + n_taxa * np.log(2 * np.pi))
    return ll

# ------------------------------------------------------------------------------------------------ #

# Data Loading #

# ------------------------------------------------------------------------------------------------ #

print("--- Loading data and phylogenetic tree ---")

df = pd.read_csv(os.path.join(BASE_PATH, CSV_FILE))
df.set_index('Species', inplace=True)

tree = dendropy.Tree.get(
    path=os.path.join(BASE_PATH, TREE_FILE),
    schema="newick",
    preserve_underscores=True
)

tree_taxa = [t.label for t in tree.taxon_namespace]
common_species = list(set(tree_taxa) & set(df.index))
tree.retain_taxa_with_labels(common_species)
df = df.loc[tree_taxa]

print(f"Data aligned for {len(df)} species.")
traits_columns = df.columns.tolist()

# ------------------------------------------------------------------------------------------------ #

# VCV Matrix Preparation #

# ------------------------------------------------------------------------------------------------ #

print("Calculating phylogenetic matrix...")
pdm = tree.phylogenetic_distance_matrix()
vcv_matrix = np.zeros((len(df), len(df)))

for i, sp1 in enumerate(df.index):
    for j, sp2 in enumerate(df.index):
        taxon1 = tree.taxon_namespace.get_taxon(label=sp1)
        taxon2 = tree.taxon_namespace.get_taxon(label=sp2)
        
        mrca = pdm.mrca(taxon1, taxon2)
        dist = mrca.distance_from_root()
        
        if dist is None:
            dist = 0.0
            
        vcv_matrix[i, j] = dist

vcv_inv = np.linalg.inv(vcv_matrix)
vcv_det = np.linalg.det(vcv_matrix)
print("VCV matrix calculated successfully.")

# ------------------------------------------------------------------------------------------------ #

# MCMC Execution #

# ------------------------------------------------------------------------------------------------ #

num_traits = len(traits_columns)
samples_stored = int(TOTAL_GENERATIONS / THINNING)
chain_sigma = np.zeros((samples_stored, num_traits))
current_sigmas = np.ones(num_traits) * 0.5
current_lls = np.zeros(num_traits)
data_matrix = df.values 
n_taxa = len(df)

for t in range(num_traits):
    current_lls[t] = log_likelihood_bm(data_matrix[:, t], current_sigmas[t], vcv_inv, vcv_det, n_taxa)

print(f"\nStarting MCMC ({TOTAL_GENERATIONS} generations)...")

for step in tqdm(range(TOTAL_GENERATIONS)):
    proposal_sigmas = current_sigmas + np.random.normal(0, 0.05, num_traits)
    proposal_sigmas = np.abs(proposal_sigmas)
    
    for t in range(num_traits):
        prop_ll = log_likelihood_bm(data_matrix[:, t], proposal_sigmas[t], vcv_inv, vcv_det, n_taxa)
        if prop_ll - current_lls[t] > np.log(np.random.rand()):
            current_sigmas[t] = proposal_sigmas[t]
            current_lls[t] = prop_ll
            
    if step % THINNING == 0:
        idx = int(step / THINNING)
        if idx < samples_stored:
            chain_sigma[idx, :] = current_sigmas

# ------------------------------------------------------------------------------------------------ #

# Post-Sampling Analysis #

# ------------------------------------------------------------------------------------------------ #

print("\n--- Processing Results ---")

df_results = pd.DataFrame(chain_sigma, columns=traits_columns)

burn_in_idx = int(samples_stored * BURN_IN_PERCENT)
df_posterior = df_results.iloc[burn_in_idx:].copy()

ess_values = []
for col in traits_columns:
    ess = calculate_ess_manual(df_posterior[col].values)
    ess_values.append(ess)

summary_table = pd.DataFrame({
    'Mean_Sigma_sq': df_posterior.mean(),
    'SD': df_posterior.std(),
    'HPD_95_low': df_posterior.quantile(0.025),
    'HPD_95_high': df_posterior.quantile(0.975),
    'ESS': ess_values
})

print("\n=== RESULTS SUMMARY TABLE ===")
print(summary_table)

output_file = os.path.join(BASE_PATH, "mcmc_results_summary.csv")
summary_table.to_csv(output_file)
print(f"\nResults saved to: {output_file}")


# ------------------------------------------------------------------------------------------------ #

# Visualization #

# ------------------------------------------------------------------------------------------------ #

fig, axes = plt.subplots(num_traits, 1, figsize=(10, 15), sharex=True)
if num_traits == 1: axes = [axes]

for i, col in enumerate(traits_columns):
    ax = axes[i]
    ax.plot(df_results[col], color='gray', alpha=0.3, label='Burn-in')
    ax.plot(np.arange(burn_in_idx, samples_stored), df_posterior[col], label='Posterior', color='blue')
    ax.set_ylabel(col)
    ax.axvline(x=burn_in_idx, color='red', linestyle='--')
    ess_val = summary_table.loc[col, 'ESS']
    ax.text(0.02, 0.9, f"ESS: {ess_val:.1f}", transform=ax.transAxes, bbox=dict(facecolor='white', alpha=0.8))

axes[0].legend()
plt.tight_layout()
plt.savefig(os.path.join(BASE_PATH, "trace_plots.png"))
print("Trace plots generated.")
plt.show()

# ------------------------------------------------------------------------------------------------ #

### End Not Run


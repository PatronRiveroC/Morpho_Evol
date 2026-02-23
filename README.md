# Morphological evolution of size, scales and shape in *Porthidium* neotropical hognose pitvipers

**Authors:**  
Carlos Patron-Rivero*, Carlos Yañez-Arenas*, Xavier Chiappa-Carrara, Octavio Rojas-Soto & Sara Ruane  
\**Corresponding authors

---

## 📖 Abstract
The evolutionary mechanisms driving morphological diversification remain poorly understood in many lineages. We investigated the morphological evolution of the hognose pitviper genus *Porthidium* (Viperidae: Crotalinae). Using a comprehensive dataset of size, scales, and shape from 489 specimens, we tested patterns of modularity versus integration through standardized major axes correlations and evolutionary rates (σ²; estimated with Monte Carlo Markov Chains) and assessed trait evolutionary trajectories with three models: Brownian motion (BM), Ornstein-Uhlenbeck (OU), and early burst (EB). We tested if morphological traits were constrained by their phylogenetic relationships while accounting for ontogenetic shifts between juveniles and adults.  

Our analyses revealed a dichotomy: while direct trait correlations indicated a current integration between size and scale, our analysis of evolutionary rates supported modularity. The BM model was consistently supported across all traits, indicating a pattern resembling neutral evolution, but we also found an optimum directed evolution fit for adult snakes. Finally, we found no significant phylogenetic signal for any trait or stage.  

Our findings suggest that variable selective pressures on complex traits likely drive rapid phenotypic differentiation. This highlights that the forces shaping morphological diversity in these pitvipers are dynamic and context-dependent, emphasizing the need to consider complex evolutionary trajectories when studying morphological diversification.  

---

## 🔑 Keywords
- Morphological Evolution  
- Modularity  
- Integration  
- Evolutionary Rates  
- Evolutionary Models  
- Phylogenetic Signal  
- Morphology  
- Evolution  
- Reptiles  
- Viperidae  

---

## 🛠️ Methods Overview
- Morphological dataset of 489 specimens (*Porthidium*).  
- Standardized Major Axis (SMA) correlations to test integration.  
- Evolutionary rates (σ²) estimated via Monte Carlo Markov Chains.  
- Model fitting: Brownian Motion (BM), Ornstein-Uhlenbeck (OU), Early Burst (EB).  
- Phylogenetic signal tested across traits and ontogenetic stages.  

---

## 📂 Data & Scripts

- adu_dnc.csv # Distance to niche centroid (adults)
- adu_NVH.csv # Allometric and scale traits dataset (adults)
- adu_pc.csv # PCA scores for adult specimens (PC1, PC2, PC3…)
- adu_pc2.csv # Secondary PCA projection for adults (subset or validation)
- adu_traits_mean_pca_spp.csv # Species-level mean PCA scores (adults)
- juv_dnc.csv # Distance to niche centroid (juveniles)
- juv_NVH.csv # Allometric and scale traits dataset (juveniles)
- juv_pc.csv # PCA scores for juvenile specimens
- juv_pc2.csv # Secondary PCA projection for juveniles
- juv_traits_mean_pca_spp.csv # Species-level mean PCA scores (juveniles)
- stat_test.csv # Results of statistical tests (correlations, significance)

- 1_sma_analysis.r # Standardized Major Axis (SMA) analysis for trait integration
- 2_rates_MCMC.r # Estimation of evolutionary rates (σ²) with MCMC
- 3_evol_ps_models.r # Evolutionary models fitting (BM, OU, EB) + phylogenetic signal

---

## 📌 Citation
If you use this repository, please cite:  

**Patron-Rivero C., Yañez-Arenas C., Chiappa-Carrara X., Rojas-Soto O., & Ruane S. (2025).**


*Morphological evolution of size, scales and shape in Porthidium neotropical hognose pitvipers.*  
Current Zoology. (In review)

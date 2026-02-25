
# ------------------------------------------------------------------------------------------------ #

### Title: ### Evolutionary models and phylogenetic signals ####
### Author: Patrón-Rivero, C. ####
### Date: 25/02/2026 ###
### Project: "Integration and modular patterns shaping the morphological evolution in Neotropical hognose pitvipers Porthidium" ###
 
# ------------------------------------------------------------------------------------------------ #

# Libraries #

# ------------------------------------------------------------------------------------------------ #

library(dplyr)
library(tidyr)
library(phytools)
library(geiger)

# ------------------------------------------------------------------------------------------------ #

# Inputs #

# ------------------------------------------------------------------------------------------------ #

tr <- read.tree("D:/1_morpho_evol/P_tree.newick")
adu <- read.csv("D:/1_morpho_evol/adu_traits_mean_pca_spp.csv")
order_indices <- c(8, 1, 4, 7, 5, 6, 3, 2)
rownames(adu) <- adu[, 1]
adu <- adu[-1]
adu <- adu[order_indices, , drop = FALSE]
juv <- read.csv("D:/1_morpho_evol/juv_traits_mean_pca_spp.csv")
j_tr <- drop.tip(tr, setdiff(tr$tip.label, c("Porthidium_dunni", "Porthidium_lansbergii", "Porthidium_nasutum", 
											 "Porthidium_ophryomegas", "Porthidium_porrasi", "Porthidium_yucatanicum")))
row.names(juv) <- juv[, 1]
juv <- juv[-1]
order_j <- c(6, 2, 5, 3, 4, 1)
juv <- juv[order_j, , drop = FALSE]

# ------------------------------------------------------------------------------------------------ #

# Phylogenetic signal #

# ------------------------------------------------------------------------------------------------ #

adu_ps <- data.frame(Trait = character(),
                     Blombergs_K = numeric(),
                     P_value = numeric(),
                     Pagel_lambda = numeric(),
                     p_value = numeric(),
                     stringsAsFactors = FALSE)

for (i in 1:ncol(adu)) {
  trait_name <- colnames(adu)[i]
  tryCatch({
    k_result <- phylosig(tr, adu[, i], method = "K", test = TRUE, nsim = 10000)
    val_k <- k_result$K
    p_K <- k_result$P

    l_result <- phylosig(tr, adu[, i], method = "lambda", test = TRUE, nsim = 10000)
    val_l <- l_result$lambda
    p_l <- l_result$P

    adu_ps <- rbind(adu_ps, data.frame(Trait = trait_name,
                                       Blombergs_K = val_k,
                                       P_value = p_K,
                                       Pagel_lambda = val_l,
                                       p_value = p_l))
  }, error = function(e) {
    message(paste("Error in trait:", trait_name, "-", e$message))
  })
}

adu_ps[, 1] <- c("Morphometry", "Scales", "Body morphometry", 
			 "Head morphometry", "Body scales", "Head scales", "Shape")

juv_ps <- data.frame(Trait = character(),
                     Blombergs_K = numeric(),
                     P_value = numeric(),
                     Pagel_lambda = numeric(),
                     p_value = numeric(),
                     stringsAsFactors = FALSE)

for (i in 1:ncol(juv)) {
  trait_name <- colnames(juv)[i]
  tryCatch({
    k_result <- phylosig(j_tr, juv[, i], method = "K", test = TRUE, nsim = 10000)
    val_k <- k_result$K
    p_K <- k_result$P

    l_result <- phylosig(j_tr, juv[, i], method = "lambda", test = TRUE, nsim = 10000)
    val_l <- l_result$lambda
    p_l <- l_result$P

    juv_ps <- rbind(juv_ps, data.frame(Trait = trait_name,
                                       Blombergs_K = val_k,
                                       P_value = p_K,
                                       Pagel_lambda = val_l,
                                       p_value = p_l))
  }, error = function(e) {
    message(paste("Error in trait:", trait_name, "-", e$message))
  })
}

juv_ps[, 1] <- c("Morphometry", "Scales", "Body morphometry", 
			 "Head morphometry", "Body scales", "Head scales", "Shape")

# ------------------------------------------------------------------------------------------------ #

# Evolutionary models #

# ------------------------------------------------------------------------------------------------ #

aic_adu <- data.frame(Trait = character(), BM = numeric(), BM_w = numeric(), OU = numeric(), OU_w = numeric(), 
					  EB = numeric(), EB = numeric(), stringsAsFactors = FALSE)

for (i in 1:ncol(adu)) {
  trait_name <- colnames(adu)[i]
  conBM <- fitContinuous(tr, adu[i], model = "BM", control = list(niter = 1000))
  conOU <- fitContinuous(tr, adu[i], model = "OU", control = list(niter = 1000))
  conEB <- fitContinuous(tr, adu[i], model = "EB", control = list(niter = 1000))
  aic.con <- setNames(c(conBM$opt$aicc, conOU$opt$aicc, conEB$opt$aicc), c("BM", "OU", "EB"))
  aic_w <- aic.w(aic.con)
  aic_adu <- rbind(aic_adu, data.frame(Trait = trait_name, BM = aic.con["BM"], BM_w = aic_w["BM"],
									   OU = aic.con["OU"], OU_w = aic_w["OU"], EB = aic.con["EB"], EB_w = aic_w["EB"]))
}

t_adu <- cbind(adu_ps, aic_adu[-1])
rownames(t_adu) <- NULL

aic_juv <- data.frame(Trait = character(), BM = numeric(), OU = numeric(), EB = numeric(), stringsAsFactors = FALSE)

for (i in 1:ncol(juv)) {
  trait_name <- colnames(juv)[i]
  conBM <- fitContinuous(j_tr, juv[i], model = "BM", control = list(niter = 1000))
  conOU <- fitContinuous(j_tr, juv[i], model = "OU", control = list(niter = 1000))
  conEB <- fitContinuous(j_tr, juv[i], model = "EB", control = list(niter = 1000))
  aic.con <- setNames(c(conBM$opt$aicc, conOU$opt$aicc, conEB$opt$aicc), c("BM", "OU", "EB"))
  aic_w <- aic.w(aic.con)
  aic_juv <- rbind(aic_juv, data.frame(Trait = trait_name, BM = aic.con["BM"], BM_w = aic_w["BM"],
									   OU = aic.con["OU"], OU_w = aic_w["OU"], EB = aic.con["EB"], EB_w = aic_w["EB"]))
}

t_juv <- cbind(juv_ps, aic_juv[-1])
rownames(t_juv) <- NULL

# ------------------------------------------------------------------------------------------------ #

### EndNotRun

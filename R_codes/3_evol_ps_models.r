
# ------------------------------------------------------------------------------------------------ #

### Title: Evolutionary models and phylogenetic signals ####
### Author: Patrón-Rivero, C. ####
### Date: 07/16/2025 ###
### Project: "Morphological evolution of size, scales and shape in Porthidium neotropical hognose pitvipers" ###

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

tr <- read.tree("E:/1_Morpho_System/P_tree.newick")
adu <- read.csv("E:/1_Morpho_System/1_Paper/data_github/adu_traits_mean_pca_spp.csv")
order_indices <- c(8, 1, 4, 7, 5, 6, 3, 2)
rownames(adu) <- adu[, 1]
adu <- adu[-1]
adu <- adu[order_indices, , drop = FALSE]
juv <- read.csv("E:/1_Morpho_System/1_Paper/data_github/juv_traits_mean_pca_spp.csv")
j_tr <- drop.tip(tr, setdiff(tr$tip.label, c("P_dunni", "P_lansbergii", "P_nasutum", "P_ophryomegas", "P_porrasi", "P_yucatanicum")))
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

aic_adu <- data.frame(Trait = character(), BM = numeric(), OU = numeric(), EB = numeric(), stringsAsFactors = FALSE)

for (i in 1:ncol(adu)) {
  trait_name <- colnames(adu)[i]
  conBM <- fitContinuous(tr, adu[i], model = "BM", control = list(niter = 1000))
  conOU <- fitContinuous(tr, adu[i], model = "OU", control = list(niter = 1000))
  conEB <- fitContinuous(tr, adu[i], model = "EB", control = list(niter = 1000))
  aic.con <- setNames(c(conBM$opt$aicc, conOU$opt$aicc, conEB$opt$aicc), c("BM", "OU", "EB"))
  aic_adu <- rbind(aic_adu, data.frame(Trait = trait_name, BM = aic.con["BM"], OU = aic.con["OU"], EB = aic.con["EB"]))
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
  aic_juv <- rbind(aic_juv, data.frame(Trait = trait_name, BM = aic.con["BM"], OU = aic.con["OU"], EB = aic.con["EB"]))
}

t_juv <- cbind(juv_ps, aic_juv[-1])
rownames(t_juv) <- NULL

# ------------------------------------------------------------------------------------------------ #

### EndNotRun

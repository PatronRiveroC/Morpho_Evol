
# ------------------------------------------------------------------------------------------------ #

### Title: Standarized Major-Axis analysis ####
### Author: Patrón-Rivero, C. ####
### Date: 07/16/2025 ###
### Project: "Morphological evolution of size, scales and shape in Porthidium neotropical hognose pitvipers" ###

# ------------------------------------------------------------------------------------------------ #

# Libraries #

# ------------------------------------------------------------------------------------------------ #

library(smatr)

# ------------------------------------------------------------------------------------------------ #

# Inputs #

# ------------------------------------------------------------------------------------------------ #

adu_pc <- read.csv("E:/1_Morpho_System/1_Paper/data_github/adu_pc.csv")[-1]
juv_pc <- read.csv("E:/1_Morpho_System/1_Paper/data_github/juv_pc.csv")[-1]

# ------------------------------------------------------------------------------------------------ #

# SMA models #

# ------------------------------------------------------------------------------------------------ #

compare_pc_models <- function(all_pc) {
  models <- list()
  count <- 1
  
  for (i in 1:(ncol(all_pc) - 1)) {
    for (j in (i + 1):ncol(all_pc)) {
      model <- sma(all_pc[, i] ~ all_pc[, j], na.action = "na.omit")
      m <- vector()
      m[1] <- coef(model)[1]
      m[2] <- coef(model)[2]
      m[3] <- model$r
      m[4] <- model$p
      models[[count]] <- m
      count <- count + 1
    }
  }
  
  results_df <- do.call(rbind, models)
  colnames(results_df) <- c("Elevation", "Slope", "R_squared", "P_value")
  
  comb <- combn(colnames(all_pc), 2, simplify = FALSE)
  comparisons <- data.frame(do.call(rbind, comb), stringsAsFactors = FALSE)
  colnames(comparisons) <- c("Var1", "Var2")
  
  final_results <- cbind(comparisons, results_df)
  f <- data.frame(lapply(final_results, function(x) if(is.list(x)) unlist(x) else x))
  
  return(f)
}

sma_adu <- compare_pc_models(adu_pc)
sma_juv <- compare_pc_models(juv_pc)

# ------------------------------------------------------------------------------------------------ #

### EndNotRun

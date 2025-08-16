
# ------------------------------------------------------------------------------------------------ #

### Title: Rates of trait evolution ####
### Author: Patrón-Rivero, C. ####
### Date: 16/07/2025 ###
### Project: "Morphological evolution of size, scales and shape in Porthidium neotropical hognose pitvipers" ###

# ------------------------------------------------------------------------------------------------ #

# Libraries #

# ------------------------------------------------------------------------------------------------ #

library(phytools)
library(evomap)
library(tidyr)
library(dplyr)
library(rstatix)
library(ggplot2)
library(ggpubr)
library(patchwork)

# ------------------------------------------------------------------------------------------------ #

# Inputs adults #

# ------------------------------------------------------------------------------------------------ #

data <- read.csv("E:/1_Morpho_System/1_Paper/data_github/adu_traits_mean_pca_spp.csv", row.names = 1)
tree <- read.tree("E:/1_Morpho_System/P_tree.newick")

# ------------------------------------------------------------------------------------------------ #

# Rate analysis using MCMC #

# ------------------------------------------------------------------------------------------------ #

rate_list <- list()
			
for (i in 1:7) {
  x <- data[[i]]
  BMsigma2 <- ace(x, tree, method = "REML")$sigma2[1]
  dat <- x
  names(dat) <- rownames(data)
  
  mvBMresults <- mvBM(dat, tree, BMsigma2)
  tree_mvBM <- tree
  tree_mvBM$edge.length <- mvBMresults$rBL
    
  iterations <- 500000
  model_mvBM <- anc.Bayes(tree_mvBM, dat, ngen = iterations)
  MCMC_mvBM_sigma2 <- model_mvBM$mcmc[(iterations * 0.002):nrow(model_mvBM$mcmc), 2]
  rate_list[[colnames(data)[i]]] <- MCMC_mvBM_sigma2
}

rate_df <- as.data.frame(rate_list)
rate_adu <- pivot_longer(rate_df, cols = everything(), names_to = "Variable", values_to = "Rates")

# ------------------------------------------------------------------------------------------------ #

# Inputs juv #

# ------------------------------------------------------------------------------------------------ #

data <- read.csv("E:/1_Morpho_System/1_Paper/data_github/juv_traits_mean_pca_spp.csv", row.names = 1)
tree <- drop.tip(tree, setdiff(tree$tip.label, c("P_dunni", "P_lansbergii", "P_nasutum", "P_ophryomegas", "P_porrasi", "P_yucatanicum")))

# ------------------------------------------------------------------------------------------------ #

# Rate analysis using MCMC #

# ------------------------------------------------------------------------------------------------ #

rate_list <- list()
			
for (i in 1:7) {
  x <- data[[i]]
  BMsigma2 <- ace(x, tree, method = "REML")$sigma2[1]
  dat <- x
  names(dat) <- rownames(data)
  
  mvBMresults <- mvBM(dat, tree, BMsigma2)
  tree_mvBM <- tree
  tree_mvBM$edge.length <- mvBMresults$rBL
    
  iterations <- 500000
  model_mvBM <- anc.Bayes(tree_mvBM, dat, ngen = iterations)
  MCMC_mvBM_sigma2 <- model_mvBM$mcmc[(iterations * 0.002):nrow(model_mvBM$mcmc), 2]
  rate_list[[colnames(data)[i]]] <- MCMC_mvBM_sigma2
}

rate_df <- as.data.frame(rate_list)
rate_juv <- pivot_longer(rate_df, cols = everything(), names_to = "Variable", values_to = "Rates")

# ------------------------------------------------------------------------------------------------ #

# Plot #

# ------------------------------------------------------------------------------------------------ #

data_combined <- bind_rows(
  mutate(rate_adu, stage = "adu"),
  mutate(rate_juv, stage = "juv")
)

data <- as.data.frame(data_combined)
data["Rates"] <- log(data["Rates"])
data[[1]] <- as.factor(data[[1]])

stat.test <- data %>%
  group_by(stage) %>%
  t_test(Rates ~ Variable) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test

stat.test <- stat.test %>%
  add_xy_position(x = "stage", dodge = 0.8)
  
p <- ggboxplot(data, x = "Variable", y = "Rates",
                 color = "stage", palette = "jco") +
    xlab("") +
    ylab(expression(log (sigma^2))) +
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(labels = c(
      "Alometry" = "Morphometry", 
	  "Folidosis" = "Scales", 
      "Body_alometry" = "Body morphometry", 
      "Body_folidosis" = "Body scales", 
      "Head_alometry" = "Head morphometry", 
      "Head_folidosis" = "Head scales",
	  "Shape" = "Shape"
    )) +
    stat_pvalue_manual(
  stat.test,  label = "p", tip.length = 0
  )


plot_boxplot_with_stats <- function(data) {
  data <- as.data.frame(data)
  data["Rates"] <- log(data["Rates"])
  data[[1]] <- as.factor(data[[1]])

  cmpr <- list(
    c("Alometry", "Folidosis"),
    c("Body_alometry", "Body_folidosis"),
    c("Head_alometry", "Head_folidosis"),
    c("Body_alometry", "Head_alometry"),
    c("Body_folidosis", "Head_folidosis"),
    c("Body_alometry", "Head_folidosis"),
    c("Body_folidosis", "Head_alometry"),
	c("Alometry", "Shape"),
	c("Folidosis", "Shape"),
	c("Body_alometry", "Shape"),
	c("Body_folidosis", "Shape"),
	c("Head_alometry", "Shape"),
	c("Head_folidosis", "Shape")
  )

  p <- ggboxplot(data, x = "Variable", y = "Rates",
                 color = "Variable", palette = "jco") +
    xlab("") +
    ylab(expression(log (sigma^2))) +
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(labels = c(
      "Alometry" = "Morphometry", 
	  "Folidosis" = "Scales", 
      "Body_alometry" = "Body morphometry", 
      "Body_folidosis" = "Body scales", 
      "Head_alometry" = "Head morphometry", 
      "Head_folidosis" = "Head scales",
	  "Shape" = "Shape"
    )) +
    ggpubr::stat_compare_means(comparisons = cmpr, tip.length = 0.01,
                               label = "p.signif", 
                               symnum.args = list(
                                 cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                 symbols = c("*", "*", "*", "*", "")
                               ))
  return(p)
}

a <- plot_boxplot_with_stats(rate_adu)
j <- plot_boxplot_with_stats(rate_juv)

p <- j + a

p1 <- p + plot_annotation(
  tag_levels = list(c("A)", "B)")),
  theme = theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold")))

setwd("E:/1_Morpho_System/1_Paper/Figs")
ggsave(file = "Fig_3.jpg", plot = p1, width = 20, height = 12, dpi = 600, units = "cm", device = "jpg")
ggsave(file = "Fig_3.pdf", plot = p1, width = 20, height = 12, dpi = 600, units = "cm", device = "pdf")
ggsave(file = "Fig_3.tiff", plot = p1, width = 20, height = 12, dpi = 600, units = "cm", device = "tiff")

# ------------------------------------------------------------------------------------------------ #

### EndNotRun

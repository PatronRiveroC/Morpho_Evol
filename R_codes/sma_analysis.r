
### Title: SMA analysis ####
### Author: Patrón-Rivero, C. ####
### Date: 25/02/2026 ###
### Project: "Integration and modular patterns shaping the morphological evolution in Neotropical hognose pitvipers Porthidium" ###

# ------------------------------------------------------------------------------------------------ #

# Libraries #

# ------------------------------------------------------------------------------------------------ #

library(readr)
library(tidyverse)
library(geomorph)
library(smatr)
library(ggplot2)
library(patchwork)

# ------------------------------------------------------------------------------------------------ #

# Inputs #

# ------------------------------------------------------------------------------------------------ #

base_dir <- "D:/1_morpho_evol/2round/data"

# ------------------------------------------------------------------------------------------------ #

# PCA #

# ------------------------------------------------------------------------------------------------ #

run_stage_pca <- function(data, stage_value) {

  df <- data %>%
    filter(Stage == stage_value) %>%
    select(where(is.numeric)) %>%
    drop_na()

  pca <- prcomp(df, center = TRUE, scale. = TRUE)

  return(pca)
}

run_file_pcas <- function(file_path) {

  data <- read.csv(file_path)

  list(
    juvenile = run_stage_pca(data, "juvenile"),
    adult    = run_stage_pca(data, "adult")
  )
}

pcas <- list(

  mor  = run_file_pcas(file.path(base_dir, "mor.csv")),
  lep  = run_file_pcas(file.path(base_dir, "lep.csv")),
  body_mor = run_file_pcas(file.path(base_dir, "body_mor.csv")),
  body_lep = run_file_pcas(file.path(base_dir, "body_lep.csv")),
  head_mor = run_file_pcas(file.path(base_dir, "head_mor.csv")),
  head_lep = run_file_pcas(file.path(base_dir, "head_lep.csv"))
)

get_scores_95 <- function(pca, threshold = 0.95) {

  var_exp <- pca$sdev^2 / sum(pca$sdev^2)
  cum_var <- cumsum(var_exp)

  n_pc <- which(cum_var >= threshold)[1]

  scores <- pca$x[, 1:n_pc, drop = FALSE]

  return(scores)
}

Size_juvenile <- get_scores_95(pcas$mor$juvenile)
Size_adult    <- get_scores_95(pcas$mor$adult)
Scales_juvenile <- get_scores_95(pcas$lep$juvenile)
Scales_adult    <- get_scores_95(pcas$lep$adult)
Body_size_juvenile <- get_scores_95(pcas$body_mor$juvenile)
Body_size_adult    <- get_scores_95(pcas$body_mor$adult)
Body_scales_juvenile <- get_scores_95(pcas$body_lep$juvenile)
Body_scales_adult    <- get_scores_95(pcas$body_lep$adult)
Head_size_juvenile <- get_scores_95(pcas$head_mor$juvenile)
Head_size_adult    <- get_scores_95(pcas$head_mor$adult)
Head_scales_juvenile <- get_scores_95(pcas$head_lep$juvenile)
Head_scales_adult    <- get_scores_95(pcas$head_lep$adult)

rename_pca_scores <- function(scores, prefix) {
  colnames(scores) <- paste0(prefix, "_pc", seq_len(ncol(scores)))
  return(scores)
}

Size_juvenile  <- rename_pca_scores(Size_juvenile,  "Size")
Size_adult     <- rename_pca_scores(Size_adult,     "Size")
Scales_juvenile <- rename_pca_scores(Scales_juvenile, "Scales")
Scales_adult    <- rename_pca_scores(Scales_adult,    "Scales")
Body_size_juvenile  <- rename_pca_scores(Body_size_juvenile,  "Body_size")
Body_size_adult     <- rename_pca_scores(Body_size_adult,     "Body_size")
Body_scales_juvenile <- rename_pca_scores(Body_scales_juvenile, "Body_scales")
Body_scales_adult    <- rename_pca_scores(Body_scales_adult,    "Body_scales")
Head_size_juvenile  <- rename_pca_scores(Head_size_juvenile,  "Head_size")
Head_size_adult     <- rename_pca_scores(Head_size_adult,     "Head_size")
Head_scales_juvenile <- rename_pca_scores(Head_scales_juvenile, "Head_scales")
Head_scales_adult    <- rename_pca_scores(Head_scales_adult,    "Head_scales")

tps_file <- "D:/1_morpho_evol/Porthidium_landmarks.TPS"
id_file  <- "D:/1_morpho_evol/2round/data/ID_GM.csv"
coords_raw <- readland.tps(tps_file, specID = "ID", negNA = TRUE)
coords_raw <- estimate.missing(coords_raw, method = "TPS")
ID_GM <- read.csv(id_file)
ids_coords <- dimnames(coords_raw)[[3]]
ids_coords_clean <- trimws(ids_coords)

ID_GM <- ID_GM %>%
  mutate(
    ID = trimws(ID),
    Stage = trimws(Stage)
  ) %>%
  filter(!is.na(Stage))
  
ID_GM <- ID_GM %>%
  filter(ID %in% ids_coords_clean) %>%
  arrange(match(ID, ids_coords_clean))

coords_juvenile <- coords_raw[, , ids_coords_clean %in% ID_GM$ID[ID_GM$Stage == "juvenile"]]
coords_adult    <- coords_raw[, , ids_coords_clean %in% ID_GM$ID[ID_GM$Stage == "adult"]]
gpa_juvenile <- gpagen(coords_juvenile, print.progress = TRUE)
gpa_adult    <- gpagen(coords_adult,    print.progress = TRUE)
pca_juvenile <- gm.prcomp(gpa_juvenile$coords)
pca_adult    <- gm.prcomp(gpa_adult$coords)
Shape_juvenile <- pca_juvenile$x[, 1:10]
colnames(Shape_juvenile) <- paste0("Shape_pc", 1:10)
Shape_adult <- pca_adult$x[, 1:13]
colnames(Shape_adult) <- paste0("Shape_pc", 1:13)

juvenile_df <- cbind(
  Size_juvenile[1:164, ],
  Scales_juvenile[1:164, ],
  Body_size_juvenile[1:164, ],
  Body_scales_juvenile[1:164, ],
  Head_size_juvenile[1:164, ],
  Head_scales_juvenile[1:164, ],
  Shape_juvenile[1:164, ]
)

adult_df <- cbind(
  Size_adult[1:191, ],
  Scales_adult[1:191, ],
  Body_size_adult[1:191, ],
  Body_scales_adult[1:191, ],
  Head_size_adult[1:191, ],
  Head_scales_adult[1:191, ],
  Shape_adult[1:191, ]
)

juvenile_df <- as.data.frame(juvenile_df)
adult_df    <- as.data.frame(adult_df)

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

sma_adu <- compare_pc_models(adult_df)
sma_juv <- compare_pc_models(juvenile_df)

# ------------------------------------------------------------------------------------------------ #

# Plots #

# ------------------------------------------------------------------------------------------------ #

adult_df <- read_csv("D:/1_morpho_evol/sma_adult.csv")
adult_df$sig <- ifelse(!is.na(adult_df$P_value) & adult_df$P_value <= 0.05, "*", "")

var_order <- c(
  paste0("Size_pc", 1:3),
  paste0("Scales_pc", 1:7),
  paste0("Shape_pc", 1:13),
  paste0("Body_size_pc", 1:3),
  paste0("Body_scales_pc", 1:3),
  "Head_size_pc1",
  paste0("Head_scales_pc", 1:5)
)

adult_df <- adult_df %>%
  mutate(
    Var1 = factor(Var1, levels = var_order),
    Var2 = factor(Var2, levels = var_order)
  )
  
a <-ggplot(adult_df, aes(Var1, Var2, fill = P_value)) +
  geom_tile(color = "black", linewidth = 0.3) +
  geom_text(aes(label = sig), size = 5) +
  scale_fill_gradient2(low = "#f4a582", high = "#74a9cf", limits = c(0, 1),) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = rev) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 0),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

juvenile_df <- read_csv("D:/1_morpho_evol/sma_juvenile.csv")
juvenile_df$sig <- ifelse(!is.na(juvenile_df$P_value) & juvenile_df$P_value <= 0.05, "*", "")

var_order <- c(
  paste0("Size_pc", 1:8),
  paste0("Scales_pc", 1:8),
  paste0("Shape_pc", 1:10),
  paste0("Body_size_pc", 1:5),
  paste0("Body_scales_pc", 1:3),
  paste0("Head_size_pc", 1:3),
  paste0("Head_scales_pc", 1:5)
)

juvenile_df <- juvenile_df %>%
  mutate(
    Var1 = factor(Var1, levels = var_order),
    Var2 = factor(Var2, levels = var_order)
  )
  
j <- ggplot(juvenile_df, aes(Var1, Var2, fill = P_value)) +
  geom_tile(color = "black", linewidth = 0.3) +
  geom_text(aes(label = sig), size = 5) +
  scale_fill_gradient2(low = "#f4a582", high = "#74a9cf", limits = c(0, 1),) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = rev) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 0),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

p <- a + j

p1 <- p + plot_annotation(
  tag_levels = list(c("A)", "B)")),
  theme = theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"))
)

setwd("D:/1_morpho_evol/2round/Figs/Final")
ggsave(file = "Fig_2.pdf", plot = p1, width = 30, height = 15, dpi = 1200, units = "cm", device = "pdf")

# ------------------------------------------------------------------------------------------------ #

# MOdel selection approach #

# ------------------------------------------------------------------------------------------------ #

compare_integration <- function(data, y_var, x_var) {
  
  df <- data[, c(y_var, x_var)]
  df <- na.omit(df)
  
  if (nrow(df) < 5) {
    return(NULL)
  }
  
  model_int <- lm(df[[y_var]] ~ df[[x_var]]) 
  model_mod <- lm(df[[y_var]] ~ 1)           
  aic_int <- AIC(model_int)
  aic_mod <- AIC(model_mod)
  delta_int <- aic_int - min(aic_int, aic_mod)
  delta_mod <- aic_mod - min(aic_int, aic_mod)
  w_int <- exp(-0.5 * delta_int)
  w_mod <- exp(-0.5 * delta_mod)
  w_sum <- w_int + w_mod
  w_int <- w_int / w_sum
  w_mod <- w_mod / w_sum
  
  data.frame(
    Y = y_var,
    X = x_var,
    AIC_integration = aic_int,
    AIC_modular = aic_mod,
    delta_AIC_integration = delta_int,
    delta_AIC_modular = delta_mod,
    wAIC_integration = w_int,
    wAIC_modular = w_mod,
    interpretation = interpretation
  )
}

vars <- colnames(adult_df)

adult_integration_results <- do.call(
  rbind,
  combn(
    vars,
    2,
    function(v) compare_integration(adult_df, y_var = v[1], x_var = v[2]),
    simplify = FALSE
  )
)

vars <- colnames(juvenile_df)

juv_integration_results <- do.call(
  rbind,
  combn(
    vars,
    2,
    function(v) compare_integration(juvenile_df, y_var = v[1], x_var = v[2]),
    simplify = FALSE
  )
)

write.csv(
  juv_integration_results,
  file = "D:/1_morpho_evol/juv_int_res.csv",
  row.names = FALSE
)

# ------------------------------------------------------------------------------------------------ #

### EndNotRun

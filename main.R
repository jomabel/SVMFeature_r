#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  #library(devtools)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(glue)
  library(scales)
  library(tools)

})

# ===== Helper functions =====
simulate <- function(method, data_scaled, costs, population_size, input_cols, output_col, num_features, iter_type, iter_val,
                     num_obj, clones, p_mutation, p_mut_ind, p_mut_fea, p_mut_coord, mut_coord) {
  message(glue::glue("🚀 Running NSGA-II simulation with method: {method}, population size: {population_size}, number of features: {num_features}, iteration type: {iter_type}, and iteration value: {iter_val}"))

  obj <- SVMFeature(
    data = data_scaled,
    inputs = input_cols,
    output = output_col,
    costs = costs,
    pop_size = population_size,
    num_fea = num_features,
    n_iter = iter_val,
    max_time = iter_val,
    mode = iter_type,
    objective = method,
    num_obj = num_obj,
    clones = clones,
    p_mutation = p_mutation,
    p_mut_ind = p_mut_ind,
    p_mut_fea = p_mut_fea,
    p_mut_coord = p_mut_coord,
    mut_coord = mut_coord
  )

  result <- run.SVMFeature(obj)
  message(glue::glue("✅ NSGA-II simulation completed with method: {method}"))
  return(result)
}

# ======= Parameter configuration =======
dataset        <- "data/09Colon_p.txt"
has_header     <- FALSE         # Set to TRUE if the file has a header
has_index      <- TRUE          # Set to TRUE if the file has an index column
has_costs      <- TRUE          # Set to TRUE if the dataset includes feature costs
separator      <- ";"           # Delimiter (e.g., "," or ";")
should_normalize <- TRUE        # Whether to normalize data
method         <- "distance-epsilon" # "distance-epsilon" or "confusion-matrix"
population_size <- 10
num_features   <- 5
iter_type      <- "iters"       # "iters" for iterations, "time" for seconds
iter_val       <- 5             # Number of iterations or time limit in seconds

# Advanced Population parameters
num_obj        <- 2
clones         <- 0
p_mutation     <- 0.7
p_mut_ind      <- 0.4
p_mut_fea      <- 0.4
p_mut_coord    <- 0.2
mut_coord      <- 0

# ====== Dataset loading ======
message("Loading dataset...")

data <- tryCatch({
  read_delim(dataset, delim = separator, col_names = has_header)
}, error = function(e) {
  stop(glue::glue("Error loading dataset: {e$message}"))
})

if (has_index) {
  data <- data[,-1]
}

# ====== Data preparation ======
output_col <- names(data)[ncol(data)]
input_cols <- names(data)[-ncol(data)]

# Cargar costes si existen
costs <- if (has_costs) rep(0, length(input_cols)) else NULL

# Normalización opcional
if (should_normalize) {
  # Normalization is handled internally by SVMFeature, but additional pre-processing can be done here if needed.
}

message(glue::glue("Dataset loaded: {nrow(data)} rows and {ncol(data)} columns"))

# ====== Algorithm execution ======
# Loading the library
devtools::load_all(".") 
library(SVMFeature)
library(scales)

# ===== Simulation =====
result <- simulate(method, data, costs, population_size, input_cols, output_col, num_features, iter_type, iter_val,
                   num_obj, clones, p_mutation, p_mut_ind, p_mut_fea, p_mut_coord, mut_coord)

# ===== Save results =====
dataset_name <- tools::file_path_sans_ext(basename(dataset))
save_results(result, dataset_name = dataset_name)

# ===== Visualization =====
plot(result)

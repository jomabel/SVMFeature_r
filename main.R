#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(glue)
  library(scales)
  library(tools)
})

# ===== Helper functions =====
simulate <- function(method, data_scaled, population_size, input_cols, output_col, num_features, iter_type, iter_val,
                     clones, p_mutation, p_mut_ind, p_mut_fea, p_mut_coord, mut_coord) {
  message(glue::glue("Running NSGA-II simulation with method: {method}, population size: {population_size}, number of features: {num_features}, iteration type: {iter_type}, and iteration value: {iter_val}"))

  obj <- svm_fs_meta(
    data = data_scaled,
    inputs = input_cols,
    output = output_col,
    pop_size = population_size,
    num_fea = num_features,
    n_iter = iter_val,
    max_time = iter_val,
    mode = iter_type,
    objective = method,
    clones = clones,
    p_mutation = p_mutation,
    p_mut_ind = p_mut_ind,
    p_mut_fea = p_mut_fea,
    p_mut_coord = p_mut_coord,
    mut_coord = mut_coord
  )

  result <- run.svm_fs_meta(obj)
  message(glue::glue("NSGA-II simulation completed with method: {method}"))
  return(result)
}

# ======= Parameter configuration =======
dataset         <- "Colon"
method          <- "DIST-EPS"
population_size <- 10
num_features    <- 5
iter_type       <- "iters"
iter_val        <- 5
normalize       <- TRUE

# Advanced Population parameters (using library defaults)
clones          <- 0
p_mutation      <- 0.7
p_mut_ind       <- 0.4
p_mut_fea       <- 0.4
p_mut_coord     <- 0.2
mut_coord       <- 0

# ====== Dataset loading and preprocessing ======
message("Loading dataset...")

if (!requireNamespace("plsgenomics", quietly = TRUE)) {
  stop("Install the example data package with install.packages('plsgenomics').",
       call. = FALSE)
}

# Load the example from its source package without attaching it.
dataset_env <- new.env(parent = baseenv())
utils::data(list = dataset, package = "plsgenomics", envir = dataset_env)
colon_data <- dataset_env[[dataset]]

# plsgenomics codes normal tissue as 1 and tumor tissue as 2.
# Keep the previous convention: normal = 1 and tumor = -1.
data <- data.frame(y = ifelse(colon_data$Y == 2, -1, 1), colon_data$X)
names(data) <- c("y", paste0("x", seq_len(ncol(colon_data$X))))

output_col <- names(data)[1]
input_cols <- names(data)[-1]

# Normalize input features (Min-Max scaling to [0, 1]).
# Normalization is performed outside the library.
if (normalize) {
  scaler <- function(x) {
    rng <- range(x, na.rm = TRUE)
    if (rng[1] == rng[2]) return(rep(0.5, length(x)))
    return((x - rng[1]) / (rng[2] - rng[1]))
  }
  data_scaled <- data
  data_scaled[, input_cols] <- as.data.frame(lapply(data[, input_cols], scaler))
} else {
  data_scaled <- data
}

message(glue::glue("Dataset prepared: {nrow(data_scaled)} rows and {ncol(data_scaled)} columns"))

# ===== Load svm_fs_meta package =====
# load_all() is used during development to simulate loading the installed package.
devtools::load_all(".")
library(svmfsmeta)

# ===== Algorithm execution =====
# The simulate() helper function wraps the svm_fs_meta object creation and run.svm_fs_meta() call.
result <- simulate(method, data_scaled, population_size, input_cols, output_col, num_features, iter_type, iter_val,
                   clones, p_mutation, p_mut_ind, p_mut_fea, p_mut_coord, mut_coord)

# ===== Result Export and Visualization =====
# save_results() exports the Pareto front to a CSV file in the 'results/' directory.
dataset_name <- tolower(dataset)
save_results(result, "results", dataset_name = dataset_name)

# draw_solution() provides a visual representation of the Pareto front using ggplot2.
draw_solution(result)

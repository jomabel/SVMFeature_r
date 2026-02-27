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

# ===== Funciones auxiliares =====
simulate <- function(method, data_scaled, costs, population_size, input_cols, output_col, num_features, iter_type, iter_val) {
  cat(glue::glue("🚀 Running NSGA-II simulation with method: {method}, population size: {population_size}, number of features: {num_features}, iteration type: {iter_type}, and iteration value: {iter_val}", "\n"))

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
    objective = method
  )

  result <- run.SVMFeature(obj)
  cat(glue::glue("✅ NSGA-II simulation completed with method: {method}"), "\n")
  return(result)
}

save_results <- function(result, method, dataset_name, population_size, iter_type, iter_val) {
  # Extraer resultados del frente de Pareto
  pareto_front <- result$population$df_solutions %>% dplyr::filter(FRONT == 1)

  # Convertir columnas tipo lista a texto
  pareto_front <- pareto_front %>%
    dplyr::mutate(across(where(is.list), ~ sapply(., function(x) paste(x, collapse = ","))))

  output_dir <- file.path("results", dataset_name, method)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  base_filename <- paste0(population_size, "_", iter_type, "_", iter_val, ".csv")
  output_file <- file.path(output_dir, base_filename)

  counter <- 1
  while (file.exists(output_file)) {
    output_file <- file.path(output_dir, paste0(population_size, "_", iter_type, "_", iter_val, "_", counter, ".csv"))
    counter <- counter + 1
  }

  write.csv(pareto_front, output_file, row.names = FALSE)
  cat(glue::glue("📁 Pareto front saved in: {output_file}\n"))
}


# ======= Configuración de argumentos =======
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Uso: Rscript main.R <dataset> [has_header yes/no] [has_index yes/no] [has_costs yes/no] [sep ,] [normalize yes/no] [method] [population_size] [num_features] [iter_type] [iter_val]")
}

dataset <- args[1]
has_header <- ifelse(length(args) >= 2, tolower(args[2]) %in% c("yes", "y"), FALSE)
has_index <- ifelse(length(args) >= 3, tolower(args[3]) %in% c("yes", "y"), TRUE)
has_costs <- ifelse(length(args) >= 4, tolower(args[4]) %in% c("yes", "y"), TRUE)
separator <- ifelse(length(args) >= 5, args[5], ",")
normalize <- ifelse(length(args) >= 6, tolower(args[6]) %in% c("yes", "y"), TRUE)
method <- ifelse(length(args) >= 7, args[7], "distance-epsilon")
population_size <- ifelse(length(args) >= 8, as.integer(args[8]), 10)
num_features <- ifelse(length(args) >= 9, as.integer(args[9]), 5)
iter_type <- ifelse(length(args) >= 10, args[10], "iters")
iter_val <- ifelse(length(args) >= 11, as.integer(args[11]), 5)

# ====== Lectura del dataset ======
cat("📂 Loading dataset...\n")

header_opt <- if (has_header) TRUE else FALSE
data <- read_delim(dataset, delim = separator, col_names = header_opt)

if (has_index) {
  data <- data[,-1]
}

if (!has_header) {
  num_cols <- ncol(data)
  names(data) <- c("y", paste0("x", 1:(num_cols-1)))
}

output_col <- names(data)[1]
input_cols <- names(data)[-1]

if (has_costs) {
  costs <- as.numeric(data[1, -1])
  data <- data[-1, ]
} else {
  costs <- rep(0.0, length(input_cols))
}

if (normalize) {
  data[input_cols] <- as.data.frame(lapply(data[input_cols], scales::rescale))
}

# ===== Carga del paquete SVMFeature =====
devtools::load_all(".")

# ===== Ejecución del algoritmo =====
result <- simulate(method, data, costs, population_size, input_cols, output_col, num_features, iter_type, iter_val)

# ===== Guardar resultados =====
dataset_name <- tools::file_path_sans_ext(basename(dataset))
save_results(result, method, dataset_name, population_size, iter_type, iter_val)

# Run after installing svmfsmeta and the suggested data package plsgenomics.
library(svmfsmeta)
if (!requireNamespace("plsgenomics", quietly = TRUE)) {
  stop("Install the example data with install.packages('plsgenomics').")
}
dataset_env <- new.env(parent = baseenv())
utils::data("Colon", package = "plsgenomics", envir = dataset_env)
colon <- dataset_env$Colon
data <- data.frame(y = ifelse(colon$Y == 2, -1, 1), colon$X)
names(data) <- c("y", paste0("x", seq_len(ncol(colon$X))))
inputs <- names(data)[-1]
set.seed(42)

# Stratify while allocating exactly ceiling(n / 3) observations to test.
groups <- split(seq_len(nrow(data)), data$y)
quotas <- lengths(groups) * ceiling(nrow(data) / 3) / nrow(data)
test_sizes <- floor(quotas)
remaining <- ceiling(nrow(data) / 3) - sum(test_sizes)
if (remaining > 0) {
  extra <- order(quotas - test_sizes, decreasing = TRUE)[seq_len(remaining)]
  test_sizes[extra] <- test_sizes[extra] + 1
}
test_rows <- unlist(Map(function(rows, n) rows[sample.int(length(rows), n)],
                        groups, test_sizes), use.names = FALSE)
train <- data[-test_rows, , drop = FALSE]
test <- data[test_rows, , drop = FALSE]

# Learn scaling from training only. Constants use the R example convention 0.5.
training_min <- vapply(train[inputs], min, numeric(1))
training_span <- vapply(train[inputs], max, numeric(1)) - training_min
for (feature in inputs) {
  if (training_span[[feature]] == 0) {
    train[[feature]] <- rep(0.5, nrow(train))
    test[[feature]] <- rep(0.5, nrow(test))
  } else {
    train[[feature]] <- (train[[feature]] - training_min[[feature]]) / training_span[[feature]]
    test[[feature]] <- (test[[feature]] - training_min[[feature]]) / training_span[[feature]]
  }
}
cat("Colon:", nrow(train), "training observations;", nrow(test), "test observations.\n")
print(cbind(TRAIN = table(train$y), TEST = table(test$y)))

optimizer <- svm_fs_meta(train, inputs, "y", pop_size = 10, num_fea = 5,
                         n_iter = 5, objective = "MC")
result <- run.svm_fs_meta(optimizer)
front <- Filter(function(sol) sol$front == 1, result$population$solution_list)
cat("\nPrediction for one unlabeled test observation using the first solution:\n")
new_observation <- test[1, inputs, drop = FALSE]
print(predict_solution(front[[1]], newdata = new_observation, type = "class"))
print(predict_solution(front[[1]], newdata = new_observation, type = "score"))
results <- do.call(rbind, lapply(front, function(sol) {
  held_out <- evaluate_test(sol, test, "y")
  row <- list(SOL = sol$num)
  train_fn <- unname(sol$conf_matrix["1", "-1"])
  train_fp <- unname(sol$conf_matrix["-1", "1"])
  for (split in c("TRAIN", "TEST")) {
    metrics <- if (split == "TRAIN") sol else held_out
    matrix <- if (split == "TRAIN") sol$conf_matrix else held_out$confusion_matrix
    fn <- if (split == "TRAIN") train_fn else held_out$fn
    fp <- if (split == "TRAIN") train_fp else held_out$fp
    n <- if (split == "TRAIN") nrow(train) else nrow(test)
    values <- list(FN = fn, FP = fp, TN = unname(matrix["-1", "-1"]),
                   TP = unname(matrix["1", "1"]), N_SAMPLES = n,
                   N_ERRORS = fn + fp, ERROR_RATE = (fn + fp) / n,
                   CONFUSION_MATRIX = paste0("[[", paste(matrix[1, ], collapse = ", "),
                                             "], [", paste(matrix[2, ], collapse = ", "), "]]"))
    for (metric in c("accuracy", "precision", "recall", "f1", "kappa", "auc")) {
      values[[toupper(metric)]] <- metrics[[metric]]
    }
    names(values) <- paste0(split, "_", names(values))
    row <- c(row, values)
  }
  row$TRAIN_DIST <- sol$objective[1]
  row$TRAIN_EPS <- sol$objective[2]
  as.data.frame(row, stringsAsFactors = FALSE)
}))
for (split in c("TRAIN", "TEST")) {
  cat("\n", split, "metrics for the training Pareto front:\n", sep = "")
  print(results[c("SOL", grep(paste0("^", split, "_"), names(results), value = TRUE))],
        row.names = FALSE)
}
cat("\nMC uses FN and FP; class 1 is positive. Error rates are proportions.\n")
cat("Test data do not influence optimization or the selection of the front.\n")

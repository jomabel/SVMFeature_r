# 🧠 SVM_FS_Meta: Multi-objective SVM Classifier with Feature Selection via Metaheuristics

`svmfsmeta` is an R package for binary **soft-margin Support Vector Machine (SVM)** classification with **Feature Selection**, based on the methodologies studied by Valero-Carreras et al. (2023) and Alcaraz (2024).

The optimizer returns an object containing the final population of candidate classifiers. Its first non-dominated front approximates the **Pareto front** for the selected pair of objectives. Every candidate uses exactly `num_fea` predictors; this number stays fixed during optimization.

---

## 🚀 Key Features

- ⚙️ **NSGA-II Metaheuristic**: Robust multi-objective genetic algorithm.
- 📈 **Optimization Objectives**:
  - `DIST-EPS`: maximize the distance between margin planes and minimize the sum of margin-violation distances (epsilon).
  - `MC`: minimize false negatives (FN) and false positives (FP) separately.
- 🔍 **Visualization**: Built-in `draw_solution()` to inspect the Pareto front.
- 📤 **Data Export**: Save optimal solutions directly to CSV.

Both `MC` objectives use the central classifier and the same predictions as the reported classification metrics. Class `1` is positive, and a zero decision score predicts `1`. In solution summaries and CSV exports, the historical column names `MCPOS` and `MCNEG` contain FN and FP, respectively. The confusion matrix has true classes in rows and predicted classes in columns: FN is `conf_matrix["1", "-1"]`, and FP is `conf_matrix["-1", "1"]`.

# 📦 Installation

To install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("jomabel/SVMFeature_r")
```

## 💻 Basic Usage

The library accepts data frames with numeric predictors and a binary response encoded as `-1` and `1`. This example loads `Colon` from the optional `plsgenomics` package: 62 observations and 2,000 gene-expression predictors. Install the example dependency with `install.packages("plsgenomics")` if necessary. Its original classes are recoded from `1` (normal tissue) to `1`, and from `2` (tumor tissue) to `-1`.

```r
library(svmfsmeta)

# 1. Load Colon from the optional example dependency
if (!requireNamespace("plsgenomics", quietly = TRUE)) {
  stop("Install the example dependency with install.packages('plsgenomics').")
}
data("Colon", package = "plsgenomics")
data <- data.frame(y = ifelse(Colon$Y == 2, -1, 1), Colon$X)
names(data) <- c("y", paste0("x", seq_len(ncol(Colon$X))))

# 2. Normalize input features (Min-Max scaling to [0, 1])
input_cols <- names(data)[-1]
scaler <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep(0.5, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}
data[, input_cols] <- as.data.frame(lapply(data[, input_cols], scaler))

# 3. Initialize the optimizer
set.seed(42)
obj <- svm_fs_meta(
  data = data,
  inputs = names(data)[-1],
  output = "y",
  pop_size = 10,
  num_fea = 5,
  n_iter = 10,
  mode = "iters",
  objective = "DIST-EPS"
)

# 4. Execute optimization
result <- run.svm_fs_meta(obj)

# 5. Visualize and save
draw_solution(result)
save_results(result, output_dir = tempdir(), dataset_name = "colon")
```

The dataset is described by [Alon et al. (1999)](https://doi.org/10.1073/pnas.96.12.6745). See `help("Colon", package = "plsgenomics")` for its source and documentation. Predictor names are assigned as `x1` to `x2000` in this example; the original gene identifiers remain available in `Colon$gene.names`.

---

## 🛠️ API Reference

### `svm_fs_meta` (Constructor)

Initializes the optimization object.

| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `data` | `data.frame` | - | The dataset containing features and target column. |
| `inputs` | `vector` | - | Names of the input feature columns. |
| `output` | `character` | - | Name of the target (class) column. |
| `pop_size` | `numeric` | - | Size of the population. |
| `num_fea` | `numeric` | - | Exact number of predictors selected in every candidate; a positive integer no larger than `length(inputs)`. |
| `objective` | `character` | `"DIST-EPS"` | Objectives: `"DIST-EPS"` or `"MC"`. |
| `n_iter` | `numeric` | `10` | Number of iterations (used if `mode="iters"`). |
| `max_time` | `numeric` | `300` | Approximate optimization time budget in seconds (used if `mode="time"`); initialization is excluded. |
| `mode` | `character` | `"iters"` | Execution mode: `"iters"` or `"time"`. |
| `clones` | `numeric` | `0` | `0` forbids duplicate solutions; `1` allows them. |
| `p_mutation` | `numeric` | `0.7` | Global mutation probability. |
| `p_mut_ind` | `numeric` | `0.4` | Probability of mutating an individual. |
| `p_mut_fea` | `numeric` | `0.4` | Probability of mutating a feature selection. |
| `p_mut_coord` | `numeric` | `0.2` | Probability of mutating hyperplane coordinates. |

| `mut_coord` | `numeric` | `0` | Mutation strategy on coordinates (0: Scaling, 1: Full replacement). |

In time mode, the timer starts after the initial population has been generated, sorted into fronts, and summarized. The budget is checked before each iteration, and an iteration already started runs to completion. Total runtime therefore includes initialization and may exceed `max_time`; this is not a strict timeout.

### `run.svm_fs_meta(object)`
Runs the optimization and returns an updated S3 object of class `svm_fs_meta`. It contains the original data and configuration, plus the final population in `result$population`.

The population contains a summary table, `df_solutions`, and the individual candidate lists, `solution_list`. The table includes all retained fronts; use `FRONT == 1` to extract its non-dominated candidates:

```r
solutions_table <- result$population$df_solutions
pareto_front <- solutions_table[solutions_table$FRONT == 1, , drop = FALSE]

front_solutions <- Filter(
  function(s) s$front == 1,
  result$population$solution_list
)
selected <- front_solutions[[1]]
selected$features      # Selected predictor names
selected$plane_coord   # Normal coefficients, in the same predictor order
selected$plane_term_b  # Intercepts: negative-class margin, positive-class margin, central plane
selected[c("accuracy", "precision", "recall", "f1", "kappa", "auc", "conf_matrix")]
```

Here the first front candidate is inspected. Candidate selection depends on the desired objective trade-off. The reported metrics describe the data used for optimization.

`DOMINATES_TO` and `DOMINATED_BY` contain `SOL` identifiers from the full
population. They are not row positions in the summary table or Pareto-only CSV.
After selection and renumbering, `fnds()` rebuilds these relationships; the
optimizer does this automatically. Internal `fronts` membership uses one-based
positions in `solution_list`.
After sorting, `SOL_DOM_BY` equals the number of identifiers in `DOMINATED_BY`;
the sorting algorithm uses a separate working counter.

### `predict_solution(solution, newdata, type = "class")`

Use `predict_solution()` on an evaluated solution to classify one or more
observations without their true labels:

```r
new_observations <- test[selected$features]  # Already scaled using training parameters
predict_solution(selected, newdata = new_observations, type = "class")
predict_solution(selected, newdata = new_observations, type = "score")
predict_solution(selected, newdata = new_observations[1, , drop = FALSE])
```

Predictor columns are matched by name; extra columns are ignored. The returned
vector keeps row order and uses row names as observation names. Scores are
oriented toward class `1`; a score of zero predicts `1`. Scores are not
probabilities. Reuse preprocessing fitted on training data; `predict_solution()`
does not fit or apply normalization. It accepts evaluated solution lists directly.

### `evaluate_test(solution, data, output)`

Evaluates a previously trained solution on held-out data using its stored
classifier. Returns a separate list with `error_rate`, `n_errors`, `n_samples`,
FN/FP (`fn`, `fp`), accuracy, precision, recall, F1, kappa, AUC,
`confusion_matrix`, `predictions` and `decision_scores`.
The function obtains its decision scores through `predict_solution()`.
The training metrics, objectives and Pareto ranks are unchanged. Class `1` is
positive, including scores of zero. AUC is `NA` for a single true class.
Undefined precision and recall follow the training convention (`NA`/`NaN`).

Split the raw data first, learn normalization only on training, and reuse that
transformation on test. For example, after selecting a trained solution:

```r
test_metrics <- evaluate_test(selected, test, output = "y")
test_metrics$error_rate  # (FN + FP) / number of test observations
test_metrics$confusion_matrix
```

The installed package includes a complete stratified Colon example (41 training
and 21 test observations), with training-only normalization and all training/test
classification metrics for the MC Pareto front:

```r
# install.packages("plsgenomics")
source(system.file("examples", "colon_train_test.R", package = "svmfsmeta"))
```

Keep the final test set independent of model selection. If held-out performance
is used to choose a classifier, use a validation set and a separate final test.

### `draw_solution(x)`
Generates a `ggplot2` scatter plot of the Pareto front (Front 1).

### `save_results(object, output_dir, dataset_name="dataset")`
Exports the Pareto front to a CSV file.

## 👥 Authors
*   **José Manuel Martínez Belando**
*   **Ángel de Lara**
*   **Daniel Valero Carreras**
*   **Antonio Jesús Banegas Luna**

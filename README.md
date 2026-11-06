# 🧠 SVMFsMeta: Multi-objective SVM Classifier with Feature Selection

`SVMFsMeta` is an R package that implements the algorithms proposed by Valero et al. (2023) and Alcaraz, J. (2024). It provides a metaheuristic approach to solve the **soft-margin Support Vector Machine (SVM)** problem combined with **Feature Selection**.

The result of the algorithm is a **Pareto front**, where each solution is an SVM-based classifier, allowing decision-makers to select the optimal model based on their specific trade-offs between accuracy and complexity.

---

## 🚀 Key Features

- ⚙️ **NSGA-II Metaheuristic**: Robust multi-objective genetic algorithm.
- 📈 **Optimization Objectives**:
  - Margin distance and Epsilon (error).
  - False Positives (FP) and False Negatives (FN).
- 🔍 **Visualization**: Built-in `draw_solution()` to inspect the Pareto front.
- 📤 **Data Export**: Save optimal solutions directly to CSV.

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
| `num_fea` | `numeric` | - | Maximum number of features to select. |
| `objective` | `character` | `"DIST-EPS"` | Objectives: `"DIST-EPS"` or `"MC"`. |
| `n_iter` | `numeric` | `10` | Number of iterations (used if `mode="iters"`). |
| `max_time` | `numeric` | `300` | Maximum seconds (used if `mode="time"`). |
| `mode` | `character` | `"iters"` | Execution mode: `"iters"` or `"time"`. |
| `clones` | `numeric` | `0` | `0` forbids duplicate solutions; `1` allows them. |
| `p_mutation` | `numeric` | `0.7` | Global mutation probability. |
| `p_mut_ind` | `numeric` | `0.4` | Probability of mutating an individual. |
| `p_mut_fea` | `numeric` | `0.4` | Probability of mutating a feature selection. |
| `p_mut_coord` | `numeric` | `0.2` | Probability of mutating hyperplane coordinates. |
| `mut_coord` | `numeric` | `0` | Mutation strategy on coordinates (0: Scaling, 1: Full replacement). |

### `run.svm_fs_meta(object)`
Starts the optimization process (NSGA-II).

### `draw_solution(x)`
Generates a `ggplot2` scatter plot of the Pareto front (Front 1).

### `save_results(object, output_dir, dataset_name="dataset")`
Exports the Pareto front to a CSV file.

## 👥 Authors
*   **José Manuel Martínez Belando**
*   **Ángel de Lara**
*   **Daniel Valero Carreras**
*   **Antonio Jesús Banegas Luna**

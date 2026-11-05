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

The library is designed to work with structured bioinformatics datasets (no headers, index column, and class labels in the first column).

```r
library(svmfsmeta)
library(readr)

# 1. Load the dataset (example: 09Colon_p.txt)
dataset_file <- "data/09Colon_p.txt"
data <- as.data.frame(read_delim(dataset_file, delim = ";", col_names = FALSE, show_col_types = FALSE))

# 2. Preprocess: drop index, rename columns, and filter binary labels
data <- data[,-1]
names(data) <- c("y", paste0("x", 1:(ncol(data)-1)))
data <- data[data$y %in% c(-1, 1), ]

# 3. Normalize input features (Min-Max scaling to [0, 1])
input_cols <- names(data)[-1]
scaler <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep(0.5, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}
data[, input_cols] <- as.data.frame(lapply(data[, input_cols], scaler))

# 4. Initialize and run the optimizer
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

# 5. Execute optimization
result <- run.svm_fs_meta(obj)

# 6. Visualize and Save
draw_solution(result)
save_results(result, output_dir = tempdir(), dataset_name = "09Colon_p")
```

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

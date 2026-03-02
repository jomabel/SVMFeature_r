# 🧠 SVMFeature: Multi-objective SVM Classifier with Feature Selection

`SVMFeature` is an R package that implements the algorithms proposed by Valero et al. (2023) and Alcaraz, J. (2024). It provides a metaheuristic approach to solve the **soft-margin Support Vector Machine (SVM)** problem combined with **Feature Selection**.

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
devtools::install_github("ajbanegas/SVMFeature")
```

## 💻 Basic Usage

The library is designed to work with structured bioinformatics datasets (no headers, index column, and class labels in the first column).

```r
library(SVMFeature)
library(readr)

# 1. Load the dataset (example: 09Colon_p.txt)
dataset_file <- "data/09Colon_p.txt"
data <- as.data.frame(read_delim(dataset_file, delim = ";", col_names = FALSE, show_col_types = FALSE))

# 2. Preprocess: drop index, rename columns, extract costs, and filter binary labels
data <- data[,-1]
names(data) <- c("y", paste0("x", 1:(ncol(data)-1)))
costs <- as.numeric(data[1, -1])
data  <- data[-1, ]
data  <- data[data$y %in% c(-1, 1), ]

# 3. Initialize and run the optimizer
obj <- SVMFeature(
  data = data,
  inputs = names(data)[-1],
  output = "y",
  costs = costs,
  pop_size = 10,
  num_fea = 5,
  n_iter = 10,
  mode = "iters",
  objective = "distance-epsilon"
)

# 4. Execute optimization
result <- run.SVMFeature(obj)

# 5. Visualize and Save
draw_solution(result)
save_results(result, dataset_name = "09Colon_p")
```

---

## 🛠️ API Reference

### `SVMFeature` (Constructor)

Initializes the optimization object.

| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `data` | `data.frame` | - | The dataset containing features and target column. |
| `inputs` | `vector` | - | Names of the input feature columns. |
| `output` | `character` | - | Name of the target (class) column. |
| `costs` | `numeric` | - | Costs associated with each feature. |
| `pop_size` | `numeric` | - | Size of the population. |
| `num_fea` | `numeric` | - | Maximum number of features to select. |
| `objective` | `character` | `"distance-epsilon"` | Objectives: `"distance-epsilon"` or `"confusion-matrix"`. |
| `n_iter` | `numeric` | `10` | Number of iterations (used if `mode="iters"`). |
| `max_time` | `numeric` | `300` | Maximum seconds (used if `mode="time"`). |
| `mode` | `character` | `"iters"` | Execution mode: `"iters"` or `"time"`. |
| `clones` | `numeric` | `0` | `0` or `1`. If `1`, checks for duplicates in population. |
| `p_mutation` | `numeric` | `0.7` | Global mutation probability. |
| `p_mut_ind` | `numeric` | `0.4` | Probability of mutating an individual. |
| `p_mut_fea` | `numeric` | `0.4` | Probability of mutating a feature selection. |
| `p_mut_coord` | `numeric` | `0.2` | Probability of mutating hyperplane coordinates. |

### `run.SVMFeature(object)`
Starts the optimization process (NSGA-II).

### `draw_solution(x)`
Generates a `ggplot2` scatter plot of the Pareto front (Front 1).

### `save_results(object, output_dir="results", dataset_name="dataset")`
Exports the Pareto front to a CSV file.

## 👥 Authors
*   **José Manuel Martínez Belando**
*   **Ángel de Lara**
*   **Daniel Valero Carreras**
*   **Antonio Jesús Banegas Luna**

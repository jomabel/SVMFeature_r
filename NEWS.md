# svmfsmeta 0.1.0

* Initial release.
* Provides NSGA-II-based Pareto-front approximation for binary soft-margin SVM
  classification with a fixed number of selected predictors and two alternative
  pairs of objectives.
* Supports `DIST-EPS`, maximizing the distance between margin planes and
  minimizing the sum of margin-violation distances (epsilon).
* Supports `MC`, minimizing false negatives and false positives computed using
  the central classifier. Summary and CSV columns `MCPOS` and `MCNEG` contain
  FN and FP, respectively; the MC plot uses FN on the horizontal axis and FP
  on the vertical axis.
* Computes accuracy, precision, recall, F1, Cohen's kappa, and ROC AUC using the
  data provided for optimization. Returns the final population, candidate
  coefficients, selected features, and a summary table of solutions.
* Provides first-front visualization and CSV export of non-dominated solutions.
* Supports a fixed iteration count or an approximate `max_time` budget. Timing
  excludes initialization and is checked before each iteration; an iteration
  already in progress runs to completion.
* Includes an executable vignette using the `Colon` dataset from the optional
  `plsgenomics` package, with reproducible examples and temporary CSV output.
* Requires R >= 4.1.0, dplyr >= 1.1.0, and ggplot2 >= 3.4.0. Unit tests require
  testthat >= 3.1.7.
* Includes unit tests for input validation, solution evaluation, dominance,
  population generation, crowding distance, and reduction.

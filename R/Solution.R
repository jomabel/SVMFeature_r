  #' @title Definition of the solution Class
#'
#' @description
#'   Defines the `solution` class, which represents an individual solution.
#'
#' @param num solution number
#' @param num_features Number of features
#' @param obj Objective
#'
#' @return Object of the `solution` class
#'
#' @export
#'
#' @examples
#' sol <- solution(num = 1, num_features = 2, obj = "DIST-EPS")
solution <- function(num, num_features, obj) {
  # Define instance variables
  solution <- list(
    num = num,
    num_features = num_features,
    data_sol = NULL,
    features = list(),
    vectors = list(),
    plane_coord = list(),
    plane_term_b = list(),
    obj_fn = obj,
    objective = numeric(4),  # (1)distance, (2)epsilon, (3)FP, (4)FN
    sol_dom_by = 0,  # Initialize domination counter
    dominates_list = list(),  # List of solutions that this solution dominates
    list_dominated_by = list(),  # List of solutions that dominate this solution
    front = -1,
    successful_evaluation = FALSE,
    crowding_distance = 0
  )

  return(solution)
}



#' @title Method to convert the solution to a dictionary
#'
#' @description
#'   Method to convert the solution to a dictionary
#'
#' @param solution solution class object
#'
#' @return List with the information of the solution
#'
#' @export
#'
#' @examples
#' sol <- solution(num = 1, num_features = 2, obj = "DIST-EPS")
#' sol_dict <- to_dict(sol)
to_dict <- function(solution) {
  return(list(
    SOL = if(!is.null(solution$num)) solution$num else NULL,
    VECTORS = if(!is.null(solution$vectors)) solution$vectors else NULL,
    PLANO_COOR = if(!is.null(solution$plane_coord)) solution$plane_coord else NULL,
    PLANO_COOR_B = if(!is.null(solution$plane_term_b)) solution$plane_term_b else NULL,
    FEATURES = if(!is.null(solution$features)) solution$features else NULL,
    DIST = if(!is.null(solution$objective) && length(solution$objective) >= 1) solution$objective[1] else NULL,
    EPS = if(!is.null(solution$objective) && length(solution$objective) >= 2) solution$objective[2] else NULL,
    MCPOS = if(!is.null(solution$objective) && length(solution$objective) >= 3) solution$objective[3] else NULL,
    MCNEG = if(!is.null(solution$objective) && length(solution$objective) >= 4) solution$objective[4] else NULL,
    DOMINATES_TO = if(!is.null(solution$dominates_list)) solution$dominates_list else NULL,
    DOMINATED_BY = if(!is.null(solution$list_dominated_by)) solution$list_dominated_by else NULL,
    SOL_DOM_BY = if(!is.null(solution$sol_dom_by)) solution$sol_dom_by else NULL,
    FRONT = if(!is.null(solution$front)) solution$front else NULL,
    CROW_DIST = if(!is.null(solution$crowding_distance)) solution$crowding_distance else NULL,
    ACCURACY = if(!is.null(solution$accuracy)) solution$accuracy else NULL,
    F1 = if(!is.null(solution$f1)) solution$f1 else NULL,
    PRECISION = if(!is.null(solution$precision)) solution$precision else NULL,
    RECALL = if(!is.null(solution$recall)) solution$recall else NULL,
    KAPPA = if(!is.null(solution$kappa)) solution$kappa else NULL,
    AUC = if(!is.null(solution$auc)) solution$auc else NULL
  ))
}



#' @title Method to Obtain a Random Class Vector
#'
#' @description
#'   Method to obtain a random class vector from the dataset.
#'
#' @param solution solution class object containing the dataset
#' @param clazz Class for which to obtain a vector
#' @param data A data.frame with the input features
#' @param output A vector with class labels
#'
#' @return Index of the random class vector
#' @export
#'
#' @examples
#' data <- data.frame(y = c(-1, 1, 1, -1))
#' sol <- solution(num = 1, num_features = 2, obj = "DIST-EPS")
#' idx <- get_class_vector(sol, clazz = 1, data = data, output = "y")
get_class_vector <- function(solution, clazz, data, output) {
  indices <- which(data[[output]] == clazz)
  if (length(indices) > 0) {
    return(sample(indices, 1))
  } else {
    return(NA)  # NA is a numeric value in R, indicating 'not available'
  }
}



#' @title Method to Generate a Random solution
#'
#' @description
#'   Method to generate a random solution by selecting random features and class vectors.
#'
#' @param solution solution class object to modify
#' @param data A data.frame with the input features
#' @param inputs A vector of input features
#' @param output A vector with class labels
#'
#' @importFrom stats runif
#'
#' @return solution class object with a randomly generated solution
#'
#' @export
#'
#' @examples
#' data <- data.frame(x1 = rnorm(4), x2 = rnorm(4), y = c(-1, 1, 1, -1))
#' sol <- solution(num = 1, num_features = 2, obj = "DIST-EPS")
#' sol <- generate_random_solution(sol, data, c("x1", "x2"), "y")
generate_random_solution <- function(solution, data, inputs, output) {
  # Select two points from class A and B
  classes <- sort(unique(data[[output]]))
  solution$vectors <- vector("list", length(classes))

  for (i in seq_along(classes)) {
    solution$vectors[[i]] <- get_class_vector(solution, classes[i], data, output)
  }

  # Seleccionar aleatoriamente p características y p coordenadas entre [-1, 1]
  if (solution$num_features == length(inputs)) {
    solution$features <- inputs
  } else {
    # Verify that the inputs exist as columns
    possible_features <- inputs[which(inputs %in% names(data))]
    solution$features <- sort(
      sample(possible_features, solution$num_features)
    )
  }

  # Calculate plane coordinates
  # Random coordinates for each feature
  solution$plane_coord <- runif(solution$num_features, -1, 1)

  # Adjust solution dataframe
  indices <- unlist(solution$vectors)
  if (length(indices) == 0 || length(solution$features) == 0) {
    # Returns the unmodified solution if no valid indices or selected features
    return(solution)
  }

  solution$data_sol <- data[indices, solution$features, drop = FALSE]

  return(solution)
}



#' @title Method to Construct Planes
#'
#' @description
#'   Method to construct planes by calculating the coefficients based on features.
#'
#' @param solution solution class object containing data and features to construct planes
#' @param data A data.frame with the input features
#'
#' @return solution class object with constructed planes
#'
#' @export
#'
#' @examples
#' data <- data.frame(x1 = rnorm(4), x2 = rnorm(4), y = c(-1, 1, 1, -1))
#' sol <- solution(num = 1, num_features = 2, obj = "DIST-EPS")
#' sol <- generate_random_solution(sol, data, c("x1", "x2"), "y")
#' sol <- construct_planes(sol, data)
construct_planes <- function(solution, data) {
  # Prepare data_sol dataframe by selecting only necessary rows and columns
  valid_indices <- unlist(solution$vectors)
  valid_indices <- valid_indices[valid_indices <= nrow(data) & valid_indices > 0]

  if (length(valid_indices) == 0) {
    stop("No valid index in 'vectors'.")
  }

  # Filter features to ensure they are present in the data
  valid_features <- solution$features[solution$features %in% names(data)]
  if (length(valid_features) == 0) {
    stop("No valid feature for selection.")
  }

  # Subselect the dataframe
  solution$data_sol <- data[valid_indices, valid_features, drop = FALSE]

  # Initialize the plane terms vector
  solution$plane_term_b <- numeric()

  # Calculate the independent terms of the plane equation
  for (i in seq_len(nrow(solution$data_sol))) {
    independent_term <- 0
    for (j in seq_along(valid_features)) {
      independent_term <- independent_term + solution$data_sol[i, j] * solution$plane_coord[j]
    }
    # Add the negative value of 'independent_term' to the plane terms vector
    solution$plane_term_b <- c(solution$plane_term_b, -independent_term)
  }

  # Add the intermediate plane term
  if(length(solution$plane_term_b) > 0) {
    middle <- mean(solution$plane_term_b)
    solution$plane_term_b <- c(solution$plane_term_b, middle)
  } else {
    stop("No plane terms have been calculated, so the mean can not be calculated.")
  }

  return(solution)
}



#' @title Method to Calculate the Distance Objective
#'
#' @description
#'   This method calculates the distance objective based on the plane equation.
#'
#' @param solution solution class object containing the plane equation
#'
#' @return solution class object with the distance objective calculated
#'
#' @export
#'
#' @examples
#' data <- data.frame(x1 = rnorm(4), x2 = rnorm(4), y = c(-1, 1, 1, -1))
#' sol <- solution(num = 1, num_features = 2, obj = "DIST-EPS")
#' sol <- generate_random_solution(sol, data, c("x1", "x2"), "y")
#' sol <- construct_planes(sol, data)
#' sol <- calculate_distance_objective(sol)
calculate_distance_objective <- function(solution) {
  # Initialize the distance objective to zero
  solution$objective[[1]] <- 0
  denominator <- sum(solution$plane_coord^2) ^ 0.5
  distance <- abs(solution$plane_term_b[[1]] - solution$plane_term_b[[2]])

  # Calculate normalized distance if the denominator is not zero
  if (denominator != 0) {
    solution$objective[[1]] <- distance / denominator
  } else {
    solution$objective[[1]] <- -1  # Error value or non-calculable indicator
  }
  return(solution)
}



#' @title Method to Calculate the Epsilon Objective
#'
#' @description
#'   This method calculates the epsilon objective by evaluating misclassified points.
#'
#' @param solution solution class object containing data and features for classification
#' @param data A data.frame with the input features.
#' @param output A vector with class labels.
#'
#' @return solution class object with the epsilon objective calculated
#'
#' @export
#'
#' @examples
#' data <- data.frame(x1 = rnorm(4), x2 = rnorm(4), y = c(-1, 1, 1, -1))
#' sol <- solution(num = 1, num_features = 2, obj = "DIST-EPS")
#' sol <- generate_random_solution(sol, data, c("x1", "x2"), "y")
#' sol <- construct_planes(sol, data)
#' sol <- calculate_epsilon_objective(sol, data, "y")
calculate_epsilon_objective <- function(solution, data, output) {
  # Initialize misclassification counters
  solution$objective[[2]] <- 0
  solution$objective[[3]] <- 0
  solution$objective[[4]] <- 0

  # Initialize necessary lists to calculate metrics
  y_true <- c()
  y_pred <- c()
  decision_scores <- c()

  # Support vectors are ordered by class: -1, then 1. Orient scores so
  # positive values point toward class 1 for either normal direction.
  orientation <- if (solution$plane_term_b[[1]] > solution$plane_term_b[[2]]) {
    1
  } else {
    -1
  }

  # Iterar sobre todos los datos para calcular la suma de errores
  for (i in 1:nrow(data)) {
    clazz <- data[i, output]
    y_true <- c(y_true, clazz) # Save actual class
    distance <- 0
    denominator <- 0

    # Calculate the distance of the vector to the planes d = (wx + b)/||w||
    for (index in seq_along(solution$features)) {
      fea <- solution$features[[index]]
      distance <- distance + solution$plane_coord[index] * data[i, fea]
      denominator <- denominator + (solution$plane_coord[index] ** 2)
    }

    # Use the same class-oriented score for predictions and AUC.
    f_x <- orientation * (distance + solution$plane_term_b[[3]])
    if (f_x >= 0) {
      pred <- 1
    } else {
      pred <- -1
    }
    y_pred <- c(y_pred, pred)
    decision_scores <- c(decision_scores, f_x)

    # Adjust for the independent term and count errors
    # If the class is 1, distance is calculated with plane b[2]
    if (clazz == 1) {
      distance <- distance + solution$plane_term_b[[2]]
    } else {
      distance <- distance + solution$plane_term_b[[1]]
    }

    distance <- distance / denominator ** (1/2)

    # If b[0] > b[1], class -1 is above the plane and class 1 below
    if (solution$plane_term_b[[1]] > solution$plane_term_b[[2]]) {
      # If class is -1 and distance is positive, it is misclassified
      if (clazz == -1) {
        if (distance > 0) {
          solution$objective[[4]] <- solution$objective[[4]] + 1
          solution$objective[[2]] <- solution$objective[[2]] + distance
        }
      }
      # If class is 1 and distance is negative, it is misclassified
      else {
        if (distance < 0) {
          solution$objective[[3]] <- solution$objective[[3]] + 1
          solution$objective[[2]] <- solution$objective[[2]] + abs(distance)
        }
      }
    }
    # If b[0] < b[1], class -1 is below the plane and class 1 above
    else {
      # If class is -1 and distance is negative, it is misclassified
      if (clazz == -1) {
        if (distance < 0) {
          solution$objective[[4]] <- solution$objective[[4]] + 1
          solution$objective[[2]] <- solution$objective[[2]] + abs(distance)
        }
      }
      # If class is 1 and distance is positive, it is misclassified
      else {
        if (distance > 0) {
          solution$objective[[3]] <- solution$objective[[3]] + 1
          solution$objective[[2]] <- solution$objective[[2]] + distance
        }
      }
    }
  }
  solution <- svmfsmeta::calculate_metrics(solution, y_true, y_pred, decision_scores)
  return(solution)
}


#' @title Calculate Classification Metrics
#'
#' @description
#'   This function calculates common classification metrics including accuracy, precision,
#'   recall, F1-score, Cohen's kappa, AUC, and the confusion matrix. The results are stored
#'   as fields in the provided solution object.
#'
#' @param solution A solution class object where the metrics will be stored
#' @param y_true A vector of true class labels (-1 or 1)
#' @param y_pred A vector of predicted class labels (-1 or 1)
#' @param decision_scores A numeric vector of decision scores (f(x) values)
#'
#' @return The input solution object with the following additional fields:
#'   \item{conf_matrix}{Confusion matrix}
#'   \item{accuracy}{Classification accuracy}
#'   \item{precision}{Precision}
#'   \item{recall}{Recall}
#'   \item{f1}{F1 score}
#'   \item{kappa}{Cohen's kappa statistic}
#'   \item{auc}{Area Under the ROC Curve (AUC)}
#'
#' @export
#'
#' @examples
#' sol <- solution(num = 1, num_features = 2, obj = "DIST-EPS")
#' y_true <- c(-1, 1, 1, -1)
#' y_pred <- c(-1, 1, -1, -1)
#' scores <- c(-0.5, 0.8, -0.1, -0.9)
#' sol <- svmfsmeta::calculate_metrics(sol, y_true, y_pred, scores)
calculate_metrics <- function(solution, y_true, y_pred, decision_scores){
  y_true_bin <- ifelse(y_true == 1, 1, 0)
  y_pred_bin <- ifelse(y_pred == 1, 1, 0)
  # Confusion matrix
  solution$conf_matrix <- table(factor(y_true, levels=c(-1,1)),
                                factor(y_pred, levels=c(-1,1)))

  # Accuracy
  solution$accuracy <- Metrics::accuracy(y_true_bin, y_pred_bin)
  # Precision, Recall, F1 (for positive class = 1)
  solution$precision <- Metrics::precision(y_true_bin, y_pred_bin)
  solution$recall <- Metrics::recall(y_true_bin, y_pred_bin)
  f1_denominator <- solution$precision + solution$recall
  solution$f1 <- if (is.na(f1_denominator) || f1_denominator == 0) {
    0
  } else {
    2 * solution$precision * solution$recall / f1_denominator
  }

  # Cohen Kappa
  solution$kappa <- psych::cohen.kappa(solution$conf_matrix)$kappa

  # AUC
  solution$auc <- Metrics::auc(y_true_bin, decision_scores)

  return(solution)
}


#' @title Evaluate solution
#'
#' @description
#'   This function evaluates a solution by calculating its objectives.
#'
#' @param solution solution class object to be evaluated
#' @param data A data.frame with the input features
#' @param output A vector with class labels
#'
#' @return solution class object with updated evaluation status and objectives
#'
#' @export
#'
#' @examples
#' data <- data.frame(x1 = rnorm(4), x2 = rnorm(4), y = c(-1, 1, 1, -1))
#' sol <- solution(num = 1, num_features = 2, obj = "DIST-EPS")
#' sol <- generate_random_solution(sol, data, c("x1", "x2"), "y")
#' sol <- evaluate_solution(sol, data, "y")
evaluate_solution <- function(solution, data, output) {
  # Verify if the sum of plane_coord is zero
  if (sum(solution$plane_coord) == 0) {
    message("The solution cannot be evaluated, there are no coordinates.")
    solution$successful_evaluation <- FALSE
    return(solution)
  } else {
    # Procedures to evaluate the solution if there are valid coordinates
    solution <- construct_planes(solution, data)

    solution <- svmfsmeta::calculate_distance_objective(solution)
    solution <- svmfsmeta::calculate_epsilon_objective(solution, data, output)


    # Assume evaluations were successful if this point is reached
    solution$successful_evaluation <- TRUE
    return(solution)
  }
}



#' @title Method to Determine Dominance
#'
#' @description
#'   This function determines if one solution dominates another based on objectives.
#'
#' @param solution1 First solution object
#' @param solution2 Second solution object
#'
#' @return TRUE if solution1 dominates solution2, otherwise FALSE
#'
#' @export
#'
#' @examples
#' sol1 <- solution(num = 1, num_features = 2, obj = "DIST-EPS")
#' sol1$objective <- c(1, 0, 0, 0)
#' sol2 <- solution(num = 2, num_features = 2, obj = "DIST-EPS")
#' sol2$objective <- c(0, 1, 1, 1)
#' dominate(sol1, sol2)
dominate <- function(solution1, solution2) {
  if (!is_valid_solution(solution1) || !is_valid_solution(solution2)) {
    # Return FALSE or handle the case as per your application logic

    return(FALSE)
  }

  # Proceed with comparison if all values are suitable
  if (solution1$obj_fn == 'DIST-EPS') {
    return(solution1$objective[[1]] >= solution2$objective[[1]] &&
             solution1$objective[[2]] <= solution2$objective[[2]])

  } else if (solution1$obj_fn == 'MC') {
    return(solution1$objective[[3]] <= solution2$objective[[3]] &&
             solution1$objective[[4]] <= solution2$objective[[4]])
  }else{
    return(FALSE)
  }
}


#' @title Method to assess if a solution is valid
#'
#' @description
#'   This function determines if the objectives of a given solution
#'   contain numerical values.
#'
#' @param solution A solution object
#'
#' @return TRUE if solution is valid, otherwise FALSE
#'
#' @export
#'
#' @examples
#' sol <- solution(num = 1, num_features = 2, obj = "DIST-EPS")
#' sol$objective <- c(1, 0, 0, 0)
#' is_valid_solution(sol)
is_valid_solution <- function(solution) {
  if (solution$obj_fn == 'DIST-EPS') {
    if (is.na(solution$objective[[1]]) || !is.numeric(solution$objective[[1]]) ||
        is.na(solution$objective[[2]]) || !is.numeric(solution$objective[[2]])) {
      return(FALSE)
    }
  } else if (solution$obj_fn == 'MC') {
    if (is.na(solution$objective[[3]]) || !is.numeric(solution$objective[[3]]) ||
        is.na(solution$objective[[4]]) || !is.numeric(solution$objective[[4]])) {
      return(FALSE)
    }
  }

  return(TRUE)
}




#' @title Method to Determine Dominance with Three States
#'
#' @description
#'   This function determines the dominance state between two solutions.
#'
#' @param solution1 First solution object
#' @param solution2 Second solution object
#'
#' @return 0 if the solutions are equal, 1 if solution1 dominates solution2, 2 if solution2 dominates solution1
#'
#' @export
#'
#' @examples
#' sol1 <- solution(num = 1, num_features = 2, obj = "DIST-EPS")
#' sol1$objective <- c(1, 0, 0, 0)
#' sol2 <- solution(num = 2, num_features = 2, obj = "DIST-EPS")
#' sol2$objective <- c(0, 1, 1, 1)
#' dominate2(sol1, sol2)
dominate2 <- function(solution1, solution2) {
  if (solution1$obj_fn == 'DIST-EPS') {

    # Check if objectives are equal
    if (solution1$objective[[1]] == solution2$objective[[1]] &&
        solution1$objective[[2]] == solution2$objective[[2]]) {
      return(0)  # No dominance, they are equal
    } else if (solution1$objective[[1]] >= solution2$objective[[1]] &&
               solution1$objective[[2]] <= solution2$objective[[2]]) {
      return(1)  # solution1 dominates solution2
    } else {
      return(2)  # solution2 dominates solution1
    }
  } else if (solution1$obj_fn == 'MC') {
    # Check if objectives are equal
    if (solution1$objective[[3]] == solution2$objective[[3]] &&
        solution1$objective[[4]] == solution2$objective[[4]]) {
      return(0)  # No dominance, they are equal
    } else if (solution1$objective[[3]] <= solution2$objective[[3]] &&
               solution1$objective[[4]] <= solution2$objective[[4]]) {
      return(1)  # solution1 dominates solution2
    } else {
      return(2)  # solution2 dominates solution1
    }
  }
}


#' @title Method to Compare solutions
#'
#' @description
#'   This function compares two solutions based on their front and dominance information.
#'
#' @param solution1 First solution object
#' @param solution2 Second solution object
#'
#' @return The solution object with better dominance or front
#'
#' @export
#'
#' @examples
#' sol1 <- solution(num = 1, num_features = 2, obj = "DIST-EPS")
#' sol1$front <- 1
#' sol2 <- solution(num = 2, num_features = 2, obj = "DIST-EPS")
#' sol2$front <- 2
#' compare_solutions(sol1, sol2)
compare_solutions <- function(solution1, solution2) {
  if (solution1$front == solution2$front) {
    if (length(solution1$list_dominated_by) == length(solution2$list_dominated_by)) {
      coin_flip <- sample(0:1, 1)
      if (coin_flip == 0) {
        return(solution1)
      } else {
        return(solution2)
      }
    } else if (length(solution1$list_dominated_by) < length(solution2$list_dominated_by)) {
      return(solution1)
    } else {
      return(solution2)
    }
  } else {
    if (solution1$front > solution2$front) {
      return(solution2)
    } else {
      return(solution1)
    }
  }
}



#' @title Method to Mutate Vectors
#'
#' @description
#'   This function mutates the vectors in the solution by assigning new class vectors.
#'
#' @param solution solution class object to mutate
#' @param data A data.frame with the input features.
#' @param output A vector with class labels.
#'
#' @return solution class object with updated vectors
#'
#' @export
#'
#' @examples
#' data <- data.frame(x1 = rnorm(4), x2 = rnorm(4), y = c(-1, 1, 1, -1))
#' sol <- solution(num = 1, num_features = 2, obj = "DIST-EPS")
#' sol <- generate_random_solution(sol, data, c("x1", "x2"), "y")
#' sol <- mutate_vectors(sol, data, "y")
mutate_vectors <- function(solution, data, output) {
  # Get unique classes and sort them
  classes <- sort(unique(data[[output]]))

  # Iterate over each class
  for (i in seq_along(classes)) {
    # Call get_class_vector for the current class
    vector <- get_class_vector(solution, classes[i], data, output)

    # Verify that the new vector is different from the current one; if not, find another
    while (identical(vector, solution$vectors[[i]])) {
      vector <- get_class_vector(solution, classes[i], data, output)
    }

    # Update the vector in the solution
    solution$vectors[[i]] <- vector
  }

  return(solution)
}

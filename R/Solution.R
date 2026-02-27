  library(caret)
  library(Metrics)
  library(psych) #Para cohen.kappa

#' @title Definition of the Solution Class
#'
#' @description
#'   Defines the `Solution` class, which represents an individual solution.
#'
#' @param num Solution number
#' @param num_features Number of features
#' @param obj Objective
#'
#' @return Object of the `Solution` class
#'
#' @export
Solution <- function(num, num_features, obj) {
  # Definir variables de instancia
  solution <- list(
    num = num,
    num_features = num_features,
    data_sol = NULL,
    features = list(),
    vectors = list(),
    plane_coord = list(),
    plane_term_b = list(),
    obj_fn = obj,
    objective = numeric(4),  # (1)distancia, (2)epsilon, (3)FP, (4)FN
    sol_dom_by = 0,  # Inicializar contador de dominación
    dominates_list = list(),  # Lista de soluciones que esta solución domina
    list_dominated_by = list(),  # Lista de soluciones que dominan a esta solución
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
#' @param solution Solution class object
#'
#' @return List with the information of the solution
#'
#' @export
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
    # COSTE = if(!is.null(solution$obj_coste)) solution$obj_coste else NULL,  # Descomentar si se usa
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
#' @param solution Solution class object containing the dataset
#' @param clazz Class for which to obtain a vector
#' @param data A data.frame with the input features
#' @param output A vector with class labels
#'
#' @return Index of the random class vector
#' @export
get_class_vector <- function(solution, clazz, data, output) {
  indices <- which(data[[output]] == clazz)
  if (length(indices) > 0) {
    return(sample(indices, 1))
  } else {
    return(NA)  # NA es un valor numérico en R, que indica 'no disponible'
  }
}



#' @title Method to Generate a Random Solution
#'
#' @description
#'   Method to generate a random solution by selecting random features and class vectors.
#'
#' @param solution Solution class object to modify
#' @param data A data.frame with the input features
#' @param inputs A vector of input features
#' @param output A vector with class labels
#'
#' @importFrom stats runif
#'
#' @return Solution class object with a randomly generated solution
#'
#' @export
generate_random_solution <- function(solution, data, inputs, output) {
  # Seleccionar dos puntos de la clase A y B
  classes <- sort(unique(data[[output]]))
  solution$vectors <- vector("list", length(classes))

  for (i in seq_along(classes)) {
    solution$vectors[[i]] <- get_class_vector(solution, classes[i], data, output)
  }

  # Seleccionar aleatoriamente p características y p coordenadas entre [-1, 1]
  if (solution$num_features == length(inputs)) {
    solution$features <- inputs
  } else {
    # Se verifica que las entradas existan como columnas
    possible_features <- inputs[which(inputs %in% names(data))]
    solution$features <- sort(
      sample(possible_features, solution$num_features)
    )
  }

  # Calcular las coordenadas del plano
  # Coordenadas aleatorias para cada característica
  solution$plane_coord <- runif(solution$num_features, -1, 1)

  # Ajustar dataframe de solución
  indices <- unlist(solution$vectors)
  if (length(indices) == 0 || length(solution$features) == 0) {
    # Retorna la solución sin modificar si no hay índices válidos o características seleccionadas
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
#' @param solution Solution class object containing data and features to construct planes
#' @param data A data.frame with the input features
#'
#' @return Solution class object with constructed planes
#'
#' @export
construct_planes <- function(solution, data) {
  # Preparar el dataframe data_sol seleccionando solo las filas y columnas necesarias
  valid_indices <- unlist(solution$vectors)
  valid_indices <- valid_indices[valid_indices <= nrow(data) & valid_indices > 0]

  if (length(valid_indices) == 0) {
    stop("No valid index in 'vectors'.")
  }

  # Filtrar las características para asegurarse de que están presentes en los datos
  valid_features <- solution$features[solution$features %in% names(data)]
  if (length(valid_features) == 0) {
    stop("No valid feature for selection.")
  }

  # Subseleccionar el dataframe
  solution$data_sol <- data[valid_indices, valid_features, drop = FALSE]

  # Inicializar el vector de términos del plano
  solution$plane_term_b <- numeric()

  # Calcular los términos independientes de la ecuación del plano
  for (i in seq_len(nrow(solution$data_sol))) {
    independent_term <- 0
    for (j in seq_along(valid_features)) {
      independent_term <- independent_term + solution$data_sol[i, j] * solution$plane_coord[j]
    }
    # Añadir el valor negativo de 'independent_term' al vector de términos del plano
    solution$plane_term_b <- c(solution$plane_term_b, -independent_term)
  }

  # Añadir el término del plano intermedio
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
#' @param solution Solution class object containing the plane equation
#'
#' @return Solution class object with the distance objective calculated
#'
#' @export
calculate_distance_objective <- function(solution) {
  # Inicializar el objective de distancia a cero
  solution$objective[[1]] <- 0
  denominator <- sum(solution$plane_coord^2) ^ 0.5
  distance <- abs(solution$plane_term_b[[1]] - solution$plane_term_b[[2]])

  # Calcular la distancia normalizada si el denominador no es cero
  if (denominator != 0) {
    solution$objective[[1]] <- distance / denominator
  } else {
    solution$objective[[1]] <- -1  # Valor de error o indicativo de no calculable
  }
  return(solution)
}



#' @title Method to Calculate the Epsilon Objective
#'
#' @description
#'   This method calculates the epsilon objective by evaluating misclassified points.
#'
#' @param solution Solution class object containing data and features for classification
#' @param data A data.frame with the input features.
#' @param output A vector with class labels.
#'
#' @return Solution class object with the epsilon objective calculated
#'
#' @export
calculate_epsilon_objective <- function(solution, data, output) {
  # Inicializar contadores de clasificaciones malas
  solution$objective[[2]] <- 0
  solution$objective[[3]] <- 0
  solution$objective[[4]] <- 0

  # Inicializar listas necesarias para calcular métricas
  y_true <- c()
  y_pred <- c()
  decision_scores <- c()


  # Iterar sobre todos los datos para calcular la suma de errores
  for (i in 1:nrow(data)) {
    clazz <- data[i, output]
    y_true <- c(y_true, clazz) # Guardar clase real
    distance <- 0
    denominator <- 0

    # Calculamos la distancia del vector a los planos d= wx + b/||w||
    for (index in seq_along(solution$features)) {
      fea <- solution$features[[index]]
      distance <- distance + solution$plane_coord[index] * data[i, fea]
      denominator <- denominator + (solution$plane_coord[index] ** 2)
    }

    # Función de decisión: f(x) = w·x + b
    f_x <- distance + solution$plane_term_b[3]
    if (f_x >= 0) {
      pred <- 1
    } else {
      pred <- -1
    }
    y_pred <- c(y_pred, pred)
    decision_scores <- c(decision_scores, f_x)

    # Ajustar por el término independiente y contar errores
    # Si la clase es 1, la distancia se calcula con el plano b[1]
    if (clazz == 1) {
      distance <- distance + solution$plane_term_b[[2]]
    } else {
      distance <- distance + solution$plane_term_b[[1]]
    }

    distance <- distance / denominator ** (1/2)

    # Si b[0] > b[1], la clase -1 está por encima del plano y la clase 1 por debajo
    if (solution$plane_term_b[[1]] > solution$plane_term_b[[2]]) {
      # Si la clase es -1 y la distancia es positiva está mal clasificada
      if (clazz == -1) {
        if (distance > 0) {
          solution$objective[[4]] <- solution$objective[[4]] + 1
          solution$objective[[2]] <- solution$objective[[2]] + distance
        }
      }
      # Si la clase es 1 y la distancia es negativa está mal clasificada
      else {
        if (distance < 0) {
          solution$objective[[3]] <- solution$objective[[3]] + 1
          solution$objective[[2]] <- solution$objective[[2]] + abs(distance)
        }
      }
    }
    # Si b[0] < b[1], la clase -1 está por debajo del plano y la clase 1 por encima
    else {
      # Si la clase es -1 y la distancia es negativa está mal clasificada
      if (clazz == -1) {
        if (distance < 0) {
          solution$objective[[4]] <- solution$objective[[4]] + 1
          solution$objective[[2]] <- solution$objective[[2]] + abs(distance)
        }
      }
      # Si la clase es 1 y la distancia es positiva está mal clasificada
      else {
        if (distance > 0) {
          solution$objective[[3]] <- solution$objective[[3]] + 1
          solution$objective[[2]] <- solution$objective[[2]] + distance
        }
      }
    }
  }
  solution <- calculate_metrics(solution, y_true, y_pred, decision_scores)
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
calculate_metrics <- function(solution, y_true, y_pred, decision_scores){
  y_true_bin <- ifelse(y_true == 1, 1, 0)
  y_pred_bin <- ifelse(y_pred == 1, 1, 0)
  # Matriz de confusión
  solution$conf_matrix <- table(factor(y_true, levels=c(-1,1)),
                                factor(y_pred, levels=c(-1,1)))

  # Accuracy
  solution$accuracy <- accuracy(y_true_bin, y_pred_bin)
  # Precisión, Recall, F1 (para la clase positiva = 1)
  solution$precision <- Metrics::precision(y_true_bin, y_pred_bin)
  solution$recall <- Metrics::recall(y_true_bin, y_pred_bin)
  solution$f1 <- ifelse(is.na(solution$precision) | is.na(solution$recall), 0, 2 * solution$precision * solution$recall / (solution$precision + solution$recall))

  # Cohen Kappa
  solution$kappa <- cohen.kappa(solution$conf_matrix)$kappa

  # AUC
  solution$auc <- auc(y_true_bin, decision_scores)

  return(solution)
}


#' @title Method to Calculate the FP/FN Objective
#'
#' @description
#'   This method calculates the FP/FN objective by evaluating misclassified points.
#'
#' @param solution Solution class object containing data and features for classification
#' @param data A data.frame with the input features
#' @param output A vector with class labels
#'
#' @return Solution class object with the FP/FN objective calculated
#'
#' @export
calculate_misclassified_objective <- function(solution, data, output) {
  # Inicializar contadores de clasificaciones malas
  solution$objective[[3]] <- 0
  solution$objective[[4]] <- 0

  # Iterar sobre todos los datos para calcular la suma de errores
  for (i in 1:nrow(data)) {
    clazz <- data[i, output]
    distance <- 0
    denominator <- 0

    # Calculamos la distancia del vector a los planos d= wx + b/||w||
    for (index in seq_along(solution$features)) {
      fea <- solution$features[[index]]
      distance <- distance + solution$plane_coord[index] * data[i, fea]
      denominator <- denominator + (solution$plane_coord[index] ** 2)
    }

    # Ajustar por el término independiente y contar errores
    # Si la clase es 1, la distancia se calcula con el plano b[0]
    if (clazz == 1) {
      distance <- distance + solution$plane_term_b[[1]]
    } else {
      distance <- distance + solution$plane_term_b[[2]]
    }

    distance <- distance / denominator ** (1/2)

    # Si b[0] < b[1], la clase 1 está por encima del plano y la clase -1 por debajo
    if (solution$plane_term_b[[1]] < solution$plane_term_b[[2]]) {
      # Si la clase es -1 y la distancia es positiva está mal clasificada
      if (clazz == -1) {
        if (distance > 0) {
          solution$objective[[4]] <- solution$objective[[4]] + 1
        }
      }
      # Si la clase es 1 y la distancia es negativa está mal clasificada
      else {
        if (distance < 0) {
          solution$objective[[3]] <- solution$objective[[3]] + 1
        }
      }
    }
    # Si b[0] > b[1], la clase 1 está por encima del plano y la clase -1 por debajo
    else {
      # Si la clase es -1 y la distancia es negativa está mal clasificada
      if (clazz == -1) {
        if (distance < 0) {
          solution$objective[[4]] <- solution$objective[[4]] + 1
        }
      }
      # Si la clase es 1 y la distancia es positiva está mal clasificada
      else {
        if (distance > 0) {
          solution$objective[[3]] <- solution$objective[[3]] + 1
        }
      }
    }
  }

  return(solution)
}

# calcular_objetivo_costes <- function(solucion) {
#   solucion$obj_coste <- 0
#   for (feature in solucion$features) {
#     solucion$obj_coste <- solucion$obj_coste + solucion$costs[feature]
#   }
#   return(solucion)
# }



#' @title Evaluate Solution
#'
#' @description
#'   This function evaluates a solution by calculating its objectives.
#'
#' @param solution Solution class object to be evaluated
#' @param data A data.frame with the input features
#' @param output A vector with class labels
#'
#' @return Solution class object with updated evaluation status and objectives
#'
#' @export
evaluate_solution <- function(solution, data, output) {
  # Verifica si la suma de plane_coord es cero
  if (sum(solution$plane_coord) == 0) {
    cat("The solution cannot be evaluated, there are no coordinates.\n")
    solution$successful_evaluation <- FALSE
    return(solution)
  } else {
    # Procedimientos para evaluar la solución si hay coordenadas válidas
    solution <- construct_planes(solution, data)

    solution <- calculate_distance_objective(solution)
    solution <- calculate_epsilon_objective(solution, data, output)
    #solution <- calculate_misclassified_objective(solution, data, output)


    # Asumir que las evaluaciones fueron exitosas si se alcanza este punto
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
#' @return TRUE if solucion1 dominates solucion2, otherwise FALSE
#'
#' @export
dominate <- function(solution1, solution2) {
  if (!is_valid_solution(solution1) || !is_valid_solution(solution2)) {
    # Retorna FALSE o maneja el caso de forma que se ajuste a tu lógica de aplicación

    return(FALSE)
  }

  # Procede con la comparación si todos los valores son adecuados
  if (solution1$obj_fn == 'distance-epsilon') {
    return(solution1$objective[[1]] >= solution2$objective[[1]] &&
             solution1$objective[[2]] <= solution2$objective[[2]])

  } else if (solution1$obj_fn == 'confusion-matrix') {
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
is_valid_solution <- function(solution) {
  if (solution$obj_fn == 'distance-epsilon') {
    if (is.na(solution$objective[[1]]) || !is.numeric(solution$objective[[1]]) ||
        is.na(solution$objective[[2]]) || !is.numeric(solution$objective[[2]])) {
      return(FALSE)
    }
  } else if (solution$obj_fn == 'confusion-matrix') {
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
#' @return 0 if the solutions are equal, 1 if solucion1 dominates solucion2, 2 if solucion2 dominates solucion1
#'
#' @export
dominate2 <- function(solution1, solution2) {
  if (solution1$obj_fn == 'distance-epsilon') {

    # Comprobar si los objetivos son iguales
    if (solution1$objective[[1]] == solution2$objective[[1]] &&
        solution1$objective[[2]] == solution2$objective[[2]]) {
      return(0)  # No hay dominancia, son iguales
    } else if (solution1$objective[[1]] >= solution2$objective[[1]] &&
               solution1$objective[[2]] <= solution2$objective[[2]]) {
      return(1)  # solucion1 domina a solucion2
    } else {
      return(2)  # solucion2 domina a solucion1
    }
  } else if (solution1$obj_fn == 'confusion-matrix') {
    # Comprobar si los objetivos son iguales
    if (solution1$objective[[3]] == solution2$objective[[3]] &&
        solution1$objective[[4]] == solution2$objective[[4]]) {
      return(0)  # No hay dominancia, son iguales
    } else if (solution1$objective[[3]] <= solution2$objective[[3]] &&
               solution1$objective[[4]] <= solution2$objective[[4]]) {
      return(1)  # solucion1 domina a solucion2
    } else {
      return(2)  # solucion2 domina a solucion1
    }
  }
}


#' @title Method to Compare Solutions
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
#' @param solution Solution class object to mutate
#' @param data A data.frame with the input features.
#' @param output A vector with class labels.
#'
#' @return Solution class object with updated vectors
#'
#' @export
mutate_vectors <- function(solution, data, output) {
  # Obtener las clases únicas y ordenarlas
  classes <- sort(unique(data[[output]]))

  # Iterar sobre cada clase
  for (i in seq_along(classes)) {
    # Llamar a get_class_vector para la clase actual
    vector <- get_class_vector(solution, classes[i], data, output)

    # Verificar que el nuevo vector sea diferente al actual y si no, buscar otro
    while (identical(vector, solution$vectors[[i]])) {
      vector <- get_class_vector(solution, classes[i], data, output)
    }

    # Actualizar el vector en la solución
    solution$vectors[[i]] <- vector
  }

  return(solution)
}

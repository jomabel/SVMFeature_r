#' @title Definition of the Population Class
#'
#' @description
#'   This function defines the `Population` class with its attributes and methods.
#'
#' @param pop_size Maximum population size
#' @param num_features Number of features selected
#' @param num_obj Number of objective (default: 2) //Este parámetro no debería estar es esta versión del algoritmo//
#' @param clones Number of clones (default: 0)
#' @param p_mutation Mutation probability (default: 0.7)
#' @param p_mut_ind Individual mutation probability (default: 0.4)
#' @param p_mut_fea Feature mutation probability (default: 0.4)
#' @param p_mut_coord Mutation coordinates (default: 0.2)
#' @param mut_coord Mutation coordinates (default: 0)
#' @param objective The objective (default: DIST-EPS)
#'
#' @return Object of the `Population` class
#'
#' @export
#'
#' @examples
#' pop <- Population(pop_size = 10, num_features = 3, objective = "DIST-EPS")
Population <- function(pop_size, num_features, num_obj = 2, clones = 0, p_mutation = 0.7,
                       p_mut_ind = 0.4, p_mut_fea = 0.4, p_mut_coord = 0.2, mut_coord = 0,
                       objective = "DIST-EPS")
{
  population <- list(
    num_features = num_features,
    pop_size = pop_size,
    clones = clones,
    solution_list = list(),
    num_fronts = NULL,
    num_sol_fron = NULL,
    df_solutions = NULL,
    fronts = list(),
    p_mutation = p_mutation,
    p_mut_ind = p_mut_ind,
    p_mut_fea = p_mut_fea,
    p_mut_coord = p_mut_coord,
    mut_coord = mut_coord,
    num_obj = num_obj,
    objective = objective
  )

  class(population) <- "Population"
  return(population)
}



#' @title Generate initial population
#'
#' @description
#'   This function generates the initial population.
#'
#' @param population Population class object
#' @param data A data.frame with the input features
#' @param inputs A vector of input features
#' @param output A vector with class labels
#' @param costs A vector with feature costs //Esto no debería estar en este algoritmo//
#'
#' @return Population class object updated with the generated initial population
#'
#' @export
#'
#' @examples
#' set.seed(42)
#' data <- data.frame(
#'   x1 = rnorm(20), x2 = rnorm(20),
#'   y = sample(c(-1, 1), 20, replace = TRUE)
#' )
#' pop <- Population(pop_size = 2, num_features = 2, objective = "DIST-EPS")
#' pop <- generate_initial_population(pop, data, c("x1", "x2"), "y", c(1, 1))
generate_initial_population <- function(population, data, inputs, output, costs) {
  message(sprintf("Creating initial population of size %d", population$pop_size))
  num_sol <- 0

  while (length(population$solution_list) < population$pop_size) {
    # Create and generate a random solution
    sol <- Solution(num_sol, population$num_features, population$objective)
    sol <- generate_random_solution(sol, data, inputs, output)
    sol <- evaluate_solution(sol, data, output)

    # Check if the solution was successfully evaluated
    if (sol$successful_evaluation) {
      # Additional clone check if necessary
      if (population$clones == 1 || check_clones(population, sol) == 0) {
        # Add the evaluated and updated solution to the population
        population$solution_list[[num_sol + 1]] <- sol
        num_sol <- num_sol + 1
      } else {
        message("Detected and skipped clone.")
      }
    } else {
      message("Disregarded invalid solution.")
    }
  }
  return(population)
}



#' @title Print population
#'
#' @description
#'   This function prints the population.
#'
#' @param population Population class object
#'
#' @return Prints the population in the console
#'
#' @export
#'
#' @examples
#' set.seed(42)
#' data <- data.frame(
#'   x1 = rnorm(20), x2 = rnorm(20),
#'   y = sample(c(-1, 1), 20, replace = TRUE)
#' )
#' pop <- Population(pop_size = 2, num_features = 2, objective = "DIST-EPS")
#' pop <- generate_initial_population(pop, data, c("x1", "x2"), "y", c(1, 1))
#' pop <- fnds(pop)
#' print_population(pop)
print_population <- function(population) {
  # Convert the list of solutions to a data frame
  df_solutions <- do.call(rbind, lapply(population$solution_list, function(solution) {
    solution_dict <- to_dict(solution)
    # Ensure all values are vectors and not lists
    solution_dict <- lapply(solution_dict, function(item) if(is.list(item)) unlist(item) else item)
    # Convert to data frame
    solution_df <- as.data.frame(t(solution_dict), stringsAsFactors = FALSE)
    return(solution_df)
  }))

  # Convert necessary columns to numeric
  df_solutions$FRONT <- as.numeric(as.character(df_solutions$FRONT))

  # Verify if the FRONT column is now numeric
    if(!is.numeric(df_solutions$FRONT)) {
    stop("FRONT column could not be converted to numeric.")
  }

  # Sort the data frame by the 'FRONT' column
  df_solutions <- df_solutions[order(df_solutions$FRONT), ]
}



#' @title Check for Clones
#'
#' @description
#'   This function checks if a solution already exists in the solution list to avoid duplicates.
#'
#' @param population Population
#' @param solution_eva Solution object to be evaluated for duplication
#'
#' @return 1 if a clone is found, otherwise 0
#'
#' @export
#'
#' @examples
#' set.seed(42)
#' data <- data.frame(
#'   x1 = rnorm(20), x2 = rnorm(20),
#'   y = sample(c(-1, 1), 20, replace = TRUE)
#' )
#' pop <- Population(pop_size = 2, num_features = 2, objective = "DIST-EPS")
#' pop <- generate_initial_population(pop, data, c("x1", "x2"), "y", c(1, 1))
#' sol <- pop$solution_list[[1]]
#' clone_res <- check_clones(pop, sol)
check_clones <- function(population, solution_eva) {
  clone <- 0
  for (sol in population$solution_list) {
    # Compare all relevant features
    if (population$objective == 'DIST-EPS') {
      sol_obj_a = sol$objective[[1]]
      sol_obj_b = sol$objective[[2]]
      soleva_obj_a = solution_eva$objective[[1]]
      soleva_obj_b = solution_eva$objective[[2]]

    } else if (population$objective == 'MC') {
      sol_obj_a = sol$objective[[3]]
      sol_obj_b = sol$objective[[4]]
      soleva_obj_a = solution_eva$objective[[3]]
      soleva_obj_b = solution_eva$objective[[4]]
    }

    if (all(sol$features == solution_eva$features) &&
          abs(sol_obj_a - soleva_obj_a) < 0.1 && # 1e-5
          abs(sol_obj_b - soleva_obj_b) < 0.1 && # 1e-5
          all(abs(sol$plane_coord - solution_eva$plane_coord) < 0.001) && # 1e-5
          identical(sol$vectors, solution_eva$vectors)) {
        clone <- 1
      break
    } else {
      clone <- 0
    }
  }
  return(clone)
}



#' @title Fast Non-Dominated Sorting (fnds)
#'
#' @description
#'   This function performs the Fast Non-Dominated Sorting on a given population.
#'
#' @param population Population class object
#'
#' @return Updated Population class object with sorted fronts
#'
#' @export
#'
#' @examples
#' set.seed(42)
#' data <- data.frame(
#'   x1 = rnorm(20), x2 = rnorm(20),
#'   y = sample(c(-1, 1), 20, replace = TRUE)
#' )
#' pop <- Population(pop_size = 2, num_features = 2, objective = "DIST-EPS")
#' pop <- generate_initial_population(pop, data, c("x1", "x2"), "y", c(1, 1))
#' pop <- fnds(pop)
fnds <- function(population) {
  pop_size <- length(population$solution_list)
  fronts <- list()
  front <- list(solutions = numeric())

  # Initialize solutions
  for (i in 1:pop_size) {
    population$solution_list[[i]]$sol_dom_by <- 0
    population$solution_list[[i]]$dominates_list <- list()
    population$solution_list[[i]]$list_dominated_by <- list()
    population$solution_list[[i]]$front <- -1
  }

  # Calculate dominance
  for (n in 1:(pop_size - 1)) {
    for (s in (n + 1):pop_size) {
      if (dominate(population$solution_list[[n]], population$solution_list[[s]])) {
        population$solution_list[[n]]$dominates_list <- unique(c(population$solution_list[[n]]$dominates_list, s))
        population$solution_list[[s]]$list_dominated_by <- unique(c(population$solution_list[[s]]$list_dominated_by, n))
        population$solution_list[[s]]$sol_dom_by <- population$solution_list[[s]]$sol_dom_by + 1

      } else if (dominate(population$solution_list[[s]], population$solution_list[[n]])) {
        population$solution_list[[s]]$dominates_list <- unique(c(population$solution_list[[s]]$dominates_list, n))
        population$solution_list[[n]]$list_dominated_by <- unique(c(population$solution_list[[n]]$list_dominated_by, s))
        population$solution_list[[n]]$sol_dom_by <- population$solution_list[[n]]$sol_dom_by + 1
      }
    }
  }

  # Assign front 1 to non-dominated solutions
  for (i in 1:pop_size) {
    if (population$solution_list[[i]]$sol_dom_by == 0) {
      front$solutions <- c(front$solutions, i)
      population$solution_list[[i]]$front <- 1
    }
  }
  fronts[[1]] <- front

  # Create subsequent fronts
  r <- 1
  while (length(fronts[[r]]$solutions) > 0) {
    new_front <- numeric()
    for (sol_num in fronts[[r]]$solutions) {
      dominated <- population$solution_list[[sol_num]]$dominates_list
      for (i in dominated) {
        population$solution_list[[i]]$sol_dom_by <- population$solution_list[[i]]$sol_dom_by - 1
        if (population$solution_list[[i]]$sol_dom_by == 0) {
          new_front <- c(new_front, i)
          population$solution_list[[i]]$front <- r + 1
        }
      }
    }
    r <- r + 1
    front <- list(solutions = new_front)
    fronts[[r]] <- front
  }

  population$fronts <- fronts
  return(population)
}



#' @title Tournament Selection for Parent
#'
#' @description
#'   This function selects a parent solution from the population using tournament selection.
#'
#' @param population Population class object
#'
#' @return Selected parent solution
#'
#' @export
#'
#' @examples
#' set.seed(42)
#' data <- data.frame(
#'   x1 = rnorm(20), x2 = rnorm(20),
#'   y = sample(c(-1, 1), 20, replace = TRUE)
#' )
#' pop <- Population(pop_size = 4, num_features = 2, objective = "DIST-EPS")
#' pop <- generate_initial_population(pop, data, c("x1", "x2"), "y", c(1, 1))
#' pop <- fnds(pop)
#' parent <- tournament_select_parent(pop)
tournament_select_parent <- function(population) {
  # Randomly select two candidates
  num1 <- sample(seq_len(population$pop_size), 1)
  num2 <- sample(seq_len(population$pop_size), 1)

  # Ensure they are not the same candidate
  while (num1 == num2) {
    num2 <- sample(seq_len(population$pop_size), 1)
  }

  # Extract candidates directly
  parent1 <- population$solution_list[[num1]]
  parent2 <- population$solution_list[[num2]]

  # Compare the two parents and return the "better" one
  return(compare_solutions(parent1, parent2))
}



#' @title Crossover Solutions
#'
#' @description
#'   This function performs crossover between two parent solutions to create offspring.
#'
#' @param population Population class object
#' @param parent First parent solution
#' @param mother Second parent solution
#' @param num Unique identifier for the new solutions
#' @param data A data.frame with the input features
#' @param output A vector with class labels
#' @param costs A vector with feature costs //¿Lo quitamos?
#'
#' @importFrom stats runif
#'
#' @return List of two new offspring solutions
#'
#' @export
#'
#' @examples
#' set.seed(42)
#' data <- data.frame(
#'   x1 = rnorm(20), x2 = rnorm(20),
#'   y = sample(c(-1, 1), 20, replace = TRUE)
#' )
#' pop <- Population(pop_size = 4, num_features = 2, objective = "DIST-EPS")
#' pop <- generate_initial_population(pop, data, c("x1", "x2"), "y", c(1, 1))
#' pop <- fnds(pop)
#' parent1 <- tournament_select_parent(pop)
#' parent2 <- tournament_select_parent(pop)
#' offspring <- crossover_solutions(pop, parent1, parent2, 5, data, "y", c(1, 1))
crossover_solutions <- function(population, parent, mother, num, data, output, costs) {
  # Create four offspring solutions using crossover from the parents
  # First child: father's features and coordinates, father's first vector and mother's second vector
  child <- Solution(num, population$num_features, population$objective)
  child$features <- parent$features
  child$plane_coord <- parent$plane_coord
  child$vectors <- list(parent$vectors[[1]], mother$vectors[[2]])

  # Second child: mother's features and coordinates, father's first vector and mother's second vector
  child1 <- Solution(num, population$num_features, population$objective)
  child1$features <- mother$features
  child1$plane_coord <- mother$plane_coord
  child1$vectors <- list(parent$vectors[[1]], mother$vectors[[2]])

  # Third child: father's features and coordinates, mother's first vector and father's second vector
  child2 <- Solution(num, population$num_features, population$objective)
  child2$features <- parent$features
  child2$plane_coord <- parent$plane_coord
  child2$vectors <- list(mother$vectors[[1]], parent$vectors[[2]])

  # Fourth child: mother's features and coordinates, mother's first vector and father's second vector
  child3 <- Solution(num, population$num_features, population$objective)
  child3$features <- mother$features
  child3$plane_coord <- mother$plane_coord
  child3$vectors <- list(mother$vectors[[1]], parent$vectors[[2]])

  # Evaluate the offspring solutions 
  child <- evaluate_solution(child, data, output)
  child1 <- evaluate_solution(child1, data, output)
  child2 <- evaluate_solution(child2, data, output)
  child3 <- evaluate_solution(child3, data, output)

  # First round of comparison
  res1 <- select_best_offspring(child, child1)
  winner1 <- res1$winner
  loser1  <- res1$loser

  res2 <- select_best_offspring(child2, child3)
  winner2 <- res2$winner
  loser2  <- res2$loser

  # Final between both winners
  res_final <- select_best_offspring(winner1, winner2)
  best <- res_final$winner
  second_candidate <- res_final$loser

  # Choose second best among remaining three (two losers + second best from final)
  second_pool <- list(loser1, loser2, second_candidate)
  second_best <- second_pool[[1]]
  for (i in 2:length(second_pool)) {
    res <- select_best_offspring(second_best, second_pool[[i]])
    second_best <- res$winner
  }

  # Assign solution IDs
  best$num <- num
  second_best$num <- num + 1

  return(list(best, second_best))
}


#' @title Compares two offspring and selects the best one based on dominance relations.
#'
#' @description
#'   This method uses the concept of dominance in multi-objective optimization to compare the two offspring.
#'   If one offspring dominates the other, the dominant one is selected. If neither dominates the other,
#'   the tie is broken randomly.
#'
#' @param child1 The first offspring to compare.
#' @param child2 The second offspring to compare.
#'
#' @return List[Solution][Solution]
#'   A list with two solutions, where the first one is the better offspring (the winner), and the second one
#'   is the inferior offspring (the loser).
#'
#' @export
#'
#' @examples
#' sol1 <- Solution(num = 1, num_features = 2, obj = "DIST-EPS")
#' sol1$objective <- c(1, 0, 0, 0)
#' sol2 <- Solution(num = 2, num_features = 2, obj = "DIST-EPS")
#' sol2$objective <- c(0, 1, 1, 1)
#' best_offspring <- select_best_offspring(sol1, sol2)
select_best_offspring <- function(child1, child2) {
  dominates <- dominate2(child1, child2)
  if (dominates == 1) {        # child1 dominates child2
    return(list(winner = child1, loser = child2))
  } else if (dominates == 2) { # child2 dominates child1
    return(list(winner = child2, loser = child1))
  } else {                     # Tie: decide randomly
    if (runif(1) < 0.5) {
      return(list(winner = child1, loser = child2))
    } else {
      return(list(winner = child2, loser = child1))
    }
  }
}


#' @title Create New Population
#'
#' @description
#'   This function creates a new population through crossover.
#'
#' @param population Population class object
#' @param data A data.frame with the input features
#' @param inputs A vector of input features
#' @param output A vector with class labels
#' @param costs A vector with feature costs //¿Lo quitamos?
#'
#' @importFrom stats runif
#'
#' @return Updated Population class object with the new population generated
#'
#' @export
#'
#' @examples
#' set.seed(42)
#' data <- data.frame(
#'   x1 = rnorm(20), x2 = rnorm(20),
#'   y = sample(c(-1, 1), 20, replace = TRUE)
#' )
#' pop <- Population(pop_size = 2, num_features = 2, objective = "DIST-EPS")
#' pop <- generate_initial_population(pop, data, c("x1", "x2"), "y", c(1, 1))
#' pop <- fnds(pop)
#' pop_new <- new_population(pop, data, c("x1", "x2"), "y", c(1, 1))
new_population <- function(population, data, inputs, output, costs) {
  i <- population$pop_size

  while (length(population$solution_list) < 2 * population$pop_size) {
    parent <- tournament_select_parent(population)
    mother <- tournament_select_parent(population)
    while (parent$num == mother$num) {
      mother <- tournament_select_parent(population)
    }

    offspring <- crossover_solutions(population, parent, mother, i, data, output, costs)
    for (child in offspring) {
      if (length(population$solution_list) >= 2 * population$pop_size) {
        break  # Avoid adding more children if the maximum size is reached
      }

      if (runif(1) < population$p_mutation) {
        child <- mutate_solution(child, population, data, inputs, output)
      }

      child <- evaluate_solution(child, data, output)
      while (!child$successful_evaluation) {
        child <- generate_random_solution(child, data, inputs, output)
        child <- evaluate_solution(child, data, output)
      }

      while (check_clones(population, child) == 1 || !child$successful_evaluation) {
        child <- mutate_solution(child, population, data,inputs, output)
        child <- evaluate_solution(child, data, output)
      }

      population$solution_list[[length(population$solution_list) + 1]] <- child
      i <- i + 1
    }
  }

  return(population)
}



#' @title Mutate Solution
#'
#' @description
#'   This function mutates a solution based on mutation probabilities and coordinates.
#'
#' @param solution Solution object to mutate
#' @param population Population class object that contains mutation parameters
#' @param data A data.frame with the input features
#' @param inputs A vector of input features
#' @param output A vector with class labels
#'
#' @importFrom stats runif
#'
#' @return Mutated solution object
#'
#' @export
#'
#' @examples
#' data <- data.frame(
#'   x1 = c(1, 2, 1, 5, 6, 7),
#'   x2 = c(1, 1, 2, 6, 5, 6),
#'   y = c(-1, -1, -1, 1, 1, 1)
#' )
#' pop <- Population(pop_size = 2, num_features = 2, objective = "DIST-EPS")
#' sol <- Solution(num = 1, num_features = 2, obj = "DIST-EPS")
#' sol <- generate_random_solution(sol, data, c("x1", "x2"), "y")
#' mut_sol <- mutate_solution(sol, pop, data, c("x1", "x2"), "y")
mutate_solution <- function(solution, population, data, inputs, output) {
  random <- runif(1)

  # Vector mutation
  if (random < population$p_mut_ind) {
    solution <- mutate_vectors(solution, data, output)
  }

  # Feature mutation
  if (length(inputs) > population$num_features) {
    for (i in seq_along(solution$features)) {
      rand_val <- runif(1)
      if (rand_val < population$p_mut_fea) {
        possible_features <- setdiff(inputs, solution$features)
        solution$features[i] <- sample(possible_features, 1)
        solution$plane_coord[i] <- runif(1, -1, 1)
      }
    }
  }

  # Plane mutation
  for (i in seq_along(solution$plane_coord)) {
    rand_val <- runif(1)
    if (rand_val < population$p_mut_coord) {
      if (population$mut_coord == 0) {
        percentage <- runif(1, 0, 0.25)
        rand_direction <- runif(1)
        if (rand_direction < 0.5) {  # increment
          solution$plane_coord[i] <- solution$plane_coord[i] * (1 + percentage)
          if (solution$plane_coord[i] > 1) {
            solution$plane_coord[i] <- solution$plane_coord[i] / 10
          }
        } else {  # decrement
          solution$plane_coord[i] <- solution$plane_coord[i] * (1 - percentage)
          if (solution$plane_coord[i] < -1) {
            solution$plane_coord[i] <- solution$plane_coord[i] / 10
          }
        }
      } else {
        solution$plane_coord[i] <- runif(1, -1, 1)
      }
    }
  }

  return(solution)
}



#' @title Reduce Population Size
#'
#' @description
#'   This function reduces the population size based on fronts and crowding distance.
#'
#' @param population Population class object to reduce
#'
#' @return Reduced Population class object
#'
#' @export
#'
#' @examples
#' set.seed(42)
#' data <- data.frame(
#'   x1 = rnorm(20), x2 = rnorm(20),
#'   y = sample(c(-1, 1), 20, replace = TRUE)
#' )
#' pop <- Population(pop_size = 5, num_features = 2, objective = "DIST-EPS")
#' pop <- generate_initial_population(pop, data, c("x1", "x2"), "y", c(1, 1))
#' pop <- fnds(pop)
#' # Create a larger population to reduce
#' pop <- new_population(pop, data, c("x1", "x2"), "y", c(1, 1))
#' pop <- fnds(pop)
#' pop_reduced <- reduce_population(pop)
reduce_population <- function(population) {
  # Create new empty population
  cat("Creating new population with", population$num_fea, "features\n")
  new_population <- Population(pop_size = population$pop_size,
                               num_features = population$num_fea,
                               objective = population$objective)

  # Front counter
  front_counter <- 1
  # Add solutions from the fronts to the new population until sizing is complete
  while (length(new_population$solution_list) +
         length(population$fronts[[front_counter]]$solutions) <= population$pop_size)
  {
    front_solutions <- lapply(population$solution_list, function(sol) {
      if (sol$front == front_counter) return(sol) else return(NULL)
    })
    front_solutions <- front_solutions[!sapply(front_solutions, is.null)]
    new_population$solution_list <- c(new_population$solution_list, front_solutions)

    front_counter <- front_counter + 1
  }

  # If the population is not complete, select the remaining best solutions using crowding distance
  if (length(new_population$solution_list) < population$pop_size
      && front_counter <= length(population$fronts))
  {
    front_solutions <- crowding_distance(population, front_counter)
    remaining <- population$pop_size - length(new_population$solution_list)
    new_population$solution_list <- c(new_population$solution_list, front_solutions[1:remaining])
  }

  # Rename the solutions
  for (i in seq_along(new_population$solution_list)) {
    new_population$solution_list[[i]]$num <- i - 1
  }
  return(new_population)
}



#' @title Crowding Distance Calculation
#'
#' @description
#'   This function calculates the crowding distance for solutions in a specified front.
#'
#' @param population Population class object
#' @param num_front The front number to calculate the crowding distance for
#'
#' @return List of solutions with updated crowding distances
#'
#' @export
#'
#' @examples
#' set.seed(42)
#' data <- data.frame(
#'   x1 = rnorm(20), x2 = rnorm(20),
#'   y = sample(c(-1, 1), 20, replace = TRUE)
#' )
#' pop <- Population(pop_size = 4, num_features = 2, objective = "DIST-EPS")
#' pop <- generate_initial_population(pop, data, c("x1", "x2"), "y", c(1, 1))
#' pop <- fnds(pop)
#' crowd_sols <- crowding_distance(pop, num_front = 1)
crowding_distance <- function(population, num_front) {
  # Extract solutions from the specified front
  front_solutions <- lapply(population$solution_list, function(sol) {
    if (sol$front == num_front) return(sol) else return(NULL)
  })
  front_solutions <- front_solutions[!sapply(front_solutions, is.null)]

  num_sol_front <- length(front_solutions)

  # Initialize crowding distance to zero for all solutions
  for (k in seq_along(front_solutions)) {
    front_solutions[[k]]$crowding_distance <- 0
  }

  # Process each objective if there are enough solutions to compare
  if (num_sol_front > 2) {
    for (obj in 1:population$num_obj) {
      # Sort solutions based on the current objective
      front_solutions <- front_solutions[order(sapply(front_solutions, function(x) x$objective[[obj]]))]

      # Assign infinite crowding distance to extreme solutions
      front_solutions[[1]]$crowding_distance <- Inf
      front_solutions[[num_sol_front]]$crowding_distance <- Inf

      max_obj <- front_solutions[[num_sol_front]]$objective[[obj]]
      min_obj <- front_solutions[[1]]$objective[[obj]]

      if (max_obj > min_obj) {
        for (j in 2:(num_sol_front - 1)) {
          increment <- (front_solutions[[j + 1]]$objective[[obj]] - front_solutions[[j - 1]]$objective[[obj]]) / (max_obj - min_obj)
          front_solutions[[j]]$crowding_distance <- front_solutions[[j]]$crowding_distance + increment
          }
      }
    }
  }

  # Sort solutions by crowding distance from highest to lowest
  front_solutions <- front_solutions[order(sapply(front_solutions, function(x) x$crowding_distance), decreasing = TRUE)]

  return(front_solutions)
}


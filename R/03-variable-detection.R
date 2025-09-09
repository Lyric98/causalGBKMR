#' @file 03-variable-detection.R
#' @title Variable pattern detection for g-BKMR
#' @description Automatically detects variable structures in prepared g-BKMR data.
#' Identifies exposure variables, time-dependent covariates, and their naming patterns.

#' Detect variable patterns in g-BKMR data
#'
#' @description Automatically detects the structure of exposure and
#' time-dependent covariate variables in prepared g-BKMR data. This function
#' is essential for the automated analysis pipeline.
#'
#' @param data Data frame containing the variables to analyze. Should be in
#'   g-BKMR format (wide format with proper variable naming).
#' @param T Integer. Number of time points in the study.
#'
#' @return A list containing the detected variable structure:
#' \describe{
#'   \item{p}{Integer. Number of exposure variables per time point}
#'   \item{Ldim}{Integer. Number of time-dependent covariates per time point}
#'   \item{td_covariate_names}{Character vector. Names of time-dependent covariates}
#'   \item{detected_pattern}{Character. Pattern used for detection}
#'   \item{baseline_td_vars}{Character vector. Baseline time-dependent variables}
#'   \item{td_vars_by_time}{List. Time-dependent variables for each time point}
#' }
#'
#' @details
#' The function recognizes several naming patterns for time-dependent covariates:
#' \itemize{
#'   \item **known_with_underscore**: bmi_0, bp_0, bmi_1, bp_1, etc.
#'   \item **known_ending_zero**: bmi0, bp0, bmi1, bp1, etc.
#'   \item **generated_format**: waist0_1, waist0_2, waist1_1, waist1_2, etc.
#'   \item **generic_format**: td_covariate1_0, td_covariate2_0, etc.
#' }
#'
#' For exposure variables, it looks for the pattern: logM1_0, logM2_0, logM1_1, logM2_1, etc.
#'
#' @examples
#' \dontrun{
#' # Create test data in g-BKMR format
#' test_data <- data.frame(
#'   id = 1:100,
#'   sex = rbinom(100, 1, 0.5),
#'   bmi_0 = rnorm(100, 25, 3),
#'   bp_0 = rnorm(100, 120, 15),
#'   logM1_0 = rnorm(100, 0, 1),
#'   logM2_0 = rnorm(100, 0, 1),
#'   logM1_1 = rnorm(100, 0, 1),
#'   logM2_1 = rnorm(100, 0, 1),
#'   bmi_1 = rnorm(100, 25, 3),
#'   bp_1 = rnorm(100, 120, 15),
#'   Y = rnorm(100, 0, 1)
#' )
#'
#' # Detect variable patterns
#' detection_result <- detect_variable_patterns(test_data, T = 2)
#'
#' # View results
#' print(detection_result)
#' }
#'
#' @export
detect_variable_patterns <- function(data, T) {

  # Detect number of exposures per time point (p)
  p <- length(grep("^logM\\d+_0$", names(data)))

  # Try to detect actual time-dependent covariate names from data
  # Look for baseline TD covariates with various patterns
  baseline_td_patterns <- list(
    # Pattern 1: Known covariate names with 0 (e.g., bmi_0, bp_0)
    "known_with_underscore" = "^[a-zA-Z_]+_0$",
    # Pattern 2: Known covariate names ending with 0 (e.g., bmi0, bp0)
    "known_ending_zero" = "^[a-zA-Z_]+0$",
    # Pattern 3: Generated format (waist0_1, waist0_2, etc.)
    "generated_format" = "^[a-zA-Z_]+0_\\d+$",
    # Pattern 4: Generic td_covariate format
    "generic_format" = "^td_covariate\\d+_0$"
  )

  baseline_td_vars <- character(0)
  detected_pattern <- ""
  td_covariate_names <- character(0)

  # Try each pattern to find baseline TD covariates
  for (pattern_name in names(baseline_td_patterns)) {
    pattern <- baseline_td_patterns[[pattern_name]]
    matches <- grep(pattern, names(data), value = TRUE)
    if (length(matches) > 0) {
      baseline_td_vars <- matches
      detected_pattern <- pattern_name

      # Extract base names from the matches
      if (pattern_name == "known_with_underscore") {
        # Extract names like "bmi", "bp" from "bmi_0", "bp_0"
        td_covariate_names <- gsub("_0$", "", matches)
      } else if (pattern_name == "known_ending_zero") {
        # Extract names like "bmi", "bp" from "bmi0", "bp0"
        td_covariate_names <- gsub("0$", "", matches)
      } else if (pattern_name == "generated_format") {
        # For waist0_1, waist0_2, extract "waist" and number the covariates
        base_name <- gsub("0_\\d+$", "", matches[1])
        n_covariates <- length(matches)
        td_covariate_names <- paste0(base_name, 1:n_covariates)
      } else if (pattern_name == "generic_format") {
        # For td_covariate1_0, td_covariate2_0, extract the names
        td_covariate_names <- gsub("_0$", "", matches)
      }
      break
    }
  }

  # If no baseline TD covariates found, set defaults
  if (length(baseline_td_vars) == 0) {
    baseline_td_vars <- character(0)
    Ldim <- 0
    td_covariate_names <- character(0)
    warning("No baseline time-dependent covariates detected. Setting Ldim = 0.")
  } else {
    Ldim <- length(baseline_td_vars)
  }

  # Detect time-dependent covariates for each time point t = 1, 2, ..., T-1
  td_vars_by_time <- list()

  if (Ldim > 0 && T > 1) {
    for (t in 1:(T-1)) {
      td_vars_t <- character(0)

      # Generate possible variable names based on detected pattern
      if (detected_pattern == "known_with_underscore") {
        # Look for bmi_1, bp_1, bmi_2, bp_2, etc.
        potential_vars <- paste0(td_covariate_names, "_", t)
      } else if (detected_pattern == "known_ending_zero") {
        # Look for bmi1, bp1, bmi2, bp2, etc.
        potential_vars <- paste0(td_covariate_names, t)
      } else if (detected_pattern == "generated_format") {
        # Look for waist1_1, waist1_2, waist2_1, waist2_2, etc.
        base_name <- gsub("0_\\d+$", "", baseline_td_vars[1])
        potential_vars <- paste0(base_name, t, "_", 1:Ldim)
      } else if (detected_pattern == "generic_format") {
        # Look for td_covariate1_1, td_covariate2_1, etc.
        potential_vars <- paste0(td_covariate_names, "_", t)
      }

      # Check which of these potential variables exist in the data
      existing_vars <- intersect(potential_vars, names(data))

      if (length(existing_vars) == Ldim) {
        td_vars_t <- existing_vars
      } else {
        # Fallback: try pattern matching
        time_patterns <- c(
          paste0("^.*", t, "_\\d+$"),    # any_name1_1, any_name1_2, etc.
          paste0("^.*_", t, "_\\d+$"),   # any_name_1_1, any_name_1_2, etc.
          paste0("^.*", t, "$")          # any_name1 (single TD covariate)
        )

        for (pattern in time_patterns) {
          matches <- grep(pattern, names(data), value = TRUE)
          # Exclude exposure variables and outcome/id
          matches <- matches[!grepl("^logM\\d+_", matches)]
          matches <- matches[!matches %in% c("Y", "id")]

          if (length(matches) >= Ldim) {
            td_vars_t <- matches[1:Ldim]
            break
          }
        }
      }

      if (length(td_vars_t) != Ldim) {
        stop("Cannot detect time-dependent covariates for time ", t,
             ". Expected ", Ldim, " variables but found ", length(td_vars_t),
             ". Available variables: ", paste(names(data), collapse = ", "))
      }

      td_vars_by_time[[t]] <- td_vars_t
    }
  }

  result <- list(
    p = p,
    Ldim = Ldim,
    td_covariate_names = td_covariate_names,
    detected_pattern = detected_pattern,
    baseline_td_vars = baseline_td_vars,
    td_vars_by_time = td_vars_by_time
  )

  # Add informative class for printing
  class(result) <- c("gbkmr_detection", "list")

  return(result)
}

#' Print method for gbkmr_detection objects
#'
#' @param x A gbkmr_detection object
#' @param ... Additional arguments (not used)
#'
#' @return Invisible x
#' @export
print.gbkmr_detection <- function(x, ...) {
  cat("g-BKMR Variable Detection Results\n")
  cat("================================\n")
  cat("Number of exposures per time point (p):", x$p, "\n")
  cat("Number of time-dependent covariates per time point (Ldim):", x$Ldim, "\n")

  if (length(x$td_covariate_names) > 0) {
    cat("Time-dependent covariate names:", paste(x$td_covariate_names, collapse = ", "), "\n")
  }

  cat("Detected naming pattern:", x$detected_pattern, "\n")

  if (length(x$baseline_td_vars) > 0) {
    cat("Baseline time-dependent variables:", paste(x$baseline_td_vars, collapse = ", "), "\n")
  }

  if (length(x$td_vars_by_time) > 0) {
    cat("\nTime-dependent variables by time point:\n")
    for (t in names(x$td_vars_by_time)) {
      cat("  Time", t, ":", paste(x$td_vars_by_time[[t]], collapse = ", "), "\n")
    }
  }

  invisible(x)
}

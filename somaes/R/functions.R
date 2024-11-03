# functions.R

set.seed(137)


## Functions to be implemented

### Discretization algorithms

#### *Algoritmos de discretización para un solo atributo y para un dataset completo (ambas opciones): Igual frecuencia e igual anchura*


# Helper function to create names for factors based on cut_points
create_factor_names <- function(cut_points) {
  c <- cut_points
  names <- sprintf("I1:(-Inf,%.3f)", c[1])  # First bin is special
  for (i in 2:length(c)) {
    names <- c(names, sprintf("I%d:(%.3f,%.3f)", i, c[i-1], c[i]))
  }
  # Our last bin is also special, tack it onto the end
  names <- c(names, sprintf("I%d:(%.3f,+Inf)", length(c) + 1, c[length(c)]))
  return(names)
}


# Helper function to categorize individual values based on cut points
categorize <- function(value, cut_points, factor_names) {

  # There must be a more elegant way, but I'll just use if() and for()
  if (value < cut_points[1]) {  # first bin
    return(factor_names[1])
  } else if (value >= cut_points[length(cut_points)]) {  # last bin
    return(factor_names[length(factor_names)])
  }
  for (i in 2:length(factor_names)) {  # otherwise we are in the middle
    if (cut_points[i - 1] <= value && value < cut_points[i]) {
      return(factor_names[i])
    }
  }
  return(NA)  # if no category is found
}


# ----------------------------------------------------------------------------
#   Equal-width discretization algorithm, vector version
# ----------------------------------------------------------------------------

discretizeEW <- function(x, num_bins) {
  ### Discretize numerical vector into categorical with bins of
  ###   equal width

  #    min                      max
  #     |    |    |    |    |    |
  #  I1------| I2 | I3 | ...|------I{num_bins}
  #          |              |
  #          ^cut_points[0] |
  #                         ^cut_points[num_bins-1]

  bin_width <- (max(x) - min(x)) / num_bins
  cpoints   <- min(x) + (1:(num_bins-1)) * bin_width  # R is 1-based

  # Make factor names
  fnames <- create_factor_names(cpoints)

  # Apply categorize() to each value in x
  x_discretized <- sapply(x, function(val) categorize(val, cpoints, fnames))

  return(list(discretized = x_discretized, cut_points = cpoints))
}


# ----------------------------------------------------------------------------
#   Equal-frequency discretization algorithm, vector version
# ----------------------------------------------------------------------------

discretizeEF <- function(x, num_bins) {
  ### Discretize numerical vector into categorical with bins of
  ###   equal frequency

  x_ordered <- sort(x)  # sort input into ascending order for later use

  # Create num_bins bins with an equal number of members with floor division
  bin_populations <- rep(floor(length(x) / num_bins), num_bins)

  # Unless the size of x happens to be an exact multiple of the number
  #   of bins, there will be leftovers from the floor division, we retrieve
  #   the number of leftovers with the modulo and assign then to random bins
  leftovers <- sample(1:num_bins, length(x) %% num_bins, replace = FALSE)
  bin_populations[leftovers] <- bin_populations[leftovers] + 1

  # Split the vector into bins of size bin_populations
  last_this_bin <- cumsum(bin_populations)
  first_next_bin <- c(1, last_this_bin[-length(last_this_bin)] + 1)
  end_of_bin <- x_ordered[last_this_bin[-length(last_this_bin)]]
  beginning_of_bin <- x_ordered[first_next_bin[-1]]

  # Take cut points at the arithmetic mean of the endpoints of the bins
  cpoints <- end_of_bin + (beginning_of_bin - end_of_bin) / 2

  # Make factor names
  fnames <- create_factor_names(cpoints)

  # Apply categorize() to each value in x
  x_discretized <- sapply(x, function(val) categorize(val, cpoints, fnames))

  return(list(discretized = x_discretized, cut_points = cpoints))
}


# ----------------------------------------------------------------------------
#   Equal-width discretization algorithm, multi-column dataframe version
# ----------------------------------------------------------------------------

discretize_EW_by_column <- function(df, num_bins, keep_original=FALSE) {
  ### Discretize input dataframe column by column with into num_bins bins of
  ###  equal width and return in a new dataframe. Set keep_original=True to
  ###  keep the input columns

  ret_df <- data.frame(matrix(nrow = nrow(df), ncol = 0))

  # Apply equal-width discretization to each column
  ew <- sapply(df, function(col) {
    discretizeEW(as.numeric(col), num_bins)$discretized
  })

  # Combine (optionally) the original df and discretization into a new df
  for (col in names(df)) {
    if (keep_original) {
      ret_df[[col]] <- df[[col]]  # original column
    }
    ret_df[[paste0(col, "_EW")]] <- ew[, col]  # add discretized column
  }

  return(ret_df)
}


# ----------------------------------------------------------------------------
#   Equal-frequency discretization algorithm, multi-column dataframe version
# ----------------------------------------------------------------------------

discretize_EF_by_column <- function(df, num_bins, keep_original=FALSE) {
  ### Discretize input dataframe column by column with into num_bins bins
  ###  with an approximately equal number of members and return in a new
  ###  dataframe. Set keep_original=True to keep the input columns

  ret_df <- data.frame(matrix(nrow = nrow(df), ncol = 0))

  # Apply equal-frequency discretization to each column
  ef <- lapply(df, function(col) {
    discretizeEF(as.numeric(col), num_bins)$discretized
  })

  # Combine (optionally) the original df and discretization into a new df
  for (col in names(df)) {
    if (keep_original) {
      ret_df[[col]] <- df[[col]]  # original column
    }
    ret_df[[paste0(col, "_EF")]] <- ef[[col]]  # add discretized column
  }

  return(ret_df)
}


# ----------------------------------------------------------------------------
#   Compare results of equal-width and equal-frequecy with original columns
# ----------------------------------------------------------------------------

discretize_EW_EF_by_column <- function(df, num_bins) {
  ### Discretize column by column using discretizeEF() and discretizeEW()
  ###  and combine with input dataframe into a new dataframe.

  ret_df <- data.frame(matrix(nrow = nrow(df), ncol = 0))

  # Apply equal frequency discretization column by column
  ef <- lapply(df, function(col) {
    discretizeEF(as.numeric(col), num_bins)$discretized
  })

  # Apply equal width discretization column by column
  ew <- lapply(df, function(col) {
    discretizeEW(as.numeric(col), num_bins)$discretized
  })

  for (col in names(df)) {
    ret_df[[col]] <- df[[col]]                 # original column
    ret_df[[paste0(col, "_EF")]] <- ef[[col]]  # EF discretized column
    ret_df[[paste0(col, "_EW")]] <- ew[[col]]  # EW discretized column
  }

  return(ret_df)
}


### Calculation of metrics for the attributes of a dataset

#### *Cálculo de métricas para los atributos de un dataset: Varianza y AUC para las variables contínuas y entropía para las discretas. La función deberá reconocer el tipo de atributo y actuar en consecuencia. Notese que en el caso del AUC, el dataset debe ser supervisado, es decir, es necesario especificar una variable clase binaria con la que evaluar el AUC de los atributos numéricos.*

# ----------------------------------------------------------------------------
#   Calculation of sample variance by column from a dataframe
# ----------------------------------------------------------------------------

calculate_variance <- function(df, sample=TRUE) {
  ### Calculate variance by column from a dataframe with numeric values.
  ###  If the option sample=True (the default) is passed,
  ###  the sample variance is returned, otherwise the population variance.
  ###  Returns an array with the selected variance for each column

  # Check that the input is a dataframe
  if (!is.data.frame(df)) {
    stop("error, input must be a DataFrame")
  }

  # Check that all columns are numeric
  if (!all(sapply(df, is.numeric))) {
    stop("error, all columns must be numeric")
  }

  # Convert dataframe to matrix for easier calculations
  x_cols <- as.matrix(df)
  N <- nrow(x_cols)
  ddof <- if (sample) 1 else 0  # Degrees of freedom

  # Population variance is defined as:
  #
  #                1
  #   variance =  ——— · sum_i=1_N( (x_i - xbar)^2 )
  #                N
  #
  # Sample variance is defined as:
  #
  #                1
  #   variance =  ——— · sum_i=1_N( (x_i - xbar)^2 )
  #               N-1
  #

  # Get the mean value of columns
  x_bar_cols <- colMeans(x_cols)

  # Get (x_cols - x_bar_cols) column-wise while avoiding R broadcasting
  difference_cols <- sweep(x_cols, 2, x_bar_cols, "-")

  # Get variance for each column
  variance <- (1 / (N - ddof)) * colSums(difference_cols ^ 2)

  return(variance)
}


# ----------------------------------------------------------------------------
#   ROC (AUC) calculation from two-column dataframe
# ----------------------------------------------------------------------------

calculate_roc_auc <- function(df, sample=TRUE) {
  ### Calculate receiver operating characteristic (ROC) curve and
  ###  accompanying area under curve (AUC) value for a variable (instead of
  ###  the more typical use with full classifier), given an input dataframe
  ###  with numeric (decimal) values in the first column, and booleans in
  ###  the second.

  # Ensure the input is a dataframe with two columns
  if (!is.data.frame(df) || ncol(df) != 2) {
    stop("Error: input must be a data frame with two columns")
  }

  # Ensure the first column is numeric
  if (!is.numeric(df[[1]])) {
    stop("Error: the first column must be numeric")
  }

  # Ensure the second column is boolean
  if (!is.logical(df[[2]])) {
    stop("Error: the second column must be boolean")
  }

  # Extract values
  unordered_numeric_values <- df[[1]]
  unordered_boolean_ground_truth <- df[[2]]

  # Order by the numeric values in descending order
  ordered_indices <- order(-unordered_numeric_values)
  candidate_thresholds <- unordered_numeric_values[ordered_indices]
  ground_truth <- unordered_boolean_ground_truth[ordered_indices]

  # Get total number of positives (P) and negatives (N)
  P <- sum(ground_truth)
  N <- length(ground_truth) - P

  # Initialize vectors for the ROC curve with one extra point for (1,1)
  num_curve_points <- length(ground_truth) + 1
  fpr_curve <- numeric(num_curve_points)
  tpr_curve <- numeric(num_curve_points)
  fpr_curve[num_curve_points] <- 1
  tpr_curve[num_curve_points] <- 1

  # Calculate TPR and FPR for each threshold
  for (i in seq_along(candidate_thresholds)) {
    c <- candidate_thresholds[i]
    predicted <- candidate_thresholds > c

    TP <- sum(ground_truth & predicted)
    FP <- sum(!ground_truth & predicted)

    tpr_curve[i] <- TP / P
    fpr_curve[i] <- FP / N
  }

  # Get AUC by numerical integration using trapezoidal rule
  AUC <- sum(diff(fpr_curve)*(head(tpr_curve, -1) + tail(tpr_curve, -1))/2)

  return(list(fpr_curve = fpr_curve, tpr_curve = tpr_curve, auc = AUC))
}


# ----------------------------------------------------------------------------
#   Shannon entropy calculation, dataframe version
# ----------------------------------------------------------------------------

calculate_entropy <- function(df) {
  ### Calulate Shannon entropy of a single dataframe column, note that all
  ###  data types are interpreted as categorical. Returns decimal.

  # Ensure the input is a dataframe with one column
  if (!is.data.frame(df) || ncol(df) != 1) {
    stop("Error: input must be a data frame with one column")
  }

  # Convert the column to a vector
  vector <- df[[1]]

  # Calculate the frequency of unique values
  unique_values <- unique(vector)
  if (length(unique_values) == 0) { return(NULL) } # entropy is undefined

  # Calculate probabilities for each unique value
  P <- table(vector) / length(vector)

  # Calculate entropy according to H(X) = -sum(p_i * log2(p_i))
  H <- -sum(P * log2(P[P > 0]))

  return(H)
}


# ----------------------------------------------------------------------------
#   Extract metrics of the appropriate type from any column in a dataframe
# ----------------------------------------------------------------------------

extract_dataset_metrics <- function(df) {
  ### Extract the implemented metrics (entropy, variance, auc) from
  ###  compatible columns of the input dataframe. Returns list of three
  ###  lists for each metric (entropy, variance, auc), each containing
  ###  yet another list with:
  ###        column name
  ###        metric as calculated by the helper functions
  ###
  ###  Note that for a column to have the Area Under Curve (auc) metric
  ###  calculated, it must be a numeric (float) column with a boolean column
  ###  immediately to its right. Booleans that have been used for AUC are
  ###  ignored for further metrics. Booleans that have NOT been used for
  ###  AUC are considered categorical data.

  # lists for the return list start out empty
  entropy_list <- list()
  variance_list <- list()
  auc_list <- list()

  # cycle through the columns of input the dataframe
  for (i in seq_along(df)) {

    c <- colnames(df)[i]
    col_type <- class(df[[c]])

    if (col_type %in% c("character", "factor", "integer")) {
      # We consider these categorical, get entropy
      e <- list(column=c,
                metric=as.numeric(calculate_entropy(df[c])))
      entropy_list[[length(entropy_list)+1]] <- e

      # logical columns can be used for AUC or categorical
    } else if (col_type == "logical") {
      if ((i-1) >= 1 && class(df[[i-1]]) != "numeric") {
        # If column to the left is NOT numeric, also get entropy
        e <- list(column=c,
                  metric=as.numeric(calculate_entropy(df[c])))
        entropy_list[[length(entropy_list)+1]] <- e
      }

    } else if (col_type == "numeric") {
      # Numeric, get variance
      v <- list(column=c,
                metric=as.numeric(calculate_variance(df[c])))
      variance_list[[length(variance_list)+1]] <- v

      # Check the column to the right, if it is logical...
      if ((i+1) <= ncol(df) && class(df[[i+1]]) == "logical") {
        # we can calculate the AUC
        a <- list(column=c,
                  metric=calculate_roc_auc(df[c((i), (i+1))])[[3]])
        auc_list[[length(auc_list)+1]] <- a
      }

    } else {
      # Print an informational message, but do NOT dump out
      message("Unknown column type at column: ", c)
    }
  }

  return(list(entropy=entropy_list, variance=variance_list, auc=auc_list))
}


# ----------------------------------------------------------------------------
#   Generate summary dataframe with the appropriate metrics per column
# ----------------------------------------------------------------------------

dataset_metrics_summary <- function(df, display_precision=3, append_input=FALSE) {
  ### Generate a summary dataframe with the implemented metrics (entropy,
  ###  variance, auc) from compatible columns of the input dataframe. Returns
  ###  a dataframe populated with text, with the number of significant digits
  ###  set by display_precision. Include append_input=True to include the
  ###  input dataframe in the output

  # Initialize a new dataframe with the same columns as the input dataframe
  summary_df <- data.frame(matrix(NA, nrow = 4, ncol = ncol(df)))
  colnames(summary_df) <- colnames(df)
  rownames(summary_df) <- c('(type)', 'variance', 'AUC', 'entropy')

  # Extract metrics from input dataframe using extract_dataset_metrics()
  metrics <- extract_dataset_metrics(df)
  e_list <- metrics$entropy
  v_list <- metrics$variance
  a_list <- metrics$auc

  # Populate the type row column by column
  for (c in colnames(df)) { summary_df['(type)', c] <- class(df[[c]]) }

  # Fill in the variances
  for (variance in v_list) {
    column_name <- variance$column
    variance_value <- round(variance$metric, display_precision)
    summary_df['variance', column_name] <- variance_value
  }

  # Fill in any AUC values
  for (auc in a_list) {
    column_name <- auc$column
    auc_value <- round(auc$metric, display_precision)
    summary_df['AUC', column_name] <- auc_value
  }

  # Now entropy values
  for (entropy in e_list) {
    column_name <- entropy$column
    entropy_value <- round(entropy$metric, display_precision)
    summary_df['entropy', column_name] <- entropy_value
  }

  # Clear out NAs for the sake of vanity
  summary_df[is.na(summary_df)] <- ''

  # Append original dataframe if requested
  if (append_input) {
    summary_df <- rbind(summary_df, df)
  }

  return(summary_df)
}


### Normalization and standardization of variables

#### *Normalización y estandarización de variables, tanto de manera individual como para el dataset completo. Esto solo debe ser aplicado a atributos que sean numéricos.*

# ----------------------------------------------------------------------------
#   Normalize dataframe by column
# ----------------------------------------------------------------------------

normalize_by_column <- function(df) {

  # Verify that the input is a dataframe
  if (!is.data.frame(df)) { stop('error, the input must be a dataframe') }

  # Verify that everything is numeric
  if (!all(sapply(df, is.numeric))) { stop('error, must be numeric') }

  # Apply the normalization to each column
  normalized_df <- as.data.frame(lapply(df, function(col) {
    (col - min(col, na.rm=T)) / (max(col, na.rm=T) - min(col, na.rm=T))
  }))

  return(normalized_df)
}


# ----------------------------------------------------------------------------
#   Standardize dataframe by column
# ----------------------------------------------------------------------------

standardize_by_column <- function(df) {

  # Verify that the input is a dataframe
  if (!is.data.frame(df)) { stop('error, the input must be a dataframe') }

  # Verify that everything is numeric
  if (!all(sapply(df, is.numeric))) { stop('error, must be numeric') }

  # Apply standardization to each column
  standardized_df <- as.data.frame(lapply(df, function(col) {
    (col - mean(col, na.rm=T)) / sd(col, na.rm=T)
  }))

  return(standardized_df)
}


### Filtering of variables based on the implemented metrics

#### *Filtrado de variables en base a las métricas implementadas. Es decir, partiendo de un dataset, obtener uno nuevo donde todas las variables cumplan los requisitos indicado (por ejemplo, una entropía superior a un cierto umbral).*

# ----------------------------------------------------------------------------
#   Select variables from dataframe based on rules and metrics
# ----------------------------------------------------------------------------

select_variables_by_metrics <- function(df, rules) {

  ### Create a new dataframe by selecting variables/columns/features from an
  ###  input dataframe whose metrics comply with a set of rules. Returns a
  ###  dataframe with a subset of the columns of the input. Mandatory
  ###  parameter rules is a list of lists defining rules for filtering. Each
  ###  rule contains:
  ###    metric (str):   'entropy', 'variance', 'auc'.
  ###    relation (str): '>=', '<=', '=='.
  ###    cutoff (float): comparison value for the metric selected'''

  # First extract metrics using extract_dataset_metrics()
  metrics <- extract_dataset_metrics(df)
  entropy_list <- metrics$entropy
  variance_list <- metrics$variance
  auc_list <- metrics$auc

  # Bundle metrics into a list with metric names as keys
  metrics_list <- list(
    entropy =  setNames(sapply(metrics$entropy,  `[[`, "metric"),
                        sapply(metrics$entropy,  `[[`, "column")  ),
    variance = setNames(sapply(metrics$variance, `[[`, "metric"),
                        sapply(metrics$variance, `[[`, "column")  ),
    auc =      setNames(sapply(metrics$auc,      `[[`, "metric"),
                        sapply(metrics$auc,      `[[`, "column")  )
  )

  # Initialize list of selected columns to all columns in the input df
  selected_cols <- colnames(df)  # column names should already be unique,
  # making this work like a set

  # Process each rule in our list
  for (rule in rules) {
    metric <- rule[[1]]
    relation <- rule[[2]]
    cutoff <- rule[[3]]

    metric_values <- metrics_list[[metric]]  # retrive, if present

    # New vector for columns that comply with the rule
    cols_complying_with_rule <- character(0)

    # Check each column and value according to the relation and cutoff
    for (col in names(metric_values)) {
      value <- metric_values[[col]]

      # Apply the specified relation
      if ((relation == ">=" && value >= cutoff) ||
          (relation == "<=" && value <= cutoff) ||
          (relation == "==" && value == cutoff)) {
        cols_complying_with_rule <- c(cols_complying_with_rule, col)
      }
    }

    # Update selected_cols with only those columns meeting the rule
    selected_columns <- intersect(selected_cols, cols_complying_with_rule)
  }

  # Return a dataframe, with only what we've selected
  return(df[, selected_columns, drop=FALSE])
}


### Calculation of the correlation/mutual information

#### *Cálculo de la correlación (información mutua en el caso de variables categóricas) por pares entre variables de un dataset. La función deberá considerar de que tipo es cada variable*

# ----------------------------------------------------------------------------
#   Correlation calculation, vector version
# ----------------------------------------------------------------------------

correlation <- function(x, y) {

  ### Calculate Pearson's correlation coefficient between the two numerical
  ###  input vectors x and y, passed as numpy arrays. Returns numerical
  ###  scalar.

  # All subsequent operations depend on equal length, so check to make sure
  #   the lengths of the x and y arrays match
  if (length(x) != length(y)) {
    stop("Error: Vectors x and y must be the same length.") }

  # (Pearson's) correlation coefficient (for samples or populations) is:
  #
  #       covariance(X,Y)                 sum( (X-Xbar)·(Y-Ybar) )
  # r = ——————————————————— = —————————————————————————————————————————————
  #     stddev(X)·stddev(Y)   sqrt(sum((X-Xbar)^2)) · sqrt(sum((Y-Ybar)^2))
  #
  # which for clarity and efficiency with numpy is better expressed as:
  #
  #                   sum_i=1_N( (x_i-xbar)·(y_i-ybar) )
  #   r = ——————————————————————————————————————————————————————————
  #       sqrt( sum_i=1_N( (x_i-xbar)^2) · sum_i=1_N( (y-ybar)^2 ) )
  #
  # Note: The above applies to both samples and populations--both the
  #         covariance and product of standard deviations actually have
  #         either N+1 or N in their denominators depending on whether they
  #         are samples, or entire populations, but in any case both
  #         denominators will be the same, and therefore cancel out.

  # Get mean values for x and y
  x_bar <- mean(x)
  y_bar <- mean(y)

  # Get the numerator of the covariance of x and y
  covariance_numerator <- sum((x - x_bar) * (y - y_bar))

  # Get the numerator of the product of the standard deviations of x and y
  stddevprod_numerator <- sqrt(sum((x - x_bar)^2) * sum((y - y_bar)^2))

  # Check the stddevprod_numerator (for the denominator),
  if (stddevprod_numerator == 0) {   # if 0, the variability is also 0,
    return(0)                      # and the last division would explode
  }

  # Calculate correlation coefficient
  correlation <- covariance_numerator / stddevprod_numerator

  return(correlation)
}


# ----------------------------------------------------------------------------
#   Mutual information calculation, vector version
# ----------------------------------------------------------------------------

mutual_information <- function(x, y) {

  ### Calculate mutual information coefficient between the two numerical
  ###  input vectors x and y, passed in as arrays. Returns numerical scalar.

  # Force inputs to be treated as factors
  x <- factor(x)
  y <- factor(y)

  # All subsequent operations depend on equal length, so check to make sure
  #   the lengths of the x and y arrays match
  if (length(x) != length(y)) {
    stop("Error: Vectors x and y must be of the same length.") }

  # Mutual information is defined as:
  #
  #                                  p(x,y)
  #   I = sum_X( sum_Y( p(x,y)·log(———————————) ) )
  #                                 p(x)·p(y)
  #

  # Get joint probabilities p_xy
  joint_table <- table(x, y)
  p_xy <- joint_table / sum(joint_table)

  # Calculate marginal probabilities for x and y
  p_x <- rowSums(p_xy)
  p_y <- colSums(p_xy)

  # Get expected probability matrix p_x * p_y with outer()
  p_x_X_p_y <- outer(p_x, p_y)

  # Calculate mutual information with:
  #
  #                       p_xy
  #   I = sum( p_xy·log(—————————) )
  #                      p_x·p_y
  #

  # These calculations were blowing up for any input type that was
  #  non-numeric, in spite of forcing the inputs to be factors. To avoid
  #  infinities and undefined results, we can only operation where the
  #  joint probability != 0. Let's get a list of indices where p_xy is
  #  non-zero, and use ONLY THOSE in our summation.
  nz <- p_xy > 0

  # Calculate mutual information using vectorized operations
  I <- sum(p_xy[nz] * log(p_xy[nz] / p_x_X_p_y[nz]))

  return(I)
}


# ----------------------------------------------------------------------------
#   Relationships with correlation and mutual information, dataframe version
# ----------------------------------------------------------------------------

column_relationships <- function(df) {
  ### Build tables exhibiting the relationships between variables in a
  ###  dataset, one for all numerical features using correlation, and another
  ###  for all categorical features using mutual information. Input is a
  ### data frame with arbitrary columns, returns two data frames.

  mi_vars <- c()  # List for categorical variables
  co_vars <- c()  # List for numerical variables

  # Check each column in the input, divide into columns to be included in
  #   categorical cross-comparison and numerical cross-comparison
  for (column in names(df)) {
    col_type <- class(df[[column]])
    if (col_type %in% c("character", "factor", "integer", "logical")) {
      mi_vars <- c(mi_vars, column)
    } else if (col_type == "numeric") {
      co_vars <- c(co_vars, column)
    } else {
      message(sprintf("'%s' has unsupported type %s, skipping.",
                      column, col_type))
    }
  }

  # Create the square correlation and association grid dataframes
  correlation_df <- as.data.frame(matrix(NA, nrow=length(co_vars),
                                         ncol=length(co_vars)  ))
  rownames(correlation_df) <- co_vars
  colnames(correlation_df) <- co_vars

  association_df <- as.data.frame(matrix(NA, nrow=length(mi_vars),
                                         ncol=length(mi_vars)  ))
  rownames(association_df) <- mi_vars
  colnames(association_df) <- mi_vars


  # Calculate mutual information for the relevant variables
  for (i in seq_along(mi_vars)) {
    for (j in seq(i, length(mi_vars))) {

      # Note that the mutual info H(x,x) SHOULD be the entropy,
      #  but it seems to be distinct from Shannon entropy which
      #  does not seem to give the same results

      I <- mutual_information(as.factor(df[[mi_vars[i]]]),
                              as.factor(df[[mi_vars[j]]]))
      association_df[i, j] <- I
      association_df[j, i] <- I
    }
  }

  # Calculate correlation for the relevant variables
  for (i in seq_along(co_vars)) {
    for (j in seq(i, length(co_vars))) {  # self-correlation is always 1
      r <- correlation(df[[co_vars[i]]],
                       df[[co_vars[j]]])
      correlation_df[i, j] <- r
      correlation_df[j, i] <- r
    }
  }

  return(list(correlation_df = correlation_df, association_df = association_df))
}


### Plots for the AUC and correlation/mutual information matrices

#### *Plots para el AUC y para las matrices de correlación/información mutua.*

# ----------------------------------------------------------------------------
#   Plot ROC (receiver operating characteristic) curve and display AUC
# ----------------------------------------------------------------------------

plot_roc_auc <- function(roc_fpr, roc_tpr, auc, title = 'ROC Curve') {

  graphics::plot.new()

  # Create a dataframe for plotting
  roc_data <- data.frame(FPR=roc_fpr, TPR=roc_tpr)

  # Plot the ROC curve with ggplot
  p <- ggplot(roc_data, aes(x=FPR, y=TPR)) +   # ROC curve itself
    geom_line(color="mediumblue", size=1) +
    geom_point(color="mediumblue") +
    geom_area(aes(y=TPR), fill="lightblue", alpha=0.3) + # the filling

    # Now we draw the typical straight line of random chance
    geom_abline(intercept=0, slope=1,
                color="red", linetype="dashed") +
    annotate("text", x=0.1, y=0.05,
             label = "Random Chance", color="red", size=4, hjust=0) +

    # Some labeling for the axes
    labs(x='False Positive Rate', y='True Positive Rate',
         title = paste0(title, ', AUC = ', auc)) +

    # We extend the limits to create space for points at 0
    xlim(-0.05, 1.05) + ylim(-0.05, 1.05) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))

  # Draw the plot
  print(p)
}


# ----------------------------------------------------------------------------
#   Visualizing relationships with correlation and/or mutual information
# ----------------------------------------------------------------------------

plot_relationships <- function(df, title='Heatmap') {
  ### Draw triangular heatmap to exhibit relationships between variables in
  ###  a dataset. Input is a square dataframe. Nothing is returned.

  if (nrow(df) == 0) { return() }  # Do nothing no data was passed to us

  df_as_mat <- as.matrix(df)  # Triangles are easier with matrices

  # If the map is symmetric about the diagonal, we hide the bottom right
  if(isTRUE(all.equal(df_as_mat, t(df_as_mat), tolerance = 1e-6))) {
    df_as_mat[lower.tri(df_as_mat)] <- NA
  }

  # Build some labels for the cells themselves
  cell_numbers <- as.data.frame(round(df_as_mat,3))
  cell_numbers[is.na(cell_numbers)] <- ''

  # A basic reasonable relative color pallete
  color_palette <- colorRampPalette(c('white', 'blue'))(100)

  # Plot a heatmap of the relationship using pheatmap
  pheatmap(df_as_mat,
           main = title,
           display_numbers = cell_numbers,
           fontsize = 10,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           color = color_palette,
           labels_row = colnames(df),
           labels_col = colnames(df),
           na_col = "white"
  )
}

#' Run wild clustered bootstrap for one bootby variable
#'
#' @param data Random data to perform wild clustered bootstrap on
#' @param model lm class model for bootstrapping on
#' @param x_interest X variable of interest, defaults to 'X'
#' @param clusterby String or vector of strings indicating variables to cluster by
#' @param boot_dist Vector of bootstrap distribution
#' @param boot_reps Integer for the number of bootstrap replications
#' @param bootby String or vector of strings indicating variables to bootstrap by
#' @param H0 Number indicating the null hypothesis, default is 0
#' @param cores Number of cores for parallel bootstrap
#' @return Vector with p-value for clusterby variable
#'
wild_cluster_boot_bootby <- function(data, model, x_interest, clusterby, boot_dist, boot_reps, bootby = clusterby, H0 = 0, cores = 1){

  #Run bootstrap function and save boot.out object
  boot.out <- do.call(what = wildclusterboot::wild_cluster_boot, args = as.list(environment()))

  #Use boot_p_val function to get vector of p-values and return
  p_vec <- wildclusterboot::boot_p_val(boot.out)

  #Turn p vector into list and add names
  p_list <- as.list(p_vec)
  names(p_list) <- clusterby

  return(p_list)

}

#' Run one replication of wild clustered bootstrap simulation
#'
#' @param dgp Function for generating random data
#' @param dgp_args Named list of arguments to be passed to dgp function
#' @param formula Formula for generating model from random data
#' @param x_interest X variable of interest, defaults to 'X'
#' @param clusterby String or vector of strings indicating variables to cluster by
#' @param boot_dist Vector of bootstrap distribution
#' @param boot_reps Integer for the number of bootstrap replications
#' @param bootby String or vector of strings indicating variables to bootstrap by
#' @param H0 Number indicating the null hypothesis, default is 0
#' @param cores Number of cores for parallel bootstrap
#' @param n Repition number
#' @return Named vector of p-values, with names representing combination of clusterby and bootby
#'
wild_cluster_boot_rep <- function(dgp, dgp_args, formula,
                                  x_interest, clusterby, boot_dist, boot_reps,
                                  bootby = clusterby, H0 = 0, cores = 1, progbar = NULL, n){

  #Create data using DGP and DGP arguments
  data <-  do.call(what = dgp, args = dgp_args)

  #Create model using data and formula
  model <- lm(data = data, formula = formula)

  clusterby_df <- data[clusterby]

  #Generate asymptotic t_values
  se_values <- sapply(X = clusterby_df, FUN = wildclusterboot::clustered_se, model = model)
  t_values <- coef(model)[x_interest]/se_values[x_interest,]

  #Generate asymptotic p_values
  df <- sapply(X = clusterby_df, FUN = function(x) length(unique(x)) - 1)
  asym_p_values <- 2*pt(t_values, df, lower.tail = FALSE)
  names(asym_p_values) <- paste('asym', names(asym_p_values), sep = '_')

  #Use wild_cluster_boot to get p-values for each bootby variable
  p_values <- sapply(X = bootby,
                   FUN = wild_cluster_boot_bootby,
                   data = data,
                   model = model,
                   clusterby = clusterby,
                   boot_dist = boot_dist,
                   boot_reps = boot_reps,
                   x_interest = x_interest,
                   H0 = H0,
                   cores = cores,
                   USE.NAMES = TRUE,
                   simplify = FALSE)

  if(!is.null(FALSE)){
    setTxtProgressBar(progbar, n)
  }

  #Flatten out list of p-values
  p_vec <- c(p_values, recursive = TRUE)

  #Replace periods with underscores in combined names
  names(p_vec) <- gsub('\\.', '_', names(p_vec))

  #Combine p-value vectors and return
  combined_p_vec <- append(p_vec, asym_p_values)
  return(combined_p_vec)

}

#' Run wild clustered bootstrap Monte Carlo and write to file
#'
#' @param dgp Function for generating random data
#' @param dgp_args Named list of arguments to be passed to dgp function
#' @param formula Formula for generating model from random data
#' @param Reps Number of total Monte Carlo repetitions
#' @param x_interest X variable of interest, defaults to 'X'
#' @param clusterby String or vector of strings indicating variables to cluster by
#' @param boot_dist Vector of bootstrap distribution
#' @param boot_reps Integer for the number of bootstrap replications
#' @param bootby String or vector of strings indicating variables to bootstrap by
#' @param H0 Number indicating the null hypothesis, default is 0
#' @param cores Number of cores for parallel bootstrap
#' @param progbar Boolean indicating whether to use progress bar
#' @param save_type Method for saving data, currently only supports db for SQLite
#' @param save_file String indicating file name to save to
#' @return Named matrix of p-values with row for eahch rep and columns for each combination of clusterby and bootby
#' @export
wild_cluster_boot_mc <- function(dgp, dgp_args, formula, reps,
                                 x_interest, clusterby, boot_dist, boot_reps,
                                 bootby,  H0 = 0, cores = 1, progbar = FALSE,
                                 save_type = 'db', save_file = 'wildclusterbootsim.db'){

  param_string <- paste(names(dgp_args), dgp_args, sep = ': ', collapse = ', ')
  print(paste('   Running replication on parameters:', param_string), quote = FALSE)

  #If progbar is true, create progress bar object
  if(progbar){

    progbar <- txtProgressBar(min = 0, max = reps, style = 3)

  }

  #Apply wild_cluster_boot_rep for rep replications
  p_values <- sapply(X = 1:reps, FUN = wild_cluster_boot_rep,
                     dgp = dgp, dgp_args = dgp_args, formula = formula,
                     x_interest = x_interest, clusterby = clusterby, boot_dist = boot_dist, boot_reps = boot_reps,
                     bootby = bootby, H0 = H0, cores = cores, progbar = progbar)

  p_df <- data.frame(t(p_values))

  all_df <- merge(p_df, dgp_args)

  #If saving to DB, save out to SQLite db
  if(save_type == 'db'){

    table <- 'p_values'
    con <-DBI::dbConnect(RSQLite::SQLite(), save_file)
    DBI::dbWriteTable(conn = con, name = table, value = data.frame(all_df), append = TRUE)
    DBI::dbDisconnect(conn = con)

  }

  return(t(p_values))

}


#' Run wild clustered bootstrap Monte Carlo on multiple combinations of parameters
#'
#' @param dgp Function for generating random data
#' @param dgp_args Named list of arguments to be passed to dgp function
#' @param formula Formula for generating model from random data
#' @param Reps Number of total Monte Carlo repetitions
#' @param x_interest X variable of interest, defaults to 'X'
#' @param clusterby String or vector of strings indicating variables to cluster by
#' @param boot_dist Vector of bootstrap distribution
#' @param boot_reps Integer for the number of bootstrap replications
#' @param bootby String or vector of strings indicating variables to bootstrap by
#' @param H0 Number indicating the null hypothesis, default is 0
#' @param cores Number of cores for parallel bootstrap
#' @param progbar Boolean indicating whether to use progress bar
#' @param save_type Method for saving data, currently only supports db for SQLite
#' @param save_file String indicating file name to save to
#' @return List of named matrices of p-values, with row for eahch rep and columns for each combination of clusterby and bootby
#' @export
muli_wild_cluster_boot_mc <- function(dgp, dgp_args, formula, reps,
                                      x_interest, clusterby, boot_dist, boot_reps,
                                      bootby = clusterby,  H0 = 0, cores = 1, progbar = FALSE,
                                      save_type = 'db', save_file = 'wildclusterbootsim.db'){

  #Expand dgp_args into all possible combinations
  parameter_df <- expand.grid(dgp_args)

  #Create list of lists out of dataframe
  parameter_lists <- split(x = parameter_df, f = seq(nrow(parameter_df)))

  #Run wild_cluster_boot_mc over
  lapply(X = parameter_lists, FUN = wild_cluster_boot_mc,
         dgp = dgp, formula = formula, reps = reps,
         x_interest = x_interest, clusterby = clusterby, boot_dist = boot_dist, boot_reps = boot_reps,
         bootby = bootby, H0 = H0, cores = cores, progbar = progbar)

}

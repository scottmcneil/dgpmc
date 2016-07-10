#' Run wild clustered bootstrap for one bootby variable
#'
#' @param bootby String or vector of strings indicating variables to bootstrap by
#' @param data Random data to perform wild clustered bootstrap on
#' @param model lm class model for bootstrapping on
#' @param clusterby String or vector of strings indicating variables to cluster by
#' @param boot_dist Vector of bootstrap distribution
#' @param boot_reps Integer for the number of bootstrap replications
#' @param H0 Number indicating the null hypothesis, default is 0
#' @param x_interest X variable of interest, defaults to 'X'
#' @return Vector with p-value for clusterby variable
#'
wildclustboot <- function(bootby, data, model, clusterby, boot_dist, boot_reps, x_interest, H0 = 0, cores = 1){

  #Set parallel option based on number of cores
  if(cores == 1){
    parallel <- 'no'
  } else {
    parallel <- 'snow'
  }

  #Create ran.gen and statistic functions
  ran.gen <- wildclusterboot::wild_clust_ran(model = model, x_interest = x_interest, bootby = bootby, boot_dist = boot_dist, H0 = H0)
  statistic <- wildclusterboot::wild_clust_statistic(model = model, x_interest = x_interest, clusterby = clusterby, H0 = H0)

  #Run bootstrap function and save boot.out object
  boot.out <- boot::boot(data = data, statistic = statistic, R = boot_reps, sim = 'parametric', ran.gen = ran.gen, parallel = parallel, ncpus = cores)

  #Use boot_p_val function to get vector of p-values and return
  p_vec <- wildclusterboot::boot_p_val(boot.out)
  p_list <- as.list(p_vec)
  names(p_list) <- clusterby
  return(p_list)

}

#' Run one replication of wild clustered bootstrap simulation
#'
#' @param dgp Function for generating random data
#' @param dgp_args Named list of arguments to be passed to dgp function
#' @param formula Formula for generating model from random data
#' @param bootby String or vector of strings indicating variables to bootstrap by
#' @param clusterby String or vector of strings indicating variables to cluster by
#' @param boot_dist Vector of bootstrap distribution
#' @param boot_reps Integer for the number of bootstrap replications
#' @param H0 Number indicating the null hypothesis, default is 0
#' @param x_interest X variable of interest, defaults to 'X'
#' @return Named vector of p-values, with names representing combination of clusterby and bootby
#'
wildclustboot_rep <- function(dgp, dgp_args, formula, bootby, clusterby, boot_dist, boot_reps, x_interest, H0 = 0, cores = 1){

  #Create data using DGP and DGP arguments
  data <-  do.call(what = dgp, args = dgp_args)

  #Create model using data and formula
  model <- lm(data = data, formula = formula)

  #Use wildclustboot to get p-values for each bootby variable
  p_vals <- sapply(X = bootby,
                   FUN = wildclustboot,
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

  #Flatten out list of p-values
  p_vec <- c(p_vals, recursive = TRUE)

  #Replace periods with underscores in combined names
  names(p_vec) <- gsub('\\.', '_', names(p_vec))
  return(p_vec)

}

#' Run one replication of wild clustered bootstrap simulation
#'
#' @param dgp Function for generating random data
#' @param dgp_args Named list of arguments to be passed to dgp function
#' @param formula Formula for generating model from random data
#' @param bootby String or vector of strings indicating variables to bootstrap by
#' @param clusterby String or vector of strings indicating variables to cluster by
#' @param boot_dist Vector of bootstrap distribution
#' @param boot_reps Integer for the number of bootstrap replications
#' @param H0 Number indicating the null hypothesis, default is 0
#' @param x_interest X variable of interest, defaults to 'X'
#' @return Named matrix of p-values, with row for each rep and names representing combination of clusterby and bootby
#' @export
wildclustboot_mc <- function(reps, dgp, dgp_args, formula, bootby, clusterby, boot_dist, boot_reps, x_interest, H0 = 0, cores = 1){

  p_values <- replicate(n = reps, expr = wildclustboot_rep(dgp = dgp,
                                                           dgp_args = dgp_args,
                                                           formula = formula,
                                                           bootby = bootby,
                                                           clusterby = clusterby,
                                                           boot_dist = boot_dist,
                                                           boot_reps = boot_reps,
                                                           x_interest = x_interest,
                                                           H0 = H0,
                                                           cores = cores))
  return(t(p_values))

}

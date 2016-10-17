#' Run one Monte Carlo repetition
#'
#' @param stat_func Function that takes a dataframe from rand_func and returns a vector of statistics
#' @param stat_args List of arguments to be passed to stat_func
#' @param rand_func Function that generates a dataframe of random data given rand_args
#' @param rand_args List of arguments to be passed to stat_func
#' @param rep Integer for working with multicore
#' @return Vector of numerical statistics
#'
monte_carlo_stat <- function(stat_func, stat_args, rand_func, rand_args, rep){

  #Generate set of random data, formatted as named list
  data <- list(data = do.call(what = rand_func, args = rand_args))

  #Run statistical function using concat of data and stat_args
  stat_vec <- do.call(what = stat_func, args = c(data, stat_args))

  #Switch to matrix for maintaining dimensions
  stat_mat <- matrix(stat_vec)

  return(stat_mat)

}

#' Run block of Monte Carlo repetitions
#'
#' @param block Integer indicating number of repetitions for block
#' @param cumulative Integer indicating total repetitions completed including this block
#' @param stat_func Function that takes a dataframe from rand_func and returns a vector of statistics
#' @param stat_args List of arguments to be passed to stat_func
#' @param rand_func Function that generates a dataframe of random data given rand_args
#' @param rand_args List of arguments to be passed to stat_func
#' @param pb Progress bar object to be incremented
#' @param seed List of seeds to be passed to each replication, will be NULL if lecuyer is set to FALSE
#' @param cores Integer indicating number of cores to use
#' @return Vector of numerical statistics
#'
monte_carlo_block <- function(block, cumulative, stat_func, stat_args, rand_func, rand_args, pb, seed, cores){

  if(cores == 1){

    #If just running on one core, run non-paralell process
    stat_mat <- monte_carlo_stat(stat_func = stat_func, stat_args = stat_args, rand_func = rand_func, rand_args = rand_args)

  } else {

    # If seed is passed from above, set it using L'Ecuyer
    if(!is.null(seed)){

      .GlobalEnv$.Random.seed <- if(block > 1) seed else parallel::nextRNGStream(seed)

    }

    #If running in paralell, run and create list of statistics
    stat_mats <- parallel::mcmapply(FUN = monte_carlo_stat, rep = 1:block,  mc.cores = cores, SIMPLIFY = FALSE,
                                    mc.set.seed = TRUE, MoreArgs = list(stat_func = stat_func, stat_args = stat_args,
                                                                         rand_func = rand_func, rand_args = rand_args))

    stat_mat <- do.call(what = cbind, args = stat_mats)
  }

  #Increment pb if not NULL
  if(class(pb) == 'txtProgressBar'){
    setTxtProgressBar(pb = pb, value = cumulative)
  }

  return(stat_mat)

}

#' Helper function for L'Ecuyer seeds
#'
#' @param blocks List of integers indicating number of repetitions for each block
#' @param reps Integer indicating the total Monte Carlo repetitions
#' @param seed Initial seed either from current state or passed explicitely
#' @param lecuyer Boolean indicating whether to use L'Ecuyer-style seed for multicore reproducibility
#' @return Vector or matrix of numerical statistics
#'
multicore_seeds <- function(blocks, reps, seed, lecuyer){

  seeds <- if(lecuyer){

    #Create reps long list of seeds
    seeds <- multicore_seeds_recur(n = 1, reps = reps, seed = seed)

    #Create vector of block ids
    block_ids <- rep(1:length(blocks), times = blocks)

    #Create list of lists of L'Ecuyer seeds
    all_seeds <- split(x = seeds, f = block_ids)
    lapply(all_seeds, function(seed_list) seed_list[[1]])

  } else{

    #Create list of lists of NULL values to be passed as seeds when L'Ecuyer is FALSE
    rep(list(NULL), length(blocks))

  }

  return(seeds)

}

#' Recursive function for creating list of needed L'Ecuyer seeds
#'
#' @param n Integer indicating the current repetition
#' @param reps Integer indicating the total number of reps
#' @param seed Initial seed either from current state or passed explicitely
#' @param lecuyer Boolean indicating whether to use L'Ecuyer-style seed for multicore reproducibility
#' @return Vector or matrix of numerical statistics
#'
multicore_seeds_recur <- function(n, reps, seed){

  if(n >= reps){
    return(list(seed))
  } else{
    return(c(list(seed), multicore_seeds_recur(n + 1, reps, parallel::nextRNGStream(seed))))
  }

}

#' Repeat Monte Carlo repetitions
#'
#' @param reps Integer indicating the total Monte Carlo repetitions
#' @param stat_func Function that takes a dataframe from rand_func and returns a vector of statistics
#' @param stat_args List of arguments to be passed to stat_func
#' @param rand_func Function that generates a dataframe of random data given rand_args
#' @param rand_args List of arguments to be passed to stat_func
#' @param progbar Progress bar object to be incremented
#' @param lecuyer Boolean indicating whether to use L'Ecuyer-style seed for multicore reproducibility
#' @param cores Integer indicating number of cores to use
#' @return Vector or matrix of numerical statistics
#'
monte_carlo_recur <- function(reps, stat_func, stat_args, rand_func, rand_args, progbar, lecuyer, cores){

  #Create pb object if TRUE
  if(progbar == TRUE){
    pb <- txtProgressBar(min = 0, max = reps, style = 3)
  } else {
    pb <- NULL
  }

  if(class(pb) == 'txtProgressBar' | cores > 1){

    #If using progress bar, chunk up repetitions by core
    minus_remain <- rep(cores, reps %/% cores)
    blocks <- if(reps %% cores == 0) minus_remain else c(minus_remain, reps %% cores)
    cumulative <- cumsum(blocks)

    #Create list of lists of seeds
    seed <- .GlobalEnv$.Random.seed
    seeds <- multicore_seeds(blocks = blocks, reps = reps, seed = seed, lecuyer = lecuyer)

    #Then run monte_carlo_block
    stat_mats <- mapply(monte_carlo_block, block = blocks, seed = seeds,
                        cumulative = cumulative, SIMPLIFY = FALSE, MoreArgs = list(stat_func = stat_func, stat_args = stat_args,
                                                                                   rand_func = rand_func, rand_args = rand_args,
                                                                                   pb = pb, cores = cores))

    #Combine into matrix
    stat_mat <- do.call(cbind, stat_mats)

  } else{

    #Otherwise run replications sequentially
    stat_mat <- mapply(FUN = monte_carlo_stat, rep = 1:reps, MoreArgs = list(stat_func = stat_func, stat_args = stat_args,
                                                                             rand_func = rand_func, rand_args = rand_args))

  }

  return(stat_mat)

}

#' Run random data Monte Carlo
#'
#' @param reps Integer indicating the total Monte Carlo repetitions
#' @param stat_func Function that takes a dataframe from rand_func and returns a vector of statistics
#' @param stat_args List of arguments to be passed to stat_func
#' @param rand_func Function that generates a dataframe of random data given rand_args
#' @param rand_args List of arguments to be passed to stat_func
#' @param progbar Boolean for whether progress bar should be used
#' @param cores Integer indicating number of cores to use
#' @param seed Integer indicating the initial seed to use, default is Mersenne Twister
#' @param lecuyer Boolean indicating whether to use L'Ecuyer-style seed for multicore reproducibility
#' @return Vector or matrix of numerical statistics
#' @export
#'
dgpmc <- function(reps, stat_func, stat_args, rand_func, rand_args, names = NULL, progbar = FALSE, seed = NULL, lecuyer = FALSE, cores = 1){

  if(!is.null(seed)){

    seed_type <- if(lecuyer) "L'Ecuyer-CMRG" else 'Mersenne-Twister'

    set.seed(seed, seed_type)

  } else if(lecuyer){

    stop('You must choose a seed if setting lecuyer to TRUE')

  }

  #Generate stat_mat from monte_carlo_recur
  stat_mat <- monte_carlo_recur(reps = reps, stat_func = stat_func, stat_args = stat_args,
                                rand_func = rand_func, rand_args = rand_args,
                                progbar = progbar, lecuyer = lecuyer, cores = cores)

  #Check if vector or matrix and convert to dataframe
  if(class(stat_mat) == 'matrix'){

    stat_df <- data.frame(t(stat_mat))

  } else{

    stat_df <- data.frame(stat_mat)

  }

  #Set names of dataframe
  names(stat_df) <- names

  return(list(stat = stat_df, stat_func = stat_func, stat_args = stat_args, rand_func = rand_func, rand_args = rand_args))

}

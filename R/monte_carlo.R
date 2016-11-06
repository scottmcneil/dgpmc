#' Run one Monte Carlo repetition
#'
#' @param stat_func Function that takes a dataframe from rand_func and returns a vector of statistics
#' @param stat_args List of arguments to be passed to stat_func
#' @param rand_func Function that generates a dataframe of random data given rand_args
#' @param rand_args List of arguments to be passed to stat_func
#' @param seed Seed which is either null or sets the L'Ecuyer seed for multicore
#' @return Vector of numerical statistics
#'
monte_carlo_stat <- function(stat_func, stat_args, rand_func, rand_args, seed = NULL, rep){

  #If using L'Ecuyer seeds, set seed:
  if(!is.null(seed)){
    .GlobalEnv$.Random.seed <- seed
  }

  #Generate set of random data, formatted as named list
  data <- list(data = do.call(what = rand_func, args = rand_args))

  #Run statistical function using concat of data and stat_args
  stat_vec <- do.call(what = stat_func, args = c(data, stat_args))

  #Switch to matrix for maintaining dimensions
  stat_mat <- matrix(stat_vec, ncol = 1)

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
#' @param seed_ffwd Integer indicating the number of L'Ecuyer streams to move ahead, for reproducibility purposes
#' @param cores Integer indicating number of cores to use
#' @return Vector of numerical statistics
#'
monte_carlo_block <- function(reps, cumulative, seed, stat_func, stat_args, rand_func, rand_args, pb, cores, seed_ffwd){


  if(seed_ffwd > 0){

    #Set ffwd block
    block <- min(seed_ffwd, ceiling(sqrt(reps)))

    #Generate block seeds
    seeds <- multicore_seeds_recur(1, block, seed)

    #Pass seeds to next block
    return(monte_carlo_block(reps = reps, cumulative = cumulative, seed = seeds[[length(seeds)]],
                             stat_func = stat_func, stat_args = stat_args, rand_func = rand_func, rand_args = rand_args,
                             pb = pb, cores = cores, seed_ffwd = seed_ffwd - block))

  } else{

    #Set block number
    block <- min(reps - cumulative, max(cores, ceiling(reps/100)))

    #Create vector of seeds
    seeds <- if(!is.null(seed)) multicore_seeds_recur(1, block, seed) else rep(list(NULL), block)

    stat_mat <- if(cumulative < seed_ffwd){

      NULL

    } else if(cores == 1 & is.null(seed)){

      #If just running on one core with no passed in seeds, run non-paralell process
      stat_mats <- mapply(FUN = monte_carlo_stat, seed = seeds,  SIMPLIFY = FALSE, MoreArgs = list(stat_func = stat_func, stat_args = stat_args,
                                                                                                    rand_func = rand_func, rand_args = rand_args))

      do.call(what = cbind, args = stat_mats)

    } else {

      #If running in paralell, run and create list of statistics
      stat_mats <- parallel::mcmapply(FUN = monte_carlo_stat, seed = seeds,  mc.cores = cores, SIMPLIFY = FALSE,
                                      mc.set.seed = FALSE, MoreArgs = list(stat_func = stat_func, stat_args = stat_args,
                                                                           rand_func = rand_func, rand_args = rand_args))

      #Column bind all stat_mats
      do.call(what = cbind, args = stat_mats)

    }

    #Increment pb if not NULL
    if(class(pb) == 'txtProgressBar'){
      setTxtProgressBar(pb = pb, value = cumulative + block)
    }

    if(block + cumulative == reps){

      #If at end of reps, return current stat_mat
      return(stat_mat)

    } else{

      #Return column bind of current stat_mat and recursion
      return(cbind(stat_mat, monte_carlo_block(reps = reps, cumulative = cumulative + block, seed = seeds[[length(seeds)]],
                                               stat_func = stat_func, stat_args = stat_args, rand_func = rand_func, rand_args = rand_args,
                                               pb = pb, cores = cores, seed_ffwd = seed_ffwd)))
    }
  }

}

#' Recursive function for creating list of needed L'Ecuyer seeds
#'
#' @param n Integer indicating the current repetition
#' @param block Integer indicating the total number of reps
#' @param seed Initial seed either from current state or passed explicitely
#' @return Seed for stream advanced block number times
#'
multicore_seeds_recur <- function(n, block, seed){

  if(n >= block){
    return(seed)
  } else{
    return(multicore_seeds_recur(n + 1, block, parallel::nextRNGStream(seed)))
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
#' @param seed_ffwd Integer indicating the number of L'Ecuyer streams to move ahead, for reproducibility purposes
#' @param cores Integer indicating number of cores to use
#' @return Vector or matrix of numerical statistics
#'
monte_carlo_recur <- function(reps, stat_func, stat_args, rand_func, rand_args, progbar, lecuyer, cores, seed_ffwd){

  #Create pb object if TRUE
  if(progbar == TRUE){
    pb <- txtProgressBar(min = 0, max = reps, style = 3)
  } else {
    pb <- NULL
  }

  stat_mat <- if(class(pb) == 'txtProgressBar' | cores > 1){

    #If using L'Ecuyer seeds, save current random seed
    seed <- if(lecuyer) .GlobalEnv$.Random.seed else NULL
    cumulative <- 0

    #Then run monte_carlo_block
    monte_carlo_block(reps = reps, cumulative = cumulative, seed = seed, stat_func = stat_func, stat_args = stat_args,
                                   rand_func = rand_func, rand_args = rand_args, pb = pb, cores = cores, seed_ffwd = seed_ffwd)

  } else{

    #Otherwise run replications sequentially
    mapply(FUN = monte_carlo_stat, rep = 1:reps, MoreArgs = list(stat_func = stat_func, stat_args = stat_args,
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
#' @param seed Integer indicating the initial seed to use, default is Mersenne Twister
#' @param lecuyer Boolean indicating whether to use L'Ecuyer-style seed for multicore reproducibility
#' @param seed_ffwd Integer indicating the number of L'Ecuyer streams to move ahead, for reproducibility purposes
#' @param cores Integer indicating number of cores to use
#' @param email List of options to pass to mailR send.mail function to send email alert when simulation completed
#' @return Vector or matrix of numerical statistics
#' @export
#'
dgpmc <- function(reps, stat_func, stat_args, rand_func, rand_args, names = NULL, progbar = FALSE,
                  seed = NULL, lecuyer = FALSE, seed_ffwd = 0, cores = 1, email = NULL){

  if(!is.null(seed)){

    seed_type <- if(lecuyer) "L'Ecuyer-CMRG" else 'Mersenne-Twister'

    set.seed(seed, seed_type)

  } else if(lecuyer){

    stop('You must choose a seed if setting lecuyer to TRUE')

  }

  #Generate stat_mat from monte_carlo_recur
  stat_mat <- monte_carlo_recur(reps = reps, stat_func = stat_func, stat_args = stat_args,
                                rand_func = rand_func, rand_args = rand_args,
                                progbar = progbar, lecuyer = lecuyer, cores = cores, seed_ffwd = seed_ffwd)

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

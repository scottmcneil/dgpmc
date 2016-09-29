#' Run one Monte Carlo repetition
#'
#' @param rep Integer indicating number of repetitions
#' @param stat_func Function that takes a dataframe from rand_func and returns a vector of statistics
#' @param stat_args List of arguments to be passed to stat_func
#' @param rand_func Function that generates a dataframe of random data given rand_args
#' @param rand_args List of arguments to be passed to stat_func
#' @return Vector of numerical statistics
#'
monte_carlo_stat <- function(rep, stat_func, stat_args, rand_func, rand_args){

  #Generate set of random data
  data <- do.call(what = rand_func, args = rand_args)

  stat_args[['data']] <- data

  #Run statistical function
  stat_vec <- do.call(what = stat_func, args = stat_args)

  return(stat_vec)

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
#' @param cores Integer indicating number of cores to use
#' @return Vector of numerical statistics
#'
monte_carlo_block <- function(block, cumulative, stat_func, stat_args, rand_func, rand_args, pb, cores){

  stat_mat <- parallel::mcmapply(FUN = monte_carlo_stat, rep = 1:block, mc.cores = cores,
                                 MoreArgs = list(stat_func = stat_func, stat_args = stat_args,
                                                 rand_func = rand_func, rand_args = rand_args))

  #Increment pb if not NULL
  if(class(pb) == 'txtProgressBar'){
    setTxtProgressBar(pb = pb, value = cumulative)
  }

  return(stat_mat)

}

#' Repeat Monte Carlo repetitions
#'
#' @param reps Integer indicating the total Monte Carlo repetitions
#' @param stat_func Function that takes a dataframe from rand_func and returns a vector of statistics
#' @param stat_args List of arguments to be passed to stat_func
#' @param rand_func Function that generates a dataframe of random data given rand_args
#' @param rand_args List of arguments to be passed to stat_func
#' @param progbar Progress bar object to be incremented
#' @param cores Integer indicating number of cores to use
#' @return Vector or matrix of numerical statistics
#'
monte_carlo_recur <- function(reps, stat_func, stat_args, rand_func, rand_args, progbar, cores){

  #Create pb object if TRUE
  if(progbar == TRUE){
    pb <- txtProgressBar(min = 0, max = reps, style = 3)
  } else {
    pb <- NULL
  }

  if(class(pb) == 'txtProgressBar'){

    #If using progress bar, chunk up repetitions by core
    minus_remain <- rep(cores, reps %/% cores)
    blocks <- if(reps %% cores == 0) minus_remain else c(minus_remain, reps %% cores)
    cumulative <- cumsum(blocks)

    #Then run monte_carlo_block
    stat_mat <- mapply(monte_carlo_block, block = blocks, cumulative = cumulative, MoreArgs = list(stat_func = stat_func, stat_args = stat_args,
                                                                                                   rand_func = rand_func, rand_args = rand_args,
                                                                                                   pb = pb))
  } else{

    #Otherwise run replications based on total cores selected
    stat_mat <- parallel::mcmapply(FUN = monte_carlo_stat, rep = 1:block, mc.cores = cores,
                                   MoreArgs = list(stat_func = stat_func, stat_args = stat_args,
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
#' @return Vector or matrix of numerical statistics
#' @export
#'
dgpmc <- function(reps, stat_func, stat_args, rand_func, rand_args, names = NULL, progbar = FALSE, cores = 1){

  #Generate stat_mat from monte_carlo_recur
  stat_mat <- monte_carlo_recur(reps = reps, stat_func = stat_func, stat_args = stat_args,
                                rand_func = rand_func, rand_args = rand_args, progbar = progbar)

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

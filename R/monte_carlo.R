#' Run one Monte Carlo repetition
#'
#' @param rep Integer indicating the current Monte Carlo repetition
#' @param stat_func Function that takes a dataframe from rand_func and returns a vector of statistics
#' @param stat_args List of arguments to be passed to stat_func
#' @param rand_func Function that generates a dataframe of random data given rand_args
#' @param rand_args List of arguments to be passed to stat_func
#' @param progbar Progress bar object to be incremented
#' @return Vector of numerical statistics
#'
monte_carlo_stat <- function(rep, stat_func, stat_args, rand_func, rand_args, progbar){

  #Generate set of random data
  data <- do.call(what = rand_func, args = rand_args)

  stat_args[['data']] <- data

  #Run statistical function
  stat_vec <- do.call(what = stat_func, args = stat_args)

  #Increment progbar if not NULL
  if(class(progbar) == 'txtProgressBar'){
    setTxtProgressBar(pb = progbar, value = rep)
  }

  return(stat_vec)

}

#' Repeat Monte Carlo repetitions
#'
#' @param reps Integer indicating the total Monte Carlo repetitions
#' @param stat_func Function that takes a dataframe from rand_func and returns a vector of statistics
#' @param stat_args List of arguments to be passed to stat_func
#' @param rand_func Function that generates a dataframe of random data given rand_args
#' @param rand_args List of arguments to be passed to stat_func
#' @param progbar Progress bar object to be incremented
#' @return Vector or matrix of numerical statistics
#'
monte_carlo_recur <- function(reps, stat_func, stat_args, rand_func, rand_args, progbar){

  #Run monte_carlo_stat reps times
  stat_mat <- sapply(1:reps, monte_carlo_stat, stat_func = stat_func, stat_args = stat_args,
                     rand_func = rand_func, rand_args = rand_args, progbar = progbar)

  return(stat_mat)

}

#' Run random data Monte Carlo
#'
#' @param reps Integer indicating the total Monte Carlo repetitions
#' @param stat_func Function that takes a dataframe from rand_func and returns a vector of statistics
#' @param stat_args List of arguments to be passed to stat_func
#' @param rand_func Function that generates a dataframe of random data given rand_args
#' @param rand_args List of arguments to be passed to stat_func
#' @param progbar Progress bar object to be incremented
#' @return Vector or matrix of numerical statistics
#' @export
#'
dgpmc <- function(reps, stat_func, stat_args, rand_func, rand_args, names = NULL, progbar = FALSE){

  #Create progbar object if TRUE
  if(progbar){
    progbar <- txtProgressBar(min = 0, max = reps, style = 3)
  }

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

#' Checks the number of rows already in a table with given id
#'
#' @param db Database to check values in
#' @param table Table to check variable/value in
#' @param id ID value to look for
#' @param id_name Name of ID variable, default is "id"
#' @return Number of rows in table with id
#'
check_completed <- function(db, table, id, id_name = 'id'){

  #Create SQLite connection
  con <- con <- DBI::dbConnect(RSQLite::SQLite(), db)

  #If table exists, read table to check number of rows with id
  if(dbExistsTable(conn = con, name = table) == 1){

    query <- paste('SELECT COUNT(*) FROM', table, 'WHERE', id_name, '=', id)
    completed <- DBI::dbFetch(DBI::dbSendQuery(conn = con, statement = query))

  } else{

    #Otherwise set to 0
    completed <- 0

  }

  #Disconnect from DB
  DBI::dbDisconnect(conn = con)

  return(completed)

}

#' Writes matrix of p-values to database
#'
#' @param p_vals Named matrix or dataframe of p-values
#' @param db Name of SQLite database
#' @param table Name of SQLite database
#' @param id ID value to be added to each row
#' @param id_name Name of ID variable, default is "id"
#'
db_write_pvals <- function(p_vals, db, table, id, id_name = 'id'){

  #Coerce p-values to dataframe and add ID column
  data <- data.frame(p_vals)
  data[id_name] <- id

  #Create connection and write table
  con <-DBI::dbConnect(RSQLite::SQLite(), db)
  DBI::dbWriteTable(conn = con, name = table, value = data, append = TRUE)

  #Close connection
  DBI::dbDisconnect(conn = con)

}

#' Runs monte carlo and write results to database
#'
#' @param db Name of SQLite database
#' @param table Name of SQLite database
#' @param id ID value to be added to each row
#' @param id_name Name of ID variable, default is "id"
#' @param parameters List of DGP and bootstrap parameters
#' @param reps Number of Monte Carlo replication
#' @param dgp DGP function
#' @param formula Formula for generating model from random data
#' @param bootby String or vector of strings indicating variables to bootstrap by
#' @param clusterby String or vector of strings indicating variables to cluster by
#' @param H0 Number indicating the null hypothesis, default is 0
#' @param x_interest X variable of interest, defaults to 'X'
#' @param cores Number of cores for parallel bootstrap
#' @export
#'
db_write_mc <- function(parameters, db, table, id_name = 'id', reps, dgp, formula, bootby, clusterby, x_interest, H0 = 0, cores = 1){

  boot_dist <- parameters['boot_dist']
  boot_reps <- parameters['boot_reps']
  id <- parameters[id_name]
  dgp_args <- parameters[!parameters %in% c('boot_dist', 'boot_reps', id_name)]

  completed <- check_completed(db = db, table = table, id = id, id_name = id_name)

  if(reps < completed){

    reps <- reps - completed

    p_vals <- wildclustboot_mc(reps = reps,
                               dgp = dgp,
                               dgp_args = dgp_args,
                               formula = formula,
                               bootby = bootby,
                               clusterby = clusterby,
                               boot_dist = boot_dist,
                               boot_reps = boot_reps,
                               x_interest = x_interest,
                               H0 = H0,
                               cores = cores)

    db_write_pvals(p_vals = p_vals,
                   db = db,
                   table = table,
                   id = id,
                   id_name = id_name)

  }


}

#' Writes parameter table to DB
#'
#' @param db Name of SQLite database
#' @param table Name of SQLite database
#' @param parameters_list Named list of parameters
#' @export
#'
db_write_param <- function (db, table, parameters_list){

  #Cartesian join parameters and add id
  parameters <- expand.grid(parameter_list)
  parameters['id'] <- rownames(parameters)

  #Open connection and write table
  con <- con <- DBI::dbConnect(RSQLite::SQLite(), db)
  DBI::dbWriteTable(conn = con, name = 'parameters', value = parameters, overwrite = TRUE)
  DBI::dbDisconnect(conn = con)

}

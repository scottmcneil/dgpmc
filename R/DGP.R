#' Generate autocorrelated numbers by group
#'
#' @param t The number of time periods
#' @param G The number of groups
#' @param rho Autocorrelation displacement parameter
#' @param mu Mean for random normal distribution
#' @param sigma Standard deviation for ranfom normal distribution
#' @return Matrix with autocorrelated, t and G columns
#'
autocorr_group <- function(t, G, rho = 1, mu = 0, sigma = 1){

  #Generate matrix of random normal variables
  rand_mat <- matrix(data = rnorm(n = t*G), nrow = t, ncol = G)

  #Generate vector of rho^t and repeat t times
  rho_vec <- rep(x = rho ^ as.numeric(0:t), times = t)

  #Roll up rho vec into matrix and 0 upper triangle
  rho_mat <- matrix(data = rho_vec, nrow = t, ncol = t)
  rho_mat[upper.tri(x = rho_mat)] <- 0

  #Create matrix of u multiplying rho and ran mats, then unroll unto vector
  u_vec <- as.vector(rho_mat %*% rand_mat)

  #Create t and G indices
  G_vec <- rep(x = 1:G, each = t)
  t_vec <- rep(x = 1:t, times = G)

  #Cbind vectors and return matrix
  u_mat <- cbind(t = t_vec, G = G_vec, u = u_vec)
  return(u_mat)

}

#'Generates t*sum(ng) or t*G*ng random normal variables indexed by t and G
#'
#' @param t The number of time periods for all G
#' @param G The number of groups in each time period
#' @param ng Obs per group, either int or G-length vector
#' @param mu Mean for random normal distribution
#' @param sigma Standard deviation for ranfom normal distribution
#' @return Matrix with random variable, t and G columns
#' @export
autocorr_dgp <- function(t, G, ng, rho = 1, lambda = 1, gamma = 1, mu = 0, sigma = 1){

  #Return error if ng not integer or G length vector
  if (length(ng) != 1 & length(ng) != G){
    stop('ng must be length 1 or length G')
  }

  time_level_data <- cbind(t = 1:t, v = rnorm(n = t, mean = mu, sd = sigma))
  group_level_data <- cbind(G = 1:G, w = rnorm(n = G, mean = mu, sd = sigma))
  auto_corr_data <- autocorr_group(t = t, G = G, mu = mu, sigma = sigma, rho = rho)

  G_levels <- rep(x = 1:G, each = t)
  t_levels <- rep(x = 1:t, times = G)

  #Check if ng integer or vector
  if (length(ng) > 1){
    #Expand ng vector to be t*G length
    ng_vec <- rep(x = ng, each = t)

    #Expand G, t and autocorr vectors by ng vector
    t_ind <- rep(x = t_levels, times = ng_vec)
    G_ind <- rep(x = G_levels, times = ng_vec)
    i <- t * sum(ng)

  } else {
    #Expand G, t and autocorr vectors by ng
    t_ind <- rep(x = t_levels, each = ng)
    G_ind <- rep(x = G_levels, each = ng)
    i <- t * G * ng
  }

  ind_level_data <- cbind(t = t_ind, G = G_ind, e = rnorm(n = i, mean = mu, sd = sigma))
  full_data <- data.frame(ind_level_data)
  full_data['tG'] <- paste(full_data[,'t'], full_data[,'G'], sep = '-')
  full_data['i'] <- row.names(full_data)

  full_data <- merge(x = full_data, y = time_level_data, all.x = TRUE)
  full_data <- merge(x = full_data, y = group_level_data, all.x = TRUE)
  full_data <- merge(x = full_data, y = auto_corr_data, all.x = TRUE)

  full_data['Y'] <- full_data[, 'v'] + full_data[, 'w'] + full_data[, 'u'] + full_data[, 'e']
  full_data['X'] <- full_data[, 'v'] + full_data[, 'w']
  final_data <- full_data[c('i', 'G', 'tG', 'Y', 'X')]
  return(final_data)

}

#'Generates G*ng random normal variables indexed by G
#'
#' @param G The number of groups in each time period
#' @param ng Obs per group, either int or G-length vector
#' @return Dataframe with X, Y and clusterby random variables
#' @export
simple_DGP <- function(ng,G){

  #Create clusterby variable
  clusterby <- as.integer(0:(ng*G-1)/ng+1)

  #Create individual level random data
  zig <- rnorm(n = G*ng,mean = 0,sd = 1)
  eig <- rnorm(n = G*ng,mean = 0,sd = 1)

  #Create group level random data
  zg <- rep(rnorm(n = G,mean = 0,sd = 1),each=ng)
  eg <- rep(rnorm(n = G,mean = 0,sd = 1),each=ng)

  #Combine data
  X <- zg + zig
  Y <- zg + zig + eg + eig

  return(data.frame(Y = Y, X = X, clusterby = clusterby))
}

#' Generate multi-way grouped random data
#'
#' @param dims Integer indicating the number of group dimensions
#' @param groups Vector of dims length with integers indicating total groups for each group dimension
#' @param rho Float indicating portion of first dimension groups for which cells wil have theta observations instead of 1
#' @param theta Integer indicating number of members of larger groups set out in rho
#' @param dim_names Character vector of names for dimensions, default is starting with H and working backward alphabetically
#' @param heterosked Boolean indicating whether to include heteroskedasticity (currently calculated by making sd(e_i) = 9*W^2)
#' @return dataframe of group-correlated random values
#' @export
multiway_DGP <- function(num_dims, groups, rho = 0, theta = 1, dim_names = rev(LETTERS)[1:num_dims+18], heterosked = FALSE){

  groups_list <- mapply(function(dim_name, groups) paste0(dim_name, 1:groups),
                        dim_name = dim_names, groups = groups, SIMPLIFY = FALSE)

  if(rho != 0){
    all_groups <- multiway_group_num(groups_list = groups_list, rho = rho, theta = theta)
  } else {
    all_groups <- groups_list
  }

  #Get combinatorial dataframe of groups from different dimensions
  group_combinations <- expand.grid(all_groups)

  #Generate group level x-values
  group_level_x <- mapply(function(groups, name) setNames(data.frame(groups, rnorm(length(groups))), paste0(c('', 'x_'), name)),
                          groups = groups_list, name = names(groups_list), SIMPLIFY = FALSE)

  #Generate individual level x-values
  ind_level_x <- rnorm(nrow(group_combinations))

  #Combine x-values
  x_data <- Reduce(f = merge, x = group_level_x, init = cbind(group_combinations, x_i = ind_level_x))

  #Sum up x-values across rows
  W <- rowSums(x_data[grepl(pattern = 'x_', names(x_data))])

  #Create group-level e-values
  group_level_e <- mapply(function(groups, name) setNames(data.frame(groups, rnorm(length(groups))), paste0(c('', 'e_'), name)),
                          groups = groups_list, name = names(groups_list), SIMPLIFY = FALSE)

  #Add heteroskedasticity if TRUE
  sd <- if(heterosked) 9*W^2 else 1

  #Generate individual level e-values (potentially with heteroskedasticity)
  ind_level_e <- rnorm(nrow(group_combinations), sd = sd)

  #Combine e-values
  e_data <- Reduce(f = merge, x = group_level_e, init = cbind(group_combinations, e_i = ind_level_e))

  #Sum up all e values and W
  Y <- W + rowSums(e_data[grepl(pattern = 'e_', names(e_data))])

  return(cbind(group_combinations, Y, W))
}


multiway_group_num <- function(groups_list, rho, theta){

  theta_groups <- ceiling(rho*length(groups_list[[1]]))
  one_groups <- length(groups_list[[1]]) - theta_groups

  dim_one_vec <- rep(c(theta, 1), times = c(theta_groups, one_groups))
  dim_one_groups <- list(rep(groups_list[[1]], times= dim_one_vec))
  names(dim_one_groups) <- names(groups_list)[1]

  all_groups <- c(dim_one_groups, groups_list[2:length(groups_list)])

}

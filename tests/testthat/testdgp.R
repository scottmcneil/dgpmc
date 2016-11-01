library(testthat)
library(dgpmc)

test_that('multiway_DGP works with two groups',{

  test_data <- data.frame(H = c('H1', 'H2', 'H1', 'H2'),
                          G = c('G1' ,'G1' ,'G2' ,'G2'),
                          Y = c(4.0727878, 0.6556725, 7.6870908, 2.8333871),
                          W = c(2.1383552, -0.3076943, 3.5153430, -0.0264946))


  num_dims <- 2
  groups <- c(2, 2)

  set.seed(42, 'Mersenne')
  function_data <- multiway_DGP(num_dims, groups)

  #Check if equal
  expect_equal(test_data, function_data, tolerance= 0.00001)

})

test_that('multiway_DGP works with unbalanced clusters',{

  test_data <- data.frame(H = c('H1', 'H1', 'H2', 'H1', 'H1', 'H2'),
                          G = c('G1', 'G1', 'G1', 'G2', 'G2', 'G2'),
                          Y = c(1.9210428, 2.1799217, 1.9234840, 0.2787875, 2.6078587, 3.3334203),
                          W = c(2.138355182, 1.627962342, 1.309952237, 1.909162014, 4.022244766, 0.005450335))


  num_dims <- 2
  groups <- c(2, 2)
  rho <- 0.5
  theta <- 2

  set.seed(42, 'Mersenne')
  function_data <- multiway_DGP(num_dims, groups, rho = rho, theta = theta)


  #Check if equal
  expect_equal(test_data, function_data, tolerance= 0.00001)

})

test_that('multiway_DGP works with heteroskedasticity',{

  test_data <- data.frame(H = c('H1', 'H2', 'H3', 'H1', 'H2', 'H3', 'H1', 'H2', 'H3'),
                          G = c('G1' ,'G1', 'G1' ,'G2' ,'G2', 'G2', 'G3', 'G3', 'G3'),
                          Y = c(-2.7178469, -2.85362704, 0.96973011, 4.9354339, 1.89093725, 0.56823349, -0.0831612, -1.90517057, -2.73148324),
                          W = c(3.515343, -0.02649464, 3.01441471, 1.7125126, 1.14443975, 3.05404209, -0.1240268, -0.94961147, 0.12368256))


  num_dims <- 2
  groups <- c(3, 3)
  heterosked <- TRUE

  set.seed(42, 'Mersenne')
  function_data <- multiway_DGP(num_dims, groups, heterosked = heterosked)
  function_data <- function_data[order(function_data['G'], function_data['H']),]
  row.names(function_data) <- 1:nrow(function_data)

  #Check if equal
  expect_equal(test_data, function_data, tolerance= 0.00001)

})

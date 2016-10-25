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

  test_data <- data.frame(H = c('H1', 'H2', 'H1', 'H2'),
                          G = c('G1' ,'G1' ,'G2' ,'G2'),
                          Y = c(-3.4479839, 0.6771162, 6.4144015, 2.2479845),
                          W = c(2.1383552, -0.3076943, 3.5153430, -0.0264946))


  num_dims <- 2
  groups <- c(2, 2)
  heterosked <- TRUE

  set.seed(42, 'Mersenne')
  function_data <- multiway_DGP(num_dims, groups, heterosked = heterosked)

  #Check if equal
  expect_equal(test_data, function_data, tolerance= 0.00001)

})

test_that('cgm_multiway_DGP works with two groups',{

  test_data <- data.frame(H = c('H1', 'H2', 'H1', 'H2'),
                          G = c('G1' ,'G1' ,'G2' ,'G2'),
                          Y = c(2.350344, -3.438974, 1.626835, 3.937426),
                          w_H = c(1.7752268, -0.6708227, 2.8824804, -0.6593572),
                          w_G = c(2.3815521, 0.3004143, 1.9377323, 2.9195080))


  num_dims <- 2
  groups <- c(2, 2)

  set.seed(42, 'Mersenne')
  function_data <- cgm_multiway_DGP(num_dims, groups)

  #Check if equal
  expect_equal(test_data, function_data, tolerance= 0.00001)

})

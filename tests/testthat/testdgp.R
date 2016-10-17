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

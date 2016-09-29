library(testthat)
library(dgpmc)

#Create dummy stat function
stat_func <- function(data, formula, x_interest, H0){

  model <- lm(data = data, formula = formula)

  x_coef <- coef(model)[[x_interest]]

  x_se <- sqrt(diag(vcov(model)))[[x_interest]]

  2*pt(abs(x_coef - H0)/x_se, df = df.residual(model), lower.tail = FALSE)

}

#Create dummy rand function
rand_func <- function(n, coef){

  x <- rnorm(n)
  y <- rnorm(n) + coef*x

  data.frame(y, x)

}

test_that('monte_carlo_stat works properly',{

  #Set out stat and rand arguments for single returned stat value
  stat_args <- list(formula = y ~ x, x_interest = 'x', H0 = 1)
  rand_args <- list(n = 10, coef = 1)

  #Create test p-value
  set.seed(42)
  test_p <- monte_carlo_stat(stat_func = stat_func, stat_args = stat_args, rand_func = rand_func, rand_args = rand_args)

  #Check if equal
  expect_equal(test_p, matrix(0.2857), tolerance = 0.0001)

  #Set out stat and rand arguments for multiple returned stat values
  stat_args <- list(formula = y ~ x, x_interest = 'x', H0 = c(1, 1.2))
  rand_args <- list(n = 10, coef = 1)

  #Create test p-value
  set.seed(42)
  test_p <- monte_carlo_stat(stat_func = stat_func, stat_args = stat_args, rand_func = rand_func, rand_args = rand_args)

  #Check if equal
  expect_equal(test_p, matrix(c(0.2857, 0.1833)), tolerance = 0.0001)

})

test_that('monte_carlo_block works properly',{

  #Set out arguments for single returned statistic
  stat_args <- list(formula = y ~ x, x_interest = 'x', H0 = 1)
  rand_args <- list(n = 10, coef = 1)

  cores <- 1
  block <- 1
  cumulative <- 3

  #Create progbar object
  pb <- txtProgressBar(min = 0, max = 5, style = 3)

  single_p <- matrix(0.2528999)

  #Set seed and generate
  set.seed(42, "L'Ecuyer")
  single_test_p <- monte_carlo_block(block = block, cumulative = cumulative, stat_func = stat_func, stat_args = stat_args,
                                     rand_func = rand_func, rand_args = rand_args, pb = pb, cores = cores)

  expect_equal(single_test_p, single_p, tolerance = 0.0001)
  expect_equal(getTxtProgressBar(pb), 3)

  #Try with two test statistics
  stat_args <- list(formula = y ~ x, x_interest = 'x', H0 = c(1, 1.2))

  #Set out two p-values
  double_p <- matrix(c(0.2528999, 0.5870307))

  #Set seed and generate values
  set.seed(42, "L'Ecuyer")
  double_test_p <- monte_carlo_block(block = block, cumulative = cumulative, stat_func = stat_func, stat_args = stat_args,
                                     rand_func = rand_func, rand_args = rand_args, pb = pb, cores = cores)

  #Test if equal
  expect_equal(double_test_p, double_p, tolerance = 0.0001)


  #Set out arguments for single returned statistic with multicore
  stat_args <- list(formula = y ~ x, x_interest = 'x', H0 = 1)
  rand_args <- list(n = 10, coef = 1)

  #Set block and core to 2
  cores <- 2
  block <- 2
  cumulative <- 4

  #Create progbar object
  pb <- txtProgressBar(min = 0, max = 5, style = 3)

  #Set seed and generate
  single_test_p <- monte_carlo_block(block = block, cumulative = cumulative, stat_func = stat_func, stat_args = stat_args,
                                     rand_func = rand_func, rand_args = rand_args, pb = pb, cores = cores)

  expect_equal(length(single_test_p), 2)
  expect_equal(getTxtProgressBar(pb), 4)


  #Set out arguments for double returned statistic with multicore
  stat_args <- list(formula = y ~ x, x_interest = 'x', H0 = c(1, 1.2))

  #Set block and core to 2
  cores <- 2
  block <- 2
  cumulative <- 4

  #Create progbar object
  pb <- txtProgressBar(min = 0, max = 5, style = 3)

  #Set seed and generate
  single_test_p <- monte_carlo_block(block = block, cumulative = cumulative, stat_func = stat_func, stat_args = stat_args,
                                     rand_func = rand_func, rand_args = rand_args, pb = pb, cores = cores)

  expect_equal(nrow(single_test_p), 2)
  expect_equal(ncol(single_test_p), 2)
  expect_equal(getTxtProgressBar(pb), 4)

})

test_that('monte_carlo_recur works properly',{

  #Set out arguments for one p-value
  stat_args <- list(formula = y ~ x, x_interest = 'x', H0 = 1)
  rand_args <- list(n = 10, coef = 1)
  reps <- 3
  cores <- 1
  progbar <- TRUE

  #Set out expected p-values
  single_p <- matrix(c(0.2528999, 0.7366565, 0.3029334), nrow = 1, ncol = 3)

  #Set seed and generate
  set.seed(42, "L'Ecuyer")
  single_test_p <- monte_carlo_recur(reps, stat_func = stat_func, stat_args = stat_args, rand_func = rand_func, rand_args = rand_args,
                                     progbar = progbar, cores = cores)

  #Test if equal
  expect_equal(single_test_p, single_p, tolerance = 0.001)

  #Set out stat args for two values returned
  stat_args <- list(formula = y ~ x, x_interest = 'x', H0 = c(1, 1.2))

  #Set out expected returned values
  double_p <- matrix(c(0.2528999, 0.7366565, 0.3029334, 0.5870307, 0.5618727, 0.1071590), nrow = 2, ncol = 3, byrow = TRUE)

  #Set seed and generate values
  set.seed(42, "L'Ecuyer")
  double_test_p <- monte_carlo_recur(reps, stat_func = stat_func, stat_args = stat_args, rand_func = rand_func, rand_args = rand_args,
                                     progbar = progbar, cores = cores)

  #Check if equal
  expect_equal(double_test_p, double_p, tolerance = 0.001)

  #Set out args for 1 value and 2 cores
  stat_args <- list(formula = y ~ x, x_interest = 'x', H0 = 1)
  cores <- 2

  #Generate data
  single_test_p <- monte_carlo_recur(reps, stat_func = stat_func, stat_args = stat_args, rand_func = rand_func, rand_args = rand_args,
                                     progbar = progbar, cores = cores)

  #Check if dimensions correct
  expect_equal(nrow(single_test_p), 1)
  expect_equal(ncol(single_test_p), 3)


  #Set out args for 2 values and 2 cores
  stat_args <- list(formula = y ~ x, x_interest = 'x', H0 = c(1, 1.2))
  cores <- 2

  #Generate data
  single_test_p <- monte_carlo_recur(reps, stat_func = stat_func, stat_args = stat_args, rand_func = rand_func, rand_args = rand_args,
                                     progbar = progbar, cores = cores)

  #Check if dimensions correct
  expect_equal(nrow(single_test_p), 2)
  expect_equal(ncol(single_test_p), 3)

})


test_that('dgpmc works properly',{

  #Set out arguments
  stat_args <- list(formula = y ~ x, x_interest = 'x', H0 = 1)
  rand_args <- list(n = 10, coef = 1)
  reps <- 3
  cores <- 1

  single_p <- data.frame(First = c(0.2528999, 0.7366565, 0.3029334))

  set.seed(42)
  single_test_p <- dgpmc(reps, stat_func = stat_func, stat_args = stat_args, rand_func = rand_func, rand_args = rand_args, names = 'First')

  expect_equal(single_test_p$stat, single_p, tolerance = 0.001)
  expect_equal(single_test_p$stat_func, stat_func)
  expect_equal(single_test_p$stat_args, stat_args)
  expect_equal(single_test_p$rand_func, rand_func)
  expect_equal(single_test_p$rand_args, rand_args)

  cores <- 2
  single_test_p <- dgpmc(reps, stat_func = stat_func, stat_args = stat_args, rand_func = rand_func, rand_args = rand_args, names = 'First',
                         cores = cores, progbar = TRUE)

  expect_equal(nrow(single_test_p$stat), 3)
  expect_equal(ncol(single_test_p$stat), 1)

  stat_args <- list(formula = y ~ x, x_interest = 'x', H0 = c(1, 1.2))
  rand_args <- list(n = 10, coef = 1)
  reps <- 3

  double_p <- data.frame(First = c(0.2528999, 0.7366565, 0.3029334), Second = c(0.5870307, 0.5618727, 0.1071590))

  set.seed(42)
  double_test_p <- dgpmc(reps, stat_func = stat_func, stat_args = stat_args, rand_func = rand_func, rand_args = rand_args, names = c('First','Second'))

  expect_equal(double_test_p$stat, double_p, tolerance = 0.001)
  expect_equal(double_test_p$stat_func, stat_func)
  expect_equal(double_test_p$stat_args, stat_args)
  expect_equal(double_test_p$rand_func, rand_func)
  expect_equal(double_test_p$rand_args, rand_args)

})

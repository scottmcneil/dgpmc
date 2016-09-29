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
  expect_equal(test_p, 0.2857, tolerance = 0.0001)

  #Set out stat and rand arguments for multiple returned stat values
  stat_args <- list(formula = y ~ x, x_interest = 'x', H0 = c(1, 1.2))
  rand_args <- list(n = 10, coef = 1)

  #Create test p-value
  set.seed(42)
  test_p <- monte_carlo_stat(stat_func = stat_func, stat_args = stat_args, rand_func = rand_func, rand_args = rand_args)

  #Check if equal
  expect_equal(test_p, c(0.2857, 0.1833), tolerance = 0.0001)

})

test_that('monte_carlo_block works properly',{

  #Set out arguments for single returned statistic
  stat_args <- list(formula = y ~ x, x_interest = 'x', H0 = 1)
  rand_args <- list(n = 10, coef = 1, seed = 42)

  cores <- 1
  block <- 1
  cumulative <- 3

  #Create progbar object
  pb <- txtProgressBar(min = 0, max = 5, style = 3)

  single_p <- c(0.2857, 0.8692, 0.5286)

  #Set seed and generate
  set.seed(42, "L'Ecuyer")
  single_test_p <- monte_carlo_block(block = block, cumulative = cumulative, stat_func = stat_func, stat_args = stat_args,
                                     rand_func = rand_func, rand_args = rand_args, pb = pb, cores = cores)

  expect_equal(single_test_p, single_p, tolerance = 0.001)
  expect_equal(getTxtProgressBar(progbar), 3)

#
#   stat_args <- list(formula = y ~ x, x_interest = 'x', H0 = c(1, 1.2))
#   rand_args <- list(n = 10, coef = 1)
#   reps <- 3
#
#   progbar <- txtProgressBar(min = 0, max = reps, style = 3)
#
#   double_p <- matrix(c(0.2857, 0.1833, 0.8692, 0.4705, 0.5286, 0.8119), nrow = 2, ncol = 3)
#
#   set.seed(42)
#   double_test_p <- monte_carlo_recur(reps, stat_func = stat_func, stat_args = stat_args, rand_func = rand_func, rand_args = rand_args, progbar = progbar)
#
#   expect_equal(double_test_p, double_p, tolerance = 0.001)
#
#   stat_args <- list(formula = y ~ x, x_interest = 'x', H0 = c(1, 1.2))
#   rand_args <- list(n = 10, coef = 1)
#   reps <- 3
#
#   progbar <- txtProgressBar(min = 0, max = reps, style = 3)
#
#   double_p <- matrix(c(0.2857, 0.1833, 0.8692, 0.4705, 0.5286, 0.8119), nrow = 2, ncol = 3)
#
#   set.seed(42)
#   double_test_p <- monte_carlo_recur(reps, stat_func = stat_func, stat_args = stat_args, rand_func = rand_func, rand_args = rand_args,
#                                      progbar = progbar, cores = 2)
#
#   expect_equal(double_test_p, double_p, tolerance = 0.001)

})

# test_that('monte_carlo_recur works properly',{
#
#   stat_args <- list(formula = y ~ x, x_interest = 'x', H0 = 1)
#   rand_args <- list(n = 10, coef = 1)
#   reps <- 3
#
#   progbar <- txtProgressBar(min = 0, max = reps, style = 3)
#
#   single_p <- c(0.2857, 0.8692, 0.5286)
#
#   set.seed(42)
#   single_test_p <- monte_carlo_recur(reps, stat_func = stat_func, stat_args = stat_args, rand_func = rand_func, rand_args = rand_args, progbar = progbar)
#
#   expect_equal(single_test_p, single_p, tolerance = 0.001)
#   expect_equal(getTxtProgressBar(progbar), 3)
#
#
#   stat_args <- list(formula = y ~ x, x_interest = 'x', H0 = c(1, 1.2))
#   rand_args <- list(n = 10, coef = 1)
#   reps <- 3
#
#   progbar <- txtProgressBar(min = 0, max = reps, style = 3)
#
#   double_p <- matrix(c(0.2857, 0.1833, 0.8692, 0.4705, 0.5286, 0.8119), nrow = 2, ncol = 3)
#
#   set.seed(42)
#   double_test_p <- monte_carlo_recur(reps, stat_func = stat_func, stat_args = stat_args, rand_func = rand_func, rand_args = rand_args, progbar = progbar)
#
#   expect_equal(double_test_p, double_p, tolerance = 0.001)
#
#   stat_args <- list(formula = y ~ x, x_interest = 'x', H0 = c(1, 1.2))
#   rand_args <- list(n = 10, coef = 1)
#   reps <- 3
#
#   progbar <- txtProgressBar(min = 0, max = reps, style = 3)
#
#   double_p <- matrix(c(0.2857, 0.1833, 0.8692, 0.4705, 0.5286, 0.8119), nrow = 2, ncol = 3)
#
#   set.seed(42)
#   double_test_p <- monte_carlo_recur(reps, stat_func = stat_func, stat_args = stat_args, rand_func = rand_func, rand_args = rand_args,
#                                      progbar = progbar, cores = 2)
#
#   expect_equal(double_test_p, double_p, tolerance = 0.001)
#
# })
#
# test_that('dgpmc works properly',{
#
#   stat_func <- function(data, formula, x_interest, H0){
#
#     model <- lm(data = data, formula = formula)
#
#     x_coef <- coef(model)[[x_interest]]
#
#     x_se <- sqrt(diag(vcov(model)))[[x_interest]]
#
#     2*pt(abs(x_coef - H0)/x_se, df = df.residual(model), lower.tail = FALSE)
#
#   }
#
#   rand_func <- function(n, coef){
#     x <- rnorm(n)
#     y <- rnorm(n) + coef*x
#     data.frame(y, x)
#   }
#
#   stat_args <- list(formula = y ~ x, x_interest = 'x', H0 = 1)
#   rand_args <- list(n = 10, coef = 1)
#   reps <- 3
#
#   single_p <- data.frame(First = c(0.2857, 0.8692, 0.5286))
#
#   set.seed(42)
#   single_test_p <- dgpmc(reps, stat_func = stat_func, stat_args = stat_args, rand_func = rand_func, rand_args = rand_args, names = 'First')
#
#   expect_equal(single_test_p$stat, single_p, tolerance = 0.001)
#   expect_equal(single_test_p$stat_func, stat_func)
#   expect_equal(single_test_p$stat_args, stat_args)
#   expect_equal(single_test_p$rand_func, rand_func)
#   expect_equal(single_test_p$rand_args, rand_args)
#
#   stat_args <- list(formula = y ~ x, x_interest = 'x', H0 = c(1, 1.2))
#   rand_args <- list(n = 10, coef = 1)
#   reps <- 3
#
#   double_p <- data.frame(First = c(0.2857, 0.8692, 0.5286), Second = c(0.1833, 0.4705, 0.8119))
#
#   set.seed(42)
#   double_test_p <- dgpmc(reps, stat_func = stat_func, stat_args = stat_args, rand_func = rand_func, rand_args = rand_args, names = c('First','Second'))
#
#   expect_equal(double_test_p$stat, double_p, tolerance = 0.001)
#   expect_equal(double_test_p$stat_func, stat_func)
#   expect_equal(double_test_p$stat_args, stat_args)
#   expect_equal(double_test_p$rand_func, rand_func)
#   expect_equal(double_test_p$rand_args, rand_args)
#
# })

library(testthat)
library(dgpmc)

test_that('monte_carlo_stat works properly',{

  stat_func <- function(data, formula, x_interest, H0){

    model <- lm(data = data, formula = formula)

    x_coef <- coef(model)[[x_interest]]

    x_se <- sqrt(diag(vcov(model)))[[x_interest]]

    2*pt(abs(x_coef - H0)/x_se, df = df.residual(model), lower.tail = FALSE)

  }

  rand_func <- function(n, coef, seed){
    set.seed(seed)
    x <- rnorm(n)
    y <- rnorm(n) + coef*x
    data.frame(y, x)
  }

  stat_args <- list(formula = y ~ x, x_interest = 'x', H0 = c(1, 1.2))
  rand_args <- list(n = 10, coef = 1, seed = 42)
  reps <- 10
  rep <- 5

  progbar <- txtProgressBar(min = 0, max = reps, style = 3)

  test_p <- monte_carlo_stat(rep, stat_func = stat_func, stat_args = stat_args, rand_func = rand_func, rand_args = rand_args, progbar = progbar)

  expect_equal(test_p, c(0.2857, 0.1833), tolerance = 0.0001)
  expect_equal(getTxtProgressBar(progbar), 5)

  stat_args <- list(formula = y ~ x, x_interest = 'x', H0 = 1)
  rand_args <- list(n = 10, coef = 1, seed = 42)

  progbar <- txtProgressBar(min = 0, max = reps, style = 3)

  test_p <- monte_carlo_stat(rep, stat_func = stat_func, stat_args = stat_args, rand_func = rand_func, rand_args = rand_args, progbar = progbar)
  expect_equal(test_p, 0.2857, tolerance = 0.0001)

})

test_that('monte_carlo_recur works properly',{

  stat_func <- function(data, formula, x_interest, H0){

    model <- lm(data = data, formula = formula)

    x_coef <- coef(model)[[x_interest]]

    x_se <- sqrt(diag(vcov(model)))[[x_interest]]

    2*pt(abs(x_coef - H0)/x_se, df = df.residual(model), lower.tail = FALSE)

  }

  rand_func <- function(n, coef){
    x <- rnorm(n)
    y <- rnorm(n) + coef*x
    data.frame(y, x)
  }

  stat_args <- list(formula = y ~ x, x_interest = 'x', H0 = 1)
  rand_args <- list(n = 10, coef = 1)
  reps <- 3

  progbar <- txtProgressBar(min = 0, max = reps, style = 3)

  single_p <- c(0.2857, 0.8692, 0.5286)

  set.seed(42)
  single_test_p <- monte_carlo_recur(reps, stat_func = stat_func, stat_args = stat_args, rand_func = rand_func, rand_args = rand_args, progbar = progbar)

  expect_equal(single_test_p, single_p, tolerance = 0.001)
  expect_equal(getTxtProgressBar(progbar), 3)


  stat_args <- list(formula = y ~ x, x_interest = 'x', H0 = c(1, 1.2))
  rand_args <- list(n = 10, coef = 1)
  reps <- 3

  progbar <- txtProgressBar(min = 0, max = reps, style = 3)

  double_p <- matrix(c(0.2857, 0.1833, 0.8692, 0.4705, 0.5286, 0.8119), nrow = 2, ncol = 3)

  set.seed(42)
  double_test_p <- monte_carlo_recur(reps, stat_func = stat_func, stat_args = stat_args, rand_func = rand_func, rand_args = rand_args, progbar = progbar)

  expect_equal(double_test_p, double_p, tolerance = 0.001)

})

test_that('dgpmc works properly',{

  stat_func <- function(data, formula, x_interest, H0){

    model <- lm(data = data, formula = formula)

    x_coef <- coef(model)[[x_interest]]

    x_se <- sqrt(diag(vcov(model)))[[x_interest]]

    2*pt(abs(x_coef - H0)/x_se, df = df.residual(model), lower.tail = FALSE)

  }

  rand_func <- function(n, coef){
    x <- rnorm(n)
    y <- rnorm(n) + coef*x
    data.frame(y, x)
  }

  stat_args <- list(formula = y ~ x, x_interest = 'x', H0 = 1)
  rand_args <- list(n = 10, coef = 1)
  reps <- 3

  single_p <- data.frame(First = c(0.2857, 0.8692, 0.5286))

  set.seed(42)
  single_test_p <- dgpmc(reps, stat_func = stat_func, stat_args = stat_args, rand_func = rand_func, rand_args = rand_args, names = 'First')

  expect_equal(single_test_p$stat, single_p, tolerance = 0.001)
  expect_equal(single_test_p$stat_func, stat_func)
  expect_equal(single_test_p$stat_args, stat_args)
  expect_equal(single_test_p$rand_func, rand_func)
  expect_equal(single_test_p$rand_args, rand_args)

  stat_args <- list(formula = y ~ x, x_interest = 'x', H0 = c(1, 1.2))
  rand_args <- list(n = 10, coef = 1)
  reps <- 3

  double_p <- data.frame(First = c(0.2857, 0.8692, 0.5286), Second = c(0.1833, 0.4705, 0.8119))

  set.seed(42)
  double_test_p <- dgpmc(reps, stat_func = stat_func, stat_args = stat_args, rand_func = rand_func, rand_args = rand_args, names = c('First','Second'))

  expect_equal(double_test_p$stat, double_p, tolerance = 0.001)
  expect_equal(double_test_p$stat_func, stat_func)
  expect_equal(double_test_p$stat_args, stat_args)
  expect_equal(double_test_p$rand_func, rand_func)
  expect_equal(double_test_p$rand_args, rand_args)

})

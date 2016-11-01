
<!-- README.md is generated from README.Rmd. Please edit that file -->
dgpmc
=====

dgpmc (Data Generation Process Monte Carlo) is a package designed as a thin wrapper for performing Monte Carlo simulations on statistical functions and random data.

Installation
------------

The easiest way to install the package is via devtools. First, run the following code in the RStudio console:

    install.packages('devtools')

Once devtools is installed, go get a personal API key [here](https://github.com/settings/tokens). The page will look like the screenshots below. Select "Generate new token" then you may need to input your password.

<kbd>![](screenshots/API%20Key%201.png)</kbd>

When you're on the screen below, you'll need to give your token a description (it does not need to be "dgpmc"). Also make sure the "repo" option is selected, like below.

<kbd>![](screenshots/API%20Key%202.png)</kbd>

Finally, make sure you copy the new API token before closing the page, or else you'll need to re-generate it.

Once you've got your API token, run the code below in the RStudio console, replacing YOUR\_AUTH\_TOKEN with the token you just generated.

    devtools::install_github(repo = 'scottmcneil/dgpmc', auth_token = 'YOUR_AUTH_TOKEN')

Usage
-----

This package contains just one function for running Monte Carlo simulations, `dgpmc`. The function has five required arguments and three optional arguments.

The five required arguments are:

<table style="width:72%;">
<colgroup>
<col width="12%" />
<col width="59%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">Argument</th>
<th align="left">Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left"><code>reps</code></td>
<td align="left">The number of total Monte Carlo replications</td>
</tr>
<tr class="even">
<td align="left"><code>stat_func</code></td>
<td align="left">A function that will take data and return a vector of statistics</td>
</tr>
<tr class="odd">
<td align="left"><code>stat_args</code></td>
<td align="left">A list of additional arguments to pass to <code>stat_func</code></td>
</tr>
<tr class="even">
<td align="left"><code>rand_func</code></td>
<td align="left">A function that will generate random data to pass to <code>stat_func</code></td>
</tr>
<tr class="odd">
<td align="left"><code>rand_args</code></td>
<td align="left">A list of arguments to pass to <code>rand_func</code></td>
</tr>
</tbody>
</table>

The two optional arguments are:

| Argument  | Description                                                          |
|-----------|----------------------------------------------------------------------|
| `names`   | Names to identify each element of the vector returned by `stat_func` |
| `progbar` | A boolean indicating whether or not to display a progress bar        |
| `cores`   | An integer indicating the number of cores you'd like to use          |

For each replication, `dgpmc` will run `rand_func` with the elements of `rand_args` as parameters. It will then pass the resulting data and the elements of `stat_args` to `stat_func`. The results of all the replications will be combined into an R datframe with `rep` rows and a column for each statistic and columns named using `names`.

The function returns a list, where the `stat` element contains the dataframe of results. The other elements in the list are the supplied `stat_func`, `stat_args`, `rand_func` and `rand_args`.

Example
-------

Below is an example of usage with very simple examples for `stat_func` and `rand_func`.


    stat_func <- function(data, formula, x_interest, H0){

      #Create linear model based on data and formula
      model <- lm(data = data, formula = formula)

      #Extract coefficient of x of interest
      x_coef <- coef(model)[[x_interest]]

      #Extract standard error of x of interest
      x_se <- sqrt(diag(vcov(model)))[[x_interest]]

      #Return vectorized hypothesis test on vector of H0
      2*pt(abs(x_coef - H0)/x_se, df = df.residual(model), lower.tail = FALSE)

    }

    rand_func <- function(n, coef){

      #Create x vector of random normal variables
      x <- rnorm(n)
      
      #Create y vector correlated to x vecto
      y <- rnorm(n) + coef*x
      
      #Return dataframe of combined observations
      data.frame(y, x)

    }

    stat_args <- list(formula = y ~ x, x_interest = 'x', H0 = c(0.8, 1, 1.2))
    rand_args <- list(n = 10, coef = 1)
    reps <- 8
    stat_names <- c('H01', 'H02', 'H03')

    mc_results <- dgpmc(reps = reps, stat_func = stat_func, stat_args = stat_args, rand_func = rand_func, rand_args = rand_args, names = statnames)

    print(mc_results$stat)

Will produce the following table:

|        H01|        H02|        H03|
|----------:|----------:|----------:|
|  0.4299363|  0.2857353|  0.1833407|
|  0.6873946|  0.8691702|  0.4704641|
|  0.3146273|  0.5281681|  0.8119407|
|  0.3351245|  0.6254144|  0.9920525|
|  0.0468053|  0.0988551|  0.2032164|
|  0.5201904|  0.8142137|  0.8564774|
|  0.9401863|  0.6682665|  0.3618710|
|  0.4594450|  0.2550542|  0.1324828|

Built-in DGP functions
----------------------

This library currently has two built-in DGP functions:

-   `simple_dgp` creates a random normal data set, clusterd by group
-   `autocorr_dgp` creates auto-correlated random data, clustered by both time and group
-   `multiway_DGP` creates data correlated by an arbitrary number of group dimensions

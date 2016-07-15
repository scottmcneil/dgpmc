
<!-- README.md is generated from README.Rmd. Please edit that file -->
wildclusterbootsim
==================

wildclusterboot is a package designed to implement Monte Carlo simulations of clustered standard errors and the wild cluster bootstrap and save the results.

Installation
------------

The easiest way to install this package is via devtools. First, run the following code:

    install.packages('devtools')

Once devtools is installed, go get a [personal API key](https://github.com/settings/tokens) and then run the following code:

    devtools::install_github(repo = 'scottmcneil/wildclusterbootsim', auth_token = [auth token you just generated])

Usage
-----

The main two functions in this library are `wild_cluster_boot_mc` and `muli_wild_cluster_boot_mc`.

`wild_cluster_boot_mc` allows you to specify a DGP function, a set of DGP arguments, a set of bootstrap arguments and where to save the results. It will then run a set number of Monte Carlo replications and save the results along with the corresponding DGP parameters.

For example:

    #Specify DGP function
    dgp <- test_dgp

    #Specify DGP parameters
    t <- 5
    G <- 5
    ng <- 30

    #Create list of dgp_args
    dgp_args <- list(t = t, G = G, ng = ng, rho = rho, lambda = lambda, gamma = gamma)

    #Specify bootstrap parameters
    formula <- 'Y ~ X'
    reps <- 1000
    x_interest <- 'X'
    clusterby <- c('G', 'tG')
    boot_dist <- 'six_pt'
    boot_reps <- 399
    bootby <- c('G', 'tG', 'i')
    H0 <- 1
    cores <- 3
    progbar = TRUE
    save_type = 'db'
    save_file <- 'wildclusterbootsim.db'

    #Run simulation
    wild_cluster_boot_mc(dgp = dgp, dgp_args = dgp_args, formula = formula, reps = reps,
                              x_interest = x_interest, clusterby = clusterby, boot_dist = boot_dist, boot_reps = boot_reps,
                              bootby = bootby,  H0 = H0, cores = cores, progbar = progbar,
                              save_type = save_type, save_file = save_file)

Currently the only save type that's supported is `'db'` but I can add others.

`multi_wild_cluster_boot_mc` is essentially the same as `wild_cluster_boot_mc`, but it allows you to specify vectors for DGP arguments, and will run a Monte Carlo on each possible combination of those arguments. It also allows you to turn on a progress bar.

For example:

    #Specify several group and time counts
    t <- c(5, 6, 7)
    G <- (5, 10, 15)
    ng <- 30

    #Create your list of DGP arguments
    dgp_args <- list(t = t, G = G, ng = ng, rho = rho, lambda = lambda, gamma = gamma)

    #Then run the simulation as before
    muli_wild_cluster_boot_mc(dgp = dgp, dgp_args = dgp_args, formula = formula, reps = reps,
                              x_interest = x_interest, clusterby = clusterby, boot_dist = boot_dist, boot_reps = boot_reps,
                              bootby = bootby,  H0 = H0, cores = cores, progbar = progbar,
                              save_type = save_type, save_file = save_file)

The above code would run nine total simulations, one for each combination of `t` and `G`.

Built-in DGP functions
----------------------

This library currently has two built-in DGP functions:

-   `simple_dgp` creates a random normal data clusterd by group
-   `autocorr_dgp` creates auto-correlated random data, clustered by both time and group

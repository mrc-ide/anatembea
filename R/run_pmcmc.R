#' Run a pMCMC
#'
#' This function sets up and runs a particle MCMC that uses Dust, Odin and MCState
#'
#' @param data_raw Time series data to fit model
#' @param n_particles Number of particles to be used in pMCMC (default = 200)
#' @param init_EIR A single value or a dataframe with two columns (time and EIR)
#'                  to specify historical malaria transmission levels before
#'                  data collection began.
#' @param target_prev Return an initial EIR value (from the equilibrium solution),
#'                    given a target prevalence in under 5yos
#' @param proposal_matrix Proposal matrix for MCMC parameters
#' @param max_param Ceiling for proposed stochastic parameter (either EIR or betaa) values (default = 1000)
#' @param prop_treated Proportion of clinical cases that receive effective treatment (default = 40%)
#' @param n_steps Number of MCMC steps in a single chain (default = 500)
#' @param n_threads Number of processing threads (default = 4)
#' @param n_chains Number of chains (default = 1)
#' @param n_workers Number of workers (default = 4)
#' @param state_check If state_check = TRUE, returns expected deriv values which should equal 0 and sets stochastic model to have EIR constant at init_EIR
#'                    If state_check = TRUE and seasonality_on = 1, then the deterministic seasonal model is still run, but theta2 is forced to 1, forcing a constant seasonality profile
#'                    If state_check = FALSE, no values are printed
#' @param country Name of country (needed if using seasonality model)
#' @param admin_unit Name of administrative unit (needed if using seasonality model)
#' @param preyears Length of time in years the deterministic seasonal model should run before Jan 1 of the year observations began (default = 2)
#' @param seasonality_on Toggle seasonality model run before observed time period (default = 1)
#' @param seasonality_check Toggle saving values of seasonality equilibrium (default = 1)
#' @param check_flexibility Toggle saving values of flexibility period
#' @param seed Allows user to specify a seed (default = 1L)
#' @param start_pf_time Number of days before first observation that particle filter will start (default = 30)
#' @param particle_tune Logical to determine if tuning the number of particles should be performed.
#' @param initial Is the initial equilibrium state informed by the user ('user-informed') or by the observed data ('fitted')?
#' @param comparison The comparison function to be used. Either 'u5' which
#'          equates the observed prevalence to prevalence under 5 years old in
#'          the model or 'pgmg' which calculates prevalence in primigravid and
#'          multigravid pregnant women for comparison with observed ANC data.
#' @export
run_pmcmc <- function(data_raw=NULL,
                      data_raw_pg=NULL,
                      data_raw_mg=NULL,
                      init_EIR = 10,
                      target_prev=NULL,
                      target_prev_group='u5',
                      n_particles=200,
                      proposal_matrix,
                      max_param=1000,
                      prop_treated = 0.4,
                      n_steps = 500,
                      n_threads = 4,
                      n_chains = 1,
                      n_workers = 1,
                      state_check = 0,## Run equilibrium checks
                      # If state_check = TRUE, returns expected deriv values which should equal 0 and sets stochastic model to have EIR constant at init_EIR
                      # If state_check = TRUE and seasonality_on = 1, then the deterministic seasonal model is still run, but theta2 is forced to 1, forcing a constant seasonality profile
                      # If state_check = FALSE, no values are printed
                      country = NULL,
                      admin_unit = NULL,
                      seasonality_on = TRUE,  ## state_check = TRUE runs a deterministic seasonal model before running the stochastic model to get more realistic immunity levels
                      ## If seasonality_on = FALSE, runs the stochastic model based on the standard equilibrium solution
                      seasonality_check = FALSE,##If TRUE, saves values of seasonality equilibrium
                      check_flexibility = FALSE,
                      seed = 1L,
                      start_pf_time = 30,
                      particle_tune = FALSE,
                      comparison = c('u5','pg','sg','mg','pgmg','pgsg','ancall'),
                      initial = 'informed'){
  ##Merge primigrav and multigrav datasets if necessary.
  if(comparison=='pgmg' | comparison=='pgsg'){
    data_raw_pg <- format_na(data_raw_pg)
    data_raw_mg <- format_na(data_raw_mg)

    data_raw <- dplyr::left_join(data_raw_pg,data_raw_mg,by=c('month'),suffix = c('.pg','.mg'))
    pg_avg_prev <- sum(data_raw[1:12,'positive.pg'],na.rm=TRUE)/sum(data_raw[1:12,'tested.pg'],na.rm=TRUE)
    mg_avg_prev <- sum(data_raw[1:12,'positive.mg'],na.rm=TRUE)/sum(data_raw[1:12,'tested.mg'],na.rm=TRUE)
    avg_prev <- c(pg_avg_prev,mg_avg_prev)
  } else {
    data_raw <- format_na(data_raw)
    avg_prev <- sum(data_raw[1:12,'positive'],na.rm=TRUE)/sum(data_raw[1:12,'tested'],na.rm=TRUE)
  }
  Sys.setenv("MC_CORES"=n_threads)
  # ## Modify dates from data
  data_proc <- data_process(data_raw=data_raw,start_pf_time=start_pf_time,check_flexibility = check_flexibility)
  data <- data_proc$data
  stochastic_schedule <- data_proc$stochastic_schedule

  #Provide age categories, and number of heterogeneity brackets
  init_age <- c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3.5, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80)
  het_brackets <- 5
  lag_rates <- 10
  max_steps <- 1e7
  atol <- 1e-6 # Probably alright to bring up a bit, maybe to 1e-4 (same with rtol)
  rtol <- 1e-6
  preyears <- 5 #Length of time in years the deterministic seasonal model should run before Jan 1 of the year observations began


  #Create model parameter list. Also loads seasonality profile data file to match to desired admin_unit and country
  mpl_pf <- mamasante::model_param_list_create(init_age = init_age,
                                    prop_treated = prop_treated,
                                    het_brackets = het_brackets,
                                    max_param = max_param,
                                    target_prev = target_prev,
                                    target_prev_group = target_prev_group,
                                    state_check = state_check,
                                    lag_rates = lag_rates,
                                    country = country,
                                    admin_unit = admin_unit,
                                    start_stoch = data_proc$start_stoch,
                                    time_origin = data_proc$time_origin,
                                    seasonality_on = seasonality_on,
                                    preyears = preyears,
                                    particle_tune = particle_tune,
                                    comparison = comparison,
                                    avg_prev = avg_prev)

  ## If a deterministic seasonal model is needed prior to the stochastic model, this loads the deterministic odin model
  det_model <- NULL
  if(seasonality_on & !is.data.frame(init_EIR)){
    odin_det <- system.file("odin", "odin_model_stripped_seasonal.R", package = "mamasante")
    det_model <- suppressMessages(odin::odin(odin_det))
  } else if(!(seasonality_on) & is.data.frame(init_EIR)){
    odin_det <- system.file("odin", "odin_model_stripped_matched.R", package = "mamasante")
    det_model <- suppressMessages(odin::odin(odin_det))
  } else if(seasonality_on & is.data.frame(init_EIR)){
    print('Seasonality not supported with multiple EIR values.')
    print('Reverting to piece-wise constant EIR.')
    odin_det <- system.file("odin", "odin_model_stripped_matched.R", package = "mamasante")
    det_model <- suppressMessages(odin::odin(odin_det))
  }

  ## Load stochastic model in odin.dust
  ##Switch between EIR and mosquito emergence models
  stoch_file <- 'odinmodelmatchedstoch_mozemerg.R'
  odin_stoch <- system.file("odin", stoch_file, package = "mamasante")

  model <- suppressMessages(odin.dust::odin_dust(odin_stoch))

  if(!model$public_methods$has_openmp()) warning('openmp must be enabled to run particle filter in parallel')

  volatility <- mcstate::pmcmc_parameter("volatility", rgamma(1,shape = 3.4, rate = 3.1), min = 0,
                                         prior = function(p) dgamma(p, shape = 3.4, rate = 3.1, log = TRUE))
  if(initial == 'informed'){
    ## Set initial state based on a user-given equilibrium EIR
    init_state <- initialise(init_EIR=init_EIR,mpl=mpl_pf,det_model=det_model)
    ### Set pmcmc parameters
    # init_betaa <- mcstate::pmcmc_parameter("init_betaa", rgamma(1,shape = 0.64, rate = 0.057), min = 0,
    #                                        prior = function(p) dgamma(p, shape = 0.64, rate = 0.057, log = TRUE))

    # pars = list(init_betaa = init_betaa, volatility = volatility) ## Put pmcmc parameters into a list
    pars = list(volatility = volatility) ## Put pmcmc parameters into a list


    mcmc_pars <- mcstate::pmcmc_parameters$new(pars,
                                               proposal_matrix,
                                               transform = user_informed(init_state)) ## Calls transformation function based on pmcmc parameters
  }else if(initial == 'fitted'){
    log_init_EIR <- mcstate::pmcmc_parameter("log_init_EIR", rnorm(1, mean = 4, sd = 3),
                                             prior = function(p) dnorm(p, mean = 4, sd = 3, log = TRUE) + p) #Add p to adjust for sampling on log scale
    pars = list(log_init_EIR = log_init_EIR, volatility = volatility) ## Put pmcmc parameters into a list
    mcmc_pars <- mcstate::pmcmc_parameters$new(pars,
                                               proposal_matrix,
                                               transform = data_informed(mpl_pf,det_model)) ## Calls transformation function based on pmcmc parameters
  }


  n_threads <- dust::dust_openmp_threads(n_threads, action = "fix")
  if(comparison=='u5'){
    ##Output from particle filter
    ##    run: output used for likelihood calculation
    ##    state: output used for visualization
    index <- function(info) {
      list(run = c(prev = info$index$prev),
           state = c(prev_05 = info$index$prev,
                     EIR = info$index$EIR_out,
                     betaa = info$index$betaa_out,
                     clininc_all = info$index$inc,
                     prev_all = info$index$prevall,
                     clininc_05 = info$index$inc05,
                     Dout = info$index$Dout,
                     Aout = info$index$Aout,
                     Uout = info$index$Uout,
                     p_det_out = info$index$p_det_out,
                     phi_out = info$index$phi_out,
                     b_out = info$index$b_out,
                     IC_out = info$index$IC_out,
                     IB_out = info$index$IB_out,
                     ID_out = info$index$ID_out,
                     spz_rate = info$index$spz_rate,
                     eff_moz_pop = info$index$eff_moz_pop,
                     moz2human_ratio = info$index$moz2human_ratio))
    }
    compare_funct <- compare_u5
 } else if(comparison=='pg'){
    ##Output from particle filter
    ##    run: output used for likelihood calculation
    ##    state: output used for visualization
    index <- function(info) {
      list(run = c(prev_pg = info$index$prev_preg_pg),
           state = c(prev_05 = info$index$prev,
                     prev_pg = info$index$prev_preg_pg,
                     prev_mg = info$index$prev_preg_mg,
                     EIR = info$index$EIR_out,
                     betaa = info$index$betaa_out,
                     clininc_all = info$index$inc,
                     prev_all = info$index$prevall,
                     clininc_05 = info$index$inc05,
                     Dout = info$index$Dout,
                     Aout = info$index$Aout,
                     Uout = info$index$Uout,
                     p_det_out = info$index$p_det_out,
                     phi_out = info$index$phi_out,
                     b_out = info$index$b_out,
                     IC_out = info$index$IC_out,
                     IB_out = info$index$IB_out,
                     ID_out = info$index$ID_out,
                     spz_rate = info$index$spz_rate,
                     eff_moz_pop = info$index$eff_moz_pop,
                     moz2human_ratio = info$index$moz2human_ratio))
    }
    compare_funct <- compare_pg

  }else if(comparison=='pgmg'){
    ##Output from particle filter
    ##    run: output used for likelihood calculation
    ##    state: output used for visualization
    index <- function(info) {
      list(run = c(prev_pg = info$index$prev_preg_pg,
                   prev_mg = info$index$prev_preg_mg),
           state = c(prev_05 = info$index$prev,
                     prev_pg = info$index$prev_preg_pg,
                     prev_mg = info$index$prev_preg_mg,
                     EIR = info$index$EIR_out,
                     betaa = info$index$betaa_out,
                     clininc_all = info$index$inc,
                     prev_all = info$index$prevall,
                     clininc_05 = info$index$inc05,
                     Dout = info$index$Dout,
                     Aout = info$index$Aout,
                     Uout = info$index$Uout,
                     p_det_out = info$index$p_det_out,
                     phi_out = info$index$phi_out,
                     b_out = info$index$b_out,
                     IC_out = info$index$IC_out,
                     IB_out = info$index$IB_out,
                     ID_out = info$index$ID_out,
                     spz_rate = info$index$spz_rate,
                     eff_moz_pop = info$index$eff_moz_pop,
                     moz2human_ratio = info$index$moz2human_ratio))
    }
    compare_funct <- compare_pgmg

  } else if(comparison=='pgsg'){
    index <- function(info) {
      list(run = c(prev_pg = info$index$prev_preg_pg,
                   prev_mg = info$index$prev_preg_sg),
           state = c(prev_05 = info$index$prev,
                     prev_pg = info$index$prev_preg_pg,
                     prev_sg = info$index$prev_preg_sg,
                     EIR = info$index$EIR_out,
                     betaa = info$index$betaa_out,
                     clininc_all = info$index$inc,
                     prev_all = info$index$prevall,
                     clininc_05 = info$index$inc05,
                     Dout = info$index$Dout,
                     Aout = info$index$Aout,
                     Uout = info$index$Uout,
                     p_det_out = info$index$p_det_out,
                     phi_out = info$index$phi_out,
                     b_out = info$index$b_out,
                     IC_out = info$index$IC_out,
                     IB_out = info$index$IB_out,
                     ID_out = info$index$ID_out,
                     spz_rate = info$index$spz_rate,
                     eff_moz_pop = info$index$eff_moz_pop,
                     moz2human_ratio = info$index$moz2human_ratio))
    }
    compare_funct <- compare_pgmg
  } else if(comparison=='sg'){
    index <- function(info) {
      list(run = c(prev_mg = info$index$prev_preg_sg),
           state = c(prev_05 = info$index$prev,
                     prev_pg = info$index$prev_preg_pg,
                     prev_sg = info$index$prev_preg_sg,
                     EIR = info$index$EIR_out,
                     betaa = info$index$betaa_out,
                     clininc_all = info$index$inc,
                     prev_all = info$index$prevall,
                     clininc_05 = info$index$inc05,
                     Dout = info$index$Dout,
                     Aout = info$index$Aout,
                     Uout = info$index$Uout,
                     p_det_out = info$index$p_det_out,
                     phi_out = info$index$phi_out,
                     b_out = info$index$b_out,
                     IC_out = info$index$IC_out,
                     IB_out = info$index$IB_out,
                     ID_out = info$index$ID_out,
                     spz_rate = info$index$spz_rate,
                     eff_moz_pop = info$index$eff_moz_pop,
                     moz2human_ratio = info$index$moz2human_ratio))
    }
    compare_funct <- compare_mg
  } else if(comparison=='mg'){
    index <- function(info) {
      list(run = c(prev_mg = info$index$prev_preg_mg),
           state = c(prev_05 = info$index$prev,
                     prev_pg = info$index$prev_preg_pg,
                     prev_mg = info$index$prev_preg_mg,
                     EIR = info$index$EIR_out,
                     betaa = info$index$betaa_out,
                     clininc_all = info$index$inc,
                     prev_all = info$index$prevall,
                     clininc_05 = info$index$inc05,
                     Dout = info$index$Dout,
                     Aout = info$index$Aout,
                     Uout = info$index$Uout,
                     p_det_out = info$index$p_det_out,
                     phi_out = info$index$phi_out,
                     b_out = info$index$b_out,
                     IC_out = info$index$IC_out,
                     IB_out = info$index$IB_out,
                     ID_out = info$index$ID_out,
                     spz_rate = info$index$spz_rate,
                     eff_moz_pop = info$index$eff_moz_pop,
                     moz2human_ratio = info$index$moz2human_ratio))
    }
    compare_funct <- compare_mg
  } else if(comparison=='ancall'){
    index <- function(info) {
      list(run = c(prev_anc_all = info$index$prev_preg_all),
           state = c(prev_05 = info$index$prev,
                     prev_pg = info$index$prev_preg_pg,
                     prev_sg = info$index$prev_preg_sg,
                     prev_anc_all = info$index$prev_preg_all,
                     EIR = info$index$EIR_out,
                     betaa = info$index$betaa_out,
                     clininc_all = info$index$inc,
                     prev_all = info$index$prevall,
                     clininc_05 = info$index$inc05,
                     Dout = info$index$Dout,
                     Aout = info$index$Aout,
                     Uout = info$index$Uout,
                     p_det_out = info$index$p_det_out,
                     phi_out = info$index$phi_out,
                     b_out = info$index$b_out,
                     IC_out = info$index$IC_out,
                     IB_out = info$index$IB_out,
                     ID_out = info$index$ID_out,
                     spz_rate = info$index$spz_rate,
                     eff_moz_pop = info$index$eff_moz_pop,
                     moz2human_ratio = info$index$moz2human_ratio))
    }
    compare_funct <- compare_ancall
  }

  ##Tune particles
  if(particle_tune){
    control_tune <- mcstate::pmcmc_control(
      n_steps=100,
      save_state = FALSE,
      save_trajectories = FALSE,
      progress = TRUE,
      n_chains = n_chains, #TO DO: Make parameter to easily change
      n_workers = n_workers, #TO DO: Make parameter to easily change
      n_threads_total = n_threads,
      rerun_every = 1,
      rerun_random = FALSE)

    init_state$fixed_vol <- 1
    init_state$fixed_betaa <- init_state$betaa_eq

    min <- 10
    max <- 200
    test <- union(min * 2**(seq(0, floor(log2(max / min)))), max)
    found_good <- FALSE
    id <- 0
    target_variance <- 1

    while (!found_good && (id + 1) < length(test)) {
      id <- id + 1
      pf_tune <- mcstate::particle_filter$new(data, model, n_particles=test[id], compare_funct,
                                                    index = index,
                                                    stochastic_schedule = stochastic_schedule,
                                                    ode_control = dust::dust_ode_control(max_steps = max_steps, atol = atol, rtol = rtol),
                                                    n_threads = n_threads)
      pmcmc_run <- mcstate::pmcmc(mcmc_pars, pf_tune, control = control_tune)

      var_loglik <- stats::var(pmcmc_run$probabilities[,'log_likelihood'])
      print(pmcmc_run$probabilities)
      print(pmcmc_run$probabilities[,'log_likelihood'])
      print(var_loglik)
      if (is.na(var_loglik)) var_loglik <- 0

      message(
        date(), " ", test[id], " particles, loglikelihod variance: ",
        var_loglik
      )

      if (var_loglik < target_variance) {
        ## choose smallest var-loglikelihood < target_variance
        found_good <- TRUE
      }
    }
    if (!found_good) id <- length(test) #Take largest number of particles if still not hitting variance target
    n_particles <- test[id]
  }

  set.seed(seed) #To reproduce pMCMC results

  ### Set particle filter
  # print('about to set up particle filter')
  # print('set up particle filter')
  pf <- mcstate::particle_filter$new(data, model, n_particles, compare_funct,
                                     index = index, seed = seed,
                                     stochastic_schedule = stochastic_schedule,
                                     ode_control = dust::dust_ode_control(max_steps = max_steps, atol = atol, rtol = rtol,debug_record_step_times=FALSE),
                                     n_threads = n_threads)

  # print('about to set up pmcmc control')
  ### Set pmcmc control
  control <- mcstate::pmcmc_control(
    n_steps,
    save_state = TRUE,
    save_trajectories = TRUE,
    progress = TRUE,
    n_chains = n_chains, #TO DO: Make parameter to easily change
    n_workers = n_workers, #TO DO: Make parameter to easily change
    n_threads_total = n_threads,
    rerun_every = 50, #Re-runs particle filter about every 50 steps (random distribution, mean=50)
    rerun_random = TRUE)
  # print('set up pmcmc control')

  # print('parameters set')
  # print('starting pmcmc run')
  # print(mpl_pf)
  ### Run pMCMC
  start.time <- Sys.time()
  pmcmc_run <- mcstate::pmcmc(mcmc_pars, pf, control = control)
  run_time <- difftime(Sys.time(),start.time,units = 'secs')
  print(run_time)

  # if(n_chains > 1) pmcmc_run <- mcstate::pmcmc_combine(pmcmc_run)

  pars <- pmcmc_run$pars
  probs <- pmcmc_run$probabilities
  mcmc <- coda::as.mcmc(cbind(probs, pars))

  ##Save seasonality equilibrium trajectories if checking equilibrium
  seas_pretime <- NULL
  if(seasonality_on & seasonality_check){
    print('Saving seasonality equilibrium trajectories')
    # Create list of seasonality trajectories for each set of sampled parameters in the posterior
    seas_pretime <- parallel::mclapply(1:nrow(pars), function(x) check_seasonality(theta=pars[x,],mpl_pf=mpl_pf,det_model=det_model))
  }
  to_return <- list(threads = n_threads,
                    particles = n_particles,
                    run_time = run_time,
                    mcmc = as.data.frame(mcmc),
                    pars = as.data.frame(pars),
                    probs = as.data.frame(probs),
                    times = pmcmc_run$trajectories$time,
                    history = pmcmc_run$trajectories$state,
                    seas_history = seas_pretime)

  return(to_return)
}

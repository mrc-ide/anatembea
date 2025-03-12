#' Test the particle filter
#'
#' This function sets up and runs a particle MCMC that uses Dust, Odin and MCState
#'
#' @param data_raw Time series data to fit model
#' @param data_raw_pg Time series of primigravidae ANC data to fit model
#' @param data_raw_mg Time series of multigravidae ANC data to fit model
#' @param n_particles Number of particles to be used in pMCMC (default = 200)
#' @param proposal_matrix Proposal matrix for MCMC parameters
#' @param max_param Ceiling for proposed stochastic parameter (either EIR or betaa) values (default = 1000)
#' @param max_steps Maximum steps for particle filter (default = 1e7)
#' @param atol atol for particle filter (default = 1e-3)
#' @param rtol rtol for particle filter (default = 1e-6)
#' @param n_steps Number of MCMC steps in a single chain (default = 500)
#' @param n_threads Number of processing threads (default = 4)
#' @param lag_rates Number of delay compartments (default = 10)
#' @param state_check Run equilibrium checks, if state_check = 1, returns expected deriv values which should equal 0 and sets stochastic model to have EIR constant at init_EIR
#'                    If state_check = 1 and seasonality_on = 1, then the deterministic seasonal model is still run, but theta2 is forced to 1, forcing a constant seasonality profile
#'                    If state_check = 0, no values are printed
#' @param country Name of country (needed if using seasonality model)
#' @param admin_unit Name of administrative unit (needed if using seasonality model)
#' @param preyears Length of time in years the deterministic seasonal model should run before Jan 1 of the year observations began (default = 2)
#' @param seasonality_on Toggle seasonality model run before observed time period (default = 1)
#' @param seasonality_check Toggle saving values of seasonality equilibrium (default = 1)
#' @param seed Allows user to specify a seed (default = 1L)
#' @param start_pf_time Number of days before first observation that particle filter will start (default = 30)
#' @param comparison The comparison function to be used. Either 'u5' which
#'          equates the observed prevalence to prevalence under 5 years old in
#'          the model or 'pgmg' which calculates prevalence in primigravid and
#'          multigravid pregnant women for comparison with observed ANC data. c('u5','pg','sg','mg','pgmg','pgsg','ancall')
#' @param stoch_param Which parameter to test

#' @export
test_pf <- function(data_raw=NULL,
                    data_raw_pg=NULL,
                    data_raw_mg=NULL,
                    n_particles=200,
                    proposal_matrix,
                    max_param=1000,
                    # EIR_vol,
                    # proposal_dist,
                    # init_EIR = 100,
                    max_steps = 1e7,
                    atol = 1e-6,
                    rtol = 1e-6,
                    n_steps = 500,
                    n_threads = 4,
                    lag_rates = 10,
                    state_check = 0,## Run equilibrium checks
                    # If state_check = 1, returns expected deriv values which should equal 0 and sets stochastic model to have EIR constant at init_EIR
                    # If state_check = 1 and seasonality_on = 1, then the deterministic seasonal model is still run, but theta2 is forced to 1, forcing a constant seasonality profile
                    # If state_check = 0, no values are printed
                    country = NULL,
                    admin_unit = NULL,
                    preyears = 2, #Length of time in years the deterministic seasonal model should run before Jan 1 of the year observations began
                    seasonality_on = 1,  ## state_check = 1 runs a deterministic seasonal model before running the stochastic model to get more realistic immunity levels
                    ## If seasonality_on = 0, runs the stochastic model based on the standard equilibrium solution
                    seasonality_check = 0,##If 1, saves values of seasonality equilibrium
                    seed = 1L,
                    start_pf_time = 30,
                    stoch_param = c('EIR','betaa'),
                    comparison = c('u5','pgmg')){
  ## Modify dates from data
  # str(data_raw$month)
  # print(data_raw$month)
  ##Merge primigrav and multigrav datasets if necessary.
  if(comparison=='pgmg'){
    data_raw <- dplyr::left_join(data_raw_pg,data_raw_mg,by=c('month','t'),suffix = c('.pg','.mg'))
  }
  start_obs <- min(zoo::as.Date(zoo::as.yearmon(data_raw$month)))#Month of first observation (in Date format)
  time_origin <- zoo::as.Date(paste0(lubridate::year(start_obs)-1,'-01-01')) #January 1 of year before observation (in Date format)
  data_raw_time <- data_raw
  data_raw_time$date <- zoo::as.Date(zoo::as.yearmon(data_raw$month), frac = 0.5) #Convert dates to middle of month
  data_raw_time$t <- as.integer(difftime(data_raw_time$date,time_origin,units="days")) #Calculate date as number of days since January 1 of year before observation
  initial_time <- min(data_raw_time$t) - start_pf_time #Start particle filter a given time (default = 30d) before first observation
  #Create daily sequence from initial_time to end of observations
  #This is so the trajector histories return daily values (otherwise it returns
  #model values only at the dates of observation)
  time_list <- data.frame(t=initial_time:max(data_raw_time$t))
  data_raw_time <- dplyr::left_join(time_list,data_raw_time,by='t')
  start_stoch <- zoo::as.Date(start_obs - start_pf_time) #Start of stochastic schedule; needs to start when particle filter starts
  data <- mcstate::particle_filter_data(data_raw_time, time = "t", rate = NULL, initial_time = initial_time) #Declares data to be used for particle filter fitting
  # print('Data processed')

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
                   spz_rate = info$index$spz_rate))
  }

  ## Provide schedule for changes in stochastic process (in this case EIR)
  ## Converts a sequence of dates (from start_stoch to 1 month after last observation point) to days since January 1 of the year before observation
  stoch_update_dates <- seq.Date(start_stoch,max(as.Date(data_raw_time$date+30),na.rm = TRUE),by='month')
  stochastic_schedule <- as.integer(difftime(stoch_update_dates,time_origin,units="days"))
  # print('stochastic_schedule assigned')

  #Provide age categories, proportion treated, and number of heterogeneity brackets
  init_age <- c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3.5, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80)
  prop_treated <- 0.4
  het_brackets <- 5

  #Create model parameter list. Also loads seasonality profile data file to match to desired admin_unit and country
  mpl_pf <- model_param_list_create(init_age = init_age,
                                    prop_treated = prop_treated,
                                    het_brackets = het_brackets,
                                    max_param = max_param,
                                    state_check = state_check,
                                    lag_rates = lag_rates,
                                    country = country,
                                    admin_unit = admin_unit,
                                    start_stoch = start_stoch,
                                    time_origin = time_origin,
                                    seasonality_on = seasonality_on,
                                    preyears = preyears)
  # print('model parameter list created')
  # print(mpl_pf$state_check)
  # print(mpl_pf$ssa0)

  ## If a deterministic seasonal model is needed prior to the stochastic model, this loads the deterministic odin model
  if(seasonality_on == 1){
    odin_det <- system.file("odin", "odin_model_stripped_seasonal.R", package = "anatambea")
    season_model <- odin::odin(odin_det)
  }


  ## Load stochastic model in odin.dust
  # print('about to load stochastic model')
  ##Switch between EIR and mosquito emergence models
  stoch_file <- ifelse(stoch_param=='betaa','odinmodelmatchedstoch_mozemerg.R','odinmodelmatchedstoch.R')
  odin_stoch <- system.file("odin", stoch_file, package = "anatambea")
  model <- odin.dust::odin_dust(odin_stoch)
  # print('loaded stochastic model')

  set.seed(seed) #To reproduce pMCMC results

  if(comparison=='u5'){
    pf <- mcstate::particle_filter$new(data, model, n_particles, compare_u5,
                                       index = index, seed = seed,
                                       stochastic_schedule = stochastic_schedule,
                                       ode_control = dust::dust_ode_control(max_steps = max_steps, atol = atol, rtol = rtol),
                                       n_threads = n_threads)
  }
  else if(comparison=='pgmg'){
    mpl_pf <- append(mpl_pf,list(coefs_pg_df = as.data.frame(readRDS('./inst/extdata/pg_corr_sample.RDS')),
                                 coefs_mg_df = as.data.frame(readRDS('./inst/extdata/mg_corr_sample.RDS'))))
    pf <- mcstate::particle_filter$new(data, model, n_particles, compare_pgmg,
                                       index = index, seed = seed,
                                       stochastic_schedule = stochastic_schedule,
                                       ode_control = dust::dust_ode_control(max_steps = max_steps, atol = atol, rtol = rtol, debug_record_step_times=TRUE),
                                       n_threads = n_threads)

  }
  # print('set up particle filter')
   # print('set up particle filter')

  # # print('about to set up pmcmc control')
  # ### Set pmcmc control
  # control <- mcstate::pmcmc_control(
  #   n_steps,
  #   save_state = TRUE,
  #   save_trajectories = TRUE,
  #   progress = TRUE,
  #   n_chains = 1, #TO DO: Make parameter to easily change
  #   n_workers = 1, #TO DO: Make parameter to easily change
  #   n_threads_total = n_threads,
  #   rerun_every = 50, #Re-runs particle filter about every 50 steps (random distribution, mean=50)
  #   rerun_random = TRUE)
  # # print('set up pmcmc control')

  ### Set pmcmc parameters
  proposals <- data.frame(volatility = rlnorm(1, meanlog = -.2, sdlog = 0.5),
                      log_init_EIR = rnorm(1, mean = 0, sd = 10))

  transform_pars <- transform(mpl_pf,season_model)
  # print('parameters set')
  # print('starting pmcmc run')

  ### Run pMCMC
  run <- lapply(1, function(x,proposal){
    equil <- transform_pars(proposal[x,])
    pf$run(pars=equil,save_history = TRUE)
  # Our simulated trajectories, with the "real" data superimposed
    return(list(history=pf$history(),stats=pf$ode_statistics()))
  },proposal=proposals)
  # start.time <- Sys.time()
  # pmcmc_run <- mcstate::pmcmc(mcmc_pars, pf, control = control)
  # run_time <- difftime(Sys.time(),start.time,units = 'secs')
  # print(run_time)
  # pars <- pmcmc_run$pars
  # probs <- pmcmc_run$probabilities
  # mcmc <- coda::as.mcmc(cbind(probs, pars))
  #

  return(run)
}

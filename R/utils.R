#' load_file
#'
#' \code{load_file} loads package file
#'
#' @description Load a file from within the inst/extdata folder of the
#'   mamasante package. File extension must be one of .csv, .txt, or .rds.
#'
#' @param name the name of a file within the inst/extdata folder.
#'
#' @importFrom utils read.csv
#'
#' @export

load_file <- function(name) {

  # check that valid file extension
  ext <- strsplit(name, "\\.")[[1]]
  ext <- ext[length(ext)]
  if(is.element(ext, c("csv", "rds", "RDS")) == FALSE){
    stop("file extension not valid")
  }

  # get full file path
  name_full <- system.file("extdata/", name, package="mamasante", mustWork = TRUE)
  message(name_full)

  # read in file
  if (ext == "rds" | ext == "RDS") {
    ret <- readRDS(name_full)
  } else {
    ret <-  read.csv(file=name_full, header=TRUE, sep=",")
  }

  return(ret)
}

#------------------------------------------------
#' match clean
#'
#' @param a First string to compare. Default = NULL
#' @param b Second string to compare. Default = NULL
#'
#' @importFrom RecordLinkage levenshteinSim
#'
#' @export

match_clean <- function(a, b){

  a <- gsub("[[:punct:][:space:]]", "", tolower(stringi::stri_trans_general(a, "latin-ascii")))
  b <- gsub("[[:punct:][:space:]]", "", tolower(stringi::stri_trans_general(b, "latin-ascii")))

  ret <- which(b %in% a)

  if(length(ret) == 0){
    distance <- levenshteinSim(a, b)
    ret <- which.max(distance)
  }
  return(ret)
}

#------------------------------------------------
#' match admin region
#'
#' \code{admin_match} Matches the user input admin unit and country with data
#'
#' @param country Character for country within which admin unit is in.
#'   Default = NULL
#' @param admin_unit Character for admin region. Some fuzzy logic will be used to
#'   match. If not provided then no seasonality is introduced. Default = NULL
#' @param admin_units_seasonal Dataframe of seasonality data for country and admin unit
#'
#' @export

admin_match <- function(admin_unit = NULL, country = NULL,
                        admin_units_seasonal){

  # intialise admin match as no match
  admin_matches <- 0

  if (!is.null(admin_unit)) {

    # if there is no country given then search for the admin unit
    if (is.null(country)) {

      # find exact match
      admin_matches <- which(tolower(admin_units_seasonal$admin1) %in% tolower(admin_unit))

      # if exact does not match try closest match
      if (length(admin_matches) < 1) {
        admin_matches <- match_clean(admin_unit, admin_units_seasonal$admin1)
      } else if (length(admin_matches) > 1){
        stop('Please specify the country of admin unit.  There are multiple with same name.')
      }

      # if we do have a country though find that match first and then find admin
    } else {

      # first find an exact match
      country_matches <- which(tolower(admin_units_seasonal$country) %in% tolower(country))

      if (length(country_matches) < 1) {
        country_name <- admin_units_seasonal$country[match_clean(country, admin_units_seasonal$country)]
        country_matches <- which(tolower(admin_units_seasonal$country) %in% tolower(country_name))
      }

      sub_admin_units_seasonal <- admin_units_seasonal[country_matches,]

      # find exact match
      admin_sub_matches <- which(tolower(sub_admin_units_seasonal$admin1) %in% tolower(admin_unit))

      # if exact does not match try closest match
      if (length(admin_sub_matches) != 1) {
        admin_sub_matches <- match_clean(admin_unit,
                                         sub_admin_units_seasonal$admin1)
      } else if (length(admin_sub_matches) > 1){
        stop('There are multiple admins with same name - check the data file!')
      }

      admin_matches <- country_matches[admin_sub_matches]
    }

    message("Requested: ", admin_unit,
            "\nReturned: ", admin_units_seasonal$admin1[admin_matches], ", ",
            admin_units_seasonal$country[admin_matches])
  }

  return(admin_matches)
}

#------------------------------------------------
#' Change missing data to 0s so dbinom function returns 0 (instead of NA)
#'
#' \code{format_na} Change missing data to 0s so dbinom function returns 0 (instead of NA)
#'
#' @param df Data input from user
#'   Default = NULL
#' @export
format_na <- function(df){
  df$tested <- ifelse(is.na(df$tested),0,df$tested)
  df$positive <- ifelse(df$tested==0,0,df$positive)
  return(df)
}
#------------------------------------------------
#' Process raw data for use in pMCMC
#'
#' \code{data_process} Transform data file
#'
#' @param data_raw Data input from user
#'   Default = NULL
#'@param start_pf_time Amount of time before first observation to start particle filter
#'@param check_flexibility Save pre-observation run
#' @export

data_process <- function(data_raw = NULL, start_pf_time = NULL, check_flexibility){
  ## Modify dates from data
  start_obs <- min(zoo::as.Date(zoo::as.yearmon(data_raw$month)))#Month of first observation (in Date format)
  time_origin <- zoo::as.Date(paste0(lubridate::year(start_obs)-1,'-01-01')) #January 1 of year before observation (in Date format)
  data_raw_time <- data_raw
  data_raw_time$date <- zoo::as.Date(zoo::as.yearmon(data_raw$month), frac = 0.5) #Convert dates to middle of month
  data_raw_time$t <- as.integer(difftime(data_raw_time$date,time_origin,units="days")) #Calculate date as number of days since January 1 of year before observation
  initial_time <- min(data_raw_time$t) - start_pf_time #Start particle filter a given time (default = 30d) before first observation

  #Create daily sequence from initial_time to end of observations
  #This is so the trajector histories return daily values (otherwise it returns
  #model values only at the dates of observation)
  start_stoch <- lubridate::floor_date(zoo::as.Date(min(data_raw_time$date) - start_pf_time), unit='month') #Start of stochastic schedule; needs to start when particle filter starts
  ## Provide schedule for changes in stochastic process (in this case EIR)
  ## Converts a sequence of dates (from start_stoch to 1 month after last observation point) to days since January 1 of the year before observation
  stoch_update_dates <- c(seq.Date(start_stoch,start_obs,by='month'),seq.Date(start_obs,max(as.Date(data_raw_time$date+30),na.rm = TRUE),by='month')[-1])
  stochastic_schedule <- as.integer(difftime(stoch_update_dates,time_origin,units="days"))
  if(check_flexibility){
    flex_time_dates <- zoo::as.Date(zoo::as.yearmon(seq.Date(start_stoch,start_obs,by='month')),frac=0.5)
    dummy_dates <- data.frame(date=c(start_stoch,flex_time_dates[-length(flex_time_dates)]))
    data_raw_time <- dplyr::bind_rows(dummy_dates,data_raw_time)
    data_raw_time$t <- as.integer(difftime(data_raw_time$date,time_origin,units="days")) #Calculate date as number of days since January 1 of year before observation
    initial_time <- min(data_raw_time$t)
  }
  data <- mcstate::particle_filter_data(data_raw_time, time = "t", rate = NULL, initial_time = initial_time) #Declares data to be used for particle filter fitting




  return(list(data=data,stochastic_schedule=stochastic_schedule,start_stoch=start_stoch,time_origin=time_origin))
}

#------------------------------------------------
#' transform final state of a seasonal model into the initial state of a stochastic model
#'
#' \code{transform_init} Transform final state of a seasonal model into the initial state of a stochastic model
#'
#' @param final_state Final time point of a seasonal model
#'   Default = NULL
#'
#' @export

transform_init <- function(final_state = NULL){
  f <- final_state
  last <- length(f$FOI_init[,1,1])
  list(FOI_eq = f$FOI_init[last,,],
       init_EIR = f$EIR_init[last,1,1],
       init_S = f$S_init[last,,],
       init_T = f$T_init[last,,],
       init_D = f$D_init[last,,],
       init_A = f$A_init[last,,],
       init_U = f$U_init[last,,],
       init_P = f$P_init[last,,],
       init_IB = f$IB_init[last,,],
       init_ID = f$ID_init[last,,],
       init_ICA = f$ICA_init[last,,],
       ICM_age = f$ICM_age_init[last,],
       age_rate = f$age_rate_init[last,],
       het_wt = f$het_wt_init[last,],
       foi_age = f$foi_age_init[last,],
       rel_foi = f$rel_foi_init[last,],
       na = f$na_init[last],
       nh = f$nh_init[last],
       x_I = f$x_I_init[last,],
       omega = f$omega_init[last],
       den = f$den_init[last,],
       age59 = f$age59_init[last],
       age05 = f$age05_init[last],
       age = f$age_init[last,],
       ft = f$ft_init[last],
       age20l = f$age20l_init[last],
       age20u = f$age20u_init[last],
       age_20_factor = f$age_20_factor_init[last],
       eta = f$eta_init[last],
       rA = f$rA_init[last],
       rT = f$rT_init[last],
       rD = f$rD_init[last],
       rU = f$rU_init[last],
       rP = f$rP_init[last],
       uCA = f$uCA_init[last],
       dCA = f$dCA_init[last],
       dB = f$dB_init[last],
       uB = f$uB_init[last],
       dID = f$dID_init[last],
       uD = f$uD_init[last],
       PM = f$PM_init[last],
       phi0 = f$phi0_init[last],
       phi1 = f$phi1_init[last],
       IC0 = f$IC0_init[last],
       kC = f$kC_init[last],
       b0 = f$b0_init[last],
       b1 = f$b1_init[last],
       kB = f$kB_init[last],
       IB0 = f$IB0_init[last],
       aD = f$aD_init[last],
       fD0 = f$fD0_init[last],
       gammaD = f$gammaD_init[last],
       d1 = f$d1_init[last],
       ID0 = f$ID0_init[last],
       kD = f$kD_init[last],
       dE = f$dE_init[last],
       DY = f$DY_init[last],
       lag_rates = f$lag_rates_init[last],
       Q0 = f$Q0_init[last],
       state_check = f$state_check_init[last],
       tau1 = f$tau1_init[last],
       tau2 = f$tau2_init[last],
       init_Sv = f$Sv_init[last],
       init_Ev = f$Ev_init[last],
       init_Iv = f$Iv_init[last],
       mv0 = f$mv_init[last],
       betaa_eq = f$betaa_eq[last],
       FOIv_eq = f$FOIv_init[last],
       prev = f$prev[last],
       init_EL = f$EL_init[last],
       init_LL = f$LL_init[last],
       init_PL = f$PL_init[last],
       EIR_eq = f$EIR_eq[last,,],
       pi = f$pi_eq[last],
       prev05 = f$prev[last],
       inc = f$inc[last],
       inc05 = f$inc05[last],
       delayGam = f$delayGam_eq[last]


  )
}

#------------------------------------------------
#' Compare function to calculate likelihood
#'
#' \code{compare_u5} Compare function that compares observed data with model estimate
#'  to calculate likelihood for the particle filter. Equates observed prevalence
#'  to prevalence under 5 year olds in the model.
#'
#' @param state Model output. Default = NULL
#' @param observed Oberved data. Default = NULL
#' @param pars Parameters, optional.
#'
#' @export
compare_u5 <- function(state, observed, pars = NULL) {
  # print('in compare function')
  # cat('Observed positive:\n',observed$positive,'\n')
  # cat('Observed tested:\n',observed$tested,'\n')
  # cat('Estimated prevalence:\n',state[1,],'\n')
  ll <- dbinom(x = observed$positive,
         size = observed$tested,
         prob = state[1,],
         log = TRUE)
  # print(ll)
  return(ll)
}
#------------------------------------------------
#' Compare function to calculate likelihood
#'
#' \code{compare_pg} Compare function that compares observed data with model estimate
#'  to calculate likelihood for the particle filter. Converts under 5 year old
#'  prevalence in the model to prevalence among primigravid
#'  women based on coefficients from a separate regression analysis.
#'
#' @param state Model output. Default = NULL
#' @param observed Oberved data. Default = NULL
#' @param pars Parameters, optional.
#'
#' @export
compare_pg <- function(state, observed, pars = NULL) {

  ll_pg <- dbinom(x = observed$positive,
                  size = observed$tested,
                  prob = state['prev_pg',],
                  log = TRUE)

  return(ll_pg)
}

#------------------------------------------------
#' Compare function to calculate likelihood
#'
#' \code{compare_mg} Compare function that compares observed data with model estimate
#'  to calculate likelihood for the particle filter. Converts under 5 year old
#'  prevalence in the model to prevalence among multigravid
#'  women based on coefficients from a separate regression analysis.
#'
#' @param state Model output. Default = NULL
#' @param observed Oberved data. Default = NULL
#' @param pars Parameters, optional.
#'
#' @export
compare_mg <- function(state, observed, pars = NULL) {

  ll_mg <- dbinom(x = observed$positive,
                  size = observed$tested,
                  prob = state['prev_mg',],
                  log = TRUE)

  return(ll_mg)
}
#------------------------------------------------
#' Compare function to calculate likelihood
#'
#' \code{compare_pgmg} Compare function that compares observed data with model estimate
#'  to calculate likelihood for the particle filter. Converts under 5 year old
#'  prevalence in the model to prevalence among primigravid and multigravid
#'  women based on coefficients from a separate regression analysis.
#'
#' @param state Model output. Default = NULL
#' @param observed Oberved data. Default = NULL
#' @param pars Parameters, optional.
#'
#' @export
compare_pgmg <- function(state, observed, pars = NULL) {

  ll_pg <- dbinom(x = observed$positive.pg,
                          size = observed$tested.pg,
                          prob = state['prev_pg',],
                          log = TRUE)

  ll_mg <- dbinom(x = observed$positive.mg,
                          size = observed$tested.mg,
                          prob = state['prev_mg',],
                          log = TRUE)

  return(ll_pg+ll_mg)
}
#------------------------------------------------
#' Compare function to calculate likelihood
#'
#' \code{compare_ancall} Compare function that compares observed data with model estimate
#'  to calculate likelihood for the particle filter. Converts under 5 year old
#'  prevalence in the model to prevalence among all
#'  women based on coefficients from a separate regression analysis (TZ data only).
#'
#' @param state Model output. Default = NULL
#' @param observed Oberved data. Default = NULL
#' @param pars Parameters, optional.
#'
#' @export
compare_ancall <- function(state, observed, pars = NULL) {
  #skip comparison if data is missing
  if(is.na(observed$positive)) {
    return(numeric(length(state[1,])))}
  ll <- dbinom(x = observed$positive,
               size = observed$tested,
               prob = state['prev_anc_all',],
               log = TRUE)
  return(ll)
}
#------------------------------------------------
#' Estimate the initial state given user inputs
#'
#' \code{initialise} Calculates the model equilibrium based on an initial EIR values,
#'    then optionally runs a deterministic seasonal model and returns initial
#'    values to be used for the stochastic model fitting.
#'
#' @param init_EIR User supplied initial EIR
#' @param mpl Model parameter list. Default = NULL
#' @param det_model Seasonality model to be used for the optional deterministic
#'    model. Default = NULL
#'
#' @export
initialise <- function(init_EIR,mpl,det_model){
  EIR_vals <- NULL
  EIR_times <- NULL
  ##Convert EIR dataframe into vectors for use in the piecewise deterministic models.
  if(is.data.frame(init_EIR)){
    # print('Input EIR: ')
    # print(init_EIR)
    EIR_vals <- init_EIR$EIR
    EIR_times <- init_EIR$time
    init_EIR <- EIR_vals[1]
    # print('Output EIR: ')
    # print(EIR_vals)
    # cat('New initial EIR: ',init_EIR,'\n')
  }
  if(!is.null(mpl$target_prev)){
    message('Optimizing initial EIR based on target prevalence.')
    opt_EIR <- suppressMessages(stats::optim(1,fn=mamasante::get_init_EIR,mpl=mpl,method='Brent',lower=0,upper=2000))
    init_EIR <- opt_EIR$par
  }
  # else if(!is.null(mpl$target_prev)&mpl$target_prev_group=='u5'){
  #   message('Optimizing initial EIR based on target prevalence in under 5 year olds.')
  #   opt_EIR <- suppressMessages(stats::optim(1,fn=mamasante::get_init_EIR_u5,mpl=mpl,age_group='u5',method='Brent',lower=0,upper=2000))
  #   init_EIR <- opt_EIR$par
  # }
  message(paste0('Initial EIR set to ',round(init_EIR,digits=1),'.'))

  ## Run equilibrium function
  state <- equilibrium_init_create_stripped(init_EIR = init_EIR,
                                            model_param_list = mpl,
                                            age_vector = mpl$init_age,
                                            ft = mpl$prop_treated,
                                            het_brackets = mpl$het_brackets,
                                            state_check = mpl$state_check)
  state <- append(state,list(EIR_times=EIR_times,EIR_vals=EIR_vals))
  # print('equilibrium state calculated')
  # print(state)

  ##run deterministic seasonality model first if seasonality_on == 1
  if(!is.null(det_model)){
    # print('creating seasonality equilirium')
    #Keep only necessary parameters
    state_use <- state[names(state) %in% coef(det_model)$name]
    # create model with initial values
    mod <- det_model$new(user = state_use, use_dde = TRUE)

    # Define time length of the deterministic model run (preyears = how many years you want the deterministic model to run before the particle filter)
    # deterministic_stop: defines the day the seasonality model should stop so that
    #                      the particle filter begins at the right time of year
    deterministic_stop <- as.integer(difftime(mpl$start_stoch,mpl$time_origin,units="days"))
    # print(mpl$preyears*365+deterministic_stop)
    tt <- seq(deterministic_stop, mpl$preyears*365+deterministic_stop,length.out=100)
    # cat('pre-model time:\n',tt,'\n')
    # run seasonality model
    mod_run <- mod$run(tt, verbose=FALSE,step_size_max=9)

    # shape output
    out <- mod$transform_variables(mod_run)
    # windows(10,8)
    # plot(out$t,out$prev,type='l')
    # windows(10,8)
    # plot(out$t,out$mv_init,type='l')
    # windows(10,8)
    # plot(out$t,out$betaa_eq,type='l')
    # windows(10,8)
    # plot(out$t,out$theta2_eq,type='l')
    # windows(10,8)
    # plot(out$t,out$lag_seas_eq,type='l')

    # View(out)

    # Transform seasonality model output to match expected input of the stochastic model
    init4pmcmc <- transform_init(out)
    # print(init4pmcmc)
    # cat('prev equilibrium: ',state_use$prev,'\n')
    # cat('prev seasonal: ',init4pmcmc$prev,'\n')
    #Print some equilibrium checks if state_check==1
    if(mpl$state_check){
      print('running equilibrium checks')
      H <- sum(init4pmcmc$init_S) + sum(init4pmcmc$init_T) + sum(init4pmcmc$init_D) + sum(init4pmcmc$init_A) + sum(init4pmcmc$init_U) + sum(init4pmcmc$init_P)

      deriv_S11 <- -init4pmcmc$FOI_eq[1,1]*init4pmcmc$init_S[1,1] + init4pmcmc$rP*init4pmcmc$init_P[1,1] + init4pmcmc$rU*init4pmcmc$init_U[1,1] +
        init4pmcmc$eta*H*init4pmcmc$het_wt[1] - (init4pmcmc$eta+init4pmcmc$age_rate[1])*init4pmcmc$init_S[1,1]
      print('Seasonal equilibrium')
      cat('deriv S check: ',deriv_S11,'\n')
      b <- init4pmcmc$b0 * ((1-init4pmcmc$b1)/(1+(init4pmcmc$init_IB[1,1]/init4pmcmc$IB0)^init4pmcmc$kB)+init4pmcmc$b1)
      cat('b check: ',b,'\n')
      EIR_eq11 <- init4pmcmc$init_EIR/365 * init4pmcmc$rel_foi[1] * init4pmcmc$foi_age[1]
      FOI_lag <- EIR_eq11 * (if(init4pmcmc$init_IB[1,1]==0) init4pmcmc$b0 else b)
      deriv_FOI111 <- (init4pmcmc$lag_rates/init4pmcmc$dE)*FOI_lag - (init4pmcmc$lag_rates/init4pmcmc$dE)*init4pmcmc$FOI_eq[1,1]
      cat('deriv FOI check: ',deriv_FOI111,'\n')

      print('Compare equilibrium state to seasonal equilibrium')
      print('The following values should be 0 or very close.')
      cat('S check: ',init4pmcmc$init_S-state$init_S,'\n')
      cat('T check: ',init4pmcmc$init_T-state$init_T,'\n')
      cat('D check: ',init4pmcmc$init_D-state$init_D,'\n')
      cat('A check: ',init4pmcmc$init_A-state$init_A,'\n')
      cat('U check: ',init4pmcmc$init_U-state$init_U,'\n')
      cat('P check: ',init4pmcmc$init_P-state$init_P,'\n')
      cat('Iv check: ',init4pmcmc$init_Iv-state$init_Iv,'\n')
      cat('init_EIR check: ',state$init_EIR-init4pmcmc$init_EIR,'\n')
      cat('prev check: ',state$prev-init4pmcmc$prev,'\n')

      print('Helpful values for reference:')
      cat('init_EIR: ',state$init_EIR,'\n')
      cat('Equilibrium prev: ',state$prev,'\n')
      cat('Seasonal equilibrium prev: ',init4pmcmc$prev,'\n')
      saveRDS(append(mpl,init4pmcmc),'seasonal_equil_values.RDS')
      saveRDS(state,'original_equil_values.RDS')

      # mpl['init_EIR'] <- NULL
      # View(init4pmcmc)
      # View(state_use)
    }
    return(append(mpl,init4pmcmc)) #Append all parameters from model parameter list for stochastic model
  }
  else{
    return(state)
  }
}
#------------------------------------------------
#' Transformation function that calculates initial values for stochastic model
#'
#' \code{user_informed} Calculates the model equilibrium based on an initial EIR values,
#'    then optionally runs a deterministic seasonal model and returns initial
#'    values to be used for the stochastic model fitting.
#'
#' @param init_state Initial state
#'
#' @export
user_informed <- function(init_state){ ## Wraps transformation function in a 'closure' environment so you can pass other parameters that you aren't fitting with the pMCMC
  function(theta) {
    ## theta: particle filter parameters that are being fitted (and so are changing at each MCMC step)
    # print('in transform function')
    # init_EIR <- exp(theta[["log_init_EIR"]]) ## Exponentiate EIR since MCMC samples on the log scale for EIR
    init_state$volatility <- theta[["volatility"]]
    # print(init_state$betaa_eq)
    # init_state$betaa_eq <- theta[["init_betaa"]]

    if(init_state$particle_tune) {
      init_state$betaa_eq <- init_state$fixed_betaa
      init_state$volatility <- init_state$fixed_vol
    }
    # print(init_state$betaa_eq)
    return(init_state)
  }
}
#------------------------------------------------
#' Transformation function that calculates initial values for stochastic model
#'
#' \code{data_informed} Calculates the model equilibrium based on an initial EIR values,
#'    then optionally runs a deterministic seasonal model and returns initial
#'    values to be used for the stochastic model fitting.
#'
#' @param mpl_pf Model parameter list. Default = NULL
#' @param season_model Seasonality model to be used for the optional deterministic
#'    model. Default = NULL
#'
#' @export
data_informed <- function(mpl_pf,season_model){ ## Wraps transformation function in a 'closure' environment so you can pass other parameters that you aren't fitting with the pMCMC
  function(theta) {
    ## theta: particle filter parameters that are being fitted (and so are changing at each MCMC step)
    # print('in transform function')
    init_EIR <- exp(theta[["log_init_EIR"]]) ## Exponentiate EIR since MCMC samples on the log scale for EIR
    vol <- theta[["volatility"]]
    mpl <- append(mpl_pf,list(volatility = vol)) ## Add extra MCMC parameters to model parameter list that aren't needed to calculate equilibrium

    ## Run equilibrium function
    state <- equilibrium_init_create_stripped(init_EIR = init_EIR,
                                              model_param_list = mpl,
                                              age_vector = mpl$init_age,
                                              ft = mpl$prop_treated,
                                              het_brackets = mpl$het_brackets,
                                              state_check = mpl$state_check)
    # print('equilibrium state calculated')
    # print(state)

    ##run deterministic seasonality model first if seasonality_on == 1
    if(mpl$seasonality_on){
      # print('creating seasonality equilirium')
      #Keep only necessary parameters
      state_use <- state[names(state) %in% coef(season_model)$name]

      # create model with initial values
      mod <- season_model$new(user = state_use, use_dde = TRUE)

      # Define time length of the deterministic model run (preyears = how many years you want the deterministic model to run before the particle filter)
      # deterministic_stop: defines the day the seasonality model should stop so that
      #                      the particle filter begins at the right time of year
      deterministic_stop <- as.integer(difftime(mpl$start_stoch,mpl$time_origin,units="days"))
      # print(mpl$preyears*365+deterministic_stop)
      tt <- seq(0, mpl$preyears*365+deterministic_stop,length.out=100)

      # run seasonality model
      mod_run <- mod$run(tt, verbose=FALSE,step_size_max=9)

      # shape output
      out <- mod$transform_variables(mod_run)
      # windows(10,8)
      # plot(out$t,out$prev,type='l')
      # View(out)

      # Transform seasonality model output to match expected input of the stochastic model
      init4pmcmc <- transform_init(out)
      # print(init4pmcmc)
      # cat('prev equilibrium: ',state_use$prev,'\n')
      # cat('prev seasonal: ',init4pmcmc$prev,'\n')

      #Print some equilibrium checks if state_check==1
      if(mpl$state_check){
        print('running equilibrium checks')
        H <- sum(init4pmcmc$init_S) + sum(init4pmcmc$init_T) + sum(init4pmcmc$init_D) + sum(init4pmcmc$init_A) + sum(init4pmcmc$init_U) + sum(init4pmcmc$init_P)

        deriv_S11 <- -init4pmcmc$FOI_eq[1,1]*init4pmcmc$init_S[1,1] + init4pmcmc$rP*init4pmcmc$init_P[1,1] + init4pmcmc$rU*init4pmcmc$init_U[1,1] +
          init4pmcmc$eta*H*init4pmcmc$het_wt[1] - (init4pmcmc$eta+init4pmcmc$age_rate[1])*init4pmcmc$init_S[1,1]
        print('Seasonal equilibrium')
        cat('deriv S check: ',deriv_S11,'\n')
        b <- init4pmcmc$b0 * ((1-init4pmcmc$b1)/(1+(init4pmcmc$init_IB[1,1]/init4pmcmc$IB0)^init4pmcmc$kB)+init4pmcmc$b1)
        cat('b check: ',b,'\n')
        EIR_eq11 <- init4pmcmc$init_EIR/365 * init4pmcmc$rel_foi[1] * init4pmcmc$foi_age[1]
        FOI_lag <- EIR_eq11 * (if(init4pmcmc$init_IB[1,1]==0) init4pmcmc$b0 else b)
        deriv_FOI111 <- (init4pmcmc$lag_rates/init4pmcmc$dE)*FOI_lag - (init4pmcmc$lag_rates/init4pmcmc$dE)*init4pmcmc$FOI_eq[1,1]
        cat('deriv FOI check: ',deriv_FOI111,'\n')

        print('Compare equilibrium state to seasonal equilibrium')
        print('The following values should be 0 or very close.')
        cat('S check: ',init4pmcmc$init_S-state$init_S,'\n')
        cat('T check: ',init4pmcmc$init_T-state$init_T,'\n')
        cat('D check: ',init4pmcmc$init_D-state$init_D,'\n')
        cat('A check: ',init4pmcmc$init_A-state$init_A,'\n')
        cat('U check: ',init4pmcmc$init_U-state$init_U,'\n')
        cat('P check: ',init4pmcmc$init_P-state$init_P,'\n')
        cat('Iv check: ',init4pmcmc$init_Iv-state$init_Iv,'\n')
        cat('init_EIR check: ',state$init_EIR-init4pmcmc$init_EIR,'\n')
        cat('prev check: ',state$prev-init4pmcmc$prev,'\n')

        print('Helpful values for reference:')
        cat('init_EIR: ',state$init_EIR,'\n')
        cat('Equilibrium prev: ',state$prev,'\n')
        cat('Seasonal equilibrium prev: ',init4pmcmc$prev,'\n')
        saveRDS(append(mpl,init4pmcmc),'seasonal_equil_values.RDS')
        saveRDS(state,'original_equil_values.RDS')

        # mpl['init_EIR'] <- NULL
        # View(init4pmcmc)
        # View(state_use)
      }
      return(append(mpl,init4pmcmc)) #Append all parameters from model parameter list for stochastic model
    }
    else{
      return(state)
    }
  }
}
#------------------------------------------------
#' Function that returns values from optional seasonal deterministic model
#'
#' \code{check_seasonality} Calculates time series of the prefix seasonal deterministic
#'    model given the posterior distribution of pMCMC parameters
#'
#' @param theta Posterior distribution of pMCMC parameters. Default = NULL
#' @param mpl_pf Model parameter list. Default = NULL
#' @param season_model Seasonality model to be used for the optional deterministic
#'    model. Default = NULL
#'
#' @export
check_seasonality <- function(theta,mpl_pf,season_model){
  init_EIR <- exp(theta[["log_init_EIR"]]) ## Exponentiate EIR since MCMC samples on the log scale for EIR
  EIR_vol <- theta[["volatility"]]
  mpl <- append(mpl_pf,list(volatility = EIR_vol)) ## Add MCMC parameters to model parameter list

  ## Run equilibrium function
  state <- equilibrium_init_create_stripped(age_vector = mpl$init_age,
                                            init_EIR = init_EIR,
                                            ft = mpl$prop_treated,
                                            model_param_list = mpl,
                                            het_brackets = mpl$het_brackets,
                                            state_check = mpl$state_check)
  # print(state)
  ##run seasonality model first if seasonality_on == TRUE
  state_use <- state[names(state) %in% coef(season_model)$name]

  # create model with initial values
  mod <- season_model$new(user = state_use, use_dde = TRUE)

  # tt <- c(0, preyears*365+as.integer(difftime(mpl$start_stoch,mpl$time_origin,units="days")))
  tt <- seq(0, mpl$preyears*365+as.integer(difftime(mpl$start_stoch,mpl$time_origin,units="days")),length.out=100)

  # run seasonality model
  mod_run <- mod$run(tt, verbose=FALSE,step_size_max=9)

  # shape output
  out <- mod$transform_variables(mod_run)
  out.df <- data.frame(t=out$t,
                       prev05 = out$prev,
                       prev_all = out$prev_all,
                       clininc_05 = out$inc05,
                       clininc_all = out$inc)
  return(out.df)
}

#------------------------------------------------
#' Converts log-odds to probability (or prevalence)
#'
#' \code{get_prev_from_log_odds} Converts log-odds to probability (or prevalence)
#'
#' @param log_odds A log-odds value.
#'
#' @export
get_prev_from_log_odds<-function(log_odds){
  return(exp(log_odds)/(1+exp(log_odds)))
}
#------------------------------------------------
#' Converts probability (prevalence) to log-odds
#'
#' \code{get_odds_from_prev} Converts probability (prevalence) to log-odds
#'
#' @param prev A probability or prevalence value. Must be between 0 and 1.
#'
#' @export
get_odds_from_prev<-function(prev){
  return(prev/(1-prev))
}
#' Return an initial EIR based on a user-given target prevalence in <5 yead old children
#'
#' \code{get_init_EIR} Return an initial EIR based on a user-given target prevalence in <5 yead old children
#'
#' @param par An initial EIR value from optim function
#' @param mpl Model parameter list
#'
#' @export
get_init_EIR <- function(par,mpl){
  init_EIR <- par[1]
  equil <- mamasante::equilibrium_init_create_stripped(age_vector = mpl$init_age,
                                                    ft = mpl$prop_treated,
                                                    het_brackets = mpl$het_brackets,
                                                    init_EIR = init_EIR,
                                                    model_param_list = mpl)
  prev_age <- ifelse(mpl$target_prev_group =='u5','prev05','prev2.10')


  return((equil[[prev_age]] - mpl$target_prev)^2)
}
#' Return an initial EIR based on a user-given target prevalence in <5 yead old children
#'
#' \code{get_init_EIR_u5} Return an initial EIR based on a user-given target prevalence in <5 yead old children
#'
#' @param par An initial EIR value from optim function
#' @param mpl Model parameter list
#'
#' @export
get_init_EIR_u5 <- function(par,mpl){
  init_EIR <- par[1]
  equil <- mamasante::equilibrium_init_create_stripped(age_vector = mpl$init_age,
                                                    ft = mpl$prop_treated,
                                                    het_brackets = mpl$het_brackets,
                                                    init_EIR = init_EIR,
                                                    model_param_list = mpl)

  return((equil$prev05 - mpl$target_prev)^2)
}

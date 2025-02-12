#------------------------------------------------
#' Equilibrium initialisation list creation
#'
#' \code{equilibrium_init_create_stripped} creates an equilibrium initialisation state to be
#' used within later model runs
#' @param age_vector Vector of age brackets.
#' @param het_brackets Integer number of biting heteogenity compartments.
#' @param country String for country of interest. If NULL the seasonal parameters
#' will attempt to be loaded using just the admin unit, however if there is ambiguity
#' in the admin unit an error will be thrown. If both NULL then no seasonality is
#' assumed. Default = NULL.
#' @param admin_unit String for admin unit with country for loading seasonal
#' parameters. If country is NULL, the admin unit will attempt to be located,however
#' if there is ambiguity in the admin unit an error will be thrown. If both country
#' and admin_unit are NULL then no seasonality is assumed. Default = NULL.
#' @param ft Numeric for the frequency of people seeking treatment.
#' @param init_EIR Numeric for desired annual EIR.
#' @param model_param_list List of epidemiological parameters created by
#' @param state_check If state_check = TRUE, returns expected deriv values which should equal 0 and sets stochastic model to have EIR constant at init_EIR
#'                    If state_check = TRUE and seasonality_on = 1, then the deterministic seasonal model is still run, but theta2 is forced to 1, forcing a constant seasonality profile
#'                    If state_check = FALSE, no values are printed
#'
#' @importFrom stringi stri_trans_general
#' @importFrom statmod gauss.quad.prob
#'
#'
#' @export


equilibrium_init_create_stripped <- function(age_vector, het_brackets,
                                             country = NULL, admin_unit = NULL, ft,
                                             init_EIR, model_param_list, state_check = 0)
{

  # mpl is shorter :)
  mpl <- model_param_list

  ## Check Parameters
  if(!is.numeric(age_vector)) stop("age_vector provided is not numeric")
  if(!is.numeric(het_brackets)) stop("het_brackets provided is not numeric")
  if(!(is.null(country) | is.character(country))) stop("country specified is not character string")
  if(!(is.null(admin_unit) | is.character(admin_unit))) stop("admin_unit specified is not character string")
  if(!is.numeric(ft)) stop("ft provided is not numeric")
  if(!is.numeric(init_EIR)) stop("EIR provided is not numeric")

  ## Handle parameters
  # database for admin units is all in Latin-ASCII for CRAN reasons so must
  # encode parameters accordingly
  if(!is.null(country)) country <- stringi::stri_trans_general(country,"Latin-ASCII")
  if(!is.null(admin_unit)) admin_unit <- stringi::stri_trans_general(admin_unit, "Latin-ASCII")

  ## population demographics
  age <- age_vector * mpl$DY
  na <- as.integer(length(age))  # number of age groups
  nh <- as.integer(het_brackets)  # number of heterogeneity groups
  h <- statmod::gauss.quad.prob(nh, dist = "normal")

  age_rate <- age_width <- age_mid_point <- den <- c()
  for (i in 1:(na-1))
  {
    age_width[i] <- age[i+1] - age[i]
    age_rate[i] <- 1/(age[i + 1] - age[i])  # vector of rates at which people leave each age group (1/age group width)
    age_mid_point[i] <- 0.5 * (age[i] + age[i + 1])  # set age group vector to the midpoint of the group

  }
  age_rate[na] = 0


  den <- 1/(1 + age_rate[1]/mpl$eta)
  for (i in 1:(na-1))
  {
    den[i+1] <- age_rate[i] * den[i]/(age_rate[i+1] + mpl$eta)  # proportion in each age_vector group
  }

  age59 <- which(age_vector * 12 > 59)[1] - 1  # index of age vector before age is >59 months
  age05 <- which(age_vector > 5)[1] - 1  # index of age vector before age is 5 years
  age02 <- which(age_vector > 2)[1] - 1  # index of age vector before age is 5 years
  age10 <- which(age_vector > 10)[1] - 1  # index of age vector before age is 5 years

  age0 <- which(age_vector > mpl$age_min)[1] - 1 #Allow a user-defined age group
  age1 <- which(age_vector > mpl$age_max)[1] - 1

  ## force of infection
  foi_age <- c()
  # for (i in 1:na){
  #    foi_age[i] <- 1 - (mpl$rho * exp(-age[i]/mpl$a0))  #force of infection for each age group
  #  }

  for (i in 1:(na-1))
  {
    foi_age[i] <- 1 - (mpl$rho * exp(-(age[i]+age[i+1])/2/mpl$a0))  #force of infection for each age group
  }
  foi_age[na]<-1 - (mpl$rho * exp(-(age[i])/mpl$a0))
  fden <- foi_age * den
  omega <- sum(fden)  #normalising constant
  ## heterogeneity
  het_x <- h$nodes
  het_wt <- h$weights
  den_het <- outer(den, het_wt)
  # rel_foi <- exp(-mpl$sigma2/2 + sqrt(mpl$sigma2) * het_x)/sum(het_wt * exp(-mpl$sigma2/2 + sqrt(mpl$sigma2) * het_x))
  rel_foi <- exp(-mpl$sigma2/2 + sqrt(mpl$sigma2) * het_x)

  ## EIR
  EIRY_eq <- init_EIR  # initial annual EIR
  EIRd_eq <- EIRY_eq/mpl$DY
  EIR_eq <- outer(foi_age, rel_foi) * EIRd_eq

  ## Immunity and FOI
  x_I <- den[1]/mpl$eta
  for (i in 2:na)
  {
    x_I[i] <- den[i]/(den[i - 1] * age_rate[i - 1])  #temporary variables
  }
  # fd <- 1 - (1 - mpl$fD0)/(1 + (age/mpl$aD)^mpl$gammaD)
  fd <- c()
  for (i in 1:(na-1))
  {
    fd[i] <- 1 - (1 - mpl$fD0)/(1 + ((age[i]+age[i+1])/2/mpl$aD)^mpl$gammaD)
  }
  fd[na]<-1 - (1 - mpl$fD0)/(1 + (age[na]/mpl$aD)^mpl$gammaD)

  # maternal immunity begins at a level proportional to the clinical
  # immunity of a 20 year old, this code finds that level
  age20i <- rep(0, na)
  for (i in 2:na)
  {
    age20i[i] <- ifelse(age[i] >= (20 * mpl$DY) & age[i - 1] < (20 * mpl$DY),
                        i, age20i[i - 1])
  }
  age20u <- as.integer(age20i[na])
  age20l <- as.integer(age20u - 1)
  age_20_factor <- (20 * mpl$DY - age[age20l] - 0.5 * age_width[age20l]) *
    2/(age_width[age20l] + age_width[age20u])

  # finding initial values for all immunity states
  IB_eq <- matrix(0, na, nh)
  FOI_eq <- matrix(0, na, nh)
  ID_eq <- matrix(0, na, nh)
  ICA_eq <- matrix(0, na, nh)
  IC_20 <- matrix(0, 1, nh)
  ICM_age <- matrix(0, na,1)
  cA_eq <- matrix(0, na, nh)
  FOIvij_eq <- matrix(0, na, nh)
  p_det_eq <- matrix(0, na, nh)
  for (j in 1:nh)
  {
    for (i in 1:na)
    {
      IB_eq[i, j] <- (ifelse(i == 1, 0, IB_eq[i - 1, j]) +
                        EIR_eq[i,j]/(EIR_eq[i, j] * mpl$uB + 1) * x_I[i])/(1 + x_I[i]/mpl$dB)
      FOI_eq[i, j] <- EIR_eq[i, j] * ifelse(IB_eq[i, j] == 0, mpl$b0,
                                            mpl$b0 * ((1 - mpl$b1)/(1 + (IB_eq[i, j]/mpl$IB0)^mpl$kB) + mpl$b1))
      ID_eq[i, j] <- (ifelse(i == 1, 0, ID_eq[i - 1, j]) +
                        FOI_eq[i, j]/(FOI_eq[i, j] * mpl$uD + 1) * x_I[i])/(1 + x_I[i]/mpl$dID)
      ICA_eq[i, j] <- (ifelse(i == 1, 0, ICA_eq[i - 1, j]) +
                         FOI_eq[i,j]/(FOI_eq[i, j] * mpl$uCA + 1) * x_I[i])/(1 + x_I[i]/mpl$dCA)
      p_det_eq[i, j] <- mpl$d1 + (1 - mpl$d1)/(1 + fd[i] * (ID_eq[i, j]/mpl$ID0)^mpl$kD)
      cA_eq[i, j] <- mpl$cU + (mpl$cD - mpl$cU) * p_det_eq[i, j]^mpl$gamma1
    }
  }
  # needs to be calculated after because it references ICA


  for (j in 1:nh)
    IC_20[j] <- mpl$PM * (ICA_eq[age20l, j] + age_20_factor *
                            (ICA_eq[age20u, j] - ICA_eq[age20l, j]))
  for (i in 1:(na-1))
    ICM_age[i]<- mpl$dCM/(age[i+1]-age[i])*(exp(-age[i]/mpl$dCM)-exp(-age[i+1]/mpl$dCM))
  ICM_age[na]<-0
  ICM_eq<-ICM_age%*%IC_20
  IC_eq <- ICA_eq + ICM_eq
  phi_eq <- mpl$phi0 * ((1 - mpl$phi1)/(1 + (IC_eq/mpl$IC0)^mpl$kC) + mpl$phi1)

  # human states
  gamma <- mpl$eta + c(age_rate[1:(na - 1)], 0)
  delta <- c(mpl$eta, age_rate[1:(na - 1)])

  betaT <- matrix(rep(mpl$rT + gamma, rep(nh, na)), ncol = nh, byrow = TRUE)
  betaD <- matrix(rep(mpl$rD + gamma, rep(nh, na)), ncol = nh, byrow = TRUE)
  betaP <- matrix(rep(mpl$rP + gamma, rep(nh, na)), ncol = nh, byrow = TRUE)

  aT <- FOI_eq * phi_eq * ft/betaT
  aD <- FOI_eq * phi_eq * (1 - ft)/betaD
  aP <- mpl$rT * aT/betaP

  Z_eq <- array(dim = c(na, nh, 4))
  Z_eq[1, , 1] <- den_het[1, ]/(1 + aT[1, ] + aD[1, ] + aP[1, ])
  Z_eq[1, , 2] <- aT[1, ] * Z_eq[1, , 1]
  Z_eq[1, , 3] <- aD[1, ] * Z_eq[1, , 1]
  Z_eq[1, , 4] <- aP[1, ] * Z_eq[1, , 1]

  for (j in 1:nh)
  {
    for (i in 2:na)
    {
      Z_eq[i, j, 1] <- (den_het[i, j] - delta[i] * (Z_eq[i - 1, j, 2]/betaT[i, j] +
                                                      Z_eq[i - 1, j, 3]/betaD[i, j] +
                                                      (mpl$rT *  Z_eq[i - 1, j, 2]/betaT[i, j]
                                                       + Z_eq[i - 1, j, 4])/betaP[i, j]))/(1 + aT[i, j] + aD[i, j] + aP[i, j])
      Z_eq[i, j, 2] <- aT[i, j] * Z_eq[i, j, 1] + delta[i] * Z_eq[i -
                                                                    1, j, 2]/betaT[i, j]
      Z_eq[i, j, 3] <- aD[i, j] * Z_eq[i, j, 1] + delta[i] * Z_eq[i -
                                                                    1, j, 3]/betaD[i, j]
      Z_eq[i, j, 4] <- aP[i, j] * Z_eq[i, j, 1] + delta[i] * (mpl$rT *
                                                                Z_eq[i - 1, j, 2]/betaT[i, j] + Z_eq[i - 1, j, 4])/betaP[i,j]

    }
  }

  Y_eq <- matrix(Z_eq[, , 1], nrow = na, ncol=nh)
  T_eq <- matrix(Z_eq[, , 2], nrow = na, ncol=nh)
  D_eq <- matrix(Z_eq[, , 3], nrow = na, ncol=nh)
  P_eq <- matrix(Z_eq[, , 4], nrow = na, ncol=nh)

  betaS <- apply(FOI_eq, MARGIN = 2, FUN = function(x, y)
  {
    x + y
  }, y = gamma)
  betaA <- apply(FOI_eq * phi_eq + mpl$rA, MARGIN = 2, FUN = function(x, y)
  {
    x + y
  }, y = gamma)
  betaU <- apply(FOI_eq + mpl$rU, MARGIN = 2, FUN = function(x, y)
  {
    x + y
  }, y = gamma)

  A_eq <- matrix(ncol = nh, nrow = na)
  U_eq <- matrix(ncol = nh, nrow = na)
  S_eq <- matrix(ncol = nh, nrow = na)
  prev_eq <- matrix(ncol = nh, nrow = na)

  for (i in 1:na)
  {
    for (j in 1:nh)
    {
      A_eq[i, j] <- (delta[i] * ifelse(i == 1, 0, A_eq[i - 1, j]) +
                       FOI_eq[i, j] * (1 - phi_eq[i, j]) * Y_eq[i, j] +
                       mpl$rD * D_eq[i,j])/(betaA[i, j] + FOI_eq[i, j] * (1 - phi_eq[i, j]))
      U_eq[i, j] <- (mpl$rA * A_eq[i, j] + delta[i] * ifelse(i == 1,
                                                             0, U_eq[i - 1, j]))/betaU[i, j]
      S_eq[i, j] <- Y_eq[i, j] - A_eq[i, j] - U_eq[i, j]
      FOIvij_eq[i, j] <- foi_age[i] * mpl$av0 * (mpl$cT * T_eq[i, j] + mpl$cD *
                                                 D_eq[i, j] + cA_eq[i, j] * A_eq[i, j] + mpl$cU * U_eq[i, j]) * rel_foi[j]/omega
      prev_eq[i,j] <- T_eq[i,j] + D_eq[i,j]  + A_eq[i,j]*p_det_eq[i,j]
      # print(FOIvij_eq[i,j])
    }
  }
  prev <- sum(prev_eq[1:age59,])/sum(den[1:age59])
  prev2.10 <- sum(prev_eq[age02:age10,])/sum(den[age02:age10])

  age_flex_length <- age1 - age0 + 1

  prev_flex <- sum(prev_eq[age0:age1,])/sum(den[age0:age1])

  # print(FOIvij_eq)
  # mosquito states
  FOIv_eq <- sum(FOIvij_eq)
  # print(FOIv_eq)
  Iv_eq <- FOIv_eq * mpl$Surv0/(FOIv_eq + mpl$mu0)
  Sv_eq <- mpl$mu0 * Iv_eq/(FOIv_eq * mpl$Surv0)
  Ev_eq <- 1 - Sv_eq - Iv_eq

  # mosquito density needed to give this EIR
  mv0 <- omega * EIRd_eq/(Iv_eq * mpl$av0)
  betaa_eq <- (FOIv_eq + mpl$mu0)*Sv_eq*mv0

  # larval states
  K0 <- 2 * mv0 * mpl$dLL * mpl$mu0 * (1 + mpl$dPL * mpl$muPL) * mpl$gammaL * (mpl$lambda + 1)/(mpl$lambda/(mpl$muLL *
  mpl$dEL) - 1/(mpl$muLL * mpl$dLL) - 1)
  PL_eq <- 2 * mpl$dPL * mpl$mu0 * mv0
  LL_eq <- mpl$dLL * (mpl$muPL + 1/mpl$dPL) * PL_eq
  EL_eq <- (LL_eq/mpl$dLL + mpl$muLL* LL_eq * (1 + mpl$gammaL * LL_eq/K0))/(1/mpl$dEL - mpl$muLL * mpl$gammaL * LL_eq/K0)

  IB_eq = array(IB_eq, c(na, nh))
  ID_eq = array(ID_eq, c(na, nh))
  ICA_eq = array(ICA_eq, c(na, nh))
  ICM_eq = array(ICM_eq, c(na, nh))
  #ssa0 <- ssa1 <- ssa2 <- ssa3 <- ssb1 <- ssb2 <- ssb3 <- theta_c <- 0
  # better het bounds for equilbirum initialisation in individual model
  zetas <- rlnorm(n = 1e5,meanlog = -mpl$sigma2/2, sdlog = sqrt(mpl$sigma2))
  while(sum(zetas>100)>0){
    zetas[zetas>100] <- rlnorm(n = sum(zetas>100),meanlog = -mpl$sigma2/2, sdlog = sqrt(mpl$sigma2))
  }

  wt_cuts <- round(cumsum(het_wt)*1e5)
  zeros <- which(wt_cuts==0)
  wt_cuts[zeros] <- 1:length(zeros)
  larges <- which(wt_cuts==1e5)
  wt_cuts[larges] <- (1e5 - (length(larges)-1)):1e5
  wt_cuts <- c(0,wt_cuts)
  het_bounds <- sort(zetas)[wt_cuts]
  het_bounds[length(het_bounds)] <- (mpl$max_age/365)+1
  b <- mpl$b0 * ((1-mpl$b1)/(1+(IB_eq/mpl$IB0)^mpl$kB)+mpl$b1)
  phi <- mpl$phi0*((1-mpl$phi1)/(1+((ICM_eq+ICA_eq)/mpl$IC0)^mpl$kC) + mpl$phi1)
  phi <- array(phi, c(na, nh))
  clin_inc <- phi*FOI_eq*Y_eq
  clin_inc <- array(clin_inc, c(na, nh))

  inc <- sum(clin_inc)
  inc05 <- sum(clin_inc[1:age05,])/sum(den[1:age59])
  inc_flex <- sum(clin_inc[age0:age1,])/sum(den[age0:age1])


  ## collate init
  res <- list(init_S = S_eq, init_T = T_eq, init_D = D_eq, init_A = A_eq, init_U = U_eq,
              init_P = P_eq, init_IB = IB_eq, init_ID = ID_eq, init_ICA = ICA_eq,
              ICM_age = as.vector(ICM_age),
              init_Iv = Iv_eq, init_Sv = Sv_eq,init_Ev = Ev_eq, init_PL = PL_eq, init_LL = LL_eq, init_EL = EL_eq,
              age_rate = age_rate, het_wt = het_wt,
              foi_age = foi_age, rel_foi = rel_foi,
              na = na, nh = nh, x_I = x_I,
              omega = omega, mv0 = mv0, betaa_eq = betaa_eq,
              FOIv_eq = FOIv_eq,
              FOI_eq = FOI_eq, EIR_eq = EIR_eq, init_EIR = init_EIR,
              den = den, age59 = age59, age05 = age05,age10=age10,
              age0 = age0, age1 = age1,
              age_flex_length = age_flex_length,
              pi = pi,
              prev05 = prev,inc = inc, inc05=inc05,
              prev_flex = prev_flex, inc_flex = inc_flex,
              prev2.10 = prev2.10,
              age = age_vector*mpl$DY, ft = ft,
              age20l = age20l, age20u = age20u, age_20_factor = age_20_factor
              ##extras for checking:
              # ,ssa0 = mpl$ssa0, ssa1 = mpl$ssa1,ssa2 = mpl$ssa2, ssa3 = mpl$ssa3, ssb1 = mpl$ssb1, ssb2 = mpl$ssb2, ssb3 = mpl$ssb3,theta_c = mpl$theta_c,
              # init_Y = Y_eq, init_ICM = ICM_eq,IC_20=IC_20,age_width = age_width,
              # het_x = het_x,K0 = K0,cA_eq = cA_eq,
              # betaS = betaS, betaA = betaA, betaU = betaU,FOIvij_eq=FOIvij_eq,
              # age_mid_point = age_mid_point, het_bounds = het_bounds,
              # p_det_eq = p_det_eq, b_eq = b, phi_eq = phi


              )

  ##Check that equilibrium solution produces an equilibrium for
  ##the desired model
  if(state_check==1){
    if(mpl$seasonality_on==0){
      H <- sum(S_eq) + sum(T_eq) + sum(D_eq) + sum(A_eq) + sum(U_eq) + sum(P_eq)

      deriv_S1 <- -FOI_eq[1,]*S_eq[1,] + mpl$rP*P_eq[1,] + mpl$rU*U_eq[1,] +
        mpl$eta*H*het_wt - (mpl$eta+age_rate[1])*S_eq[1,]
      deriv_S2 <- -FOI_eq[2,]*S_eq[2,] + mpl$rP*P_eq[2,] + mpl$rU*U_eq[2,] - (mpl$eta+age_rate[2])*S_eq[2,] + age_rate[1]*S_eq[1,]
      # deriv(S[2:na, 1:nh]) <- -FOI[i,j]*S[i,j] + rP*P[i,j] + rU*U[i,j] -
      #   (eta+age_rate[i])*S[i,j] + age_rate[i-1]*S[i-1,j]
      phi <- mpl$phi0*((1-mpl$phi1)/(1+((ICM_eq+ICA_eq)/mpl$IC0)^mpl$kC) + mpl$phi1)
      phi <- array(phi, c(na, nh))

      clin_inc <- phi*FOI_eq*Y_eq
      clin_inc <- array(clin_inc, c(na, nh))

      deriv_D1 <- (1-ft)*clin_inc[1,] - mpl$rD*D_eq[1,] -
        (mpl$eta+age_rate[1])*D_eq[1,]
      # deriv(D[1, 1:nh]) <- (1-ft)*clin_inc[i,j] - rD*D[i,j] -
      #   (eta+age_rate[i])*D[i,j]
      deriv_A1 <- (1-phi[1,])*FOI_eq[1,]*Y_eq[1,] - FOI_eq[1,]*A_eq[1,] +
        mpl$rD*D_eq[1,] - mpl$rA*A_eq[1,] - (mpl$eta+age_rate[1])*A_eq[1,]
      # deriv(A[2:na, 1:nh]) <- (1-phi[i,j])*FOI[i,j]*Y[i,j] - FOI[i,j]*A[i,j] +
      #   rD*D[i,j] - rA*A[i,j] - (eta+age_rate[i])*A[i,j] + age_rate[i-1]*A[i-1,j]
      deriv_ICA2 <- FOI_eq[2,]/(FOI_eq[2,] * mpl$uCA + 1) - 1/mpl$dCA*ICA_eq[2,] - (ICA_eq[2,]-ICA_eq[1,])/x_I[2]
      # deriv(ICA[2:na, 1:nh]) <- FOI[i,j]/(FOI[i,j] * uCA + 1) - 1/dCA*ICA[i,j] - (ICA[i,j]-ICA[i-1,j])/x_I[i]
      deriv_IB3 <- EIR_eq[3,]/(EIR_eq[3,]* mpl$uB + 1) - IB_eq[3,]/mpl$dB - (IB_eq[3,]-IB_eq[2,])/x_I[3]
      # deriv(IB[2:na, 1:nh]) <- EIR[i,j]/(EIR[i,j]* uB + 1) - IB[i,j]/dB - (IB[i,j]-IB[i-1,j])/x_I[i]
      deriv_ID4 <- FOI_eq[4,]/(FOI_eq[4,]*mpl$uD + 1) - ID_eq[4,]/mpl$dID - (ID_eq[4,]-ID_eq[3,])/x_I[4]
      # deriv(ID[2:na, 1:nh]) <- FOI[i,j]/(FOI[i,j]*uD + 1) - ID[i,j]/dID - (ID[i,j]-ID[i-1,j])/x_I[i]

      ##FOI delay derivs
      # deriv(FOI[,,1]) <- (lag_rates/dE)*FOI_lag[i,j] - (lag_rates/dE)*FOI[i,j,1]
      # deriv(FOI[,,2:lag_rates]) <- (lag_rates/dE)*FOI[i,j,k-1] - (lag_rates/dE)*FOI[i,j,k]
      b <- mpl$b0 * ((1-mpl$b1)/(1+(IB_eq/mpl$IB0)^mpl$kB)+mpl$b1)

      FOI_lag<- matrix(ncol = nh, nrow = na)
      for(i in 1:na){
        for(j in 1:nh){
          FOI_lag[i,j] <- EIR_eq[i,j] * (if(IB_eq[i,j]==0) mpl$b0 else b[i,j])
        }
      }
      init_FOI_delay <- FOI_eq
      deriv_FOI1 <- (mpl$lag_rates/mpl$dE)*FOI_lag - (mpl$lag_rates/mpl$dE)*init_FOI_delay
      # deriv(FOI[,,2:lag_rates]) <- (lag_rates/dE)*FOI[i,j,k-1] - (lag_rates/dE)*FOI[i,j,k]

      cat('S[1,] derivative = ',deriv_S1,'\n')
      cat('S[2,] derivative = ',deriv_S2,'\n')
      cat('D[1,] derivative = ',deriv_D1,'\n')
      cat('A[1,] derivative = ',deriv_A1,'\n')
      cat('ICA[2,] derivative = ',deriv_ICA2,'\n')
      cat('IB[3,] derivative = ',deriv_IB3,'\n')
      cat('ID[4,] derivative = ',deriv_ID4,'\n')
      cat('FOI[,,1] derivative = ',deriv_FOI1,'\n')
    }
    else{
      H <- sum(S_eq) + sum(T_eq) + sum(D_eq) + sum(A_eq) + sum(U_eq) + sum(P_eq)

      deriv_S1 <- -FOI_eq[1,]*S_eq[1,] + mpl$rP*P_eq[1,] + mpl$rU*U_eq[1,] +
        mpl$eta*H*het_wt - (mpl$eta+age_rate[1])*S_eq[1,]
      deriv_S2 <- -FOI_eq[2,]*S_eq[2,] + mpl$rP*P_eq[2,] + mpl$rU*U_eq[2,] - (mpl$eta+age_rate[2])*S_eq[2,] + age_rate[1]*S_eq[1,]
      # deriv(S[2:na, 1:nh]) <- -FOI[i,j]*S[i,j] + rP*P[i,j] + rU*U[i,j] -
      #   (eta+age_rate[i])*S[i,j] + age_rate[i-1]*S[i-1,j]
      phi <- mpl$phi0*((1-mpl$phi1)/(1+((ICM_eq+ICA_eq)/mpl$IC0)^mpl$kC) + mpl$phi1)
      phi <- array(phi, c(na, nh))

      clin_inc <- phi*FOI_eq*Y_eq
      clin_inc <- array(clin_inc, c(na, nh))

      deriv_D1 <- (1-ft)*clin_inc[1,] - mpl$rD*D_eq[1,] -
        (mpl$eta+age_rate[1])*D_eq[1,]
      # deriv(D[1, 1:nh]) <- (1-ft)*clin_inc[i,j] - rD*D[i,j] -
      #   (eta+age_rate[i])*D[i,j]
      deriv_A1 <- (1-phi[1,])*FOI_eq[1,]*Y_eq[1,] - FOI_eq[1,]*A_eq[1,] +
        mpl$rD*D_eq[1,] - mpl$rA*A_eq[1,] - (mpl$eta+age_rate[1])*A_eq[1,]
      # deriv(A[2:na, 1:nh]) <- (1-phi[i,j])*FOI[i,j]*Y[i,j] - FOI[i,j]*A[i,j] +
      #   rD*D[i,j] - rA*A[i,j] - (eta+age_rate[i])*A[i,j] + age_rate[i-1]*A[i-1,j]
      # deriv(IB[1, 1:nh]) <- EIR[i,j]/(EIR[i,j]* uB + 1) - IB[i,j]/dB - IB[i,j]/x_I[i]
      deriv_IB1 <- EIR_eq[1,]/(EIR_eq[1,]* mpl$uB + 1) - IB_eq[1,]/mpl$dB - IB_eq[1,]/x_I[1]
      deriv_ICA2 <- FOI_eq[2,]/(FOI_eq[2,] * mpl$uCA + 1) - 1/mpl$dCA*ICA_eq[2,] - (ICA_eq[2,]-ICA_eq[1,])/x_I[2]
      # deriv(ICA[2:na, 1:nh]) <- FOI[i,j]/(FOI[i,j] * uCA + 1) - 1/dCA*ICA[i,j] - (ICA[i,j]-ICA[i-1,j])/x_I[i]
      deriv_IB3 <- EIR_eq[3,]/(EIR_eq[3,]* mpl$uB + 1) - IB_eq[3,]/mpl$dB - (IB_eq[3,]-IB_eq[2,])/x_I[3]
      # deriv(IB[2:na, 1:nh]) <- EIR[i,j]/(EIR[i,j]* uB + 1) - IB[i,j]/dB - (IB[i,j]-IB[i-1,j])/x_I[i]
      deriv_ID4 <- FOI_eq[4,]/(FOI_eq[4,]*mpl$uD + 1) - ID_eq[4,]/mpl$dID - (ID_eq[4,]-ID_eq[3,])/x_I[4]
      # deriv(ID[2:na, 1:nh]) <- FOI[i,j]/(FOI[i,j]*uD + 1) - ID[i,j]/dID - (ID[i,j]-ID[i-1,j])/x_I[i]

      ##FOI delay derivs
      # deriv(FOI[,,1]) <- (lag_rates/dE)*FOI_lag[i,j] - (lag_rates/dE)*FOI[i,j,1]
      # deriv(FOI[,,2:lag_rates]) <- (lag_rates/dE)*FOI[i,j,k-1] - (lag_rates/dE)*FOI[i,j,k]
      # FOI[,] <- EIR[i,j] * (if(IB[i,j]==0) b0 else b[i,j])
      b <- mpl$b0 * ((1-mpl$b1)/(1+(IB_eq/mpl$IB0)^mpl$kB)+mpl$b1)
      FOI_eq11 <- EIR_eq[1,1] * (if(IB_eq[1,1]==0) b0 else b[1,1])
      FOI_eq23 <- EIR_eq[2,3] * (if(IB_eq[2,3]==0) b0 else b[2,3])

      ince <- FOIv_eq*Sv_eq*mv0
      betaa <- 0.5*PL_eq/mpl$dPL
      # betaa <- mpl$mu0 * mv0
      deriv_Sv <- -ince - mpl$mu0*Sv_eq*mv0 + betaa
      # deriv(Sv) <- -ince - mu*Sv + betaa
      surv <- exp(-mpl$mu0*mpl$delayMos)
      incv <- ince * surv

      deriv_Ev <- ince - incv - mpl$mu0*Ev_eq*mv0
      deriv_Iv <- incv - mpl$mu0*Iv_eq*mv0
      # deriv(Ev) <- ince - incv - mu*Ev
      # deriv(Iv) <- incv - mu*Iv

      K0_eq <- 2*mv0*mpl$dLL*mpl$mu0*(1+mpl$dPL*mpl$muPL)*mpl$gammaL*(mpl$lambda+1)/(mpl$lambda/(mpl$muLL*mpl$dEL)-1/(mpl$muLL*mpl$dLL)-1)

      deriv_PL <- LL_eq/mpl$dLL - mpl$muPL*PL_eq - PL_eq/mpl$dPL
      # deriv(PL) <- LL/dLL - muPL*PL - PL/dPL
      deriv_LL <- EL_eq/mpl$dEL - mpl$muLL*(1+mpl$gammaL*(EL_eq + LL_eq)/K0_eq)*LL_eq - LL_eq/mpl$dLL
      # deriv(LL) <- EL/dEL - muLL*(1+gammaL*(EL + LL)/KL)*LL - LL/dLL
      fv <- 1/( mpl$tau1 + mpl$tau2 ) # mosquito feeding rate (zbar from intervention param.)
      eov <- mpl$betaL/mpl$mu0*(exp(mpl$mu0/fv)-1)
      beta_larval <- eov*mpl$mu0*exp(-mpl$mu0/fv)/(1-exp(-mpl$mu0/fv)) # Number of eggs laid per day
      deriv_EL <- beta_larval*mv0-mpl$muEL*(1+(EL_eq+LL_eq)/K0_eq)*EL_eq - EL_eq/mpl$dEL
      # deriv(EL) <- beta_larval*mv-muEL*(1+(EL+LL)/KL)*EL - EL/dEL

      fv <- 1/( mpl$tau1 + mpl$tau2 )
      av <- mpl$Q0*fv
      EIR <- rel_foi[1] * foi_age[1] * av * Iv_eq*mv0/omega
      t <- c(0:10)
      theta2 <- (mpl$ssa0+mpl$ssa1*cos(2*pi*t/365)+mpl$ssa2*cos(2*2*pi*t/365)+mpl$ssa3*cos(3*2*pi*t/365)+mpl$ssb1*sin(2*pi*t/365)+mpl$ssb2*sin(2*2*pi*t/365)+ mpl$ssb3*sin(3*2*pi*t/365) ) /mpl$theta_c

      cat('S[1,] derivative = ',deriv_S1,'\n')
      cat('S[2,] derivative = ',deriv_S2,'\n')
      cat('D[1,] derivative = ',deriv_D1,'\n')
      cat('A[1,] derivative = ',deriv_A1,'\n')
      cat('ICA[2,] derivative = ',deriv_ICA2,'\n')
      cat('IB[1,] derivative = ',deriv_IB1,'\n')
      cat('IB[3,] derivative = ',deriv_IB3,'\n')
      cat('ID[4,] derivative = ',deriv_ID4,'\n')
      cat('Sv derivative = ',deriv_Sv,'\n')
      cat('Ev derivative = ',deriv_Ev,'\n')
      cat('Iv derivative = ',deriv_Iv,'\n')
      cat('Sv initial = ',Sv_eq,'\n')
      cat('Ev initial = ',Ev_eq,'\n')
      cat('Iv initial = ',Iv_eq,'\n')
      cat('EL derivative = ',deriv_EL,'\n')
      cat('LL derivative = ',deriv_LL,'\n')
      cat('PL derivative = ',deriv_PL,'\n')
      cat('FOI[1,] = ',FOI_eq11,'\n')
      cat('FOI[2,3] = ',FOI_eq23,'\n')
      cat('betaa = ',betaa,'\n')
      cat('betaa_eq = ',betaa_eq,'\n')
      cat('ince = ',ince,'\n')
      cat('FOIv = ',FOIv_eq,'\n')
      cat('mv = ',mv0,'\n')
      cat('K0 = ',K0,'\n')
      cat('K0_eq = ',K0_eq,'\n')
      cat('av = ',av,'\n')
      cat('av0 = ',mpl$av0,'\n')
      cat('EIR[1,1] = ',EIR,'\n')
      cat('EIR_eq[1,1] = ',EIR_eq[1,1],'\n')
      cat('theta2 = ',theta2,'\n')

    }
  }
  res <- append(res,mpl)

  return(res)
}

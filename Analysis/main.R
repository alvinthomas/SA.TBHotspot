###### Tuberculosis Deterministic Model ######
### Made for Infectious Disease Dynamics
### Based on "Heterogeneity in tuberculosis transmission and the
### role of geographic hotspots in propagating epidemics"
### Dowdy et. al
###
### Quarter 4, 2015
### Johns Hopkins Univeristy
### Alvin G. Thomas
### Created: May 6, 2015
### Last Modified: October 17, 2015

## Required Packages ##
# install.packages("ggplot2")
# install.packages("ggthemes")
library("ggplot2")
library("ggthemes")

## Support Scripts ##
source("~/Documents/R Workspace/SA.TBHotspot/R/DeterministicSIR.r")

##### Variable Definitions #####

# ### Parameter Definitions for Rio de Janerio
# ### Based on Dowdy et. al
# # Number of transmissions per active TB per year, community
# DEF.BETA0 = 3.71
# # Number of transmissions per active TB per year, hotspot
# DEF.BETA1 = 9.74
# # Relative rate of transmission, hotspot-to-hotspot vs. hotspot-to-community
# DEF.RT = 0.03
# # Proportion of total population residing in the hotspot
# DEF.Z1 = 0.06
# # Proportion of total population residing in the community
# DEF.Z0 = 1 - DEF.Z1
# # Rate of rapid progression after recent infection, HIV-positive, per year
# DEF.ZETA1 = 0.31
# # Rate of slow progression after remote infection, HIV-positive, per year
# DEF.NU1 = 0.08
# # TB mortality rate, HIV-negative, per year
# DEF.MUTB0 = 0.031
# # TB mortality rate, HIV-positive, per year
# DEF.MUTB1 = 0.074
# # TB detection/treatment rate, HIV-negative, per year
# DEF.RHO0 = 0.87
# # TB detection/treatment rate, HIV-positive, per year
# DEF.RHO1 = 1.74
# # HIV incidence, per year
# DEF.PHI = 0.00015
# # HIV mortality rate (non-TB), per year
# DEF.MU1 = 0.026
# # TB relapse rate, per year
# DEF.PSI = 0.0083
# # Relative infectivity of HIV/TB cases
# DEF.RI = 0.68
# # Partial immunity to reinfection if latently infected
# DEF.P = 0.56
# # Inverse Duration of `recent infection' phase
# DEF.ETA = 1/5
# # Rate of rapid progression during this phase, HIV-negative, per year
# DEF.ZETA0 = 0.03
# # Rate of slow progression of remote TB infection, HIV-negative, per year
# DEF.NU0 = 0.0005
# # Inverse Life expectancy
# DEF.MU0 = 1/73
# # Relative HIV incidence in hotspot vs. community
# DEF.RH = 2.13
# # One millenium
# MAX.TIME = 1000

### Parameter Definitions for South Africa, Gauteng Province, Carletonville
# Number of transmissions per active TB per year, community
DEF.BETA0 = 1.77
# Number of transmissions per active TB per year, hotspot
DEF.BETA1 = 2.427
# Relative rate of transmission, hotspot-to-hotspot vs. hotspot-to-community
DEF.RT = 0.03
# Proportion of total population residing in the hotspot
DEF.Z1 = 0.0015
# Proportion of total population residing in the community
DEF.Z0 = 1 - DEF.Z1
# Rate of rapid progression after recent infection, HIV-positive, per year
DEF.ZETA1 = 0.31
# Rate of slow progression after remote infection, HIV-positive, per year
DEF.NU1 = 0.08
# TB mortality rate, HIV-negative, per year
DEF.MUTB0 = 48/100000
# TB mortality rate, HIV-positive, per year
DEF.MUTB1 = 121/100000
# TB detection/treatment rate, HIV-negative, per year
DEF.RHO0 = 0.53
# TB detection/treatment rate, HIV-positive, per year
DEF.RHO1 = 1.06
# HIV incidence, per year
DEF.PHI = 0.01
# HIV mortality rate (non-TB), per year
DEF.MU1 = 0.007
# TB relapse rate, per year
DEF.PSI = 0.0083
# Relative infectivity of HIV/TB cases
DEF.RI = 0.68
# Partial immunity to reinfection if latently infected
DEF.P = 0.56
# Inverse Duration of `recent infection' phase
DEF.ETA = 1/5
# Rate of rapid progression during this phase, HIV-negative, per year
DEF.ZETA0 = 0.03
# Rate of slow progression of remote TB infection, HIV-negative, per year
DEF.NU0 = 0.0005
# Inverse Life expectancy
DEF.MU0 = 1/56
# Relative HIV incidence in hotspot vs. community
DEF.RH = 2.65
# One millenium
MAX.TIME = 1000

dx.dt.SL2ARHOT <- function(t, y, param) {
  # There are five primary compartments
  # Each compartment is subdivided by HIV and Hotspot
  # h0=HIV uninfected
  # h1=HIV infected
  # g0=general population
  # g1=hotspot population

  if (t > 500){
    paramList <- param
  } else {
    paramList <- c(beta0=DEF.BETA0,
                   beta1=DEF.BETA1,
                   rt=DEF.RT,
                   z1=DEF.Z1,
                   z0=DEF.Z0,
                   zeta1=DEF.ZETA1,
                   nu1=DEF.NU1,
                   mutb0=DEF.MUTB0,
                   mutb1=DEF.MUTB1,
                   rho0=DEF.RHO0,
                   rho1=DEF.RHO1,
                   phi=DEF.PHI,
                   mu1=DEF.MU1,
                   psi=DEF.PSI,
                   ri=DEF.RI,
                   p=DEF.P,
                   eta=DEF.ETA,
                   zeta0=DEF.ZETA0,
                   nu0=DEF.NU0,
                   mu0=DEF.MU0,
                   rh=DEF.RH)
  }

  # MORTALITY
  mort_g0 <- as.numeric(paramList["mu0"] * (y["S_h0_g0"] + y["L1_h0_g0"] + y["L2_h0_g0"] +
                                              y["A_h0_g0"] + y["R_h0_g0"] + y["S_h1_g0"] +
                                              y["L1_h1_g0"] + y["L2_h1_g0"] +
                                              y["A_h1_g0"] + y["R_h1_g0"]) +
                          paramList["mu1"] * (y["S_h1_g0"] +
                                                y["L1_h1_g0"] + y["L2_h1_g0"] +
                                                y["A_h1_g0"] + y["R_h1_g0"]) +
                          paramList["mutb0"] * y["A_h0_g0"] +
                          paramList["mutb1"] * y["A_h1_g0"])

  mort_g1 <- as.numeric(paramList["mu0"] * (y["S_h0_g1"] + y["L1_h0_g1"] + y["L2_h0_g1"] +
                                              y["A_h0_g1"] + y["R_h0_g1"] + y["S_h1_g1"] +
                                              y["L1_h1_g1"] + y["L2_h1_g1"] +
                                              y["A_h1_g1"] + y["R_h1_g1"]) +
                          paramList["mu1"] * (y["S_h1_g1"] +
                                                y["L1_h1_g1"] + y["L2_h1_g1"] +
                                                y["A_h1_g1"] + y["R_h1_g1"]) +
                          paramList["mutb0"] * y["A_h0_g1"] +
                          paramList["mutb1"] * y["A_h1_g1"])

  # FORCE OF INFECTION
  lambda_g0 <- as.numeric((paramList["beta0"] / paramList["z0"]) *
                            (y["A_h0_g0"] + paramList["ri"] * y["A_h1_g0"]) +
                            (paramList["beta1"] / paramList["z1"]) * paramList["rt"] *
                            (y["A_h0_g1"] + paramList["ri"] * y["A_h1_g1"]))
  lambda_g1 <- as.numeric((paramList["beta1"] / paramList["z1"]) *
                            (y["A_h0_g1"] + paramList["ri"] * y["A_h1_g1"]) +
                            (paramList["beta0"] / paramList["z0"]) * paramList["rt"] *
                            (y["A_h0_g0"] + paramList["ri"] * y["A_h1_g0"]))

  # HIV
  # Corrected from Dowdy et al.
  hiv_h0_g0 <- as.numeric(-paramList["phi"])
  hiv_h0_g1 <- as.numeric(-paramList["rh"] * paramList["phi"])
  hiv_h1_g0 <- as.numeric(paramList["phi"])
  hiv_h1_g1 <- as.numeric(paramList["rh"] * paramList["phi"])

  paramList2 <- setNames(c(mort_g0,mort_g1,
                           lambda_g0,lambda_g1,
                           hiv_h0_g0,hiv_h0_g1,hiv_h1_g0,hiv_h1_g1),
                         c("mort_g0","mort_g1",
                           "lambda_g0","lambda_g1",
                           "hiv_h0_g0","hiv_h0_g1","hiv_h1_g0","hiv_h1_g1"))

  # SUSCEPTIBLES
  dS_h0_g0 <- as.numeric(paramList2["mort_g0"] -
                           (paramList2["lambda_g0"]+ paramList["mu0"]) * y["S_h0_g0"] +
                           (paramList2["hiv_h0_g0"] * y["S_h0_g0"]))
  dS_h0_g1 <- as.numeric(paramList2["mort_g1"] -
                           (paramList2["lambda_g1"] + paramList["mu0"]) * y["S_h0_g1"] +
                           (paramList2["hiv_h0_g1"] * y["S_h0_g1"]))
  dS_h1_g0 <- as.numeric(-(paramList2["lambda_g0"] + paramList["mu0"] + paramList["mu1"]) * y["S_h1_g0"] +
                           (paramList2["hiv_h1_g0"] * y["S_h0_g0"]))
  dS_h1_g1 <- as.numeric(-(paramList2["lambda_g1"] + paramList["mu0"] + paramList["mu1"]) * y["S_h1_g1"] +
                           (paramList2["hiv_h1_g1"] * y["S_h0_g1"]))

  # LATENT, RECENTLY INFECTED
  dL1_h0_g0 <- as.numeric(paramList2["lambda_g0"] * (y["S_h0_g0"] +
                                                       y["L2_h0_g0"] * (1-paramList["p"]) + y["R_h0_g0"]) -
                            (paramList["eta"] + paramList["zeta0"] + paramList["mu0"] + 0) * y["L1_h0_g0"] +
                            (paramList2["hiv_h0_g0"] * y["L1_h0_g0"]))
  dL1_h0_g1 <- as.numeric(paramList2["lambda_g1"] * (y["S_h0_g1"] +
                                                       y["L2_h0_g1"] * (1-paramList["p"]) + y["R_h0_g1"]) -
                            (paramList["eta"] + paramList["zeta0"] + paramList["mu0"] + 0) * y["L1_h0_g1"] +
                            (paramList2["hiv_h0_g1"] * y["L1_h0_g1"]))
  dL1_h1_g0 <- as.numeric(paramList2["lambda_g0"] * (y["S_h1_g0"] +
                                                       y["L2_h1_g0"] * (1-paramList["p"]) + y["R_h1_g0"]) -
                            (paramList["eta"] + paramList["zeta1"] + paramList["mu0"] + paramList["mu1"]) *
                            y["L1_h1_g0"] +(paramList2["hiv_h1_g0"] * y["L1_h0_g0"]))
  dL1_h1_g1 <- as.numeric(paramList2["lambda_g1"] * (y["S_h1_g1"] +
                                                       y["L2_h1_g1"] * (1-paramList["p"]) + y["R_h1_g1"]) -
                            (paramList["eta"] + paramList["zeta1"] + paramList["mu0"] + paramList["mu1"]) *
                            y["L1_h1_g1"] + (paramList2["hiv_h1_g1"] * y["L1_h0_g1"]))

  # LATENT, REMOTELY INFECTED
  dL2_h0_g0 <- as.numeric((paramList["eta"] * y["L1_h0_g0"]) -
                            (paramList2["lambda_g0"] * (1 - paramList["p"]) +
                               paramList["nu0"] + paramList["mu0"] + 0) * y["L2_h0_g0"] +
                            (paramList2["hiv_h0_g0"] * y["L2_h0_g0"]))
  dL2_h0_g1 <- as.numeric((paramList["eta"] * y["L1_h0_g1"]) -
                            (paramList2["lambda_g1"] * (1 - paramList["p"]) + paramList["nu0"] + paramList["mu0"] + 0) *
                            y["L2_h0_g1"] + (paramList2["hiv_h0_g1"] * y["L2_h0_g1"]))
  dL2_h1_g0 <- as.numeric((paramList["eta"] * y["L1_h1_g0"]) -
                            (paramList2["lambda_g0"] * (1 - paramList["p"]) +
                               paramList["nu1"] + paramList["mu0"] + paramList["mu1"]) * y["L2_h1_g0"] +
                            (paramList2["hiv_h1_g0"] * y["L2_h0_g0"]))
  dL2_h1_g1 <- as.numeric((paramList["eta"] * y["L1_h1_g1"]) -
                            (paramList2["lambda_g1"] * (1 - paramList["p"]) +
                               paramList["nu1"] + paramList["mu0"] + paramList["mu1"]) * y["L2_h1_g1"] +
                            (paramList2["hiv_h1_g1"] * y["L2_h0_g1"]))

  # ACTIVE TB
  dA_h0_g0 <- as.numeric((paramList["zeta0"] * y["L1_h0_g0"]) +
                           (paramList["nu0"] * y["L2_h0_g0"]) +
                           (paramList["psi"] * y["R_h0_g0"]) -
                           (paramList["rho0"] + paramList["mu0"] + 0 + paramList["mutb0"]) * y["A_h0_g0"] +
                           (paramList2["hiv_h0_g0"] * y["A_h0_g0"]))
  dA_h0_g1 <- as.numeric((paramList["zeta0"] * y["L1_h0_g1"]) +
                           (paramList["nu0"] * y["L2_h0_g1"]) +
                           (paramList["psi"] * y["R_h0_g1"]) -
                           (paramList["rho0"] + paramList["mu0"] + 0 + paramList["mutb0"]) * y["A_h0_g1"] +
                           (paramList2["hiv_h0_g1"] * y["A_h0_g1"]))
  dA_h1_g0 <- as.numeric((paramList["zeta1"] * y["L1_h1_g0"]) +
                           (paramList["nu1"] * y["L2_h1_g0"]) +
                           (paramList["psi"] * y["R_h1_g0"]) -
                           (paramList["rho1"] + paramList["mu0"] + paramList["mu1"] + paramList["mutb1"]) * y["A_h1_g0"] +
                           (paramList2["hiv_h1_g0"] * y["A_h0_g0"]))
  dA_h1_g1 <- as.numeric((paramList["zeta1"] * y["L1_h1_g1"]) +
                           (paramList["nu1"] * y["L2_h1_g1"]) +
                           (paramList["psi"] * y["R_h1_g1"]) -
                           (paramList["rho1"] + paramList["mu0"] + paramList["mu1"] + paramList["mutb1"]) * y["A_h1_g1"] +
                           (paramList2["hiv_h1_g1"] * y["A_h0_g1"]))

  # RECOVERED/TREATED
  dR_h0_g0 <- as.numeric((paramList["rho0"] * y["A_h0_g0"]) -
                           (paramList2["lambda_g0"] + paramList["psi"] + paramList["mu0"] + 0) * y["R_h0_g0"] +
                           (paramList2["hiv_h0_g0"] * y["R_h0_g0"]))
  dR_h0_g1 <- as.numeric((paramList["rho0"] * y["A_h0_g1"]) -
                           (paramList2["lambda_g1"] + paramList["psi"] + paramList["mu0"] + 0) * y["R_h0_g1"] +
                           (paramList2["hiv_h0_g1"] * y["R_h0_g1"]))
  dR_h1_g0 <- as.numeric((paramList["rho1"] * y["A_h1_g0"]) -
                           (paramList2["lambda_g0"] + paramList["psi"] + paramList["mu0"] + paramList["mu1"]) * y["R_h1_g0"] +
                           (paramList2["hiv_h1_g0"] * y["R_h0_g0"]))
  dR_h1_g1 <- as.numeric((paramList["rho1"] * y["A_h1_g1"]) -
                           (paramList2["lambda_g1"] + paramList["psi"] + paramList["mu0"] + paramList["mu1"]) * y["R_h1_g1"] +
                           (paramList2["hiv_h1_g1"] * y["R_h0_g1"]))

  dCI_h0_g0 <- as.numeric((paramList["zeta0"] * y["L1_h0_g0"]) +
                            (paramList["nu0"] * y["L2_h0_g0"]) +
                            (paramList["psi"] * y["R_h0_g0"]))

  dCI_h0_g1 <- as.numeric((paramList["zeta0"] * y["L1_h0_g1"]) +
                            (paramList["nu0"] * y["L2_h0_g1"]) +
                            (paramList["psi"] * y["R_h0_g1"]))

  dCI_h1_g0 <-  as.numeric((paramList["zeta1"] * y["L1_h1_g0"]) +
                             (paramList["nu1"] * y["L2_h1_g0"]) +
                             (paramList["psi"] * y["R_h1_g0"]))

  dCI_h1_g1 <- as.numeric((paramList["zeta1"] * y["L1_h1_g1"]) +
                            (paramList["nu1"] * y["L2_h1_g1"]) +
                            (paramList["psi"] * y["R_h1_g1"]))

  # Return the results as a list
  return(list(c(dS_h0_g0, dS_h0_g1, dS_h1_g0, dS_h1_g1,
                dL1_h0_g0, dL1_h0_g1, dL1_h1_g0, dL1_h1_g1,
                dL2_h0_g0, dL2_h0_g1, dL2_h1_g0, dL2_h1_g1,
                dA_h0_g0, dA_h0_g1, dA_h1_g0, dA_h1_g1,
                dR_h0_g0, dR_h0_g1, dR_h1_g0, dR_h1_g1,
                dCI_h0_g0, dCI_h0_g1, dCI_h1_g0, dCI_h1_g1)))
}

runSL2ARHOT<-function(beta0=DEF.BETA0,
                      beta1=DEF.BETA1,
                      rt=DEF.RT,
                      z1=DEF.Z1,
                      z0=DEF.Z0,
                      zeta1=DEF.ZETA1,
                      nu1=DEF.NU1,
                      mutb0=DEF.MUTB0,
                      mutb1=DEF.MUTB1,
                      rho0=DEF.RHO0,
                      rho1=DEF.RHO1,
                      phi=DEF.PHI,
                      mu1=DEF.MU1,
                      psi=DEF.PSI,
                      ri=DEF.RI,
                      p=DEF.P,
                      eta=DEF.ETA,
                      zeta0=DEF.ZETA0,
                      nu0=DEF.NU0,
                      mu0=DEF.MU0,
                      rh=DEF.RH,
                      initial.state= c(S_h0_g0=DEF.S_H0_G0, S_h0_g1=DEF.S_H0_G1,
                                       S_h1_g0=DEF.S_H1_G0, S_h1_g1=DEF.S_H1_G1,
                                       L1_h0_g0=DEF.L1_H0_G0, L1_h0_g1=DEF.L1_H0_G1,
                                       L1_h1_g0=DEF.L1_H1_G0, L1_h1_g1=DEF.L1_H1_G1,
                                       L2_h0_g0=DEF.L2_H0_G0, L2_h0_g1=DEF.L2_H0_G1,
                                       L2_h1_g0=DEF.L2_H1_G0, L2_h1_g1=DEF.L2_H1_G1,
                                       A_h0_g0=DEF.A_H0_G0, A_h0_g1=DEF.A_H0_G1,
                                       A_h1_g0=DEF.A_H1_G0, A_h1_g1=DEF.A_H1_G1,
                                       R_h0_g0=DEF.R_H0_G0, R_h0_g1=DEF.R_H0_G1,
                                       R_h1_g0=DEF.R_H1_G0, R_h1_g1=DEF.R_H1_G1,
                                       CI_h0_g0=0, CI_h0_g1=0, CI_h1_g0=0, CI_h1_g1=0),
                      max.time = MAX.TIME) {

  # Parameter Matrix
  param <- c(beta0=beta0,
             beta1=beta1,
             rt=rt,
             z1=z1,
             z0=z0,
             zeta1=zeta1,
             nu1=nu1,
             mutb0=mutb0,
             mutb1=mutb1,
             rho0=rho0,
             rho1=rho1,
             phi=phi,
             mu1=mu1,
             psi=psi,
             ri=ri,
             p=p,
             eta=eta,
             zeta0=zeta0,
             nu0=nu0,
             mu0=mu0,
             rh=rh)

  #Sequence of times at which we want estimates..
  times <- seq(0,max.time,1)

  #Run the ODE-Solver
  sir.output <- lsoda(initial.state, times, dx.dt.SL2ARHOT, param)

  #now we will return the output.
  return(sir.output)
}

addCompartments <-function(run) {
  run <- as.data.frame(run)
  run$S <- run$S_h0_g0 + run$S_h1_g0 + run$S_h0_g1 + run$S_h1_g1
  run$L1 <- run$L1_h0_g0 + run$L1_h1_g0 + run$L1_h0_g1 + run$L1_h1_g1
  run$L2 <- run$L2_h0_g0 + run$L2_h1_g0 + run$L2_h0_g1 + run$L2_h1_g1
  run$A <- run$A_h0_g0 + run$A_h1_g0 + run$A_h0_g1 + run$A_h1_g1
  run$R <- run$R_h0_g0 + run$R_h1_g0 + run$R_h0_g1 + run$R_h1_g1
  run$CI <- run$CI_h0_g0 + run$CI_h1_g0 + run$CI_h0_g1 + run$CI_h1_g1
  run$CI0 <- run$CI_h0_g0 + run$CI_h1_g0
  run$CI1 <- run$CI_h0_g1 + run$CI_h1_g1
  return(run)
}

MAXPOP = 1
a=0.04
b=0.0003
DEF.S_H0_G0 <- (0.6985-a)*MAXPOP
DEF.S_H0_G1 <- (0.0010-b)*MAXPOP
DEF.S_H1_G0 <- (0.30-a)*MAXPOP
DEF.S_H1_G1 <- (0.0005-b)*MAXPOP
DEF.L1_H0_G0 <- 0*MAXPOP
DEF.L1_H0_G1 <- 0*MAXPOP
DEF.L1_H1_G0 <- 0*MAXPOP
DEF.L1_H1_G1 <- 0*MAXPOP
DEF.L2_H0_G0 <- 0*MAXPOP
DEF.L2_H0_G1 <- 0*MAXPOP
DEF.L2_H1_G0 <- 0*MAXPOP
DEF.L2_H1_G1 <- 0*MAXPOP
DEF.A_H0_G0 <- a*MAXPOP
DEF.A_H0_G1 <- b*MAXPOP
DEF.A_H1_G0 <- a*MAXPOP
DEF.A_H1_G1 <- b*MAXPOP
DEF.R_H0_G0 <- 0.00*MAXPOP
DEF.R_H0_G1 <- 0.00*MAXPOP
DEF.R_H1_G0 <- 0.00*MAXPOP
DEF.R_H1_G1 <- 0.00*MAXPOP

run0 <- runSL2ARHOT()
run0 <- addCompartments(run0)
r0_incid_5y_diff <- (run0$CI[505] - run0$CI[500]) * 100000
r0_hotincid_5y_diff <- (run0$CI1[505] - run0$CI1[500])/DEF.Z1 * 100000
g1 <- ggplot(data=run0, aes(x=time, color=State)) +
  geom_point(aes(y=S, color="Susceptible")) +
  geom_point(aes(y=L1, color="Latent, recent")) +
  geom_point(aes(y=L2, color= "Latent, remote")) +
  geom_point(aes(y=A, color="Active")) +
  geom_point(aes(y=R, color="Recovered")) +
  #scale_color_manual(values = cbPalette[4:8]) +
  scale_color_brewer(palette = "Set1") +
  #geom_point(aes(y=CI, color="Cumulative Incidence")) +
  #theme_tufte() +
  theme(legend.position="bottom") +
  labs(x = "Time (years)", y = "Proportion of Population", title = "TB Hotspot/Community Model: Gauteng \n") +
  theme(text = element_text(size=25))
g1
r0_incid_5y_diff/5
r0_hotincid_5y_diff/5


run1 <- runSL2ARHOT(nu1=0.19*DEF.NU1)
run1 <- addCompartments(run1)
r1_incid_5y_diff <- (run1$CI[505] - run1$CI[500]) * 100000
run2 <- runSL2ARHOT(beta1 = DEF.BETA0)
run2 <- addCompartments(run2)
r2_incid_5y_diff <- (run2$CI[505] - run2$CI[500]) * 100000


#### SENSITIVITY ANALYSES ####
# Create data frame of percent differences
varList <- c("beta0", "beta1", "rt", "zeta1", "mutb0", "mutb1",
             "rho0", "phi", "mu1", "psi","ri", "p", "eta", "zeta0", "nu0", "mu0", "rh")
vList <- 1:17
sensList <- data.frame(Variable=varList, Low=vList, High=vList)

# Baseline Scenario
years <- 501:550
time <- 1:50
baseline <- runSL2ARHOT()
baseline <- addCompartments(baseline)
b_incid_5y_diff <- (baseline$CI[505] - baseline$CI[500]) * 100000
b_incid_1y <- (baseline$CI[years]-baseline$CI[years-1])*100000

# Intervention 1
intervention <- runSL2ARHOT(beta1=DEF.BETA0)
intervention <- addCompartments(intervention)
i_incid_5y_diff <- (intervention$CI[505] - intervention$CI[500]) * 100000
i_diff = (b_incid_5y_diff-i_incid_5y_diff)

i_incid_1y <- (intervention$CI[years]-intervention$CI[years-1])*100000
i_diff1 <- (b_incid_1y-i_incid_1y)
dat = data.frame(time=time,community=b_incid_1y,intervention1=i_incid_1y,i_diff1=i_diff1)

# Intervention 2
intervention <- runSL2ARHOT(beta1=DEF.BETA0,rh=1)
intervention <- addCompartments(intervention)
i_incid_5y_diff <- (intervention$CI[505] - intervention$CI[500]) * 100000
i_diff = (b_incid_5y_diff-i_incid_5y_diff)

i_incid_1y <- (intervention$CI[years]-intervention$CI[years-1])*100000
i_diff1 <- (b_incid_1y-i_incid_1y)
dat = cbind(dat,intervention2=i_incid_1y,i_diff2=i_diff1)

# Intervention 3
intervention <- runSL2ARHOT(rho0=DEF.RHO0*1.25, rho1=DEF.RHO1*1.25)
intervention <- addCompartments(intervention)
i_incid_5y_diff <- (intervention$CI[505] - intervention$CI[500]) * 100000
i_diff = (b_incid_5y_diff-i_incid_5y_diff)

i_incid_1y <- (intervention$CI[years]-intervention$CI[years-1])*100000
i_diff1 <- (b_incid_1y-i_incid_1y)
dat = cbind(dat,intervention3=i_incid_1y,i_diff3=i_diff1)

# Intervention 4
intervention <- runSL2ARHOT(beta0=DEF.BETA0*0.90)
intervention <- addCompartments(intervention)
i_incid_5y_diff <- (intervention$CI[505] - intervention$CI[500]) * 100000
i_diff = (b_incid_5y_diff-i_incid_5y_diff)

i_incid_1y <- (intervention$CI[years]-intervention$CI[years-1])*100000
i_diff1 <- (b_incid_1y-i_incid_1y)
dat = cbind(dat,intervention4=i_incid_1y,i_diff4=i_diff1)

g2 <- ggplot(data=dat, aes(x=time, color=Package)) +
  geom_line(aes(y=community, color="Baseline"),size=1, linetype="solid") +
  geom_line(aes(y=intervention1, color="Reduce B1"),size=1, linetype="dashed") +
  geom_line(aes(y=intervention2, color="Reduce B1 + Reduce RH"),size=1, linetype="dotted") +
  geom_line(aes(y=intervention3, color="25% Increase of Screen+Treat"),size=1, linetype="dotdash") +
  geom_line(aes(y=intervention4, color="10% Reduction in B0"),size=1, linetype="longdash") +
  #scale_color_manual(values = cbPalette[4:8]) +
  scale_color_brewer(palette = "Set1") +
  #theme_few() +
  #geom_point(aes(y=CI, color="Cumulative Incidence")) +
  theme(legend.position="bottom") +
  ylim(0,400) +
  labs(x = "Time (years)", y = "TB Incidence per 100000",
       title = "Reductions in TB Incidence after Successful Interventions, Gauteng \n") +
  theme(text = element_text(size=25))
g2


# beta0 value=1.77 range 1.33-2.21
DEF.BETA0=1.33
baseline <- runSL2ARHOT()
baseline <- addCompartments(baseline)
b_incid_5y_diff <- (baseline$CI[505] - baseline$CI[500]) * 100000

intervention <- runSL2ARHOT(beta1=DEF.BETA0)
intervention <- addCompartments(intervention)
i_incid_5y_diff <- (intervention$CI[505] - intervention$CI[500]) * 100000

i_diff = (b_incid_5y_diff-i_incid_5y_diff)
sensList[1,2] <- i_diff

DEF.BETA0=2.21
baseline <- runSL2ARHOT()
baseline <- addCompartments(baseline)
b_incid_5y_diff <- (baseline$CI[505] - baseline$CI[500]) * 100000

intervention <- runSL2ARHOT(beta1=DEF.BETA0)
intervention <- addCompartments(intervention)
i_incid_5y_diff <- (intervention$CI[505] - intervention$CI[500]) * 100000

i_diff = (b_incid_5y_diff-i_incid_5y_diff)
sensList[1,3] <- i_diff
DEF.BETA0=1.77

# beta1 value=2.427 range 1.82-3.03
DEF.BETA1=1.82
baseline <- runSL2ARHOT()
baseline <- addCompartments(baseline)
b_incid_5y_diff <- (baseline$CI[505] - baseline$CI[500]) * 100000

intervention <- runSL2ARHOT(beta1=DEF.BETA0)
intervention <- addCompartments(intervention)
i_incid_5y_diff <- (intervention$CI[505] - intervention$CI[500]) * 100000

i_diff = (b_incid_5y_diff-i_incid_5y_diff)
sensList[2,2] <- i_diff

DEF.BETA1=3.03
baseline <- runSL2ARHOT()
baseline <- addCompartments(baseline)
b_incid_5y_diff <- (baseline$CI[505] - baseline$CI[500]) * 100000

intervention <- runSL2ARHOT(beta1=DEF.BETA0)
intervention <- addCompartments(intervention)
i_incid_5y_diff <- (intervention$CI[505] - intervention$CI[500]) * 100000

i_diff = (b_incid_5y_diff-i_incid_5y_diff)
sensList[2,3] <- i_diff
DEF.BETA1=2.427

# rt value=0.03 range 0.01-0.063
DEF.RT = 0.01
baseline <- runSL2ARHOT()
baseline <- addCompartments(baseline)
b_incid_5y_diff <- (baseline$CI[505] - baseline$CI[500]) * 100000

intervention <- runSL2ARHOT(beta1=DEF.BETA0)
intervention <- addCompartments(intervention)
i_incid_5y_diff <- (intervention$CI[505] - intervention$CI[500]) * 100000

i_diff = (b_incid_5y_diff-i_incid_5y_diff)
sensList[3,2] <- i_diff

DEF.RT = 0.063
baseline <- runSL2ARHOT()
baseline <- addCompartments(baseline)
b_incid_5y_diff <- (baseline$CI[505] - baseline$CI[500]) * 100000

intervention <- runSL2ARHOT(beta1=DEF.BETA0)
intervention <- addCompartments(intervention)
i_incid_5y_diff <- (intervention$CI[505] - intervention$CI[500]) * 100000

i_diff = (b_incid_5y_diff-i_incid_5y_diff)
sensList[3,3] <- i_diff
DEF.RT = 0.03

# zeta1 value=0.31 range 0.15-0.82 / nu1 value=0.08 range 0.02-0.21
DEF.ZETA1 = 0.15
DEF.NU1 = 0.02
baseline <- runSL2ARHOT()
baseline <- addCompartments(baseline)
b_incid_5y_diff <- (baseline$CI[505] - baseline$CI[500]) * 100000

intervention <- runSL2ARHOT(beta1=DEF.BETA0)
intervention <- addCompartments(intervention)
i_incid_5y_diff <- (intervention$CI[505] - intervention$CI[500]) * 100000

i_diff = (b_incid_5y_diff-i_incid_5y_diff)
sensList[4,2] <- i_diff

DEF.ZETA1 = 0.82
DEF.NU1 = 0.21
baseline <- runSL2ARHOT()
baseline <- addCompartments(baseline)
b_incid_5y_diff <- (baseline$CI[505] - baseline$CI[500]) * 100000

intervention <- runSL2ARHOT(beta1=DEF.BETA0)
intervention <- addCompartments(intervention)
i_incid_5y_diff <- (intervention$CI[505] - intervention$CI[500]) * 100000

i_diff = (b_incid_5y_diff-i_incid_5y_diff)
sensList[4,3] <- i_diff
DEF.ZETA1 = 0.31
DEF.NU1 = 0.08

# mutb0 value=48/100000 range 28-73
DEF.MUTB0 = 28/100000
baseline <- runSL2ARHOT()
baseline <- addCompartments(baseline)
b_incid_5y_diff <- (baseline$CI[505] - baseline$CI[500]) * 100000

intervention <- runSL2ARHOT(beta1=DEF.BETA0)
intervention <- addCompartments(intervention)
i_incid_5y_diff <- (intervention$CI[505] - intervention$CI[500]) * 100000

i_diff = (b_incid_5y_diff-i_incid_5y_diff)
sensList[5,2] <- i_diff

DEF.MUTB0 = 73/100000
baseline <- runSL2ARHOT()
baseline <- addCompartments(baseline)
b_incid_5y_diff <- (baseline$CI[505] - baseline$CI[500]) * 100000

intervention <- runSL2ARHOT(beta1=DEF.BETA0)
intervention <- addCompartments(intervention)
i_incid_5y_diff <- (intervention$CI[505] - intervention$CI[500]) * 100000

i_diff = (b_incid_5y_diff-i_incid_5y_diff)
sensList[5,3] <- i_diff
DEF.MUTB0 = 48/100000

# mutb1 value=121/100000 range 90-158
DEF.MUTB1 = 90/100000
baseline <- runSL2ARHOT()
baseline <- addCompartments(baseline)
b_incid_5y_diff <- (baseline$CI[505] - baseline$CI[500]) * 100000

intervention <- runSL2ARHOT(beta1=DEF.BETA0)
intervention <- addCompartments(intervention)
i_incid_5y_diff <- (intervention$CI[505] - intervention$CI[500]) * 100000

i_diff = (b_incid_5y_diff-i_incid_5y_diff)
sensList[6,2] <- i_diff

DEF.MUTB1 = 158/100000
baseline <- runSL2ARHOT()
baseline <- addCompartments(baseline)
b_incid_5y_diff <- (baseline$CI[505] - baseline$CI[500]) * 100000

intervention <- runSL2ARHOT(beta1=DEF.BETA0)
intervention <- addCompartments(intervention)
i_incid_5y_diff <- (intervention$CI[505] - intervention$CI[500]) * 100000

i_diff = (b_incid_5y_diff-i_incid_5y_diff)
sensList[6,3] <- i_diff
DEF.MUTB1 = 121/100000

# rho0 value=0.53 range 0.46-0.59 / rho1 value=1.06 range 0.92-1.18
DEF.RHO0 = 0.46
REF.RHO1 = 1.06
baseline <- runSL2ARHOT()
baseline <- addCompartments(baseline)
b_incid_5y_diff <- (baseline$CI[505] - baseline$CI[500]) * 100000

intervention <- runSL2ARHOT(beta1=DEF.BETA0)
intervention <- addCompartments(intervention)
i_incid_5y_diff <- (intervention$CI[505] - intervention$CI[500]) * 100000

i_diff = (b_incid_5y_diff-i_incid_5y_diff)
sensList[7,2] <- i_diff

DEF.RHO0 = 0.59
REF.RHO1 = 1.18
baseline <- runSL2ARHOT()
baseline <- addCompartments(baseline)
b_incid_5y_diff <- (baseline$CI[505] - baseline$CI[500]) * 100000

intervention <- runSL2ARHOT(beta1=DEF.BETA0)
intervention <- addCompartments(intervention)
i_incid_5y_diff <- (intervention$CI[505] - intervention$CI[500]) * 100000

i_diff = (b_incid_5y_diff-i_incid_5y_diff)
sensList[7,3] <- i_diff
DEF.RHO0 = 0.53
REF.RHO1 = 1.06

# phi value=0.01 range 0.0075-0.0125
DEF.PHI = 0.0075
baseline <- runSL2ARHOT()
baseline <- addCompartments(baseline)
b_incid_5y_diff <- (baseline$CI[505] - baseline$CI[500]) * 100000

intervention <- runSL2ARHOT(beta1=DEF.BETA0)
intervention <- addCompartments(intervention)
i_incid_5y_diff <- (intervention$CI[505] - intervention$CI[500]) * 100000

i_diff = (b_incid_5y_diff-i_incid_5y_diff)
sensList[8,2] <- i_diff

DEF.PHI = 0.0125
baseline <- runSL2ARHOT()
baseline <- addCompartments(baseline)
b_incid_5y_diff <- (baseline$CI[505] - baseline$CI[500]) * 100000

intervention <- runSL2ARHOT(beta1=DEF.BETA0)
intervention <- addCompartments(intervention)
i_incid_5y_diff <- (intervention$CI[505] - intervention$CI[500]) * 100000

i_diff = (b_incid_5y_diff-i_incid_5y_diff)
sensList[8,3] <- i_diff
DEF.PHI = 0.01

# mu1 value=0.007 range 0.00525-0.00875
DEF.MU1=0.00525
baseline <- runSL2ARHOT()
baseline <- addCompartments(baseline)
b_incid_5y_diff <- (baseline$CI[505] - baseline$CI[500]) * 100000

intervention <- runSL2ARHOT(beta1=DEF.BETA0)
intervention <- addCompartments(intervention)
i_incid_5y_diff <- (intervention$CI[505] - intervention$CI[500]) * 100000

i_diff = (b_incid_5y_diff-i_incid_5y_diff)
sensList[9,2] <- i_diff

DEF.MU1=0.00875
baseline <- runSL2ARHOT()
baseline <- addCompartments(baseline)
b_incid_5y_diff <- (baseline$CI[505] - baseline$CI[500]) * 100000

intervention <- runSL2ARHOT(beta1=DEF.BETA0)
intervention <- addCompartments(intervention)
i_incid_5y_diff <- (intervention$CI[505] - intervention$CI[500]) * 100000

i_diff = (b_incid_5y_diff-i_incid_5y_diff)
sensList[9,3] <- i_diff
DEF.MU1=0.007

# psi value=0.0083 range 0.0063-0.0187
DEF.PSI = 0.0063
baseline <- runSL2ARHOT()
baseline <- addCompartments(baseline)
b_incid_5y_diff <- (baseline$CI[505] - baseline$CI[500]) * 100000

intervention <- runSL2ARHOT(beta1=DEF.BETA0)
intervention <- addCompartments(intervention)
i_incid_5y_diff <- (intervention$CI[505] - intervention$CI[500]) * 100000

i_diff = (b_incid_5y_diff-i_incid_5y_diff)
sensList[10,2] <- i_diff

DEF.PSI = 0.0187
baseline <- runSL2ARHOT()
baseline <- addCompartments(baseline)
b_incid_5y_diff <- (baseline$CI[505] - baseline$CI[500]) * 100000

intervention <- runSL2ARHOT(beta1=DEF.BETA0)
intervention <- addCompartments(intervention)
i_incid_5y_diff <- (intervention$CI[505] - intervention$CI[500]) * 100000

i_diff = (b_incid_5y_diff-i_incid_5y_diff)
sensList[10,3] <- i_diff
DEF.PSI = 0.0083

# ri value=0.68 range 0-1
DEF.RI=0
baseline <- runSL2ARHOT()
baseline <- addCompartments(baseline)
b_incid_5y_diff <- (baseline$CI[505] - baseline$CI[500]) * 100000

intervention <- runSL2ARHOT(beta1=DEF.BETA0)
intervention <- addCompartments(intervention)
i_incid_5y_diff <- (intervention$CI[505] - intervention$CI[500]) * 100000

i_diff = (b_incid_5y_diff-i_incid_5y_diff)
sensList[11,2] <- i_diff

DEF.RI=1
baseline <- runSL2ARHOT()
baseline <- addCompartments(baseline)
b_incid_5y_diff <- (baseline$CI[505] - baseline$CI[500]) * 100000

intervention <- runSL2ARHOT(beta1=DEF.BETA0)
intervention <- addCompartments(intervention)
i_incid_5y_diff <- (intervention$CI[505] - intervention$CI[500]) * 100000

i_diff = (b_incid_5y_diff-i_incid_5y_diff)
sensList[11,3] <- i_diff
DEF.RI=0.68

# p value=0.56 range 0.42-0.70
DEF.P=0.42
baseline <- runSL2ARHOT()
baseline <- addCompartments(baseline)
b_incid_5y_diff <- (baseline$CI[505] - baseline$CI[500]) * 100000

intervention <- runSL2ARHOT(beta1=DEF.BETA0)
intervention <- addCompartments(intervention)
i_incid_5y_diff <- (intervention$CI[505] - intervention$CI[500]) * 100000

i_diff = (b_incid_5y_diff-i_incid_5y_diff)
sensList[12,2] <- i_diff

DEF.P=0.70
baseline <- runSL2ARHOT()
baseline <- addCompartments(baseline)
b_incid_5y_diff <- (baseline$CI[505] - baseline$CI[500]) * 100000

intervention <- runSL2ARHOT(beta1=DEF.BETA0)
intervention <- addCompartments(intervention)
i_incid_5y_diff <- (intervention$CI[505] - intervention$CI[500]) * 100000

i_diff = (b_incid_5y_diff-i_incid_5y_diff)
sensList[12,3] <- i_diff
DEF.P=0.56

# eta value=1/5 range 1/4.2-1/6.2
DEF.ETA=1/4.2
baseline <- runSL2ARHOT()
baseline <- addCompartments(baseline)
b_incid_5y_diff <- (baseline$CI[505] - baseline$CI[500]) * 100000

intervention <- runSL2ARHOT(beta1=DEF.BETA0)
intervention <- addCompartments(intervention)
i_incid_5y_diff <- (intervention$CI[505] - intervention$CI[500]) * 100000

i_diff = (b_incid_5y_diff-i_incid_5y_diff)
sensList[13,2] <- i_diff

DEF.ETA=1/6.2
baseline <- runSL2ARHOT()
baseline <- addCompartments(baseline)
b_incid_5y_diff <- (baseline$CI[505] - baseline$CI[500]) * 100000

intervention <- runSL2ARHOT(beta1=DEF.BETA0)
intervention <- addCompartments(intervention)
i_incid_5y_diff <- (intervention$CI[505] - intervention$CI[500]) * 100000

i_diff = (b_incid_5y_diff-i_incid_5y_diff)
sensList[13,3] <- i_diff
DEF.ETA=1/5

# zeta0 value=0.03 range 0.026-0.036
DEF.ZETA0=0.026
baseline <- runSL2ARHOT()
baseline <- addCompartments(baseline)
b_incid_5y_diff <- (baseline$CI[505] - baseline$CI[500]) * 100000

intervention <- runSL2ARHOT(beta1=DEF.BETA0)
intervention <- addCompartments(intervention)
i_incid_5y_diff <- (intervention$CI[505] - intervention$CI[500]) * 100000

i_diff = (b_incid_5y_diff-i_incid_5y_diff)
sensList[14,2] <- i_diff

DEF.ZETA0=0.036
baseline <- runSL2ARHOT()
baseline <- addCompartments(baseline)
b_incid_5y_diff <- (baseline$CI[505] - baseline$CI[500]) * 100000

intervention <- runSL2ARHOT(beta1=DEF.BETA0)
intervention <- addCompartments(intervention)
i_incid_5y_diff <- (intervention$CI[505] - intervention$CI[500]) * 100000

i_diff = (b_incid_5y_diff-i_incid_5y_diff)
sensList[14,3] <- i_diff
DEF.ZETA0=0.03

# nu0 value=0.0005 range 0.0002-0.0011
DEF.NU0=0.0002
baseline <- runSL2ARHOT()
baseline <- addCompartments(baseline)
b_incid_5y_diff <- (baseline$CI[505] - baseline$CI[500]) * 100000

intervention <- runSL2ARHOT(beta1=DEF.BETA0)
intervention <- addCompartments(intervention)
i_incid_5y_diff <- (intervention$CI[505] - intervention$CI[500]) * 100000

i_diff = (b_incid_5y_diff-i_incid_5y_diff)
sensList[15,2] <- i_diff

DEF.NU0=0.0011
baseline <- runSL2ARHOT()
baseline <- addCompartments(baseline)
b_incid_5y_diff <- (baseline$CI[505] - baseline$CI[500]) * 100000

intervention <- runSL2ARHOT(beta1=DEF.BETA0)
intervention <- addCompartments(intervention)
i_incid_5y_diff <- (intervention$CI[505] - intervention$CI[500]) * 100000

i_diff = (b_incid_5y_diff-i_incid_5y_diff)
sensList[15,3] <- i_diff
DEF.NU0=0.0005

# mu0 value=1/56 range 1/27-1/100
DEF.MU0=1/27
baseline <- runSL2ARHOT()
baseline <- addCompartments(baseline)
b_incid_5y_diff <- (baseline$CI[505] - baseline$CI[500]) * 100000

intervention <- runSL2ARHOT(beta1=DEF.BETA0)
intervention <- addCompartments(intervention)
i_incid_5y_diff <- (intervention$CI[505] - intervention$CI[500]) * 100000

i_diff = (b_incid_5y_diff-i_incid_5y_diff)
sensList[16,2] <- i_diff

DEF.MU0=1/100
baseline <- runSL2ARHOT()
baseline <- addCompartments(baseline)
b_incid_5y_diff <- (baseline$CI[505] - baseline$CI[500]) * 100000

intervention <- runSL2ARHOT(beta1=DEF.BETA0)
intervention <- addCompartments(intervention)
i_incid_5y_diff <- (intervention$CI[505] - intervention$CI[500]) * 100000

i_diff = (b_incid_5y_diff-i_incid_5y_diff)
sensList[16,3] <- i_diff
DEF.MU0=1/56

# rh value=2.65 range 1.0-5.0
DEF.RH=1.0
baseline <- runSL2ARHOT()
baseline <- addCompartments(baseline)
b_incid_5y_diff <- (baseline$CI[505] - baseline$CI[500]) * 100000

intervention <- runSL2ARHOT(beta1=DEF.BETA0)
intervention <- addCompartments(intervention)
i_incid_5y_diff <- (intervention$CI[505] - intervention$CI[500]) * 100000

i_diff = (b_incid_5y_diff-i_incid_5y_diff)
sensList[17,2] <- i_diff

DEF.RH=5.0
baseline <- runSL2ARHOT()
baseline <- addCompartments(baseline)
b_incid_5y_diff <- (baseline$CI[505] - baseline$CI[500]) * 100000

intervention <- runSL2ARHOT(beta1=DEF.BETA0)
intervention <- addCompartments(intervention)
i_incid_5y_diff <- (intervention$CI[505] - intervention$CI[500]) * 100000

i_diff = (b_incid_5y_diff-i_incid_5y_diff)
sensList[17,3] <- i_diff
DEF.RH=2.13

sensList

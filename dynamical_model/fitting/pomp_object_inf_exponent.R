library(pomp)



simul_pf <- pomp::Csnippet(r"{
  // compute transmission rate
  double betaIN = exp(b1*season1+b2*season2+b3*season3+b4*season4+
    b5*season5+b6*season6+season4*bH*covariate);

  // gamma white noise
  double dW = rgammawn(sigPRO,dt);

  // force of infection
  double foi = (betaOUT+pow((I1+q0*I2)/pop,alpha)*betaIN)*(dW/dt);

  double dBS1 = (delta*pop+dpopdt)*dt; // FIXME: POSSIBLY MODEL POPULATION SIZE AS STOCHASTIC?
  double dS2S1 = muS2S1*S2*dt;
  double dS1E = F*S1*dt;

  double dEI1 = muEI1*E*dt;
  double dI1S2 = muI1S2*I1*dt;
  double dS2I2 = F*S2*dt;
  double dI2S2 = muI2S2*I2*dt;
  double dS1D = delta*S1*dt;
  double dED = delta*E*dt;
  double dI1D = delta*I1*dt;
  double dS2D = delta*S2*dt;
  double dI2D = delta*I2*dt;
  double dKK = ((foi-K)/(tau/2.0))*dt;
  double dFF = ((K-F)/(tau/2.0))*dt;

  // compute equations
  S1 += dBS1 + dS2S1 - dS1E - dS1D;
  E  += dS1E - dEI1 - dED;
  I1 += dEI1 - dI1S2 - dI1D;
  S2 += dI1S2 - dS2S1 - dS2I2 +dI2S2 - dS2D;
  I2 += dS2I2 - dI2S2 - dI2D;
  K  += dKK;
  F  += dFF;
  cases += rho*dEI1;
  W += (dW-dt)/sigPRO;
}")

############ rmeas #################
rmeas_pf <- pomp::Csnippet("
  double size = 1.0/sigOBS/sigOBS;
  PF = rnbinom_mu(size,cases);
")

############ dmeas #################
dmeas_pf <- pomp::Csnippet("
  double size = 1.0/sigOBS/sigOBS;
  if(ISNA(PF)){
    if(!give_log){
      lik = 1;
    } else {
      lik = 0;
    }
  } else{
  lik = dnbinom_mu(PF, size, cases + 0.1, 1);
  }
  if (!give_log) lik = exp(lik);
")
############ initlz #################

initlz_pf <- pomp::Csnippet("
  double m = pop/(S1_0 + E_0 + I1_0 + S2_0 + I2_0);

  S1 = nearbyint(m*S1_0);
  E = nearbyint(m*E_0);
  I1 = nearbyint(m*I1_0);
  S2 = nearbyint(m*S2_0);
  I2 = nearbyint(m*I2_0);

  K = K_0;
  F = F_0;
  cases = 0;
  W = 0;
")

par_names <- c("sigOBS", "sigPRO", "muS2S1", "muEI1", "muI1S2", "muI2S2", "betaOUT", "delta", "rho", "tau", "q0", "alpha", "b1", "b2", "b3", "b4", "b5", "b6", "bH")
vp_names <- c("S1_0", "E_0", "I1_0", "I2_0", "S2_0", "K_0", "F_0")
log_transf <- c("muS2S1", "muEI1", "muI1S2", "muI2S2", "betaOUT", "sigOBS", "sigPRO", "tau", "K_0", "F_0")
logit_transf <- c("rho", "q0", "alpha")
barycentric_transf <- c("S1_0", "E_0", "I1_0", "S2_0", "I2_0")
state_names <- c("cases", "S1", "E", "I1", "S2", "I2", "K", "F", "W")

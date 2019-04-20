library(RNetLogo)

#setup path to load NetLogo model
nl.path <- "C:/Program Files/NetLogo 6.0.3/app"
nl.jarname <- "netlogo-6.0.3.jar"
model.path <- paste(getwd(), "Bacteria-Rust model.nlogo", sep = "/")

#set initial conditions for NetLogo model
NLStart(nl.path = nl.path,
        gui = FALSE,
        nl.jarname = nl.jarname)
NLLoadModel(model.path)

totalTime <- 168
alphaGlobal <- 0.9
betaGlobal <- 0.5
hGlobal <- 0.15
gammaGlobal <- 0.5
muGlobal <- 0.5
psiGlobal <- 0.1

alphaLocal <- 0.5
betaLocal <- 0.8
hLocal <- 0.07
gammaLocal <- 0.3
muLocal <- 0.6
deltaLocal <- 0.3
psiLocal <- 1

n <- 100

NLCommand("setup")
initialSetup <- NLReport("[proportion] of turtles")

#run one simulation of NetLogo model
runSimulation <-
  function(ticks = totalTime,
           alpha,
           beta,
           h,
           gamma,
           mu = mu0,
           psi,
           delta,
           global) {
    NLCommand("setup")
    #set initial conditions
    NLCommand(paste("set global-control ", global))
    NLCommand(paste("set immigration-rate ", alpha / 24))
    NLCommand(paste("set emmigration-rate ", beta / 24))
    NLCommand(paste("set growth-rate-rust ", h))
    NLCommand(paste("set conversion-rate ", gamma / 24))
    NLCommand(paste("set prop-prob ", psi))
    NLCommand(paste("set main-irrigation ", mu))
    NLCommand(paste("set secondary-irrigation ", delta))
    
    totalTurtles <- NLReport("count turtles")
    for (i in 1:totalTurtles) {
      NLCommand(paste("ask turtle", i - 1, "[set proportion", initialSetup[i], "]"))
    }
    
    #generate time series of proportion of rust at each plant
    turtleSeries <- matrix(0L, nrow = ticks, ncol = totalTurtles)
    for (i in 1:ticks) {
      NLCommand("go")
      turtleSeries[i, ] <- NLReport("[proportion] of turtles")
    }
    
    return(turtleSeries)
  }

observationsGlobal <-
  runSimulation(
    alpha = alphaGlobal,
    beta = betaGlobal,
    h = hGlobal,
    gamma = gammaGlobal,
    mu = muGlobal,
    psi = psiGlobal,
    delta = 0,
    global = TRUE
  )
observationsLocal <-
  runSimulation(
    alpha = alphaLocal,
    beta = betaLocal,
    h = hLocal,
    gamma = gammaLocal,
    mu = muLocal,
    psi = psiLocal,
    delta = deltaLocal,
    global = FALSE
  )

#calculate Equation (14)
calculateNorm <-
  function(x,
           global,
           alpha,
           beta,
           h,
           gamma,
           mu,
           psi,
           delta) {
    simulation <- runSimulation(
      alpha = alpha,
      beta = beta,
      h = h,
      gamma = gamma,
      mu = mu,
      psi = psi,
      delta = delta,
      global = global
    )
    
    if (global) {
      val <- norm(observationsGlobal - simulation, "f")
    } else{
      val <- norm(observationsLocal - simulation, "f")
    }
    return(val)
  }

#run ABC procedure for simulation
runProcedure <- function(param, global, expected, tol) {
  library(EasyABC)
  
  model <- function(x) {
    if (global) {
      simulation <- switch(
        param,
        "alpha" = calculateNorm(
          global = TRUE,
          alpha = x,
          beta = betaGlobal,
          h = hGlobal,
          gamma = gammaGlobal,
          mu = muGlobal,
          psi = psiGlobal,
          delta = 0
        ),
        "beta" = calculateNorm(
          global = TRUE,
          alpha = alphaGlobal,
          beta = x,
          h = hGlobal,
          gamma = gammaGlobal,
          mu = muGlobal,
          psi = psiGlobal,
          delta = 0
        ),
        "h" = calculateNorm(
          global = TRUE,
          alpha = alphaGlobal,
          beta = betaGlobal,
          h = x,
          gamma = gammaGlobal,
          mu = muGlobal,
          psi = psiGlobal,
          delta = 0
        ),
        "gamma" = calculateNorm(
          global = TRUE,
          alpha = alphaGlobal,
          beta = betaGlobal,
          h = hGlobal,
          gamma = x,
          mu = muGlobal,
          psi = psiGlobal,
          delta = 0
        ),
        "mu" = calculateNorm(
          global = TRUE,
          alpha = alphaGlobal,
          beta = betaGlobal,
          h = hGlobal,
          gamma = gammaGlobal,
          mu = x,
          psi = psiGlobal,
          delta = 0
        ),
        "psi" = calculateNorm(
          global = TRUE,
          alpha = alphaGlobal,
          beta = betaGlobal,
          h = hGlobal,
          gamma = gammaGlobal,
          mu = muGlobal,
          psi = x,
          delta = 0
        )
      )
    } else{
      simulation <- switch(
        param,
        "alpha" = calculateNorm(
          global = FALSE,
          alpha = x,
          beta = betaLocal,
          h = hLocal,
          gamma = gammaLocal,
          mu = muLocal,
          psi = psiLocal,
          delta = deltaLocal
        ),
        "beta" = calculateNorm(
          global = FALSE,
          alpha = alphaLocal,
          beta = x,
          h = hLocal,
          gamma = gammaLocal,
          mu = muLocal,
          psi = psiLocal,
          delta = deltaLocal
        ),
        "h" = calculateNorm(
          global = FALSE,
          alpha = alphaLocal,
          beta = betaLocal,
          h = x,
          gamma = gammaLocal,
          mu = muLocal,
          psi = psiLocal,
          delta = deltaLocal
        ),
        "gamma" = calculateNorm(
          global = FALSE,
          alpha = alphaLocal,
          beta = betaLocal,
          h = hLocal,
          gamma = x,
          mu = muLocal,
          psi = psiLocal,
          delta = deltaLocal
        ),
        "mu" = calculateNorm(
          global = FALSE,
          alpha = alphaLocal,
          beta = betaLocal,
          h = hLocal,
          gamma = gammaLocal,
          mu = x,
          psi = psiLocal,
          delta = deltaLocal
        ),
        "psi" = calculateNorm(
          global = FALSE,
          alpha = alphaLocal,
          beta = betaLocal,
          h = hLocal,
          gamma = gammaLocal,
          mu = muLocal,
          psi = x,
          delta = deltaLocal
        ),
        "delta" = calculateNorm(
          global = FALSE,
          alpha = alphaLocal,
          beta = betaLocal,
          h = hLocal,
          gamma = gammaLocal,
          mu = muLocal,
          psi = psiLocal,
          delta = x
        )
      )
    }
    return(simulation)
  }
  
  switch(
    param,
    "alpha" = ABC_sequential(
      method = "Beaumont",
      model = model,
      tolerance_tab = tol,
      summary_stat_target = expected,
      nb_simul = n,
      prior = list(c("unif", 0, 1))
    ),
    "beta" = ABC_sequential(
      method = "Beaumont",
      model = model,
      tolerance_tab = tol,
      summary_stat_target = expected,
      nb_simul = n,
      prior = list(c("unif", 0, 1))
    ),
    "h" = ABC_sequential(
      method = "Beaumont",
      model = model,
      tolerance_tab = tol,
      summary_stat_target = expected,
      nb_simul = n,
      prior = list(c("unif", 0, 0.3))
    ),
    "gamma" = ABC_sequential(
      method = "Beaumont",
      model = model,
      tolerance_tab = tol,
      summary_stat_target = expected,
      nb_simul = n,
      prior = list(c("unif", 0, 1))
    ),
    "mu" = ABC_sequential(
      method = "Beaumont",
      model = model,
      tolerance_tab = tol,
      summary_stat_target = expected,
      nb_simul = n,
      prior = list(c("unif", 0, 1))
    ),
    "psi" = ABC_sequential(
      method = "Beaumont",
      model = model,
      tolerance_tab = tol,
      summary_stat_target = expected,
      nb_simul = n,
      prior = list(c("unif", 0, 20))
    ),
    "delta" = ABC_sequential(
      method = "Beaumont",
      model = model,
      tolerance_tab = tol,
      summary_stat_target = expected,
      nb_simul = n,
      prior = list(c("unif", 0, 1))
    )
  )
}
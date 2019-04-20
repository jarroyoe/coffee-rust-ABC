library(parallel)

#create connection of NetLogo to R
RNetlogoSetup <- function(dummy,
                          gui,
                          nl.path,
                          nl.jarname,
                          model.path) {
  library(RNetLogo)
  NLStart(nl.path = nl.path,
          gui = gui,
          nl.jarname = nl.jarname)
  NLLoadModel(model.path)
}

NetLogoSimulation <-
  function(dataSetPath) {
    #generate cluster and load with variables
    #NOTICE: NetLogo Connection.R has to be loaded first
    #in order for the Global Environment to feed each core
    load(dataSetPath)
    cores <- detectCores()
    cl <- makeCluster(cores)
    clusterExport(cl, ls())
    clusterEvalQ(cl, as.vector(lsf.str(.GlobalEnv)))
    
    #setup Netlogo
    parLapply(
      cl,
      1:cores,
      RNetlogoSetup,
      gui = FALSE,
      nl.path = nl.path,
      nl.jarname = nl.jarname,
      model.path = model.path
    )
    
    #procedure
    paramsGlobal <- c("alpha", "beta", "gamma", "psi", "h", "mu")
    paramsLocal <-
      c("alpha", "beta", "gamma", "psi", "h", "mu", "delta")
    
    #generate prior mean and variance for global and local control strategies
    testGlobal <-
      parsapply(
        cl,
        1:n,
        calculateNorm,
        TRUE,
        alphaGlobal,
        betaGlobal,
        hGlobal,
        gammaGlobal,
        muGlobal,
        psiGlobal,
        0
      )
    expectedGlobal <- mean(testGlobal)
    tolGlobal <- c(2 * sd(testGlobal), sd(testGlobal))
    testLocal <-
      parsapply(
        cl,
        1:n,
        calculateNorm,
        FALSE,
        alphaLocal,
        betaLocal,
        hLocal,
        gammaLocal,
        muLocal,
        psiLocal,
        deltaLocal
      )
    expectedLocal <- mean(testLocal)
    tolLocal <- c(2 * sd(testLocal), sd(testLocal))
    
    resultsGlobal <-
      parSapply(cl,
                paramsGlobal,
                runProcedure(),
                TRUE,
                expectedGlobal,
                tolGlobal)
    resultsLocal <-
      parSapply(cl,
                paramsLocal,
                runProcedure(),
                FALSE,
                expectedLocal,
                tolLocal)
    
    #close Netlogo and report results
    parLapply(cl, RNetLogoQuit, 1:cores)
    return(list(Global = resultsGlobal, Local = resultsLocal))
  }

RNetLogoQuit <- function(dummy) {
  NLQuit()
}
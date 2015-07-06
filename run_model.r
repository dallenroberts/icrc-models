##################################################
## Allen Roberts
## July 1, 2015
##################################################

rm(list = ls())

library(data.table)

## Global variables
year_start <- 1970
year_end <- 2020
tstep <- 0.25 # years
nsteps <- (year_end - year_start) * (1 / tstep)

## Attribute values
hiv <- c(0, 1)
age <- seq(1, 12)
male <- c(0, 1)
risk <- seq(1, 3)
cd4 <- seq(0, 5) ## For now CD4/VL = 0 means not applicable (susceptible or on treatment)
vl <- seq(0, 5)
circ <- c(0, 1)
prep <- c(0, 1)
condom <- c(0, 1)

## Variable for all of the dimensions in our model.  Usefully for switching the data.table keys
all_keys <- c("hiv", "age", "male", "risk", "cd4", "vl", "circ", "prep", "condom")

## Source functions
source("interpolate.r")
source("demography_functions.r")
source("progressDisease.r")
source("seedInfections.r")
source("riskAdjust.r")
# source("lambda_functions.r)

## Load input epidemiological parameters
source("load_parameters.r")

## Initialize population matrix
pop <- as.data.table(expand.grid(hiv, age, male, risk, cd4, vl, circ, prep, condom))
setattr(pop, 'names', c("hiv", "age", "male", "risk", "cd4", "vl", "circ", "prep", "condom"))
pop$count <- 0
pop$diff <- 0

## Run model
for(tt in 1:nsteps) {
  
  ## Seed infections for first year
  if(tt == 1) {
    ## Add intial populations.  Initially all are susceptible. 
    setkey(pop, hiv, age, male, cd4, vl, circ, prep, condom)
    pop[init_pop, count := pop]
    
    ## Distribute by risk group
    setkey(pop, age, male, risk)
    pop[risk_props, count := count * prop]
    
    ## Seed infections - this is currently adding 0.1% of total population to infected groups, but not subtracting them from the susceptible pool.  Need to confirm with Roger
    seedInfections(pop, 0.001) 
  }
  
  ## Calculate calendar year
  year <- floor(year_start + (tt - 1) * tstep)
  
  ## Demography
  addBirths(pop)
  subtractDeaths(pop)
  agePop(pop)
  
  ## Disease progression
  progressDisease(pop)
  
  ## Transmission
  # calcMixMat(pop) ## Sets up the mixing matrix
  # calcLambda(pop)
  
  
  # Compute end-of-year population and set difference back to zero for next iteration of loop
  pop[, c("count", "diff") := list(count + diff, 0)]
  
  # Adjust population to match risk prevalence
  riskAdjust(pop)
  
  ## Calculate some statistics
  
  # Increment time step
  tt <- tt + 1
  
}

## Make some plots


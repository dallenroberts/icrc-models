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

all_keys <- c("hiv", "age", "male", "risk", "cd4", "vl", "circ", "prep", "condom")

## Load parameters
## Initial population

init_pop <- fread("data/initial_populations.csv")
init_pop[, pop := floor(pop * 0.8)] ## Remnant from Roger's model
init_pop[, c("hiv", "cd4", "vl", "circ", "prep", "condom") := 0]
setkey(init_pop, hiv, age, male, cd4, vl, circ, prep, condom)

## Proportion of population in each risk group (by age)
risk_props <- fread("data/risk_proportions.csv")
setkey(risk_props, age, male, risk)

## Fertility
fert <- fread("data/base_fertility_rate.csv")

## Add effect modification by CD4 count
fert <- fert[, .(age, male, gamma, cd4 = rep(0:5, each = 12))]

## Add these coefficients - need to confirm this with Roger because the models and the supplementary are contradictory
fert_coeffs <- data.table(cd4 = seq(0, 5), coeff = c(1, 1, 0.59, 0.59, 0.42, 0.42))
setkey(fert_coeffs, cd4)
setkey(fert, cd4)
fert[fert_coeffs, gamma := gamma * coeff]

## Background mortality (non-HIV) by age and sex
back_mort <- fread("data/background_mortality.csv")
setkey(back_mort, age, male)

## HIV-relative mortality
hiv_mort <- fread("data/hiv_mortality.csv")
hiv_mort$hiv <- 1
setkey(hiv_mort, hiv, age, cd4)

## Disease progression - cd4_duration and vl_duration are average duration spent in that CD4 or VL category (respectively) in years
dis_prog <- fread("data/disease_progression.csv")
dis_prog$hiv <- 1

## Initialize population matrix
pop <- as.data.table(expand.grid(hiv, age, male, risk, cd4, vl, circ, prep, condom))
setattr(pop, 'names', c("hiv", "age", "male", "risk", "cd4", "vl", "circ", "prep", "condom"))
pop$count <- 0
pop$diff <- 0

## Demography functions
addBirths <- function(dt) {
  
  ## Parameters (move these outside)
  nncirc_prop <- 0.1 ## Neonatal circumcision prevalence
  
  ## Setting vertical transmission. Note that year is a global variable.  Should provide a linear interpolation once the hard-coded values can be verified
  if(year <= 2004) {
    vert_trans <- 0.34
  } else if(year > 2004 & year < 2008) {
    vert_trans <- 0.202
  } else {
    vert_trans <- 0.071
  }
  setkey(fert, age, male, cd4)
  setkey(dt, age, male, cd4)
  
  ## All births
  dt[fert, births := count * gamma]
  
  ## Calculate births from uninfected mothers
  births_from_neg <- dt[hiv == 0, sum(births, na.rm = TRUE)]
  
  ## Calculate  births from infected mothers
  births_from_pos <- dt[hiv == 1, sum(births, na.rm = TRUE)]
  
  ## Calculate number of HIV+ births
  pos_births <- births_from_pos * vert_trans
  
  ## Calculate number of HIV- births
  neg_births <- births_from_pos * (1 - vert_trans) + births_from_neg
  
  ## Add to population
  dt[, births := 0]
  
  ## Add to population by sex
  dt[hiv == 0 & age == 1 & vl == 0 & cd4 == 0 & prep == 0 & condom == 0, births := neg_births * 0.5]
  dt[hiv == 1  & age == 1 & vl == 1 & cd4 == 1 & prep == 0 & condom == 0, births := pos_births * 0.5]
  
  ## Distribute by circumcision
  dt[male == 0 & circ == 1, births := 0]
  dt[male == 1 & circ == 1, births := births * nncirc_prop]
  dt[male == 1 & circ == 0, births := births * (1 - nncirc_prop)]
  
  ## Distribute totals across risk status
  setkey(dt, age, male, risk)
  dt[risk_props, births := births * prop]
  
  ## Add births to population
  dt[, diff := diff + births]
  dt[, births := NULL]
  
}

subtractDeaths <- function(dt) {
  
  ## Subtract non-HIV deaths
  setkey(dt, age, male)
  dt[back_mort, diff := diff - count * mu]
  
  ## Subtract HIV deaths
  setkey(dt, hiv, age, cd4)
  dt[hiv_mort, diff := diff - count * alpha]
  
}



agePop <- function(dt) {

  ## Efflux - subtract 1/5 of each compartment
  dt[, diff := diff - count * 1/5]
  
  # Influx - for all but first age group, add 1/5 of corresponding compartment from previous age group
  prev_age <- copy(dt)[, age := age + 1] ## This is making a copy but still might be faster than using dcast to get the populations of the previous age group
  setnames(prev_age, "count", "count_prev")
  setkeyv(dt, all_keys)
  dt[prev_age, diff := diff + 1/5 * count_prev]

}

## Disease progression
progressDisease <- function(dt) {
  
  # dt <- copy(pop)
  
  ## Note that we want the probability of progressing in a given time step.  This can be obtained from the mean duration using the exponential decay function (assuming a constant rate).  http://hyperphysics.phy-astr.gsu.edu/hbase/nuclear/meanlif.html
  
  setkey(dis_prog, hiv, male, vl, cd4)
  setkey(dt, hiv, male, vl, cd4)
  
  ## CD4 Efflux
  dt[dis_prog[cd4 < 5], diff := diff - count * exp(1 - 1/cd4_duration)]
  l
  ## VL Efflux
  dt[dis_prog[vl < 5], diff := diff - count * exp(1 - 1/vl_duration)]
  
  ## Influx - need to get previous categories
  ## CD4 Influx
  prev_cd4 <- copy(dt[cd4 > 0 & cd4 < 5 & hiv == 1])
  prev_cd4[dis_prog[cd4 < 5], c("prob", "count_prev") := list(exp(1 - 1/cd4_duration), count)]
  prev_cd4[, cd4 := cd4 + 1]
  setkeyv(prev_cd4, all_keys)
  setkeyv(dt, all_keys)
  dt[prev_cd4, diff := diff + prob * count_prev]
  
  ## VL Influx (check that above code actually works first)
  
}

interpolate <- function(breaks, values, step_size, start = year_start, end = year_end) {
  
  smoothed <- rep(0, (year_end - year_start)/ step_size + 1)
  
  ## Checks
  if(!is.numeric(breaks) | !is.numeric(values) | !is.numeric(step_size)) stop("Error: breaks, values, and steps must both be numeric vectors")
  if(any(is.na(breaks) | is.na(values))) stop("Error: Missing breaks or values" )
  if(min(breaks) < year_start | max(breaks) > year_end) stop("Error: breaks fall outside of year limits")
  if(any(sort(breaks) != breaks)) stop("Error: breaks are not sorted")
  
  ## Starting values
  if(breaks[1] != year_start) {
    smoothed[1:((breaks[1] - year_start)/step_size + 1)] <- values[1]
  }
  
  ## Ending values
  if(breaks[length(breaks)] != year_end) {
    smoothed[((breaks[length(breaks)] - year_start) / step_size + 1):length(smoothed)] <- values[length(values)] 
    
  }
  
  ## Intermediate values
  for(ii in 1:length(values)) {
    smoothed[(breaks[ii] - year_start)/step_size + 1] <- values[ii]
  }
  
  ## Interpolated values
  for(ii in 1:(length(values) - 1)) {
    
    smoothed[(breaks[ii] - year_start)/step_size + 1] <- values[ii]
    
    slope <- (values[ii + 1] - values[ii])/((breaks[ii + 1] - breaks[ii])/step_size + 1)
    smoothed[((breaks[ii] - year_start)/step_size + 2):((breaks[ii + 1] - year_start)/step_size )] <- values[ii] + slope * seq(1, ((breaks[ii + 1] - breaks[ii])/step_size - 1))
      
    
  }
  
}
## Add intial populations.  Initially all are susceptible. 
setkey(pop, hiv, age, male, cd4, vl, circ, prep, condom)
pop[init_pop, count := pop]

## Distribute by risk group
setkey(pop, age, male, risk)
pop[risk_props, count := count * prop]

## Seed infections - this is currently adding 0.1% of total population to infected groups, but not subtracting them from the susceptible pool.  Need to confirm with Roger.

seedInfections <- function(dt, prop) {
  
  ## Requires variable naming conventions
  ## prop is the proportion of the population that will be added to the infectious groups
  
  ## Total population
  N <- sum(dt$count)
  
  ## Distribute evenly among males and females
  dt[hiv == 1 & risk == 2 & cd4 == 1 & vl == 1 & circ == 0 & prep == 0 & condom == 0 & ((male == 1 & age == 6 ) | (male == 0 & age == 5)), count := N * prop / 2]

}

riskAdjust <- function(dt) {
  
  ## Sums the population across risk categories for each compartment
  dt[, c("risk", "sum") := list(risk, sum(count)), by = .(hiv, age, male, cd4, vl, circ, prep, condom)]
  
  ## Multiplies the summed population by the risk proportions defined in the initial parameters
  setkey(dt, age, male, risk)
  dt[risk_props, count := sum * prop]
  
  dt[, sum := NULL]
}

## Run model
for(tt in 1:nsteps) {
  
  ## Seed infections for first year
  if(tt == 1) {
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





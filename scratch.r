##################################################
## Allen Roberts
## July 1, 2015
##################################################

rm(list = ls())

library(data.table)


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
if(t == 0) {
  seedInfections(pop, 0.001) 
  pop[hiv == 0, sum(count)] * 0.001 == pop[hiv == 1, sum(count)]
}

# addBirths(pop)
subtractDeaths(pop)
agePop(pop)




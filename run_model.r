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
nsteps <- (year_end - year_start) / tstep + 1

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
source("lambda_functions.r")
source("transmit.r")

## Load input epidemiological parameters
source("load_parameters.r")

## Initialize population matrix
pop <- as.data.table(expand.grid(sapply(all_keys, get)))
setattr(pop, 'names', c("hiv", "age", "male", "risk", "cd4", "vl", "circ", "prep", "condom"))
pop$count <- 0
pop$diff <- 0

## Data tables for statistics
population <- as.data.table(expand.grid(hiv, age, male, risk, cd4, vl, circ, prep, condom, seq(year_start, year_end)))
setattr(population, 'names', c("hiv", "age", "male", "risk", "cd4", "vl", "circ", "prep", "condom", "yy"))
population[, pop_size := 0]


# births <- as.data.table(expand.grid(age, male, seq(1, nsteps)))
# setattr(births, 'names', c("age", "male", "time"))
# 
deaths <- as.data.table(expand.grid(hiv, age, male, risk, seq(1, nsteps)))
setattr(deaths, 'names', c("hiv", "age", "male", "risk", "time"))
deaths[, c("non_aids_deaths", "aids_deaths") := 0]

incidence <- as.data.table(expand.grid(age, male, risk, seq(1, nsteps)))
setattr(incidence, 'names', c("age", "male", "risk", "time"))
incidence[, c("vert_infections", "horiz_infections", "rate") := 0]

# art <- c(0, 1)
# interventions <- as.data.table(expand.grid(age, male, circ, prep, condom, art, seq(1, nsteps)))

## Run model
for(tt in 1:nsteps) {
  
  # tt <- 1
  print(tt)

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
  progressDisease(pop, tstep)
  
  ## Transmission
  ## Calculate the mixing matrix
  calcMixMat(pop, mixing_matrix, tt)
  
  ## Calculate adjusted partnerships per year
  adjustPartnerships(pop, mixing_matrix) 
  
  ## Calculate lambda
  calcLambda(pop, mixing_matrix, adjusted_partners)
  
  ## Transmit infections
  transmit(pop, lambda_mat)
  
  # Compute end-of-year population and set difference back to zero for next iteration of loop
  pop[, c("count", "diff") := list(count + diff, 0)]
  
  # Adjust population to match risk prevalence
  riskAdjust(pop)
  
  ## Calculate some statistics
  ## Populations
  if((year_start + (tt - 1) * tstep) %% 1 == 0) {
    setkeyv(population, c("yy", all_keys))
    setkeyv(pop, c("yy", all_keys))
    population[pop, pop_size := count]
  }
  
  # Increment time step
  tt <- tt + 1
  
}

## Make some plots
population_size <- population[, list(total_size = sum(pop_size)), by = list(male, age, yy)]
prevalence <- population[hiv == 1, list(size = sum(pop_size)), by = list(male, age, yy)]
setkey(population_size, male, age, yy)
setkey(prevalence, male, age, yy)
prevalence[population_size, prev := size/total_size]
prevalence[is.na(prev), prev := 0]

risks <- population[, list(size = sum(pop_size)), by = list(risk, yy)]

## Population size
pop_plot <- ggplot(population_size, aes(x = yy, y = total_size)) +
    geom_line(aes(colour = factor(male))) +
    facet_wrap(~age)

## HIV prevalence
prev_plot <- ggplot(prevalence, aes(x = yy, y = prev)) +
  geom_line(aes(group = male, colour = factor(male))) +
  facet_wrap(~age) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))

## HIV incidence
incidence[, c("year", "year_exact") := list(floor(year_start + (time - 1) * tstep), floor(year_start + (time - 1) * tstep))]
incidence[, infections := vert_infections + horiz_infections]

inc_rate_plot <- ggplot(data = incidence, aes(x = year_exact, y = rate * 100)) +
    geom_line(aes(colour = factor(male), linetype = factor(risk))) +
    facet_wrap(~age)

infections <- incidence[, list(yearly_infections = sum(infections)), by = list(age, male, year)]
inc_plot <- ggplot(data = infections, aes(x = year, y = yearly_infections)) +
    geom_line(aes(colour = factor(male))) +
    facet_wrap(~age)

## Deaths
deaths[is.na(non_aids_deaths), non_aids_deaths := 0]
deaths[is.na(aids_deaths), aids_deaths := 0]

deaths[, c("year", "year_exact") := list(floor(year_start + (time - 1) * tstep), floor(year_start + (time - 1) * tstep))]
deaths_data <- deaths[, list(background = sum(non_aids_deaths), hiv_relative = sum(aids_deaths), total = sum(non_aids_deaths + aids_deaths)), by = list(hiv, age, male, year)]

deaths_data <- melt(deaths_data, id_vars = c("hiv", "age", "male", "year"), measure.vars = c("background", "hiv_relative", "total"), variable.name = "death_type", value.name = "deaths")
deaths_plot <- ggplot(data = deaths_data, aes(x = year, y = deaths)) +
    geom_line(aes(colour = factor(male), linetype = death_type)) +
    facet_wrap(~age + hiv)

hiv_deaths_data <- deaths_data[hiv == 1 & death_type == "hiv_relative"]
hiv_deaths_plot <- ggplot(data = hiv_deaths_data, aes(x = year, y = deaths)) +
  geom_line(aes(colour = factor(male))) +
  facet_wrap(~age)

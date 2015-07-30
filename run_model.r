##################################################
## Allen Roberts
## July 1, 2015
##################################################

rm(list = ls())

library(data.table)
library(ggplot2)
library(reshape2)

## Global variables
year_start <- 1970
year_end <- 2020
tstep <- 0.1 # years
nsteps <- (year_end - year_start) / tstep + 1

## Attribute values
hiv <- c(0, 1)
age <- seq(1, 12)
male <- c(0, 1)
risk <- seq(1, 3)
cd4 <- seq(0, 5) ## CD4/VL = 0 right now means that they are uninfected
vl <- seq(0, 5)
circ <- c(0, 1)
prep <- c(0, 1)
condom <- c(0, 1)
art <- c(0, 1)

## Variable for all of the dimensions in our model.  Usefully for switching the data.table keys
all_keys <- c("hiv", "age", "male", "risk", "cd4", "vl", "circ", "prep", "condom", "art")

## Source functions
source("interpolate.r")
source("demography_functions.r")
source("progressDisease.r")
source("seedInfections.r")
source("riskAdjust.r")
source("lambda_functions.r")
source("transmit.r")
source("distributeART.r")
source("distributeCondoms.r")

## Load input epidemiological parameters
source("load_parameters.r")

## Initialize population matrix
pop <- as.data.table(expand.grid(sapply(all_keys, get)))
setattr(pop, 'names', c("hiv", "age", "male", "risk", "cd4", "vl", "circ", "prep", "condom", "art"))
pop$count <- 0
pop$diff <- 0

## Data tables for statistics
population <- as.data.table(expand.grid(hiv, age, male, risk, cd4, vl, circ, prep, condom, art, seq(year_start, year_end)))
setattr(population, 'names', c("hiv", "age", "male", "risk", "cd4", "vl", "circ", "prep", "condom", "art", "yy"))
population[, pop_size := 0]


# births <- as.data.table(expand.grid(age, male, seq(1, nsteps)))
# setattr(births, 'names', c("age", "male", "time"))
# 
deaths <- as.data.table(expand.grid(hiv, age, male, risk, seq(1, nsteps)))
setattr(deaths, 'names', c("hiv", "age", "male", "risk", "time"))
deaths[, c("non_aids_deaths", "aids_deaths") := 0]

incidence <- as.data.table(expand.grid(age, male, risk, seq(1, nsteps)))
setattr(incidence, 'names', c("age", "male", "risk", "time"))
incidence[, c("vert_infections", "horiz_infections") := 0]

## Distribution of CD4 and VL
dis_dist <- as.data.table(expand.grid(age, male, cd4, vl, seq(1, nsteps)))
setattr(dis_dist, 'names', c("age", "male", "cd4", "vl", "time"))

# art <- c(0, 1)
# interventions <- as.data.table(expand.grid(age, male, circ, prep, condom, art, seq(1, nsteps)))

## Run model
for(tt in 1:nsteps) {
  
  # tt <- 1
  print(tt)

  if(tt == 1) {
    ## Add intial populations.  Initially all are susceptible. 
    setkey(pop, hiv, age, male, cd4, vl, circ, prep, condom, art)
    pop[init_pop, count := pop]
    
    ## Distribute by risk group
    setkey(pop, age, male, risk)
    pop[risk_props, count := count * prop]
    
    ## Seed infections - this is currently adding 0.1% of total population to infected groups, but not subtracting them from the susceptible pool.  Need to confirm with Roger
    seedInfections(pop, 0.001) 
  }
  
  ## Calculate calendar year
  year <- floor(year_start + (tt - 1) * tstep)
  
  ## Distribute ART coverage
  distributeART(pop, tt)
  
  ## Distribute condom coverage
  distributeCondoms(pop, tt)
  
  ## Optional reduction in background mortality
  if(year >= 1990) {
    back_mort[, mu := mu * (100 - 2 * tstep)/100]
  }
  
  ## Optional decrease in fertility
  # fert[, gamma := gamma * (100 - tstep)/ 100]
  
  ## Demography
  addBirths(pop)
  subtractDeaths(pop)
  agePop(pop, tstep)
  pop[, yy := year]
  
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

  # ## HIV prevalence
  # print(year)
  # print(pop[age == 5, sum(count)/sum(pop$count[pop$age == 5]), by = hiv])
  
}

## Make some plots
population_size <- population[, list(total_size = sum(pop_size)), by = list(male, age, yy)]
age_prevalence <- population[hiv == 1, list(size = sum(pop_size)), by = list(male, age, yy)]
setkey(population_size, male, age, yy)
setkey(age_prevalence, male, age, yy)
age_prevalence[population_size, prev := size/total_size]
age_prevalence[is.na(prev), prev := 0]

risks <- population[, list(size = sum(pop_size)), by = list(risk, yy)]

## Population size
kzn_pop <- fread("data/kzn_population.csv")
zaf_pop <- fread("data/zaf_population.csv")
pop_plot <- ggplot(data = population_size, aes(x = yy, y = total_size)) +
    geom_line(aes(colour = factor(male))) +
    geom_point(data = kzn_pop, aes(x = year, y = pop, colour = factor(male))) +
    facet_wrap(~age)

## HIV prevalence
total_prev <- fread("data/total_prevalence.csv")
age_prev <- fread("data/age_specific_prevalence.csv")
prev_plot <- ggplot(age_prevalence, aes(x = yy, y = prev)) +
  geom_line(aes(group = male, colour = factor(male))) +
  geom_point(data = age_prev, aes(x = year, y = prevalence, colour = factor(male))) +
  facet_wrap(~age) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))

overall_prev <- age_prevalence[age > 3 & age <= 10, list(total_prev = sum(size) / sum(size / prev)), by = list(yy)]
total_prev_plot <- ggplot(overall_prev, aes(x = yy, y = total_prev)) +
  geom_line() +
  geom_point(data = total_prev, aes(x = year, y = prevalence, colour = factor(location))) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))

## HIV incidence
incidence[, c("year", "year_exact") := list(floor(year_start + (time - 1) * tstep), floor(year_start + (time - 1) * tstep))]
incidence[, infections := vert_infections + horiz_infections]

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

## Distribution of cd4 and vl
cd4_dist <- population[hiv == 1, list(size = sum(pop_size)), by = list(age, cd4, yy)]
cd4_dist[, pct := size/sum(size), by = list(age, yy)]
vl_dist <- population[hiv == 1, list(size = sum(pop_size)), by = list(age, vl, yy)]
vl_dist[, pct := size/sum(size), by = list(age, yy)]

cd4_plot <- ggplot(data = cd4_dist, aes(x = yy, y = pct, group = factor(cd4))) +
  geom_line(aes(x = yy, y = pct, colour = factor(cd4), position = 'stack')) +
  facet_wrap(~age)

vl_plot <- ggplot(data = vl_dist, aes(x = yy, y = pct, group = factor(vl))) +
  geom_line(aes(x = yy, y = pct, colour = factor(vl), position = 'stack')) +
  facet_wrap(~age)


##################################################
## Allen Roberts
## July 7, 2015
## Description: Load demography functions for compartmental model
##################################################

## This file contains the demography functions for the compartmental model.  "addBirths" augments the diff variable by multiplying fertility rates by counts of women. "subtractDeaths" decreases the diff variable by multiplying deaths rates by counts.  "agePop" keeps track (in "diff") of the counts entering/leaving each age group due to aging.

## Note that "dt" stands for "data.table" and "time_index" is the iteration of the loop (corresponding to the global variable tt), which represents the discrete time point. Functions that have "time_index" as an argument are time-dependent.

## Demography functions
addBirths <- function(dt, time_index = tt) {
  
  ## Parameters (move these outside)
  nncirc_prop <- 0.1 ## Neonatal circumcision prevalence
  
  setkey(fert, age, male, cd4, art)
  setkey(dt, age, male, cd4, art)
  
  ## All births
  dt[fert, births := count * gamma]
  
  ## Keep track of birth statistics
  birth_stats <- dt[male == 0, list(time = tt, num = sum(births)), by = list(hiv, age)]
  setkey(birth_stats, time, hiv, age)
  setkey(births, time, hiv, age)
  births[birth_stats, num_births := num]
  
  ## Calculate births from uninfected mothers. Count mothers on ART as "negatives"
  births_from_neg <- dt[hiv == 0 | art == 1, sum(births, na.rm = TRUE)]
  
  ## Calculate  births from infected mothers
  births_from_pos <- dt[hiv == 1 & art == 0, sum(births, na.rm = TRUE)]
  
  ## Calculate number of HIV+ births
  pos_births <- births_from_pos * vert_trans[time_index]
  
  ## Calculate number of HIV- births
  neg_births <- births_from_pos * (1 - vert_trans[time_index]) + births_from_neg
  
  ## Initialize births added
  dt[, births := 0]
  
  ## Distribute births added by sex
  dt[hiv == 0 & age == 1 & vl == 0 & cd4 == 0 & prep == 0 & condom == 0 & art == 0, births := neg_births * 0.5]
  dt[hiv == 1  & age == 1 & vl == 1 & cd4 == 1 & prep == 0 & condom == 0 & art == 0, births := pos_births * 0.5]
  
  ## Distribute births added by circumcision
  dt[male == 0 & circ == 1, births := 0]
  dt[male == 1 & circ == 1, births := births * nncirc_prop]
  dt[male == 1 & circ == 0, births := births * (1 - nncirc_prop)]
  
  ## Distribute total births added across risk status
  setkey(dt, age, male, risk)
  dt[risk_props, births := births * prop]
  
  ## Add births to population
  dt[, diff := diff + births]
  
  ## Keep track of new infections
  hiv_births <- dt[hiv == 1 & age == 1, list(inf_births = sum(births), time = tt), by = list(age, male)]
  
  setkey(hiv_births, time, age, male)
  setkey(incidence, time, age, male)
  
  incidence[hiv_births, vert_infections := inf_births]
  
  ## Clean up
  dt[, births := NULL]

  
}

subtractDeaths <- function(dt) {
  
  ## Subtract non-HIV deaths
  setkey(dt, age, male)
  dt[back_mort, back_deaths := count * mu]
  
  ## Subtract HIV deaths
  setkey(dt, hiv, age, cd4, art)
  dt[hiv_mort, hiv_deaths := count * alpha]
  dt[hiv == 0, hiv_deaths := 0]
  
  dt[, diff := diff - back_deaths - hiv_deaths]
  
  ## Keep track
  death_stats <- dt[, list(aids_deaths = sum(hiv_deaths), non_aids_deaths = sum(back_deaths), time = tt), by = list(hiv, age, male)]
  setkey(death_stats, time, hiv, age, male)
  setkey(deaths, time, hiv, age, male)
  deaths[death_stats, c("hiv_deaths", "back_deaths") := list(aids_deaths, non_aids_deaths)]
  
  ## Clean up
  dt[, c("back_deaths", "hiv_deaths") := NULL]
#   non_hiv_deaths <- dt[back_mort, list(hiv = hiv, age = age, male = male, risk = risk, deaths_count = count * mu)]
#   non_hiv_deaths <- non_hiv_deaths[, list(total_deaths = sum(deaths_count), time = tt), by = list(hiv, age, male, risk)]
#   setkey(non_hiv_deaths, time, hiv, age, male, risk)
#   setkey(deaths, time, hiv, age, male, risk)
#   deaths[non_hiv_deaths, non_aids_deaths := total_deaths]
#   
#   ## Keep track
#   hiv_deaths <- dt[hiv_mort, list(hiv = hiv, age = age, male = male, risk = risk, deaths_count = count * alpha)]
#   hiv_deaths <- hiv_deaths[, list(total_deaths = sum(deaths_count), time = tt), by = list(hiv, age, male, risk)]
#   setkey(hiv_deaths, time, hiv, age, male, risk)
#   setkey(deaths, time, hiv, age, male, risk)
#   deaths[hiv_deaths, aids_deaths := total_deaths]
}


agePop <- function(dt, time_step) {
  
  ## Efflux - subtract 1/5 of each compartment for each year
  dt[, diff := diff - count * time_step/5]
  
  # Influx - for all but first age group, add 1/5 of corresponding compartment from previous age group
  prev_age <- copy(dt)[, age := age + 1] ## This is making a copy but still might be faster than using dcast to get the populations of the previous age group
  setnames(prev_age, "count", "count_prev")
  setkeyv(dt, all_keys)
  dt[prev_age, diff := diff + time_step/5 * count_prev]
  
}

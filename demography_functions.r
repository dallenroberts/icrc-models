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
  
  setkey(fert, age, male, cd4)
  setkey(dt, age, male, cd4)
  
  ## All births
  dt[fert, births := count * gamma]
  
  ## Calculate births from uninfected mothers
  births_from_neg <- dt[hiv == 0, sum(births, na.rm = TRUE)]
  
  ## Calculate  births from infected mothers
  births_from_pos <- dt[hiv == 1, sum(births, na.rm = TRUE)]
  
  ## Calculate number of HIV+ births
  pos_births <- births_from_pos * vert_trans[time_index]
  
  ## Calculate number of HIV- births
  neg_births <- births_from_pos * (1 - vert_trans[time_index]) + births_from_neg
  
  ## Initialize births added
  dt[, births := 0]
  
  ## Distribute births added by sex
  dt[hiv == 0 & age == 1 & vl == 0 & cd4 == 0 & prep == 0 & condom == 0, births := neg_births * 0.5]
  dt[hiv == 1  & age == 1 & vl == 1 & cd4 == 1 & prep == 0 & condom == 0, births := pos_births * 0.5]
  
  ## Distribute births added by circumcision
  dt[male == 0 & circ == 1, births := 0]
  dt[male == 1 & circ == 1, births := births * nncirc_prop]
  dt[male == 1 & circ == 0, births := births * (1 - nncirc_prop)]
  
  ## Distribute total births added across risk status
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

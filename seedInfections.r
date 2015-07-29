##################################################
## Allen Roberts
## July 7, 2015
## Description: Function to seed initial HIV+ infections into the model
##################################################

## This function adds HIV+ infections to the initial population.  Currently it does not subtract those infections from the HIV- population.  "prop" is the proportion of the TOTAL population that will be infected.  These are distributed equally across the two different groups.

seedInfections <- function(dt, prop) {
  
  ## Requires variable naming conventions
  ## prop is the proportion of the population that will be added to the infectious group
  
  ## Total population
  N <- sum(dt$count)
  
  ## Distribute evenly among males and females
  dt[hiv == 1 & risk == 2 & cd4 == 1 & vl == 1 & circ == 0 & prep == 0 & condom == 0 & art == 0 & ((male == 1 & age == 6 ) | (male == 0 & age == 5)), count := N * prop / 2]
  
}

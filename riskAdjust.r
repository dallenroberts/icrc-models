##################################################
## Allen Roberts
## July 7, 2015
## Description: Function to redistribute population according to risk distribution specified in parameters.
##################################################

## This function adjusts the population such that the distribution across risk categories is constant at each time point.  The risk distribution is specified in the initial parameters.

riskAdjust <- function(dt) {
  
  ## Sums the population across risk categories for each compartment
  dt[, c("risk", "sum") := list(risk, sum(count)), by = .(hiv, age, male, cd4, vl, circ, prep, condom, art)]
  
  ## Multiplies the summed population by the risk proportions defined in the initial parameters
  setkey(dt, age, male, risk)
  dt[risk_props, count := sum * prop]
  
  dt[, sum := NULL]
}

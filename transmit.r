##################################################
## Allen Roberts
## July 15, 2015
## Description: Transmits infections based on lambda (risk of infection) and effect of interventions such as PrEP, circumcision, and condoms.
##################################################

transmit <- function(dt, lambdas) {
  
  ## Merge on risk reduction interventions.
  dt[, psi := 1]
  dt[condom == 1, psi := psi * (1 - risk_reduction[intervention == "condom", psi])]
  dt[circ == 1 & male == 1, psi := psi * (1 - risk_reduction[intervention == "circ", psi])]
  dt[prep == 1, psi := psi * (1 - risk_reduction[intervention == "prep", psi])]
  
  ## Merge on lambda
  setkey(dt, male, age, risk)
  setkey(lambdas, male, age, risk)
  dt[lambdas, lambda := lambda]
  
  ## Subtract from HIV negative population
  dt[hiv == 0, diff := diff - count * lambda * psi]

  ## Keep track of incidence
  new_infections <- dt[hiv == 0, list(time = tt, new_inf = sum(count * lambda * psi)), by = list(male, age)]
  setkey(new_infections, male, age, time)
  setkey(incidence, male, age, time)
  incidence[new_infections, horiz_infections := new_inf]
  
  ## Add to HIV positive population
  new_hiv <- copy(dt)
  new_hiv <- new_hiv[hiv == 0]
  new_hiv[, c("hiv", "new_infections") := list(1, count * lambda * psi)]
  
  ## Seed new HIV infections in vl = 1 and cd4 = 1
  new_hiv <- new_hiv[, list(new_infections = sum(new_infections)), by = list(hiv, age, male, risk, circ, prep, condom, art)]
  new_hiv[, c("cd4", "vl") := 1]
  setkeyv(new_hiv, all_keys)
  setkeyv(dt, all_keys)
  dt[new_hiv, diff := diff + new_infections]
  
  ## Clean up
  dt[, c("psi", "lambda") := NULL]
  
}
## Mixing matrix - sets the proportion of partnerships that come from each age and sexual risk group
calcMixMat <- function(dt) {
  
  dt <- copy(pop)
  ## Set epsilons - these parameters govern the extent to which mixing is random or assortative. Consider smoothing this parameter.
  epsilons <- fread("data/epsilons.csv")
  
  ## Set deltas - these parameters set the mixing distributions for completely assortative mixing.
  #   delta_risk <- ifelse(mix_mat$risk == mix_mat$risk_p, 1, 0)
  #   delta_age <- rep(0, nrow(mix_mat)) 
  #   
  #   ## Delta age is allowed to vary with time
  #   if(year <= 2004) {
  #     
  #     delta_age[mix_mat$age == mix_mat$age_p] <- 0.3
  #     delta_age[male == 0 & mix_mat$age == mix_mat$age_p + 1] <- 0.7
  #     delta_age[male == 1 & mix_mat$age == mix_mat$age_p - 1] <- 0.7
  #     
  #   } else {
  #     
  #     delta_age[mix_mat$age == mix_mat$age_p] <- 0.7
  #     delta_age[male == 0 & mix_mat$age == mix_mat$age_p + 1] <- 0.3
  #     delta_age[male == 1 & mix_mat$age == mix_mat$age_p - 1] <- 0.3
  #     
  #   }
  
  ## Calculate number of partners
  
  ## Calculate mixing matrix vector based off of random mixing
  dt[, sex_total := sum(count), by = list(male)]
  random_mix_risk <- dt[, list(prop = sum(count/sex_total)), by = list(male, risk)]
  random_mix_age <- dt[, list(prop = sum(count/sex_total)), by = list(male, age)]
  
  ## Calculate net mixing matrix as weighted average of assortative and random mixing
  
  
  
}

## Calculate per-act probability of HIV transmission by viral load of partner
baseline <- 0.0006 ## Baseline probability of HIV transmission (VL < 1000 copies/mL)
trans_probs <- fread("data/transmission_probabilities.csv")
trans_probs[, beta := baseline * scalar]

## Calculate risk reduction for HIV-negative partner based on intervention usage
risk_reduction <- fread("data/risk_reduction.csv")

## Number of coital acts per partnership
acts_per_year <- fread("data/acts_per_year.csv")

## Calculate per-partnership probability of HIV transmission per year.  Depends on risk group of HIV-negative partner and viral load of HIV positive partner.

## Calculate lambda (force of infection) for each individual. 

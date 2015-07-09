##################################################
## Allen Roberts
## July 7, 2015
## Description: Load functions to calculate force of infection for compartmental model
##################################################

require(reshape2)

## calcMixMat - sets the proportion of partnerships that come from each age and sexual risk group. It defines two scenarios - random mixing and assortative mixing. Assortative mixing is governed by delta (age and risk).  The weights governing the balance between assortative and random mixing are determined by epsilon. The mixing matrix contains, for each sex, age, and risk category, the proportion of partnerships that come from each sex, age, and risk category of the partner. The sum of the "prop" column in the mixing matrix by sex, age, and risk category should equal 1.

## adjustPartnerships - this function adjusts the annual number of sexual partnerships such that those reported by men and women are equal.

## Note that "dt" stands for "data.table" and "time_index" is the iteration of the loop (corresponding to the global variable tt), which represents the discrete time point. Functions that have "time_index" as an argument are time-dependent.

calcMixMat <- function(dt, mix_mat, time_index = tt) {
  
  ## Reset mixing matrix to zero
  mix_mat[, prop := 0]
  
  ## Set deltas - these parameters set the mixing distributions for completely assortative mixing. Eventually should move these parameters outside of the function.
  mix_mat$delta_risk <- ifelse(mix_mat$risk == mix_mat$risk_p, 1, 0)
  mix_mat$delta_age <- 0
  
  ## Delta age is allowed to vary with time
  if(year <= 2004) {
    
    mix_mat[age == age_p, delta_age := 0.3]
    mix_mat[male == 0 & age == age_p - 1, delta_age := 0.7]
    mix_mat[male == 1 & age == age_p + 1, delta_age := 0.7]
    
  } else {
    
    mix_mat[age == age_p, delta_age := 0.7]
    mix_mat[male == 0 & age == age_p - 1, delta_age := 0.3]
    mix_mat[male == 1 & age == age_p + 1, delta_age := 0.3]
    
  }
  
  ## Ensure that probabilities sum to one for youngest males/oldest females
  mix_mat[male == 0 & age == 12 & age_p == 12, delta_age := 1]
  mix_mat[male == 1 & age == 1 & age_p == 1, delta_age := 1]
  
  
  ## Calculate mixing matrix based off of random mixing
  ## First, merge on partnerships/year
  setkey(dt, age, male, risk)
  setkey(partners, age, male, risk)
 dt[partners, partners := partners] 
   ## This is really confusing because partners has three different environments here.  In order, I'm referring to the data.table in the parent frame, the new column in dt, and the column in the partners data.table. 
  
  ## Next, calculate distribution by age and sex assuming random mixing.
  dt[, sex_total := sum(partners * count), by = list(male)]
  random_mix_age <- dt[, list(prop_age = sum(partners * count/sex_total)), by = list(male, age)]
  setnames(random_mix_age, paste(names(random_mix_age), "p", sep = "_"))  ## Note that these proportions are used to find the probability of selecting the a partner with those characterstics. We add the "_p" here so it matches up with the mixing matrix notation.
  setkey(random_mix_age, male_p, age_p)

  ## Then, calculate distribution by risk given age (ie, the conditional probability of having a partnership from risk category r given the partner's age is a). Pr(a,r) = Pr(r|a) * Pr(a). Note that this really shouldn't change the way the model is currently structured, since the risk distribution is set automatically.  But this allows it to change should we choose to relax that assumption later.
  dt[, sex_age_total := sum(partners * count), by = list(male, age)]
  random_mix_risk_age <- dt[, list(prop_risk_by_age = sum(partners * count/sex_age_total)), by = list(male, risk, age)]
  setnames(random_mix_risk_age, paste(names(random_mix_risk_age), "p", sep = "_"))
  setkey(random_mix_risk_age, male_p, age_p, risk_p)
  
  ## Calculate net mixing matrix as weighted average of assortative and random mixing
  ## Calculate Pr(a), which is the first term in the product.
  setkey(mix_mat, male_p, age_p)
  mix_mat <- mix_mat[random_mix_age, pa := prop_age_p * epsilons[time_index] + (1 - epsilons[time_index]) * delta_age]
  
  ## Calcualte Pr(r|a), which is the second term in the product.
  setkey(mix_mat, male_p, age_p, risk_p)
  mix_mat <- mix_mat[random_mix_risk_age, pr_a:= prop_risk_by_age_p * epsilons[time_index] + (1 - epsilons[time_index]) * delta_risk]
  
  ## Multiply together to get the mixing matrix
  mix_mat[, prop := pa * pr_a]

  ## Clean up the workspace
  mix_mat[, c("pa", "pr_a", "delta_age", "delta_risk") := NULL]
  rm(random_mix_age, random_mix_risk_age)
  
  ## Checks - note that for Roger's paper the males in the youngest age group and females in the oldest age group have proportions less than 1.  This is because the delta assignment is based off of partnerships with older males and younger females.  Roger's paper assumed the effect was negligible since the age groups aren't as relevant for HIV.  We get around this with the last lines of code in the delta section.
  mix_mat[, sum(prop), by = list(male, age, risk)]
  
}

## Adjust number of partnerships. The number of partnerships is adjusted for differences in reported number of sexual partners between males and females. The global parameter "theta" determines whether the difference is driven by males or females, where theta = 1 indicates entirely male-driven and theta = 0 indicates entirely female-driven.
adjustPartnerships <- function(dt, mix_mat) {
  
  ## Sum number of partnerships by age, sex, and risk
  sums <- dt[, list(count = sum(count)), by = list(age, male, risk)]
  sums[partners, partners_count := count * partners]
  
  ## Multiply number of partnerships by mixing matrix
  setkey(mix_mat, male, age, risk)
  setkey(sums, male, age, risk)
  mix_mat[sums, partners := partners_count * prop]
  
  ## Calculate male and female discrepancy in reported partners
  disc <- dcast.data.table(mix_mat, age + risk + age_p + risk_p ~ male, value.var = "partners")
  setnames(disc, "0", "female")
  setnames(disc, "1", "male")
  disc[, discrepancy := male/female]
  disc[, c("female", "male") := NULL]
  disc <- disc[rep(1:nrow(disc), times = 2)]
  disc[, male := rep(c(0, 1), each = nrow(disc) / 2)]
  
  ## Calculate adjusted partnerships/year.  Note that the adjusted partnerships/year, unlike partnerships/year, seems to take into account the age and risk category of the partner.
  ## Merge on the number of partners (unadjusted)
  setkey(disc, male, age, risk)
  setkey(partners, male, age, risk)
  disc[partners, partners := partners]
  
  ## Calculated adjusted partners based on global variable theta
  disc[male == 1, adjusted_partners := partners * discrepancy ^ -(1 - theta)]
  disc[male == 0, adjusted_partners := partners * discrepancy ^ theta]
  
  return(disc)
}
  



## Calculate per-partnership probability of HIV transmission per year.  Depends on risk group of HIV-negative partner and viral load of HIV positive partner.

## Calculate lambda (force of infection) for each individual. 

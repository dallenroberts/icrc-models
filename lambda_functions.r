##################################################
## Allen Roberts
## July 7, 2015
## Description: Load functions to calculate force of infection for compartmental model
##################################################

## calcMixMat - sets the proportion of partnerships that come from each age and sexual risk group. It defines two scenarios - random mixing and assortative mixing. Assortative mixing is governed by delta (age and risk).  The weights governing the balance between assortative and random mixing are determined by epsilon. The mixing matrix contains, for each sex, age, and risk category, the proportion of partnerships that come from each sex, age, and risk category of the partner. The sum of the "prop" column in the mixing matrix by sex, age, and risk category should equal 1.

## adjustPartnerships - this function adjusts the annual number of sexual partnerships such that those reported by men and women are equal.

## calcLambda - this function calculates the force of infection, which is the probability of HIV acquisition for each HIV negative person in the population. To do so, we first must expand the mixing matrix by viral loads of the partners, since transmission probabilities are a function of the sex and risk status of the HIV-negative partner and the viral load of the HIV-positive partner. We then merge on the transmission risks (betas) to the expanded mixing matrix (hereafter referred to as "lambda_mat"). The expanded mixing matrix is collapsed by weighting the transmission risk by the probability of encountering an HIV-positive partner in that viral load category (note that this is relevant for vl = 0, which includes both HIV negative and HIV+ on ART) within each age/sex/risk partnership, which is just the proportion of the total population in that age/sex/risk category that is HIV + with viral load v'. This is the per-partnership per year risk of transmission (from the partner to the individual). We then calculate the per-individual per-year risk of transmission by including the number of partnerships that individual has. We then collapse across partnerships to get each individual's risk of infection. This does not include the effect of preventive interventions for the HIV-negative individual, such as circumcision, condom usage, or PrEP. 

## Note that "dt" stands for "data.table" and "time_index" is the iteration of the loop (corresponding to the global variable tt), which represents the discrete time point. Functions that have "time_index" as an argument are time-dependent.

calcMixMat <- function(dt, mix_mat, time_index = tt) {
  
  ## Reset mixing matrix to zero
  mix_mat[, prop := 0]
  
  ## Set deltas - these parameters set the mixing distributions for completely assortative mixing. Eventually should move these parameters outside of the function.
  mix_mat[, delta_risk := ifelse(risk == risk_p, 1, 0)]
  mix_mat[, delta_age := 0]
  
  ## Delta age is allowed to vary with time
  mix_mat[age == age_p, delta_age := deltas[time_index]]
  mix_mat[male == 0 & age == age_p - 1, delta_age := 1 - deltas[time_index]]
  mix_mat[male == 1 & age == age_p + 1, delta_age := 1 - deltas[time_index]]

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
  mix_mat[random_mix_age, pa := prop_age_p * epsilons[time_index] + (1 - epsilons[time_index]) * delta_age]
  
  ## Calcualte Pr(r|a), which is the second term in the product.
  setkey(mix_mat, male_p, age_p, risk_p)
  mix_mat[random_mix_risk_age, pr_a:= prop_risk_by_age_p * epsilons[time_index] + (1 - epsilons[time_index]) * delta_risk]
  
  ## Multiply together to get the mixing matrix
  mix_mat[, prop := pa * pr_a]

  ## Clean up the workspace
  dt[, c("partners", "sex_total", "sex_age_total") := NULL]
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
  
  ## Calculate male and female discrepancy in reported partners.  Could probably figure out a way to reshape to clean this up, but it is somewhat non-trivial.
  male_reports <- mix_mat[male == 1, .(age, risk, age_p, risk_p, partners)]
  setnames(male_reports, c("age", "risk", "age_p", "risk_p", "partners"), c("age_male", "risk_male", "age_female", "risk_female", "partners_male"))
  
  female_reports <- mix_mat[male == 0, .(age, risk, age_p, risk_p, partners)]
  setnames(female_reports, c("age", "risk", "age_p", "risk_p", "partners"), c("age_female", "risk_female", "age_male", "risk_male", "partners_female"))
  
  setkey(male_reports, age_male, risk_male, age_female, risk_female)
  setkey(female_reports, age_male, risk_male, age_female, risk_female)
  male_reports[female_reports, c("partners_female", "discrepancy") := list(partners_female, partners_male / partners_female)]
  female_reports[male_reports, c("partners_male", "discrepancy") := list(partners_male, partners_male / partners_female)]
  
  ## Reformat
  setnames(male_reports, c("age_male", "risk_male", "age_female", "risk_female"), c("age", "risk", "age_p", "risk_p"))
  male_reports[, c("male", "partners_male", "partners_female") := list(1, NULL, NULL)]

  setnames(female_reports, c("age_female", "risk_female", "age_male", "risk_male"), c("age", "risk", "age_p", "risk_p"))
  female_reports[, c("male", "partners_male", "partners_female") := list(0, NULL, NULL)]
  
  disc <- rbindlist(list(male_reports, female_reports))
  
  ## Calculate adjusted partnerships/year.  Note that the adjusted partnerships/year, unlike partnerships/year, seems to take into account the age and risk category of the partner.
  ## Merge on the number of partners (unadjusted)
  setkey(disc, male, age, risk)
  setkey(partners, male, age, risk)
  disc[partners, partners := partners]
  
  ## Calculated adjusted partners based on global variable theta
  disc[male == 1, adjusted_partners := partners * discrepancy ^ -(1 - theta)]
  disc[male == 0, adjusted_partners := partners * discrepancy ^ theta]
  
  ## Note that superassigment operator here - for now (for debugging) we want adjusted_partners in the parent environment.  We can figure out later a way to merge some of the lambda_functions together so we don't need to keep track of disc in the run_model environment.
  adjusted_partners <<- disc[, .(male, age, risk, age_p, risk_p, adjusted_partners)]
  
  ## Clean up
  mix_mat[, partners := NULL]
}

## Calculate lambda (force of infection) for each individual. 
calcLambda <- function(dt, mix_mat, adj_parts) {
  
  ## Multiply mixing matrix by adjusted partnership matrix to get number of partners per person per year in each possible partnership type
  setkey(mix_mat, male, age, risk, age_p, risk_p)
  setkey(adj_parts, male, age, risk, age_p, risk_p)
  mix_mat[adj_parts, adjusted_partners := adjusted_partners * prop]
  
  ## Expand mixing matrix by viral load of partner
  lambda_mat <- rbindlist(lapply(0:5, function(x, d) data.table(d, vl_p = x), d = mix_mat))
  ## Expand mixing matrix by ART status of partner
  lambda_mat <- rbindlist(lapply(0:1, function(x, d) data.table(d, art_p = x), d = lambda_mat))
  
  ## Merge on transmission probabilities for each partnership
  setkey(lambda_mat, male, risk, vl_p, art_p)
  setkey(betas, male, risk, vl_p, art_p)
  lambda_mat[betas, transmission_risk := transmission_risk]
  
  ## Calculate number HIV + people in each ART and viral load category by sex, age, and risk
  art_vl_prev <- dt[hiv == 1, .(art_vl_count = sum(count)), by = list(male, age, risk, vl, art)]
  
  ## Calculate total number of people in each sex, age, and risk category
  total_counts <- dt[, .(total = sum(count)), by = list(male, age, risk)]
  setkey(total_counts, male, age, risk)
  setkey(art_vl_prev, male, age, risk)
  art_vl_prev[total_counts, total := total]
  setnames(art_vl_prev, c("male", "age", "risk", "vl", "art"), c("male_p", "age_p", "risk_p", "vl_p", "art_p"))
  setkey(art_vl_prev, male_p, age_p, risk_p, vl_p, art_p)
  
  ## Merge on counts of people in each partnership/vl category to mixing matrix
  setkey(lambda_mat, male_p, age_p, risk_p, vl_p, art_p)
  lambda_mat[art_vl_prev, c("art_vl_count", "total") := list(art_vl_count, total)]
  
  ## Calculate per-partnership per year risk - weighted average of transmission risk based on counts of HIV+ in each viral load category in each partnership divided by total population (HIV+ and HIV-) in each age/sex/risk category
  lambda_mat <- lambda_mat[, list(pp_risk = sum(art_vl_count * transmission_risk / total), adjusted_partners = median(adjusted_partners)), by = list(male, age, risk, male_p, age_p, risk_p)]

  ## Multiply number of partners per person per year in each possible partnership type by the per-partnership per year transmission risk. Note that this formula differs from Roger's supplemental since we're explicitly using risks here. 
  lambda_mat[, total_risk := 1 - (1 - pp_risk) ^ adjusted_partners]
  
  ## Calculate the risk over all partnerships
  lambda_mat <<- lambda_mat[, list(lambda = 1 - prod(1 - total_risk)), by = list(male, age, risk)]

  ## Clean up
  mix_mat[, adjusted_partners := NULL]
  
}



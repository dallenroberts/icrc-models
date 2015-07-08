##################################################
## Allen Roberts
## July 7, 2015
## Description: Load model input parameters
##################################################

## This loads the epidemiological parameters that govern the compartmental model functions.  Ideally these numbers should be carefully referenced, with labeled units, and the hardcoded values should be moved out to csvs.

require(data.table)

init_pop <- fread("data/initial_populations.csv")
init_pop[, pop := floor(pop * 0.8)] ## Remnant from Roger's model
init_pop[, c("hiv", "cd4", "vl", "circ", "prep", "condom") := 0]
setkey(init_pop, hiv, age, male, cd4, vl, circ, prep, condom)

## Proportion of population in each risk group (by age)
risk_props <- fread("data/risk_proportions.csv")
setkey(risk_props, age, male, risk)

## Fertility
fert <- fread("data/base_fertility_rate.csv")

## Add effect modification by CD4 count
fert <- fert[, .(age, male, gamma, cd4 = rep(0:5, each = 12))]

## Add these fertility coefficients - need to confirm this with Roger because the models and the supplementary are contradictory
fert_coeffs <- data.table(cd4 = seq(0, 5), coeff = c(1, 1, 0.59, 0.59, 0.42, 0.42))
setkey(fert_coeffs, cd4)
setkey(fert, cd4)
fert[fert_coeffs, gamma := gamma * coeff]
rm(fert_coeffs)

## Vertical transmission
vert_trans <- fread("data/vertical_transmission.csv")
vert_trans <- interpolate(breaks = vert_trans$year, values = vert_trans$vert)

## Background mortality (non-HIV) by age and sex
back_mort <- fread("data/background_mortality.csv")
setkey(back_mort, age, male)

## HIV-relative mortality
hiv_mort <- fread("data/hiv_mortality.csv")
hiv_mort$hiv <- 1
setkey(hiv_mort, hiv, age, cd4)

## Disease progression - cd4_duration and vl_duration are average duration spent in that CD4 or VL category (respectively) in years
dis_prog <- fread("data/disease_progression.csv")
dis_prog$hiv <- 1

## Initialize mixing matrix.
mix_mat <- as.data.table(expand.grid(male, age, risk, male, age, risk))

## The "_p" indicates that the attribute corresponds to the individual's partner
setattr(mix_mat, 'names', c("male", "age", "risk", "male_p", "age_p", "risk_p"))

## Assuming only heterosexual transmission
mix_mat <- mix_mat[mix_mat$male != mix_mat$male_p]

mix_mat$prop <- 0

## Set epsilons - these parameters govern the extent to which mixing is random or assortative. Note that we could break these out separately by for epsilon_age and epsilon_risk.
epsilons <- fread("data/epsilons.csv")
epsilons <- interpolate(breaks= epsilons$year, values = epsilons$epsilon)

## Number of partners per year by age, sex, and risk
partners <- fread("data/partners_per_year.csv")

## Theta - parameter that governs the extent to which differences in reported number of sexual partners between males and females is male (1) or female (0) driven
theta <- 0.5

## Number of coital acts per partnership
acts <- fread("data/acts_per_year.csv")

## Per-act probability of HIV transmission by viral load of partner
baseline <- 0.0006 ## Baseline probability of HIV transmission (VL < 1000 copies/mL)
trans_probs <- fread("data/transmission_probabilities.csv")
trans_probs[, beta := baseline * scalar]

## Risk reduction for HIV-negative partner based on intervention usage
risk_reduction <- fread("data/risk_reduction.csv")




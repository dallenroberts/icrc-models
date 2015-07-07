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

## Vertical transmission
vert_trans <- fread("data/vertical_transmission.csv")
vert_trans <- interpolate(breaks = vert_trans$year, values = vert_trans$vert)

## Add effect modification by CD4 count
fert <- fert[, .(age, male, gamma, cd4 = rep(0:5, each = 12))]

## Add these fertility coefficients - need to confirm this with Roger because the models and the supplementary are contradictory
fert_coeffs <- data.table(cd4 = seq(0, 5), coeff = c(1, 1, 0.59, 0.59, 0.42, 0.42))
setkey(fert_coeffs, cd4)
setkey(fert, cd4)
fert[fert_coeffs, gamma := gamma * coeff]
rm(fert_coeffs)

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

## Initialize mixing matrix - can we do this just once or does it need to be sex specific?.  Male is only entered in once because all partnerships are assumed to be heterosexual, so male only identifies the individual (the partner is abs(male - 1))
mix_mat <- as.data.table(expand.grid(male, age, risk, age, risk))
## The "_p" indicates that the attribute corresponds to the individual's partner
setattr(mix_mat, 'names', c("male", "age", "risk", "age_p", "risk_p"))
mix_mat$prop <- 0

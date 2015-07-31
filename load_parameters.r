##################################################
## Allen Roberts
## July 7, 2015
## Description: Load model input parameters
##################################################

## This loads the epidemiological parameters that govern the compartmental model functions.  Ideally these numbers should be carefully referenced, with labeled units, and the hardcoded values should be moved out to csvs.

require(data.table)

init_pop <- fread("data/initial_populations.csv")
init_pop[, pop := floor(pop * 0.8)] ## Remnant from Roger's model - I think this is saying that the population in 1970 should be about 80% of that in 1985
init_pop[, c("hiv", "cd4", "vl", "circ", "prep", "condom", "art") := 0]
setkey(init_pop, hiv, age, male, cd4, vl, circ, prep, condom, art)

## ART coverage
art_coverage <- fread("data/art_coverage.csv")

art_cov <- lapply(sort(unique(art_coverage$cd4)), function(cd4_ind) {
  
  df <- art_coverage[cd4 == cd4_ind]
  interpolate(breaks = df$year, values = df$prop_art)
})
rm(art_coverage)

## Condom usage
condom_coverage <- fread("data/condom_usage.csv")
# condom_coverage[, usage := usage * 0.73] ## Optional bias correction
condom_coverage <- condom_coverage[order(condom_coverage$year), ]
condom_cov <- lapply(sort(unique(condom_coverage$age)), function(age_cat) {
  
  df <- condom_coverage[age == age_cat]
  interpolate(breaks = df$year, values = df$usage)
})
rm(condom_coverage)

## Proportion of population in each risk group (by age)
risk_props <- fread("data/risk_proportions.csv")
setkey(risk_props, age, male, risk)

## Fertility
fert <- fread("data/base_fertility_rate.csv")
fert[, gamma := gamma * 1.1]

## Add effect modification by CD4 count
fert <- fert[, .(art, age, male, gamma, cd4 = rep(0:5, each = 12))]

## Add these fertility coefficients
fert_coeffs <- data.table(cd4 = seq(0, 5), coeff = c(1, 1, 0.59, 0.59, 0.42, 0.42))
fert_coeffs <- rbindlist(lapply(0:1, function(x, d) data.table(d, art = x), d = fert_coeffs))
fert_coeffs[art == 1, coeff := 1] ## No reduction in fertility for mothers on ART
setkey(fert_coeffs, cd4, art)
setkey(fert, cd4, art)
fert[fert_coeffs, gamma := gamma * coeff]
rm(fert_coeffs)

## Adjust fertility by time step and convert to risk
fert[, gamma := 1 - exp(-gamma * tstep)]

## Vertical transmission
vert_trans <- fread("data/vertical_transmission.csv")
vert_trans <- interpolate(breaks = vert_trans$year, values = vert_trans$vert)

## Background mortality (non-HIV) by age and sex. These are per-year mortality rates.
back_mort <- fread("data/background_mortality.csv")
back_mort[, mu := 1 - exp(-mu * tstep)] # Adjust for time step and convert to risk
setkey(back_mort, age, male)

## HIV-relative mortality
hiv_mort <- fread("data/hiv_mortality.csv")
hiv_mort <- rbindlist(lapply(0:1, function(x, d) data.table(d, art = x), d = hiv_mort))
hiv_mort[art == 1, alpha := 0]
hiv_mort[, c("hiv", "alpha") := list(1, 1 - exp(-alpha * tstep))]
setkey(hiv_mort, hiv, age, cd4, art)

## Disease progression - cd4_duration and vl_duration are average duration spent in that CD4 or VL category (respectively) in years
dis_prog <- fread("data/disease_progression.csv")
dis_prog$hiv <- 1
dis_prog$art <- 0

## Initialize mixing matrix.
mixing_matrix <- as.data.table(expand.grid(male, age, risk, male, age, risk))

## The "_p" indicates that the attribute corresponds to the individual's partner
setattr(mixing_matrix, 'names', c("male", "age", "risk", "male_p", "age_p", "risk_p"))

## Assuming only heterosexual transmission
mixing_matrix <- mixing_matrix[mixing_matrix$male != mixing_matrix$male_p]

mixing_matrix$prop <- 0

## Set deltas - these parameters govern the mixing pattern for "completely assortative" mixing, which isn't truly "completely assortative" but is allowed to vary. The delta that gets loaded here is the proportion of partnerships that are with the same age group. 1-delta is the proportion that are with the age group one above (females) or below (males).
deltas <- fread("data/deltas.csv")
deltas <- interpolate(breaks = deltas$year, values = deltas$delta)

## Set epsilons - these parameters govern the extent to which mixing is random or assortative. Note that we could break these out separately by for epsilon_age and epsilon_risk.
epsilons <- fread("data/epsilons.csv")
epsilons <- interpolate(breaks= epsilons$year, values = epsilons$epsilon)

## Number of partners per year by age, sex, and risk
partners <- fread("data/partners_per_year.csv")
## Adjust by time-step
partners[, partners := partners * tstep]

# Test adjustment by factor
# partners[age <= 6, partners := partners * 1.05]
# partners[age == 7, partners := partners * 0.6]
# partners[age > 7, partners := partners * 0.4]

# ## Test redistribution of partners by age
# partners[age < 6, partners := partners * 1.2]
# partners[age > 6, partners := partners * 0.8]


## Theta - parameter that governs the extent to which differences in reported number of sexual partners between males and females is male (1) or female (0) driven
theta <- 0.5

## Number of coital acts per partnership 
acts <- fread("data/acts_per_partnership.csv")
## For testing purposes - reduce number of acts by some factor to make prevalence data seem more reasonable
# acts[, acts := acts/5]

## Per-act probability of HIV transmission by viral load of partner. Note that ART risk reduction is built in here. We might want to consider moving it 
baseline <- 0.0006 ## Baseline probability of HIV transmission (VL < 1000 copies/mL)
trans_probs <- fread("data/transmission_probabilities.csv")
trans_probs[, chi := baseline * scalar]

## Calculate per-partnership probability of HIV transmission per year.  Depends on risk group of HIV-negative partner and viral load of HIV positive partner.
betas <- data.table(expand.grid(male, risk, vl, art))
setattr(betas, 'names', c("male", "risk", "vl", "art"))

## Join transmission probabilities per act
setkey(betas, vl)
setkey(trans_probs, vl)
betas[trans_probs, chi := chi]

## Join sexual acts per year per partnership
setkey(acts, male, risk)
setkey(betas, male, risk)
betas[acts, acts := acts]

## Note that the vl and art variables actually correspond to the partner
setnames(betas, c("vl", "art"), c("vl_p", "art_p"))

## Calculate transmission probability per partnership per time step
betas[, transmission_risk := 1 - (1 - chi) ^ acts]

## Risk reduction for ART usage in HIV positive partner
betas[art_p == 1, transmission_risk := transmission_risk * 0.08]


## Risk reduction for HIV-negative partner based on intervention usage
risk_reduction <- fread("data/risk_reduction.csv")

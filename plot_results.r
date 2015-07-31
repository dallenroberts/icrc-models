##################################################
## Allen Roberts
## July 30, 2015
##################################################

load("output/2015-07-31/defaults.RData")

library(ggplot2)

## Formatting
population[, c("year", "year_exact") := list(floor(year_start + (time - 1) * tstep), year_start + (time - 1) * tstep)] 
population$male <- factor(population$male, levels = c(0, 1), labels = c("Female", "Male"))
population$age <- factor(population$age, levels = seq(1, 12), labels = c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59"))
population$hiv <- factor(population$hiv, levels = c(0, 1), labels = c("Negative", "Positive"))

births[, c("year", "year_exact") := list(floor(year_start + (time - 1) * tstep), year_start + (time - 1) * tstep)] 
births$age <- factor(births$age, levels = seq(1, 12), labels = c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59"))
births$hiv <- factor(births$hiv, levels = c(0, 1), labels = c("Negative", "Positive"))

deaths[, c("year", "year_exact") := list(floor(year_start + (time - 1) * tstep), year_start + (time - 1) * tstep)] 
deaths$male <- factor(deaths$male, levels = c(0, 1), labels = c("Female", "Male"))
deaths$age <- factor(deaths$age, levels = seq(1, 12), labels = c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59"))

incidence[, c("year", "year_exact") := list(floor(year_start + (time - 1) * tstep), year_start + (time - 1) * tstep)] 
incidence$male <- factor(incidence$male, levels = c(0, 1), labels = c("Female", "Male"))
incidence$age <- factor(incidence$age, levels = seq(1, 12), labels = c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59"))


interventions[, c("year", "year_exact") := list(floor(year_start + (time - 1) * tstep), year_start + (time - 1) * tstep)] 
interventions$male <- factor(interventions$male, levels = c(0, 1), labels = c("Female", "Male"))
interventions$age <- factor(interventions$age, levels = seq(1, 12), labels = c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59"))
interventions$hiv <- factor(interventions$hiv, levels = c(0, 1), labels = c("Negative", "Positive"))

## Plot control
theme_set(theme_bw())
colors <- c("blue4", "green4")
names(colors) <- c("Female", "Male")
sexColors <- scale_colour_manual(name = "Sex", values = colors)

## Population data
pop_data <- fread("data/kzn_population.csv")
pop_data$male <- factor(pop_data$male, levels = c(0, 1), labels = c("Female", "Male"))
pop_data$age <- factor(pop_data$age, levels = seq(1, 12), labels = c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59"))

## Total population
total_pop <- population[, list(total = sum(pop_size)), by= list(male, year_exact, year)]
total_pop_data <- pop_data[, list(total = sum(pop)), by = list(male, year)]

total_pop_plot <- ggplot(data = total_pop, aes(x = year_exact, y = total / 1000000)) +
  geom_line(aes(colour = male)) +
  geom_point(data = total_pop_data, aes(x = year, y = total / 1000000, colour = male)) +
  sexColors + 
  scale_x_continuous(limits = c(1970, 2020), breaks = seq(1970, 2020, by = 10)) +
  xlab("Year") + ylab("KZN Population (Millions)")


## Age-specific population
age_pop <- population[, list(total = sum(pop_size)), by= list(age, male, year_exact, year)]

age_pop_plot <- ggplot(data = age_pop, aes(x = year_exact, y = total / 1000000)) +
  geom_line(aes(colour = male)) +
  geom_point(data = pop_data, aes(x = year, y = pop / 1000000, colour = male)) +
  sexColors + 
  scale_x_continuous(limits = c(1970, 2020), breaks = seq(1970, 2020, by = 10)) +
  xlab("Year") + ylab("KZN Population (Millions)") +
  facet_wrap(~age)

## Birth rates
births <- births[, list(births = sum(num_births)), by = list(age, year_exact, year)]
births$male <- "Female"

setkey(births, age, year_exact, male)
setkey(age_pop, age, year_exact, male)

births[age_pop, denom := total]
births <- births[, list(birth_rate = sum(births)/sum(denom * tstep)), by = list(year, age, male)]

birth_rate_age_plot <- ggplot(data = births, aes(x = year, y = birth_rate)) +
    geom_line() +
    scale_x_continuous(limits = c(1970, 2020), breaks = seq(1970, 2020, by = 10)) +
    xlab("Year") + ylab("Annual birth rate") +
    facet_wrap(~age)

## Death rates
deaths <- deaths[, list(non_aids_deaths = sum(back_deaths), aids_deaths = sum(hiv_deaths)), by = list(age, male, year_exact, year)]

setkey(deaths, age, year_exact, male)
setkey(age_pop, age, year_exact, male)

deaths[age_pop, denom := total]
deaths <- deaths[, list(hiv_death_rate = sum(aids_deaths)/sum(denom * tstep), back_death_rate = sum(non_aids_deaths)/sum(denom * tstep)), by = list(year, age, male)]

deaths <- melt(deaths, id.vars = c("age", "year", "male"), measure.vars = c("hiv_death_rate", "back_death_rate"), variable.name = "type", value.name = "rate")

death_rate_age_plot <- ggplot(data = deaths, aes(x = year, y = rate)) +
  geom_line(aes(colour = male, linetype = type)) +
  sexColors + 
  scale_x_continuous(limits = c(1970, 2020), breaks = seq(1970, 2020, by = 10)) +
  xlab("Year") + ylab("Annual death rate") +
  facet_wrap(~age)

## HIV prevalence 15-49
total_adult_pop <- population[as.numeric(age) > 3 & as.numeric(age) <= 10, list(total = sum(pop_size)), by= list(male, year_exact)] 
total_prev <- population[as.numeric(age) > 3 & as.numeric(age) <= 10 , list(hiv_total = sum(pop_size)), by = list(hiv, male, year_exact)]
setkey(total_adult_pop, male, year_exact)
setkey(total_prev, male, year_exact)
total_prev[total_adult_pop, prev := 100 * hiv_total/total]
total_prev <- total_prev[hiv == "Positive"]

total_prev_data <- fread("data/total_prevalence.csv")
total_prev_data[, prev := prevalence * 100]

total_prev_plot <- ggplot(data = total_prev, aes(x = year_exact, y = prev)) +
    geom_line(aes(colour = male)) +
    geom_point(data = total_prev_data, aes(x = year, y = prev)) +
    sexColors + 
    scale_x_continuous(limits = c(1970, 2020), breaks = seq(1970, 2020, by = 10)) +
    xlab("Year") + ylab("HIV Prevalence (%)") +
    ggtitle("HIV Prevalence (15-49)")

## Age-specific prevalence
age_prev <- population[, list(hiv_total = sum(pop_size)), by = list(hiv, male, age, year_exact)]
setkey(age_prev, male, age, year_exact)
setkey(age_pop, male, age, year_exact)
age_prev[age_pop, prev := 100 * hiv_total/total]
age_prev <- age_prev[hiv == "Positive"]

age_prev_data <- fread("data/age_specific_prevalence.csv")
age_prev_data <- age_prev_data[male  < 2]
age_prev_data$male <- factor(age_prev_data$male, levels = c(0, 1), labels = c("Female", "Male"))
age_prev_data$age <- factor(age_prev_data$age, levels = seq(1, 12), labels = c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59"))
age_prev_data[, prev := prevalence * 100]

age_prev_plot <- ggplot(data = age_prev, aes(x = year_exact, y = prev)) +
  geom_line(aes(colour = male)) +
  geom_point(data = age_prev_data, aes(x = year, y = prev, colour = male)) +
  sexColors + 
  scale_x_continuous(limits = c(1970, 2020), breaks = seq(1970, 2020, by = 10)) +
  xlab("Year") + ylab("HIV Prevalence (%)") +
  facet_wrap(~age)

## HIV Incidence
## Total
## By Age
incidence <- incidence[, list(infections = sum(horiz_infections + vert_infections)), by = list(age, male, year_exact, year)]

setkey(incidence, age, year_exact, male)
setkey(age_pop, age, year_exact, male)

incidence[age_pop, denom := total]
incidence <- incidence[, list(incidence_rate = sum(infections)/sum(denom * tstep)), by = list(year, age, male)]

incidence_rate_age_plot <- ggplot(data = incidence, aes(x = year, y = incidence_rate * 100)) +
  geom_line(aes(colour = male)) +
  sexColors + 
  scale_x_continuous(limits = c(1970, 2020), breaks = seq(1970, 2020, by = 10)) +
  xlab("Year") + ylab("Incidence rate (%)") +
  facet_wrap(~age)

## Intervention coverage
## ART
art_age <- interventions[hiv == "Positive", list(size = sum(total)), by = list(age, male, year_exact, art)]
art_age <- art_age[art == 1]
hiv_age_pop <- population[hiv == "Positive", list(denom = sum(pop_size)), by = list(age, male, year_exact)]
setkey(art_age, age, male, year_exact)
setkey(hiv_age_pop, age, male, year_exact)
art_age[hiv_age_pop, cov := 100 * size / denom]
art_age[is.na(cov), cov := 0]

art_age_plot <- ggplot(data = art_age, aes(x = year_exact, y = cov)) +
  geom_line(aes(colour = male)) +
  sexColors + 
  scale_x_continuous(limits = c(1970, 2020), breaks = seq(1970, 2020, by = 10)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  xlab("Year") + ylab("ART Coverage (%)") +
  facet_wrap(~age)

## Condoms
condom_age <- interventions[, list(size = sum(total)), by = list(age, male, year_exact, condom)]
condom_age <- condom_age[condom == 1]
denom_age_pop <- population[, list(denom = sum(pop_size)), by = list(age, male, year_exact)]
setkey(condom_age, age, male, year_exact)
setkey(denom_age_pop, age, male, year_exact)
condom_age[denom_age_pop, cov := 100 * size / denom]
condom_age[is.na(cov), cov := 0]

condom_age_plot <- ggplot(data = condom_age, aes(x = year_exact, y = cov)) +
  geom_line(aes(colour = male)) +
  sexColors + 
  scale_x_continuous(limits = c(1970, 2020), breaks = seq(1970, 2020, by = 10)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  xlab("Year") + ylab("Condom Coverage (%)") +
  facet_wrap(~age)


## Circumcision
circ_age <- interventions[male == "Male", list(size = sum(total)), by = list(age, year_exact, circ)]
circ_age <- circ_age[circ == 1]
male_age_pop <- population[male == "Male", list(denom = sum(pop_size)), by = list(age, year_exact)]
setkey(circ_age, age, year_exact)
setkey(male_age_pop, age,  year_exact)
circ_age[male_age_pop, cov := 100 * size / denom]
circ_age[is.na(cov), cov := 0]

circ_age_plot <- ggplot(data = circ_age, aes(x = year_exact, y = cov)) +
  geom_line() +
  scale_x_continuous(limits = c(1970, 2020), breaks = seq(1970, 2020, by = 10)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  xlab("Year") + ylab("Circumcision Coverage (%)") +
  facet_wrap(~age)


# ## HIV incidence
# incidence[, c("year", "year_exact") := list(floor(year_start + (time - 1) * tstep), floor(year_start + (time - 1) * tstep))]
# incidence[, infections := vert_infections + horiz_infections]
# 
# infections <- incidence[, list(yearly_infections = sum(infections)), by = list(age, male, year)]
# inc_plot <- ggplot(data = infections, aes(x = year, y = yearly_infections)) +
#   geom_line(aes(colour = factor(male))) +
#   facet_wrap(~age)
# 

# 
# ## Distribution of cd4 and vl
# cd4_dist <- population[hiv == 1, list(size = sum(pop_size)), by = list(age, cd4, yy)]
# cd4_dist[, pct := size/sum(size), by = list(age, yy)]
# vl_dist <- population[hiv == 1, list(size = sum(pop_size)), by = list(age, vl, yy)]
# vl_dist[, pct := size/sum(size), by = list(age, yy)]
# 
# cd4_plot <- ggplot(data = cd4_dist, aes(x = yy, y = pct, group = factor(cd4))) +
#   geom_line(aes(x = yy, y = pct, colour = factor(cd4), position = 'stack')) +
#   facet_wrap(~age)
# 
# vl_plot <- ggplot(data = vl_dist, aes(x = yy, y = pct, group = factor(vl))) +
#   geom_line(aes(x = yy, y = pct, colour = factor(vl), position = 'stack')) +
#   facet_wrap(~age)
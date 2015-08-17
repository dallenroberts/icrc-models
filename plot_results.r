##################################################
## Allen Roberts
## July 30, 2015
##################################################

library(reshape2)
library(data.table)
library(ggplot2)
library(RColorBrewer)

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

stopifnot(dis_dist[cd4 == 0 | vl == 0, sum(total)] == 0) ## No HIV-infected individual should have VL or CD4 = 0 - that's reserved for uninfected (hiv == 0) individuals
dis_dist <- dis_dist[cd4 != 0 & vl != 0]
dis_dist[, c("year", "year_exact") := list(floor(year_start + (time - 1) * tstep), year_start + (time - 1) * tstep)]
dis_dist$male <- factor(dis_dist$male, levels = c(0, 1), labels = c("Female", "Male"))
dis_dist$age <- factor(dis_dist$age, levels = seq(1, 12), labels = c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59"))
dis_dist$cd4 <- factor(dis_dist$cd4, levels = seq(1, 5), labels = c("Acute", "> 500", "350-500", "200-349", "< 200"))
dis_dist$vl <- factor(dis_dist$vl, levels = seq(1, 5), labels = c("Acute", "< 1,000", "1,000-10,000", "10,000-50,000", "> 50,000"))


## Plot control
theme_set(theme_bw())
colors <- c("blue4", "green4", "red4")
names(colors) <- c("Female", "Male", "Both")
sexColors <- scale_colour_manual(name = "Sex", values = colors)

vl_colors <- rev(brewer.pal(6, "Spectral"))
names(vl_colors) <- c("On ART", "Acute", "< 1,000", "1,000-10,000", "10,000-50,000", "> 50,000")
vlColors <- scale_fill_manual(name = "Viral Load", values = vl_colors)

cd4_colors <- rev(brewer.pal(6, "Spectral"))
names(cd4_colors) <- c("On ART", "Acute", "> 500", "350-500", "200-349", "< 200")
cd4Colors <- scale_fill_manual(name = "CD4", values = cd4_colors)


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
  scale_x_continuous(limits = c(year_start, year_end), breaks = seq(year_start, year_end, by = 10)) +
  xlab("Year") + ylab("KZN Population (Millions)") +
  ggtitle("Total population")


## Age-specific population
age_pop <- population[, list(total = sum(pop_size)), by= list(age, male, year_exact, year)]

age_pop_plot <- ggplot(data = age_pop, aes(x = year_exact, y = total / 1000000)) +
  geom_line(aes(colour = male)) +
  geom_point(data = pop_data, aes(x = year, y = pop / 1000000, colour = male)) +
  sexColors + 
  scale_x_continuous(limits = c(year_start, year_end), breaks = seq(year_start, year_end, by = 10)) +
  xlab("Year") + ylab("KZN Population (Millions)") +
  facet_wrap(~age) +
  ggtitle("Age-specific population")

## Birth rates
births <- births[, list(births = sum(num_births)), by = list(age, year_exact, year)]
births$male <- "Female"

setkey(births, age, year_exact, male)
setkey(age_pop, age, year_exact, male)

births[age_pop, denom := total]
births <- births[, list(birth_rate = sum(births)/sum(denom * tstep)), by = list(year, age, male)]

birth_rate_age_plot <- ggplot(data = births, aes(x = year, y = birth_rate)) +
    geom_line() +
    scale_x_continuous(limits = c(year_start, year_end), breaks = seq(year_start, year_end, by = 10)) +
    xlab("Year") + ylab("Annual birth rate") +
    facet_wrap(~age) +
  ggtitle("Age-specific fertility rates")

## Death rates
deaths <- deaths[, list(non_aids_deaths = sum(back_deaths), aids_deaths = sum(hiv_deaths)), by = list(age, male, year_exact, year)]

setkey(deaths, age, year_exact, male)
setkey(age_pop, age, year_exact, male)

deaths[age_pop, denom := total]
deaths <- deaths[, list(hiv_death_rate = sum(aids_deaths)/sum(denom * tstep), back_death_rate = sum(non_aids_deaths)/sum(denom * tstep)), by = list(year, age, male)]

deaths <- melt(deaths, id.vars = c("age", "year", "male"), measure.vars = c("hiv_death_rate", "back_death_rate"), variable.name = "type", value.name = "rate")
deaths$type <- factor(deaths$type, levels = c("hiv_death_rate", "back_death_rate"), labels = c("HIV", "Background"))
death_rate_age_plot <- ggplot(data = deaths, aes(x = year, y = rate)) +
  geom_line(aes(colour = male, linetype = type)) +
  sexColors + 
  scale_x_continuous(limits = c(year_start, year_end), breaks = seq(year_start, year_end, by = 10)) +
  xlab("Year") + ylab("Annual death rate") +
  facet_wrap(~age) +
  ggtitle("HIV and background age-specific mortality rates")

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
    geom_point(data = total_prev_data, aes(x = year, y = prev, shape = location)) +
    sexColors + 
    scale_x_continuous(limits = c(year_start, year_end), breaks = seq(year_start, year_end, by = 10)) +
    xlab("Year") + ylab("HIV Prevalence (%)") +
    ggtitle("HIV Prevalence (15-49)")

## Age-specific prevalence
age_prev <- population[, list(hiv_total = sum(pop_size)), by = list(hiv, male, age, year_exact)]
setkey(age_prev, male, age, year_exact)
setkey(age_pop, male, age, year_exact)
age_prev[age_pop, prev := 100 * hiv_total/total]
age_prev <- age_prev[hiv == "Positive"]

age_prev_data <- fread("data/age_specific_prevalence.csv")
# age_prev_data <- age_prev_data[male  < 2]
age_prev_data$male <- factor(age_prev_data$male, levels = c(0, 1, 2), labels = c("Female", "Male", "Both"))
age_prev_data$age <- factor(age_prev_data$age, levels = seq(1, 12), labels = c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59"))
age_prev_data[, prev := prevalence * 100]

age_prev_plot <- ggplot(data = age_prev, aes(x = year_exact, y = prev)) +
  geom_line(aes(colour = male)) +
  geom_point(data = age_prev_data, aes(x = year, y = prev, colour = male, shape = location)) +
  sexColors + 
  scale_x_continuous(limits = c(year_start, year_end), breaks = seq(year_start, year_end, by = 10)) +
  xlab("Year") + ylab("HIV Prevalence (%)") +
  facet_wrap(~age) +
  ggtitle("Age-specific HIV prevalence")

## HIV Incidence
## Total
total_incidence <- incidence[, list(infections = sum(horiz_infections + vert_infections)), by = list(year, year_exact, male)]
setkey(total_incidence, male, year_exact, year)
setkey(total_pop, male, year_exact, year)
total_incidence[total_pop, denom := total]
total_incidence <- total_incidence[, list(incidence_rate = sum(infections)/sum(denom * tstep)), by = list(year, male)]

incidence_rate_plot <- ggplot(data = total_incidence, aes(x = year, y = incidence_rate * 100)) +
  geom_line(aes(colour = male)) +
  sexColors + 
  scale_x_continuous(limits = c(year_start, year_end), breaks = seq(year_start, year_end, by = 10)) +
  xlab("Year") + ylab("Incidence rate (%)") +
  ggtitle("Total HIV incidence")

## Incidence in women ages 15-49 and men ages 15-54
adult_incidence <- incidence[as.numeric(age) > 3 & ((male == "Female" & as.numeric(age) <= 10) | male == "Male" & as.numeric(age) <= 11), list(infections = sum(horiz_infections + vert_infections)), by = list(year, year_exact, male)]
setkey(adult_incidence, male, year_exact)
setkey(total_adult_pop, male, year_exact)
adult_incidence[total_adult_pop, denom := total]
adult_incidence <- adult_incidence[, list(incidence_rate = sum(infections)/sum(denom * tstep)), by = list(year, male)]

adult_incidence_data <- fread("data/adult_incidence.csv")
adult_incidence_data$male <- factor(adult_incidence_data$male, levels = c(0, 1, 2), labels = c("Female", "Male", "Both"))
adult_incidence_data[, incidence := incidence * 100]



adult_incidence_rate_plot <- ggplot(data = adult_incidence, aes(x = year, y = incidence_rate * 100)) +
  geom_line(aes(colour = male)) +
  geom_point(data = adult_incidence_data, aes(x = year, y = incidence, colour = male, shape = location)) + 
  sexColors + 
  scale_x_continuous(limits = c(year_start, year_end), breaks = seq(year_start, year_end, by = 10)) +
  xlab("Year") + ylab("Incidence rate (%)") +
  ggtitle("Adult  HIV incidence")


## By Age
age_incidence <- incidence[, list(infections = sum(horiz_infections + vert_infections)), by = list(age, male, year_exact, year)]

setkey(age_incidence, age, year_exact, male)
setkey(age_pop, age, year_exact, male)

age_incidence[age_pop, denom := total]
age_incidence <- age_incidence[, list(incidence_rate = sum(infections)/sum(denom * tstep)), by = list(year, age, male)]

age_incidence_data <- fread("data/age_specific_incidence.csv")
age_incidence_data$male <- factor(age_incidence_data$male, levels = c(0, 1, 2), labels = c("Female", "Male", "Both"))
age_incidence_data$age <- factor(age_incidence_data$age, levels = seq(1, 12), labels = c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59"))
age_incidence_data[, incidence := incidence * 100]
incidence_rate_age_plot <- ggplot(data = age_incidence, aes(x = year, y = incidence_rate * 100)) +
  geom_line(aes(colour = male)) +
  geom_point(data = age_incidence_data, aes(x = year, y = incidence, colour = male, shape = location)) +
  sexColors + 
  scale_x_continuous(limits = c(year_start, year_end), breaks = seq(year_start, year_end, by = 10)) +
  xlab("Year") + ylab("Annual incidence (%)") +
  facet_wrap(~age) +
  ggtitle("Age-specific HIV incidence")

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
  scale_x_continuous(limits = c(year_start, year_end), breaks = seq(year_start, year_end, by = 10)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  xlab("Year") + ylab("ART Coverage (%)") +
  facet_wrap(~age) +
  ggtitle("Age-specific ART coverage")

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
  scale_x_continuous(limits = c(year_start, year_end), breaks = seq(year_start, year_end, by = 10)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  xlab("Year") + ylab("Condom Coverage (%)") +
  facet_wrap(~age) +
  ggtitle("Age-specific condom usage")


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
  scale_x_continuous(limits = c(year_start, year_end), breaks = seq(year_start, year_end, by = 10)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  xlab("Year") + ylab("Circumcision Coverage (%)") +
  facet_wrap(~age) +
  ggtitle("Age-specific circumcision coverage")

## Distribution of disease
cd4_dist <- dis_dist
cd4_dist[art == 1, cd4 := "On ART"]
cd4_dist <- cd4_dist[, list(size = sum(total)), by = list(age, cd4, year_exact)]
cd4_dist[, pct := size / sum(size), by = list(age, year_exact)]
vl_dist <- dis_dist
vl_dist[art == 1, vl := "On ART"]
vl_dist <- vl_dist[, list(size = sum(total)), by = list(age, vl, year_exact)]
vl_dist[, pct := size / sum(size), by = list(age, year_exact)]

cd4_plot <- ggplot(data = cd4_dist, aes(x = year_exact, y = pct * 100)) +
  geom_area(aes(x = year_exact, y = pct * 100, fill = cd4, position = 'stack')) +
  cd4Colors +
  scale_x_continuous(limits = c(year_start, year_end), breaks = seq(year_start, year_end, by = 10)) +
  xlab("Year") + ylab("Percentage") +
  facet_wrap(~age) + ggtitle("CD4 Distribution among HIV+")

vl_plot <- ggplot(data = vl_dist, aes(x = year_exact, y = pct * 100)) +
  geom_area(aes(fill = vl, position = 'stack')) +
  vlColors + 
  scale_x_continuous(limits = c(year_start, year_end), breaks = seq(year_start, year_end, by = 10)) +
  xlab("Year") + ylab("Percentage") +
  facet_wrap(~age) + ggtitle("Viral Load Distribution among HIV+")

## Save as PDF
pdf(file = paste0("output/", date, "/", name, ".pdf"), width = 10, height = 8)
print(total_pop_plot)
print(age_pop_plot)
print(birth_rate_age_plot)
print(death_rate_age_plot)
print(total_prev_plot)
print(age_prev_plot)
print(incidence_rate_plot)
print(adult_incidence_rate_plot)
print(incidence_rate_age_plot)
print(art_age_plot)
print(circ_age_plot)
print(condom_age_plot)
print(cd4_plot)
print(vl_plot)
dev.off()

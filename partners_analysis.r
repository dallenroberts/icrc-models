rm(list = ls())

library(ggplot2)
library(reshape2)
library(MASS)

partners <- read.csv("data/assa_partners_dist.csv", stringsAsFactors = FALSE)
partners$age <- seq(from = 12.5, to = 87.5, by = 5)

ggplot(data = partners, aes(x = age, y = male)) +
    geom_point()



icrc <- read.csv("data/partners_per_year.csv", stringsAsFactors = FALSE)
icrc$age <- icrc$age*5 - 2.5


ggplot(data = icrc, aes(x = age, y = partners)) +
    geom_line(aes(colour = factor(male), linetype = factor(risk)))

wide <- dcast(icrc, age + male ~ risk, value.var = 'partners')
wide$coeff1 <- wide[, "2"]/wide[, "1"]
wide$coeff2 <- wide[, "3"]/wide[, "2"]

assa_partners <- melt(partners, id.var = 'age', measure.vars = c('male', 'female'), variable.name = 'sex', value.name = 'assa_partners')
assa_partners$male <- ifelse(assa_partners$sex == "male", 1, 0)
assa_partners$risk <- 3

new_partners <- merge(assa_partners[, !names(assa_partners) %in% c("risk", "sex")], wide, by = c("age", "male"), all = TRUE)
new_partners$assa_partners[is.na(new_partners$assa_partners)] <- new_partners[is.na(new_partners$assa_partners), "3"]
new_partners[, "3"] <- new_partners$assa_partners * 5.5
new_partners[, "2"] <- new_partners[, "3"]/new_partners$coeff2
new_partners[, "1"] <- new_partners[, "2"]/new_partners$coeff1

new_partners <- melt(new_partners, id.vars = c("age", "male"), measure.vars = c("1", "2", "3"), variable.name = "risk", value.name = "partners")
new_partners <- new_partners[new_partners$age < 60, ]
ggplot(data = new_partners, aes(x = age, y = partners)) +
  geom_line(aes(colour = factor(male), linetype = factor(risk)))

new_partners$age <- floor((new_partners$age + 2.5)/5)
new_partners$risk <- as.numeric(as.character(new_partners$risk))
write.csv(new_partners, file = "data/new_partners_per_year.csv", row.names = FALSE)

## Male comparison
ggplot(data = icrc[icrc$risk == 3 & icrc$male == 1, ], aes(x = age, y = partners)) +
  geom_point(colour = 'blue') +
  geom_point(data = partners, aes(x = age, y = male * 6), colour = 'red')


## Female comparison
ggplot(data = icrc[icrc$risk == 3 & icrc$male == 0, ], aes(x = age, y = partners)) +
  geom_point(colour = 'blue') +
  geom_point(data = partners, aes(x = age, y = female * 5.5), colour = 'red')



gm <- glm(male ~ age, data = partners, family = Gamma)
partners$male_preds <- predict(gm)
ggplot(data = partners, aes(x = age, y = male)) +
  geom_point() +
  geom_line(aes(x = age, y = male_preds))

## Males
qplot(partners$age, partners$male)
male_fit <- fitdistr(partners$male, dgamma, start = list(shape = 2, rate = 0.5), lower = 0.0001)
curve(dweibull(x, scale=7.4, shape=3),from=0, to=100, main="Weibull
distribution")

test <- fitdistr(partners$male/850, 'gamma')
test <- fitdistr(partners$male, dweibull, start = list(scale = 1, shape = 1), lower = 0.01)


ggplot(data = partners, aes(x = age, y = male)) +
    geom_point() +
  geom_line(aes(x = age, y = 850*dgamma(partners$age, scale = 8, shape = 4.25)))


curve(dgamma(x, scale=10, shape=2),from=0, to=100, main="Gamma
distribution")


fitdistr(partners$male/100, 'gamma')
ages <- seq(1, 100, by = 1)


qplot(ages, dgamma(ages, shape = male_fit$estimate["shape"], rate = male_fit$estimate["rate"]))

## FEmales
qplot(partners$age, partners$female)
female_fit <- fitdistr(partners$female, "gamma")

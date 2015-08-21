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
new_partners[, "3"] <- new_partners$assa_partners * 5.6
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
  geom_point(data = partners, aes(x = age, y = male * 5.6), colour = 'red')


## Female comparison
ggplot(data = icrc[icrc$risk == 3 & icrc$male == 0, ], aes(x = age, y = partners)) +
  geom_point(colour = 'blue') +
  geom_point(data = partners, aes(x = age, y = female *5.6), colour = 'red')


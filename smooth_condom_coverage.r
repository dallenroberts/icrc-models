library(boot)
library(ggplot2)

rm(list = ls())
usage <- read.csv("data/condom_usage.csv", stringsAsFactors = FALSE)

ages <- seq(1, 12)
years <- seq(min(usage$year), max(usage$year))
# years <- seq(year_start, year_end)
lm <- loess(usage ~ year + factor(age), data = usage, degree = 1, span = 0.9)

pred_frame <- expand.grid(factor(ages), years)
names(pred_frame) <- c("age", "year")
preds <- as.data.frame(predict(lm, pred_frame))
names(preds) <- as.character(seq(min(usage$year), max(usage$year)))
preds$age <- seq(1, 12)
preds <- melt(preds, id.vars = c("age"), variable.name = "year", value.name = "usage")
preds$year <- as.numeric(as.character(preds$year))
preds$usage[preds$usage < 0] <- 0

ggplot(data = preds, aes(x = year, y = usage)) +
  geom_line() + 
  facet_wrap(~age)

write.csv(preds, file = "data/smoothed_condom_usage.csv", row.names = FALSE)

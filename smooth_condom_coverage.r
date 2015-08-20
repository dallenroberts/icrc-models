library(boot)
library(ggplot2)

rm(list = ls())
usage <- read.csv("data/condom_usage.csv", stringsAsFactors = FALSE)

ages <- seq(1, 12)
years <- seq(min(usage$year), max(usage$year))
# years <- seq(year_start, year_end)

count <- 1
for(aa in unique(usage$age)) {
  sub <- usage[usage$age == aa, ]
  lm <- loess(usage ~ year, data = sub, degree = 1, span = 2)
  
  preds <- data.frame("year" = years)
  preds <- as.data.frame(predict(lm, preds))
  preds$year <- years
  names(preds) <- c("usage", "year")
  preds$age <- aa
  preds$usage[preds$usage < 0] <- 0
  
  if(count == 1) {
    pred_frame <- preds
  } else {
    pred_frame <- rbind(pred_frame, preds)
  }
  
  count <- count  + 1
}


ggplot(data = pred_frame, aes(x = year, y = usage)) +
  geom_line() + 
  facet_wrap(~age)

write.csv(pred_frame, file = "data/smoothed_condom_usage.csv", row.names = FALSE)

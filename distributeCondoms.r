##################################################
## Allen Roberts
## July 30, 2015
## Description: Distributes people into condom users to match survey reports. Stratifies coverage by age
##################################################

distributeCondoms <- function(dt, time_index) {
  
  ## ART coverage by CD4 count
  coverage <- as.data.table(data.frame("age" = 1:12, prop =  
      sapply(seq(1, length(condom_cov)), function(x) {
        condom_cov[[x]][time_index]
      })
  ))
  props <- rbindlist(lapply(0:1, function(x, d) data.table(d, condom = x), d = coverage))
  props[condom == 0, prop := 1 - prop]
  setkey(props, age, condom)
  
  ## Sums the population across CD4 categories
  dt[, c("condom", "sum") := list(condom, sum(count)), by = .(hiv, age, male, risk, cd4, vl, circ, prep, art)]
  
  ## Multiplies the summed population by the CD4 proportions
  setkey(dt, age, condom)
  dt[props, count := sum * prop]
  
  dt[, sum := NULL]
  
}
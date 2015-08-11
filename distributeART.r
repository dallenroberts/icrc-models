##################################################
## Allen Roberts
## July 27, 2015
## Description: Distributes HIV+ people onto ART to match expected coverage. Stratifies coverage by CD4 count
##################################################

distributeART <- function(dt, time_index) {
  
  ## ART coverage by CD4 count
  coverage <- as.data.table(data.frame("cd4" = 1:5, prop =  
    sapply(seq(1, length(art_cov)), function(x) {
      art_cov[[x]][time_index]
    })
  ))
  props <- rbindlist(lapply(0:1, function(x, d) data.table(d, art = x), d = coverage))
  props[art == 0, prop := 1 - prop]
  props$hiv <- 1
  setkey(props, hiv, cd4, art)
 
  ## Sums the population across CD4 categories
  dt[, c("art", "sum") := list(art, sum(count)), by = .(hiv, age, male, risk, cd4, vl, circ, prep, condom)]
  
  ## Multiplies the summed population by the CD4 proportions
  setkey(dt, hiv, cd4, art)
  dt[props, count := sum * prop]
  
  dt[, sum := NULL]
  
}
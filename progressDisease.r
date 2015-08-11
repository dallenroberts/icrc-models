##################################################
## Allen Roberts
## July 7, 2015
## Description: Disease progression function (CD4 count and viral load)
##################################################

## This function takes the population data.table as input and progresses HIV+ persons (not on ART) through CD4 and viral load categories. Note that it modifies dt by reference, which points to the data.table that was originally passed to it.  This means that it's modifying that data.table in the parent environment.

progressDisease <- function(dt, time_step) {
  
  ## Note that we want the probability of progressing in a given time step.  This can be obtained from the mean duration using the exponential decay function (assuming a constant rate).  http://hyperphysics.phy-astr.gsu.edu/hbase/nuclear/meanlif.html
  
  setkey(dis_prog, hiv, male, vl, cd4, art)
  setkey(dt, hiv, male, vl, cd4, art)
  dt[dis_prog, c("cd4_duration", "vl_duration") := list(cd4_duration, vl_duration)]
  
  # dt <- merge(dt, dis_prog, all.x = TRUE)
  
  ## CD4 Efflux
  dt[!is.na(cd4_duration), diff := diff - count * (1 - exp(-time_step/cd4_duration))]
  
  ## VL Efflux
  dt[!is.na(vl_duration), diff := diff - count * (1 - exp(-time_step/vl_duration))]
  
  ## Influx - need to get counts from previous CD4 and VL categories, calculate how many people leave the category, increment the category by 1, and merge back to the original data frame and add that number of people.  Probably a cleaner way of doing this
  ## CD4 Influx
  prev_cd4 <- copy(dt[vl != 0 & cd4 >= 1 & cd4 < 5 & hiv == 1 & art == 0])
  prev_cd4[, cd4_influx := count * (1 - exp(-time_step/cd4_duration))]
  prev_cd4[, cd4 := cd4 + 1]
  prev_cd4[, c("count", "cd4_duration", "vl_duration", "diff") := NULL]
  setkeyv(prev_cd4, all_keys)
  setkeyv(dt, all_keys)
  dt[prev_cd4, cd4_influx := cd4_influx]
  # dt <- merge(dt, prev_cd4, all.x = TRUE)
  dt[!is.na(cd4_influx), diff := diff + cd4_influx]
  dt[, cd4_influx := NULL]
  
  ## VL Influx
  prev_vl <- copy(dt[cd4 != 0 & vl >= 1 & vl < 5 & hiv == 1 & art == 0])
  prev_vl[, vl_influx := count * (1 - exp(-time_step/vl_duration))]
  prev_vl[, vl := vl + 1]
  prev_vl[, c("count", "cd4_duration", "vl_duration", "diff") := NULL]
  setkeyv(prev_vl, all_keys)
  setkeyv(dt, all_keys)
  dt[prev_vl, vl_influx := vl_influx ]
  # dt <- merge(dt, prev_vl, all.x = TRUE)
  dt[!is.na(vl_influx), diff := diff + vl_influx]
  dt[, c("vl_influx", "cd4_duration", "vl_duration") := NULL]
  
}

##################################################
## Allen Roberts
## July 7, 2015
## Description: Function for linear interpolation
##################################################

## This function is used to interpolate linearly between observed values to create a smoothed time series.  The "breaks" argument contains the sorted year breaks (ie, 1990, 2000) and the "values" argument contains the values corresponding to the breaks.  "step_size" specifies the steps in years (ie, 0.25 = quarter-year step).  Default arguments for the range of the time series are the global variables year_start and year_end.

## The function will return a vector of values.  The length of the vector will be the total number of time steps between (and including) the start and end years.  The values at the breaks will be the values passed in the "values" argument to the function.  The values between the breaks will be linear interpolations between the passed values.  If the minimum and maximum breaks are greater than or less than (respectively) the year_start and year_end, the function will return repeated corresponding values for the years between year_start and the minimum break or year_end and the maximum break.

interpolate <- function(breaks, values, step_size = tstep, start = year_start, end = year_end) {
  
#   breaks <- deltas$year
#   values <- deltas$delta
#   step_size = tstep
#   start = year_start
#   end = year_end
#   
  smoothed <- rep(0, (year_end - year_start + 1)/ step_size)
  
  ## Checks
  if(!is.numeric(breaks) | !is.numeric(values) | !is.numeric(step_size)) stop("Error: breaks, values, and steps must both be numeric vectors")
  if(any(is.na(breaks) | is.na(values))) stop("Error: Missing breaks or values" )
  if(min(breaks) < year_start | max(breaks) > year_end) stop("Error: breaks fall outside of year limits")
  if(any(sort(breaks) != breaks)) stop("Error: breaks are not sorted")
  
  ## Starting values
  if(breaks[1] != year_start) {
    smoothed[1:((breaks[1] - year_start)/step_size)] <- values[1]
  }
  
  ## Ending values
  if(breaks[length(breaks)] != year_end) {
    smoothed[((breaks[length(breaks)] - year_start) / step_size + 1):length(smoothed)] <- values[length(values)] 
    
  } else {
    smoothed[(length(smoothed) - 1/step_size + 1):length(smoothed)] <- values[length(values)]
  }
  
  ## Interpolated values
  for(ii in 1:(length(values) - 1)) {
    
    ## Add specified break values
    smoothed[((breaks[ii] - year_start)/step_size + 1)] <- values[ii]
    
    ## Interpolate between specified break values
    slope <- (values[ii + 1] - values[ii])/((breaks[ii + 1] - breaks[ii])/step_size - 1)
    smoothed[((breaks[ii] - year_start)/step_size + 2):((breaks[ii + 1] - year_start)/step_size)] <- values[ii] + slope * seq(1, ((breaks[ii + 1] - breaks[ii])/step_size - 1))
  }
  
  return(smoothed)
  
}

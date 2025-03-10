# this model cacluates the hydrometeor temperature accounting t.1 and rh.1
# which is used in logistic regression mapping
compute_Ti <- function(df) {
  
  # Compute saturated vapor pressure (Tetens formula)
  saturated_vapor_pressure <- function(T) {
    return(6.112 * exp((17.62 * T) / (243.12 + T)))  # In hPa
  }
  
  # Compute actual vapor pressure from RH
  actual_vapor_pressure <- function(T, RH) {
    return(saturated_vapor_pressure(T) * (RH / 100))
  }
  
  # Function to compute hydrometeor temperature (Ti) iteratively
  compute_Ti <- function(Ta, RH) {
    L <- 2.5e6  # Latent heat of vaporization (J/kg)
    lt <- 0.024  # Thermal conductivity of air (W/m·K)
    D <- 2.5e-5  # Diffusivity of water vapor in air (m^2/s)
    
    Ti <- Ta  # Initial guess (starting with air temperature)
    for (i in 1:10) {  # Iterative approach
      r_Ta <- actual_vapor_pressure(Ta, RH) / (461.5 * (Ta + 273.15))  # Water vapor density
      r_sat_Ti <- saturated_vapor_pressure(Ti) / (461.5 * (Ti + 273.15))  # Saturated vapor density at Ti
      
      # Update Ti using the Pomeroy equation
      Ti_new <- round(Ta + (D / lt) * L * (r_Ta - r_sat_Ti), 3)
    }
    return(Ti)
  }
  
  # Apply compute_Ti function to each row of the dataframe
  df$Ti <- mapply(compute_Ti, df$t.1, df$rh.1)
  
  return(df)  # Return updated dataframe
}

# Define the derivative function
dydx <- function(x, y) {
  return(x^2 + y^2)
}

# Define the Runge-Kutta method
rk4 <- function(y0, x0, h, n) {
  # Initialize arrays
  x <- rep(0, n+1)
  y <- rep(0, n+1)
  
  # Set initial values
  x[1] <- x0
  y[1] <- y0
  
  # Iterate over n steps
  for (i in 1:n) {
    k1 <- h * dydx(x[i], y[i])
    k2 <- h * dydx(x[i] + h/2, y[i] + k1/2)
    k3 <- h * dydx(x[i] + h/2, y[i] + k2/2)
    k4 <- h * dydx(x[i] + h, y[i] + k3)
    y[i+1] <- y[i] + (k1 + 2*k2 + 2*k3 + k4) / 6
    x[i+1] <- x[i] + h
  }
  # Return the arrays
  return(data.frame(x=x, y=y))
}

# Example usage
result=rk4(y0=0, x0=0, h=0.05, n=20)
plot(result$x, result$y, type='b', col='red', main = '4th-order Runge-Kutta Method');grid()


print(result)

rm( list = ls() )
fun.reaction <- function (x, dt, t.end, rt.col = 1:3){
  f1 <- function(sigma, x){
    sigma * x[2] - sigma * x[1]
  }
  f2 <- function(rho, x){
    rho * x[1] - x[1] * x[3] - x[2]
    # rho * x[1] -  x[2]
  }
  f3 <- function(beta, x){
    x[1] * x[2] - beta * x[3]
    # - beta * x[3]
  }
  NT = t.end / dt
  mat= matrix(0, NT,3)
  for( i in 1:NT){
    x[1] = x[1] + f1(sigma, x) * dt
    x[2] = x[2] + f2(rho, x) * dt
    x[3] = x[3] + f3(beta, x) * dt
    mat[i, ]= x
  }
  ret = cbind(1:NT * dt,  mat[, rt.col])
  colnames(ret) = c('Time', 'x', 'y', 'z')
  ret= data.frame(ret)
  ret
}

sigma = 10; beta = 8/3; rho = 28

x= c(1, 1, 1) # IC
t.end = 50
dt = 0.01

x0 = c(10, 2, 1)
x1 = fun.reaction(x = x0, dt, t.end)

#' ==================================================================
#' ==================================================================
x = x1
par(mfrow=c(2,2))
plot(x$x, x$y, type = 'l')
plot(x$x, x$z, type = 'l')
plot(x$y, x$z, type = 'l')
par(mfrow=c(1,1))

# 
# rgl::plot3d(x[, 2:4])
# stop()
#' ==================================================================
#' ==================================================================
icol=1
epsilon = c(0,1,0) * 10^(-14)

x2 = fun.reaction(x = x0 + epsilon, dt, t.end)
# x2 = fun.reaction(x = x0+ c(0, 10^-13, 0), dt, t.end)
tr = (1:nrow(x1))[x1[, 1] > 40]
tr = (1:nrow(x1))[]
par(mfrow=c(3,1), mar=c(2, 4, 1, 1))
plot(x1$Time[tr], x1$x[tr], type='l', col=1, xlab='Time', ylab='X'); grid()
lines(x2$Time[tr], x2$x[tr], col=2)

plot(x1$Time[tr], x1$y[tr], type='l', col=1, xlab='Time', ylab='Y'); grid()
lines(x2$Time[tr], x2$y[tr], col=2)

plot(x1$Time[tr], x1$z[tr], type='l', col=1, xlab='Time', ylab='Z'); grid()
lines(x2$Time[tr], x2$z[tr], col=2)
#' ==================================================================
#' ==================================================================
library(plot3D)
x=x1
par(mfrow=c(1,1), mar=c(1, 1, 1, 1))
scatter3D (x=x$x, y=x$y, z=x$z, type = "l", theta=45, phi=10, bty='b2')

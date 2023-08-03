rm(list=ls())
library(plot3D)
source("Code/ch03/ch3_ConfAq2DEx.R")
source("Code/ch03/ch3_ConfAq2DIm.R")
Tmax = 1e4
dx = dy = 50
t0 = Sys.time()
x1 = CA.Explicit(Tmax = Tmax, DT = 25, 
                 Lxy = c(1000, 1000), Dxy = c(dx, dy),
                 ignore.cfl = TRUE)
t1 = Sys.time()
tu.ex = t1 - t0

t0 = Sys.time()
x2 = CA.Implicit(Tmax = Tmax, DT = 25, 
                 Lxy = c(1000, 1000), Dxy = c(dx, dy),
                 ignore.cfl = TRUE)
t1 = Sys.time()
tu.im = t1 - t0

df = data.frame( 'Type' = c('Explicit', 'Implicit'),
                 'NT' = c(x1$NT, x2$NT), 
                 'DT' = c(x1$DT, x2$DT),
                 'Time' = c(tu.ex, tu.im),
                 'CFL.X' = c(x1$CFL[1], x2$CFL[1]),
                 'CFL.Y' = c(x1$CFL[2], x2$CFL[2])
)
print(df)

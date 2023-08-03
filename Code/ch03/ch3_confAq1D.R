#' Problem: 1D Confined Aquifer
#'  governing Eqn: du/dt = k * B * (dd u / d x^2)
#' BC U0 = 10, UL=5
#' IC uic = 2 
#'  X = c(0, 1000)
#'  D = 23 mm2/s = 2.3e-5 m^2/s
#'  DX = 1 m
#'  DT = 10 s
#'  Time = 0 to 1e6 s

rm(list=ls())

ConfAq.explicit <- function(
  U0=10, UL=5, uic = 2,
  X = 1000,
  DX= 10,
  DT = 1,
  DD = 2e-2,
  T1 = 1e5,
  epsilon = 1e-4, bc2 = NULL,  ignore.cfl = FALSE, plot = TRUE
){
  T0 = 0
  tt = seq(T0, T1, DT)
  NT = length(tt)
  xx = seq(0, X, DX)
  NX = X / DX + 1
  
  alpha = DD  * DT / (DX * DX)
  beta = 1-2*alpha
  CFL = DD  * DT / (DX * DX)
  print(CFL)  
  if(!ignore.cfl){
    if(CFL >=1 ){
      warning('CFL  = ', CFL)
      stop()
    }
  }
  #' =========================================
  x0 = rep(uic, NX)
  #' =========================================
  ylim = c(min(x0), max(U0, UL))
  xlim=c(0, X)
  if(plot){
  plot(xx, x0, type='b', col=2,  lwd=3, ylim=ylim, xlim=xlim, xlab=xlab, ylab=ylab); grid()
  lines(x=c(1,1) * 0, y=c(min(x0), U0), lwd=2, col=3, type='b')
  lines(x=c(1,1) * X, y=c(min(x0), UL), lwd=2, col=3, type='b')
  text(x=X/2 , y = uic + diff(ylim)*.051, 'Initial condition', font=2)
  text(x=X * 0.05 , y = U0, 'BC 1', font=2)
  text(x=X * 0.95 , y = UL, 'BC 2', font=2)
  }
  #' =========================================
  mat = matrix(0, nrow = NX, ncol = NX)
  for(i in 1:NX){
    for(j in 1:NX){
      if(i==j){
        mat[i, j] = beta
      }
      if(i+1 == j | i-1 == j ){
        mat[i, j] = alpha
      }
    }
  }
  mat[1,]=c(1, rep(0, NX-1))
  mat[NX,]=c(rep(0, NX-1), 1)
  
  xm = matrix(NA, nrow=NX, ncol=NT)
  vs = cbind(rep(0, NX))
  if(!is.null(bc2)){
    vs[bc2$id] = bc2$val
  }
  
  b=bx = cbind(x0)
  for(i in 1:NT){
    if(i>1){
      bx = mat %*% b + vs * DT
      if(any(is.nan(bx))) { break }
      if(mean(abs(b-bx)) < epsilon ) { break }
    }
    bx[1] = U0
    bx[NX] = UL
    xm[, i]=bx
    b = bx
  }
  NT = i
  xm=xm[, 1:NT]
  message('CFL value = ', CFL)  
  message('Total Timesteps (dt * nt)= ', DT, ' * ', NT)
  ret = list('x' = xx, 
             'time' = tt,
             'u' = xm,
             'CFL' = CFL,
             'DT' = DT,
             'NT' = NT,
             'xlim' = xlim,
             'ylim' = ylim)
  return(ret)
}
plot1 <- function(x, nout = 20){
  NT = x$NT
  # id = sort(unique(c(1, 5, 10, 20, round(seq(1, NT, length.out = min(10, NT))))))
  id=10^seq(0, log10(x$NT), length.out = nout)
  col=colorspace::diverge_hcl(n=length(id)); lty=1
  matplot(x=x$x, y=x$u[, id], type='l', ylim=x$ylim, xlim=x$xlim, 
          xlab=xlab, ylab=ylab, col=col, lty=lty)
  legend('topright', paste0('T=', x$time[id]+1), col=col, lty=lty, bg='transparent')
}
plot2 <- function(x){
  NT = x$NT
  id = c(2, 4, 6, 8); nid=length(id)
  lty=1:nid; col=lty
  matplot(t(x$u[id, ]), type='l', xlab=xlab, ylab=ylab, col=col, lty=lty); grid()
  legend('bottomright',col=col, lty=lty, paste('Node', id))
}
#' =========================================
xlab ='Distance (m)'
ylab = 'Temperature (C)'
x = ConfAq.explicit(DX= 25, DT = 5, U0=10, UL=5, uic = 2.5,
                    DD=50, epsilon = 1e-5, ignore.cfl = FALSE,
                    bc2 = list(id = 20:30, val= +0.001))
plot1(x)
# plot2(x)


#' Problem: 1D Heat Transfer
#'  governing Eqn: du/dt = k/r/c * (dd u / d x^2)
#'  wiki: https://en.wikipedia.org/wiki/Thermal_diffusivity
#'  wiki: https://en.wikipedia.org/wiki/Heat_equation
#' BC U0 = 100, UL=50
#' IC uic = 25 
#'  X = c(0, 1)
#'  D = 23 mm2/s = 2.3e-5 m^2/s
#'  DX = 0.01 m
#'  DT = 10 s
#'  Time = 0 to 1e6 s

HT.implicit <- function(
  U0=100, UL=50, uic = 25,
  X = 1,
  DX= 0.1,
  DT = 1,
  DD = 2.3e-5,
  Tmax = 1e5,
  epsilon = 1e-4, bc2 = NULL,  ignore.cfl = FALSE, plot = TRUE
){
  T0 = 0
  tt = seq(T0, Tmax, DT)
  NT = length(tt)
  xx = seq(0, X, DX)
  NX = X / DX + 1
  
  alpha = -DD  * DT / (DX * DX)
  beta = 1 + 2 * DD  * DT / (DX * DX)
  CFL = DD  * DT / (DX * DX)
  print(CFL)  
  if(!ignore.cfl){
    if(CFL >=1 ){
      stop()
    }
  }
  #' =========================================
  x0 = rep(uic, NX)
  #' =========================================
  ylim = c(min(x0), max(U0, UL))
  xlim=c(0, X)
  if(plot){
    plot(xx, x0, type='b', col=2,  lwd=3, ylim=ylim, xlim=xlim, xlab=xlab, ylab=ylab)
    grid()
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
  # vs[NX/2] = ss
  b=bx=cbind(x0)
  for(i in 1:NT){
    if(i>1){
      bx = solve(mat, b)
      if(any(is.nan(bx))) { break }
      if(mean(abs(b-bx)) < epsilon) { break }
    }
    bx[1] = U0
    bx[NX] = UL
    xm[, i]=bx
    b = bx
  }
  NT = i
  xm=xm[, 1:NT]
  # message('CFL value = ', CFL)  
  # message('Total Timesteps (dt * nt)= ', DT, ' * ', NT)
  yy = xm; yy[abs(yy)>1e20] = NA
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
  id=10^seq(0, log10(x$NT), length.out = nout)
  col=colorspace::diverge_hcl(n=length(id)); lty=1
  matplot(x=x$x, y=x$u[, id], type='l', ylim=x$ylim, xlim=x$xlim, 
          xlab=xlab, ylab=ylab, col=col, lty=lty)
  legend('topright', paste0('T=', x$time[id]+1), col=col, lty=lty, bg='transparent')
  mtext(text = paste('CFL =', x$CFL ), side=3, cex=1.5)
}
plot2 <- function(x){
  NT = x$NT
  id = c(2, 4, 6, 8); nid=length(id)
  lty=1:nid; col=lty
  matplot(t(x$u[id, ]), type='l', xlab=xlab, ylab=ylab, col=col, lty=lty); grid()
  legend('bottomright',col=col, lty=lty, paste('Node', id))
}


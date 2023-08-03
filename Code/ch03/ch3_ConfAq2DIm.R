
#' Problem: 1D Confined Aquifer
#'  governing Eqn: du/dt = DDx * (dd u / d x^2) + DDy * (dd u / d y^2)
#'  
diag.matrix <- function(id = c(-1, 0, 1), x = rep(1, length(id)), n = 3, def.val = 0){
  val = matrix(x, ncol=length(id), nrow=1)
  mat = matrix(def.val, n, n)
  nid = length(id)
  for(i in 1:n){
    for(j in 1:n){
      for(k in 1:nid){
        if(i + id[k] == j){
          mat[i,j] = val[k]
        }
      }
    }
  }
  mat
}

toBC <- function(idl, x, val){
  nbc = length(idl)
  for(i in 1:nbc){
    x[idl[[i]]] = val[i]
  }
  x
}

CA.Implicit <- function(  bc1 = c(0, 0, 0,0),
                          bc2 = NULL,
                          uic = 25,
                          Lxy = c(1000, 1000),
                          Dxy = c(50,50),
                          DD = rep(23 ,2),
                          epsilon = 0.001,
                          DT = 25,
                          Tmax = 1e5,
                          ignore.cfl = FALSE, plot = TRUE){
  DX=Dxy[1]; DY=Dxy[2];
  tt = seq(0, Tmax, DT); 
  NT = length(tt)
  xx = seq(0, Lxy[1], DX); NX = length(xx)
  yy = seq(0, Lxy[2], DY); NY = length(yy)
  # NX = Lxy[1] / DX + 1;   NY = Lxy[2] / DY + 1
  N = NX * NY
  alpha = -DD[1]  * DT / (DX * DX)
  beta  = -DD[2]  * DT / (DY * DY)
  CFL.x =  DD[1]  * DT / (DX * DX)
  CFL.y =  DD[2]  * DT / (DY * DY)
  gamma = 1 + 2 * DD[1]  * DT / (DX * DX) + 2 * DD[2]  * DT / (DY * DY)
  message('CFL value = (', CFL.x, '\t', CFL.y, ')')  
  if(!ignore.cfl){
    if(CFL.x >=.5 | CFL.y >=.5){
      stop()
    }
  }
  #' =========================================
  x0 = rep(uic, NX)
  mat = diag.matrix(id = c(-NY, -1, 0, 1, NY), n=N,
                    x=c(alpha, beta, gamma, beta, alpha), def.val = 0)
  dmat = diag.matrix(id=0, x=1, n=N, def.val = 0)
  idl = list(1:NY,
             1+(1:NX - 1)*(NY),
             (1:NX) * NY,
             (NX-1)*(NY)+1:NY) 
  nbc = length(idl)
  id.bc = sort(unique(unlist(idl)))
  mat[id.bc, ] = dmat[id.bc,]
  
  arr = array(NA, dim=c(NY,NX,NT))
  vs = cbind(rep(0, N))
  b=bx=cbind(rep(uic, N))
  b = toBC(idl = idl, x=b, val=bc1)
  
  for(i in 1:NT){
    if(i>1){
      bx = solve(mat, b + vs * DT)
      if(any(is.nan(bx))) { break }
      if(mean(abs(b-bx)) < epsilon) { break }
    }
    bx = toBC(idl = idl, x=bx, val=bc1)
    arr[, , i] = matrix(bx, nrow = NY, ncol = NX)
    b = bx
  }
  NT = i
  arr=arr[,, 1:NT]
  # message('Total Timesteps (dt * nt)= ', DT, ' * ', NT)
  # yy = arr; yy[abs(yy)>1e20] = NA
  ret = list('x' = xx,
             'y' = yy,
             'z' = arr,
             'time' = tt,
             'CFL' = c(CFL.x, CFL.y),
             'DT' = DT,
             'NT' = NT
  )
  return(ret)
}
plot.3d <- function(x, nr=3, nc=4, clim=NULL){
  par(mfrow=c(nr, nc), mar=c(1,1,1,1))
  idx = round(10^seq(0, log10(x$NT), length.out = nc*nr))
  z=x$z; 
  z[ is.infinite(abs(z)) ] = NA
  if(is.null(clim)){
    clim = range(z, na.rm = TRUE)
  }
  for(i in idx ){
    plot3D::persp3D(z=z[, , i], clim=clim,
            colvar=z[, , i])
    mtext(paste('T =', x$DT * i), side= 3, line=-1)
  }
}



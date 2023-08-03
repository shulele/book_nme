rm(list=ls())

library(plot3D)
library(grid)
logo <- function(col = jet.col(), resfac=1){
  # par(mfrow = c(1, 1), mar=c(1,1,1,1) * 0)
  AA <- plot3D::Hypsometry$z;  AA[AA<=0] <- NA
  dd = 55
  ZZ = rbind(AA[(dd+1):359, ], AA[1:dd, ])
  ZZ = ZZ ^ 0.75
  plot3D::spheresurf3D(ZZ, NAcol = rgb(0, 94/255, 220/255, .25), 
                       resfac = resfac, # output resolution
                       box = T,bty = "g", colkey = FALSE, 
                       # xlim = lim, ylim = lim, zlim = lim, 
                       theta = 90+dd, phi = 25, 
                       col = col, border='gold',
                       lwd=0.1,
                       lighting = TRUE, ltheta = 135+dd, lphi = 35, # Light
                       shade = TRUE)
}

plot.cover <- function(wd = 53, 
                       ht = 81, resfac = 1,
                       fts = c('Kai', 'Hei', 'Arial', 'Hei', 'Hei') 
){
  X = wd * resfac
  Y = ht * resfac
  xlim=c(0, X)
  ylim=c(0, Y)
  #' ==== 背景颜色调色板 =====
  colfunc <- colorRampPalette(c(rep('#3C293C', 9), "#342D35","#2D2B2E", "#292929"), bias=1)
  # colfunc <- colorRampPalette(c("#292929", "#2D2B2E", "#342D35", '#3C293C', '#3C293C'), bias=10)
  g <- matrix(colfunc(ht),nrow=ht, ncol=wd)
  
  
  par(new=FALSE, mar=rep(0, 4))
  #' ==== 紫色背景 =====
  plot(0, type='n', axes=F, xlim=xlim, ylim=ylim, xlab='', ylab='')
  grid::grid.raster(g)
  
  #' ==== 白色圆角矩形 =====
  grid::grid.roundrect(x=.5, y=0.47, width=.85, height=.42, name="rr", just='center')
  #' ==== 左上角白色条 =====
  # grid::grid.roundrect(x=1, y=0.933, width=.39, height=0.04,  name="rr", just='right',)
  
  par(new=TRUE, family=fts[1])
  # plot(0, type='n', axes=F, xlim=xlim, ylim=ylim, xlab='', ylab='')
  #' ==== 文字 =====
  xy = cbind(c(0.648, rep(.03, 4)) * X,  
             c(1.0, .86, .76, .18, .1)*Y )
  col.text = c('gray40', 'white', 'white', 'white', 'white')
  txt=c('   高性能科学计算', ## '数值计算方法'
        '数值地球','Numeric Methods in Earth Sciences',  
        '舒乐乐，孟宪红','中国科学院西北生态环境资源研究院')
  adj=c(0, 0, 0, 0, 0)
  cex = c(1.5, 7, 2.5, 2.5, 2.5) * resfac
  font=c(2, 2, 2, 2, 2)
  
  for(i in 1:5){
    if(i==1){
      legend(x=xy[i, 1], y=xy[i, 2], legend = txt[1], adj=adj[i],  cex=cex[i], 
             text.font=font[i],text.col = col.text[i], text.width = 800*resfac,
             # box.col = 'white',
             bg = "white", box.lwd = 0)
    }else{
      text(x=xy[i, 1], y=xy[i, 2], labels = txt[i],
           col=col.text[i], adj=adj[i], cex=cex[i], font=font[i], family=fts[i])
    }
  }
  #' ==== 地球logo =====
  par(new=TRUE, fig=c(0.1, 0.9, c(0.40, 0.80) - 0.13), mar=rep(0, 4))
  logo(resfac = .5)
  
  #' ==== 曲线 =====
  fx <- function(n, npi = 4){
    y=0; x = seq(-1 * npi * pi, npi * pi, length.out=100*npi)
    for(i in 1:n){ y = y+ sin(x/i)/i }
    y
  }
  y = fx(6, 4)
  par(new=TRUE, fig=c(0.1, 0.9, c(0.40, 0.60) - 0.13), mar=rep(0, 4))
  plot( y, col='gray80', lwd=6*resfac, type='l', xlab='', ylab='', axes=F)
  lines(y, col=2, lwd=1.5*resfac)
  #' ====  =====
}
# 
library(plot3D)
library(grid)
dir.create('image', showWarnings = FALSE, recursive = TRUE)
resfac = .5
jpeg('image/Cover.jpg', width=8.48*resfac, height=12.96*resfac, res=200, unit='in')
plot.cover(resfac = resfac,
           fts = c('Kai', 'Hei', 'Arial', 'Songti SC', 'Songti SC') )
dev.off()

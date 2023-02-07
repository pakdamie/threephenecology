library(RColorBrewer)
image.mvcwt2_CC1_Modified_Supp <-  function (x, z.fun = "Re", bound = 1, reset.par = TRUE, maxCol = 1, minCol = 0, family,...)
{
  z.fun = match.fun(z.fun)
  opar = par(no.readonly = TRUE)
  if (reset.par)
    on.exit(par(opar))
  pal = colorRampPalette(rev(brewer.pal(11, "Spectral")))(1024)[round(1+ minCol*1024):round(maxCol*1024)]
  with(x, {
    nvar = ifelse(length(dim(z)) == 3, dim(z)[3], 1)
    par(mfrow = c(nvar, 1), mar = rep(0.2, 4), oma = rep(5,
                                                         4))
    for (i in 1:nvar) {
      image(x, y, z.fun(z[, , i]), log = "y", col = pal,
            axes = FALSE, ...,)
      if (i%%2)
        axis(2)
      else axis(4)
      if (exists("z.boot") && !is.null(z.boot)) {
        z.boot = 1 - abs(1 - 2 * z.boot)
          contour(x, y, z.boot[, , i], levels = 0.05, lty = 3,
            add = TRUE, drawlabels = FALSE)
        zb = p.adjust(as.vector(z.boot), method = "BY")
        dim(zb) = dim(z.boot)
           contour(x, y, zb[, , i], levels = 0.05, lwd = 2,
            add = TRUE, drawlabels = FALSE)
      }
      if (is.finite(bound)) {
        lines(min(x) + bound * y, y, lty = 2, lwd = 2,
              col = "darkgrey")
        lines(max(x) - bound * y, y, lty = 2, lwd = 2,
              col = "darkgrey")
      }
      box()
    }
    axis(1, at= c(2,4,6,8),labels=c(2004,2006,2008,2010))
    mtext("", 1, 3, outer = TRUE)
    mtext("Scale", 2, 3, outer = TRUE)
    title(family, outer=TRUE)
  })
  return(invisible(x))
}

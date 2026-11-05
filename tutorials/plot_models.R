# Path diagrams for the two-trait tutorial models, drawn from their R objects.
# Rectangles denote traits, ellipses denote sources, and arrow labels are
# parameter names or fixed values. No fitted coefficients are displayed.
plot_causal_model <- function(model, main = "") {
  if (!inherits(model, "gcmsem_model") || length(model$observed) != 2L ||
      length(model$latent) > 1L)
    stop("This tutorial helper expects two observed traits and at most one latent source")
  old <- graphics::par(mar = c(.4, .4, 3, .4), xpd = NA)
  on.exit(graphics::par(old), add = TRUE)
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1), asp = NA)
  graphics::title(main = main, col.main = "#203447", cex.main = 1.25)
  observed <- model$observed
  pos <- rbind(c(.22, .48), c(.78, .48))
  rownames(pos) <- observed
  if (length(model$latent)) pos <- rbind(pos, c(.5, .84))
  pos <- as.matrix(pos)
  rownames(pos) <- c(observed, model$latent)
  ellipse <- function(x, y, rx, ry, label, fill = "#EDF3F7", cex = 1) {
    theta <- seq(0, 2*pi, length.out = 150)
    graphics::polygon(x+rx*cos(theta),y+ry*sin(theta),col=fill,border="#416078",lwd=1.7)
    graphics::text(x,y,label,cex=cex,col="#203447")
  }
  arrows <- Filter(function(r) r$kind == "path", model$rows)
  both <- all(vapply(list(c(observed[1],observed[2]),c(observed[2],observed[1])),
    function(pair) any(vapply(arrows,function(r)identical(model$nodes[r$tuple],pair),logical(1))),logical(1)))
  for (r in arrows) {
    to <- model$nodes[r$tuple[1]]; from <- model$nodes[r$tuple[2]]
    label <- if(r$parameter > 0) names(model$start)[r$parameter] else as.character(r$value)
    a <- pos[from, ]; b <- pos[to, ]
    if (from %in% observed && to %in% observed) {
      direction <- sign(b[1]-a[1]); offset <- if(both) .065*direction else 0
      a <- a+c(.14*direction,offset); b <- b+c(-.14*direction,offset)
      graphics::arrows(a[1],a[2],b[1],b[2],length=.12,lwd=2,col="#2374AB")
      graphics::text(mean(c(a[1],b[1])),a[2]+.04*direction,label,col="#205A80",cex=.98)
    } else {
      a <- a+c(sign(b[1]-a[1])*.07,-.085); b <- b+c(0,.075)
      graphics::arrows(a[1],a[2],b[1],b[2],length=.12,lwd=1.8,col="#416078")
      graphics::text(mean(c(a[1],b[1]))+sign(b[1]-a[1])*.035,
                     mean(c(a[2],b[2]))+.01,label,col="#416078")
    }
  }
  for (v in observed) {
    xy <- pos[v, ]
    graphics::rect(xy[1]-.14,xy[2]-.075,xy[1]+.14,xy[2]+.075,
                   col="#DCEEF8",border="#2374AB",lwd=1.8)
    graphics::text(xy[1],xy[2],v,cex=1.18,col="#163D56")
    ellipse(xy[1],.22,.13,.065,paste0("u_",v),fill="#F3EFE6",cex=.95)
    graphics::arrows(xy[1],.285,xy[1],xy[2]-.075,length=.11,lwd=1.6,col="#7A6A4A")
    graphics::text(xy[1],.1,"Free variance, k3, k4",cex=.83,col="#665B48")
  }
  if (length(model$latent)) ellipse(.5,.84,.16,.10,paste0(model$latent,"\nGaussian source"),cex=.95) else
    graphics::text(.5,.84,"No common source\nStable branch requested in fitting",cex=1.05,col="#56616A")
  graphics::mtext("Arrow labels are parameters or fixed values; diagrams show model assumptions.",
                   side=1,line=-1.1,cex=.72,col="#56616A")
  invisible(model)
}

gen.sim.corr <- function(x, sigma = 0.001) {
     library(expm)
     r <- nrow(x)
     c <- ncol(x)
     r.c <- diag(r)
     c.c <- diag(c)
     r.m <- rowMeans(x)
     c.m <- colMeans(x)
     d.r.c <- 1
     d.c.c <- 1
     k <- 1
     while (d.r.c > sigma | d.c.c > sigma) {
          p.r.c <- r.c
          p.c.c <- c.c
          for (i in 1: r) {
               for (j in 1:r) {
                    if (i != j) {
                         r.x <- x[i, ] - r.m[i]
                         r.y <- x[j, ] - r.m[j]
                         r.xy <- r.x %*% c.c * t(r.y)
                         r.xx <- r.x %*% c.c * t(r.x)
                         r.yy <- r.y %*% c.c * t(r.y) 
                         r.num <- sum(r.xy)
                         r.den <- sqrt(sum(r.xx)) * sqrt(sum(r.yy))
                         r.c[i, j] <- round(r.num / r.den, 2)
                    }
               }
          }
          for (i in 1: c) {
               for (j in 1:c) {
                    if (i != j) {
                         c.x <- x[, i] - c.m[i]
                         c.y <- x[, j] - c.m[j]
                         c.xy <- c.x %*% r.c * t(c.y) 
                         c.xx <- c.x %*% r.c * t(c.x) 
                         c.yy <- c.y %*% r.c * t(c.y) 
                         c.num <- sum(c.xy)
                         c.den <- sqrt(sum(c.xx)) * sqrt(sum(c.yy))
                         c.c[i, j] <- round(c.num / c.den, 2)
                    }
               }
          }
          d.r.c <- sum(rowSums(p.r.c - r.c))
          d.c.c <- sum(rowSums(p.c.c - c.c))
          k <- k + 1
     }
     rownames(r.c) <- rownames(x)
     colnames(r.c) <- rownames(x)
     rownames(c.c) <- colnames(x)
     colnames(c.c) <- colnames(x)
     return(list(row.sims = r.c, col.sims = c.c, n.iter = k))
}
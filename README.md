# radish

Fast gradient-based optimization of resistance surfaces.

Requires [corMLPE](https://github.com/nspope/corMLPE): `devtools::install_github("nspope/corMLPE")`

```r
library(radish)
library(raster)

data(melip)

covariates <- raster::stack(melip.altitude, melip.forestcover)

# currently, it's crucial to scale variables
scale_to_0_1 <- function(x) (x - min(x, na.rm = TRUE))/(max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
values(melip.altitude) <- scale_to_0_1(values(melip.altitude))
values(melip.forestcover) <- scale_to_0_1(values(melip.forestcover))

covariates <- raster::stack(list(altitude=melip.altitude, forestcover=melip.forestcover))

plot(covariates[["altitude"]])
points(melip.coords)

surface <- radish_conductance_surface(covariates, melip.coords, directions = 8)
fit_mlpe <- radish_optimize(loglinear_conductance, mlpe, surface, melip.Fst)

fit_mlpe$ztable #MLE
fit_mlpe$loglik #loglikelihood
fit_mlpe$AIC    #AIC
fit_mlpe$vcor   #correlation in asympotitic distribution of MLE
plot(dist_from_cov(as.matrix(fit_mlpe$fit$covariance)), melip.Fst, pch = 19,
     xlab = "Optimized resistance distance", ylab = "Fst") #visualize fit

# visualise likelihood surface across grid (takes awhile)
theta <- as.matrix(expand.grid(x=seq(-6,6,length.out=21), y=seq(-6,6,length.out=21)))
grid <- radish_grid(loglinear_conductance, mlpe, surface, melip.Fst, theta, covariance=FALSE)

library(ggplot2)
ggplot(data.frame(loglik=grid$loglik, grid$theta)) + 
  geom_tile(aes(x=x,y=y,fill=-loglik)) + theme_bw() +
  geom_contour(aes(x=x,y=y,z=-loglik), color="black") +
  xlab(expression(theta[altitude])) +
  ylab(expression(theta[forestcover]))

# calculate resistance distances across grid
distances <- radish_distance(loglinear_conductance, surface, theta)

plot(distances$distance[,,which(theta[,1] == 0 & theta[,2] == 0)],
     melip.Fst, pch = 19, xlab = "IBD resistance distance", ylab = "Fst")


```
 
This is still very much a work-in-progress. Contact at nspope@utexas.edu.

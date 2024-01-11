rm(list = ls())
graphics.off()
cat("\14")

library(LakeEnsemblR)
library(ggplot2)
library(ggpubr)
library(lubridate)
library(reshape2)
library(tidyr)



qual_fun <- function(O, P){
  
  # Arguments:
  #^^^^^^^^^^
  # O: observed values
  # P: predicted values
  #
  # Return Value:
  #^^^^^^^^^^^^^^
  # qual: A data.frame containing the six quality estimates
  # set of both O and P where both have no NAs
  id <- !((is.na(O) | is.na(P)) | (is.na(O) & is.na(P)))
  O <- O[id]
  P <- P[id]
  
  # nash sutcliff
  nse <- 1 - sum((O - P)^2, na.rm = TRUE)/sum((O - mean(O, na.rm=TRUE))^2, na.rm = TRUE)
  
  # pearson corelation coef
  r <- sum((O - mean(O, na.rm = TRUE))*(P - mean(P, na.rm = TRUE)),
           na.rm = TRUE)/sqrt(sum((O - mean(O, na.rm = TRUE))^2, na.rm = TRUE)*
                                sum((P - mean(P, na.rm = TRUE))^2, na.rm = TRUE))
  
  # bias
  bias <- mean((P - O), na.rm = TRUE)
  
  # mean absolute error
  mae <- mean(abs(P - O), na.rm = TRUE)
  
  # mean square error
  mse <- mean((P - O)^2, na.rm = TRUE)
  
  # root mean square error
  rmse <- sqrt(mse)
  
  # Absolute Maximum Error
  ame <- max(abs(P - O))
  
  # Mean Relative Mean Error or Bias
  mrme <- mean((P - O)/O, na.rm = TRUE)
  
  # normalised mean absolute error
  nmae <- mean(abs((P - O)/O), na.rm = TRUE)
  
  # normalised mean error
  nme <- mean((P - O), na.rm = TRUE)/mean(O, na.rm = TRUE)  
  
  qual <- data.frame(rmse = rmse, nse = nse, r = r, bias = bias, mae = mae, nmae = nmae, mse=mse,ame=ame,mrme=mrme,nme=nme)
  
  return(qual)
}

# Plot depth and time-specific results
thm <- theme_pubr(base_size = 17) + grids()

config_file <- "LakeEnsemblR.yaml"

model <- c("GLM", "GOTM", "Simstrat", "FLake")

# Import the LER output into your workspace
ens_out <- paste0("output/ensemble_output.nc")

# load water temp for the cali and vali periods
wtemp <- load_var(ens_out, "temp",dim = "model")

# for the rviewers sake calculate metrics at each depth and then average
depth_qual <- list()
cols_c <- which(apply(wtemp$Obs, 2, function(x)sum(!is.na(x))) != 0)[-1]
rows_c <- which(apply(wtemp$Obs, 1, function(x)sum(!is.na(x))) > 1)

for (m in model) {
  depth_qual[[m]]$calibration <- melt(lapply(cols_c, function(c) qual_fun(wtemp$Obs[rows_c, c],
                                                                          wtemp[[m]][rows_c, c])))
  depth_qual[[m]]$calibration$L1 <- as.numeric(gsub("wtr_", "",depth_qual[[m]]$calibration$L1))
  colnames(depth_qual[[m]]$calibration) <- c("Metric", "Value", "Depth")

}

depth_qual <- reshape2::melt(depth_qual, id.vars = c("Depth", "Metric"))
depth_qual <- depth_qual[, -3]
colnames(depth_qual) <- c("Depth", "Metric", "Value", "Period", "Model")

palm_depth <- ggplot(depth_qual) + geom_line(aes(x = Depth, y = Value, col = Period)) +
  coord_flip() + facet_grid(Model ~ Metric, scales = "free_x")  + scale_x_reverse() + thm +
  xlab("Depth (m)")

write.csv(depth_qual,'depth_qual.csv')


tot_qual <- list()
for (m in model) {
  OO <- wtemp$Obs[rows_c,cols_c]
  OO<-melt(OO)
  PP <- wtemp[[m]][rows_c,cols_c]
  PP<-melt(PP)
  tot_qual[[m]] <- t(qual_fun(OO$value,PP$value ))
}

write.csv(tot_qual,'tot_qual.csv')


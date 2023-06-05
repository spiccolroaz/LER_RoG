#initial clean up
rm(list = ls())
graphics.off()
cat("\f")

# set working directory
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# load library
library(LakeEnsemblR)

# Set config file
config_file <- "LakeEnsemblR.yaml"
model <- c("FLake", "GLM", "GOTM", "Simstrat")
# 1. Example - creates directories with all model setup
export_config(config_file = config_file, model = c("FLake", "GLM", "GOTM", "Simstrat"),
              folder = ".")
# 2. Run ensemble lake models
run_ensemble(config_file = config_file,
             model = c("FLake", "GLM", "GOTM", "Simstrat"),
             return_list = FALSE, parallel = FALSE)

# path of the output netcdf file
ncdf <- "output/ensemble_output.nc"
# plot heatmap
plot_heatmap(ncdf)

### calibration
cali_res <- cali_ensemble(config_file = config_file, num = 10000, cmethod = "MCMC",
                          parallel = TRUE, model = model)
# get best parameters
best_par <- setNames(lapply(model, function(m)cali_res[[m]]$bestpar), model)
print(best_par)

# > print(best_par) #2023-5-31 3:16 PM CDT
# $FLake
# wind_speed         swr   c_relax_C
# 1.244507421 1.099784772 0.008672142
#
# $GLM
# wind_speed                 swr mixing/coef_mix_hyp
# 1.040228            1.073248            0.826737
#
# $GOTM
# wind_speed              swr turb_param/k_min
# 1.247361e+00     1.099835e+00     2.869447e-06
#
# $Simstrat
# wind_speed          swr     a_seiche
# 1.0862895980 0.9847077065 0.0008157647

# Load libraries for post-processing
library(gotmtools)
library(ggplot2)

export_config(config_file = config_file, model = c("FLake", "GLM", "GOTM", "Simstrat"),
              folder = ".")

run_ensemble(config_file = config_file,
             model = c("FLake", "GLM", "GOTM", "Simstrat"),
             return_list = FALSE, parallel = FALSE)

# path of the output netcdf file
ncdf <- "output/ensemble_output.nc"
# plot heatmap
plot_heatmap(ncdf)

## Plot model output using gotmtools/ggplot2
# Extract names of all the variables in netCDF
ncdf <- "output/ensemble_output.nc"
vars <- gotmtools::list_vars(ncdf)
vars # Print variables

p1 <- plot_heatmap(ncdf)
p1
# Change the theme and increase text size for saving
p1 <- p1 +
  theme_classic(base_size = 24) +
  scale_colour_gradientn(limits = c(0, 21),
                         colours = rev(RColorBrewer::brewer.pal(11, "Spectral")))
# Save as a png file
ggsave("output/ensemble_heatmap.png", p1,  dpi = 300,width = 384,height = 280, units = "mm")

calc_fit(ncdf = "output/ensemble_output.nc",
         model = c("FLake", "GLM",  "GOTM", "Simstrat"),
         var = "temp")

plist <- plot_resid(ncdf = "output/ensemble_output.nc",var = "temp",
                    model = c('FLake', 'GLM',  'GOTM', 'Simstrat'))

# Johannes:
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

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
  # rmse
  rmse <- sqrt(mean((O - P)^2, na.rm = TRUE))

  # nash sutcliff
  nse <- 1 - sum((O - P)^2, na.rm = TRUE)/sum((O - mean(O, na.rm=TRUE))^2, na.rm = TRUE)

  # pearson corelation coef
  r <- sum((O - mean(O, na.rm = TRUE))*(P - mean(P, na.rm = TRUE)),
           na.rm = TRUE)/sqrt(sum((O - mean(O, na.rm = TRUE))^2, na.rm = TRUE)*
                                sum((P - mean(P, na.rm = TRUE))^2, na.rm = TRUE))

  # bias
  bias <- mean((P - O), na.rm = TRUE)

  # mean absolute error
  mae <- mean(abs(O - P), na.rm = TRUE)

  # normalised mean absolute error
  nmae <- mean(abs((O - P)/O), na.rm = TRUE)

  qual <- data.frame(rmse = rmse, nse = nse, r = r, bias = bias, mae = mae, nmae = nmae)

  return(qual)
}

# Plot depth and time-specific results
thm <- theme_pubr(base_size = 17) + grids()

config_file <- "LakeEnsemblR.yaml"

model <- c("FLake","GLM", "GOTM", "Simstrat")#, "FLake")

# Import the LER output into your workspace
ens_out <- paste0("output/", get_yaml_value(config_file, "output", "file"),
                    ".nc")

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

palm_depth <- ggplot(depth_qual) + geom_line(aes(x = Depth, y = Value), size = 1.5) +
  coord_flip() + facet_grid(Model ~ Metric, scales = "free_x")  + scale_x_reverse() + thm +
  xlab("Depth (m)")

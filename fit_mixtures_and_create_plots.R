source("utils.R")

library(mixdist)
library(tidyverse)

files <- dir(path="data/", pattern="*.CSV", full.names=TRUE)
basenames <- dir(path="data/", pattern="*.CSV", full.names=FALSE)

fittings <- tibble(
  name = character(),
  norm.pi1 = numeric(),
  norm.pi2 = numeric(),
  norm.chisq = numeric(),
  lnorm.pi1 = numeric(),
  lnorm.pi2 = numeric(),
  lnorm.chisq = numeric(),
  lnorm.mu1 = numeric(),
  lnorm.mu2 = numeric(),
  lnorm.sigma1 = numeric(),
  lnorm.sigma2 = numeric()
)

# Calculate fittings
for (i in 1:length(files)) {
  data <- read.coulter.data(files[i])
  params <- get.fitting.parameters(data)
  fittings <- add_row(fittings, name=basenames[i], norm.pi1=params$norm.pi1, norm.pi2=params$norm.pi2, norm.chisq=params$norm.chisq,
                      lnorm.pi1=params$lnorm.pi1, lnorm.pi2=params$lnorm.pi2, lnorm.chisq=params$lnorm.chisq,
                      lnorm.mu1=params$lnorm.mu1, lnorm.mu2=params$lnorm.mu2, lnorm.sigma1=params$lnorm.sigma1, lnorm.sigma2=params$lnorm.sigma2)
}

for (i in 1:length(files)) {
  basename <- tools::file_path_sans_ext(basename(basenames[i]))
  data <- read.coulter.data(files[i])
  plot <- fit.and.compare(data, basename)
  fname <- paste("output/", tools::file_path_sans_ext(basename(basenames[i])), '.pdf', sep='')
  ggsave(fname, plot)
}

write_csv2(fittings, "output/all_fittings.csv")

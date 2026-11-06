

library(readxl)
setwd("C:/Users/lekshmi.venugopal/OneDrive - Southwest District Health/Documents/Outbreak")

cases <- read_xlsx("data.xlsx")
saveRDS(cases, "data.rds", compress = "xz")




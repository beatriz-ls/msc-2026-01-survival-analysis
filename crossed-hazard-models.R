library(tidyverse)
library(survstan)
library(GGally)


# Load ata -------------------------------------------------------
data(ipass)
glimpse(ipass)

ipass <- ipass %>%
  mutate(
    arm = as.factor(ifelse(arm == 1, "gefitinib", "carboplatin/paclitaxel"))
  )

# KM curves ------------------------------------------------------
km <- survfit(Surv(time, status) ~ arm, data = ipass)

## Graphic
ggsurv(km)

# ajustando os modelos PH, PO e YP:

m <- ceiling(sqrt(nrow(ipass)))
dist <- "loglogistic"

aft <- aftreg(Surv(time, status) ~ arm, data = ipass, dist = dist)
ah <- ahreg(Surv(time, status) ~ arm, data = ipass, dist = dist)
eh <- ehreg(Surv(time, status) ~ arm, data = ipass, dist = dist)
ph <- phreg(Surv(time, status) ~ arm, data = ipass, dist = dist)
po <- poreg(Surv(time, status) ~ arm, data = ipass, dist = dist)
yp <- ypreg(Surv(time, status) ~ arm, data = ipass, dist = dist)

newdata <- data.frame(
  arm = as.factor(factor(c("gefitinib", "carboplatin/paclitaxel")))
)

surv <- survfit(ah, newdata = newdata)

# Plotando as curvas estimadas:
ggsurv(km) +
  geom_line(data = surv, aes(x = time, y = surv, color = arm, group = arm))


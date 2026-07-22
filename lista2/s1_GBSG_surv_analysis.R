#install.packages("shrink")
#install.packages(c('bellreg', 'LambertW', 'peppm', 'YPPE'))
#install.packages(file.choose(), repos = NULL, type="source") # instalando survcurv

library(survival)
library(survstan)
library(tidyverse)
library(shrink)
library(gtsummary)

# Script 1: Modelagem do dataset GBSG
# [x] Curvas de KM
# [x] Testes de comparação de curvas para váriáveis categóricas
# [x] Teste de tendecia para variável ordinal
# [] Avaliando proporcional hazards
# [] 

# Load data --------------------------------------------------------------------

data(GBSG) # from shrink

# Preprocessing data -----------------------------------------------------------

data(GBSG)

gbsg <- GBSG %>%
  mutate(
    htreat = factor(htreat,
                    levels = c(0, 1),
                    labels = c("Não", "Sim")),
    menostat = factor(menostat,
                      levels = c(1, 2),
                      labels = c("Pré-menopausa",
                                 "Pós-menopausa")),
    tumgrad = factor(tumgrad,
                     levels = c(1, 2, 3),
                     labels = c("Grau I",
                                "Grau II",
                                "Grau III"))
  ) %>%
  select(
    id,
    htreat,
    age,
    menostat,
    tumsize,
    tumgrad,
    posnodal,
    prm,
    esm,
    rfst,
    cens
  )

rm(GBSG)

# Analise descritiva -----------------------------------------------------------

gbsg <- gbsg %>%
  mutate(status = factor(cens,
                         levels = c(0, 1),
                         labels = c("Censurado", "Evento")))

tab_status <- gbsg %>%
  select(-id, -cens) %>%
  tbl_summary(
    by = status,
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    missing = "ifany"
  ) %>%
  add_p() %>%   # testa diferença entre grupos
  bold_labels()

tab_status

gbsg %>%
  group_by(status) %>%
  summarise(
    age_mean = mean(age, na.rm = TRUE),
    age_sd = sd(age, na.rm = TRUE),
    
    tumsize_mean = mean(tumsize, na.rm = TRUE),
    tumsize_sd = sd(tumsize, na.rm = TRUE),
    
    posnodal_mean = mean(posnodal, na.rm = TRUE),
    posnodal_sd = sd(posnodal, na.rm = TRUE),
    
    rfst_mean = mean(rfst, na.rm = TRUE),
    rfst_sd = sd(rfst, na.rm = TRUE)
  )

gbsg$age <- scale(gbsg$age)
gbsg$tumsize <- scale(gbsg$tumsize)
gbsg$posnodal <- scale(gbsg$posnodal)
gbsg$prm <- scale(gbsg$prm)
gbsg$esm <- scale(gbsg$esm)

# Kaplan Meyer -----------------------------------------------------------------

surv_gbsg <- with(
  gbsg,
  Surv(rfst, cens)
)

## KM general

fit <- survfit(surv_gbsg ~ 1, data = gbsg)

summary(fit, times = c(365, 730, 1095)) # estimativas pomtuais

km_plot <- ggsurvplot(
  fit,
  data = gbsg,
  risk.table = TRUE,
  
  title = "Curva de Kaplan-Meier",
  xlab = "Tempo até recorrência (dias)",
  ylab = "Probabilidade de sobrevivência livre de recorrência",
  
  risk.table.title = "Número de pacientes em risco",
  risk.table.y.text = FALSE
)

ggsave(
  filename = "lista2/plot/gbsg/km.png",
  plot = km_plot$plot,
  width = 8,
  height = 6,
  dpi = 300
)

## KM hormonal tratment

fit_htreat <- survfit(surv_gbsg ~ htreat, data = gbsg)

km_htreat_plot <- ggsurvplot(
  fit_htreat,
  data = gbsg,
  pval = TRUE,
  risk.table = TRUE,
  title = "Sobrevivência livre de recorrência por tratamento hormonal",
  xlab = "Tempo (dias)",
  ylab = "Probabilidade de sobrevivência",
  legend.title = "Tratamento",
  
  risk.table.title = "Número de pacientes em risco",
  risk.table.y.text = FALSE
)

ggsave(
  filename = "lista2/plot/gbsg/km_htreat.png",
  plot = km_htreat_plot$plot,
  width = 8,
  height = 6,
  dpi = 300
)

## KM Menopause

fit_meno <- survfit(surv_gbsg ~ menostat, data = gbsg)

km_meno_plot <- ggsurvplot(
  fit_meno,
  data = gbsg,
  pval = TRUE,
  risk.table = TRUE,
  title = "Sobrevivência livre de recorrência por estado menopausal",
  xlab = "Tempo (dias)",
  ylab = "Probabilidade de sobrevivência",
  legend.title = "Estado menopausal",
  
  risk.table.title = "Número de pacientes em risco",
  risk.table.y.text = FALSE
)

ggsave(
  filename = "lista2/plot/gbsg/km_meno.png",
  plot = km_meno_plot$plot,
  width = 8,
  height = 6,
  dpi = 300
)

## Tumour degree

fit_grad <- survfit(surv_gbsg ~ tumgrad, data = gbsg)

km_grad_plot <- ggsurvplot(
  fit_grad,
  data = gbsg,
  pval = TRUE,
  risk.table = TRUE,
  title = "Sobrevivência livre de recorrência por grau tumoral",
  xlab = "Tempo (dias)",
  ylab = "Probabilidade de sobrevivência",
  legend.title = "Grau",
  
  risk.table.title = "Número de pacientes em risco",
  risk.table.y.text = FALSE
)

ggsave(
  filename = "lista2/plot/gbsg/km_grad.png",
  plot = km_meno_plot$plot,
  width = 8,
  height = 6,
  dpi = 300
)


# Testes de curvas -------------------------------------------------------------

# htreat ------------------------------------------------------------

teste_htreat_0  <- survdiff(Surv(rfst, cens) ~ htreat, data = gbsg, rho = 0)
teste_htreat_05 <- survdiff(Surv(rfst, cens) ~ htreat, data = gbsg, rho = 0.5)
teste_htreat_1  <- survdiff(Surv(rfst, cens) ~ htreat, data = gbsg, rho = 1)
teste_htreat_m1 <- survdiff(Surv(rfst, cens) ~ htreat, data = gbsg, rho = -1)

# menostat ----------------------------------------------------------

teste_meno_0  <- survdiff(Surv(rfst, cens) ~ menostat, data = gbsg, rho = 0)
teste_meno_05 <- survdiff(Surv(rfst, cens) ~ menostat, data = gbsg, rho = 0.5)
teste_meno_1  <- survdiff(Surv(rfst, cens) ~ menostat, data = gbsg, rho = 1)
teste_meno_m1 <- survdiff(Surv(rfst, cens) ~ menostat, data = gbsg, rho = -1)

# tumgrad -----------------------------------------------------------

teste_grad_0  <- survdiff(Surv(rfst, cens) ~ tumgrad, data = gbsg, rho = 0)
teste_grad_05 <- survdiff(Surv(rfst, cens) ~ tumgrad, data = gbsg, rho = 0.5)
teste_grad_1  <- survdiff(Surv(rfst, cens) ~ tumgrad, data = gbsg, rho = 1)
teste_grad_m1 <- survdiff(Surv(rfst, cens) ~ tumgrad, data = gbsg, rho = -1)

# tabela resumo -----------------------------------------------------

resultado <- tibble(
  variavel = c(
    rep("htreat", 4),
    rep("menostat", 4),
    rep("tumgrad", 4)
  ),
  rho = rep(c(0, 0.5, 1, -1), 3),
  chisq = c(
    teste_htreat_0$chisq,
    teste_htreat_05$chisq,
    teste_htreat_1$chisq,
    teste_htreat_m1$chisq,
    
    teste_meno_0$chisq,
    teste_meno_05$chisq,
    teste_meno_1$chisq,
    teste_meno_m1$chisq,
    
    teste_grad_0$chisq,
    teste_grad_05$chisq,
    teste_grad_1$chisq,
    teste_grad_m1$chisq
  ),
  gl = c(
    rep(length(teste_htreat_0$n) - 1, 4),
    rep(length(teste_meno_0$n) - 1, 4),
    rep(length(teste_grad_0$n) - 1, 4)
  )
) %>%
  mutate(
    pvalor = pchisq(chisq, df = gl, lower.tail = FALSE)
  )

resultado


# Teste de tendencia para grau tumoral

logrank <- survdiff(
  Surv(rfst,cens) ~ tumgrad,
  data = gbsg
)

u <- with(logrank, obs-exp)

V <- logrank$var

a <- c(1,2,3)

c <- crossprod(a,u)

C <- a %*% V %*% a

z <- c/sqrt(C)

2*pnorm(abs(z), lower.tail = FALSE)

# Avaliando PH -----------------------------------------------------------------

## Criando Modelo de cox geral do zero

cox0 <- coxph(
  Surv(rfst, cens) ~
    htreat +
    age +
    menostat +
    tumsize +
    tumgrad +
    posnodal +
    prm +
    esm,
  data = gbsg
)

## Avaliando residuos

# PH assumption
ph_test <- cox.zph(cox0)
print(ph_test)

teste_ph <- cox.zph(cox0)

teste_ph

if (!dir.exists("lista2/plot/gbsg")) {
  dir.create("lista2/plot/gbsg", recursive = TRUE)
}

png(
  filename = "lista2/plot/gbsg/diagnostico_ph.png",
  width = 3600,
  height = 2400,
  res = 300
)

par(
  mfrow = c(2, 3),
  mar = c(5, 5, 4, 2),
  cex.axis = 1.2,
  cex.lab = 1.3,
  cex.main = 1.4
)

# usa índice numérico (forma correta do cox.zph)
for (i in seq_len(nrow(teste_ph$table))) {
  plot(teste_ph[i],
       main = rownames(teste_ph$table)[i],
       xlab = "Tempo",
       ylab = "Schoenfeld residuals")
  abline(h = 0, col = "red")
}

par(mfrow = c(1, 1))

dev.off()

# resíduos
mart_cox0 <- residuals(cox0, type = "martingale")
residuals(cox0, type = "deviance")

# Ajuste Yang and Prentice -----------------------------------------------------

m <- ceiling(sqrt(nrow(gbsg)))

yp <- ypreg(
  Surv(rfst, cens) ~ htreat + age + menostat +
    tumsize + tumgrad + posnodal + prm + esm,
  data = gbsg,
  #dist = piecewise(m = m)
)


po <- poreg(
  Surv(rfst, cens) ~ htreat + age + menostat +
    tumsize + tumgrad + posnodal + prm + esm,
  data = gbsg,
  #dist = piecewise(m = m)
)

# comparação PO vs YP
AIC(po, yp)
anova(po, yp)

ggresiduals(yp)
ggresiduals(po)

# Selecionando váriaveis e melhorando modelo geral -----------------------------


form <- Surv(rfst, cens) ~
  htreat + age + menostat +
  tumsize + tumgrad +
  posnodal + prm + esm

yp_weib <- ypreg(form, data = gbsg, baseline = "weibull")

yp_exp <- ypreg(form, data = gbsg, baseline = "exponential")

yp_logn <- ypreg(form, data = gbsg, baseline = "lognormal")

yp_llog <- ypreg(form, data = gbsg, baseline = "loglogistic")

yp_fat  <- ypreg(form, data = gbsg, baseline = "fatigue")

## Comparação

comparacao <- data.frame(
  Modelo = c("Weibull", "Fatigue"),
  
  logLik = c(
    as.numeric(logLik(yp_weib)),
    as.numeric(logLik(yp_fat))
  ),
  
  AIC = c(
    as.numeric(AIC(yp_weib)),
    as.numeric(AIC(yp_fat))
  )
)

comparacao$DeltaAIC <- comparacao$AIC - min(comparacao$AIC)

comparacao

## Avaliando modelo de fadiga

ggresiduals(yp_fat, type = "coxsnell")
ggresiduals(yp_fat, type = "martingale")
ggresiduals(yp_fat, type = "deviance")


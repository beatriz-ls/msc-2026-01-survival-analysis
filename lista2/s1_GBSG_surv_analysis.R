#install.packages("shrink")
#install.packages(c('bellreg', 'LambertW', 'peppm', 'YPPE'))
#install.packages(file.choose(), repos = NULL, type="source") # instalando survcurv

library(shrink)
library(survcure)
library(dplyr)
library(survival)
library(survminer)
library(sobrevivencia)

# Script 1: Modelagem do dataset GBSG
# [x] Curvas de KM
# [] Testes de comparação de curvas para váriáveis categóricas
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

# Kaplan Meyer -----------------------------------------------------------------

surv_gbsg <- with(
  gbsg,
  Surv(rfst, cens)
)

## KM general

km_plot <- ggsurvplot(
  km,
  data = gbsg,
  risk.table = TRUE,
  
  title = "Curva de Kaplan-Meier",
  xlab = "Tempo até recorrência (dias)",
  ylab = "Probabilidade de sobrevivência livre de recorrência",
  
  risk.table.title = "Número de pacientes em risco",
  risk.table.y.text = FALSE
)

ggsave(
  filename = "lista2/plot/km.png",
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
  filename = "lista2/plot/km_htreat.png",
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
  filename = "lista2/plot/km_meno.png",
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
  filename = "lista2/plot/km_grad.png",
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

ggresiduals(cox0) # cox snell

ggresiduals( # martingale
  cox0,
  type = "martingale"
)

ggresiduals( # deviance
  cox0,
  type = "deviance"
)

## Teste de proporcionalidade

teste_ph <- cox.zph(
  cox0,
  transform = "km"
)

teste_ph # Modelo não funciona bem para PH



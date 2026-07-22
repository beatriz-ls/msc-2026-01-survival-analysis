#install.packages("shrink")
#install.packages(c('bellreg', 'LambertW', 'peppm', 'YPPE'))
#install.packages(file.choose(), repos = NULL, type="source") # instalando survcurv

library(dplyr)
library(gtsummary)
library(survival)
library(survminer)
library(ggplot2)
library(tibble)
library(survcure)

# Script 2: Modelagem do dataset e1690
# [x] Curvas de KM
# [x] Testes de comparação de curvas para váriáveis categóricas
# [x] Teste de tendecia para variável ordinal
# [] Ajustando um modelo de fração de cura
# [] 

# Load data --------------------------------------------------------------------

data(e1690)

# Preprocessing data -----------------------------------------------------------

data <- e1690 %>%
  mutate(
    # Tratamento
    tratamento = factor(trt,
                        levels = c(1, 2),
                        labels = c("Tratamento 1", "Tratamento 2")),
    
    # EVENTO PRINCIPAL: reincidência
    reincidencia = factor(failcens,
                          levels = c(0, 1),
                          labels = c("Censurado", "Reincidência")),
    
    # (opcional) sobrevida geral - não é o foco agora
    status_sobrevida = factor(survcens,
                              levels = c(0, 1),
                              labels = c("Censurado", "Evento")),
    
    sexo = factor(sex,
                  levels = c(1, 2),
                  labels = c("Masculino", "Feminino")),
    
    performance = factor(ps,
                         levels = c(0, 1),
                         labels = c("Bom", "Ruim")),
    
    idade = age,
    nodos = node,
    breslow = as.numeric(breslow)
  ) %>%
  select(
    tratamento,
    tempo_reincidencia = failtime,
    reincidencia,
    idade,
    nodos,
    sexo,
    performance,
    breslow
  )

rm(e1690)

# Analise descritiva -----------------------------------------------------------

tab_status <- data %>%
  tbl_summary(
    by = reincidencia,
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    missing = "ifany"
  ) %>%
  add_p() %>%
  bold_labels()

tab_status

# Médias e desvios-padrão

data %>%
  group_by(reincidencia) %>%
  summarise(
    idade_media   = mean(idade, na.rm = TRUE),
    idade_dp      = sd(idade, na.rm = TRUE),
    
    nodos_media   = mean(nodos, na.rm = TRUE),
    nodos_dp      = sd(nodos, na.rm = TRUE),
    
    breslow_media = mean(breslow, na.rm = TRUE),
    breslow_dp    = sd(breslow, na.rm = TRUE)
  )

# Padronização das variáveis contínuas

data <- data %>%
  mutate(
    idade = scale(idade)[,1],
    nodos = scale(nodos)[,1],
    breslow = scale(breslow)[,1]
  )

# Kaplan Meyer -----------------------------------------------------------------

surv_e1690 <- with(
  data,
  Surv(tempo_reincidencia,
       reincidencia == "Reincidência")
)

## KM geral

#===============================================================================
# Objeto de sobrevivência
#===============================================================================

surv_e1690 <- with(
  data,
  Surv(tempo_reincidencia,
       reincidencia == "Reincidência")
)

#===============================================================================
# Curva geral
#===============================================================================

km <- survfit(surv_e1690 ~ 1)

km_plot <- ggsurvplot(
  km,
  data = data,
  risk.table = TRUE,
  
  title = "Curva de Kaplan-Meier",
  xlab = "Tempo até recorrência",
  ylab = "Probabilidade livre de recorrência",
  
  risk.table.title = "Número de pacientes em risco",
  risk.table.y.text = FALSE
)

km_plot

ggsave(
  filename = "lista2/plot/e1690/km_global.png",
  plot = km_plot$plot,
  width = 8,
  height = 6,
  dpi = 300
)

## KM tratamento

fit_trat <- survfit(surv_e1690 ~ tratamento, data = data)

km_trat_plot <- ggsurvplot(
  fit_trat,
  data = data,
  pval = TRUE,
  risk.table = TRUE,
  
  title = "Sobrevivência livre de recorrência por tratamento",
  xlab = "Tempo",
  ylab = "Probabilidade de sobrevivência",
  
  legend.title = "Tratamento",
  
  risk.table.title = "Número de pacientes em risco",
  risk.table.y.text = FALSE
)

km_trat_plot

ggsave(
  filename = "lista2/plot/e1690/km_tratamento.png",
  plot = km_trat_plot$plot,
  width = 8,
  height = 6,
  dpi = 300
)

fit_sexo <- survfit(surv_e1690 ~ sexo, data = data)

km_sexo_plot <- ggsurvplot(
  fit_sexo,
  data = data,
  pval = TRUE,
  risk.table = TRUE,
  
  title = "Sobrevivência livre de recorrência por sexo",
  xlab = "Tempo",
  ylab = "Probabilidade de sobrevivência",
  
  legend.title = "Sexo",
  
  risk.table.title = "Número de pacientes em risco",
  risk.table.y.text = FALSE
)

km_sexo_plot

ggsave(
  filename = "lista2/plot/e1690/km_sexo.png",
  plot = km_sexo_plot$plot,
  width = 8,
  height = 6,
  dpi = 300
)

fit_ps <- survfit(surv_e1690 ~ performance, data = data)

km_ps_plot <- ggsurvplot(
  fit_ps,
  data = data,
  pval = TRUE,
  risk.table = TRUE,
  
  title = "Sobrevivência livre de recorrência por Performance Status",
  xlab = "Tempo",
  ylab = "Probabilidade de sobrevivência",
  
  legend.title = "Performance",
  
  risk.table.title = "Número de pacientes em risco",
  risk.table.y.text = FALSE
)

km_ps_plot

ggsave(
  filename = "lista2/plot/e1690/km_performance.png",
  plot = km_ps_plot$plot,
  width = 8,
  height = 6,
  dpi = 300
)

#===============================================================================
# Tratamento
#===============================================================================

teste_trat_0  <- survdiff(surv_e1690 ~ tratamento, data = data, rho = 0)
teste_trat_05 <- survdiff(surv_e1690 ~ tratamento, data = data, rho = 0.5)
teste_trat_1  <- survdiff(surv_e1690 ~ tratamento, data = data, rho = 1)
teste_trat_m1 <- survdiff(surv_e1690 ~ tratamento, data = data, rho = -1)

#===============================================================================
# Sexo
#===============================================================================

teste_sexo_0  <- survdiff(surv_e1690 ~ sexo, data = data, rho = 0)
teste_sexo_05 <- survdiff(surv_e1690 ~ sexo, data = data, rho = 0.5)
teste_sexo_1  <- survdiff(surv_e1690 ~ sexo, data = data, rho = 1)
teste_sexo_m1 <- survdiff(surv_e1690 ~ sexo, data = data, rho = -1)

#===============================================================================
# Performance
#===============================================================================

teste_ps_0  <- survdiff(surv_e1690 ~ performance, data = data, rho = 0)
teste_ps_05 <- survdiff(surv_e1690 ~ performance, data = data, rho = 0.5)
teste_ps_1  <- survdiff(surv_e1690 ~ performance, data = data, rho = 1)
teste_ps_m1 <- survdiff(surv_e1690 ~ performance, data = data, rho = -1)

resultado <- tibble(
  variavel = c(
    rep("Tratamento", 4),
    rep("Sexo", 4),
    rep("Performance", 4)
  ),
  
  rho = rep(c(0, 0.5, 1, -1), 3),
  
  chisq = c(
    teste_trat_0$chisq,
    teste_trat_05$chisq,
    teste_trat_1$chisq,
    teste_trat_m1$chisq,
    
    teste_sexo_0$chisq,
    teste_sexo_05$chisq,
    teste_sexo_1$chisq,
    teste_sexo_m1$chisq,
    
    teste_ps_0$chisq,
    teste_ps_05$chisq,
    teste_ps_1$chisq,
    teste_ps_m1$chisq
  ),
  
  gl = c(
    rep(length(teste_trat_0$n)-1,4),
    rep(length(teste_sexo_0$n)-1,4),
    rep(length(teste_ps_0$n)-1,4)
  )
) %>%
  mutate(
    pvalor = pchisq(chisq,
                    df = gl,
                    lower.tail = FALSE)
  )

resultado

# Teste de tendencia

logrank <- survdiff(
  Surv(tempo_reincidencia,
       reincidencia == "Reincidência") ~ node,
  data = e1690
)

u <- with(logrank, obs - exp)
V <- logrank$var

a <- c(1, 2, 3, 4)

c <- crossprod(a, u)
C <- a %*% V %*% a

z <- c / sqrt(C)

2 * pnorm(abs(z), lower.tail = FALSE)


# Ajuste modelo survcure -------------------------------------------------------

dados_cura <- data %>%
  mutate(
    status = ifelse(reincidencia == "Reincidência", 1, 0)
  ) %>%
  filter(tempo_reincidencia > 0)

inc_bern <- survcure(
  Surv(tempo_reincidencia, status) ~ tratamento + idade + nodos +
    sexo + performance + breslow,
  data = dados_cura,
  incidence = "bernoulli",
  latency = "weibull"
)

inc_pois <- survcure(
  Surv(tempo_reincidencia, status) ~ tratamento + idade + nodos +
    sexo + performance + breslow,
  data = dados_cura,
  incidence = "poisson",
  latency = "weibull"
)

inc_bell <- survcure(
  Surv(tempo_reincidencia, status) ~ tratamento + idade + nodos +
    sexo + performance + breslow,
  data = dados_cura,
  incidence = "bell",
  latency = "weibull"
)

inc_nb <- survcure(
  Surv(tempo_reincidencia, status) ~ tratamento + idade + nodos +
    sexo + performance + breslow,
  data = dados_cura,
  incidence = "negbin",
  latency = "weibull"
)


AIC(
  inc_bern,
  inc_pois,
  inc_bell,
  inc_nb
)

weib <- survcure(
  Surv(tempo_reincidencia, status) ~
    tratamento + idade + nodos +
    sexo + performance + breslow,
  
  data = dados_cura,
  
  incidence = "negbin",
  
  latency = "weibull",
  
  init = 0
)

pe <- survcure(
  Surv(tempo_reincidencia, status) ~
    tratamento + idade + nodos +
    sexo + performance + breslow,
  
  data = dados_cura,
  
  incidence = "negbin",
  
  latency = "pe",
  
  init = 0
)

bp <- survcure(
  Surv(tempo_reincidencia, status) ~
    tratamento + idade + nodos +
    sexo + performance + breslow,
  
  data = dados_cura,
  
  incidence = "negbin",
  
  latency = "bp",
  
  init = 0
)

AIC(weib, pe, bp)


#===============================================================================
# Estruturas de regressão para a latência Weibull
#===============================================================================

modelo_base <- survcure(
  Surv(tempo_reincidencia, status) ~
    tratamento + idade + nodos +
    sexo + performance + breslow,
  
  data = dados_cura,
  
  incidence = "negbin",
  
  latency = "weibull",
  
  init = 0
)

modelo_ph <- survcure(
  Surv(tempo_reincidencia, status) ~
    tratamento + idade + nodos +
    sexo + performance + breslow |
    tratamento + idade + nodos +
    sexo + performance + breslow,
  
  data = dados_cura,
  
  incidence = "negbin",
  
  latency = weibull("ph"),
  
  init = 0
)

modelo_po <- survcure(
  Surv(tempo_reincidencia, status) ~
    tratamento + idade + nodos +
    sexo + performance + breslow |
    tratamento + idade + nodos +
    sexo + performance + breslow,
  
  data = dados_cura,
  
  incidence = "negbin",
  
  latency = weibull("po"),
  
  init = 0
)

#===============================================================================
# Testes da razão de verossimilhança
#===============================================================================

anova(modelo_base, modelo_ph)

anova(modelo_base, modelo_po)


# Avaliando PH -----------------------------------------------------------------

#===============================================================================
# Modelo de Cox
#===============================================================================

cox_fit <- coxph(
  Surv(tempo_reincidencia,
       status) ~
    tratamento +
    idade +
    nodos +
    sexo +
    performance +
    breslow,
  data = dados_cura
)

summary(cox_fit)

#===============================================================================
# Teste da hipótese de riscos proporcionais
#===============================================================================

teste_ph <- cox.zph(cox_fit)

teste_ph


if (!dir.exists("plots")) {
  dir.create("plots")
}

png(
  filename = "lista2/plot/e1690/diagnostico_ph.png",
  width = 3600,
  height = 2400,
  res = 300
)

par(
  mfrow = c(2, 3),
  mar = c(5, 5, 4, 2),
  cex.axis = 1.4,
  cex.lab = 1.5,
  cex.main = 1.6
)

plot(teste_ph)

par(mfrow = c(1, 1))

dev.off()

#===============================================================================
# Resumo do modelo final
#===============================================================================

summary(modelo_ph)

coef(modelo_ph)

confint(modelo_ph)



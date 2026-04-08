library(survival)
library(tidyverse)

# Script X: Modelos de tempo acelerados
# - Seguindo o cap 5.1. do livro do professor
# - Exemplo 5.1 Sobrevida de pacientes com leucemia aguda: 
# - Considere os tempos de sobrevivência, em semanas, de 17 pacientes com leucemia aguda.

# Criando dataset de exemplo ----------------------------------------------------

leucemia <- data.frame(
  tempo = c(65,156,100,134,16,108,121,4,39,143,56,26,22,1,1,5,65),
  evento = rep(1,17),
  lwbc = c(3.36,2.88,3.63,3.41,3.78,4.02,4.00,4.23,3.73,3.85,3.97,4.51,4.54,5.00,5.00,4.72,5.00)  
)

exponencial1 <- survreg(Surv(tempo, evento)~lwbc, data = leucemia, dist = "exponential")
weibull1 <- survreg(Surv(tempo, evento)~lwbc, data = leucemia, dist = "weibull")

anova(exponencial1, weibull1)

# =====================================================================
# Parcial I - Econometría II (Universidad del Quindío) 
# Juan José Garcés Pineda - juanj.garcesp@uqvirtual.edu.co
# Determinantes del salario de jugadores de la NBA (base nbasal)
# Corte transversal | wage en miles de dólares
# =====================================================================

# ---- 0. Paquetes y datos --------------------------------------------
paquetes <- c("wooldridge", "tidyverse", "moments", "tseries",
              "lmtest", "sandwich", "car", "broom")
nuevos <- paquetes[!(paquetes %in% installed.packages()[, "Package"])]
if (length(nuevos)) install.packages(nuevos)
invisible(lapply(paquetes, library, character.only = TRUE))

data("nbasal")
df <- nbasal %>% drop_na(wage, points, rebounds, assists, exper, guard, forward)
# Categoría base de posición: center

# ---- (a) Especificación ---------------------------------------------
# lwage = b0 + b1*points + b2*rebounds + b3*assists + b4*exper
#         + b5*expersq + b6*guard + b7*forward + u
# Signos esperados: b1>0, b2>0, b3>0, b4>0, b5<0, b6 y b7: indeterminados
formula_m <- lwage ~ points + rebounds + assists + exper + expersq + guard + forward

# ---- (b) Análisis descriptivo con tests -----------------------------
desc <- df %>%
  select(lwage, wage, points, rebounds, assists, exper) %>%
  pivot_longer(everything(), names_to = "variable") %>%
  group_by(variable) %>%
  summarise(n = n(), media = mean(value), mediana = median(value),
            sd = sd(value), min = min(value), max = max(value),
            asimetria = skewness(value), curtosis = kurtosis(value))
print(desc)

# Normalidad: salario vs. log del salario (justifica usar logaritmo)
jarque.bera.test(df$wage)       # H0: normalidad
jarque.bera.test(df$lwage)

# Diferencia de salarios entre guards y no guards (Welch)
t.test(lwage ~ guard, data = df, var.equal = FALSE)

# ANOVA de una vía: salario por posición
df <- df %>% mutate(posicion = case_when(guard == 1 ~ "Guard",
                                         forward == 1 ~ "Forward",
                                         TRUE ~ "Center"))
summary(aov(lwage ~ posicion, data = df))

# Correlaciones
print(cor(df %>% select(lwage, points, rebounds, assists, exper)))

# Gráficos
ggplot(df, aes(wage)) + geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  labs(title = "Histograma del salario", x = "Salario (miles USD)")
ggplot(df, aes(lwage)) + geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  labs(title = "Histograma de ln(salario)", x = "ln(salario)")
ggplot(df, aes(points, lwage)) + geom_point(alpha = .4) + geom_smooth(method = "lm") +
  labs(title = "ln(salario) vs. puntos por partido")
ggplot(df, aes(posicion, lwage, fill = posicion)) + geom_boxplot() +
  labs(title = "ln(salario) por posición", x = "", y = "ln(salario)")

# ---- (c) Estimación por MCO -----------------------------------------
m <- lm(formula_m, data = df)
summary(m)
confint(m)

b <- coef(m)
cat("Un punto más por partido aumenta el salario en (%):", 100 * b["points"], "\n")
cat("Experiencia óptima (años):", -b["exper"] / (2 * b["expersq"]), "\n")

# ---- (d) Tabla ANOVA ------------------------------------------------
n <- nobs(m); k <- length(coef(m)) - 1
SCT <- sum((df$lwage - mean(df$lwage))^2)
SCR <- sum(resid(m)^2)
SCE <- SCT - SCR
CME <- SCE / k
CMR <- SCR / (n - k - 1)

anova_tab <- tibble(
  Fuente = c("Regresión", "Residuos", "Total"),
  SC     = c(SCE, SCR, SCT),
  gl     = c(k, n - k - 1, n - 1),
  CM     = c(CME, CMR, NA),
  F      = c(CME / CMR, NA, NA)
)
print(anova_tab)
cat("R2 =", SCE / SCT, "| R2 ajustado =", summary(m)$adj.r.squared, "\n")

# ---- (e) Significancia individual (t) -------------------------------
tidy(m) %>%
  mutate(decision = ifelse(p.value < 0.05, "Rechaza H0 (significativa)",
                           "No rechaza H0")) %>%
  print()
cat("t crítico (5%):", qt(0.975, n - k - 1), "\n")

# ---- (f) Significancia global (F) -----------------------------------
Fest <- summary(m)$fstatistic
cat("F =", Fest[1], "| p-valor =",
    pf(Fest[1], Fest[2], Fest[3], lower.tail = FALSE), "\n")
cat("F crítico (5%):", qf(0.95, k, n - k - 1), "\n")

# Prueba conjunta de posición (guard = forward = 0)
linearHypothesis(m, c("guard = 0", "forward = 0"))

# ---- (g) Verificación de supuestos ----------------------------------
# 1. Normalidad de los errores
jarque.bera.test(resid(m))
shapiro.test(resid(m))
qqnorm(resid(m)); qqline(resid(m))

# 2. Homocedasticidad
bptest(m)                                           # Breusch-Pagan
bptest(m, ~ fitted(m) + I(fitted(m)^2))             # White reducido
plot(fitted(m), resid(m), xlab = "Ajustados", ylab = "Residuos"); abline(h = 0, col = "red")

# 3. Multicolinealidad (exper y expersq: VIF alto por construcción)
vif(m)

# 4. Forma funcional
resettest(m, power = 2:3, type = "fitted")

# 5. Media de los errores
mean(resid(m))

# 6. Observaciones influyentes
plot(m, which = 4)
which(cooks.distance(m) > 4 / n)

# ---- Corrección si hay heterocedasticidad ---------------------------
coeftest(m, vcov = vcovHC(m, type = "HC3"))

# Diagnóstico general
par(mfrow = c(2, 2)); plot(m); par(mfrow = c(1, 1))

# =====================================================================
# Análisis exploratorio y dispersión del rendimiento
# =====================================================================

## 1. Distribución del salario en nivel (wage)
hist(df$wage, breaks = 20, col = "#AEDFF7", border = "black",
     main = "1. Distribución del Salario\nen Nivel (wage)",
     xlab = "Salario (miles de USD)", ylab = "Frecuencia")

## 2. Distribución del logaritmo del salario (lwage)
hist(df$lwage, breaks = 20, col = "#B7F3C0", border = "black",
     main = "2. Distribución del Log(Salario)\n(lwage)",
     xlab = "Logaritmo del salario", ylab = "Frecuencia")

## 3. Log(Wage) vs. Puntos
plot(df$points, df$lwage, col = "darkgreen", pch = 16,
     main = "3. Log(Wage) vs. Puntos", xlab = "Puntos por Partido",
     ylab = "Logaritmo del Salario (lwage)")
abline(lm(lwage ~ points, data = df), col = "red", lwd = 2)

## 4. Log(Wage) vs. Rebotes
plot(df$rebounds, df$lwage, col = "darkorange", pch = 16,
     main = "4. Log(Wage) vs. Rebotes", xlab = "Rebotes por Partido",
     ylab = "Logaritmo del Salario (lwage)")
abline(lm(lwage ~ rebounds, data = df), col = "red", lwd = 2)

## 5. Log(Wage) vs. Asistencias
plot(df$assists, df$lwage, col = "purple", pch = 16,
     main = "5. Log(Wage) vs. Asistencias", xlab = "Asistencias por Partido",
     ylab = "Logaritmo del Salario (lwage)")
abline(lm(lwage ~ assists, data = df), col = "red", lwd = 2)

## 6. Log(Wage) vs. Experiencia
plot(df$exper, df$lwage, col = "steelblue", pch = 16,
     main = "6. Log(Wage) vs. Experiencia", xlab = "Años de Experiencia",
     ylab = "Logaritmo del Salario (lwage)")
abline(lm(lwage ~ exper, data = df), col = "red", lwd = 2)

## 7. Comparación de retornos marginales
anios <- 0:20
ret_points <- b["points"]   * anios
ret_reb    <- b["rebounds"] * anios
ret_ast    <- b["assists"]  * anios
ret_exp    <- b["exper"] * anios + b["expersq"] * anios^2

plot(anios, ret_points, type = "l", col = "red", lwd = 2,
     ylim = range(c(ret_points, ret_reb, ret_ast, ret_exp)),
     main = "Comparación de Retornos Marginales\nal Rendimiento y la Experiencia",
     xlab = "Unidades adicionales", ylab = "Efecto acumulado sobre Log(Salario)")
lines(anios, ret_reb, col = "darkgreen", lwd = 2, lty = 2)
lines(anios, ret_ast, col = "purple", lwd = 2, lty = 3)
lines(anios, ret_exp, col = "steelblue", lwd = 2, lty = 4)
legend("topleft",
       legend = c("Puntos por partido", "Rebotes por partido",
                  "Asistencias por partido", "Experiencia (años)"),
       col = c("red", "darkgreen", "purple", "steelblue"),
       lty = c(1, 2, 3, 4), lwd = 2, bty = "n")

## 8. Histograma de residuos con curva de densidad
res <- resid(m)
hist(res, breaks = 14, freq = FALSE, col = "#4C79A8", border = "black",
     main = "8. Distribución de Residuos\n(Prueba Jarque-Bera)",
     xlab = "Residuos del modelo (u)", ylab = "Densidad")
lines(density(res), col = "darkred", lwd = 2)

## 9. Q-Q plot de residuos
qqnorm(res, main = "9. Q-Q Plot de Residuos", col = "navy", pch = 16)
qqline(res, col = "red", lwd = 2)

## 10. Residuos vs. valores ajustados
plot(fitted(m), res, col = "#6699CC", pch = 16,
     main = "10. Residuos vs. Valores Ajustados\n(Prueba Breusch-Pagan)",
     xlab = "Salario Estimado (valores ajustados)", ylab = "Residuos (u)")
abline(h = 0, col = "red", lty = 2, lwd = 2)

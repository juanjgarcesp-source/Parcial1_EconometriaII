# =====================================================================
# Parcial I - Econometría II (Universidad del Quindío)
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

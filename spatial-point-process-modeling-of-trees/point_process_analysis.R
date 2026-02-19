# ----------------------------------
# Curso: Estadística Espacial
# Alumna: Sofia Dextre
# Descripción: Este trabajo es para la PC3 del curso, contiene el análisis
#              puntual espacial de la data bei y bei.extra 
# -----------------------------------

# ---------------------- L I B R E R I A S ----------------------
install.packages("spatstat.data")
install.packages("spatstat.geom")
install.packages("spatstat.explore")
install.packages("spatstat.model")

library(spatstat.data)
library(spatstat.geom)
library(spatstat.explore)
library(spatstat.model)

# ----------------- LEER DATA ----------------------------------
data(bei)

bei.extra <- spatstat.data::bei.extra 
str(bei.extra)

# ------- EXPLORACIÓN INICIAL: DISTRIBUCIONES -----------------

plot(bei, main = "Distribución espacial de los árboles")
#covariabels - ver como la elevacion cambia
plot(bei.extra$elev, main = "Elevación del terreno")
plot(bei.extra$grad, main = "Gradiente de elevación (pendiente)")

# -------- TEST DE ALEATORIEDAD ESPACIAL --------------------------------
# Test de cuadrantes
# Dividir el área en 6x6 cuadrantes
test_cuadrante <- quadrat.test(bei, nx = 5, ny = 5)
print(test_cuadrante)
plot(test_cuadrante)

# duncion k de ripley
K <- Kest(bei)
plot(K, main = "Función K de Ripley")
# agregando bandas de simulación
set.seed(123)  # Para reproducibilidad
K_env <- envelope(bei, Kest, nsim = 99)
plot(K_env, main = "Función K de Ripley con bandas de simulación")

#--------- ESTIMACIÓN DE LA FUNCIÓN DE INTENSIDAD ------------------------
# Estimación de intensidad con kernel smoothing
lambda_kernel <- density(bei)
# Visualizar el mapa de intensidad
plot(lambda_kernel, main = "Intensidad estimada (kernel smoothing)")

# Estimación no paramétrica
# Modelo inhomogéneo con elevación y gradiente
modelo_log <- ppm(bei ~ elev + grad, covariates = bei.extra)
summary(modelo_log)
# visualizacion de la intensidad estimada
plot(predict(modelo_log), main = "Intensidad estimada - modelo paramétrico")


#---------  DIAGNÓSTICO DEL MODELO ------------------------
#revisar residuos espaciales o diagnóstico del modelo
diagnose.ppm(modelo_log)

#evaluar la bondad de ajuste
Kres <- Kres(modelo_log)
env_Kres <- envelope(modelo_log, Kres)
plot(env_Kres)

#Compación entre modelos:
modelo_1 <- ppm(bei ~ elev, covariates = bei.extra)
modelo_2 <- ppm(bei ~ grad, covariates = bei.extra)
anova(modelo_1, modelo_log, test = "Chi")

#validación cruzada o simulación
simulacion <- rpoispp(modelo_log)
plot(simulacion)

# ----------- INTERACCIONES Y TRANSFORMACIONES PARA EL MODELO----------------
# MODELO CON INTERACCIÓN
modelo_int <- ppm(bei ~ elev * grad, covariates = bei.extra)
summary(modelo_int)
# MODELO CON TRANSFORMACIÓN
modelo_transf <- ppm(bei ~ I(elev^2) + sqrt(grad), covariates = bei.extra)
summary(modelo_transf)
# MODELO CON TRANSFORMACIÓN + INTERACCIÓN
modelo_full <- ppm(bei ~ I(elev^2) * sqrt(grad), covariates = bei.extra)
summary(modelo_full)

# --------- comprar los 3 modelos ------------------------------------------
#AIC
AIC(modelo_int)
AIC(modelo_transf)
AIC(modelo_full)
#comparar con diagnose ppm
par(mfrow = c(1, 3))  # Para ver los 3 plots en una fila
diagnose.ppm(modelo_int, main = "Modelo Interacción")
diagnose.ppm(modelo_transf, main = "Modelo Transformado")
diagnose.ppm(modelo_full, main = "Modelo Full")

#comparar predicciones de intensidad
par(mfrow = c(1, 3))
plot(predict(modelo_int), main = "Intensidad - Interacción")
plot(predict(modelo_transf), main = "Intensidad - Transformación")
plot(predict(modelo_full), main = "Intensidad - Full")

#loglikehood
logLik(modelo_int)
logLik(modelo_transf)
logLik(modelo_full)

# ------------- SIMULACIÓN PARA EL MODELO FULL ----------------------
# Simular patrón desde el modelo ajustado
set.seed(123)  # Para reproducibilidad
# Simulate Toma el modelo ajustado (modelo_full) y genera una simulación 
#del proceso de Poisson inhomogéneo definido por ese modelo.
#nsim=1 porque solo queremosuna simulacion
simulado <- simulate(modelo_full, nsim = 1) 

# Comparar con el patrón real
par(mfrow = c(1, 2))
plot(bei, main = "Patrón observado")
plot(simulado[[1]], main = "Simulación del modelo final")
par(mfrow = c(1, 1))





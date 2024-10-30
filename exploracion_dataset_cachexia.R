# Archivo de RStudio de la exploracion_cachexia.R

#Cargamos las librerías necesarias
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

library(ggplot2)
library(dplyr)

#Paso 1: Cargamos el dataset
dataset <- dataset <- read.csv("https://raw.githubusercontent.com/aramonarr/RAMON-ARRUFAT-ALBERT-PEC1/main/human_cachexia.csv", row.names = 1)

#Seleccionamos solo las columnas numéricas
numeric_data <- dataset[, sapply(dataset, is.numeric)]

#Una vez cargado el dataset, creo un objeto SummarizedExperiment utilizando los datos del estudio y sus metadatos:
  
library(SummarizedExperiment)
colData <- DataFrame(SampleID = colnames(numeric_data), Condition = "CondiciónEjemplo")
colData
rowData <- DataFrame(FeatureID = rownames(numeric_data), Description = "DescripciónEjemplo")
se <- SummarizedExperiment(assays = list(counts = as.matrix(numeric_data)), rowData = rowData, colData = colData)
se
#Paso 2: Resumen estadístico
summary_statistics <- summary(numeric_data)
print(summary_statistics)

#Paso 3: Realizamos la transformación logarítmica para reducir el sesgo. También evitar valores cero o negativos sumando 1 antes de aplicar el logaritmo
transformed_data <- log(numeric_data + 1)

#Paso 4: Análisis de Componentes Principales (PCA)
pca <- prcomp(transformed_data, scale. = TRUE)
pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], Condition = dataset$Muscle.loss)

# Gráfico de PCA
ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point() +
  labs(title = "PCA del Dataset Human Cachexia", x = "Componente Principal 1", y = "Componente Principal 2") +
  theme_minimal()

# Paso 5: Histogramas de metabolitos seleccionados
selected_vars <- c("X1.6.Anhydro.beta.D.glucose", "X2.Aminobutyrate", "Citrate", "Lactate")

# Creamos histogramas para cada variable seleccionada
par(mfrow = c(2, 2)) # Configuramos el área de gráficos para 4 histogramas
for (var in selected_vars) {
  hist(numeric_data[[var]], main = paste("Distribución de", var), xlab = var, col = "skyblue", breaks = 15)
}
par(mfrow = c(1, 1)) # Resetear el área de gráficos

#Paso 6: Análisis de correlación
correlation_matrix <- cor(transformed_data)
print(correlation_matrix)

#Paso 7: Gráficos de caja para metabolitos seleccionados
install.packages("tidyr")
library(tidyr)
numeric_data_long <- numeric_data %>%
  mutate(Condition = dataset$Muscle.loss) %>%
  pivot_longer(cols = -Condition, names_to = "Metabolite", values_to = "Value")

#Gráfico de caja para algunos metabolitos
ggplot(numeric_data_long %>% filter(Metabolite %in% selected_vars), aes(x = Metabolite, y = Value, fill = Condition)) +
  geom_boxplot() +
  labs(title = "Distribución de Metabolitos Seleccionados", x = "Metabolito", y = "Concentración") +
  theme_minimal()


# dvsd
clase de diversidad
---
title: "Análisis de Diversidad de Especies por Montserrat Cervantes Espinoza" 
"Montserrat Cervantes Espinoza"
output: html_document
date: "2024-09-22"
---
# Instalar y cargar los paquetes necesarios
```{r setup, include=TRUE, warning=FALSE}

library(vegan)
library(betapart)
library(ggplot2)
```
## Establecer una semilla para reproducibilidad
```{r, include=TRUE, warning=FALSE}
set.seed(123)

abundancia <- matrix(sample(0:20, 30, replace = TRUE), nrow = 5, ncol = 6)
rownames(abundancia) <- paste("Sitio", 1:5)
colnames(abundancia) <- paste("Especie", 1:6)
```
# 3.1. Riqueza
```{r , include=TRUE, warning=FALSE}
riqueza <- rowSums(abundancia > 0)
print("Riqueza de especies por sitio:")
print(riqueza)
```
# 3.2. Abundancia
## Abundancia Absoluta
```{r , include=TRUE, warning=FALSE}
abundancia_absoluta <- rowSums(abundancia)
print("Abundancia absoluta por sitio:")
print(abundancia_absoluta)
```
## Abundancia Relativa
```{r , include=TRUE, warning=FALSE}
abundancia_relativa <- prop.table(abundancia, 1) # Por filas
print("Abundancia relativa por sitio:")
print(abundancia_relativa)
```
# 3.4. Equitatividad
## Índice de Shannon
```{r , include=TRUE, warning=FALSE}
equitatividad_shannon <- diversity(abundancia)
print("Índice de equitatividad de Shannon:")
print(equitatividad_shannon)
```
## Números de Hill
```{r , include=TRUE, warning=FALSE}

hill_numbers <- renyi(abundancia, c(0, 1, 2))
print("Números de Hill (q=0,1,2):")
print(hill_numbers)
```
# Curvas de acumulación

```{r , include=TRUE}
abundancia <- matrix(c(0, 1, 0, 3, 0, 0,
                       0, 0, 0, 0, 1, 2,
                       5, 0, 0, 0, 0, 0,
                       0, 0, 2, 3, 0, 4,
                       0, 1, 1, 0, 0, 0), 
                     nrow = 5, byrow = TRUE)
rownames(abundancia) <- paste("Sitio", 1:5)
colnames(abundancia) <- paste("Especie", 1:6)
```
## Análisis de Rarefacción
```{r , include=TRUE, warning=FALSE}
abn3 <- abundancia
s <- specnumber(abn3)  # Riqueza observada por sitio
raremax <- min(rowSums(abn3))  # Máximo número de individuos para rarefacción
srare <- rarefy(abn3, raremax, se = FALSE)  # Rarefacción
dr <- drarefy(abn3, raremax)  # Rarefacción con desviaciones estándar
```
### Preparar los datos para la gráfica
```{r , include=TRUE, warning=FALSE}
rare_df <- data.frame(
  Observado = s,
  Rarefaccionado = srare
)
```
## Gráfica de Rarefacción
## Gráfica de la curva de rarefacción
```{r , include=TRUE, warning=FALSE, echo=FALSE}
rarecurve(abn3, sample = raremax, col = "purple", cex = 0.6, add = TRUE)
```
# Analisis de diversidad Beta
##matrix de abundancia
```{r , include=TRUE, warning=FALSE}
abundancia <- matrix(c(0, 1, 0, 3, 0, 0,
                       0, 0, 0, 0, 1, 2,
                       5, 0, 0, 0, 0, 0,
                       0, 0, 2, 3, 0, 4,
                       0, 1, 1, 0, 0, 0), 
                     nrow = 5, byrow = TRUE)
rownames(abundancia) <- paste("Sitio", 1:5)
colnames(abundancia) <- paste("Especie", 1:6)
```
### Imprimir la matriz de abundancia
```{r , include=TRUE, warning=FALSE}
print("Matriz de Abundancia:")
print(abundancia)
```
# Crear una matriz de presencia-ausencia
### Importante considerar con que datos se esta haciendo la diversidad beta preferentemente usar datos de presencia ausencia.
```{r , include=FALSE, warning=FALSE}
presencia_ausencia <- ifelse(abundancia > 0, 1, 0)
print("Matriz de Presencia-Ausencia:")
print(presencia_ausencia)
```
# Convertir la matriz de presencia-ausencia en un objeto de betapart (paquete)
```{r , include=FALSE, warning=FALSE}
ceram.s.core <- betapart.core(presencia_ausencia)

# Medidas de diversidad beta
ceram.s.multi <- beta.multi(ceram.s.core)
print("Medidas de Diversidad Beta:")
print(ceram.s.multi)

# Muestreo a través de sitios iguales
ceram.s.samp <- beta.sample(ceram.s.core, index.family = "sor", sites = 5, samples = 150)
print("Muestreo de Diversidad Beta:")
print(ceram.s.samp)

# Distribuciones de componentes de diversidad beta
dist.s <- ceram.s.samp$sampled.values
```

```{r , include=TRUE, warning=FALSE, echo=TRUE}
# Gráficas de densidad
x11()
par(mar = c(10, 10, 10, 10), xpd = NA, cex = 0.5)

# Gráfica de la densidad de diversidad beta
plot(density(dist.s$beta.SOR), xlim = c(0, 2), ylim = c(0, 35), 
     xlab = "Diversidad Beta", main = "", lwd = 3)
lines(density(dist.s$beta.SNE), lty = 1, lwd = 2, col = "red")
lines(density(dist.s$beta.SIM), lty = 2, lwd = 2, col = "grey60")

# Gráfica alternativa
plot(density(dist.s$beta.SOR), xlim = c(0, 0.8), ylim = c(0, 19), 
     xlab = "Diversidad Beta", main = "", lwd = 3)
lines(density(dist.s$beta.SNE), lty = 1, lwd = 2)
lines(density(dist.s$beta.SIM), lty = 2, lwd = 2)
```

# Abundancia por sitio

```{r , include=TRUE, warning=FALSE, echo=TRUE}
abundancia <- matrix(c(0, 1, 0, 3, 0, 0,
                       0, 0, 0, 0, 1, 2,
                       5, 0, 0, 0, 0, 0,
                       0, 0, 2, 3, 0, 4,
                       0, 1, 1, 0, 0, 0), 
                     nrow = 5, byrow = TRUE)
rownames(abundancia) <- paste("Sitio", 1:5)
colnames(abundancia) <- paste("Especie", 1:6)

# Calcular la abundancia total por sitio
abundancia_total <- rowSums(abundancia)
abundancia_grafico <- data.frame(Sitio = rownames(abundancia), Abundancia = abundancia_total)

# Calcular la desviación estándar y error
abundancia_sd <- apply(abundancia, 1, sd)  # Desviación estándar para cada sitio
n <- ncol(abundancia)  # Número de especies
error <- qnorm(0.975) * (abundancia_sd / sqrt(n))  # Error estándar para el IC del 95%
```

```{r , include=TRUE, warning=FALSE, echo=TRUE}
ggplot(abundancia_grafico, aes(x = Sitio, y = Abundancia)) +
  geom_bar(stat = "identity", fill = "lightgreen") +
  geom_errorbar(aes(ymin = Abundancia - error, ymax = Abundancia + error), width = 0.2, color = "red") +
  labs(title = "Abundancia por sitio", y = "Abundancia", x = "Sitios") +
  theme_minimal()
```
# 5. Análisis de la Diversidad Beta con vegan
```{r , include=TRUE, warning=FALSE, echo=TRUE}
beta_diversidad <- vegdist(abundancia)
```
# Descomposición de la Beta Diversidad
```{r , include=TRUE, warning=FALSE, echo=TRUE}
anidamiento_recambio <- betadisper(beta_diversidad, group = rownames(abundancia))
```
# Escalado Multidimensional (NMDS)
```{r , include=FALSE, warning=FALSE, echo=TRUE}
mds <- metaMDS(beta_diversidad, try = 100, trymax = 500)
```
# Gráfica de NMDS
```{r , include=TRUE, warning=FALSE, echo=TRUE}
nmds_plot <- data.frame(mds$points)
ggplot(nmds_plot, aes(x = MDS1, y = MDS2)) +
  geom_point(color = "blue", size = 3) +
  geom_text(aes(label = rownames(abundancia)), vjust = -1, color = "red") +
  labs(title = "Escalado Multidimensional (NMDS)", x = "MDS1", y = "MDS2") +
  theme_minimal()
```
# Análisis de Coordenadas Principales (PCA)
```{r , include=FALSE, warning=FALSE}
abundancia <- matrix(c(0, 1, 0, 3, 0, 0,
                       0, 0, 0, 0, 1, 2,
                       5, 0, 0, 0, 0, 0,
                       0, 0, 2, 3, 0, 4,
                       0, 1, 1, 0, 0, 0), 
                     nrow = 5, byrow = TRUE)
rownames(abundancia) <- paste("Sitio", 1:5)
colnames(abundancia) <- paste("Especie", 1:6)
pca_resultados <- prcomp(abundancia, scale. = TRUE)
pca_points <- data.frame(pca_resultados$x)
```

# Gráfica de PCA con ejes y vectores
```{r , include=TRUE, warning=FALSE, echo=TRUE}
ggplot(pca_points, aes(x = PC1, y = PC2)) +
  geom_point(color = "green", size = 3) +
  geom_text(aes(label = rownames(abundancia)), vjust = -1, color = "blue") +
  geom_segment(aes(xend = PC1 * 0.5, yend = PC2 * 0.5), 
               arrow = arrow(length = unit(0.2, "inches")), color = "orange") +
  labs(title = "Análisis de Coordenadas Principales (PCA)", 
       x = "PC1", y = "PC2") +
  theme_minimal()
```

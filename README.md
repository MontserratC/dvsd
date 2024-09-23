---
title: "Análisis de Diversidad de Especies por Montserrat Cervantes Espinoza"
output:
  html_document: default
  pdf_document: default
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
```
# Crear una matrix con la siguiente estructura
### sitio = filas
### especies = columnas
```{r, include=TRUE, warning=FALSE}
abundancia <- matrix(c(0, 1, 0, 3, 0, 0,
                       0, 0, 0, 0, 1, 2,
                       5, 0, 0, 0, 0, 0,
                       0, 0, 2, 3, 0, 4,
                       0, 1, 1, 0, 0, 0), 
                     nrow = 5, byrow = TRUE)
rownames(abundancia) <- paste("Sitio", 1:5)
colnames(abundancia) <- paste("Especie", 1:6)
```
#  Riqueza
## simplemente la suma de filas (numero de especies)
```{r , include=TRUE, warning=FALSE}
riqueza <- rowSums(abundancia > 0)
print("Riqueza de especies por sitio:")
print(riqueza)
```
#  Abundancia
## Abundancia Absoluta (suma de abundancias)
```{r , include=TRUE, warning=FALSE}
abundancia_absoluta <- rowSums(abundancia)
print("Abundancia absoluta por sitio:")
print(abundancia_absoluta)
```
## Abundancia Relativa (suma abundancias entre uno)
```{r , include=TRUE, warning=FALSE}
abundancia_relativa <- prop.table(abundancia, 1) # Por filas
print("Abundancia relativa por sitio")
print(abundancia_relativa)
```
#  Equitatividad
## Índice de Shannon
```{r , include=TRUE, warning=FALSE}
equitatividad_shannon <- diversity(abundancia)
print("Índice de equitatividad de Shannon:")
print(equitatividad_shannon)
```
## Números de Hill
## recordar 
**q=0** Equivale a la riqueza de especies (S), es decir, al número total de especies presentes.

**q=1** Equivale al exponencial del índice de Shannon (exp(H')), que considera tanto la riqueza como la equitatividad.

**q=2** Equivale al inverso del índice de Simpson (1/D), que pondera más las especies dominantes.

```{r , include=TRUE, warning=FALSE}

hill_numbers <- renyi(abundancia, c(0, 1, 2))
print("Números de Hill (q=0,1,2):")
print(hill_numbers)
```
# Curvas de acumulación

```{r , include=TRUE}
#matrix de abundancia
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
### se enfoca en comparar la diversidad en diferentes tamaños de muestra.
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
# **Analisis de diversidad Beta**
## matrix de abundancia
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
### Convertir la matriz de presencia-ausencia en un objeto de betapart (paquete)
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
### Diferencia de riqueza de especies (anidamiento) **(βSNE)**. 
### Recambio o reemplazo de la disimilitud **(βSIM)** 
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

# Abundancia por sitio comparacion simple analisis exploratorio

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

# Curva Rango-abundancia
```{r , include=TRUE, warning=FALSE, echo=TRUE}
abundancia <- matrix(sample(1:50, 30 * 5, replace = TRUE), ncol = 5)
rownames(abundancia) <- paste("Sitio", 1:30, sep = "")
colnames(abundancia) <- paste("Especie", 1:5, sep = "")

# Calcular el total de abundancia por especie (suma de filas)
abundancia_total <- colSums(abundancia)

# Ordenar las especies por su abundancia de mayor a menor
abundancia_ordenada <- sort(abundancia_total, decreasing = TRUE)

# Crear un data frame para el gráfico de rango-abundancia
rango_abundancia <- data.frame(
  Rango = seq_along(abundancia_ordenada),
  Abundancia = abundancia_ordenada,
  Especie = names(abundancia_ordenada)
)

# Graficar rango-abundancia usando ggplot2
ggplot(rango_abundancia, aes(x = Rango, y = Abundancia)) +
  geom_line(size = 1, color = "blue") +
  geom_point(size = 3, color = "red") +
  geom_text(aes(label = Especie), vjust = -0.5, size = 3, angle = 45) +  # Etiquetas de especies
  labs(title = "Curva de Rango-Abundancia", x = "Rango de Especies", y = "Abundancia") +
  theme_minimal()
```

 Análisis de la Diversidad Beta con vegan
```{r , include=TRUE, warning=FALSE, echo=TRUE}
beta_diversidad <- vegdist(abundancia)
```
# Descomposición de la Beta Diversidad
```{r , include=TRUE, warning=FALSE, echo=TRUE}
anidamiento_recambio <- betadisper(beta_diversidad, group = rownames(abundancia))
```

```{r , include=TRUE, warning=FALSE, echo=TRUE}
library(vegan)
library(ggplot2)
library(reshape2)

abundancia <- matrix(c(5, 3, 0, 2,
                       0, 1, 4, 3,
                       2, 0, 5, 1,
                       1, 3, 0, 0),
                     nrow = 4, byrow = TRUE)
rownames(abundancia) <- paste("Sitio", 1:4)
colnames(abundancia) <- paste("Bigote", 1:4)

# Convertir la matriz de abundancia a presencia/ausencia
presencia_ausencia <- ifelse(abundancia > 0, 1, 0)

# Calcular la diversidad beta usando la distancia de Jaccard
beta_diversidad <- vegdist(presencia_ausencia, method = "jaccard")

# Convertir a data frame
beta_df <- as.data.frame(as.matrix(beta_diversidad))

# Agregar los nombres de los sitios como columnas
beta_df$Sitio1 <- rownames(beta_df)

# Transformar a formato largo
beta_long <- melt(beta_df, id.vars = "Sitio1", variable.name = "Sitio2", value.name = "Diversidad_Beta")

# Filtrar para evitar comparaciones consigo mismo
beta_long <- beta_long[beta_long$Sitio1 != beta_long$Sitio2, ]

# Graficar el boxplot de la diversidad beta
ggplot(beta_long, aes(x = Sitio1, y = Diversidad_Beta)) +
  geom_boxplot(fill = "skyblue", outlier.color = "red") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "darkblue", 
               position = position_dodge(0.75)) +  # Promedio
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, 
               color = "darkblue", position = position_dodge(0.75)) +  # Intervalos de confianza
  labs(title = "Boxplot de Diversidad Beta",
       x = "Sitios",
       y = "Diversidad Beta (Jaccard)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

```{r , include=TRUE, warning=FALSE, echo=TRUE}
data <- matrix(c(1, 0, 1, 0, 
                 1, 1, 0, 1, 
                 0, 1, 1, 0, 
                 1, 0, 0, 1), 
               nrow = 4, byrow = TRUE)

# Calcular los componentes de beta diversidad
beta <- beta.pair(data, index.family = "jaccard")

# Resultados
beta.sim <- beta$beta.sim  # Reemplazo de especies
beta.sne <- beta$beta.sne  # Diferencia en riqueza de especies
beta.sor <- beta$beta.sor  # Diversidad beta total

# Crear un dataframe con los resultados de beta diversidad
df <- data.frame(
  Grupo = c("anfibios", "reptiles", "mamíferos", "aves"),
  Reemplazo = c(0.7, 0.5, 0.6, 0.2), # Valores simulados, reemplázalos con tus datos
  Riqueza = c(0, 0, 0, 0.3)          # Valores simulados, reemplázalos con tus datos
)

# Instala ggplot2 si no lo tienes
if(!require(ggplot2)) install.packages("ggplot2")

# Cargar ggplot2
library(ggplot2)

# Crear la gráfica
ggplot(df, aes(x = Grupo)) +
  geom_bar(aes(y = Reemplazo + Riqueza, fill = "Anidamiento"), stat = "identity", color = "black") +
  geom_bar(aes(y = Reemplazo, fill = "Reemplazo"), stat = "identity", color = "black") +
  scale_fill_manual(name = "Componentes", values = c("Reemplazo" = "gray", "Anidamiento" = "purple")) +
  labs(y = "Diferencias en composición de especies", x = "", 
       title = "Diferencias en composición de especies por grupo o gremio") +
  theme_minimal()
```

# Análisis de la Diversidad Beta NMDS
```{r , include=TRUE, warning=FALSE, echo=TRUE}
abundancia <- matrix(c(5, 0, 3, 0, 0, 2, 5, 1, 0, 0, 4, 2, 1, 0, 3, 6), 
                     nrow = 4, byrow = TRUE)
rownames(abundancia) <- c("Sitio1", "Sitio2", "Sitio3", "Sitio4")
```

```{r , include=FALSE, warning=FALSE, echo=TRUE}
# Cálculo de la diversidad beta usando vegdist
beta_diversidad <- vegdist(abundancia)

# NMDS
mds <- metaMDS(beta_diversidad, try = 100, trymax = 500)

# Convertir los resultados en un data frame
nmds_plot <- data.frame(mds$points)

# Crear un factor de los grupos de sitios/meses para representar en la gráfica
nmds_plot$Mes <- factor(c("Septiembre", "Noviembre", "Enero", "Marzo"))
nmds_plot$Hora <- factor(c("Diurno", "Nocturno", "Diurno", "Crepuscular"))

# Simular tres conjuntos adicionales de puntos con pequeñas variaciones
set.seed(123)  # Para reproducibilidad

# Función para generar puntos adicionales con ruido
agregar_puntos <- function(df, veces) {
  for (i in 1:veces) {
    puntos_extra <- df
    puntos_extra$MDS1 <- puntos_extra$MDS1 + rnorm(nrow(puntos_extra), sd = 0.01)
    puntos_extra$MDS2 <- puntos_extra$MDS2 + rnorm(nrow(puntos_extra), sd = 0.01)
    df <- rbind(df, puntos_extra)
  }
  return(df)
}

# Agregar el triple de puntos
nmds_plot <- agregar_puntos(nmds_plot, 3)
```

```{r , include=TRUE, warning=FALSE, echo=TRUE}
# Graficar el NMDS con los puntos adicionales
ggplot(nmds_plot, aes(x = MDS1, y = MDS2, shape = Mes, color = Hora)) +
  geom_point(size = 4) + 
  stat_ellipse(aes(group = Mes), linetype = "dotted") +  # Elipses por grupo (Mes)
  scale_shape_manual(values = c(21, 22, 23, 24, 25)) +  # Diferentes formas para los meses
  scale_color_manual(values = c("black", "gray", "blue", "red")) +  # Colores para las horas
  labs(title = "Escalado Multidimensional No-Métrico (NMDS)", 
       x = "Dimensión 1", y = "Dimensión 2") +
  theme_minimal() +
  theme(legend.position = "right")

```


# Gráfica de PCA con ejes y vectores

```{r , include=TRUE, warning=FALSE, echo=TRUE}
# Análisis de Coordenadas Principales (PCA)
set.seed(123)
abundancia <- matrix(sample(1:20, 30 * 5, replace = TRUE), ncol = 5)
rownames(abundancia) <- paste("Sitio", 1:30, sep = "")
colnames(abundancia) <- paste("Especie", 1:5, sep = "")

# Añadir una variable para los grupos
grupos <- factor(rep(c("Cardonal", "Espinar", "Matorral desértico"), each = 10))

# Cálculo de disimilitud Bray-Curtis y el PCA (o escalado multidimensional)
disimilitud <- vegdist(abundancia, method = "bray")
pca_resultados <- cmdscale(disimilitud, eig = TRUE)

# Convertir los resultados en un data frame
pca_points <- data.frame(pca_resultados$points)
colnames(pca_points) <- c("PC1", "PC2")
pca_points$Sitio <- rownames(abundancia)
pca_points$Grupo <- grupos

# Calcular los centroides para cada grupo
centroides <- aggregate(cbind(PC1, PC2) ~ Grupo, data = pca_points, FUN = mean)

# Combinar los puntos con los centroides
pca_combined <- merge(pca_points, centroides, by = "Grupo", suffixes = c("", "_centroide"))

# Graficar el PCA con puntos, centroides y elipses
ggplot(pca_combined, aes(x = PC1, y = PC2, color = Grupo, shape = Grupo)) +
  geom_point(size = 3) +  # Puntos por grupo
  geom_point(aes(x = PC1_centroide, y = PC2_centroide), shape = 4, size = 4, color = "red") +  # Centroides
  geom_segment(aes(xend = PC1_centroide, yend = PC2_centroide), linetype = "dotted", color = "black") +  # Líneas a # # centroides
  stat_ellipse(aes(group = Grupo), linetype = "dashed") +  # Elipses por grupo
  labs(title = "Análisis de Coordenadas Principales (PCA) con Grupos", 
       x = "Eje 1", y = "Eje 2") +
  theme_minimal() +
  theme(legend.position = "right")
```
<span>https://github.com/MontserratC/dvsd/blob/440e5dffa0e2f3d77c01e40e5c4737b8496fedbb/coordenadasPCA.jpeg</span><span>)</span>

coordenadasPCA.jpeg

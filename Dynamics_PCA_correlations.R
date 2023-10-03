library(ggcorrplot)
library(FactoMineR)
library(heatmaply)
library(devtools)
library(factoextra)
library(corrplot)
library(ggplot)



datos <- read.table("10-20kpc_dynamic_data_1step.csv", header = TRUE,
                    sep = ",",row.names = 1)
inc_ellipsoid <- read.table("inclination_ellipsoid_10_20_1step.csv", header = TRUE,
                    sep = ",",row.names = 1)

gas_inclination <- read.table("gas_inclination.txt", header = FALSE)
head(gas_inclination)
#Adding the extra columns to the main data frame

datos$ellipsoid <- inc_ellipsoid$Inc_20
datos$gas_inc<- gas_inclination$V1

#Dropping the lookback time
data_wo_lookback <- subset(datos, select = -c(Lookback))
summary(data_wo_lookback)

#Scale data
scaled_data <- scale(data_wo_lookback)
print(ncol(scaled_data))
matriz_correlacion <- cor(scaled_data)

#Correlation plot
png(file="Figures/correlation_map.png", width=600, height= 500)
ggcorrplot(matriz_correlacion, type = "lower", lab = TRUE)
dev.off()


#PCA Analysis
data_pca <- prcomp(scaled_data)
summary(data_pca)

res.pca <- PCA(scaled_data, graph = FALSE)

print(res.pca)


#Evaluation of PCA
var <- get_pca_var(res.pca)
head(var$coord)

corrplot(var$cos2, is.corr=FALSE)
corrplot(var$contrib, is.corr=FALSE)

#Scree plot of explained variance
png(file="Figures/PCA_variance_scree_plot.png", width=600, height= 500)
fviz_eig(res.pca, addlabels=TRUE)
dev.off()


fviz_cos2(res.pca, choice="var", axes=1)
fviz_cos2(res.pca, choice="var", axes=2)

fviz_pca_var(res.pca, col.var="cos2", gradient.cols= c("black", "orange", "green"),
             repel=TRUE)

fviz_pca_var(res.pca, col.var="cos2", axes=c(3,4),
             gradient.cols= c("black", "orange", "green"),
             repel=TRUE)


#Contribution plots
fviz_contrib(res.pca, choice = "var",axes=1, top=1/ncol(scaled_data)*100)
fviz_contrib(res.pca, choice = "var",axes=2,  top=1/ncol(scaled_data)*100)
fviz_contrib(res.pca, choice = "var",axes=3, top=1/ncol(scaled_data)*100)
fviz_contrib(res.pca, choice = "var",axes=4, top=1/ncol(scaled_data)*100)



#K-Means for clustering of the variables
set.seed(12345)

#Optimization of the number of clusters
wss <- numeric(nrow(var$coord))
for (i in 1:(nrow(var$coord)-1)){
  print(i)
  kmeans_model <- kmeans(var$coord, centers=i)
  wss[i] <- kmeans_model$tot.withinss
}

k_values <- 1:ncol(scaled_data) 
wss_values <- wss 

data_kmeans <- data.frame(K = k_values, WSS = wss_values)
png(file="Figures/K-means_WSS.png", width=600, height= 400)
ggplot(data_kmeans, aes(x = K, y = WSS)) +
  geom_line() + 
  geom_point(size = 3) +  
  labs(x = "Clusters", y = "WSS", title = "k-means of PCA") +
  scale_x_continuous(breaks = seq(min(k_values), max(k_values), by = 1))
  theme_gray() 
dev.off()

  
#Grouping by clusters
res.k_means <- kmeans(var$coord, centers=6, nstart=25)
grp <- as.factor(res.k_means$cluster)


png(file="Figures/clustering_variables.png", width=650, height= 550)
fviz_pca_var(res.pca, col.var=grp, palette="viridis",
             legend.title="Cluster")
dev.off()
fviz_pca_var(res.pca, axes=c(3,4),col.var=grp, palette="viridis",
             legend.title="Cluster")


#Some interpretation of results...

#Not surprising that it groups the gas inclination + acceleration by gas.
#This further indicates that the most acceleration by gas is due to its inclination

#It also connects the inner ellipsoid inclination with the acceleration by dark matter

#The bending seems more independent from the acceleration by dark matter, further indicating 
#that it is not the sole perturber.

#Isolated Satellites (DM_out). Their acceleration are more instant, they can provide initial kick
#but are note strongly correlated with the bending
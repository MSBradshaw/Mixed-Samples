library(readr)
library(tibble)
library(ggplot2)

data <- read_csv('MixedSamples/BF_mixtures_quantilenorm_matrix.csv')
data[is.na(data)] <- 0

names <- c("Accessions","Sequence","BloodYoung1",
           "SA7","SA8","SA9",
           "B_SA_1.1","B_SA_1.100","B_SA_100.1",
           "B_SE_1.1","B_SE_1.100","B_SE_100.1",
           "Blood_WM24","Blood_WM26",
           "SE_SA_1.1","SE_SA_1.100","SE_SA_100.1",
           "SE1005","SE1022","SE1027")

groups <- c("Accessions","Sequence","Blood",
           "SA","SA","SA",
           "B_SA_1.1","B_SA_1.100","B_SA_100.1",
           "B_SE_1.1","B_SE_1.100","B_SE_100.1",
           "Blood","Blood",
           "SE_SA_1.1","SE_SA_1.100","SE_SA_100.1",
           "SE","SE","SE")

colnames(data) <- names

m <- as.matrix(data[,3:ncol(data)])
colnames(m) <- NULL

tm <- t(m)

tib <- as.tibble(tm)

pca <- prcomp(tib)
eigs <- pca$sdev^2
pc1 <- eigs[1] / sum(eigs)
pc2 <- eigs[2] / sum(eigs)

pca_tib <- as.tibble(pca$x)
pca_tib$names <- names[3:length(names)]
pca_tib$group <- groups[3:length(groups)]


#plot all
p <- ggplot(data = pca_tib,aes(x=PC1,y=PC2,color=group)) + geom_point(size = 4) +
  xlab(paste('PC1', (pc1 * 100),'%')) + 
  ylab(paste('PC2', (pc2 * 100),'%'))
p
ggsave('MixedSamples/pca-plots/pca-all.png',width = 15, height = 10, units = c("in"))

#plot just blood semen
temp = pca_tib[pca_tib$group %in% c('Blood','SE','SA',"B_SE_1.1","B_SE_1.100","B_SE_100.1"),]
p <- ggplot(data = temp,aes(x=PC1,y=PC2,color=group)) + geom_point(size = 4) +
  xlab(paste('PC1', (pc1 * 100),'%')) + 
  ylab(paste('PC2', (pc2 * 100),'%'))
p
ggsave('MixedSamples/pca-plots/pca-blood-semen.png',width = 15, height = 10, units = c("in"))

#plot just blood saliva
temp = pca_tib[pca_tib$group %in% c('Blood','SE','SA',"B_SA_1.1","B_SA_1.100","B_SA_100.1"),]
p <- ggplot(data = temp,aes(x=PC1,y=PC2,color=group)) + geom_point(size = 4) +
  xlab(paste('PC1', (pc1 * 100),'%')) + 
  ylab(paste('PC2', (pc2 * 100),'%'))
p
ggsave('MixedSamples/pca-plots/pca-blood-saliva.png',width = 15, height = 10, units = c("in"))


#plot just semen saliva
temp = pca_tib[pca_tib$group %in% c('Blood','SE','SA',"SE_SA_1.1","SE_SA_1.100","SE_SA_100.1"),]
p <- ggplot(data = temp,aes(x=PC1,y=PC2,color=group)) + geom_point(size = 4) +
  xlab(paste('PC1', (pc1 * 100),'%')) + 
  ylab(paste('PC2', (pc2 * 100),'%'))
p
ggsave('MixedSamples/pca-plots/pca-semen-saliva.png',width = 15, height = 10, units = c("in"))


#create the average of each group
data$blood_avg = (data$BloodYoung1 + data$Blood_WM24 + data$Blood_WM26) / 3
std <- apply(data,1,function(row){
  sd(c(row[3],row[13],row[14]))
})
data$blood_std <- std

data$saliva_avg = (data$SA7 + data$SA8 + data$SA9) / 3
std <- apply(data,1,function(row){
  sd(c(row[4],row[5],row[6]))
})
data$saliva_std <- std

data$semen_avg = (data$SE1005 + data$SE1022 + data$SE1027) / 3
std <- apply(data,1,function(row){
  sd(c(row[18],row[19],row[20]))
})
data$semen_std <- std

#create fake 1 to 1
#SE to SA
data$in_silico_SE_SA_1_1 <- (data$semen_avg * .5) + (data$saliva_avg * .5)
#SE Blood
data$in_silico_SE_blood_1_1 <- (data$semen_avg * .5) + (data$blood_avg * .5)
#Saliba Blood
data$in_silico_SA_blood_1_1 <- (data$saliva_avg * .5) + (data$blood_avg * .5)

write_csv(data,'MixedSamples/avrg-std-data.csv')

names <- c(names,c('blood_average','blood_std','saliva_average','saliva_std','semen_average','semen_std','in_silico_SE_SA_1_1','in_silico_SE_blood_1_1','in_silico_SA_blood_1_1'))
old_groups <- c(groups,c('blood_average','blood_std','saliva_average','saliva_std','semen_average','semen_std','in_silico_SE_SA_1_1','in_silico_SE_blood_1_1','in_silico_SA_blood_1_1'))
groups <- c(groups,c('blood_average','blood_std','saliva_average','saliva_std','semen_average','semen_std','in_silico_SE_SA_1_1','in_silico_SE_blood_1_1','in_silico_SA_blood_1_1'))

c_groups <- groups
colors <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928','#f0017f')

group.colors <- c(A = "#333BFF", B = "#CC6600", C ="#9633FF", D = "#E2FF33", E = "#E3DB71")

color_groups = c(Blood=colors[6], SA=colors[2], B_SA_1.1=colors[8], B_SA_1.100=colors[10], B_SA_100.1=colors[12], B_SE_1.1=colors[8], B_SE_1.100=colors[10], B_SE_100.1=colors[12], SE_SA_1.1=colors[8], SE_SA_1.100=colors[10], SE_SA_100.1=colors[12] ,SE=colors[4], blood_average=colors[6],blood_std='#000000',saliva_average=colors[2],saliva_std='#000000',semen_average=colors[4],semen_std='#000000',in_silico_SE_SA_1_1=colors[13], in_silico_SE_blood_1_1=colors[13], in_silico_SA_blood_1_1=colors[13])

m <- as.matrix(data[,3:ncol(data)])
colnames(m) <- NULL

tm <- t(m)

tib <- as.tibble(tm)

pca <- prcomp(tib)
eigs <- pca$sdev^2
pc1 <- eigs[1] / sum(eigs)
pc2 <- eigs[2] / sum(eigs)

pca_tib <- as.tibble(pca$x)
pca_tib$names <- names[3:length(names)]
pca_tib$group <- groups[3:length(groups)]
pca_tib$c_group <- c_groups[3:length(groups)]
#plot just blood semen
temp = pca_tib[pca_tib$group %in% c('Blood','SE','SA',"B_SE_1.1","B_SE_1.100","B_SE_100.1","in_silico_SE_blood_1_1"),]
p <- ggplot(data = temp,aes(x=PC1,y=PC2,color=group)) + geom_point(size = 4) +
  xlab(paste('PC1', (pc1 * 100),'%')) + 
  ylab(paste('PC2', (pc2 * 100),'%'))
p
ggsave('MixedSamples/pca-plots/pca-blood-semen-in-silico.png',width = 15, height = 10, units = c("in"))

#plot just blood saliva
temp = pca_tib[pca_tib$group %in% c('Blood','SE','SA',"B_SA_1.1","B_SA_1.100","B_SA_100.1","in_silico_SA_blood_1_1"),]
p <- ggplot(data = temp,aes(x=PC1,y=PC2,color=group)) + geom_point(size = 4) +
  xlab(paste('PC1', (pc1 * 100),'%')) + 
  ylab(paste('PC2', (pc2 * 100),'%'))
p
ggsave('MixedSamples/pca-plots/pca-blood-saliva-in-silico.png',width = 15, height = 10, units = c("in"))


#plot just semen saliva
temp = pca_tib[pca_tib$group %in% c('Blood','SE','SA',"SE_SA_1.1","SE_SA_1.100","SE_SA_100.1","in_silico_SE_SA_1_1"),]
p <- ggplot(data = temp,aes(x=PC1,y=PC2,color=group)) + geom_point(size = 4) +
  xlab(paste('PC1', (pc1 * 100),'%')) + 
  ylab(paste('PC2', (pc2 * 100),'%')) + 
  scale_fill_manual(values=color_groups)
p
ggsave('MixedSamples/pca-plots/pca-semen-saliva-in-silico.png',width = 15, height = 10, units = c("in"))



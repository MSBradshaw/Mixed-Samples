library(readr)

data <- read_csv('MixedSamples/avrg-std-data.csv')


#SE SA 1 to 1 plus or minus part of std
for( i in c(1:10)){
  name <- paste('in_silico_SE_SA_1_1..',i,sep='')
  data[name] <- apply(data,1,function(row){
    stdse <- as.numeric(row[26])
    avgse <- as.numeric(row[25])
    stdsa <- as.numeric(row[24])
    avgsa <- as.numeric(row[23])
    rse <- runif(1,(-1*stdse),stdse)
    rsa <- runif(1,(-1*stdsa),stdsa)

    return(((avgse * .5)  + rse) + ((avgsa * .5)  + rsa))
  })
}

#SE SA 1 to 1 plus or minus part of std
for( i in c(1:10)){
  name <- paste('in_silico_SE_SA_3_1..',i,sep='')
  data[name] <- apply(data,1,function(row){
    stdse <- as.numeric(row[26])
    avgse <- as.numeric(row[25])
    stdsa <- as.numeric(row[24])
    avgsa <- as.numeric(row[23])
    rse <- runif(1,(-1*stdse),stdse)
    rsa <- runif(1,(-1*stdsa),stdsa)

    return(((avgse * .75)  + rse) + ((avgsa * .25)  + rsa))
  })
}


#SE SA 1 to 1 plus or minus part of std
for( i in c(1:10)){
  name <- paste('in_silico_SE_SA_1_3..',i,sep='')
  data[name] <- apply(data,1,function(row){
    stdse <- as.numeric(row[26])
    avgse <- as.numeric(row[25])
    stdsa <- as.numeric(row[24])
    avgsa <- as.numeric(row[23])
    rse <- runif(1,(-1*stdse),stdse)
    rsa <- runif(1,(-1*stdsa),stdsa)
    
    return(((avgse * .25)  + rse) + ((avgsa * .75)  + rsa))
  })
}


#SE Blood 1 to 1 plus or minus part of std
for( i in c(1:10)){
  name <- paste('in_silico_SE_Blood_1_1..',i,sep='')
  data[name] <- apply(data,1,function(row){
    stdse <- as.numeric(row[26])
    avgse <- as.numeric(row[25])
    stdblood <- as.numeric(row[22])
    avgblood <- as.numeric(row[21])
    rse <- runif(1,(-1*stdse),stdse)
    rsa <- runif(1,(-1*stdblood),stdblood)
    
    return(((avgse * .5)  + rse) + ((stdblood * .5)  + rsa))
  })
}

for( i in c(1:10)){
  name <- paste('in_silico_SE_Blood_3_1..',i,sep='')
  data[name] <- apply(data,1,function(row){
    stdse <- as.numeric(row[26])
    avgse <- as.numeric(row[25])
    stdblood <- as.numeric(row[22])
    avgblood <- as.numeric(row[21])
    rse <- runif(1,(-1*stdse),stdse)
    rsa <- runif(1,(-1*stdblood),stdblood)
    
    return(((avgse * .75)  + rse) + ((stdblood * .25)  + rsa))
  })
}


for( i in c(1:10)){
  name <- paste('in_silico_SE_Blood_1_3..',i,sep='')
  data[name] <- apply(data,1,function(row){
    stdse <- as.numeric(row[26])
    avgse <- as.numeric(row[25])
    stdblood <- as.numeric(row[22])
    avgblood <- as.numeric(row[21])
    rse <- runif(1,(-1*stdse),stdse)
    rsa <- runif(1,(-1*stdblood),stdblood)
    
    return(((avgse * .25)  + rse) + ((stdblood * .75)  + rsa))
  })
}

#SA Blood 1 to 1 plus or minus part of std
for( i in c(1:10)){
  name <- paste('in_silico_SA_Blood_1_1..',i,sep='')
  data[name] <- apply(data,1,function(row){
    stdsa <- as.numeric(row[24])
    avgsa <- as.numeric(row[23])
    stdblood <- as.numeric(row[22])
    avgblood <- as.numeric(row[21])
    rse <- runif(1,(-1*stdsa),stdsa)
    rsa <- runif(1,(-1*stdblood),stdblood)
    
    return(((avgsa * .5)  + rse) + ((stdblood * .5)  + rsa))
  })
}

for( i in c(1:10)){
  name <- paste('in_silico_SA_Blood_3_1..',i,sep='')
  data[name] <- apply(data,1,function(row){
    stdsa <- as.numeric(row[24])
    avgsa <- as.numeric(row[23])
    stdblood <- as.numeric(row[22])
    avgblood <- as.numeric(row[21])
    rse <- runif(1,(-1*stdsa),stdsa)
    rsa <- runif(1,(-1*stdblood),stdblood)
    
    return(((avgsa * .75)  + rse) + ((stdblood * .25)  + rsa))
  })
}

for( i in c(1:10)){
  name <- paste('in_silico_SA_Blood_1_3..',i,sep='')
  data[name] <- apply(data,1,function(row){
    stdsa <- as.numeric(row[24])
    avgsa <- as.numeric(row[23])
    stdblood <- as.numeric(row[22])
    avgblood <- as.numeric(row[21])
    rse <- runif(1,(-1*stdsa),stdsa)
    rsa <- runif(1,(-1*stdblood),stdblood)
    
    return(((avgsa * .25)  + rse) + ((stdblood * .75)  + rsa))
  })
}

names <- colnames(data)
groups <- colnames(data)
#make all insilico names their group name
groups <- sub('\\.\\..*','',groups)
groups[c(4,5,6)] <- 'Saliva'
groups[c(18,19,20)] <- 'Semen'
groups[c(3,13,14)] <- 'Blood'


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

#plot just semen saliva
temp = pca_tib[pca_tib$group %in% c('Blood','Semen','Saliva',"SE_SA_1.1","SE_SA_1.100","SE_SA_100.1","in_silico_SE_SA_1_1",'in_silico_SE_SA_3_1','in_silico_SE_SA_1_3'),]
p <- ggplot(data = temp,aes(x=PC1,y=PC2,color=group)) + geom_point(size = 4) +
  xlab(paste('PC1', (pc1 * 100),'%')) + 
  ylab(paste('PC2', (pc2 * 100),'%')) + 
  scale_fill_manual(values=color_groups)
p
ggsave('MixedSamples/pca-plots/pca-semen-saliva-in-silico--std.png',width = 15, height = 10, units = c("in"))

#plot just semen blood
temp = pca_tib[pca_tib$group %in% c('Blood','Semen','Saliva',"B_SE_1.1","B_SE_100.1","B_SE_1.100","in_silico_SE_Blood_1_1",'in_silico_SE_Blood_3_1','in_silico_SE_Blood_1_3'),]
p <- ggplot(data = temp,aes(x=PC1,y=PC2,color=group)) + geom_point(size = 4) +
  xlab(paste('PC1', (pc1 * 100),'%')) + 
  ylab(paste('PC2', (pc2 * 100),'%')) + 
  scale_fill_manual(values=color_groups)
p
ggsave('MixedSamples/pca-plots/pca-semen-blood-in-silico--std.png',width = 15, height = 10, units = c("in"))


#plot just saliva blood
temp = pca_tib[pca_tib$group %in% c('Blood','Semen','Saliva',"B_SA_1.1","B_SA_100.1","B_SA_1.100","in_silico_SA_Blood_1_1",'in_silico_SA_Blood_3_1','in_silico_SA_Blood_1_3'),]
p <- ggplot(data = temp,aes(x=PC1,y=PC2,color=group)) + geom_point(size = 4) +
  xlab(paste('PC1', (pc1 * 100),'%')) + 
  ylab(paste('PC2', (pc2 * 100),'%')) + 
  scale_fill_manual(values=color_groups)
p
ggsave('MixedSamples/pca-plots/pca-blood-saliva-in-silico--std.png',width = 15, height = 10, units = c("in"))

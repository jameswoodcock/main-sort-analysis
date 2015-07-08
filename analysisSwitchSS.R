rm(list = ls())
library(DistatisR)
library(MASS)
library(vegan)
library(ggplot2)
library(pvclust)
library(dendextend)
library(gplots)
library(devtools)

material = "SS"
screePlot = FALSE
ellipsePlot = FALSE
figPath = paste("./plots/",material,"/",sep="")
pvals = FALSE

data <- read.csv(paste("./data/data",material,".csv",sep=""),row.names=1,header=T)
experience <- read.csv(paste("./data/experience.csv",sep=""),row.names=1,header=T)

experience <- experience[colnames(data),]

acousticsExp <- array(rep(experience[,1],each = nrow(data)^2), c(nrow(data),nrow(data),ncol(data)))
engineerExp <- array(rep(experience[,2],each = nrow(data)^2), c(nrow(data),nrow(data),ncol(data)))
musicExp <- array(rep(experience[,3],each = nrow(data)^2), c(nrow(data),nrow(data),ncol(data)))

Nsubs = dim(data)[2]
NsubsAcoustics = sum(experience[,1])
NsubsEngineer = sum(experience[,2])
NsubsMusic = sum(experience[,3])

distCube <- DistanceFromSort(data)
distCubeAcoustics <- distCube*acousticsExp
distCubeEngineer <- distCube*engineerExp
distCubeMusic <- distCube*musicExp
distCubeNonAcoustics <- distCube*abs(acousticsExp-1)
distCubeNonEngineer <- distCube*abs(engineerExp-1)
distCubeNonMusic <- distCube*abs(musicExp-1)
totDist <- apply(distCube,c(1,2),sum)/Nsubs
totDistAcoustics <- apply(distCubeAcoustics,c(1,2),sum)/NsubsAcoustics
totDistEngineer <- apply(distCubeEngineer,c(1,2),sum)/NsubsEngineer
totDistMusic <- apply(distCubeMusic,c(1,2),sum)/NsubsMusic
totDistNonAcoustics <- apply(distCubeNonAcoustics,c(1,2),sum)/(Nsubs-NsubsAcoustics)
totDistNonEngineer <- apply(distCubeNonEngineer,c(1,2),sum)/(Nsubs-NsubsEngineer)
totDistNonMusic <- apply(distCubeNonMusic,c(1,2),sum)/(Nsubs-NsubsMusic)

clust = hclust(dist(totDist),"ward.D2")
clustAcoustics = hclust(dist(totDistAcoustics),"ward.D2")
clustEngineer = hclust(dist(totDistEngineer),"ward.D2")
clustMusic = hclust(dist(totDistMusic),"ward.D2")
clustNonAcoustics = hclust(dist(totDistNonAcoustics),"ward.D2")
clustNonEngineer = hclust(dist(totDistNonEngineer),"ward.D2")
clustNonMusic = hclust(dist(totDistNonMusic),"ward.D2")
if (pvals == TRUE){
pvalsClust <- pvclust(totDist,method.hclust = "ward.D2",method.dist="euclidean",nboot=1000)
}

ngroups <- apply(data,2,max)
ngroupsEngineer <- ngroups*experience[,1]
ngroupsEngineer[which(ngroupsEngineer==0)] = NA
ngroupsNonEngineer <- ngroups*abs(experience[,1]-1)
ngroupsNonEngineer[which(ngroupsNonEngineer==0)] = NA

ngroupsbar <- data.frame(participant = colnames(data),ngroups = ngroups)
mdsClusters = ceiling(median(ngroups))
#mdsClusters = 8
mdsClustersEngineer = ceiling(median(ngroupsEngineer,na.rm=TRUE))
mdsClustersNonEngineer = ceiling(median(ngroupsNonEngineer,na.rm=TRUE))



pdf(paste(figPath,"ngroups.pdf",sep=""),width = 10,height = 10)
print(
ggplot(ngroupsbar,aes(x=participant,ngroups)) + geom_bar(stat="identity") + labs(y="Number of groups",x="Participant") +theme_bw()
)
dev.off()

pdf(paste(figPath,"dendro_all.pdf",sep=""),width = 11, height = 8)
#plot(clust,hang=-1)
plot(clust, xlab=NA, sub=NA, main=NA, cex = 0.5)
rect.hclust(clust,k=mdsClusters,border="red")
dev.off()
pdf(paste(figPath,"dendro_engineer.pdf",sep=""),width = 11, height = 8)
#plot(clust,hang=-1)
plot(clustEngineer, xlab=NA, sub=NA, main=NA, cex = 0.5)
rect.hclust(clustEngineer,k=mdsClustersEngineer,border="red")
dev.off()
pdf(paste(figPath,"dendro_non_engineer.pdf",sep=""),width = 11, height = 8)
#plot(clust,hang=-1)
plot(clustNonEngineer, xlab=NA, sub=NA, main=NA, cex = 0.5)
rect.hclust(clustNonEngineer,k=mdsClustersNonEngineer,border="red")
dev.off()


for (i in 2:11) {
#pdf(paste(figPath,"dendro_",i,"groups.pdf",sep=""),width = 11, height = 8)
#plot(clust, xlab=NA, sub=NA, main=NA, cex = 0.5)
#rect.hclust(clust,k=i,border="red")
#dev.off()
if (i==mdsClusters){
pdf(paste(figPath,"dendro_median_groups.pdf",sep=""),width = 11, height = 8)
#plot(clust,hang=-1)
plot(clust, xlab=NA, sub=NA, main=NA, cex = 0.5)
rect.hclust(clust,k=i,border="red")
dev.off()
}
}

#Do some MDS!!!!
mdsRes <- mmds(totDist)

if (screePlot == TRUE){
nmmds <- metaMDS(totDist,k=2,zerodist="add")
stress <- nmmds$stress
nmmds <- metaMDS(totDist,k=3,zerodist="add")
stress <- append(stress,nmmds$stress)
nmmds <- metaMDS(totDist,k=4,zerodist="add")
stress <- append(stress,nmmds$stress)
nmmds <- metaMDS(totDist,k=5,zerodist="add")
stress <- append(stress,nmmds$stress)
nmmds <- metaMDS(totDist,k=6,zerodist="add")
stress <- append(stress,nmmds$stress)
nmmds <- metaMDS(totDist,k=7,zerodist="add")
stress <- append(stress,nmmds$stress)
nmmds <- metaMDS(totDist,k=8,zerodist="add")
stress <- append(stress,nmmds$stress)
nmmds <- metaMDS(totDist,k=9,zerodist="add")
stress <- append(stress,nmmds$stress)
dims <- 2:9

scree <- data.frame(dims,stress)

pdf(paste(figPath,"scree.pdf",sep=""),width = 10,height = 10)
print(
ggplot(scree, aes(x = dims, y = stress)) + geom_line(linetype=2) + geom_point(shape=1) + theme_bw() + labs(x = "Number of dimensions",y = "Stress")
)
dev.off()
}

mdsDims = 3
nmmds <- metaMDS(totDist,k=mdsDims,zerodist="ignore")
nmmdsEngineer <- metaMDS(totDistEngineer,k=mdsDims,zerodist="ignore")
nmmdsNonEngineer <- metaMDS(totDistNonEngineer,k=mdsDims,zerodist="ignore")

#dim1nm <- mdsRes$FactorScore[,1]
#dim2nm <- mdsRes$FactorScore[,2]
#dim3nm <- mdsRes$FactorScore[,3]

groups = cutree(clust,k=mdsClusters)
groupsEngineer = cutree(clustEngineer,k=mdsClustersEngineer)
groupsNonEngineer = cutree(clustNonEngineer,k=mdsClustersNonEngineer)

df3 <- data.frame(nmmds$points,groups)
df3Engineer <- data.frame(nmmdsEngineer$points,groupsEngineer)
df3NonEngineer <- data.frame(nmmdsNonEngineer$points,groupsNonEngineer)

groupLabs <- "1"
for (i in 2:10){
groupLabs <- append(groupLabs,toString(i))
}

textSize = 3

pdf(paste(figPath,"mdsdim1dim2text.pdf",sep=""),width = 10, height = 10)
print(
ggplot(df3,aes(x = MDS1,y = MDS2,label=rownames(df3),color=factor(groups))) + geom_text(size = textSize) + theme_bw() + xlim(-0.6,0.6) + ylim(-0.6,0.6) + labs(x = "Dimension I",y = "Dimension II") + scale_color_discrete(name = "Group",labels = groupLabs) + theme(legend.position = c(.9,.9)) + ggtitle(material)
)
dev.off()

pdf(paste(figPath,"mdsdim1dim3text.pdf",sep=""),width = 10, height = 10)
print(
ggplot(df3,aes(x = MDS1,y = MDS3,label=rownames(df3),color=factor(groups))) + geom_text(size = textSize) + theme_bw() + xlim(-0.6,0.6) + ylim(-0.6,0.6) + labs(x = "Dimension I",y = "Dimension III") + scale_color_discrete(name = "Group",labels = groupLabs) + theme(legend.position = c(.9,.9)) + ggtitle(material)
)
dev.off()

pdf(paste(figPath,"mdsdim2dim3text.pdf",sep=""),width = 10, height = 10)
print(
ggplot(df3,aes(x = MDS2,y = MDS3,label=rownames(df3),color=factor(groups))) + geom_text(size = textSize) + theme_bw() + xlim(-0.6,0.6) + ylim(-0.6,0.6) + labs(x = "Dimension II",y = "Dimension III") + scale_color_discrete(name = "Group",labels = groupLabs) + theme(legend.position = c(.9,.9)) + ggtitle(material)
)
dev.off()

if (mdsDims > 3){
pdf(paste(figPath,"mdsdim1dim4text.pdf",sep=""),width = 10, height = 10)
print(
ggplot(df3,aes(x = MDS1,y = MDS4,label=rownames(df3),color=factor(groups))) + geom_text(size = textSize) + theme_bw() + xlim(-0.6,0.6) + ylim(-0.6,0.6) + labs(x = "Dimension I",y = "Dimension IV") + scale_color_discrete(name = "Group",labels = groupLabs) + theme(legend.position = c(.9,.9)) + ggtitle(material)
)
dev.off()

pdf(paste(figPath,"mdsdim3dim4text.pdf",sep=""),width = 10, height = 10)
print(
ggplot(df3,aes(x = MDS3,y = MDS4,label=rownames(df3),color=factor(groups))) + geom_text(size = textSize) + theme_bw() + xlim(-0.6,0.6) + ylim(-0.6,0.6) + labs(x = "Dimension I",y = "Dimension IV") + scale_color_discrete(name = "Group",labels = groupLabs) + theme(legend.position = c(.9,.9)) + ggtitle(material)
)
dev.off()
}

##DISTATIS analysis

#testDistatis <- distatis(distCube)
#BootF <- BootFactorScores(testDistatis$res4Splus$PartialF,niter=1000)
#LeF	<- testDistatis$res4Splus$F
#PartialFS <- testDistatis$res4Splus$PartialF
#GraphDistatisBoot(LeF,BootF,PartialFS,ZeTitle="Bootstrap on Factors")

#df4 <- data.frame(testDistatis$res4Cmat$G,experience)

#ggplot(df4,aes(x = dim.1,y = dim.2,label=rownames(df4),color=factor(engineering))) + geom_point(aes(shape=factor(engineering)),size=5) + theme_bw()

##Procrustes analysis

pro <- procrustes(nmmdsEngineer,nmmdsNonEngineer,scale=TRUE)
pdf(paste(figPath,"procrustes.pdf",sep=""),width = 10, height = 10)
plot(pro)
dev.off()

proNonEngineer <- procrustes(nmmds,nmmdsNonEngineer,scale=TRUE)
df3NonEngineer["MDS1"] <- proNonEngineer$Yrot[,1]
df3NonEngineer["MDS2"] <- proNonEngineer$Yrot[,2]
df3NonEngineer["MDS3"] <- proNonEngineer$Yrot[,3]
if (mdsDims >3){
df3NonEngineer["MDS4"] <- proNonEngineer$Yrot[,4]
}

proEngineer <- procrustes(nmmds,nmmdsEngineer,scale=TRUE)
df3Engineer["MDS1"] <- proEngineer$Yrot[,1]
df3Engineer["MDS2"] <- proEngineer$Yrot[,2]
df3Engineer["MDS3"] <- proEngineer$Yrot[,3]
if (mdsDims >3){
df3Engineer["MDS4"] <- proEngineer$Yrot[,4]
}

pdf(paste(figPath,"mdsdim1dim2textengineer.pdf",sep=""),width = 10, height = 10)
#plot(dim1nm,dim2nm,type="n")
#text(dim1nm, dim2nm, labels = row.names(data), cex=.7)
print(
ggplot(df3Engineer,aes(x = MDS1,y = MDS2,label=rownames(df3),color=factor(groupsEngineer))) + geom_text(size = textSize) + theme_bw() + xlim(-0.6,0.6) + ylim(-0.6,0.6) + labs(x = "Dimension I",y = "Dimension II") + scale_color_discrete(name = "Group",labels = groupLabs) + theme(legend.position = c(.9,.9)) + ggtitle(material)
)
#plot(nmmds,type="t")
dev.off()

pdf(paste(figPath,"mdsdim1dim2textnonengineer.pdf",sep=""),width = 10, height = 10)
#plot(dim1nm,dim2nm,type="n")
#text(dim1nm, dim2nm, labels = row.names(data), cex=.7)
print(
ggplot(df3NonEngineer,aes(x = MDS1,y = MDS2,label=rownames(df3),color=factor(groupsNonEngineer))) + geom_text(size = textSize) + theme_bw() + xlim(-0.6,0.6) + ylim(-0.6,0.6) + labs(x = "Dimension I",y = "Dimension II") + scale_color_discrete(name = "Group",labels = groupLabs) + theme(legend.position = c(.9,.9)) + ggtitle(material)
)
#plot(nmmds,type="t")
dev.off()

#####CUT#######

#pdf(paste(figPath,"mdsdim1dim3.pdf",sep=""),width = 10,height = 10)
#plot(dim1nm,dim2nm,type="n")
#text(dim1nm, dim2nm, labels = row.names(data), cex=.7)
#print(
#ggplot(df3,aes(x = MDS1,y = MDS3,label=rownames(df3),color=factor(groups))) + geom_point(aes(shape=factor(groups)),size=5) + theme_bw() + xlim(-0.6,0.6) + ylim(-0.6,0.6) + labs(x = "Dimension I",y = "Dimension III")+ scale_color_discrete(name = "Group",labels = groupLabs) + scale_shape_discrete(name = "Group",labels = groupLabs) + theme(legend.position = c(.9,.9)) + ggtitle(material)
#)
#plot(nmmds,type="t")
#dev.off()

#pdf(paste(figPath,"mdsdim1dim2.pdf",sep=""),width = 10,height = 10)
#plot(dim1nm,dim2nm,type="n")
#text(dim1nm, dim2nm, labels = row.names(data), cex=.7)
#print(
#ggplot(df3,aes(x = MDS1,y = MDS2,label=rownames(df3),color=factor(groups))) + geom_point(aes(shape=factor(groups)),size=5) + theme_bw() + xlim(-0.6,0.6) + ylim(-0.6,0.6) + labs(x = "Dimension I",y = "Dimension II")+ scale_color_discrete(name = "Group",labels = groupLabs) + scale_shape_discrete(name = "Group",labels = groupLabs) + theme(legend.position = c(.9,.9)) + ggtitle(material)
#)
#plot(nmmds,type="t")
#dev.off()


#dim1nm <- nmmds$points[,1]
#dim2nm <- nmmds$points[,2]
#dim3nm <- nmmds$points[,3]
#dim4nm <- nmmds$points[,4]
#dim1nmEngineer <- nmmdsEngineer$points[,1]
#dim2nmEngineer <- nmmdsEngineer$points[,2]
#dim3nmEngineer <- nmmdsEngineer$points[,3]
#dim4nmEngineer <- nmmdsEngineer$points[,4]
#dim1nmNonEngineer <- nmmdsNonEngineer$points[,1]
#dim2nmNonEngineer <- nmmdsNonEngineer$points[,2]
#dim3nmNonEngineer <- nmmdsNonEngineer$points[,3]
#dim4nmNonEngineer <- nmmdsNonEngineer$points[,4]

#Ellipses

if (ellipsePlot == TRUE){
library(devtools)
source_url("https://raw.github.com/JoFrhwld/FAAV/master/r/stat-ellipse.R")
#ggplot(df3Engineer,aes(x = MDS1,y = MDS2,label=rownames(df3),color=factor(groups))) + stat_ellipse(geom = "polygon", alpha = 1/2, aes(fill = factor(groups))) + geom_point(aes(shape=factor(groups)),size=1)

df3centroidsEngineer <- aggregate(df3Engineer,by=list(groupsEngineer),FUN="mean")
pdf(paste(figPath,"mdsdim1dim2ellipseengineer.pdf",sep=""),width = 10, height = 10)
#plot(dim1nm,dim2nm,type="n")
#text(dim1nm, dim2nm, labels = row.names(data), cex=.7)
print(
ggplot(df3Engineer,aes(x = MDS1,y = MDS2,label=rownames(df3Engineer),color=factor(groupsEngineer))) + stat_ellipse(geom = "polygon", alpha = 1/2, aes(fill = factor(groupsEngineer))) + geom_point(size=1) + theme_bw()  + labs(x = "Dimension I",y = "Dimension II") + scale_color_discrete(name = "Group",labels = groupLabs)  + ggtitle(material) + geom_text(data=df3centroidsEngineer,aes(x = MDS1,y = MDS2,label=c("1","2","3","4","5"),color=factor(df3centroidsEngineer$groups)),fontface="bold") + theme(legend.position = "none")
)# 
#)
#plot(nmmds,type="t") + xlim(-2,2) + ylim(-2,2)
dev.off()

df3centroidsNonEngineer <- aggregate(df3NonEngineer,by=list(groupsNonEngineer),FUN="mean")
pdf(paste(figPath,"mdsdim1dim2ellipsenonengineer.pdf",sep=""),width = 10, height = 10)
#plot(dim1nm,dim2nm,type="n")
#text(dim1nm, dim2nm, labels = row.names(data), cex=.7)
print(
ggplot(df3NonEngineer,aes(x = MDS1,y = MDS2,label=rownames(df3NonEngineer),color=factor(groupsNonEngineer))) + stat_ellipse(geom = "polygon", alpha = 1/2, aes(fill = factor(groupsNonEngineer))) + geom_point(size=1) + theme_bw() + xlim(-0.6,0.6) + ylim(-0.6,0.6) + labs(x = "Dimension I",y = "Dimension II") + scale_color_discrete(name = "Group",labels = groupLabs)  + ggtitle(material)+ geom_text(data=df3centroidsNonEngineer,aes(x = MDS1,y = MDS2,label=c("Individually\n identifiable\n background","Sounds related\n to action\n and movement","Dialogue","NO CLEAR\n GROUP","Non-diegetic music\n and diffuse background"),color=factor(groupsNonEngineer))) + theme(legend.position = "none")
)# 
#)
#plot(nmmds,type="t")
dev.off()

df3centroids <- aggregate(df3,by=list(groups),FUN="mean")
pdf(paste(figPath,"mdsdim1dim2ellipseall.pdf",sep=""),width = 10, height = 10)
#plot(dim1nm,dim2nm,type="n")
#text(dim1nm, dim2nm, labels = row.names(data), cex=.7)
print(
ggplot(df3,aes(x = MDS1,y = MDS2,label=rownames(df3),color=factor(groups))) + stat_ellipse(geom = "polygon", alpha = 1/2, aes(fill = factor(groups))) + geom_point(size=1) + theme_bw() + xlim(-0.8,0.8) + ylim(-0.8,0.8) + labs(x = "Dimension I",y = "Dimension II") + scale_color_discrete(name = "Group",labels = groupLabs)  + ggtitle(material)+ geom_text(data=df3centroids,aes(x = MDS1,y = MDS2,label=c("Music in\n scene","Distant\nbackground","Vehicle\nsounds","Prominent\nevent\nsounds","Diffuse\ncontinuous\nbackground","Low level\nevent sounds","Clear\nspeech","Unitelligible\nvoices"),color=factor(groups))) + theme(legend.position = "none")
)# 
#)
#plot(nmmds,type="t")
dev.off()
pdf(paste(figPath,"mdsdim1dim3ellipseall.pdf",sep=""),width = 10, height = 10)
#plot(dim1nm,dim2nm,type="n")
#text(dim1nm, dim2nm, labels = row.names(data), cex=.7)
print(
ggplot(df3,aes(x = MDS1,y = MDS3,label=rownames(df3),color=factor(groups))) + stat_ellipse(geom = "polygon", alpha = 1/2, aes(fill = factor(groups))) + geom_point(size=1) + theme_bw() + xlim(-0.8,0.8) + ylim(-0.8,0.8) + labs(x = "Dimension I",y = "Dimension III") + scale_color_discrete(name = "Group",labels = groupLabs)  + ggtitle(material)+ geom_text(data=df3centroids,aes(x = MDS1,y = MDS3,label=c("Music in\n scene","Distant\nbackground","Vehicle\nsounds","Prominent\nevent\nsounds","Diffuse\ncontinuous\nbackground","Low level\nevent sounds","Clear\nspeech","Unitelligible\nvoices"),color=factor(groups))) + theme(legend.position = "none")
)# 
#)
#plot(nmmds,type="t")
dev.off()
pdf(paste(figPath,"mdsdim2dim3ellipseall.pdf",sep=""),width = 10, height = 10)
#plot(dim1nm,dim2nm,type="n")
#text(dim1nm, dim2nm, labels = row.names(data), cex=.7)
print(
ggplot(df3,aes(x = MDS2,y = MDS3,label=rownames(df3),color=factor(groups))) + stat_ellipse(geom = "polygon", alpha = 1/2, aes(fill = factor(groups))) + geom_point(size=1) + theme_bw() + xlim(-0.8,0.8) + ylim(-0.8,0.8) + labs(x = "Dimension II",y = "Dimension III") + scale_color_discrete(name = "Group",labels = groupLabs)  + ggtitle(material)+ geom_text(data=df3centroids,aes(x = MDS2,y = MDS3,label=c("Music in\n scene","Distant\nbackground","Vehicle\nsounds","Prominent\nevent\nsounds","Diffuse\ncontinuous\nbackground","Low level\nevent sounds","Clear\nspeech","Unitelligible\nvoices"),color=factor(groups))) + theme(legend.position = "none")
)# 
#)
#plot(nmmds,type="t")
dev.off()
}

#Analysis of labels

dataExp <- read.csv(paste("./data/labels_data/labels",material,"Exp.csv",sep=""),row.names=1,header=T)
dataNonExp <- read.csv(paste("./data/labels_data/labels",material,"NonExp.csv",sep=""),row.names=1,header=T)
dataAll <- read.csv(paste("./data/labels_data/labels",material,"All.csv",sep=""),row.names=1,header=T)

pdf(paste(figPath,"clustergramExp.pdf",sep=""),width = 16, height = 18)
heatmap(as.matrix(dataExp),hclustfun=function(d) hclust(d,method='ward.D2'),margins=c(20,20),scale='none',col=grey(seq(0.9
,0,-0.01)),main=paste(material,"experienced listerners"))
dev.off()
pdf(paste(figPath,"clustergramNonExp.pdf",sep=""),width = 16, height = 18)
heatmap(as.matrix(dataNonExp),hclustfun=function(d) hclust(d,method='ward.D2'),margins=c(20,20),scale='none',col=grey(seq(0.9
,0,-0.01)),main=paste(material,"naive listerners"))
dev.off()
pdf(paste(figPath,"clustergramAll.pdf",sep=""),width = 16, height = 18)
heatmap(as.matrix(dataAll),hclustfun=function(d) hclust(d,method='ward.D2'),margins=c(20,20),scale='none',col=grey(seq(0.9
,0,-0.01)),main=paste(material,"all listerners"))
dev.off()

rowClust <- hclust(dist(dataAll),"ward.D2")
cutMat <- as.matrix(cutree(rowClust,k=mdsClusters))
clusterTable <- split(row.names(cutMat),cutMat)

pdf(paste(figPath,"dendro_labels_all.pdf",sep=""),width = 11, height = 8)
#plot(clust,hang=-1)
plot(rowClust, xlab=NA, sub=NA, main=NA, cex = 0.5)
rect.hclust(rowClust,k=mdsClusters,border="red")
dev.off()

rowClustNon <- hclust(dist(dataNonExp),"ward.D2")
cutMatNon <- as.matrix(cutree(rowClustNon,k=mdsClustersNonEngineer))
clusterTableNon <- split(row.names(cutMatNon),cutMatNon)

pdf(paste(figPath,"dendro_labels_non_exp.pdf",sep=""),width = 11, height = 8)
#plot(clust,hang=-1)
plot(rowClustNon, xlab=NA, sub=NA, main=NA, cex = 0.5)
rect.hclust(rowClustNon,k=mdsClustersNonEngineer,border="red")
dev.off()

rowClustExp <- hclust(dist(dataExp),"ward.D2")
cutMatExp <- as.matrix(cutree(rowClustExp,k=mdsClustersEngineer))
clusterTableExp <- split(row.names(cutMatExp),cutMatExp)

pdf(paste(figPath,"dendro_labels_exp.pdf",sep=""),width = 11, height = 8)
#plot(clust,hang=-1)
plot(rowClustExp, xlab=NA, sub=NA, main=NA, cex = 0.5)
rect.hclust(rowClustExp,k=mdsClustersEngineer,border="red")
#rect.dendrogram(as.dendrogram(rowClustExp),k=4,border = 8, lty = 5, lwd = 2)
dev.off()

for (i in 1:mdsClustersNonEngineer)
{
write.table(clusterTableNon[[paste(i)]],paste(figPath,"nonExpGroups",i,".csv",sep=""),sep=",")
}
for (i in 1:mdsClustersEngineer)
{
write.table(clusterTableExp[[paste(i)]],paste(figPath,"expGroups",i,".csv",sep=""),sep=",")
}
for (i in 1:mdsClusters)
{
write.table(clusterTable[[paste(i)]],paste(figPath,"allGroups",i,".csv",sep=""),sep=",")
}

colClust <- hclust(dist(t(dataAll)),"ward.D2")
pdf(paste(figPath,paste("labels_and_object_clusters_",material,".pdf",sep=""),sep=""),width = 11, height = 13)
par(mfrow=c(2,1))
plot(colClust, xlab=NA, sub=NA, main=NA, cex = 0.5)
rect.hclust(colClust,k=mdsClusters,border="red")
plot(rowClust, xlab=NA, sub=NA, main=NA, cex = 0.5)
rect.hclust(rowClust,k=mdsClusters,border="red")
dev.off()


##Prune dendrograms for paper figures

prunedClust <- as.dendrogram(rowClust)
for (i in 1:mdsClusters)
{
prunedClust <- prune(prunedClust,clusterTable[[i]][seq(length(clusterTable[[i]])-1)])
}

pdf(paste(figPath,"dendro_labels_pruned.pdf",sep=""),width = 11, height = 8)
par(mar=c(15,5,5,5))
plot(set(prunedClust,"labels",c("Low level event sounds","Continuous background","Music\nVehicle sounds\nLow level impact sounds","Human voices and vocalisations","High level event sounds")), xlab=NA,sub=NA, main=NA, cex = 0.5)
dev.off()






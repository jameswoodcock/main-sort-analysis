rm(list = ls())
library(DistatisR)
library(MASS)
library(vegan)
library(ggplot2)
library(pvclust)

material = "FF"
screePlot = FALSE
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
rect.hclust(clust,k=6,border="red")
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
pdf(paste(figPath,"dendro_",i,"groups.pdf",sep=""),width = 11, height = 8)
#plot(clust,hang=-1)
plot(clust, xlab=NA, sub=NA, main=NA, cex = 0.5)
rect.hclust(clust,k=i,border="red")
dev.off()
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

nmmds <- metaMDS(totDist,k=4,zerodist="ignore")

#dim1nm <- mdsRes$FactorScore[,1]
#dim2nm <- mdsRes$FactorScore[,2]
#dim3nm <- mdsRes$FactorScore[,3]

groups = cutree(clust,k=mdsClusters)
dim1nm <- nmmds$points[,1]
dim2nm <- nmmds$points[,2]
dim3nm <- nmmds$points[,3]
dim4nm <- nmmds$points[,4]
df3 <- data.frame(dim1 = nmmds$points[,1],dim2 = nmmds$points[,2],groups)

groupLabs <- "1"
for (i in 2:mdsClusters){
groupLabs <- append(groupLabs,toString(i))
}

textSize = 3

pdf(paste(figPath,"mdsdim1dim2.pdf",sep=""),width = 10,height = 10)
#plot(dim1nm,dim2nm,type="n")
#text(dim1nm, dim2nm, labels = row.names(data), cex=.7)
print(
ggplot(df3,aes(x = dim1nm,y = dim2nm,label=rownames(df3),color=factor(groups))) + geom_point(aes(shape=factor(groups)),size=5) + theme_bw() + xlim(-0.6,0.6) + ylim(-0.6,0.6) + labs(x = "Dimension I",y = "Dimension II")+ scale_color_discrete(name = "Group",labels = groupLabs) + scale_shape_discrete(name = "Group",labels = groupLabs) + theme(legend.position = c(.9,.9)) + ggtitle(material)
)
#plot(nmmds,type="t")
dev.off()

pdf(paste(figPath,"mdsdim1dim2text.pdf",sep=""),width = 10, height = 10)
#plot(dim1nm,dim2nm,type="n")
#text(dim1nm, dim2nm, labels = row.names(data), cex=.7)
print(
ggplot(df3,aes(x = dim1nm,y = dim2nm,label=rownames(df3),color=factor(groups))) + geom_text(size = textSize) + theme_bw() + xlim(-0.6,0.6) + ylim(-0.6,0.6) + labs(x = "Dimension I",y = "Dimension II") + scale_color_discrete(name = "Group",labels = groupLabs) + theme(legend.position = c(.9,.9)) + ggtitle(material)
)
#plot(nmmds,type="t")
dev.off()

pdf(paste(figPath,"mdsdim1dim3.pdf",sep=""),width = 10,height = 10)
#plot(dim1nm,dim2nm,type="n")
#text(dim1nm, dim2nm, labels = row.names(data), cex=.7)
print(
ggplot(df3,aes(x = dim1nm,y = dim3nm,label=rownames(df3),color=factor(groups))) + geom_point(aes(shape=factor(groups)),size=5) + theme_bw() + xlim(-0.6,0.6) + ylim(-0.6,0.6) + labs(x = "Dimension I",y = "Dimension III")+ scale_color_discrete(name = "Group",labels = groupLabs) + scale_shape_discrete(name = "Group",labels = groupLabs) + theme(legend.position = c(.9,.9)) + ggtitle(material)
)
#plot(nmmds,type="t")
dev.off()

pdf(paste(figPath,"mdsdim1dim3text.pdf",sep=""),width = 10, height = 10)
#plot(dim1nm,dim2nm,type="n")
#text(dim1nm, dim2nm, labels = row.names(data), cex=.7)
print(
ggplot(df3,aes(x = dim1nm,y = dim3nm,label=rownames(df3),color=factor(groups))) + geom_text(size = textSize) + theme_bw() + xlim(-0.6,0.6) + ylim(-0.6,0.6) + labs(x = "Dimension I",y = "Dimension III") + scale_color_discrete(name = "Group",labels = groupLabs) + theme(legend.position = c(.9,.9)) + ggtitle(material)
)
#plot(nmmds,type="t")
dev.off()

pdf(paste(figPath,"mdsdim1dim4text.pdf",sep=""),width = 10, height = 10)
#plot(dim1nm,dim2nm,type="n")
#text(dim1nm, dim2nm, labels = row.names(data), cex=.7)
print(
ggplot(df3,aes(x = dim1nm,y = dim4nm,label=rownames(df3),color=factor(groups))) + geom_text(size = textSize) + theme_bw() + xlim(-0.6,0.6) + ylim(-0.6,0.6) + labs(x = "Dimension I",y = "Dimension IV") + scale_color_discrete(name = "Group",labels = groupLabs) + theme(legend.position = c(.9,.9)) + ggtitle(material)
)
#plot(nmmds,type="t")
dev.off()

pdf(paste(figPath,"mdsdim2dim3text.pdf",sep=""),width = 10, height = 10)
#plot(dim1nm,dim2nm,type="n")
#text(dim1nm, dim2nm, labels = row.names(data), cex=.7)
print(
ggplot(df3,aes(x = dim2nm,y = dim3nm,label=rownames(df3),color=factor(groups))) + geom_text(size = textSize) + theme_bw() + xlim(-0.6,0.6) + ylim(-0.6,0.6) + labs(x = "Dimension II",y = "Dimension III") + scale_color_discrete(name = "Group",labels = groupLabs) + theme(legend.position = c(.9,.9)) + ggtitle(material)
)
#plot(nmmds,type="t")
dev.off()

pdf(paste(figPath,"mdsdim3dim4text.pdf",sep=""),width = 10, height = 10)
#plot(dim1nm,dim2nm,type="n")
#text(dim1nm, dim2nm, labels = row.names(data), cex=.7)
print(
ggplot(df3,aes(x = dim3nm,y = dim4nm,label=rownames(df3),color=factor(groups))) + geom_text(size = textSize) + theme_bw() + xlim(-0.6,0.6) + ylim(-0.6,0.6) + labs(x = "Dimension I",y = "Dimension IV") + scale_color_discrete(name = "Group",labels = groupLabs) + theme(legend.position = c(.9,.9)) + ggtitle(material)
)
#plot(nmmds,type="t")
dev.off()

##DISTATIS analysis

#testDistatis <- distatis(distCube)
#BootF <- BootFactorScores(testDistatis$res4Splus$PartialF,niter=1000)
#LeF	<- testDistatis$res4Splus$F
#PartialFS <- testDistatis$res4Splus$PartialF
#GraphDistatisBoot(LeF,BootF,PartialFS,ZeTitle="Bootstrap on Factors")

#df4 <- data.frame(testDistatis$res4Cmat$G,experience)

#ggplot(df4,aes(x = dim.1,y = dim.2,label=rownames(df4),color=factor(engineering))) + geom_point(aes(shape=factor(engineering)),size=5) + theme_bw()



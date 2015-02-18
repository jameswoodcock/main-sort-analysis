rm(list = ls())
library(DistatisR)
library(MASS)
library(vegan)
library(ggplot2)
library(gplots)
library(dendextend)

material = "ALL"
screePlot = FALSE
figPath = paste("./plots/",material,"/",sep="")

data <- read.csv(paste("./data/data",material,".csv",sep=""),row.names=1,header=T)

distCube <- DistanceFromSort(data)
totDist <- apply(distCube,c(1,2),sum)/11
#diag(totDist) <- 1
#dissMat <- abs(totDist-1)
clust = hclust(dist(totDist)^2,"ward")

ngroups <- apply(data,2,max)
ngroupsbar <- data.frame(participant = colnames(data),ngroups = ngroups)
mdsClusters = median(ngroups)

pdf(paste(figPath,"ngroups.pdf",sep=""),width = 10,height = 10)
print(
ggplot(ngroupsbar,aes(x=participant,ngroups)) + geom_bar(stat="identity") + labs(y="Number of groups",x="Participant") +theme_bw()
)
dev.off()

for (i in 2:11) {
pdf(paste(figPath,"dendro_",i,"groups.pdf",sep=""),width = 55, height = 8)
#plot(clust,hang=-1)
plot(clust, xlab=NA, sub=NA, main=NA, cex = 0.5)
rect.hclust(clust,k=i,border="red")
dev.off()
if (i==mdsClusters){
pdf(paste(figPath,"dendro_median_groups.pdf",sep=""),width = 55, height = 8)
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


pdf(paste(figPath,"mdsdim1dim2.pdf",sep=""),width = 10,height = 10)
#plot(dim1nm,dim2nm,type="n")
#text(dim1nm, dim2nm, labels = row.names(data), cex=.7)
print(
ggplot(df3,aes(x = dim1nm,y = dim2nm,label=rownames(df3),color=factor(groups))) + geom_point(aes(shape=factor(groups)),size=5) + theme_bw() + xlim(-0.6,0.6) + ylim(-0.6,0.6) + labs(x = "Dimension I",y = "Dimension II")+ scale_color_discrete(name = "Group",labels = groupLabs) + scale_shape_discrete(name = "Group",labels = groupLabs) + theme(legend.position = c(.9,.9))
)
#plot(nmmds,type="t")
dev.off()

pdf(paste(figPath,"mdsdim1dim2text.pdf",sep=""),width = 40, height = 40)
#plot(dim1nm,dim2nm,type="n")
#text(dim1nm, dim2nm, labels = row.names(data), cex=.7)
print(
ggplot(df3,aes(x = dim1nm,y = dim2nm,label=rownames(df3),color=factor(groups))) + geom_text(aes(shape=factor(groups))) + theme_bw() + xlim(-0.6,0.6) + ylim(-0.6,0.6) + labs(x = "Dimension I",y = "Dimension II") + scale_color_discrete(name = "Group",labels = groupLabs) + theme(legend.position = c(.9,.9))
)
#plot(nmmds,type="t")
dev.off()

pdf(paste(figPath,"mdsdim1dim3.pdf",sep=""),width = 10,height = 10)
#plot(dim1nm,dim2nm,type="n")
#text(dim1nm, dim2nm, labels = row.names(data), cex=.7)
print(
ggplot(df3,aes(x = dim1nm,y = dim3nm,label=rownames(df3),color=factor(groups))) + geom_point(aes(shape=factor(groups)),size=5) + theme_bw() + xlim(-0.6,0.6) + ylim(-0.6,0.6) + labs(x = "Dimension I",y = "Dimension III")+ scale_color_discrete(name = "Group",labels = groupLabs) + scale_shape_discrete(name = "Group",labels = groupLabs) + theme(legend.position = c(.9,.9))
)
#plot(nmmds,type="t")
dev.off()

pdf(paste(figPath,"mdsdim1dim3text.pdf",sep=""),width = 40, height = 40)
#plot(dim1nm,dim2nm,type="n")
#text(dim1nm, dim2nm, labels = row.names(data), cex=.7)
print(
ggplot(df3,aes(x = dim1nm,y = dim3nm,label=rownames(df3),color=factor(groups))) + geom_text(aes(shape=factor(groups))) + theme_bw() + xlim(-0.6,0.6) + ylim(-0.6,0.6) + labs(x = "Dimension I",y = "Dimension III") + scale_color_discrete(name = "Group",labels = groupLabs) + theme(legend.position = c(.9,.9))
)
#plot(nmmds,type="t")
dev.off()

pdf(paste(figPath,"mdsdim1dim4text.pdf",sep=""),width = 40, height = 40)
#plot(dim1nm,dim2nm,type="n")
#text(dim1nm, dim2nm, labels = row.names(data), cex=.7)
print(
ggplot(df3,aes(x = dim1nm,y = dim4nm,label=rownames(df3),color=factor(groups))) + geom_text(aes(shape=factor(groups))) + theme_bw() + xlim(-0.6,0.6) + ylim(-0.6,0.6) + labs(x = "Dimension I",y = "Dimension IV") + scale_color_discrete(name = "Group",labels = groupLabs) + theme(legend.position = c(.9,.9))
)
#plot(nmmds,type="t")
dev.off()

pdf(paste(figPath,"mdsdim2dim3text.pdf",sep=""),width = 40, height = 40)
#plot(dim1nm,dim2nm,type="n")
#text(dim1nm, dim2nm, labels = row.names(data), cex=.7)
print(
ggplot(df3,aes(x = dim2nm,y = dim3nm,label=rownames(df3),color=factor(groups))) + geom_text(aes(shape=factor(groups))) + theme_bw() + xlim(-0.6,0.6) + ylim(-0.6,0.6) + labs(x = "Dimension II",y = "Dimension III") + scale_color_discrete(name = "Group",labels = groupLabs) + theme(legend.position = c(.9,.9))
)
#plot(nmmds,type="t")
dev.off()

pdf(paste(figPath,"mdsdim3dim4text.pdf",sep=""),width = 40, height = 40)
#plot(dim1nm,dim2nm,type="n")
#text(dim1nm, dim2nm, labels = row.names(data), cex=.7)
print(
ggplot(df3,aes(x = dim3nm,y = dim4nm,label=rownames(df3),color=factor(groups))) + geom_text(aes(shape=factor(groups))) + theme_bw() + xlim(-0.6,0.6) + ylim(-0.6,0.6) + labs(x = "Dimension I",y = "Dimension IV") + scale_color_discrete(name = "Group",labels = groupLabs) + theme(legend.position = c(.9,.9))
)
#plot(nmmds,type="t")
dev.off()

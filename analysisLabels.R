rm(list = ls())
library(DistatisR)
library(MASS)
library(vegan)
library(ggplot2)
library(gplots)
library(dendextend)

material = "FF"
figPath = paste("./plots/",material,"/",sep="")

dataExp <- read.csv(paste("./data/labels_data/labels",material,"Exp.csv",sep=""),row.names=1,header=T)
dataNonExp <- read.csv(paste("./data/labels_data/labels",material,"NonExp.csv",sep=""),row.names=1,header=T)
dataAll <- read.csv(paste("./data/labels_data/labels",material,"All.csv",sep=""),row.names=1,header=T)


pdf(paste(figPath,"clustergramExp.pdf",sep=""),width = 16, height = 22)
heatmap(as.matrix(dataExp),hclustfun=function(d) hclust(d,method='ward.D2'),margins=c(20,20),col=grey(seq(0.9
,0,-0.01)),main=paste(material,"experienced listerners"))
dev.off()
pdf(paste(figPath,"clustergramNonExp.pdf",sep=""),width = 16, height = 22)
heatmap(as.matrix(dataNonExp),hclustfun=function(d) hclust(d,method='ward.D2'),margins=c(20,20),col=grey(seq(0.9
,0,-0.01)),main=paste(material,"naive listerners"))
dev.off()
pdf(paste(figPath,"clustergramAll.pdf",sep=""),width = 16, height = 22)
heatmap(as.matrix(dataAll),hclustfun=function(d) hclust(d,method='ward.D2'),margins=c(20,20),col=grey(seq(0.9
,0,-0.01)),main=paste(material,"all listerners"))
dev.off()


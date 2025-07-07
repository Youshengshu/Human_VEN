# The code is adatped from this nice paper: Brian E. Kalmbach et al., 2021, Neuron
# This code performs clustering based on ephys, and contructs various tsne plots, heatmaps featured in Figure3 

# FIRST NEED TO SET THE WORKING DIRECTORY TO THE MAIN FOLDER
# This is the folder where "Clustering_tSNE.R" is located
rm(list = ls())

###Load these libraries
library(Rtsne)
library(ggplot2)
library(dplyr)
library(sigclust)
library(gplots)
library(umap)

#load the data and list of features to be used in clustering
datadir="G:/VEN project_KW/Data analysis/Patch-seq mutimodal/tsne_classification/Data_analysis_202411/" ####set to folder folder containing data files
setwd(datadir)
datafile="Ephys.csv"
L5cells<-read.csv(paste0(datadir,datafile))

datafile="features.csv"
features<-read.csv(paste0(datadir,datafile))



#select features for clustering from csv file
L5.subset<-subset(L5cells, select=c(names(features)))
# 
##perform PCA
L5.pca<-prcomp(L5.subset,center=TRUE,scale.=TRUE)
L5PCAs<-as.matrix(L5.pca$x)
variance.explained<-summary(L5.pca)$importance[2,]
variance.explained<-subset(variance.explained, variance.explained>=0.02)

##selects only PCs that explain at least 1% of variance
L5PCAsSelect<-subset(L5PCAs, select=c(names(variance.explained)))##subset of PCAs with variance explained greater than .01


##Clustering based on PCAs和层次聚类
L5Clusters<-as.data.frame(L5PCAsSelect)
L5Clusters<-hclust(dist(L5Clusters),method="ward.D2")
dendrogram<-as.dendrogram(L5Clusters,edge.root=TRUE,leaf=FALSE)
plot(dendrogram)

cluster.p<-data.frame(p=integer())
#statistical test of clustering, cuts is the number of clusters to test, data is the df to test, outcut is the df to append the p value to
sig.cluster<-function(data){
  sigClustAll<-as.data.frame(scale(data))
  sig.result<-sigclust(sigClustAll, nsim=1000, nrep=1, icovest=1)
  plot(sig.result)
  pval<-attributes(sig.result)$pval
  return(pval)
}
L5cut<-cutree(L5Clusters,k=2)
cluster.p[1,]<-sig.cluster(data=L5.subset)#Test for sig clustering at highest level

L5cut<-cutree(L5Clusters,k=3)
TempData<-L5.subset
TempData$etype=L5cut
TempDatav2<-subset(TempData, etype=="1"|etype=="2")#Test whether IT1vET cut is significant
TempDatav2<-subset(TempDatav2, select=-c(etype))
cluster.p[2,]<-sig.cluster(data=TempDatav2)
L5cells$etype=L5cut

L5cut<-cutree(L5Clusters,k=4)
TempData<-L5.subset
TempData$etype=L5cut
TempDatav2<-subset(TempData, etype=="1"|etype=="4")#Test whether there are two types of ET neurons
TempDatav2<-subset(TempDatav2, select=-c(etype))
cluster.p[3,]<-sig.cluster(data=TempDatav2)

L5cut<-cutree(L5Clusters,k=5)
TempData<-L5.subset
TempData$etype=L5cut
TempDatav2<-subset(TempData, etype=="2"|etype=="4")#Test whether there are two types of IT1 neurons
TempData2v2<-subset(TempDatav2, select=-c(etype))
cluster.p[4,]<-sig.cluster(data=TempDatav2)
rownames(cluster.p)<-c("2 clusters", "3 clusters", "4 clusters", "5 clusters")
print(cluster.p)

L5cut<-cutree(L5Clusters,k=3) 
L5cells$etype=L5cut


##plot tsne for Figure 2 and Fig S3########
#Calculate tsne
set.seed(40)

tsneout<-Rtsne(L5PCAsSelect, perplexity=15,step=1000,pca=FALSE,theta = 0.2)
L5cells$tsne1.1=(tsneout$Y[,1])/max(tsneout$Y[,1])
L5cells$tsne1.2=(tsneout$Y[,2])/max(tsneout$Y[,2])

#tsne plot function
tsne.plot.function<-function(FileName,data,x,y,var,colors,labels){
  .e <- environment()
  print(output.plot<-ggplot(data, aes(x,y),environment = .e)+
          geom_point(size=3,shape=16,colour="transparent")+geom_point(aes(colour=factor(var)))+
          theme_classic()+scale_colour_manual(values=colors,name="e-type",labels=labels)+theme(aspect.ratio=1)+
          coord_cartesian(xlim=c(-1.5,1.5),ylim=c(-1.5,1.5)))
  ggsave(FileName,plot=output.plot,width=5,height=5)
  return(output.plot)
}

# #Plot e-type tsne
FileName<-paste0(datadir,"/etype_tsne_test.pdf")
var1<-L5cells$etype
colors<-c("green","red","cyan",'gray','yellow')
labels<-c("1","2","3","4")

x<-L5cells$tsne1.1
y<-L5cells$tsne1.2
etype.tsne<-tsne.plot.function (FileName=FileName,data=L5cells,x=x,y=y,var=var1,colors=colors,labels = labels)


# #Plot M-type tsne
datafile="label.csv"
label<-read.csv(paste0(datadir,datafile))
label<-subset(label, select=c(Label))
label_mtype<-as.matrix(label)

FileName<-paste0(datadir,"/Mtype_tsne_test.pdf")
var2<-label_mtype
colors<-c("green","red","cyan",'gray','yellow')
labels<-c("VEN","ETPClike","TRI","ITPClike")
x<-L5cells$tsne1.1
y<-L5cells$tsne1.2
etype.tsne<-tsne.plot.function (FileName=FileName,data=L5cells,x=x,y=y,var=var2,colors=colors,labels = labels)
etype.tsne



#heatmap for showing parameter, no need to do cluster
my_palette<-colorRampPalette(c("blue","white","red"))(n=75)
L5cells$etype <- as.numeric(L5cells$etype)
L5.subset$etype=L5cells$etype

L5.subset<-L5.subset[order(L5.subset[,ncol(L5.subset)],decreasing=T),]
scale.L5.subset<-scale(L5.subset)
scale.L5.subset.matrix<-as.matrix(scale.L5.subset)
L5heat<-heatmap.2(scale.L5.subset.matrix,Rowv=FALSE,Colv=TRUE, dendrogram = "none",col=my_palette, trace="none",key=TRUE,keysize=2,margins=c(10,10))
#L5heat<-heatmap.2(scale.L5.subset.matrix,Rowv=TRUE,Colv=FALSE, dendrogram = "row",col=my_palette, trace="none",key=TRUE,keysize=2, hclustfun=function(x) hclust(x,method="ward.D2"),margins=c(10,10))
write.csv(L5cells,paste0(datadir,"/clusteringOutput.csv"))

















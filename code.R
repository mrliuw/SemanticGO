setwd("/raid1/LSA/LSA-Arabi")
options(stringsAsFactors = FALSE)
dat=scan("tair.gaf",what=character(),sep="\t",quote="")
library(stringi)
newdat=data.frame(gene=sub("\\|.*", "", dat[c(seq(11,4247212,17))]),go=dat[c(seq(5,4247212,17))])
dim(newdat)
newdat2=newdat[-grep("RNA",newdat[,1]),]
write.csv(newdat2,"gene-go.csv",row.names=F)
dim(newdat2)
names(newdat)
test=data.frame(newdat2,count=(rep(0,dim(newdat2)[1])))
dim(test)
ngene=length(unique(test[,1]))
ngo=length(unique(test[,2]))
ngene;ngo
testma=matrix(0,ngene,ngo)
dim(testma)
dimnames(testma) <- list(unique(test[,1]), unique(test[,2]))
str(testma)
genes=unique(test[,1])
gos=unique(test[,2])
result=acast(test, gene~go, value.var="count")
for (n in 1:dim(newdat2)[1]) {write.table(data.frame(gene=test[n,1],go=get_parent_nodes(test[n,2])[2],row.names=NULL),"childgo.csv",row.names=F,col.names=F,append=T)}
#above codes too slow
library(doParallel)
m=detectCores()
cl <- makeCluster(20)
registerDoParallel(cl)
foreach(i=1:20,.packages = c("GOfuncR")) %dopar% {
  for (n in (dim(newdat2)[1]-12491*(i-1)):(dim(newdat2)[1]-12491*i)) {write.table(data.frame(gene=newdat2[n,1],go=get_parent_nodes(newdat2[n,2])[2],row.names=NULL),paste(i,"-childgo.csv",sep=""),row.names=F,col.names=F,append=T,sep=",")}
}
stopCluster(cl)
alldat=0
num=0
for (i in 1:20) { dat=read.csv(paste(i,"-childgo.csv",sep=""),sep=",",header=F);alldat=rbind(alldat,dat,make.row.names=F);num=num+length(dat[,1])}
print(num)
colnames(alldat)=c("gene","go")
finaldat=rbind(alldat,newdat2,make.row.names=F)
fin=unique(finaldat)
final=data.frame(fin[-1,],count=rep(1,dim(fin)[1]-1))
library(reshape2)
finalresult=acast(final, gene~go, value.var="count") #the original matrix
finalresult[is.na(finalresult)] <- 0
write.csv(finalresult,"final-gene-go.csv")

##calculate MAE(=abs(original-predicted)) for k=50,100,150...250
finalresult=read.csv("final-gene-go.csv",row.names=1)
library(lsa)
space=lsa(finalresult,dims=300)
plot(space$sk,pch=20,cex=0.2,ylab="Singular value",xlab="Indices of singular value",
     main="Singular values of Arabidopsis gene-GO matrix")
D <- diag(space$sk)
##add a colname for the column of rownames before melting
finalresult2=cbind(gene = rownames(finalresult), finalresult)
rownames(finalresult) <- NULL
originalGFcolums=melt(finalresult2,value.name = "count", varnames=c('gene', 'go'))
for (k in seq(50,250,50)){
newGF=space$tk[,1:k] %*% D[1:k,1:k] %*% t(space$dk[,1:k]) #the new matrix, we can check the largest value for each go, maybe the biggest gene
MAE<-sum(abs(finalresult-newGF))/dim(finalresult)[1]/dim(finalresult)[2]
print(MAE)
newGFcolums=melt(newGF,value.name = "count", varnames=c('gene', 'go'))
combine=cbind(originalGFcolums,newcount=newGFcolums[,-c(1,2)])
tp=0;fp=0;fn=0;tn=0;error=0 #too slow, use dopar, it comsume time,debug before running 
cl <- makeCluster(4) #it comsume memory, donnot use too much cores
registerDoParallel(cl)
points=100
foreach(T=1:points,.combine=rbind) %dopar% {
  threshold=T/100 #1000 will take longer time
  tp=length(which(combine$count==1 & combine$newcount>threshold)); 
  fp=length(which(combine$count==1 & combine$newcount<threshold));  
  fn=length(which(combine$count==0 & combine$newcount>threshold));
  tn=length(which(combine$count==0 & combine$newcount<threshold));
  error=fp+fn;
  write.table(data.frame(tp=tp,fp=fp,fn=fn,tn=tn,error=error,threshold=threshold,tnr=tn/(fp+tn),
              tpr=tp/(tp+fn),f1=2*tp*tp/(tp+fp)/(tp+fn)/(tp/(tp+fp)+tp/(tp+fn))),paste("k",k,"4ROC.csv",sep="-"),row.names=F,col.names=F,sep=",",append=T)
}
stopCluster(cl)
}
#ploting ROC curves#
k=200
roc=read.csv(paste("k",k,"4ROC.csv",sep="-"),header=F,sep=",")
names(roc)=c("tp","fp","fn","tn","error","threshold","tnr","tpr","fmeasure")
roc=roc[order(roc$threshold,decreasing = F),]
plot(seq(1,points,1)/100,roc$error,type="l",col="red",xlab="Threshold",ylab="Num. of errors")
lines(seq(1,points,1)/100,roc$fp,col="green")
lines(seq(1,points,1)/100,roc$fn,col="black")
sens=roc$tp/(roc$tp+roc$fn)
spec=roc$tn/(roc$fp+roc$tn)
plot(1-spec,sens,type="l",xlab="1-Specificity",ylab="Sensitivity")
height = (sens[-1]+sens[-length(sens)])/2
width = -diff(spec) # = diff(rev(omspec))
sum(height*width)+(1-sum(width))*1
choseThreshold=roc[which.min(roc$error),]$threshold

setwd("/raid1/LSA/LSA-Arabi/revision")
fmeasure=read.csv("Fmeasure.csv")
mae=read.csv("MAE.csv")
plot(mae,type="b",lty=2,pch=16,ylab="Mean absolute error (MAE)",xlab="Rank k",
     main="Model-free mean absolute error at different k")
plot(fmeasure,type="b",lty=2,pch=16,ylab="F-measure",xlab="Rank k",
     main="Model-free performance at different k")
##########10-fold cross validation based k=200###########
finalresult=read.csv("final-gene-go.csv",row.names=1)
library(lsa)
k <- 10
indices <- sample(1:nrow(finalresult))
folds <- cut(indices, breaks = k, labels = FALSE)
for (i in 1:k) {
  cat("processing fold #", i, "\n")
  val_indices <- which(folds == i, arr.ind = TRUE)
  val_data <- finalresult[val_indices,]
  partial_train_data <- finalresult[-val_indices,]
  space=lsa(t(partial_train_data),dims=200) #fold_in new gene should be vectors
  pre_val=fold_in(t(as.matrix(val_data)),space) #transpose and comfirm val_data is.numeric,sample is row need transpose
  class(pre_val)="matrix"
  val_columns=melt(t(as.matrix(val_data)),value.name = "count", varnames=c('go', 'gene')) #transpose?
  pre_columns=melt(pre_val,value.name = "count", varnames=c('go', 'gene')) #check before melt
  combine=cbind(val_columns,newcount=pre_columns[,-c(1,2)])
  tp=0;fp=0;fn=0;tn=0;error=0 #too slow, use dopar, it comsume time,debug before running 
  cl <- makeCluster(4) #it comsume memory, donnot use too much cores
  registerDoParallel(cl)
  points=100
  foreach(T=1:points,.combine=rbind) %dopar% {
    threshold=T/100 #1000 will take longer time
    tp=length(which(combine$count==1 & combine$newcount>threshold)); 
    fp=length(which(combine$count==1 & combine$newcount<threshold));  
    fn=length(which(combine$count==0 & combine$newcount>threshold));
    tn=length(which(combine$count==0 & combine$newcount<threshold));
    error=fp+fn;
    write.table(data.frame(tp=tp,fp=fp,fn=fn,tn=tn,error=error,threshold=threshold,tnr=tn/(fp+tn),
                           tpr=tp/(tp+fn),f1=2*tp*tp/(tp+fp)/(tp+fn)/(tp/(tp+fp)+tp/(tp+fn))),paste("fold",i,"100-4ROC.csv",sep="-"),row.names=F,col.names=F,sep=",",append=T)
  }
  stopCluster(cl)
}
###foldin rice and calculate performance statistics##
a=read.csv("/raid1/LSA/LSA-rice/rice4foldin.csv",header=T,row.names = 1)
space=lsa(t(partial_train_data),dims=200)
pre_rice=fold_in(t(as.matrix(a)),space) 
class(pre_rice)="matrix"
## rice_genevec=as.matrix(a) %*% space$tk;write.csv(rice_genevec,"rice_genevec.csv",header=T)
rice_columns=melt(t(as.matrix(a)),value.name = "count", varnames=c('go', 'gene')) #transpose?
pre_columns=melt(pre_rice,value.name = "count", varnames=c('go', 'gene')) #check before melt
combine=cbind(rice_columns,newcount=pre_columns[,-c(1,2)])
choseThreshold=0.73
novelgo-rice=combine[which(combine$count==0 & combine$newcount>choseThreshold),]
write.csv(novelgo,"novelgo-rice.csv",row.names=F)
tp=0;fp=0;fn=0;tn=0;error=0 #too slow, use dopar, it comsume time,debug before running 
cl <- makeCluster(4) #it comsume memory, donnot use too much cores
registerDoParallel(cl)
points=100
foreach(T=1:points,.combine=rbind) %dopar% {
  threshold=T/100 #1000 will take longer time
  tp=length(which(combine$count==1 & combine$newcount>threshold)); 
  fp=length(which(combine$count==1 & combine$newcount<threshold));  
  fn=length(which(combine$count==0 & combine$newcount>threshold));
  tn=length(which(combine$count==0 & combine$newcount<threshold));
  error=fp+fn;
  write.table(data.frame(tp=tp,fp=fp,fn=fn,tn=tn,error=error,threshold=threshold,tnr=tn/(fp+tn),
                         tpr=tp/(tp+fn),f1=2*tp*tp/(tp+fp)/(tp+fn)/(tp/(tp+fp)+tp/(tp+fn))),paste("rice","4ROC.csv",sep="-"),row.names=F,col.names=F,sep=",",append=T)
}
stopCluster(cl)

###rice pathway analysis###
osapathway=read.csv("osapathgenes.csv",header=F)
osagenevec=read.csv("rice_genevec.csv",header=T,row.names=1)
osapathwaynames=unique(osapathway[,2]) #delet the empty one
osafinalsumpathway=data.frame()
for (p in osapathwaynames){
  osapathwaygene=osapathway[osapathway$V2==p,1]
  osapathwayvec=osagenevec[marray::rm.na(match(unlist(as.list(t(osapathwaygene))),rownames(osagenevec))),]
  osasumpathway=apply(osapathwayvec,2,sum)
  osafinalsumpathway=rbind(osafinalsumpathway,osasumpathway) #rbind的参数先后顺序有差别
}
osapathwaycos=cosine(as.matrix(t(osafinalsumpathway)))
dimnames(osapathwaycos)=list(osapathwaynames,osapathwaynames)
write.csv(osapathwaycos,"osapathwaycos.csv")
tree=hclust(as.dist(1-osapathwaycos),method="average")
plot(tree, main = "Clustering of metabolism pathways",xlab = "", sub = "",cex=0.5)
row.names(osafinalsumpathway) =osapathwaynames
write.table(data.frame(osafinalsumpathway),"osapathwayvec.csv",sep=",",col.names=F,append=F)
heatmap(osapathwaycos,labRow=osapathwaynames,symm = T,labCol=F,distfun=dist,
        hclustfun = hclust,margins = c(1, 4),cexRow=0.5,col = heat.colors(12)) 


genevec=space$tk[,1:k] %*% D[1:k,1:k]
govec=D[1:k,1:k] %*% t(space$dk[,1:k])
#####ploting ROC curves###
roc=read.csv("4ROC.csv",header=F,sep=",")
names(roc)=c("tp","fp","fn","tn","error","threshold","tnr","tpr")
roc=roc[order(roc$threshold,decreasing = F),]
plot(seq(1,points,1)/100,roc$error,type="l",col="red",xlab="Threshold",ylab="Num. of errors")
lines(seq(1,points,1)/100,roc$fp,col="green")
lines(seq(1,points,1)/100,roc$fn,col="black")
sens=roc$tp/(roc$tp+roc$fn)
spec=roc$tn/(roc$fp+roc$tn)
plot(1-spec,sens,type="l",xlab="1-Specificity",ylab="Sensitivity")
height = (sens[-1]+sens[-length(sens)])/2
width = -diff(spec) # = diff(rev(omspec))
sum(height*width)+(1-sum(width))*0.9178
choseThreshold=roc[which.min(roc$error),]$threshold
novelgo=combine[which(combine$count==0 & combine$newcount>choseThreshold),]
write.csv(novelgo,"novelgo.csv",row.names=F)
###########test the normality
library("ggpubr")
ggdensity(novelgo$newcount, 
          main = "Density plot of GO",
          xlab = "GO value")
######narrow down by threshold
novelgo8=combine[which(combine$count==0 & combine$newcount>0.8),];write.csv(novelgo8,"novelgo8.csv",row.names=F)
# Create a cosine similarity between two Terms
library(LSAfun)
#genevec=read.table("genevec.csv",header=T,row.names=1,sep=",")
#myCo <- costring('ACL1','ACL2', tvectors= genevec) 
#cosinemat=cosine(t(as.matrix(genevec))) #too slow using the dopar

cl <- makeCluster(3) #set up as many as cores,n can be zhengchu?too large error shows up
registerDoParallel(cl)  
n=dim(genevec)[1]
out=0
foreach(i=1:n,.combine=rbind,.packages = c("LSAfun")) %dopar%  {
  for (j in 1:n){out[j]=as.matrix(round(cosine(unlist(genevec[i,]),unlist(genevec[j,])),2))
  }
  write.table(t(out),"out.csv",row.names=i,col.names=F,sep=",",append=T) #need to be ordered
}
cosinemat=read.csv("out.csv",header=F)
stopCluster(cl) 

da=read.csv("hsp20cosine.csv",header=T,row.names=1)
hclust(as.dist(1-da),method="average")
tree=hclust(as.dist(1-da),method="average")
plot(tree, main = "Clustering of Hsp20 genes",xlab = "", sub = "")

options(stringsAsFactors = FALSE)
dat0=read.csv("genevec.csv",header=T)
datSummary=dat0[,1:1];dim(dat0)

datExpr = t(dat0[,2: ncol(dat0)]);no.samples = dim(datExpr)[[1]];dim(datExpr)
ArrayName= names(data.frame(dat0[,-1]))
GeneName= dat0$geneid
datExpr=data.frame(t(dat0[,-1]))
names(datExpr)=dat0[,1]
dimnames(datExpr)[[1]]=names(data.frame(dat0[,-1]))
table(dimnames(datExpr)[[1]]==datTraits$sample)
y=datTraits$gender
z=datTraits$age
sizeGrWindow(9, 5)
pdf("ClusterTreeSamples.pdf")
plotClusterTreeSamples(datExpr=datExpr, y=y)
dev.off()
rm(dat0);gc()

memory.size=100000
powers=c(seq(1,10,by=1),seq(12,14,by=2));sft=pickSoftThreshold(datExpr, powerVector=powers,networkType = "signed")
RpowerTable=sft[[2]]
sizeGrWindow(9, 5);pdf('choosing power.pdf');par(mfrow = c(1,2));cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red");dev.off()
# Mean connectivity as a function of the soft-thresholding power
sizeGrWindow(9, 5);pdf('mean connectivity.pdf');plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red");dev.off()

softPower =12
Connectivity=softConnectivity(datExpr,corFnc = "cor", corOptions = "use ='p'",power=softPower,type="signed")
pdf("scale-free.pdf");scaleFreePlot(Connectivity,nBreaks = 10,truncated = FALSE,removeFirst = FALSE, main = "");dev.off()
adjacency = adjacency(datExpr,corFnc = "cor", corOptions = "use ='p'",type = "signed", power = softPower)
TOM = TOMsimilarity(adjacency,TOMType="signed");dissTOM = 1-TOM
#method="complete"  ?
geneTree = hclust(as.dist(dissTOM), method = "average")#高版本已经用hclust
minModuleSize =30;
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = 1, pamRespectsDendro = FALSE,minClusterSize = minModuleSize,cutHeight=0.99);
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

MEList = moduleEigengenes(datExpr, colors = dynamicMods)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average");#
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")
MEDissThres = 0.2
abline(h=MEDissThres, col = "red")
merge = mergeCloseModules(datExpr, dynamicMods, cutHeight = MEDissThres, verbose = 3);
mergedColors = merge$colors;
mergedMEs = merge$newMEs;
sizeGrWindow(12, 9)
pdf("DendroAndColors.pdf")
plotDendroAndColors(geneTree, cbind(dynamicMods, mergedColors),c("Dynamic Tree Cut", "Merged dynamic"),dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleColors = mergedColors
colorOrder = c("grey", standardColors(unique(moduleColors)));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average");#
pdf("METree.pdf")
plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")
dev.off()
MEList = moduleEigengenes(datExpr, colors = dynamicMods)
nSamples=nrow(datExpr)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = cbind.data.frame(datSummary,corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
write.table(data.frame(ArrayName,MEs),"MEs.csv",row.name=F)
kMEdat=data.frame(geneModuleMembership,MMPvalue)
write.table(data.frame(datSummary,kMEdat),"kME-MMPvalue.csv",row.names=FALSE)
k.in=intramodularConnectivity(adjacency(datExpr,corFnc = "cor", corOptions = "use ='p'", type = "signed", power = softPower), moduleColors,scaleByMax = FALSE)
datout=data.frame(datSummary, colorNEW=moduleColors, k.in)
write.table(datout, file="OutputCancerNetwork.csv", sep=",", row.names=F)
hubs    = chooseTopHubInEachModule(datExpr, moduleColors)
write.csv(data.frame(module=names(hubs),moduleColor=labels2colors(names(hubs)),hub=hubs),"num2color.csv",row.names=F)
###compare hsp similarity with other 18 sampled gene similarity
hsp=c("AT5G47600","AT4G21870","AT5G37670","AT1G54050","AT5G12020","AT1G59860","AT2G29500","AT1G53540","AT5G12030","AT1G07400",
      "AT5G59720","AT2G19310","AT5G54660","AT4G10250","AT5G51440","AT4G25200","AT4G27670","AT1G52560")
hspcos=cosine(t(as.matrix(genevec[match(toupper(hsp[1:18]),rownames(genevec)),])))
hspcos=hspcos[lower.tri(hspcos)]
samplevec=0
res=0
for (i in 1:1000) {
  samplevec=genevec[sample(nrow(genevec),size=18,replace=F),]
  samplecos=cosine(t(as.matrix(samplevec)))
  samplecos=samplecos[lower.tri(samplecos)]
  res=rbind(res,data.frame(pvalue=t.test(hspcos,samplecos)$p.value,x=t.test(hspcos,samplecos)$estimate[1],y=t.test(hspcos,samplecos)$estimate[2])) #p-value < 2.2e-16,mean of x mean of y 0.7369939 0.3771999 
}
count=res[2:1001,]$pvalue>0.01
length(count[count==TRUE])/1000

##KEGG pathway representations analysis,replace ath01212
setwd("/raid1/LSA/LSA-Arabi/revision")
genevec=read.csv("genevec-200.csv",header=T,row.names = 1)
pathwayvec=read.csv("pathwayvec-200.csv",header=F,row.names = 1)
pathway=read.csv("/raid1/LSA/LSA-Arabi/pathgenes.csv",header=F)
pathwaygenecos=sort(cosine(unlist(pathwayvec["ath00920",]),as.matrix(t(genevec))),decreasing = T)
commongene=intersect(names(pathwaygenecos)[1:30],pathway[pathway$V2=="ath00920",3])
uniquegene=names(pathwaygenecos)[1:30][!names(pathwaygenecos)[1:30] %in% commongene]
randompathwaygenecos=cosine(unlist(pathwayvec["ath00920",]),as.matrix(t(genevec)))
p=lapply(1:1000,function(x) wilcox.test(pathwaygenecos[uniquegene],sample(randompathwaygenecos,20))$p.value)
sum(p<0.000001)

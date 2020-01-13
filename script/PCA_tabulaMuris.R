
## PCA 

cats<-c("known risk genes","candidate risk genes","random genes", "other genes")
pcols<-c("red","dodgerblue",rgb(0,0,0,0.4),rgb(220/255,220/255,220/255,0.5)) #"sienna2"
allg<-read.table("~/Dropbox (CGC)/PAH/US_UK/result/geneExpression/MouseHumanIDmapper.txt",
                 header = 1,comment.char = "",sep="\t",stringsAsFactors = F)

knowngenes<-"SOX17|BMPR2|TBX4|KCNK3|ATP13A3|GDF2|AQP1$|KLK1$|GGCX|ACVRL1|^KDR$|^ENG$|EIF2AK4|SMAD4|^CAV1$|NOTCH1$|SMAD4|SMAD9"
confirmgenes<-"SOX17|BMPR2|TBX4|KCNK3|ACVRL1|^ENG$|EIF2AK4|^CAV1$|NOTCH1$|SMAD4|SMAD9"


#risks<-"col6a5|fbln2|pdgfd|hn1l|jpt2|^kdr$" #chd9|mrrf|ccdc28A|vps51|mmp10|chd9|CHKB|myf5|irak4|npffr2|myo5b|ipo4|hn1l|jpt2"
risks<-"fbln2|pdgfd|^kdr$" #chd9|mrrf|ccdc28A|vps51|mmp10|chd9|CHKB|myf5|irak4|npffr2|myo5b|ipo4|hn1l|jpt2"
denovos<-read.csv("~/Dropbox (CGC)/PAH/US_UK_SPARK/de novo data/ALL_PAH_denovos.csv",header=1,stringsAsFactors = F,comment.char = "",check.names = F)
denovogenes<-denovos$GeneName[grep("^SYN",denovos$Mutation_Type)]
denovo_mouse<-allg$MouseGene[which(allg$HumanGene%in%denovogenes)]
#cands<-read.table("~/Dropbox (CGC)/PAH/US_UK/result/geneExpression/candidates.txt",header = F,comment.char = "",stringsAsFactors = F,fill = T)
#cands<-cands$V1

all_negatives<-allg[intersect(which(!allg$HumanGene%in%c(knowngenes,risks)==T),
                              grep(paste(knowngenes,risks,sep="|"),allg$HumanGene,ignore.case = T,invert = T)),"MouseGene"]



dat_all<-read.csv("~/Dropbox (CGC)/PAH/US_UK/result/geneExpression/FACS_tissues_genes.summary.csv",
              header=1,stringsAsFactors = F,comment.char = "",check.names = F)
dat_all<-dat_all[which(dat_all$gene%in%allg$MouseGene),]
#dat_all<-dat_all[grep(".y$",names(dat_all),invert = T),]
#dat<-dat[which(dat$gene%in%allg$V4),]
#ss<-unlist(lapply(1:dim(dat_all)[1],FUN=function(i){sum(dat_all[i,2:dim(dat_all)[2]])}))
#dat_all<-dat_all[which(ss>0),]

all_negatives<-all_negatives[which(all_negatives%in%dat_all$gene)]
negatives<- read.table("~/Dropbox (CGC)/PAH/US_UK/result/geneExpression/negative_selected.csv",header=1,stringsAsFactors = F,comment.char = "") #sample(all_negatives,100,replace = F)
negatives<-negatives$gene
includes<-unique(c(which(dat_all$gene%in%c(negatives)),grep(risks,dat_all$gene,ignore.case = T),
                   grep(knowngenes,dat_all$gene,ignore.case = T)))
indexs<-intersect(grep("heart|lung|aorta",names(dat_all),ignore.case = T) ,
                  grep("Lung_leukocyte|lung_$|heart_$|killer|_T|_B|smooth|antigen|trache|lung_monocyte|Aprta_fib",names(dat_all),ignore.case = T,invert = T))
dat<-dat_all[includes,]

dat_project<-dat_all[-includes,]

pcs<-prcomp(as.matrix(dat[,indexs]),scale. = T,center = T)
proj<-scale(dat_project[,indexs],pcs$center,pcs$scale) %*% pcs$rotation

n=5
cols<-colors()
cols<-sample(cols,size = n,replace = F)
#cols[c(7,8,2)]<-"white"

cells<-rownames(pcs$rotation)
cells<-gsub("rate_","",cells)
cells<-gsub("_"," ",cells)
cells<-gsub("Lung lung","Lung",cells)
cells<-gsub("of lung","",cells)

#km<-kmeans(rbind(pcs$x,proj),n)
matri<-as.matrix(dat[,indexs])
row.names(matri)<-dat$gene
colnames(matri)<-cells
matri<-matri[c(grep(confirmgenes,dat$gene,ignore.case = T),
           grep(risks,dat$gene,ignore.case = T),which(dat$gene%in%c(negatives))),]
#colnames(matri)<-gsub("rate_| cell","",colnames(matri))
pdf("~/Dropbox (CGC)/PAH/US_UK/result/geneExpression/Heatmap_no_jpt2_col6A5_2.pdf",width=10,height=10)
par(mar=c(2,5,15,2))
library(RColorBrewer)

colMain <- colorRampPalette(brewer.pal(8, "Blues"))(25)
heatmap(matri[c(grep(confirmgenes,rownames(matri),ignore.case = T),
              grep(risks,rownames(matri),ignore.case = T)),],
        cexRow = 1.5,cexCol = 1.3,col=colMain,labCol = '')
par(oma=c(0,0,0,0))
heatmap(matri,cexRow = 0.5,cexCol = 1.3,col=colMain)
dev.off()
mat<-data.frame(gene=c(dat$gene,dat_project$gene),PC2=c(pcs$x[,2],proj[,2]))
mat$human_gene<-unlist(lapply(1:dim(mat)[1],FUN=function(i){id<-which(allg$MouseGene==mat$gene[i]); if(length(id)>0){return(allg$HumanGene[id[1]])}else{return(NA)}}))
mat$transID<-unlist(lapply(1:dim(mat)[1],FUN=function(i){id<-which(allg$MouseGene==mat$gene[i]); if(length(id)>0){return(allg$transID[id[1]])}else{return(NA)}}))

mat$type=cats[4]
mat$type[which(mat$gene %in%negatives)]<-cats[3]
mat$type[which(mat$gene %in%denovo_mouse)]<-"de novo"
mat$type[grep(confirmgenes,mat$human_gene,ignore.case = T)]<-cats[1]

mat$type[grep(risks,mat$human_gene,ignore.case = T)]<-cats[2]
if(mat$PC2[which(mat$gene=="Kdr")] <0){
  mat<-mat[order(mat$PC2,decreasing = F),]
}else{
  mat<-mat[order(mat$PC2,decreasing = T),]
}
mat$rank<-1:dim(mat)[1]
mat<-merge(mat,dat_all[,c(1,indexs)],by.x="gene",by.y="gene",all.x=T)
write.csv(mat,"~/Dropbox (CGC)/PAH/US_UK/result/geneExpression/PC2_pergene_no_jpt2_col6A5.csv")
#cols<-sample(colors(),size = 15,replace = F)
cex.tex=1.3
pdf("~/Dropbox (CGC)/PAH/US_UK/result/geneExpression/PCA_projection_no_jpt2_col6A5.pdf",width=7,height=7)

## PCA
i=1;j=1;

plot(proj[,j],proj[,i+1],xlab="",
     ylab="",col=pcols[4],main="",
     pch=20,cex=1,xlim=c(min(proj[,j])*1.5,max(proj[,j])*1.1),xaxt='n',yaxt='n')
axis(1,at = seq(-3,15,3),labels = seq(-3,15,3),cex.axis=cex.tex)
axis(2,at = seq(-4,6,2),labels = seq(-4,6,2),cex.axis=cex.tex)
title(xlab="PC1",ylab="PC2",line=2.5,cex.lab=cex.tex)
#paste(format(100*pcs$sdev[i+1]/sum(pcs$sdev),trim = 2,digits = 2),"% ",
#format(100*pcs$sdev[j]/sum(pcs$sdev),trim = 2,digits = 2),"% ",
points(pcs$x[,j],pcs$x[,i+1],xlab=paste("PC",j,sep=""),ylab=paste("PC",i+1,sep=""),main="",
       pch=20,cex=1,col=pcols[3])
points(pcs$x[grep(confirmgenes,dat$gene,ignore.case = T),j],pcs$x[grep(confirmgenes,dat$gene,ignore.case = T),
                                                                  i+1],
       col=pcols[1],
       pch=20,cex=1.2)
points(pcs$x[grep(risks,dat$gene,ignore.case = T),j],pcs$x[grep(risks,dat$gene,ignore.case = T),i+1],
       col=pcols[2],
       pch=20,cex=1.2)
# points(pcs$x[which(dat$gene%in%denovo_mouse),j],pcs$x[which(dat$gene%in%denovo_mouse),i+1],
#        col="orange",
#        pch=20,cex=1)

ids<-intersect(grep(confirmgenes,dat$gene,ignore.case = T),which(pcs$x[,j] >0))
text(pcs$x[ids,j],
     pcs$x[ids,i+1],col=pcols[1],
     dat$gene[ids],cex = cex.tex,pos=4,font=4)

ids<-intersect(grep(confirmgenes,dat$gene,ignore.case = T),which(pcs$x[,j] <0))
text(pcs$x[ids,j],
     pcs$x[ids,i+1],col=pcols[1],
     dat$gene[ids],cex = cex.tex,pos=c(4,2,1,3,4,4),font=4)

ids<-intersect(grep(risks,dat$gene,ignore.case = T),which(pcs$x[,j] >0))
text(pcs$x[ids,j],pcs$x[ids,i+1],
     col=pcols[2],
     dat$gene[ids],cex = cex.tex,pos=4,font=4)


ids<-intersect(grep(risks,dat$gene,ignore.case = T),which(pcs$x[,j] <0))
text(pcs$x[ids,j],pcs$x[ids,i+1],
     col=pcols[2],
     dat$gene[ids],cex = cex.tex,pos=c(3,4,4),font=4)

legend("bottomright",legend = cats,fill=pcols,bty='n',cex=1.2)
dev.off()

## loadings

pdf("~/Dropbox (CGC)/PAH/US_UK/result/geneExpression/PCA_loadings_no_jpt2_col6A5.pdf")
i=2
par(mar=c(2.5,4,1,1))
plot(1:dim(pcs$rotation)[1],pcs$rotation[,i],
     xlim=c(-0.7,dim(pcs$rotation)[1]+2),ylim=c(-0.45,0.5),ylab="",xaxt='n',yaxt='n',xlab="", pch=20)
axis(2,seq(-0.4,0.4,0.2),seq(-0.4,0.4,0.2),cex.axis=cex.tex)
title(xlab="Cells",line=1,cex.lab=cex.tex)
title(ylab="PC2 loadings",line=2.5,cex.lab=cex.tex)
text((dim(pcs$rotation)[1]-4):dim(pcs$rotation)[1],pcs$rotation[(dim(pcs$rotation)[1]-4):dim(pcs$rotation)[1],i], 
     cells[(dim(pcs$rotation)[1]-4):dim(pcs$rotation)[1]],cex=1.2,pos=1,srt=0,pch=20)

text(4:(dim(pcs$rotation)[1]-5),pcs$rotation[4:(dim(pcs$rotation)[1]-5),i], cells[4:(dim(pcs$rotation)[1]-5)],cex=1.2,pos=3,srt=0,pch=20)

text(1:3,pcs$rotation[1:3,i], cells[1:3],cex=1.2,pos=1,srt=0,pch=20)
par(mar=c(4.5,4.5,1,1))
## density
plot(density(pcs$x[,i],from = -5,to=5),main="",xlab="",col="gray",ylim=c(0,0.9),ylab="",yaxt='n',xaxt='n')
title(xlab="PC2",ylab="Density",line=2.5,cex.lab=cex.tex)
axis(1,at=seq(-4,4,2),labels = seq(-4,4,2),cex.axis=cex.tex)
axis(2,at=seq(0,0.8,0.2),labels = seq(0,0.8,0.2),cex.axis=cex.tex)
lines(density(pcs$x[grep(confirmgenes,dat$gene,ignore.case = T),i],from = -5,to=5),col=pcols[1])
lines(density(pcs$x[grep(risks,dat$gene,ignore.case = T),i],from = -5,to=5),col=pcols[2])
# lines(density(pcs$x[which(dat$gene%in%denovo_mouse),i]),col="orange")
#lines(density(pcs$x[which(dat$gene%in%negatives),i],from = -5,to=5),col="black")
legend("topleft",legend = cats[c(1,2,4)],cex=cex.tex,
       fill = pcols[c(1,2,4)],bty='n')



#,grep(risks,mat$gene,ignore.case = T)
hist(mat$rank[c(grep(confirmgenes,mat$gene,ignore.case = T))],
     xlab="Rank of known risk genes",main="Rank distribution",ylim=c(0,6),cex.lab=cex.tex,cex.axis=cex.tex,
     breaks = c(0,1000,2000,4000,6000,8000,10000,12000,13907),freq = T)

par(mar=c(1,5,1,1))
xx<-mat$rank[c(grep(confirmgenes,mat$gene,ignore.case = T),grep(risks,mat$gene,ignore.case = T))]
boxplot(mat$rank,yaxt='n')
axis(2,seq(0,14000,2000),seq(0,14000,2000),cex.axis=cex.tex)
title(ylab="Rank of known and candidate risk genes",cex.lab=cex.tex)
xx<-mat$rank[c(grep(confirmgenes,mat$gene,ignore.case = T))]
points(c(1.015,0.985,1,0.97,1.025,1,1,rep(1,4)), sort(xx),pch=20,cex=1.2,col=pcols[1])
xx<-mat$rank[c(grep(risks,mat$gene,ignore.case = T))]
#points(c(0.94,1,0.94,1,1,1), sort(xx),pch=15,cex=1,col=pcols[2])
points(c(0.96,0.985,0.96), sort(xx),pch=15,cex=1,col=pcols[2])
legend("topleft",legend = cats[1:2],
       col= pcols[1:2],pch=c(20,15),bty='n',cex=cex.tex)

dev.off()

pdf("~/Dropbox (CGC)/PAH/US_UK/result/geneExpression/PCA_all_no_jpt2_col6A5.pdf")
 for(i in 1:10){
 plot(1:dim(pcs$rotation)[1],pcs$rotation[,i],xlim=c(0,dim(pcs$rotation)[1]),ylim=c(-0.7,0.7),ylab=paste("PC",i),xlab="cells")
 text(1:dim(pcs$rotation)[1],pcs$rotation[,i], cells,cex=0.75,pos=3,srt=45,pch=20)
 }
 
par(oma=c(1,1,1,1))
for(i in 1:5){
plot(density(pcs$x[,i]),main="",xlab=paste("PC",i),col=pcols[4])
lines(density(pcs$x[grep(confirmgenes,dat$gene,ignore.case = T),i]),col=pcols[1])
lines(density(pcs$x[grep(risks,dat$gene,ignore.case = T),i]),col=pcols[2])
#lines(density(pcs$x[which(dat$gene%in%denovo_mouse),i]),col="orange")
lines(density(pcs$x[which(dat$gene%in%negatives),i]),col=pcols[3])
legend("topright",legend = cats,col=pcols, pch=rep(20,4),bty='n')
}
i=1
for(i in c(0:3)){
  j=2
  plot(proj[,j],proj[,i+1],xlab=paste("PC",j,": ",format(pcs$sdev[j],trim = 2,digits = 2),sep=""),
       ylab=paste("PC",i+1,": ",format(pcs$sdev[i+1],trim = 2,digits = 2),sep=""),col=pcols[4],main="FACS-4tissues",
       pch=20,cex=1)
  
  points(pcs$x[,j],pcs$x[,i+1],xlab=paste("PC",j,sep=""),ylab=paste("PC",i+1,sep=""),col=pcols[3],main="FACS-4tissues",
       pch=20,cex=1)
  points(pcs$x[grep(confirmgenes,dat$gene,ignore.case = T),j],pcs$x[grep(confirmgenes,dat$gene,ignore.case = T),
                                i+1],
         col=pcols[1],
         pch=20,cex=1)
  points(pcs$x[grep(risks,dat$gene,ignore.case = T),j],pcs$x[grep(risks,dat$gene,ignore.case = T),i+1],
         col=pcols[2],
         pch=20,cex=1)
 # points(pcs$x[which(dat$gene%in%denovo_mouse),j],pcs$x[which(dat$gene%in%denovo_mouse),
  #                                                                  i+1],
 #        col="orange",
  #       pch=20,cex=1)
  
  # text(pcs$x[which(dat$gene%in%denovo_mouse),j],
  #      pcs$x[which(dat$gene%in%denovo_mouse),i+1],col="orange",
  #      dat$gene[which(dat$gene%in%denovo_mouse)],cex = 1,pos=1)
  
  text(pcs$x[grep(confirmgenes,dat$gene,ignore.case = T),j],
       pcs$x[grep(confirmgenes,dat$gene,ignore.case = T),i+1],col=pcols[1],
       dat$gene[grep(confirmgenes,dat$gene,ignore.case = T)],cex = 1,pos=1)

  text(pcs$x[grep(risks,dat$gene,ignore.case = T),j],pcs$x[grep(risks,dat$gene,ignore.case = T),i+1],
       col=pcols[2],
       dat$gene[grep(risks,dat$gene,ignore.case = T)],cex = 1,pos=1)
#  legend("topright",legend = c("known","candidate","de novo","negative"),fill=c("red","blue","orange","black"),bty='n')
  legend("topright",legend =cats,col= pcols,pch=20,bty='n')
}
dev.off()
library("elasticnet")
y<-rep(0,dim(pcs$x)[1])
y[grep(confirmgenes,dat$gene,ignore.case = T)]<-1
ob<-enet(x=pcs$x,y=y,lambda = 0.5)
pdf("~/Dropbox (CGC)/PAH/US_UK/result/geneExpression/Enet_no_jpt2_col6A5.pdf")
plot(ob,main="PCs",las=2)
write.csv(attr(ob$beta.pure,"scaled:scale")*ob$beta.pure,file="~/Dropbox (CGC)/PAH/US_UK/result/geneExpression/Enet_pcs.result.csv")

fits<-predict.enet(ob,rbind(pcs$x,proj),mode="step",s=10,type="fit")

ob<-enet(x=as.matrix(dat[,indexs]),y=y,lambda = 0.5)
plot(ob,main="rawdata",las=2)
write.csv(attr(ob$beta.pure,"scaled:scale")*ob$beta.pure,file="~/Dropbox (CGC)/PAH/US_UK/result/geneExpression/Enet_rawdata.result.csv")
fits<-predict.enet(ob,as.matrix(dat_all[,indexs]),mode="step",s=10,type="fit")


allgenes<-c(dat$gene,dat_project$gene)
outfit<- data.frame(gene=allgenes,predictor=fits$fit,stringsAsFactors = F)
outfit$known=0;
outfit[grep(confirmgenes,outfit$gene,ignore.case = T),"known"]<-"known"

outfit[grep(risks,outfit$gene,ignore.case = T),"known"]<-"candidate"
outfit[which(outfit$gene%in%denovo_mouse),"known"]<-"denovo"
outfit<-outfit[order(outfit$predictor),]
outfit$rank<-1:dim(outfit)[1]
write.csv(outfit,file="~/Dropbox (CGC)/PAH/US_UK/result/geneExpression/predict_genes.csv")
plot(fits$fit)
plot(density(fits$fit))
keyg="fbln2|bmpr2|sox17|^kdr$|pdgfd|gdf2|tbx4|acvrl1|col6a5|hn1l|smad9|notch1"
plot(density(fits$fit))
abline(v=fits$fit[grep(keyg,allgenes,ignore.case = T)],col="gray",lty=2)
text(fits$fit[grep(keyg,allgenes,ignore.case = T)],seq(1,50,by = 5),allgenes[grep(keyg,allgenes,ignore.case = T)],col="gray")
dev.off()



test<-function(){
  
  pdf("/Users/nazhu/Dropbox (CGC)/PAH/US_UK/result/geneExpression/TabulaMutis.pdf") #("~/terra/Resources/GeneExpression/single_cell/TabulaMuris/TabulaMuris_PCA.pdf")
  knowngenes<-"SOX17|BMPR2|TBX4|KCNK3|ATP13A3|GDF2|AQP1$|KLK1$|GGCX|ACVRL1|KDR|^ENG$|EIF2AK4|SMAD4|^CAV1$|NOTCH1$|SMAD4|SMAD9|Fbln2"
  risks<-"col6a5|fbln2|pdgfd|chd9|hn1l|jpt2" #|mrrf|ccdc28A|vps51|mmp10|chd9|CHKB|myf5|irak4|npffr2|myo5b|ipo4|hn1l|jpt2"
  negatives<-
    dat1<-read.csv("~/terra/Resources/GeneExpression/single_cell/TabulaMuris/Droplet_tissues_genes.summary.csv",
                   header=1,stringsAsFactors = F,comment.char = "",check.names = F)
  dat1<-dat1[which(dat1$gene%in%allg$V4),]
  ss<-unlist(lapply(1:dim(dat1)[1],FUN=function(i){sum(dat1[i,2:dim(dat1)[2]])}))
  dat1<-dat1[which(ss>0),]
  indexs<-2:dim(dat1)[2]
  indexs<-grep("heart|lung|Trachea|aorta",names(dat1),ignore.case = T)
  pcs<-prcomp(as.matrix(dat1[,indexs]),scale. = T,center = T)
  n=8
  km<-kmeans(pcs$x,n)
  cols<-colors()
  cols<-sample(cols,size = n,replace = F)
  
  for(i in c(1:3)){
    j=1;
    plot(pcs$x[,j],pcs$x[,i+1],xlab=paste("PC",j,sep=""),ylab=paste("PC",i+1,sep=""),col=cols[km$cluster],main="droplet-4tissues",pch=20,cex=0.5)
    points(pcs$x[grep(knowngenes,dat1$gene,ignore.case = T),j],pcs$x[grep(knowngenes,dat1$gene,ignore.case = T),i+1],
           col=cols[km$cluster[grep(knowngenes,dat1$gene,ignore.case = T)]],
           pch=20,cex=3)
    points(pcs$x[grep(knowngenes,dat1$gene,ignore.case = T),j],pcs$x[grep(knowngenes,dat1$gene,ignore.case = T),i+1],col="black",
           pch=1,cex=3)
    # points(pcs$x[grep("Fbln2",dat$gene,ignore.case = T),j],pcs$x[grep("Fbln2",dat$gene,ignore.case = T),i+1],
    #        col=cols[km$cluster[grep("Fbln2",dat$gene,ignore.case = T)]],pch=20,cex=3)
    text(pcs$x[grep(knowngenes,dat1$gene,ignore.case = T),j],pcs$x[grep(knowngenes,dat1$gene,ignore.case = T),i+1],col="blue",
         dat1$gene[grep(knowngenes,dat1$gene,ignore.case = T)],cex = 0.5)
    points(pcs$x[grep(risks,dat1$gene,ignore.case = T),j],pcs$x[grep(risks,dat1$gene,ignore.case = T),i+1],col="black",
           pch=1,cex=3)
    text(pcs$x[grep(risks,dat1$gene,ignore.case = T),j],pcs$x[grep(risks,dat1$gene,ignore.case = T),i+1],col="darkgreen",
         dat1$gene[grep(risks,dat1$gene,ignore.case = T)],cex = 0.5)
  }
  ## K-means
  
  dev.off()
  
  dat<-read.csv("~/terra/Resources/GeneExpression/single_cell/TabulaMuris/FACS_tissues_genes.summary.csv",
                header=1,stringsAsFactors = F,comment.char = "",check.names = F)
  dat<-dat[which(dat$gene%in%allg$V4),]
  ss<-unlist(lapply(1:dim(dat)[1],FUN=function(i){sum(dat[i,2:dim(dat)[2]])}))
  dat<-dat[which(ss>0),]
  indexs<-2:dim(dat)[2]
  indexs<-intersect(grep("heart|lung|aorta",names(dat),ignore.case = T),grep("strom|endothelial",names(dat),ignore.case = T))
  pcs<-prcomp(as.matrix(dat[,indexs]),scale. = T,center = T)
  n=6
  km<-kmeans(pcs$x,n)
  xx<-pcs$x[grep(confirmgenes,dat$gene,ignore.case = T),c(26,30)]
  yy<-pcs$x[,c(26,30)]
  distmat<-as.matrix(pdist(xx,yy,indices.A = 1:dim(xx)[1],indices.B = 1:dim(yy)[1]))
  dscore<-unlist(lapply(1:dim(distmat)[2],FUN = function(i){min(distmat[,i]) }))
  table(km$cluster[grep(knowngenes,dat$gene,ignore.case = T)])
  
  for(c in 1:n){
    print(c)
    print(dat$gene[intersect(which(km$cluster==c),grep(knowngenes,dat$gene,ignore.case = T))])
  }
  table(km$cluster)
  
  cols<-colors()
  cols<-sample(cols,size = n,replace = F)
  for(i in c(1:3)){
    j=1;
    plot(pcs$x[,j],pcs$x[,i+1],xlab=paste("PC",j,sep=""),ylab=paste("PC",i+1,sep=""),col=cols[km$cluster],main="FACS-4tissues",
         pch=20,cex=0.5)
    points(pcs$x[grep(knowngenes,dat$gene,ignore.case = T),j],
           pcs$x[grep(knowngenes,dat$gene,ignore.case = T),i+1],
           col=cols[km$cluster[grep(knowngenes,dat$gene,ignore.case = T)]],
           pch=20,cex=3)
    
    points(pcs$x[grep(risks,dat$gene,ignore.case = T),j],
           pcs$x[grep(risks,dat$gene,ignore.case = T),i+1],
           col=cols[km$cluster[grep(risks,dat$gene,ignore.case = T)]],
           pch=20,cex=3)
    
    #points(pcs$x[grep("Fbln2",dat$gene,ignore.case = T),j],pcs$x[grep("Fbln2",dat$gene,ignore.case = T),i+1],
    #       col=cols[km$cluster[grep("Fbln2",dat$gene,ignore.case = T)]],pch=20,cex=3)
    # points(pcs$x[grep(knowngenes,dat$gene,ignore.case = T),j],pcs$x[grep(knowngenes,dat$gene,ignore.case = T),i+1],col="black",
    #       pch=1,cex=3)
    text(pcs$x[grep(knowngenes,dat$gene,ignore.case = T),j],pcs$x[grep(knowngenes,dat$gene,ignore.case = T),i+1],col="green",
         dat$gene[grep(knowngenes,dat$gene,ignore.case = T)],cex = 1)
    # points(pcs$x[grep(risks,dat$gene,ignore.case = T),j],pcs$x[grep(risks,dat$gene,ignore.case = T),i+1],col="black",
    #         pch=1,cex=3)
    text(pcs$x[grep(risks,dat$gene,ignore.case = T),j],pcs$x[grep(risks,dat$gene,ignore.case = T),i+1],col="green",
         dat$gene[grep(risks,dat$gene,ignore.case = T)],cex = 1)
    points(pcs$x[,j],pcs$x[,i+1],xlab=paste("PC",j,sep=""),ylab=paste("PC",i+1,sep=""),col=cols[as.integer(dscore/0.5)],main="FACS-4tissues",pch=20,cex=0.5)
  }
  
  dev.off()
  
  n=4
  km<-kmeans(xx,n)
  
  table(km$cluster[grep(knowngenes,testgenes,ignore.case = T)])
  
  for(c in 1:n){
    print(c)
    print(testgenes[intersect(which(km$cluster==c),grep(knowngenes,testgenes,ignore.case = T))])
  }
  table(km$cluster)
  
  
  
  
  for(i in c(1:3)){
    j=1;
    plot(xx[,j],xx[,i+1],xlab=paste("PC",j,sep=""),ylab=paste("PC",i+1,sep=""),col=cols[km$cluster],main="FACS-4tissues",
         pch=20,cex=0.5)
    # points(proj[,j],proj[,i+1],xlab=paste("PC",j,sep=""),ylab=paste("PC",i+1,sep=""),col="gray",
    #      pch=20,cex=0.5)
    points(xx[grep(knowngenes,testgenes,ignore.case = T),j],
           xx[grep(knowngenes,testgenes,ignore.case = T),i+1],
           col=cols[km$cluster[grep(knowngenes,testgenes,ignore.case = T)]],
           pch=20,cex=5)
    
    points(xx[grep(risks,testgenes,ignore.case = T),j],
           xx[grep(risks,testgenes,ignore.case = T),i+1],
           col=cols[km$cluster[grep(risks,testgenes,ignore.case = T)]],
           pch=20,cex=5)
    
    #points(pcs$x[grep("Fbln2",dat$gene,ignore.case = T),j],pcs$x[grep("Fbln2",dat$gene,ignore.case = T),i+1],
    #       col=cols[km$cluster[grep("Fbln2",dat$gene,ignore.case = T)]],pch=20,cex=3)
    #points(pcs$x[grep(knowngenes,dat$gene,ignore.case = T),j],pcs$x[grep(knowngenes,dat$gene,ignore.case = T),i+1],col="black",
    #        pch=1,cex=3)
    text(xx[grep(knowngenes,testgenes,ignore.case = T),j],xx[grep(knowngenes,testgenes,ignore.case = T),i+1],col="blue",
         testgenes[grep(knowngenes,testgenes,ignore.case = T)],cex = 1)
    #points(pcs$x[grep(risks,dat$gene,ignore.case = T),j],pcs$x[grep(risks,dat$gene,ignore.case = T),i+1],col="black",
    #       pch=1,cex=3)
    text(xx[grep(risks,testgenes,ignore.case = T),j],xx[grep(risks,testgenes,ignore.case = T),i+1],col="darkgreen",
         testgenes[grep(risks,testgenes,ignore.case = T)],cex = 1)
    # legend("topright",legend = unique(km$cluster),fill=cols[unique(km$cluster)])
    
  }
  
  dev.off()
  
  write.csv(dat_t,file="~/Dropbox (CGC)/PAH/US_UK/result/geneExpression/located_relevantcluster.txt",row.names = F)
  dat1<-read.csv("~/Dropbox (CGC)/PAH/US_UK/result/geneExpression/",
                 header=1,stringsAsFactors = F,comment.char = "",check.names = F)
  dat1<-dat1[unique(c(which(dat1$gene%in%c(negatives)),grep(risks,dat1$gene,ignore.case = T),grep(knowngenes,dat1$gene,ignore.case = T))),]
  ss<-unlist(lapply(1:dim(dat1)[1],FUN=function(i){sum(dat1[i,2:dim(dat1)[2]])}))
  dat1<-dat1[which(ss>0),]
  indexs<-2:dim(dat1)[2]
  #indexs<-grep("heart|lung|Trachea|aorta",names(dat1),ignore.case = T)
  pcs<-prcomp(as.matrix(dat1[,indexs]),scale. = T,center = T)
  n=8
  km<-kmeans(pcs$x,n)
  cols<-colors()
  cols<-sample(cols,size = n,replace = F)
  
  for(i in c(1:3)){
    j=1;
    plot(pcs$x[,j],pcs$x[,i+1],xlab=paste("PC",j,sep=""),ylab=paste("PC",i+1,sep=""),col=cols[km$cluster],main="droplet-4tissues",pch=20,cex=0.5)
    points(pcs$x[grep(knowngenes,dat1$gene,ignore.case = T),j],pcs$x[grep(knowngenes,dat1$gene,ignore.case = T),i+1],
           col=cols[km$cluster[grep(knowngenes,dat1$gene,ignore.case = T)]],
           pch=20,cex=3)
    points(pcs$x[grep(knowngenes,dat1$gene,ignore.case = T),j],pcs$x[grep(knowngenes,dat1$gene,ignore.case = T),i+1],col="black",
           pch=1,cex=3)
    # points(pcs$x[grep("Fbln2",dat$gene,ignore.case = T),j],pcs$x[grep("Fbln2",dat$gene,ignore.case = T),i+1],
    #        col=cols[km$cluster[grep("Fbln2",dat$gene,ignore.case = T)]],pch=20,cex=3)
    text(pcs$x[grep(knowngenes,dat1$gene,ignore.case = T),j],pcs$x[grep(knowngenes,dat1$gene,ignore.case = T),i+1],col="blue",
         dat1$gene[grep(knowngenes,dat1$gene,ignore.case = T)],cex = 0.5)
    points(pcs$x[grep(risks,dat1$gene,ignore.case = T),j],pcs$x[grep(risks,dat1$gene,ignore.case = T),i+1],col="black",
           pch=1,cex=3)
    text(pcs$x[grep(risks,dat1$gene,ignore.case = T),j],pcs$x[grep(risks,dat1$gene,ignore.case = T),i+1],col="darkgreen",
         dat1$gene[grep(risks,dat1$gene,ignore.case = T)],cex = 0.5)
  }
  ## K-means
  
  
  
}
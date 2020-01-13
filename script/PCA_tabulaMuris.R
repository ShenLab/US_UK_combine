setwd("yourfolder/")
## PCA 

cats<-c("known risk genes","candidate risk genes","random genes", "other genes")
pcols<-c("red","dodgerblue",rgb(0,0,0,0.4),rgb(220/255,220/255,220/255,0.5)) #"sienna2"
allg<-read.table("data/MouseHumanIDmapper.txt",
                 header = 1,comment.char = "",sep="\t",stringsAsFactors = F)

knowngenes<-"SOX17|BMPR2|TBX4|KCNK3|ATP13A3|GDF2|AQP1$|KLK1$|GGCX|ACVRL1|^KDR$|^ENG$|EIF2AK4|SMAD4|^CAV1$|NOTCH1$|SMAD4|SMAD9"
confirmgenes<-"SOX17|BMPR2|TBX4|KCNK3|ACVRL1|^ENG$|EIF2AK4|^CAV1$|NOTCH1$|SMAD4|SMAD9"

risks<-"fbln2|pdgfd|^kdr$" 

all_negatives<-allg[intersect(which(!allg$HumanGene%in%c(knowngenes,risks)==T),
                              grep(paste(knowngenes,risks,sep="|"),allg$HumanGene,ignore.case = T,invert = T)),"MouseGene"]



dat_all<-read.csv("data/FACS_tissues_genes.summary.csv",
              header=1,stringsAsFactors = F,comment.char = "",check.names = F)
dat_all<-dat_all[which(dat_all$gene%in%allg$MouseGene),]

all_negatives<-all_negatives[which(all_negatives%in%dat_all$gene)]
negatives<-sample(all_negatives,100,replace = F)
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

matri<-as.matrix(dat[,indexs])
row.names(matri)<-dat$gene
colnames(matri)<-cells
matri<-matri[c(grep(confirmgenes,dat$gene,ignore.case = T),
           grep(risks,dat$gene,ignore.case = T),which(dat$gene%in%c(negatives))),]

pdf("result/Heatmap.pdf",width=10,height=10)
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
write.csv(mat,"result/PC2_pergene.csv")
#cols<-sample(colors(),size = 15,replace = F)
cex.tex=1.3
pdf("result/PCA_projection.pdf",width=7,height=7)

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

pdf("result/PCA_loadings.pdf")
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

par(mar=c(4,4,1,1))
hist(mat$rank[c(grep(confirmgenes,mat$gene,ignore.case = T))],col=rgb(1,0,0,alpha = 0.5),breaks= seq(0,14000,1000),
     freq = T,main="",xlab="",ylab="",xaxt='n',yaxt='n')
hist(mat$rank[c(grep(risks,mat$gene,ignore.case = T))],breaks =  seq(0,14000,1000),
     col=rgb(0,206/256,230/256,alpha=0.5),add=T,freq=T)
axis(1,at = seq(0,14000,1000),labels = paste(c(0:14),"k",sep=""),cex.axis=1.5,las=2)
axis(2,at = seq(0,6,1),labels = seq(0,6,1),cex.axis=1.5,las=2)
title(xlab="", ylab="Frequency",cex.lab=1.5,line=2.5)
title(xlab="PC2 rank", cex.lab=1.5,line=3)
legend("topright",legend = cats[1:2],fill=c(rgb(1,0,0,alpha = 0.5),rgb(0,206/256,230/256,alpha=0.5)),cex=1.5,bty='n')


dev.off()


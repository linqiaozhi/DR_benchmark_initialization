####################
## Computing embeddings for Samusik_01
####################

library(RColorBrewer)
library(rsvd)
source("../utils.R")

## Loading FItSNE - Has to be installed separetely first, see https://github.com/KlugerLab/FIt-SNE
source("../FIt-SNE/fast_tsne.R")

## Loading and subsetting the data
load("./data/Samusik_01/dataset.RData")
xp=d$xp_src
regular_channels=subset(read.csv("./annots/template_Samusik_01.csv",header=TRUE,stringsAsFactors=FALSE),Select=="y")[,2]

## Saving transformed and column-subsetted dataset for later benchmarking
xp=asinh(xp[,regular_channels])
save(xp,file="./rdata/Samusik_01/xp.RData")

load("./rdata/Samusik_01/xp.RData")

rerun=TRUE ## Set to TRUE to re-run UMAP and t-SNE
if(rerun){
    ## Running UMAP on Samusik_01
    library(reticulate)
    use_python("/data/george/miniconda3/bin/python",TRUE)
    umap_module=import("umap",convert=TRUE)
    umap=umap_module$UMAP(n_neighbors=15L,min_dist=0.2,metric="euclidean",verbose=TRUE,random_state=123L)$fit_transform(xp[,regular_channels])
    colnames(umap)=c("UMAP1","UMAP2")
    save(umap,file="./rdata/Samusik_01/umap.RData")


    umap_random=umap_module$UMAP(n_neighbors=15L,min_dist=0.2,metric="euclidean",verbose=TRUE,random_state=123L, init='random')$fit_transform(xp[,regular_channels])
    colnames(umap_random)=c("UMAP1","UMAP2")
    save(umap_random,file="./rdata/Samusik_01/umap_random.RData")

    set.seed(123)
    xp_c <- scale(xp[,regular_channels], center=T, scale=F)
    rsvdout <- rsvd(xp_c, k=2, q=10)
    xp_top2pcs <- rsvdout$u %*% diag(rsvdout$d)
    xp_top2pcs_sc <- 0.0001*(xp_top2pcs/sd(xp_top2pcs[,1]))
    fitsnefixed=fftRtsne(xp[,regular_channels],  initialization=xp_top2pcs_sc, rand_seed=123L)
    colnames(fitsnefixed)=c("FItSNE1","FItSNE2")
    save(fitsnefixed,file="./rdata/Samusik_01/fitsnefixed.RData")

    xp_c <- scale(xp[,regular_channels], center=T, scale=F)
    rsvdout <- rsvd(xp_c, k=2,q=10)
    xp_top50pcs <- rsvdout$u %*% diag(rsvdout$d)
    fitsnerandom=fftRtsne(xp_c,   rand_seed=123L)
    colnames(fitsnerandom)=c("FItSNE1","FItSNE2")
    save(fitsnerandom,file="./rdata/Samusik_01/fitsnerandom.RData")

    ### Running t-SNE on Samusik_01
    #set.seed(123)
    #library(Rtsne)
    #tsne=Rtsne(xp[,regular_channels],max_iter=1000,verbose=TRUE)$Y
    #colnames(tsne)=c("tSNE1","tSNE2")
    #save(tsne,file="./rdata/Samusik_01/tSNE.RData")
}

manual_labels=d$xp_src[,"label"]
manual_labels[is.na(manual_labels)]=0
gs=read.table(file="./annots/population_names_Samusik_and_gating.csv",sep=",",header=TRUE,stringsAsFactors=FALSE,check.names=FALSE)
codePops=c("0"="Ungated",setNames(gs[,2],gs[,1]))
populations=manual_labels
save(populations,file="./rdata/Samusik_01/populations_manual.RData")

####################
## Fig 2AB
####################

#load("./rdata/Samusik_01/xp.RData")
#load("./rdata/Samusik_01/umap.RData")
#load("./rdata/Samusik_01/tSNE.RData")
#load("./rdata/Samusik_01/populations_manual.RData")

coords=aggregate_df(umap_random,populations,median,1)
png("./graphs/Samusik_random.png",width=4000,height=2000,res=300)
par(mfrow=c(1,2),mar=c(4,4,0,0))
plot(umap_random,pch=16,cex=0.5,col=c("0"=unname(makeTransparent("gray",100)),setNames(colorRampPalette(brewer.pal(12,"Set3"))(24),1:24))[as.character(populations)],bty="n")
for(i in 1:24){
    text(
        x=coords[as.character(i),1],
        y=coords[as.character(i),2],
        labels=gs[as.character(i),"name"],
        font=2,
        cex=0.75,
        xlab="UMAP_random1",
        ylab="UMAP_random2"
    )
}
coords_fitsnerandom=aggregate_df(fitsnerandom,populations,median,1)
plot(fitsnerandom,pch=16,cex=0.5,col=c("0"=unname(makeTransparent("gray",100)),setNames(colorRampPalette(brewer.pal(12,"Set3"))(24),1:24))[as.character(populations)],bty="n")
for(i in 1:24){
    text(
        x=coords_fitsnerandom[as.character(i),1],
        y=coords_fitsnerandom[as.character(i),2],
        labels=gs[as.character(i),"name"],
        font=2,
        cex=0.75,
        xlab="t-SNE_random1",
        ylab="t-SNE_random2",
        xpd=NA
    )
}
dev.off()



coords=aggregate_df(umap,populations,median,1)
png("./graphs/Samusik_fixed.png",width=4000,height=2000,res=300)
par(mfrow=c(1,2),mar=c(4,4,0,0))
plot(umap,pch=16,cex=0.5,col=c("0"=unname(makeTransparent("gray",100)),setNames(colorRampPalette(brewer.pal(12,"Set3"))(24),1:24))[as.character(populations)],bty="n")
for(i in 1:24){
    text(
        x=coords[as.character(i),1],
        y=coords[as.character(i),2],
        labels=gs[as.character(i),"name"],
        font=2,
        cex=0.75,
        xlab="UMAP_fixed1",
        ylab="UMAP_fixed2"
    )
}
coords_fitsnefixed=aggregate_df(fitsnefixed,populations,median,1)
plot(fitsnefixed,pch=16,cex=0.5,col=c("0"=unname(makeTransparent("gray",100)),setNames(colorRampPalette(brewer.pal(12,"Set3"))(24),1:24))[as.character(populations)],bty="n")
for(i in 1:24){
    text(
        x=coords_fitsnefixed[as.character(i),1],
        y=coords_fitsnefixed[as.character(i),2],
        labels=gs[as.character(i),"name"],
        font=2,
        cex=0.75,
        xlab="t-SNE_fixed1",
        ylab="t-SNE_fixed2",
        xpd=NA
    )
}
dev.off()

#####################
### Density heatmap
#####################
#
#umap_counts=log10(1+table(cut(umap[,1],breaks=300),cut(umap[,2],breaks=300)))
#tsne_counts=log10(1+table(cut(tsne[,1],breaks=300),cut(tsne[,2],breaks=300)))
#
#png("./graphs/Density heatmap.png",height=1000,width=2000,res=300)
#par(mfrow=c(1,2))
#par(mar=c(1,1,3,3),bty="l")
#image(umap_counts,col=c("white",colorRampPalette(rev(brewer.pal(9,"RdYlBu")))(100)),breaks=seq(0,max(umap_counts),length.out=102),xaxt="n",yaxt="n")
#title(main="UMAP")
#x=seq((par("usr")[1]+4*par("usr")[2])/5,par("usr")[2],length.out=51)
#yb=rep((par("usr")[4]-par("usr")[3])*0.85+par("usr")[3],50)
#yt=rep((par("usr")[4]-par("usr")[3])*0.9,+par("usr")[3],50)
#rect(x[-length(x)],yb,x[-1],yt,col=c("white",colorRampPalette(rev(brewer.pal(9,"RdYlBu")))(49)),border=c("white",colorRampPalette(rev(brewer.pal(9,"RdYlBu")))(49)))
#text(x=mean(c(x[1],tail(x,1))),y=yt[1],pos=3,labels="Number of events\nper surface unit",xpd=NA,cex=0.5)
#text(x=x[1],y=yb[1],pos=1,labels="0",xpd=NA,cex=0.4)
#text(x=tail(x,1),y=yb[1],pos=1,labels="max(UMAP)",xpd=NA,cex=0.4)
#image(tsne_counts,col=c("white",colorRampPalette(rev(brewer.pal(9,"RdYlBu")))(100)),breaks=seq(0,max(tsne_counts),length.out=102),xaxt="n",yaxt="n")
#title(main="tSNE")
#rect(x[-length(x)],yb,x[-1],yt,col=c("white",colorRampPalette(rev(brewer.pal(9,"RdYlBu")))(49)),border=c("white",colorRampPalette(rev(brewer.pal(9,"RdYlBu")))(49)))
#text(x=mean(c(x[1],tail(x,1))),y=yt[1],pos=3,labels="Number of events\nper surface unit",xpd=NA,cex=0.5)
#text(x=x[1],y=yb[1],pos=1,labels="0",xpd=NA,cex=0.4)
#text(x=tail(x,1),y=yb[1],pos=1,labels="max(tSNE)",xpd=NA,cex=0.4)
#dev.off()
#
#####################
### Sup-erythrocytes
#####################
#
#library(matlab)
#png("./graphs/Erythrocytes.png",height=2000,width=2000,res=300)
#col=toColors_continuous(xp[,"Ter119"],palette=colorRampPalette(c("#0000AA","#0000FF","#0055FF","#00AAFF","#00FFFF","#55FFAA","#AAFF55","#FFFF00","#FFAA00","#FF5500"))(100))
#plot(umap,col=col,pch=16,cex=0.2,bty="l")
#x=seq((par("usr")[1]+4*par("usr")[2])/5,par("usr")[2],length.out=51)
#yb=rep((par("usr")[4]-par("usr")[3])*0.95+par("usr")[3],50)
#yt=rep((par("usr")[4]-par("usr")[3])*1,+par("usr")[3],50)
#rect(x[-length(x)],yb,x[-1],yt,col=jet.colors(50),border=colorRampPalette(c("#0000AA","#0000FF","#0055FF","#00AAFF","#00FFFF","#55FFAA","#AAFF55","#FFFF00","#FFAA00","#FF5500"))(50))
#text(x=mean(c(x[1],tail(x,1))),y=yb[1],pos=1,labels="Ter119")
#dev.off()

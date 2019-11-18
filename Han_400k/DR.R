## Loading the UMAP module
source("../utils.R")

## Loading FItSNE
source("../FIt-SNE/fast_tsne.R")

## Loading t-SNE
library(Rtsne)

####################
## Dimensionality reduction for the full Han dataset (400k scRNAseq)
####################

load("./rdata/pca.RData") ## Computed in ./data_parsing_full_dataset.R

library(reticulate)
use_python("/data/george/miniconda3/bin/python",TRUE)
umap_module=import("umap",convert=TRUE)
## UMAP
umap=umap_module$UMAP(n_neighbors=30L,min_dist=0.2,metric="euclidean",verbose=FALSE,random_state=123L,verbose=TRUE)$fit_transform(pca)
colnames(umap)=c("UMAP1","UMAP2")
save(umap,file="./rdata/umap.RData")

## UMAP
umaprandom=umap_module$UMAP(n_neighbors=30L,init='random',min_dist=0.2,metric="euclidean",verbose=FALSE,random_state=123L,verbose=TRUE)$fit_transform(pca)
colnames(umaprandom)=c("UMAP1","UMAP2")
save(umaprandom,file="./rdata/umaprandom.RData")

### tSNE
#tsne=Rtsne(pca)$Y
#colnames(tsne)=c("tSNE1","tSNE2")
#save(tsne,file="./rdata/tsne.RData")

## FItSNE
fitsne=fftRtsne(pca,rand_seed=123L)
colnames(fitsne)=c("FItSNE1","FItSNE2")
save(fitsne,file="./rdata/fitsne.RData")


pca_sc <- 0.0001*(pca/sd(pca[,1]))
fitsnefixed=fftRtsne(pca,  initialization=pca_sc[,1:2], rand_seed=123L)
colnames(fitsnefixed)=c("FItSNE1","FItSNE2")
save(fitsnefixed,file="./rdata/fitsnefixed.RData")


### FItSNE-le
#fitsnele=fftRtsne(pca,rand_seed=123L,start_late_exag_iter=800L,late_exag_coeff=4)
#colnames(fitsnele)=c("FItSNE_le1","FItSNE_le2")
#save(fitsnele,file="./rdata/fitsnele.RData")

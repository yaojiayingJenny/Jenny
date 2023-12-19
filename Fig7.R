######Fig7 plot--------------------------------------------------------
#ptest
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193816
# we download GSE193816 human asthma h5ad data and convert to rds file


```
import pandas as pd
import scanpy as sc
import numpy as np
import anndata
anndata.__version__

name='GSE193816_mnp_data'
adata=anndata.read(name+'.h5ad')
metadata=pd.DataFrame(data=adata.obs)
metadata.to_csv(name + '_metadata.tsv',sep=",") 
```
	
library(SeuratDisk)
library(Seurat)
library(rhdf5) #h5read
library(zellkonverter) #readH5AD

Convert(paste0(name,".h5ad"), dest = "h5seurat", overwrite = F)
sc <- LoadH5Seurat(paste0(name,".h5seurat"),meta.data = F,misc=F) #先不添加meta.data(如何报错的话)
meta<-read.csv(paste0(name,'_metadata.tsv'),header=T,row.names=1)
sc@meta.data<-meta
mnp<-sc #GSE193816_mnp_data.rds



#####-----------------------------------------------------------
#FindTransferAnchors predict celltypes with human MNP dataset 
library(scibet)
library(Seurat)
suppressMessages(library(ggplot2))
suppressMessages(library(tidyverse))
suppressMessages(library(scibet))
suppressMessages(library(viridis)) #scale_fill_viridis
suppressMessages(library(ggsci))
library(Seurat)
library(readr)
get_confusion_plot <- function(actual, predicted, align = 1){
  require(reshape2)
  require(ggplot2)

  if (class(actual) == "factor"){ actual_levels <- levels(actual) }
  else { actual_levels <- unique(actual) }
  
  if (class(predicted) == "factor"){ predicted_levels <- levels(predicted) }
  else { predicted_levels <- unique(predicted) }
  
  get_confusion_prop_tbl <- function(actual, predicted, align = 1){
    if (!(align %in% c(1,2))) {stop("align must be 1 or 2")}
    conf <- table(actual, predicted)
    if (align == 2) {
      conf <- t(conf)
      return (t(conf / rowSums(conf)))
    }
    return (conf / rowSums(conf))
  }
  
  conf_prop <- get_confusion_prop_tbl(actual, predicted, align)
  conf_prop_melt <- reshape2::melt(conf_prop, variable.name= predicted, id = actual)

  conf_prop_melt[,"actual"] = factor(conf_prop_melt$actual, levels=actual_levels)
  conf_prop_melt[,"predicted"] = factor(conf_prop_melt$predicted, levels=predicted_levels)
  
  p0<-theme(legend.text = element_text(face="plain", color="black",size=10),
        axis.text.x = element_text(color="black",size=10,angle=90),
        axis.text.y = element_text(color="black",size=10),
        axis.title.x = element_text(face="plain", color="black",size=15),
        axis.title.y = element_text(face="plain", color="black",size=15))
  plt <- ggplot2::ggplot(data = conf_prop_melt, aes(x=predicted,y=actual)) +
    geom_tile(aes(fill=value,colour='black')) +theme_bw()+p0+scale_fill_gradientn(colours=c('#46085B','#453480','#20918C','#5AC763','#DFE318'), guide=guide_colorbar(reverse=F))
  
  #print(plt)
  return(plt)
}


RenameGenesSeurat <- function(obj = ls.Seurat[[i]], newnames = tmp) { # Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.
  print("Run this before integration. It only changes obj@assays$RNA@counts and @data ")
  RNA <- obj@assays$RNA

  if (nrow(RNA) == length(newnames)) {
    if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- newnames
    if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- newnames
    # if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]]    <- newnames
  } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
  obj@assays$RNA <- RNA
  return(obj)
}

library(homologene)
mousegenes<-human2mouse(rownames(mnp),db = homologeneData2)
#rownames(mousegenes)<-mousegenes$mouseGene
mnp_human<-subset(mnp,features=mousegenes$humanGene)
can<-mnp_human@assays$RNA@counts@Dimnames[[1]]
map<-mousegenes[!duplicated(mousegenes$humanGene),]
rownames(map)<-map$humanGene
map<-map[can,]
mnp_mouse <- RenameGenesSeurat(mnp_human,newnames = map$mouseGene)

saveRDS(mnp_mouse,file='mnp_mouse.rds')


tmp1<-subset(mnp_mouse,cluster %in% c(seq(10:14)+9))
tmp2<-subset(mnp_mouse,cluster %in% seq(1:9))

tmp<-Mac;prefix<-'Mac_predicted';tmp0<-tmp2
tmp<-DC;prefix<-'DC_predicted';tmp0<-tmp1



anchors_ims<-20;anchors <- FindTransferAnchors(reference = tmp0, query = tmp, dims = 1:anchors_ims)	
predictions <- TransferData(anchorset = anchors, refdata = tmp0$celltype, dims = 1:anchors_ims)
pancreas.query0 <- AddMetaData(tmp, metadata = predictions)
pancreas.query<-pancreas.query0
table(pancreas.query$predicted.id)

library(scibet)
p<-get_confusion_plot(pancreas.query$annotation, pancreas.query$predicted.id)+geom_text(aes(label = round(value,2)))+coord_flip()+ggtitle(prefix)
pdf(paste0(prefix,'_FindTransferAnchors_Confusion_heatmap.pdf'),w=6,h=4)
print(p)
dev.off()

p1 <- DimPlot(pancreas.query, reduction = "umap", label = F,group.by = "annotation",label.size = 4,pt.size = 0.5,raster=FALSE)+ggtitle("raw labels")
p2 <- DimPlot(pancreas.query, reduction = "umap", label = F,group.by = "predicted.id",label.size = 4,pt.size = 0.5,raster=FALSE)+ggtitle("predicted labels")
p<-plot_grid(p1,p2,ncol = 2)
pdf(paste0(prefix,'_umap.pdf'),w=13,h=5)
print(p)
dev.off()



########---------
#Fig7E example of  Vlnplot  with wilcox.test
mnp<-readRDS('GSE193816_mnp_data.rds')
tmp<-mnp
genelist<-'CCL17,CXCL16,TREM2,GPNMB'
ptest1<-stat_compare_means(label = "p.signif",comparisons =list(c('AC Bln','AC Ag'),c('AA Bln','AA Ag')),method='wilcox.test',paired=FALSE,size=5,vjust = 0.5,hide.ns = TRUE)
dataset<-FetchData(tmp,vars='TGFB1')
meta<-tmp@meta.data[rownames(dataset),c('Channel','group','condition','Group','celltype')]
dataset<-cbind(dataset,meta)
dataset$Group<-factor(dataset$Group,levels=c('AC Bln','AC Ag','AA Bln','AA Ag'))
#genelist<-'CCL17,CXCL16,TREM2,GPNMB'
pair_groups4<-c('#69B38D','#339191','#E19056','#E6612C')
gene<-'CCL17';dataset$value<-dataset$TGFB1
p<-ggplot(dataset, aes(x = Group, y = value,fill=Group)) +geom_violin()+facet_wrap(~celltype,ncol=5)+ptest1+theme_bw()+p0+theme(legend.position='none')+xlab('')+scale_fill_manual(values=pair_groups4)+ylim(0,max(dataset$value+1))+ylab('Expression Level')+stat_summary (fun.data = function (x) data.frame (y=max(dataset$value), label = paste (round(mean (x),2))), geom="text",size=2)
pdf(paste0(gene,'_celltype_group_VlnPlot_ptest.pdf'),w=8,h=3)
print(p)
dev.off()
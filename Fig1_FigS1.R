library(Seurat)
library(pheatmap)
library(ggplot2)
library(CellChat)

######Fig1 plot--------------------------------------------------------
#26 celltypes
palette<-c(
    "Tcm CD4"="#F8961E",
    "Treg"="#7180B9",
    "Tgd"="#43aa8b",
    "Th2"="#F9C74F",
    "Tfh"="#90be6d",
    "Tc CD8"="#3423a6",
    "Tifn"="#9d4edd",
    "NKT"="#577590",
    "Tcm CD8"="#ff99ac",
    "Trm CD8"="#F94144", 
    "ILC1"="#B4F8C8",
    "ILC2"="#b68d40",
    "NK"="#A0E7E5",
    "B"="#6DA34D",
    "Neu"="#6ab8ee",
    "Basophils"="#f44336",
    "AM"="#fdd835",
    "IM"="#d84315",
    "Gpnmb+ Mac"="#478c5c",
    "CCR7+ DC"="#F433AB" ,
    "Dividing/Pclaf+ AMs"="#6f6f6f",
    "Plac8+ Mon"="#E3879E",
    "CD300e+ Mon"="#8C0000",
    "pDC"="#743481",
    "cDC1"="#f56038",
    "Mgl2+ DC"="#03378c")
	
rds<-readRDS('rdsfiles/IC.rds') 
ICs<-read.table('IC_celltypes.list',header=F,sep='\t')$V1
rds$celltype<-factor(rds$ann210320,levels=ICs)
rds$Group<-plyr::mapvalues(x =rds$group,from = c('PBS CD45+','HDM CD45+'),to = c('PBS','HDM'))

DefaultAssay(rds)<-'RNA'
Idents(rds)<-'celltype'

pdf('Fig1B.pdf',w=9,h=5)
p1 <- DimPlot(rds, reduction = "densMap", label = T,group.by = "celltype",label.size = 3,pt.size = 0,cols=palette,repel=T)+ggtitle('')
print(p1)
dev.off()

pdf('Fig1C.pdf',w=6.5,h=5)
p2 <- DimPlot(rds, reduction = "densMap", label = F,group.by = "group",label.size = 4,pt.size = 0,cols=c("PBS CD45+"="dark gray", "HDM CD45+"="dark blue"))+ggtitle('')
print(p2)
dev.off()


#marker genes
genelist<-'Cd79a,Cd19,S100a8,S100a9,Hdc,Itgax,Ccr7,Mgl2,Siglech,Irf8,Itgam,Gpnmb,Cd68,C1qa,Chil3,Plac8,Cd300e,Cd3e,Cd4,Cd8a,Il7r,Gata3,Foxp3,Id2,Ncr1,Gzma,Nkg7,Tcrg-C1'
genes<-unique(c(unlist(strsplit(genelist,split = ",",fixed=T))))
colors<-c("#6DA34D","#6ab8ee","#f44336","#743481","#f56038","#03378c","#F433AB","#E3879E","#8C0000","#fdd835","#d84315","#478c5c","#6f6f6f","#F8961E","#7180B9","#43aa8b","#F9C74F","#90be6d","#3423a6","#9d4edd","#577590","#ff99ac","#F94144","#B4F8C8","#b68d40","#A0E7E5")
p<-VlnPlot(rds, features = genes,split.by='celltype',pt.size = -1,slot="data",assay="RNA",stack = T,flip =T,cols=colors)+theme(legend.position='none',axis.text.y = element_text(color="black",size=10,angle=0),axis.text.x = element_text(color="black",size=10,angle=90))+xlab('')
pdf('Fig1D.pdf',w=7,h=8)
print(p)
dev.off()
##barplot of cell ratio
```
B
Neu
Basophils
pDC
cDC1
Mgl2+ DC
CCR7+ DC
Plac8+ Mon
CD300e+ Mon
AM
IM
Gpnmb+ Mac
Dividing/Pclaf+ AMs
Tcm CD4
Treg
Tgd
Th2
Tfh
Tc CD8
Tifn
NKT
Tcm CD8
Trm CD8
ILC1
ILC2
NK
```
per<-round(prop.table(table(rds@meta.data[,c('Group','celltype')]), margin=1),3)
per <-data.frame(per)
p0<-theme_bw()+theme(legend.position='right',panel.grid=element_blank(), legend.background = element_rect(colour = NA),
        legend.title=element_blank(),legend.text = element_text(face="plain", color="black",size=8),strip.text = element_text(size=20),
        axis.text.x = element_text(color="black",size=10,angle=0),
        axis.text.y = element_text(color="black",size=10),
        axis.title.x = element_text(face="plain", color="black",size=15),
        axis.title.y = element_text(face="plain", color="black",size=15))	
p2<-ggplot(data =per, aes(fill=celltype,y=Freq,x=Group)) +
geom_bar(stat = "identity",position = "stack",width=0.6)+p0+xlab('')+ylab('Fraction of total cells')+scale_fill_manual(values=palette)+guides(ncol=1)
pdf('Fig1E.pdf',w=6,h=6)
print(p2)
dev.off()

######FigS1 plot--------------------------------------------------------
genelist<-'Cd19,S100a8,Itgax,Itgam,Cd68,C1qa,Cd3e,Nkg7'
genes<-unique(c(unlist(strsplit(genelist,split = ",",fixed=T))))
pdf('FigS1A.pdf',w=16,h=6)
p<-FeaturePlot(rds,reduction = "densMap", features=genes,cols = c("lightgrey", "red"),label=F,pt.size=0,raster=F,label.size=3,ncol=4)
print(p)
dev.off()
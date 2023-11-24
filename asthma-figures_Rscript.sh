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


######Fig2 plot--------------------------------------------------------
DC<-readRDS('rdsfiles/DC.rds')

DC$celltype<-factor(DC$celltype,levels=c("CCR7+ DC","CD209e+ DC","cDC1","Dividing DC1.1","Dividing DC1.2","IFN response DC1","Mgl2+ DC","pDC","Plac8+ DC"))

Idents(DC)<-'celltype'
DefaultAssay(DC)<-'RNA'

#umap plot
#cluster_cols<-c('#D51F26','#208A42','#9983BD','#89288F','#89C75F','#3BBCA8','#D24B27','#8A9FD1','#FEE500','#272E6A','#0C727C','#6E4B9E','#C06CAB','#F47D2B','#E6C2DC','#D8A767','#90D5E4','#F37B7D','#7E1416')
stallion = c("1"="#D51F26","2"="#272E6A","3"="#208A42","4"="#89288F","5"="#F47D2B", "6"="#FEE500","7"="#8A9FD1","8"="#C06CAB","19"="#E6C2DC","10"="#90D5E4", "11"="#89C75F","12"="#F37B7D","13"="#9983BD","14"="#D24B27","15"="#3BBCA8", "16"="#6E4B9E","17"="#0C727C", "18"="#7E1416","9"="#D8A767","20"="#3D3D3D")			   
cluster_cols =as.vector(stallion)#ArchR::ArchRPalettes[['stallion']]
pdf('Fig2A.pdf',w=6.5,h=5)
p1 <- DimPlot(DC, reduction = "umap", label = F,group.by = "celltype",label.size = 4,pt.size = 1,cols=cluster_cols)+ggtitle('DC')
print(p1)
dev.off()


genelist<-'Ccr7,Ccl22,Cd209d,Cd209e,Itgae,Naaa,Mcm2,Mcm3,Mcm5,Dnmt1,Top2a,Mki67,Tubb5,Ifit1,Ifit2,Mgl2,Ear2,Cd209a,Siglech,Ly6c2,Ly6a,Cebpb,Fn1,Mafb,Csf1r'
genes<-unlist(strsplit(genelist,split = ",",fixed=T))

cluster.averages <- AverageExpression(object = DC, assays ='RNA',return.seurat = F)
data<-cluster.averages$RNA
da<-log2(data[genes,]+1)
filename<-'Fig2B.pdf'
color = colorRampPalette(c('blue','white','red'))(50)
pheatmap(da,filename=filename,color = color,width=5,scale='row',height=6,cluster_rows=F,cluster_cols=F,display_numbers = F,number_format = "%.0f",fontsize = 12,fontsize_col=8,border_color=NA,angle_col = "45")


# cell ratio 
per<-round(prop.table(table(DC@meta.data[,c('Group','celltype')]), margin=1),3)
per <-data.frame(per)
p0<-theme_bw()+theme(legend.position='right',panel.grid=element_blank(), legend.background = element_rect(colour = NA),
        legend.title=element_blank(),legend.text = element_text(face="plain", color="black",size=8),strip.text = element_text(size=20),
        axis.text.x = element_text(color="black",size=10,angle=0),
        axis.text.y = element_text(color="black",size=10),
        axis.title.x = element_text(face="plain", color="black",size=15),
        axis.title.y = element_text(face="plain", color="black",size=15))	
p2<-ggplot(data =per, aes(fill=celltype,y=Freq,x=Group)) +
geom_bar(stat = "identity",position = "stack",width=0.6)+p0+xlab('')+ylab('Fraction of total cells')+scale_fill_manual(values=cluster_cols)+guides(ncol=1)
pdf('Fig2C.pdf',w=6,h=6)
print(p2)
dev.off()


genelist<-'H2-Aa,Cd80,Cd86,Cd40,Il15,Il15ra,Cd274,Icam1,Il2ra,Itga4,Il12b,Il7r,Cd74,Tlr2,Il1b,Cd300a'
genes<-unlist(strsplit(genelist,split = ",",fixed=T))
tmp<-subset(DC,celltype %in% c('CCR7+ DC','CD209e+ DC','IFN response DC1','Mgl2+ DC'))
tmp$celltype<-factor(tmp$celltype,levels=c('CCR7+ DC','CD209e+ DC','IFN response DC1','Mgl2+ DC'))
Idents(tmp)<-'celltype'
p<-VlnPlot(tmp, features = genes,pt.size = -1,slot="data",assay="RNA",stack = T,flip =T)+theme(legend.position='none',axis.text.y = element_text(color="black",size=10,angle=0),axis.text.x = element_text(color="black",size=10,angle=30))+xlab('')
pdf('Fig2E.pdf',w=4,h=6)
print(p)
dev.off()


######FigS2 plot--------------------------------------------------------
pdf('FigS2A.pdf',w=12,h=5)
p1 <- DimPlot(DC, reduction = "umap", label = F,group.by = "celltype",split.by='Group',label.size = 4,pt.size = 1,cols=cluster_cols)+ggtitle('DC')
print(p1)
dev.off()


genes<-c('Itgam','Itgax','Xcr1')
p<-VlnPlot(DC, features = genes,pt.size = -1,slot="data",assay="RNA",stack = T,flip =T,group.by='celltype')+theme(legend.position='none',axis.text.y = element_text(color="black",size=10,angle=0),axis.text.x = element_text(color="black",size=10,angle=30))+xlab('')
pdf('FigS2B.pdf',w=6,h=3)
print(p)
dev.off()


DC$celltype<-factor(DC$celltype,levels=)
p<-VlnPlot(DC, features = 'Ccl17',pt.size = -1,slot="data",assay="RNA",stack = F,flip =T,cols=cluster_cols)+theme(legend.position='none',axis.text.y = element_text(color="black",size=10,angle=0),axis.text.x = element_text(color="black",size=10,angle=30))+xlab('')
pdf('FigS2C.pdf',w=6,h=2.5)
print(p)
dev.off()


######Fig3 plot--------------------------------------------------------
Mac<-readRDS('rdsfiles/Mac.rds')
Idents(Mac)<-'celltype'
DefaultAssay(Mac)<-'RNA'

pdf('FigS2D.pdf',w=12,h=5)
p1 <- DimPlot(Mac, reduction = "umap", label = F,group.by = "celltype",split.by='Group',label.size = 4,pt.size = 1,cols=cluster_cols)+ggtitle('Macrophage')
print(p1)
dev.off()


pdf('Fig3A.pdf',w=6.5,h=5)
p1 <- DimPlot(Mac, reduction = "umap", label = F,group.by = "celltype",label.size = 4,pt.size = 1,cols=cluster_cols)+ggtitle('Macrophage')
print(p1)
dev.off()


genelist<-'Plet1,Siglecf,Krt79,Chil3,Pparg,Ear2,Cd36,Top2a,Mki67,Tubb5,Mcm3,Mcm5,Mcm6,Spp1,Arg1,Fabp5,Gpnmb,Ly6c2,Irf7,Itgam,Trem2,Fcgr1,C5ar1,Ly6e,F13a1,Apoe,Cx3cr1,C1qa,C1qb,C1qc,Fn1,Tnip3,S100a4,S100a6,Cd209a'
genes<-unlist(strsplit(genelist,split = ",",fixed=T))
cluster.averages <- AverageExpression(object = Mac, assays ='RNA',return.seurat = F)
data<-cluster.averages$RNA
da<-log2(data[genes,]+1)
da<-log2(data[genes,]+1)
filename<-'Fig3B.pdf'
color = colorRampPalette(c('blue','white','red'))(50)
pheatmap(da,filename=filename,color = color,width=5,scale='row',height=8,cluster_rows=F,cluster_cols=F,display_numbers = F,number_format = "%.0f",fontsize = 12,fontsize_col=8,border_color=NA,angle_col = "45")


per<-round(prop.table(table(Mac@meta.data[,c('Group','celltype')]), margin=1),3)
per <-data.frame(per)
p0<-theme_bw()+theme(legend.position='right',panel.grid=element_blank(), legend.background = element_rect(colour = NA),
        legend.title=element_blank(),legend.text = element_text(face="plain", color="black",size=8),strip.text = element_text(size=20),
        axis.text.x = element_text(color="black",size=10,angle=0),
        axis.text.y = element_text(color="black",size=10),
        axis.title.x = element_text(face="plain", color="black",size=15),
        axis.title.y = element_text(face="plain", color="black",size=15))	
p2<-ggplot(data =per, aes(fill=celltype,y=Freq,x=Group)) +
geom_bar(stat = "identity",position = "stack",width=0.6)+p0+xlab('')+ylab('Fraction of total cells')+scale_fill_manual(values=cluster_cols)+guides(ncol=1)
pdf('Fig3C.pdf',w=6,h=6)
print(p2)
dev.off()



tmp<-subset(Mac,name %in% c('AM','Gpnmb+ Mac','IM'))
genelist<-'Arg1,Thbs1,P2rx4,Lamp1,Trem2,Cebpb,Zbp1,Vegfa,Braf'
genes<-unlist(strsplit(genelist,split = ",",fixed=T))
prefix<-'Mac_name3_celltype'
p<-VlnPlot(tmp, features = genes,pt.size = -1,slot="data",assay="RNA",stack = T,flip =T,group.by='name')+theme(legend.position='none',axis.text.y = element_text(color="black",size=10,angle=0),axis.text.x = element_text(color="black",size=10,angle=30))+xlab('')
pdf('Fig3E.pdf',w=4,h=5)
print(p)
dev.off()


```
#CCL_gene.list
Ccl2
Ccl3
Ccl4
Ccl5
Ccl6
Ccl7
Ccl8
Ccl9
Ccl12
Ccl17
Ccl24
Il1b
Il18
Cxcl2
Cxcl3
Cxcl10
Cxcl13
Cxcl16
```
tmp1<-subset(Mac,name %in% c('AM','Gpnmb+ Mac','IM'))
genes<-as.vector(read.table('CCL_gene.list')$V1)
p_tmp1<-theme(legend.text =  element_text(color="black",size=8),legend.title =  element_text(color="black",size=10),
axis.text.x = element_text(color="black",size=12,angle=90),
axis.text.y = element_text(color="black",size=12),
axis.title.x = element_text(face="plain", color="black",size=15),
axis.title.y = element_text(face="plain", color="black",size=15))
tmp1$name.g<-paste0(tmp1$name,'_',tmp1$Group)
tmp1$name.g<-factor(tmp1$name.g,levels=names(table(tmp1$name.g)))
tmp1 <- ScaleData(tmp1, features = genes, verbose = FALSE)
Idents(tmp1)<-'name.g'
pdf('Fig3F.pdf',w=8,h=4)
p<-DotPlot(tmp1,features=genes,cols =c("lightgrey", "red"))+xlab('')+p_tmp1+ylab('')
print(p)
dev.off()

######Fig4 plot--------------------------------------------------------
Trds<-readRDS('rdsfiles/T.rds')
Idents(Trds)<-'celltype'
DefaultAssay(Trds)<-'RNA'
Trds$Group<-factor(Trds$Group,levels=c('PBS','HDM'))

palette<-c(
    "Tifn"="#9d4edd",
    "Tcm CD4"="#F8961E",
	"Tfh"="#90be6d",
    "Treg"="#7180B9",
	"Th2"="#F9C74F",
    "Trm CD8"="#F94144",
    "Tcm CD8"="#ff99ac",
    "Tc CD8"="#3423a6",
    "Tgd"="#43aa8b",
    "NKT"="#577590")
	
pdf('Fig4A.pdf',w=6,h=5)
DimPlot(Trds,label=T,cols=palette,repel=F)
dev.off()

genes1<-c("Cd4","Cd8a","Cd8b1","Cd69","Ifit1","Ifit3","Isg15","Bst2","Sell","Tcf7","Ly6c1","Cd28","Lef1","Ccr7","Icos","Ctla4","Foxp3","Il2ra","Cd27","Cd44","Gata3","Tox2","Ly6c2","Eomes","Ccl5","Nkg7","Tbx21","Prf1","Gzma","Gzmb","Gzmk","Ifng","Il7r","Rorc","Trdc",'Trgv2',"Itgae")
Idents(Trds)<-'celltype'
Trds <- ScaleData(Trds, features = genes1, verbose = FALSE)
p0=DotPlot(Trds,features =genes1)#,'FTH1','FTL'
p3<-p0+theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1))+xlab('')+ylab('')+
  guides(color = guide_colorbar(title = 'Scaled expression',order = 1),size = guide_legend("Percent expressed"),order = 0)+
  scale_colour_gradientn(colours = c("dodgerblue1","#44c1f0", "lightgoldenrod","#e20612",'#cc340c'))+
  theme(legend.position = "top")+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.3, linetype="solid"))+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)
pdf(paste0('Fig4B.pdf'), w=9,h=5)
print(p3)
dev.off()

per<-round(prop.table(table(Trds@meta.data[,c('Group','celltype')]), margin=1),3)
per <-data.frame(per)
p0<-theme_bw()+theme(legend.position='right',panel.grid=element_blank(), legend.background = element_rect(colour = NA),
        legend.title=element_blank(),legend.text = element_text(face="plain", color="black",size=8),strip.text = element_text(size=20),
        axis.text.x = element_text(color="black",size=10,angle=0),
        axis.text.y = element_text(color="black",size=10),
        axis.title.x = element_text(face="plain", color="black",size=15),
        axis.title.y = element_text(face="plain", color="black",size=15))	
p2<-ggplot(data =per, aes(fill=celltype,y=Freq,x=Group)) +
geom_bar(stat = "identity",position = "stack",width=0.6)+p0+xlab('')+ylab('Fraction of total cells')+scale_fill_manual(values=palette)+guides(ncol=1)
pdf('Fig4C.pdf',w=4,h=6)
print(p2)
dev.off()

genes<-c('Ifng','il2','il4','il5','il6','il10','Tgfb1','il13','il17a','il17f','il21')
library(stringr)
genes<-str_to_title(genes)
#StackedVlnplot2
#pSSC & CPC
modify_vlnplot1<- function(obj,
                          features,
                          pt.size = 0,
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          cols=palette,
                          slot="data",
                          assay="RNA",
                          ...) {
  p<- VlnPlot(obj, features = features,pt.size = pt.size, cols=cols, slot=slot,
              assay=assay,... )  +
    xlab("") + ylab(features) + ggtitle("") +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title= element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = rel(0.7), angle = 0, vjust = 0.5,color="black"),
          plot.margin = plot.margin )
  return(p)
}
#add plot.title= element_blank(),?
## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0,
                          slot="data",
                          assay="RNA",
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot1(obj = obj, features = x, slot=slot, assay=assay,...))
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(color="black",size=10,angle=90),axis.ticks.x = element_line())
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

p<-StackedVlnPlot(obj = Trds, features =genes,slot="data",assay="RNA")
pdf('Fig4D.pdf',w=4,h=5)
print(p)
dev.off()
	

genelist<-'Cd44,Hif1a,Ctla4,Gata3,Stat1,Stat3,Id2,Rbpj,Rora,Runx2,Il10'
genes<-unique(c(unlist(strsplit(genelist,split = ",",fixed=T))))
p<-VlnPlot(Trds, features = genes,pt.size = -1,slot="data",assay="RNA",stack = T,flip =T,split.by='Group',cols=c('#00BFC4','#F8766D'))+theme(legend.position='right',axis.text.y = element_text(color="black",size=10,angle=0),axis.text.x = element_text(color="black",size=10,angle=30))+xlab('')
pdf(paste0('Fig4I.pdf'),w=8,h=6)
print(p)
dev.off()


pdf('FigS4A.pdf',w=12,h=5)
DimPlot(Trds,split.by="Group",label=T,cols=palette,repel=T)
dev.off()



######Fig5 plot--------------------------------------------------------
nonIC<-readRDS('rdsfiles/Fig5.nonIC_celltypes_Major.rds')
Idents(nonIC)<-'celltype'
DefaultAssay(nonIC)<-'RNA'
nonIC$Group<-factor(nonIC$Group,levels=c('PBS','HDM'))

palette1=c(
'EC gCap'="#9ccc9c",
'EC arterial'="#ff1100",
'EC pVein'="#187bcd",
'EC sVein'="#ff71b5",
'EC aCap'="#ff8d38",
'EC Lym'="#046307")

palette2=c(
'AT2'="#ffd319",
'Ciliated'="#5AAE61",
'Club'="#95F5F5",
'AT1'="#A50026",
'Goblet'="#B2ABD2",
'Basal'="#05709D",
'Neuroendocrine'="#D6604D")

palette3=c(
"Alveolar Fibroblast" = "#EF7D87",
"Peribronchial Fibroblast" ="#91B5D6",
"Pericyte"="#4DAF4A",
"Adventitial Fibroblast"="#F29403",
"SMC"="#C39BC9")
palette<-c(palette1,palette2,palette3)


pdf('Fig5A.pdf',w=8,h=5)
p1 <- DimPlot(nonIC, reduction = "umap", label = T,group.by = "celltype",label.size = 3,pt.size = 0,repel=T,cols=palette)+ggtitle('')
print(p1)
dev.off()



genelist1<-'Pecam1,Cdh5,Cldn5,Kdr,Car4,Tbx2,Ednrb,Fibin,Cd93,Plvap,Gpihbp1,Vwf,Ltbp4,Nrp2'
genelist2<-'Epcam,Alcam,Wfdc2,Ndnf,Cd74,H2-Ab1,Sftpc,Sftpb,Cxcl15,Krt15,Krt5,Tmem212,Ccdc153,Tuba1a,Cldn10,Retnla,Piezo2,Resp18'
genelist3<-'Col1a1,Col1a2,Col3a1,Fmo2,Npnt,Lgr5,Hhip,Aspn,Cesld,Cspg4,Pdgfrb,Pi16,Dcn,Acta2,Cnn1,Tagln'
genes1<-unique(c(unlist(strsplit(genelist1,split = ",",fixed=T))))
genes2<-unique(c(unlist(strsplit(genelist2,split = ",",fixed=T))))
genes3<-unique(c(unlist(strsplit(genelist3,split = ",",fixed=T))))
genes<-unique(c(genes1,genes2,genes3))

Idents(nonIC)<-'celltype'
nonIC <- ScaleData(nonIC, features = genes1, verbose = FALSE)
p0=DotPlot(nonIC,features =genes)#,'FTH1','FTL'
p3<-p0+theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1))+xlab('')+ylab('')+
  guides(color = guide_colorbar(title = 'Scaled expression',order = 1),size = guide_legend("Percent expressed"),order = 0)+
  scale_colour_gradientn(colours = c("dodgerblue1","#44c1f0", "lightgoldenrod","#e20612",'#cc340c'))+
  theme(legend.position = "top")+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.3, linetype="solid"))+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)
pdf('Fig5B.pdf', w=12,h=7)
print(p3)
dev.off()

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
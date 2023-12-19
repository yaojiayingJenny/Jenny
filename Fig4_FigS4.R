library(Seurat)
library(pheatmap)
library(ggplot2)
library(CellChat)


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

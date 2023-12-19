library(Seurat)
library(pheatmap)
library(ggplot2)
library(CellChat)


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
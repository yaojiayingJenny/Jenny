library(Seurat)
library(pheatmap)
library(ggplot2)
library(CellChat)



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

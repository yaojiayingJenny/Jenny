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

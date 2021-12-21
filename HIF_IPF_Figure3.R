BiocManager::install('GSVA')
library(GSVA)
library(survival)
library(survminer)

####GSE40839 ####

GSE38958<-getGEO('GSE38958',destdir = '.')
GSE38958<-GSE38958[[1]]
GSE38958<-ReadAffy()

GSE38958_data<-exprs(GSE38958)
boxplot(GSE38958_data)
boxplot(GSE38958_data,outline=FALSE)
GSE38958_data<-affy::rma(GSE38958)
GSE38958_data<-exprs(GSE38958_data)


bb<-"./GPL5175.soft"
bb_nn <- grep("^[^#!^]", readLines(bb))[1] - 1
pfinfo_bb <- read.table(bb, sep = "\t", quote = "", header = TRUE, skip = bb_nn, fill = TRUE)
colnames(pfinfo_bb)
pfinfo_bb <- pfinfo_bb[, c(1,11)]
write.csv(pfinfo_bb,'pfinfo_bb.csv')
pfinfo_bb<-read.csv('pfinfo_bb.csv')
pfinfo_bb<-pfinfo_bb[,1:2]

GSE38958_data<-as.data.frame(GSE38958_data)
GSE38958_data$ID<-rownames(GSE38958_data)
GSE38958_data<-merge(GSE38958_data,pfinfo_bb,by='ID')
GSE38958_CHICAGO<-GSE38958[,c(1:76)]
GSE38958_CHICAGO<-as.matrix(GSE38958_CHICAGO[,-1])

GSE38958_PITTSBURGH<-GSE38958[,c(1,77:121)]
GSE38958_PITTSBURGH<-as.matrix(GSE38958_PITTSBURGH[,-1])


gsva_Hypoxia <- gsva(as.matrix(GSE38958),Hypoxia_list , mx.diff=1)
HIF1A<-as.data.frame((gsva_Hypoxia))
hclust_completem <- hclust(dist(as.data.frame(t(HIF1A))), method = "ward.D")

## Figure 3B  
HIF_signature<-GSE38958_data[Hypoxia_list,]
rownames(gsva_Hypoxia)<-'GSVA_score'
HIF_signature<-HIF_signature[,hclust_completem$order]
colnames(HIF_signature)<-HIF_signature[8,]
HIF_signature<-HIF_signature[-8,]
a3<-as.matrix(HIF_signature)
f1<-function(x){
  x<-as.numeric(x)
}
a3<-apply(a3, 2, f1)
rownames(a3)<-rownames(HIF_signature)
a3 <- a3[,order(a3[1,])]
class(a3[3,1])
colorbar<-colorRampPalette(c('darkblue','grey','red'))(n=1000)
distCor <- function(a3) as.dist(1-cor(t(a3)))
hclustAvg <- function(a3) hclust(a3, method="average") #hclustfun=hclustAvg???distfun=distCor
hc<-hclust(as.dist(1-cor(values, method="pearson")), method="average")
png('heatmap_HIF.png',width = 20000,height = 10000,res = 1200)
heatmap.2(a3, trace="none", density='none',scale="row", margins = c(10,12),cexRow = 2,
          cexCol = 0.6, zlim=c(-10,10),srtCol=45,adjCol=c(1,0),Colv = F,Rowv = F,
          distfun=distCor, col=colorbar, symbreak=FALSE,key = T,keysize = 0.8,
          ColSideColors  = col_labels) 
dev.off()}#### HIF-siganture heatmap####

## Figure 3D
Oxidative_list<-list(c('ABCC1', 'CDKN2D', 'FES', 'GCLC', 'GCLM', 'GLRX2', 'HHEX', 'IPCEF1', 
                       'JUNB', 'LAMTOR5', 'LSP1', 'MBP', 'MGST1', 'MPO', 'NDUFA6', 'PFKP', 'PRDX1', 'PRDX2','PRDX4', 
                       'PRNP', 'SBNO2','SCAF4', 'SOD1', 'SOD2', 'RXN1', 'TXN', 'TXNRD1'))
gsva_Oxidative <- gsva(raw_exprSet_final_batch,Oxidative_list , mx.diff=1)

cor.test(as.numeric(gsva_Oxidative),as.numeric(gsva_Hypoxia))

### merge list HIF score and Oxidative score
df<-rbind(gsva_Oxidative,gsva_Hypoxia)
df<-t(df)
df<-as.data.frame(df)
colnames(df)<-c('Oxidative_Stree_Score','HIF_Score')
## plot
pdf('Oxidative_Stree_HIF.pdf',width = 10,height = 10)
par(mar=c(5,8,100,5))
ggscatter(df, x = "Oxidative_Stree_Score", y = "HIF_Score", 
          size = 5,
          xlab = expression(paste('HIF',' Score')),
          ylab = expression(paste('Oxidative Stress',' Score')),
          add = "reg.line",  # Add regressin line
          #cor.coef = TRUE,
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          
          #conf.int = TRUE, # Add confidence interval
          cor.coeff.args = list(method = "pearson", label.x =-0.07 , label.sep = "\n",
                                label.y = 1.1),cor.coef.size = 10)+
  annotate("text", x =-0.025, y=1.0, label = substitute(paste(italic("R"), " = 0.48")),size=10)+
  annotate("text", x = -0.025, y=0.9, label = substitute(paste(italic("P"), " = 3.5 x ",10^-8)),size=10)+
  theme(axis.title=element_text(size=40,face="bold"),axis.text=element_text(size=35,face="bold"),axis.title.x=element_text(face = 'bold'),
        axis.title.y=element_text(face = 'bold'))
dev.off()

library(GEOquery)
library(survival)
library(GSVA)
setwd('/Users/lefan/Rcode/IPF-singlecell/GSE70866_RAW_BAL/')
#### download data ####
gse70867<-getGEO('GSE70867',destdir = '.')
rm(gse70867)
## check data, it will be displayed as a list, using list[] to extract elments
str(gse70867)
#### chose the data we need ####
b =gse70867[[2]]
boxplot(b)
### get clinical information   
phe=pData(b)
rownames(phe)
rownames(phe)
phe<-phe[23:154,1:17]
### get expression matrix
raw_exprSet=exprs(b) 
b
### delete healthy samples
raw_exprSet<-raw_exprSet[,23:154]
dim(raw_exprSet)
raw_exprSet<-na.omit(raw_exprSet)

### get expression matrix
c = gse70867[[3]]
### get clinical information
phe_c=pData(c)
phe_c<-phe_c[,1:17]
rm(b,c,phe_c,phe,raw_exprSet,raw_exprSet_c,raw_exprSet_final_numeric)
### merge two clinical information
phe_final<-rbind(phe,phe_c)

write.csv(phe_final,'phe_final.csv') ## adjust in excel

phe_final<-read.csv('phe_final_IPF.csv',header = T,sep = ',')
phe_final_HDIPF<-read.csv('phe_final_HDIPF.csv',header = T,sep = ',')
phe_final_HDIPF_COPD<-read.csv('phe_final_HDIPF_COPD.csv',header = T,sep = ',')
rm(phe_final_HDIPF_COPD)

###
raw_exprSet_c=exprs(c) 
dim(raw_exprSet_c)
raw_exprSet_c<-na.omit(raw_exprSet_c)

### COPD patients and healthy volunteer
d = gse70867[[4]]
phe_d=pData(d)
phe_d<-phe_d[,1:17]
raw_exprSet_d=exprs(d) 
?cbind
raw_exprSet_d<-na.omit(raw_exprSet_d)
dd<-"./GPL570.soft"
dd_nn <- grep("^[^#!^]", readLines(dd))[1] - 1
pfinfo_dd <- read.table(dd, sep = "\t", quote = "", header = TRUE, skip = dd_nn, fill = TRUE)
pfinfo_dd <- pfinfo_dd[, c(1,11)]
write.csv(pfinfo_dd,'pfinfo_dd.csv')

for (i in 1:54676){
  i<-1
  pfinfo_dd[i,2] <-strsplit(pfinfo_dd[i,2],' /// ')[[1]][1]  
}
pfinfo_dd<-read.csv('pfinfo_dd.csv',sep = ',')
pfinfo_dd<-pfinfo_dd[,1:2]
?strsplit
raw_exprSet_d<-as.data.frame(raw_exprSet_d)
raw_exprSet_d$ID<-rownames(raw_exprSet_d)
raw_exprSet_d<-merge(raw_exprSet_d,pfinfo_dd,by='ID')
raw_exprSet_d<-raw_exprSet_d[,-1]
raw_exprSet_d<-aggregate(x=raw_exprSet_d[,1:(ncol(raw_exprSet_d)-1)],by=list(raw_exprSet_d$Gene.Symbol),FUN=mean)
rownames(raw_exprSet_d)<-raw_exprSet_d$Group.1
raw_exprSet_d<-as.matrix(raw_exprSet_d[,-1])
save(raw_exprSet_d,pfinfo_dd,file = './COPD/COPD.Rdata')
pdf('raw_exprSet_d.pdf')
boxplot(raw_exprSet_d)
dev.off()
write.csv(raw_exprSet_d,'COPD_nor_expression.csv')
gsva_Hypoxia <- gsva(raw_exprSet_d,Hypoxia_list , mx.diff=1)
gsva_Hypoxia<-as.data.frame(t(gsva_Hypoxia))
write.csv(gsva_Hypoxia,'COPD_gsva_Hypoxia.csv')
rm(dd,dd_nn,pfinfo_dd,raw_exprSet_d)

#### Probe ID transfer ####
bb <- "./GPL14550.soft"
cc<-"./GPL17077.soft"
bb_nn <- grep("^[^#!^]", readLines(bb))[1] - 1
cc_nn <- grep("^[^#!^]", readLines(cc))[1] - 1
pfinfo_bb <- read.table(bb, sep = "\t", quote = "", header = TRUE, skip = bb_nn, fill = TRUE)
pfinfo_bb <- pfinfo_bb[, c(1,7)]
pfinfo_cc <- read.table(cc, sep = "\t", quote = "", header = TRUE, skip = cc_nn, fill = TRUE)
pfinfo_cc <- pfinfo_cc[, c(1,7)]
rm(bb,cc,bb_nn,cc_nn,pfinfo_bb,pfinfo_cc)
?grep
## bb group
raw_exprSet<-as.data.frame(raw_exprSet)
raw_exprSet$ID<-rownames(raw_exprSet)
raw_exprSet<-merge(raw_exprSet,pfinfo_bb,by='ID')
raw_exprSet<-raw_exprSet[,-1]
### get average value of genes
raw_exprSet<-aggregate(x=raw_exprSet[,1:(ncol(raw_exprSet)-1)],by=list(raw_exprSet$GENE_SYMBOL),FUN=mean)
## cc group
raw_exprSet_c<-as.data.frame(raw_exprSet_c)
raw_exprSet_c$ID<-rownames(raw_exprSet_c)
colnames(raw_exprSet_c)
raw_exprSet_c<-merge(raw_exprSet_c,pfinfo_cc,by='ID')
raw_exprSet_c<-raw_exprSet_c[,-1]
raw_exprSet_c<-aggregate(x=raw_exprSet_c[,1:(ncol(raw_exprSet_c)-1)],by=list(raw_exprSet_c$GENE_SYMBOL),FUN=mean)
rm(raw_exprSet_c)

#### microarray_unique_value function ####
microarray_unique_value<-function(x,proble_id,Func){
  x<-as.data.frame(x)
  x$ID<-rownames(x)
  colnames(x)
  x<-merge(x,proble_id,by='ID')
  x<-x[,-1]
  x<-aggregate(x=x[,1:(ncol(x)-1)],by=list(x$GENE_SYMBOL),FUN=Func)
  return(x)
}
colnames(raw_exprSet)
rownames(raw_exprSet)[10]
microarray_unique_value(raw_exprSet,pfinfo_bb,mean)
microarray_unique_value(raw_exprSet_c,pfinfo_cc,mean)
### merge two expression matrix
raw_exprSet_final_numeric<- merge(raw_exprSet, raw_exprSet_c, by = "Group.1")
#raw_exprSet_final_numeric<- merge(raw_exprSet_final_numeric, raw_exprSet_d, by = "Group.1")

rownames(raw_exprSet_final_numeric)<-raw_exprSet_final_numeric$Group.1 # add rowname
raw_exprSet_final_numeric<-raw_exprSet_final_numeric[,-1] # delete 1st col and 1st row
### convert character into number
raw_exprSet_final_numeric[, c(1:length(raw_exprSet_final_numeric[1,]))] <- sapply(raw_exprSet_final_numeric[, c(1:length(raw_exprSet_final_numeric[1,]))], as.numeric)
raw_exprSet_final_numeric<-as.matrix(raw_exprSet_final_numeric)
# check
colnames(raw_exprSet_final_numeric)
class(raw_exprSet_final_numeric[1,1])
boxplot(raw_exprSet_final_numeric)

#### Batch effect ####
BiocManager::install('sva')
library(sva)
### create batch group
GSE27957_batch<-c(rep(1,82),rep(2,64))
phe_final$status<-c(rep(1,62),rep(2,50),rep(3,64))
batch<-phe_final$status
phe_final_HDIPF_COPD$status<-c(rep(1,82),rep(2,50),rep(3,64),rep(4,57))

batch<-phe_final_HDIPF_COPD$status
dim(raw_exprSet_final_numeric)
### remove batch
?ComBat
raw_exprSet_final_batch<- ComBat(dat=raw_exprSet_final_numeric, batch=batch)
?ComBat
raw_exprSet_final_batch<-as.matrix(raw_exprSet_final_batch)
save(raw_exprSet_final_batch,phe_final_HDIPF,file = 'GSE70866_HDIPF_rmBatch_matrix.Rdata')

Hypoxia_list<-list(c('NDRG1','ENO1','VEGFA','MRPS17','TPI1','CDKN3',
                     'MIF','LDHA','ALDOA','TUBB6','PGAM1','SLC2A1','P4HA1','ACOT7','ADM'))

gsva_Hypoxia <- gsva(raw_exprSet_final_batch,Hypoxia_list , mx.diff=1)
gsva_Hypoxia <- gsva(Freiburg,Hypoxia_list , mx.diff=1)


Oxidative_list<-list(c('ABCC1', 'CDKN2D', 'FES', 'GCLC', 'GCLM', 'GLRX2', 'HHEX', 'IPCEF1', 
                       'JUNB', 'LAMTOR5', 'LSP1', 'MBP', 'MGST1', 'MPO', 'NDUFA6', 'PFKP', 'PRDX1', 'PRDX2','PRDX4', 
                       'PRNP', 'SBNO2','SCAF4', 'SOD1', 'SOD2', 'RXN1', 'TXN', 'TXNRD1'))
BAL_list<-list(c("BMP6", "IBSP", 'S100A14', 'STAB1', 'LOC284751', 'ANKRD22'))
gsva_BAL<-gsva(raw_exprSet_final_batch,BAL_list , mx.diff=1)
gsva_Hypoxia <- gsva(raw_exprSet_final_batch,Hypoxia_list , mx.diff=1)
gsva_Oxidative <- gsva(raw_exprSet_final_batch,Oxidative_list , mx.diff=1)

cor.test(as.numeric(gsva_Oxidative),as.numeric(gsva_Hypoxia)) # p= 1.341e-11,R=0.4814433 BAL | p-value = 3.532e-08,R=0.477372 PBMC
# p= 7.864e-15, R=0.898941, GSE92592
cor.test(as.numeric(gsva_BAL),as.numeric(gsva_Hypoxia))
### merge list score and LOXL2_PLOD2 score
df<-rbind(gsva_Oxidative,gsva_Hypoxia)
df<-t(df)
df<-as.data.frame(df)
colnames(df)<-c('Oxidative_Stree_Score','HIF_Score')

## plot
pdf('Oxidative_Stree_HIF_GSE92592.pdf',width = 10,height = 10)
par(mar=c(5,8,100,5))
ggscatter(df, x = "Oxidative_Stree_Score", y = "HIF_Score", # p= 1.341e-11,R=0.4814433 |p-value = 3.532e-08,R=0.477372 PBMC
          size = 5,
          xlab = expression(paste('HIF',' Score')),
          #xlab = expression(paste('TGF-',beta,' Score')),
          ylab = expression(paste('Oxidative Stress',' Score')),
          #ylab = expression(italic(LOXL2)~+~italic(PLOD2)~  Score),
          add = "reg.line",  # Add regressin line
          #cor.coef = TRUE,
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          
          #conf.int = TRUE, # Add confidence interval
          cor.coeff.args = list(method = "pearson", label.x =-0.07 , label.sep = "\n",
                                label.y = 1.1),cor.coef.size = 10)+
  annotate("text", x =-0.025, y=1.0, label = substitute(paste(italic("R"), " = 0.89")),size=10)+
  annotate("text", x = -0.025, y=0.9, label = substitute(paste(italic("P"), " = 7.9 x ",10^-15)),size=10)+
  #annotate("text", x = 0.407, y=-0.75, label = substitute(paste(italic("P"), " = 0.03 ")),size=15)+
  theme(axis.title=element_text(size=40,face="bold"),axis.text=element_text(size=35,face="bold"),axis.title.x=element_text(face = 'bold'),
        axis.title.y=element_text(face = 'bold'))
dev.off()


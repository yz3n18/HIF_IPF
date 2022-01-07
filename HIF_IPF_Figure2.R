Hypxoia_score<-HIF1A
HIF_score<-HIF1A


## Figure 2A   

HIF_signature<-raw_exprSet_final_batch[Hypoxia_list,]
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

## Figure 2B
if(T){
  GSVA_survival<-function(gene_list,location,Kmean){
    if  (grepl(toupper(location),'THREEREGIONS')){
      gsva_Hypoxia <- gsva(raw_exprSet_final_batch,gene_list , mx.diff=1)
      location<-'Freiburg+Siena+Leuven'
    } else if (grepl(toupper(location),'FREIBURG')){
      gsva_Hypoxia <- gsva(Freiburg,gene_list , mx.diff=1)
      location<-'Freiburg'
    } else if (grepl(toupper(location),'SIENA')){
      gsva_Hypoxia <- gsva(Siena,gene_list , mx.diff=1)
      location<-'Siena'
    } else if (grepl(toupper(location),'LEUVEN')){
      gsva_Hypoxia <- gsva(Leuven,gene_list , mx.diff=1)
      location<-'Leuven'
    }
    HIF1A<-as.data.frame((gsva_Hypoxia))
    hclust_completem <- hclust(dist(as.data.frame(t(HIF1A))), method = "ward.D")
    hif1a_es <- cutree(hclust_completem,k=Kmean)
    hif1a_es <- hif1a_es[hclust_completem$order]
    hif1a_es<-as.data.frame(hif1a_es)
    hif1a_es$X<-rownames(hif1a_es)
    phe_final_HIF1a<-merge(hif1a_es,phe_final,by='X')
    if (Kmean >2){
      phe_final_HIF1a<-phe_final_HIF1a[phe_final_HIF1a[,2] < 2 |phe_final_HIF1a[,2] > 2 ,]
    }
    #write.csv(phe_final_HIF1a,paste0('GSVA_HIF_final_','_',Kmean,'.csv'))
    Surv_data<-Surv(time = phe_final_HIF1a$day,event =phe_final_HIF1a$event )
    a<-survdiff(Surv_data~phe_final_HIF1a$hif1a_es)
    a_pvalue<-sprintf("%.2e", pchisq(a$chisq, length(a$n)-1, lower.tail = FALSE))
    a_pvalue
    Surv_fit<-survfit(Surv_data~phe_final_HIF1a$hif1a_es)
    pdf(paste0(location,'_','_','P=',as.character(a_pvalue),'_K',
               as.character(Kmean),'_GSVA_HIF.pdf'),width = 8,height = 8)
    par(mar=c(6,6,4,2))
    ## low-1
    HIF1A_1<-HIF1A[,phe_final_HIF1a[order(phe_final_HIF1a[,2],decreasing = F)[1],1]]
    ## high-2
    HIF1A_2<-HIF1A[,phe_final_HIF1a[order(phe_final_HIF1a[,2],decreasing = T)[1],1]]
    if (HIF1A_1 > HIF1A_2){
      plot(Surv_fit,conf.int='none',col=c('red','blue'),
           lwd=2,mark.time=T,cex.lab=3,cex.axis=2.5,main=location)#,xlab='\n\nTime to death (days)')
      legend('bottomleft',c(paste("High ", "(n=", table(phe_final_HIF1a$hif1a_es)[1], ")", sep=""),
                            paste("Low ", "(n=", table(phe_final_HIF1a$hif1a_es)[2], ")", sep="")),
             col = c('red','blue'),lwd = 3,bty='n',cex=1)
    } else {
      plot(Surv_fit,conf.int='none',col=c('blue','red'),
           lwd=2,mark.time=T,cex.lab=3,cex.axis=2.5,main=location)#,xlab='\n\nTime to death (days)'
      legend('bottomleft',c(paste("Low ", "(n=", table(phe_final_HIF1a$hif1a_es)[1], ")", sep=""),
                            paste("High ", "(n=", table(phe_final_HIF1a$hif1a_es)[2], ")", sep="")),
             col = c('blue','red'),lwd = 3,bty='n',cex=1)
    }
    title(ylab="Survival", line=3,cex.lab=3,cex.axis=2.5)
    title(xlab="Time to death (days)", line=3.5,cex.lab=3,cex.axis=2.5)
    dev.off()  
  }
}
GSVA_survival(HIF_list,'THREEREGIONS',2)
# Figure 2C
library(survminer)

phe_final_HIF1a<-read.csv('phe_final.csv',sep= ',',header = T)
single_line<-Surv(time = phe_final_HIF1a$day,event = phe_final_HIF1a$event)

phe_final_HIF1a$single_line<- with(phe_final_HIF1a,single_line)
Gcox<-coxph(single_line~Gap_Stage+HIFScore,
            data = phe_final_HIF1a)
ggforest(Gcox,data = phe_final_HIF1a,cpositions = c(0.02,0.22,0.4),main = 'Hazard ratio',
         fontsize =1.8,refLabel = 'reference',noDigits = 2)

# Figure 2D
if(T){
  phe_function<-function(gene_list,Kmean){
    gsva_Hypoxia <- gsva(raw_exprSet_final_batch,gene_list,mx.diff=1)
    HIF1A<-as.data.frame((gsva_Hypoxia))
    hclust_completem <- hclust(dist(as.data.frame(t(HIF1A))), method = "ward.D")
    hif1a_es <- cutree(hclust_completem,k=Kmean)
    hif1a_es <- hif1a_es[hclust_completem$order]
    hif1a_es<-as.data.frame(hif1a_es)
    hif1a_es$X<-rownames(hif1a_es)
    phe_final_HIF1a<-merge(hif1a_es,phe_final,by='X')
    if (Kmean >2){
      phe_final_HIF1a<-phe_final_HIF1a[phe_final_HIF1a[,2] < 2 |phe_final_HIF1a[,2] > 2 ,]
    }
    return( phe_final_HIF1a)
  }
}
phe_function(HIF_list,2)


# Figure 4A
# chicago /pittsburgh

if(T){
  GSVA_survival<-function(gene_list,location,Kmean){
    
    if  (grepl(toupper(location),'TWOREGIONS')){
      gsva_Hypoxia <- gsva(as.matrix(GSE38958),gene_list , mx.diff=1)
      location<-'Chigaco+Pittsburgh'
    } else if (grepl(toupper(location),'CHICAGO')){
      gsva_Hypoxia <- gsva(as.matrix(GSE38958_CHICAGO),gene_list , mx.diff=1)
      location<-'Chigaco'
    } else if (grepl(toupper(location),'PITTSBURGH')){
      gsva_Hypoxia <- gsva(as.matrix(GSE27957_PITTSBURGH),gene_list , mx.diff=1)
      location<-'Pittsburgh'
    } 
    HIF1A<-as.data.frame((gsva_Hypoxia))
    hclust_completem <- hclust(dist(as.data.frame(t(HIF1A))), method = "ward.D")
    hclust_completem <- hclust(dist(HIF1A), method = "ward.D")
    hif1a_es <- cutree(hclust_completem,k=Kmean)
    hif1a_es <- hif1a_es[hclust_completem$order]
    hif1a_es<-as.data.frame(hif1a_es)
    hif1a_es$X<-rownames(hif1a_es)
    #phe_final_HIF1a<-merge(HIF1A,GSE38958_phe,by='X')
    phe_final_HIF1a<-merge(hif1a_es,GSE38958_phe,by='X')
    if (Kmean >2){
      phe_final_HIF1a<-phe_final_HIF1a[phe_final_HIF1a[,2] < 2 |phe_final_HIF1a[,2] > 2 ,]
    }
    Surv_data<-Surv(time = phe_final_HIF1a$Days,event =phe_final_HIF1a$Event )
    a<-survdiff(Surv_data~phe_final_HIF1a$hif1a_es)
    a
    a_pvalue<-sprintf("%.2e", pchisq(a$chisq, length(a$n)-1, lower.tail = FALSE))
    Surv_fit<-survfit(Surv_data~phe_final_HIF1a$hif1a_es)
    pdf(paste0(location,'_','_','P=',as.character(a_pvalue),'_K',
               as.character(Kmean),'_GSVA_HIFSCORE.pdf'),width = 8,height = 8)
    par(mar=c(5,5,4,2))
    ## low-1
    HIF1A_1<-HIF1A[,phe_final_HIF1a[order(phe_final_HIF1a[,2],decreasing = F)[1],1]]
    ## high-2
    HIF1A_2<-HIF1A[,phe_final_HIF1a[order(phe_final_HIF1a[,2],decreasing = T)[1],1]]
    if (HIF1A_1 > HIF1A_2){
      plot(Surv_fit,conf.int='none',col=c('red','blue'),
           lwd=2,mark.time=T,xlab='Time to death (days)',ylab='Survival',cex.lab=2,cex.axis=1.5,main=location)
      legend('bottomleft',c(paste("High ", "(n=", table(phe_final_HIF1a$hif1a_es)[1], ")", sep=""),
                            paste("Low ", "(n=", table(phe_final_HIF1a$hif1a_es)[2], ")", sep="")),
             col = c('red','blue'),lwd = 2,bty='n',cex=1.5)
    } else {
      plot(Surv_fit,conf.int='none',col=c('blue','red'),
           lwd=2,mark.time=T,xlab='Time to death (days)',ylab='Survival',cex.lab=2,cex.axis=1.5,main=location)
      legend('bottomleft',c(paste("Low ", "(n=", table(phe_final_HIF1a$hif1a_es)[1], ")", sep=""),
                            paste("High ", "(n=", table(phe_final_HIF1a$hif1a_es)[2], ")", sep="")),
             col = c('blue','red'),lwd = 2,bty='n',cex=1.5)
    }
 
    dev.off()  
    return(HIF1A)
  }
}
GSVA_survival(HIF_list,'TWOREGIONS',2)

# Figure 4B
HIF1A<-as.data.frame(t(HIF1A))
HIF1A$X<-rownames(HIF1A)
phe_final_HIF1a<-merge(HIF1A,phe_final_HIF1a,by='X')
colnames(phe_final_HIF1a)[2]<-'HIF_Score'

single_line<-Surv(time = phe_final_HIF1a$day,event = phe_final_HIF1a$event) 
Gcox<-coxph(single_line~ HIF_Score+Gender+Age,data = phe_final_HIF1a)# uni
Gcox
summary(Gcox)  #output provides HR CIs
ggforest(Gcox,data = phe_final_HIF1a,cpositions = c(0.02,0.22,0.4),main = 'Hazard ratio',
         fontsize =1.8,refLabel = 'reference',noDigits = 2)

## Figure 4-figure supplement
GSVA_survival(HIF_list,'CHICAGO',2)
GSVA_survival(HIF_list,'PITTSBURGH',2)

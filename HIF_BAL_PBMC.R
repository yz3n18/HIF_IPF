BiocManager::install('GSVA')
library(GSVA)
library(survival)
library(survminer)
raw_exprSet_final_batch<-readRDS('/Users/yihuawang//Rcode/UFPM/raw_exprSet_final_batch.rds')
load('/Users/yihuawang/IPF-singlecell/GSE70866_RAW_BAL/raw_exprSet_final_numeric.Rdata')

Freiburg<-raw_exprSet_final_batch[,1:62]
Siena<-raw_exprSet_final_batch[,63:112]
Siena<-raw_exprSet_final_numeric[,63:112]
Leuven<-raw_exprSet_final_batch[,113:176]
#Leuven<-raw_exprSet_final_numeric[,113:176]

phe_final1<-read.csv('phe_final_IPF.csv')
phe_final<-read.csv('/Users/yihuawang/Rcode/IPF-singlecell/GSE108975_macrophages/phe_final.csv')
Hypoxia_list<-list(c('NDRG1','ENO1','VEGFA','MRPS17','TPI1','CDKN3',
                     'MIF','LDHA','ALDOA','TUBB6','PGAM1','SLC2A1','P4HA1','ACOT7','ADM'))

gsva_Hypoxia <- gsva(raw_exprSet_final_batch,Hypoxia_list , mx.diff=1)
gsva_Hypoxia <- gsva(Freiburg,Hypoxia_list , mx.diff=1)
gsva_Hypoxia <- gsva(Siena,Hypoxia_list , mx.diff=1)
gsva_Hypoxia <- gsva(Leuven,Hypoxia_list , mx.diff=1)
HIF1A<-as.data.frame((gsva_Hypoxia))
hclust_completem <- hclust(dist(as.data.frame(t(HIF1A))), method = "ward.D")



hif1a_es <- cutree(hclust_completem,k=2)
hif1a_es<-as.data.frame(hif1a_es)
hif1a_es$X<-rownames(hif1a_es)
#phe_final_HIF1a<-merge(HIF1A,GSE27957_phe,by='X')
phe_final_HIF1a<-merge(hif1a_es,phe_final,by='X')# BAL

phe_final_HIF1a<-phe_final_HIF1a[grepl( 'Leuven', phe_final_HIF1a$title, fixed = TRUE),] # Freiburg/Siena/Leuven

#write.csv(phe_final_HIF1a,paste0('GSVA_final_1506','_',Kmean,'.csv'))
Surv_data<-Surv(time = phe_final_HIF1a$day,event =phe_final_HIF1a$event ) # BAL
a<-survdiff(Surv_data~phe_final_HIF1a$hif1a_es)


a_pvalue<-sprintf("%.2e", pchisq(a$chisq, length(a$n)-1, lower.tail = FALSE))


Surv_fit<-survfit(Surv_data~phe_final_HIF1a$hif1a_es)
fit<-survfit(Surv(day,event) ~ hif1a_es, data=phe_final_HIF1a)

KMsurvival_plot<-ggsurvplot(fit,data=phe_final_HIF1a,pval = TRUE, #show p-value of log-rank test，显示log-rank分析得到的P值
                            #conf.int = TRUE, #添加置信区
                            pval.size=10,
                            legend.labs =  c( "HIF score Low",'HIF score High'),
                            main = 'BAL all cohorts',
                            legend.title='', 
                            xlab = "Time to death (days)",   ###  customize X axis label.自定义x的time in years
                            #xlim=c(0,50),
                            break.x.by=500, ###改变坐标轴间距
                            ylab=paste0('Overall survival'),
                            #ylab=paste0('Overall peak day'),
                            #ylab=paste0('SRAS−CoV−2 RNA +'),
                            #surv.median.line = "hv", #添加中位生存时间的线
                            palette = c( "turquoise","indianred1"), ###  自定义颜色
                            #font.main = c(16, "bold", "darkblue"),
                            font.main = c(40, "bold"),
                            font.x = c(40,'black'), # X 轴
                            font.y = c(35,'black'), # c(14, "bold.italic", "darkred"), y 轴特
                            font.tickslab = 30,# c(12, "plain", "darkgreen"), 坐标轴大小
                            #conf.int.style = "step",  ###  customize style of confidence intervals,改变置信区间的样子
                            risk.table = "abs_pct",  ###  absolute number and percentage at risk，这里以n(%)的形式展示risk table
                            risk.table.y.text.col = T,###  colour risk table text annotations.
                            risk.table.y.text = F,###  show bars instead of names in text annotations in legend of risk table.不显示legend名字
                            tables.y.text=F,
                            tables.x.text =F,
                            #risk.table.title="My title", ## 改变 title
                            fontsize=7, ## 表格中数据大小
                            ncensor.plot = F, #我这里不显示删失的TRUE就显示
                            #tables.theme=theme_cleantable(), # 取消table边框
                            ggtheme = theme_classic()#绘图主题
)
KMsurvival_plot$plot<-KMsurvival_plot$plot+
  theme(legend.text = element_text(size = 35),plot.margin = unit(c(2,2,0,2), "cm"))

KMsurvival_plot$table<-KMsurvival_plot$table+labs(x = NULL, y = NULL)+theme(axis.text=element_text(size=35),
                                                                            axis.title=element_text(size=35),
                                                                            
                                                                            legend.text = element_text(size = 35),
                                                                            plot.title = element_text(size=35),
                                                                            plot.margin = unit(c(0,2,2,2), "cm")) # 改变numer at risk 大小)
KMsurvival_plot
# Freiburg/Siena/Leuven
ggsave(file=paste0('Freiburg','_HIF_ggsurvplot.pdf'),width = 15,height = 12,print(KMsurvival_plot),
       onefile=FALSE)

HIF1A<-as.data.frame(t(HIF1A))
HIF1A$X<-rownames(HIF1A)
phe_final_HIF1a<-merge(HIF1A,phe_final_HIF1a,by='X')
colnames(phe_final_HIF1a)[2]<-'HIF_Score'

single_line<-Surv(time = phe_final_HIF1a$day,event = phe_final_HIF1a$event) # BAL
Gcox<-coxph(single_line~ HIF_Score,data = phe_final_HIF1a)# uni
Gcox
summary(Gcox)  #output provides HR CIs


#### IPA ####
table(phe_final_HIF1a$hif1a_es)
raw_exprSet_final_batch_high<-raw_exprSet_final_batch[,phe_final_HIF1a[phe_final_HIF1a$hif1a_es==2,1]]
raw_exprSet_final_batch_low<-raw_exprSet_final_batch[,phe_final_HIF1a[phe_final_HIF1a$hif1a_es==1,1]]
BAL_batchmatrix_limma<-cbind(raw_exprSet_final_batch_low,raw_exprSet_final_batch_high)
mean(raw_exprSet_final_batch_high['ADM',])
mean(raw_exprSet_final_batch_low['ADM',])
mean(gsva_Hypoxia[,phe_final_HIF1a[phe_final_HIF1a$hif1a_es==1,1]])
mean(gsva_Hypoxia[,phe_final_HIF1a[phe_final_HIF1a$hif1a_es==2,1]])

library(limma)
design <- model.matrix(~-1+factor(c(rep(1,100),rep(2,76))))

colnames(design) <- c("Low","High") # ??????design????????????rename the colnames of design

#output the design matrix

contrastmatrix <- makeContrasts(High-Low,levels=design)##and make the contrasts
contrastmatrix # check contrastmatrix


fit <- lmFit(BAL_batchmatrix_limma, design) # Run limma to get Differentially expressed genes
fit2 <- contrasts.fit(fit, contrastmatrix)  # Run limma to get Differentially expressed genes
?lmFit
#the borrowed variance approach described in class
fit2 <- eBayes(fit2) # ??????????????????FDR????????????
myresults_new <-topTable(fit2,coef=1, # 1 for FIH, 2 for VHL, 3 for VHL_FIH
                         adjust="fdr",
                         number=nrow(BAL_batchmatrix_limma))  # ??????????????????FDR????????????, extract final result after FDR adjustment
myresults_new$Symbol<-row.names(myresults_new) # ???????????????Symbol????????????create a new col 'Symbol' by the rowname of this data frame
for (i in 1:20189){
  myresults_new$Low_mean[i]<-mean(raw_exprSet_final_batch_low[rownames(myresults_new)[i],])
  myresults_new$High_mean[i]<-mean(raw_exprSet_final_batch_high[rownames(myresults_new)[i],])
}

myresults_new_up<-subset(myresults_new,adj.P.Val< 0.05)

myresults_new_down<-subset(myresults_new_up,myresults_new_up$logFC < 0)
myresults_new_up<-subset(myresults_new_up,myresults_new_up$logFC > 0)

save(myresults_new,myresults_new_up,myresults_new_down,file='BAL_HIF_high_low_limma.RData')
write.table(myresults_new,'BAL_HIF_high_low_limma.txt',sep = '\t')

write.table(myresults_new_up,'HIF_high_low_limma_UpDEG.txt',sep = '\t')
write.table(myresults_new_down,'HIF_high_low_limma_DownDEG.txt',sep = '\t')

myresults_new[c('NDRG1','ENO1','VEGFA','MRPS17','TPI1','CDKN3',
                'MIF','LDHA','ALDOA','TUBB6','PGAM1','SLC2A1','P4HA1','ACOT7','ADM'),]



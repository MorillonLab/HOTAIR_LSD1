#!/usr/bin/env Rscript

args<-commandArgs(TRUE)

home<-args[1]

library(ggplot2)
library(plyr)
library(reshape2)


#home<-"/media/marcgabriel/eda138bc-95c2-4778-99df-8915815cb86e/Marina_lsd1_test/"

#"/media/marcgabriel/eda138bc-95c2-4778-99df-8915815cb86e/Marina_lsd1_test"
chip_seq_results_G<-read.delim(paste(home,"G_TSS_with_closest_peak.tsv",sep=""),header=T)

chip_seq_results_H<-read.delim(paste(home,"H_TSS_with_closest_peak.tsv",sep=""),header=T)

chip_seq_results_L<-read.delim(paste(home,"L_TSS_with_closest_peak.tsv",sep=""),header=T)

chip_seq_results_P<-read.delim(paste(home,"P_TSS_with_closest_peak.tsv",sep=""),header=T)

G<-c(chip_seq_results_G$dist_peak_to_transcript_start)
G<-data.frame(dist_peak_to_transcript_start=G,peak_ID=1:length(G),condition="G")

H<-c(chip_seq_results_H$dist_peak_to_transcript_start)
H<-data.frame(dist_peak_to_transcript_start=H,peak_ID=1:length(H),condition="H")

P<-c(chip_seq_results_P$dist_peak_to_transcript_start)
P<-data.frame(dist_peak_to_transcript_start=P,peak_ID=1:length(P),condition="P")

L<-c(chip_seq_results_L$dist_peak_to_transcript_start)
L<-data.frame(dist_peak_to_transcript_start=L,peak_ID=1:length(L),condition="L")



all_conds<-rbind(G,H,P,L)

all_conds$dist_peak_to_transcript_start<-all_conds$dist_peak_to_transcript_start/1000

white_background<-theme(axis.line = element_line(colour = "black"),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank())

png(filename=paste(home,"boxplot_distrib_dist_peak_to_TSS_Minus5kb_Plus5kb_G_H_P_L_model1.png",sep=""),res=1700,width=87,height=87,units ="mm")


  p<-ggplot(all_conds,aes(x=condition,y=dist_peak_to_transcript_start,fill="transparent"))+
    stat_boxplot(geom ='errorbar',aes(color =condition),lwd=0.5,width= 0.5)+
    geom_boxplot(aes(color =condition),outlier.color = "transparent",outlier.size = 0.2,lwd=0.5)+
    geom_point(data =all_conds,aes(shape=NA,colour=condition))+
    #geom_point(colour = "black", size = 1)+
    geom_jitter(position=position_dodge2(width=0.5),color="black",size=0.01)+
    xlab("")+
    ylab("Distance from TSS, kb")+
    ggtitle(paste("limits : -5kb/+5kb around TSS\n","nb TSS with G peaks : ",nrow(G),"\n","nb TSS with L peaks : ",nrow(L),"\n","nb TSS with H peaks : ",nrow(H),"\n","nb TSS with P peaks : ",nrow(P),sep=""))+
    #geom_point(data =all_conds,aes(colour=condition))+
    scale_fill_manual(values=c("G"="transparent","L"="transparent","H"="transparent","P"="transparent"),name="condition",guide=FALSE)+
    scale_color_manual(values=c("G"="#7f7f7f","L"="#e46c0a","H"="#604a7b","P"="#31859c"),name="condition",guide=FALSE)
  
  
  dat <- ggplot_build(p)$data[[1]]
  
  meds <- ddply(all_conds, .(condition,condition), summarize, med = median(dist_peak_to_transcript_start))
  
  #remark med & middle are the same..., so no need to make all this stuff if its work
  dat<-cbind(dat,meds)
  
  
  #dat <- ggplot_build(p)$data[[1]]
  
  #dat<-data.frame(dat,condition=unique(all_conds$condition))
  
  dat <- ggplot_build(p)$data[[1]]
  
  meds <- ddply(all_conds, .(condition), summarize, med = median(dist_peak_to_transcript_start))
  
  dat<-cbind(dat,meds)
  
  my_theme<-theme(axis.text.x=element_text(color="black",size=12,hjust=0.5,face="bold"),axis.text.y = element_text(color="black",size=12,face="bold"),plot.title = element_text(hjust = 0.5,size=12),axis.title=element_text(size=12,face="bold"),legend.text=element_text(size=12),legend.title=element_text(size=12,face="bold"))
  
  
  
  print(p +white_background+my_theme
        
  )

dev.off() 

cairo_ps(filename=paste(home,"boxplot_distrib_dist_peak_to_TSS_Minus5kb_Plus5kb_G_H_P_L_model1.eps",sep=""),width=3.4252,height=3.4252,fallback_resolution = 1700)

 print(p +white_background+my_theme)

dev.off()




png(filename=paste(home,"boxplot_distrib_dist_peak_to_TSS_Minus5kb_Plus5kb_G_H_P_L_model2.png",sep=""),res=1700,width=87,height=87,units ="mm")


p<-ggplot(all_conds,aes(x=condition,y=dist_peak_to_transcript_start,fill=condition))+
  stat_boxplot(geom ='errorbar',aes(color =condition),width= 0.5,lwd=0.5)+
  geom_boxplot(aes(color =condition),outlier.color = "black",outlier.size = 0.2)+
  geom_point(data =all_conds,aes(shape=NA,colour=condition))+
  #geom_point(colour = "black", size = 1)+
  xlab("")+
  ylab("Distance from TSS, kb")+
  ggtitle(paste("limits : -5kb/+5kb around TSS\n","nb TSS with G peaks : ",nrow(G),"\n","nb TSS with L peaks : ",nrow(L),"\n","nb TSS with H peaks : ",nrow(H),"\n","nb TSS with P peaks : ",nrow(P),sep=""))+
  #geom_point(data =all_conds,aes(colour=condition))+
  scale_fill_manual(values=c("G"="#7f7f7f","L"="#e46c0a","H"="#604a7b","P"="#31859c"),name="condition",guide=FALSE)+
  scale_color_manual(values=c("G"="#7f7f7f","L"="#e46c0a","H"="#604a7b","P"="#31859c"),name="condition",guide=FALSE)


dat <- ggplot_build(p)$data[[1]]

meds <- ddply(all_conds, .(condition,condition), summarize, med = median(dist_peak_to_transcript_start))

#remark med & middle are the same..., so no need to make all this stuff if its work
dat<-cbind(dat,meds)


#dat <- ggplot_build(p)$data[[1]]

#dat<-data.frame(dat,condition=unique(all_conds$condition))

dat <- ggplot_build(p)$data[[1]]

meds <- ddply(all_conds, .(condition), summarize, med = median(dist_peak_to_transcript_start))

dat<-cbind(dat,meds)

my_theme<-theme(axis.text.x=element_text(color="black",size=12,hjust=0.5,face="bold"),axis.text.y = element_text(color="black",size=12,face="bold"),plot.title = element_text(hjust = 0.5,size=12),axis.title=element_text(size=12,face="bold"),legend.text=element_text(size=12),legend.title=element_text(size=12,face="bold"))


print(p + geom_segment(data=dat, aes(x=xmin, xend=xmax, 
                                     y=med, yend=med), colour="black", size=0.5)+white_background+my_theme
      
)


dev.off()

cairo_ps(filename=paste(home,"boxplot_distrib_dist_peak_to_TSS_Minus5kb_Plus5kb_G_H_P_L_model2.eps",sep=""),width=3.4252,height=3.4252,fallback_resolution = 1700)

 print(p + geom_segment(data=dat, aes(x=xmin, xend=xmax, 
                                     y=med, yend=med), colour="black", size=0.5)+white_background+my_theme)

dev.off()

my_table<-data.frame()

my_comp<-data.frame(t(combn(unique(as.character(all_conds$condition)),2)))
names(my_comp)<-c("cond1","cond2")
for(i in 1:nrow(my_comp)){
  
  cond1<-as.character(my_comp$cond1[i])
  
  cond2<-as.character(my_comp$cond2[i])
  
  
  my_pvalue<-wilcox.test(all_conds[which(all_conds$condition==cond1),]$dist_peak_to_transcript_start,all_conds[which(all_conds$condition==cond2),]$dist_peak_to_transcript_start)$p.value
  
  my_table<-rbind(my_table,
                  data.frame(cond1=cond1,
                             cond2=cond2,
                             nb_cond1=length(all_conds[which(all_conds$condition==cond1),]$dist_peak_to_transcript_start),
                             nb_cond2=length(all_conds[which(all_conds$condition==cond2),]$dist_peak_to_transcript_start),
                             mean_dist_cond1=mean(all_conds[which(all_conds$condition==cond1),]$dist_peak_to_transcript_start),
                             mean_dist_cond2=mean(all_conds[which(all_conds$condition==cond2),]$dist_peak_to_transcript_start),
                             median_dist_cond1=median(all_conds[which(all_conds$condition==cond1),]$dist_peak_to_transcript_start),
                             median_dist_cond2=median(all_conds[which(all_conds$condition==cond2),]$dist_peak_to_transcript_start),
                             wilcoxon_test_pvalue=my_pvalue))
  
  
}

write.table(my_table,paste(home,"boxplot_distrib_dist_peak_to_TSS_Minus5kb_Plus5kb_G_H_P_L_table.tsv",sep=""),sep="\t",col.names = T,row.names = F,quote=F)


############################

############################


#"/media/marcgabriel/eda138bc-95c2-4778-99df-8915815cb86e/Marina_lsd1_test"
chip_seq_results_G<-read.delim(paste(home,"G_allPeaks_with_closest_feat.bed",sep=""),header=F)

chip_seq_results_L<-read.delim(paste(home,"L_allPeaks_with_closest_feat.bed",sep=""),header=F)

chip_seq_results_H<-read.delim(paste(home,"H_allPeaks_with_closest_feat.bed",sep=""),header=F)

chip_seq_results_P<-read.delim(paste(home,"P_allPeaks_with_closest_feat.bed",sep=""),header=F)


G<-c(chip_seq_results_G$V19)
G<-data.frame(dist_peak_to_transcript_start=G,peak_ID=1:length(G),condition="G")

H<-c(chip_seq_results_H$V19)
H<-data.frame(dist_peak_to_transcript_start=H,peak_ID=1:length(H),condition="H")

P<-c(chip_seq_results_P$V19)
P<-data.frame(dist_peak_to_transcript_start=P,peak_ID=1:length(P),condition="P")

L<-c(chip_seq_results_L$V19)
L<-data.frame(dist_peak_to_transcript_start=L,peak_ID=1:length(L),condition="L")





all_conds<-rbind(G,H,P,L)

all_conds$dist_peak_to_transcript_start<-all_conds$dist_peak_to_transcript_start/1000

white_background<-theme(axis.line = element_line(colour = "black"),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank())

png(filename=paste(home,"boxplot_distrib_dist_peak_to_TSS_alldist_G_H_P_L_model1.png",sep=""),res=1700,width=87,height=87,units ="mm")

  
  p<-ggplot(all_conds,aes(x=condition,y=dist_peak_to_transcript_start,fill="transparent"))+
    stat_boxplot(geom ='errorbar',aes(color =condition),lwd=0.5,width= 0.5)+
    geom_boxplot(aes(color =condition),outlier.color = "transparent",outlier.size = 0.2,lwd=0.5)+
    geom_point(data =all_conds,aes(shape=NA,colour=condition))+
    #geom_point(colour = "black", size = 1)+
    geom_jitter(position=position_dodge2(width=0.5),color="black",size=0.01)+
    xlab("")+
    ylab("Distance from TSS, kb")+
    ggtitle(paste("limits : no dist limits\n","nb G peaks : ",nrow(G),"\n","nb L peaks : ",nrow(L),"\n","nb H peaks : ",nrow(H),"\n","nb P peaks : ",nrow(P),sep=""))+
    #geom_point(data =all_conds,aes(colour=condition))+
    scale_fill_manual(values=c("G"="transparent","L"="transparent","H"="transparent","P"="transparent"),name="condition",guide=FALSE)+
    scale_color_manual(values=c("G"="#7f7f7f","L"="#e46c0a","H"="#604a7b","P"="#31859c"),name="condition",guide=FALSE)
  
  dat <- ggplot_build(p)$data[[1]]
  
  meds <- ddply(all_conds, .(condition,condition), summarize, med = median(dist_peak_to_transcript_start))
  
  #remark med & middle are the same..., so no need to make all this stuff if its work
  dat<-cbind(dat,meds)
  
  
  #dat <- ggplot_build(p)$data[[1]]
  
  #dat<-data.frame(dat,condition=unique(all_conds$condition))
  
  dat <- ggplot_build(p)$data[[1]]
  
  meds <- ddply(all_conds, .(condition), summarize, med = median(dist_peak_to_transcript_start))
  
  dat<-cbind(dat,meds)
  
  my_theme<-theme(axis.text.x=element_text(color="black",size=12,hjust=0.5,face="bold"),axis.text.y = element_text(color="black",size=12,face="bold"),plot.title = element_text(hjust = 0.5,size=12),axis.title=element_text(size=12,face="bold"),legend.text=element_text(size=12),legend.title=element_text(size=12,face="bold"))
  
  
  
  print(p +white_background+my_theme)

dev.off()

cairo_ps(filename=paste(home,"boxplot_distrib_dist_peak_to_TSS_alldist_G_H_P_L_model1.eps",sep=""),width=3.4252,height=3.4252,fallback_resolution = 1700)

print(p +white_background+my_theme)

dev.off()




png(filename=paste(home,"boxplot_distrib_dist_peak_to_TSS_alldist_G_H_P_L_model2.png",sep=""),res=1700,width=87,height=87,units ="mm")


p<-ggplot(all_conds,aes(x=condition,y=dist_peak_to_transcript_start,fill=condition))+
  stat_boxplot(geom ='errorbar',aes(color =condition),width = 0.5,width= 0.5)+
  geom_boxplot(aes(color =condition),outlier.color = "black",outlier.size = 0.2)+
  geom_point(data =all_conds,aes(shape=NA,colour=condition))+
  #geom_point(colour = "black", size = 1)+
  xlab("")+
  ylab("Distance from TSS, kb")+
  ggtitle(paste("limits : no dist limits\n","nb G peaks : ",nrow(G),"\n","nb L peaks : ",nrow(L),"\n","nb H peaks : ",nrow(H),"\n","nb P peaks : ",nrow(P),sep=""))+
  #geom_point(data =all_conds,aes(colour=condition))+
  scale_fill_manual(values=c("G"="#7f7f7f","L"="#e46c0a","H"="#604a7b","P"="#31859c"),name="condition",guide=FALSE)+
  scale_color_manual(values=c("G"="#7f7f7f","L"="#e46c0a","H"="#604a7b","P"="#31859c"),name="condition",guide=FALSE)


dat <- ggplot_build(p)$data[[1]]

meds <- ddply(all_conds, .(condition,condition), summarize, med = median(dist_peak_to_transcript_start))

#remark med & middle are the same..., so no need to make all this stuff if its work
dat<-cbind(dat,meds)


#dat <- ggplot_build(p)$data[[1]]

#dat<-data.frame(dat,condition=unique(all_conds$condition))

dat <- ggplot_build(p)$data[[1]]

meds <- ddply(all_conds, .(condition), summarize, med = median(dist_peak_to_transcript_start))

dat<-cbind(dat,meds)

my_theme<-theme(axis.text.x=element_text(color="black",size=12,hjust=0.5,face="bold"),axis.text.y = element_text(color="black",size=12,face="bold"),plot.title = element_text(hjust = 0.5,size=12),axis.title=element_text(size=12,face="bold"),legend.text=element_text(size=12),legend.title=element_text(size=12,face="bold"))


print(p + geom_segment(data=dat, aes(x=xmin, xend=xmax,
                                     y=med, yend=med), colour="black", size=0.5)+white_background+my_theme

)

dev.off()


cairo_ps(filename=paste(home,"boxplot_distrib_dist_peak_to_TSS_alldist_G_H_P_L_model2.eps",sep=""),width=3.4252,height=3.4252,fallback_resolution = 1700)

 print(p + geom_segment(data=dat, aes(x=xmin, xend=xmax,
                                     y=med, yend=med), colour="black", size=0.5)+white_background+my_theme
      
)
 
dev.off()

my_table<-data.frame()

my_comp<-data.frame(t(combn(unique(as.character(all_conds$condition)),2)))
names(my_comp)<-c("cond1","cond2")
for(i in 1:nrow(my_comp)){

  cond1<-as.character(my_comp$cond1[i])

  cond2<-as.character(my_comp$cond2[i])


  my_pvalue<-wilcox.test(all_conds[which(all_conds$condition==cond1),]$dist_peak_to_transcript_start,all_conds[which(all_conds$condition==cond2),]$dist_peak_to_transcript_start)$p.value

  my_table<-rbind(my_table,
                  data.frame(cond1=cond1,
                             cond2=cond2,
                             nb_cond1=length(all_conds[which(all_conds$condition==cond1),]$dist_peak_to_transcript_start),
                             nb_cond2=length(all_conds[which(all_conds$condition==cond2),]$dist_peak_to_transcript_start),
                             mean_dist_cond1=mean(all_conds[which(all_conds$condition==cond1),]$dist_peak_to_transcript_start),
                             mean_dist_cond2=mean(all_conds[which(all_conds$condition==cond2),]$dist_peak_to_transcript_start),
                             median_dist_cond1=median(all_conds[which(all_conds$condition==cond1),]$dist_peak_to_transcript_start),
                             median_dist_cond2=median(all_conds[which(all_conds$condition==cond2),]$dist_peak_to_transcript_start),
                             wilcoxon_test_pvalue=my_pvalue))


}

write.table(my_table,paste(home,"boxplot_distrib_dist_peak_to_TSS_alldist_G_H_P_L_table.tsv",sep=""),sep="\t",col.names = T,row.names = F,quote=F)



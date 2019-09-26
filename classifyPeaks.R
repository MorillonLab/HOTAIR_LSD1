#!/usr/bin/env Rscript

args <- commandArgs(TRUE)


suppressMessages(library(ggplot2))

suppressMessages(library(plotrix))

suppressMessages(library(Biostrings))

suppressMessages(library(data.table))

suppressMessages(library(RColorBrewer))

#peak attribution rules :
#exon > 5'UTR > 3'UTR > intron > promoter > intergenic


#args[1]

home<-args[1]

file<-args[2]

#file<-"/media/marcgabriel/39f07759-1445-4caf-bf52-06e8bb75d6f3/Marina_HOTAIR_lsd1/chip_seq_analysis/EPI_MES_chip_seq/all_H3K4me3_counts.tsv"
#file="/media/marcgabriel/eda138bc-95c2-4778-99df-8915815cb86e/Marina_lsd1_test/all_P_counts.tsv"

all_peaks_OneRep<-fread(file,header=T,data.table=F,sep="\t")


all_read_number<-nrow(all_peaks_OneRep)

# &
#all_peaks_OneRep$TSS==F
intronic_peaks<-all_peaks_OneRep[which(all_peaks_OneRep$intron==T &
                                         all_peaks_OneRep$exon==F &
                                         all_peaks_OneRep$five_prime_utr==F &
                                         all_peaks_OneRep$three_prime_utr==F),]

intronic<-nrow(intronic_peaks)

#all_peaks_OneRep$TSS==F &
promoter_peaks<-all_peaks_OneRep[which(all_peaks_OneRep$promoter==T &
                                       all_peaks_OneRep$five_prime_utr==F &
                                       all_peaks_OneRep$exon==F,
                                       all_peaks_OneRep$intron==F),]

promoter<-nrow(promoter_peaks)

exonic_peaks<-all_peaks_OneRep[which(all_peaks_OneRep$exon==T),]

exonic<-nrow(exonic_peaks)


intergenic_peaks<-all_peaks_OneRep[which(all_peaks_OneRep$intergenic==T &
                                           all_peaks_OneRep$intron==F &
                                           all_peaks_OneRep$promoter==F &
                                           all_peaks_OneRep$exon==F &
                                           all_peaks_OneRep$five_prime_utr==F &
                                           all_peaks_OneRep$three_prime_utr==F),]

intergenic<-nrow(intergenic_peaks)

five_prime_utr_peaks<-all_peaks_OneRep[which(all_peaks_OneRep$five_prime_utr==T & all_peaks_OneRep$exon==F),]

five_prime_utr<-nrow(five_prime_utr_peaks)

three_prime_utr_peaks<-all_peaks_OneRep[which(all_peaks_OneRep$three_prime_utr==T & all_peaks_OneRep$exon==F &
                                              all_peaks_OneRep$five_prime_utr==F & all_peaks_OneRep$promoter==F),]

three_prime_utr<-nrow(three_prime_utr_peaks)

all_colors<-rainbow(6)

mydata.frame<-data.frame(feature=c("intronic","exonic","intergenic","promoter","5' UTR","3' UTR"),
           number=c(intronic,exonic,intergenic,promoter,five_prime_utr,three_prime_utr),
           percentage=round((c(intronic,exonic,intergenic,promoter,five_prime_utr,three_prime_utr)/all_read_number)*100,2),
           label=paste(c("intronic","exonic","intergenic","promoter (<=1kb from TSS)","5' UTR","3' UTR"),
                       paste(round((c(intronic,exonic,intergenic,promoter,five_prime_utr,three_prime_utr)/all_read_number)*100,2),"%",sep=""),sep="\n"))


png(filename=paste(home,unique(all_peaks_OneRep$mark_type),"_peak_distrib.png",sep=""),width=1200,height=800)


pie(c(intronic,exonic,intergenic,promoter,five_prime_utr,three_prime_utr),mydata.frame$label,col=all_colors,
    main=paste(unique(all_peaks_OneRep$mark_type)," peak distribution across genomic features\nTotal peaks = ",format(all_read_number,big.mark=" "),sep=""))

#legend('topright', legend=paste(as.character(classified_reads$features)," ",classified_reads$percentage,"%",sep=""), fill=as.character(classified_reads$color))

dev.off()

write.table(file=paste(home,unique(all_peaks_OneRep$mark_type),"_peak_genomic_distrib.tsv",sep=""),mydata.frame[,1:(ncol(mydata.frame)-1)],row.names=F,col.names=T,quote=F,sep="\t")


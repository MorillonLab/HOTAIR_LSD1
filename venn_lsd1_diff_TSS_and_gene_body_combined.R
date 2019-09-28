#!/usr/bin/env Rscript

library(gdata)

library(eulerr)

library(lattice)

library(latticeExtra)

library(VennDiagram)

library(xlsx)
library(pdp)

args <- commandArgs(TRUE)

if (length(args)==0){
  stop("missing arguments !\nUsage :\n
       ./venn_lsd1_diff_TSS_and_gene_body_combined.R home1 home2 home3 official_annotation\n
       
       home1 : output directory for this script\n
       home2 : output directory of the script \"intersection_GPLH_Lsd1_peaks_rep_not_merged_genebody.sh\"\n
       home3 : output directory of the script \"intersection_GPLH_Lsd1_peaks_rep_not_merged.sh\"\n\n
       remarks : these 3 directories should be different!")
}

#find location of this script
LocationOfThisScript = function() # Function LocationOfThisScript returns the location of this .R script (may be needed to source other files in same dir)
{
  this.file = NULL
  # This file may be 'sourced'
  for (i in -(1:sys.nframe())) {
    if (identical(sys.function(i), base::source)) this.file = (normalizePath(sys.frame(i)$ofile))
  }
  
  if (!is.null(this.file)) return(dirname(this.file))
  
  # But it may also be called from the command line
  cmd.args = commandArgs(trailingOnly = FALSE)
  cmd.args.trailing = commandArgs(trailingOnly = TRUE)
  cmd.args = cmd.args[seq.int(from=1, length.out=length(cmd.args) - length(cmd.args.trailing))]
  res = gsub("^(?:--file=(.*)|.*)$", "\\1", cmd.args)
  
  # If multiple --file arguments are given, R uses the last one
  res = tail(res[res != ""], 1)
  if (0 < length(res)) return(dirname(res))
  
  # Both are not the case. Maybe we are in an R GUI?
  return(NULL)
}

current.dir = paste(LocationOfThisScript(),"/",sep="")

#script that allows recursive intersections
source(paste(current.dir,"getAllCombinations.R",sep=""))

# home1<-"/media/marcgabriel/Transcend/LSD1_metagenes/Marina_peaks_analysis_TSS_and_genebody/"
# home2<-"/media/marcgabriel/Transcend/LSD1_metagenes/Marina_peaks_analysis_rep_not_merged_genebody/"
# home3<-"/media/marcgabriel/Transcend/LSD1_metagenes/Marina_peaks_analysis_rep_not_merged/"
# official_annotation<-"/home/marcgabriel/Documents/gencode26lift37/gencode.v26lift37.annotation.sorted.gff3"
# gene_type<-"/home/marcgabriel/Documents/gencode26lift37/matching_gene_name_types.tsv"

home1<-args[1]
home2<-args[2]
home3<-args[3]
official_annotation<-args[4]

gene_type<-paste(current.dir,matching_gene_name_types.tsv,sep="")

#check if all suplied directories are different
if(all(c(home1,home2,home3)==home1)==T |
   length(grep("TRUE",c(home1,home2,home3)==home1))>1|
   length(grep("TRUE",c(home1,home2,home3)==home2))>1|
   length(grep("TRUE",c(home1,home2,home3)==home3))>1
   
   ){
  
  stop("you should have 3 different directories !!")
}

#check if the annotation file exists
if(length(official_annotation)==0 |file.exists(official_annotation)==F){
  
  stop("your annotation file is missing, or doesn't exit !!")
}


#process the supplied different directories
home1<-paste(home,"/",sep="")
home1<-gsub("//","/",home)

home2<-paste(home,"/",sep="")
home2<-gsub("//","/",home)

home3<-paste(home3,"/",sep="")
home3<-gsub("//","/",home3)

#create the output dir of this script
dir.create(home1,showWarnings=F,recursive=T)





#gene_body<-"/media/marcgabriel/Transcend/LSD1_metagenes/Marina_peaks_analysis_rep_not_merged_genebody/intersection_genesWithPeaks_Minus0kb_Plus0kb_G_H_P_L.tsv"
gene_body<-paste(home2,"/intersection_genesWithPeaks_Minus0kb_Plus0kb_G_H_P_L.tsv",sep="")

#TSS<-"/media/marcgabriel/Transcend/LSD1_metagenes/Marina_peaks_analysis_rep_not_merged/intersection_genesWithPeaks_Minus5kb_Plus5kb_G_H_P_L.tsv"
TSS<-paste(home3,"/intersection_genesWithPeaks_Minus5kb_Plus5kb_G_H_P_L.tsv",sep="")

setwd(home1)

gene_type<-read.delim(gene_type,sep="\t",check.names=F,header=T)

official_annotation<-read.delim(file=pipe(paste("grep -v ^# ",official_annotation," |grep -P \"\\tgene\\t\"",sep="")),
                                sep="\t",check.names=F,header=F)


names(official_annotation) <-c("chromosome","source","feature","start","end","score","strand","frame","description")
official_annotation$gene_id<-gsub("gene_id=","",grep("gene_id=",unlist(strsplit(as.character(official_annotation$description),";")),value=T))

#extract gene name (HUGO ID)
official_annotation$gene_name<-gsub("gene_name=","",grep("gene_name=",unlist(strsplit(as.character(official_annotation$description),";")),value=T))

official_annotation$gene_type<-gsub("gene_type=","",grep("gene_type",unlist(strsplit(as.character(official_annotation$description),";")),value=T))

official_annotation<-official_annotation[,c("gene_name","gene_type","chromosome","start","end","strand")]


gene_body<-read.delim(paste(gene_body,sep=""),header=T)
TSS<-read.delim(paste(TSS,sep=""),header=T)

gene_body_IDs<-gene_body[,grep("IDs",names(gene_body),value=T)]

TSS_IDs<-TSS[,grep("IDs",names(TSS),value=T)]

G_body<-c()
L_body<-c()
H_body<-c()
P_body<-c()
for(i in 1:ncol(gene_body_IDs)){
  
  tmp<-unlist(strsplit(as.character(unlist(gene_body_IDs[i]))," "))
  tmp<-tmp[!is.na(tmp)]
  tmp<-tmp[!tmp%in%"NA"]
  
  if(grepl("cond1",names(gene_body_IDs[i]))==T){
    
  G_body<-c(G_body,tmp)
  
  }
  
  if(grepl("cond2",names(gene_body_IDs[i]))==T){
    
    L_body<-c(L_body,tmp)
    
  }
  if(grepl("cond3",names(gene_body_IDs[i]))==T){
    
    H_body<-c(H_body,tmp)
    
  }
  
  if(grepl("cond4",names(gene_body_IDs[i]))==T){
    
    P_body<-c(P_body,tmp)
    
  }
  
}

G_TSS<-c()
L_TSS<-c()
H_TSS<-c()
P_TSS<-c()
for(i in 1:ncol(TSS_IDs)){
  
  tmp<-unlist(strsplit(as.character(unlist(TSS_IDs[i]))," "))
  tmp<-tmp[!is.na(tmp)]
  tmp<-tmp[!tmp%in%"NA"]
  
  if(grepl("cond1",names(TSS_IDs[i]))==T){
    
    G_TSS<-c(G_TSS,tmp)
    
  }
  
  if(grepl("cond2",names(TSS_IDs[i]))==T){
    
    L_TSS<-c(L_TSS,tmp)
    
  }
  if(grepl("cond3",names(TSS_IDs[i]))==T){
    
    H_TSS<-c(H_TSS,tmp)
    
  }
  
  if(grepl("cond4",names(TSS_IDs[i]))==T){
    
    P_TSS<-c(P_TSS,tmp)
    
  }
  
}
G_both<-unique(c(G_TSS,G_body))
L_both<-unique(c(L_TSS,L_body))
H_both<-unique(c(H_TSS,H_body))
P_both<-unique(c(P_TSS,P_body))

result3<-getAllConbinations(list(G_both,L_both,H_both,P_both),list("G","L","H","P"))

intersected_all<-intersect(intersect(intersect(G_both,H_both),P_both),L_both)

exp_names2<-c(paste("cond1\nG","\n(n = ",format(length(G_both),big.mark = " "),")",sep=""),
              paste("cond2\nL","\n(n = ",format(length(L_both),big.mark = " "),")",sep=""),
              paste("cond3\nH","\n(n = ",format(length(H_both),big.mark = " "),")",sep=""),
              paste("cond4\nP","\n(n = ",format(length(P_both),big.mark = " "),")",sep=""))

all_colors<-c("#7f7f7f","#e46c0a","#604a7b","#31859c")

#G,H,P,L
#area1,area3,area4,area2
venn.plot <- draw.quad.venn(
  area1 = length(G_both),
  area2 = length(L_both),
  area3 =length(H_both),
  area4 = length(P_both),
  n12 = length(intersect(G_both,L_both)),
  n13 = length(intersect(G_both,H_both)),
  n14 = length(intersect(G_both,P_both)),
  n23 = length(intersect(L_both,H_both)),
  n24 = length(intersect(L_both,P_both)),
  n34 = length(intersect(H_both,P_both)),
  n123 = length(intersect(intersect(G_both,L_both),H_both)),
  n124 = length(intersect(intersect(G_both,L_both),P_both)),
  n134 = length(intersect(intersect(G_both,H_both),P_both)),
  n234 = length(intersect(intersect(L_both,H_both),P_both)),
  n1234 = length(intersected_all),
  category = exp_names2,
  fill = all_colors,
  cex = 1,
  cat.cex = 0.8,
  margin = 0.04,
  alpha=c(0.7,0.7,0.7,0.7),
  cat.dist=c(0.3,0.3,0.2,0.2))



png(filename=paste(home1,"intersection_genesWithPeaks_TSS_genebody_G_H_P_L.png",sep=""),res=1700,width=87,height=87,units ="mm")

#grid.draw(venn.plot)
grid.arrange(gTree(children=venn.plot),top="Genes with peaks")

dev.off()

cairo_ps(filename=paste(home1,"intersection_genesWithPeaks_TSS_genebody_G_H_P_L.eps",sep=""),width=3.4252,height=3.4252,fallback_resolution = 1700)

grid.arrange(gTree(children=venn.plot),top="Genes with peaks")

dev.off()

venn.plot <- draw.quad.venn(
  area1 = length(G_both),
  area2 = length(L_both),
  area3 =length(H_both),
  area4 = length(P_both),
  n12 = length(intersect(G_both,L_both)),
  n13 = length(intersect(G_both,H_both)),
  n14 = length(intersect(G_both,P_both)),
  n23 = length(intersect(L_both,H_both)),
  n24 = length(intersect(L_both,P_both)),
  n34 = length(intersect(H_both,P_both)),
  n123 = length(intersect(intersect(G_both,L_both),H_both)),
  n124 = length(intersect(intersect(G_both,L_both),P_both)),
  n134 = length(intersect(intersect(G_both,H_both),P_both)),
  n234 = length(intersect(intersect(L_both,H_both),P_both)),
  n1234 = length(intersected_all),
  category = exp_names2,
  fill = all_colors,
  cex = 1,
  cat.cex = 0.8,
  label.col = rep("transparent", 15),
  margin = 0.04,
  alpha=c(0.7,0.7,0.7,0.7),
  cat.dist=c(0.3,0.3,0.2,0.2))



png(filename=paste(home1,"intersection_genesWithPeaks_TSS_genebody_G_H_P_L_NoLabels.png",sep=""),res=1700,width=87,height=87,units ="mm")

#grid.draw(venn.plot)
grid.arrange(gTree(children=venn.plot),top="Genes with peaks")

dev.off()

cairo_ps(filename=paste(home1,"intersection_genesWithPeaks_TSS_genebody_G_H_P_L_NoLabels.eps",sep=""),width=3.4252,height=3.4252,fallback_resolution = 1700)

grid.arrange(gTree(children=venn.plot),top="Genes with peaks")

dev.off()

setwd(home2)

G_genes_with_rep_peaks<-read.delim(paste(home2,"G_genes_with_peaks_in_both_rep.tsv",sep=""),header=T)
G_genes_with_rep_peaks_back<-G_genes_with_rep_peaks
G_genes_with_rep_peaks<-unique(as.character(G_genes_with_rep_peaks$gene_name))


H_genes_with_rep_peaks<-read.delim(paste(home2,"H_genes_with_peaks_in_both_rep.tsv",sep=""),header=T)
H_genes_with_rep_peaks<-unique(as.character(H_genes_with_rep_peaks$gene_name))

P_genes_with_rep_peaks<-read.delim(paste(home2,"P_genes_with_peaks_in_both_rep.tsv",sep=""),header=T)
P_genes_with_rep_peaks<-unique(as.character(P_genes_with_rep_peaks$gene_name))

L_genes_with_rep_peaks<-read.delim(paste(home2,"L_genes_with_peaks_in_both_rep.tsv",sep=""),header=T)
L_genes_with_rep_peaks_back<-L_genes_with_rep_peaks
L_genes_with_rep_peaks<-unique(as.character(L_genes_with_rep_peaks$gene_name))

all_gene_body<-c(G_genes_with_rep_peaks,H_genes_with_rep_peaks,P_genes_with_rep_peaks,L_genes_with_rep_peaks)


setwd(home3)

G_genes_with_rep_peaks2<-read.delim(paste(home3,"G_genes_with_peaks_in_both_rep.tsv",sep=""),header=T)
G_genes_with_rep_peaks2_back<-G_genes_with_rep_peaks2
G_genes_with_rep_peaks2<-unique(as.character(G_genes_with_rep_peaks2$gene_name))


H_genes_with_rep_peaks2<-read.delim(paste(home3,"H_genes_with_peaks_in_both_rep.tsv",sep=""),header=T)
H_genes_with_rep_peaks2<-unique(as.character(H_genes_with_rep_peaks2$gene_name))

P_genes_with_rep_peaks2<-read.delim(paste(home3,"P_genes_with_peaks_in_both_rep.tsv",sep=""),header=T)
P_genes_with_rep_peaks2<-unique(as.character(P_genes_with_rep_peaks2$gene_name))

L_genes_with_rep_peaks2<-read.delim(paste(home3,"L_genes_with_peaks_in_both_rep.tsv",sep=""),header=T)
L_genes_with_rep_peaks2_back<-L_genes_with_rep_peaks2
L_genes_with_rep_peaks2<-unique(as.character(L_genes_with_rep_peaks2$gene_name))


all_TSS<-c(G_genes_with_rep_peaks2,H_genes_with_rep_peaks2,P_genes_with_rep_peaks2,L_genes_with_rep_peaks2)


genes_withclosest_GL_peaks<-unlist(strsplit(result3$cond1_cond2_IDs," "))

genes_withclosest_GL_peaks_frame<-data.frame(genes_withclosest_GL_peaks)

#when "NA", check the other table (priority to tss, then genebody)
genes_withclosest_GL_peaks_frame<-merge(genes_withclosest_GL_peaks_frame,official_annotation,
                                        by.x="genes_withclosest_GL_peaks",
                                        by.y="gene_name",
                                        all.x=TRUE,
                                        all.y=FALSE)

genes_withclosest_GL_peaks_frame<-merge(genes_withclosest_GL_peaks_frame[,c("genes_withclosest_GL_peaks","chromosome","strand","gene_type")]
                                        ,G_genes_with_rep_peaks_back[,c("transcript_start","gene_name","dist_peak_to_transcript_loc","chromosome")],
                                        by.x=c("genes_withclosest_GL_peaks","chromosome"),
                                        by.y=c("gene_name","chromosome"),
                                        all.x=TRUE,
                                        all.y=FALSE)

missing_genes<-unique(as.character(genes_withclosest_GL_peaks_frame[which(is.na(genes_withclosest_GL_peaks_frame$transcript_start)),]$genes_withclosest_GL_peaks))

for(i in 1:length(missing_genes)){
  
  G_genes_with_rep_peaks2_back[which(G_genes_with_rep_peaks2_back$gene_name==missing_genes[i]),]
  
  genes_withclosest_GL_peaks_frame$transcript_start[genes_withclosest_GL_peaks_frame$genes_withclosest_GL_peaks==missing_genes[i]]<-G_genes_with_rep_peaks2_back[which(G_genes_with_rep_peaks2_back$gene_name==missing_genes[i]),]$transcript_start[1]
  
  genes_withclosest_GL_peaks_frame$dist_peak_to_transcript_loc[genes_withclosest_GL_peaks_frame$genes_withclosest_GL_peaks==missing_genes[i]]<-G_genes_with_rep_peaks2_back[which(G_genes_with_rep_peaks2_back$gene_name==missing_genes[i]),]$dist_peak_to_transcript_loc[1]
  
  
}

#genes_withclosest_GL_peaks_frame<-genes_withclosest_GL_peaks_frame[which(!is.na(genes_withclosest_GL_peaks_frame$dist_peak_to_transcript_loc)),]

#G_genes_with_rep_peaks2_back[which(G_genes_with_rep_peaks2_back$gene_name%in%missing_genes),]





colnames(genes_withclosest_GL_peaks_frame)[colnames(genes_withclosest_GL_peaks_frame)=="dist_peak_to_transcript_loc"] <- "G_dist_peak_to_transcript_loc"

colnames(genes_withclosest_GL_peaks_frame)[colnames(genes_withclosest_GL_peaks_frame)=="transcript_start"] <- "G_transcript_start"

genes_withclosest_GL_peaks_frame<-merge(genes_withclosest_GL_peaks_frame
                                        ,L_genes_with_rep_peaks_back[,c("gene_name","transcript_start","dist_peak_to_transcript_loc","chromosome")],
                                        by.x=c("genes_withclosest_GL_peaks","chromosome"),
                                        by.y=c("gene_name","chromosome"),
                                        all.x=TRUE,
                                        all.y=FALSE)

missing_genes<-unique(as.character(genes_withclosest_GL_peaks_frame[which(is.na(genes_withclosest_GL_peaks_frame$transcript_start)),]$genes_withclosest_GL_peaks))

for(i in 1:length(missing_genes)){
  
  L_genes_with_rep_peaks2_back[which(L_genes_with_rep_peaks2_back$gene_name==missing_genes[i]),]
  
  genes_withclosest_GL_peaks_frame$transcript_start[genes_withclosest_GL_peaks_frame$genes_withclosest_GL_peaks==missing_genes[i]]<-L_genes_with_rep_peaks2_back[which(L_genes_with_rep_peaks2_back$gene_name==missing_genes[i]),]$transcript_start[1]
  
  genes_withclosest_GL_peaks_frame$dist_peak_to_transcript_loc[genes_withclosest_GL_peaks_frame$genes_withclosest_GL_peaks==missing_genes[i]]<-L_genes_with_rep_peaks2_back[which(L_genes_with_rep_peaks2_back$gene_name==missing_genes[i]),]$dist_peak_to_transcript_loc[1]
  
  
}

colnames(genes_withclosest_GL_peaks_frame)[colnames(genes_withclosest_GL_peaks_frame)=="dist_peak_to_transcript_loc"] <- "L_dist_peak_to_transcript_loc"
colnames(genes_withclosest_GL_peaks_frame)[colnames(genes_withclosest_GL_peaks_frame)=="transcript_start"] <- "L_transcript_start"


genes_withclosest_GL_peaks_frame<-subset(genes_withclosest_GL_peaks_frame,select=-c(gene_type))

genes_withclosest_GL_peaks_frame<-merge(genes_withclosest_GL_peaks_frame,gene_type,
                                        by.x=c("genes_withclosest_GL_peaks"),
                                        by.y=c("ID"),
                                        all.x=TRUE,
                                        all.y=FALSE)

genes_withclosest_GL_peaks_frame<-genes_withclosest_GL_peaks_frame[!duplicated(genes_withclosest_GL_peaks_frame$genes_withclosest_GL_peaks),]

genes_withclosest_GL_peaks_frame$TSS_peak<-"no"
genes_withclosest_GL_peaks_frame$gene_body_peak<-"no"

for(i in 1:nrow(genes_withclosest_GL_peaks_frame)){
  
  if(as.character(genes_withclosest_GL_peaks_frame$genes_withclosest_GL_peaks[i])%in%all_gene_body){
    
    genes_withclosest_GL_peaks_frame$gene_body_peak[i]<-"yes"
    
  }
  
  if(as.character(genes_withclosest_GL_peaks_frame$genes_withclosest_GL_peaks[i])%in%all_TSS){
    
    genes_withclosest_GL_peaks_frame$TSS_peak[i]<-"yes"
    
  }
  
  
}

setwd(home1)

write.table(genes_withclosest_GL_peaks_frame,paste(home1,"genes_with_only_GL_peaks_TSS_and_genebody.tsv",sep=""),sep="\t",col.names = T,row.names = F,quote=F)

#it works


#########


#LMS_file<-"/home/marcgabriel/Documents/Marina_HOTAIR_lsd1/diff_expression/DESeq2_results/H_P_vs_G_L/sig_diff_upregulated_H_P_vs_G_L.tsv"

#HMS_file<-"/home/marcgabriel/Documents/Marina_HOTAIR_lsd1/diff_expression/DESeq2_results/H_P_vs_G_L/sig_diff_downregulated_H_P_vs_G_L.tsv"

#has been adaptated from the files above : by Marina (filter on foldchange upregulated : > 1.5 ; downregulated < 0.69)
#LMS_HMS_file<-"/media/marcgabriel/eda138bc-95c2-4778-99df-8915815cb86e/Marina_lsd1_test/HMS_LMS.xlsx"

LMS_HMS_file<-paste(current.dir,"HMS_LMS.xlsx",sep="")

LMS_HMS_file<-read.xls(LMS_HMS_file,
                       sheet=1,
                       method="tab",
                       header=T)


LMS_IDs<-unique(as.character(LMS_HMS_file[,2]))[unique(as.character(LMS_HMS_file[,2]))!=""]

HMS_IDs<-unique(as.character(LMS_HMS_file[,1]))[unique(as.character(LMS_HMS_file[,1]))!=""]


#chip_seq_results_G<-read.delim("/home/marcgabriel/Documents/Marina_HOTAIR_lsd1/chip_seq_analysis/G_TSS_with_closest_peak.tsv",header=T)

#chip_seq_results_L<-read.delim("/home/marcgabriel/Documents/Marina_HOTAIR_lsd1/chip_seq_analysis/L_TSS_with_closest_peak.tsv",header=T)

#genes_withclosest_GL_peaks<-unique(as.character(unlist(read.delim("/home/marcgabriel/Documents/Marina_HOTAIR_lsd1/G_L_spe_IDs.txt",header=F))))

#here we use the clean set of gene (those without any other peaks in a area, and with peaks in both replicates for the condition : result3)




result2<-getAllConbinations(list(genes_withclosest_GL_peaks,LMS_IDs,HMS_IDs),list("G&L","LMS_IDs","HMS_IDs"))

write.table(result2,paste(home1,"intersection_genesWithPeaks_G+L","_LMS_HMS.tsv",sep=""),sep="\t",col.names = T,row.names = F,quote=F)


intersected_all<-intersect(intersect(genes_withclosest_GL_peaks,LMS_IDs),HMS_IDs)


#"#c6b291"
all_colors<-c("#A1856F","#376794","#FF008C")


exp_names2<-c(paste("cond1\nG+L","\n(n = ",format(length(genes_withclosest_GL_peaks),big.mark = " "),")",sep=""),
              paste("cond2\nLMS","\n(n = ",format(length(LMS_IDs),big.mark = " "),")",sep=""),
              paste("cond3\nHMS","\n(n = ",format(length(HMS_IDs),big.mark = " "),")",sep=""))




venn.plot <- draw.triple.venn(
  area1 = length(genes_withclosest_GL_peaks),
  area2 =length(LMS_IDs),
  area3 = length(HMS_IDs),
  n12 = length(intersect(genes_withclosest_GL_peaks,LMS_IDs)),
  n13 = length(intersect(genes_withclosest_GL_peaks,HMS_IDs)),
  n23 = length(intersect(LMS_IDs,HMS_IDs)),
  n123 = length(intersected_all),
  category = exp_names2,
  fill = all_colors,
  euler.d = F,
  scaled = F,
  cex = 1,
  cat.cex = 0.8,
  margin = 0.01,
  alpha=c(0.7,0.7,0.7),
  cat.dist=c(0.2,0.2,0.2))

png(filename=paste(home1,"intersection_genesWithPeaks_G+L","_LMS_HMS.png",sep=""),res=1700,width=87,height=87,units ="mm")

#grid.draw(venn.plot,top="Genes with peaks")

grid.arrange(gTree(children=venn.plot),top="Genes with peaks")

dev.off()


cairo_ps(filename=paste(home1,"intersection_genesWithPeaks_G+L","_LMS_HMS.eps",sep=""),width=3.4252,height=3.4252,fallback_resolution = 1700)

grid.arrange(gTree(children=venn.plot),top="Genes with peaks")

dev.off()

plot.new()

venn.plot <- draw.triple.venn(
  area1 = length(genes_withclosest_GL_peaks),
  area2 =length(LMS_IDs),
  area3 = length(HMS_IDs),
  n12 = length(intersect(genes_withclosest_GL_peaks,LMS_IDs)),
  n13 = length(intersect(genes_withclosest_GL_peaks,HMS_IDs)),
  n23 = length(intersect(LMS_IDs,HMS_IDs)),
  n123 = length(intersected_all),
  category = exp_names2,
  fill = all_colors,
  euler.d = F,
  scaled = F,
  cex = 1,
  label.col = rep("transparent",7),
  cat.cex = 0.8,
  margin = 0.01,
  alpha=c(0.7,0.7,0.7),
  cat.dist=c(0.2,0.2,0.2))

png(filename=paste(home1,"intersection_genesWithPeaks_G+L","_LMS_HMS_NoLabels.png",sep=""),res=1700,width=87,height=87,units ="mm")


grid.arrange(gTree(children=venn.plot),top="Genes with peaks")

dev.off()


cairo_ps(filename=paste(home1,"intersection_genesWithPeaks_G+L","_LMS_HMS_NoLabels.eps",sep=""),width=3.4252,height=3.4252,fallback_resolution = 1700)


grid.arrange(gTree(children=venn.plot),top="Genes with peaks")

dev.off()




#HP
########################

genes_withclosest_HP_peaks<-unlist(strsplit(result3$cond3_cond4_IDs," "))

result2<-getAllConbinations(list(genes_withclosest_HP_peaks,LMS_IDs,HMS_IDs),list("H&P","LMS_IDs","HMS_IDs"))

write.table(result2,paste(home1,"intersection_genesWithPeaks_H+P","_LMS_HMS.tsv",sep=""),sep="\t",col.names = T,row.names = F,quote=F)

intersected_all<-intersect(intersect(genes_withclosest_HP_peaks,LMS_IDs),HMS_IDs)

#HMS : #FF008C ; LMS : #376794
#all_colors<-c("#415491","#FF3300","#00FFFF")
all_colors<-c("#415491","#376794","#FF008C")


exp_names2<-c(paste("cond1\nH+P","\n(n = ",format(length(genes_withclosest_HP_peaks),big.mark = " "),")",sep=""),
              paste("cond2\nLMS","\n(n = ",format(length(LMS_IDs),big.mark = " "),")",sep=""),
              paste("cond3\nHMS","\n(n = ",format(length(HMS_IDs),big.mark = " "),")",sep=""))




venn.plot <- draw.triple.venn(
  area1 = length(genes_withclosest_HP_peaks),
  area2 =length(LMS_IDs),
  area3 = length(HMS_IDs),
  n12 = length(intersect(genes_withclosest_HP_peaks,LMS_IDs)),
  n13 = length(intersect(genes_withclosest_HP_peaks,HMS_IDs)),
  n23 = length(intersect(LMS_IDs,HMS_IDs)),
  n123 = length(intersected_all),
  category = exp_names2,
  fill = all_colors,
  euler.d = F,
  scaled = F,
  cex = 1,
  cat.cex = 0.8,
  margin = 0.01,
  alpha=c(0.7,0.7,0.7),
  cat.dist=c(0.2,0.2,0.2))

png(filename=paste(home1,"intersection_genesWithPeaks_H+P","_LMS_HMS.png",sep=""),res=1700,width=87,height=87,units ="mm")

#grid.draw(venn.plot,top="Genes with peaks")

grid.arrange(gTree(children=venn.plot),top="Genes with peaks")

dev.off()


cairo_ps(filename=paste(home1,"intersection_genesWithPeaks_H+P","_LMS_HMS.eps",sep=""),width=3.4252,height=3.4252,fallback_resolution = 1700)

grid.arrange(gTree(children=venn.plot),top="Genes with peaks")

dev.off()

plot.new()

venn.plot <- draw.triple.venn(
  area1 = length(genes_withclosest_HP_peaks),
  area2 =length(LMS_IDs),
  area3 = length(HMS_IDs),
  n12 = length(intersect(genes_withclosest_HP_peaks,LMS_IDs)),
  n13 = length(intersect(genes_withclosest_HP_peaks,HMS_IDs)),
  n23 = length(intersect(LMS_IDs,HMS_IDs)),
  n123 = length(intersected_all),
  category = exp_names2,
  fill = all_colors,
  euler.d = F,
  scaled = F,
  cex = 1,
  cat.cex = 0.8,
  label.col = rep("transparent",7),
  margin = 0.01,
  alpha=c(0.7,0.7,0.7),
  cat.dist=c(0.2,0.2,0.2))

png(filename=paste(home1,"intersection_genesWithPeaks_H+P","_LMS_HMS_NoLabels.png",sep=""),res=1700,width=87,height=87,units ="mm")

#grid.draw(venn.plot,top="Genes with peaks")

grid.arrange(gTree(children=venn.plot),top="Genes with peaks")

dev.off()

cairo_ps(filename=paste(home1,"intersection_genesWithPeaks_H+P","_LMS_HMS_NoLabels.eps",sep=""),width=3.4252,height=3.4252,fallback_resolution = 1700)

grid.arrange(gTree(children=venn.plot),top="Genes with peaks")

dev.off()

plot.new()




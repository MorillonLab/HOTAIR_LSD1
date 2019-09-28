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
  stop("missing arguments !")
}


#home<-"/media/marcgabriel/Transcend/LSD1_metagenes/Marina_peaks_analysis_rep_not_merged_genebody/"
home<-args[1]

official_annotation<-args[2]

LMS_HMS_file<-args[3]

gene_type<-args[4]



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
source(paste(current.dir,"getAllCombinations.R"),sep="")



official_annotation<-read.delim(file=pipe(paste("grep -v ^# ",official_annotation," |grep -P \"\\tgene\\t\"",sep="")),
                                sep="\t",check.names=F,header=F)

gene_type<-read.delim(gene_type,sep="\t",check.names=F,header=T)

names(official_annotation) <-c("chromosome","source","feature","start","end","score","strand","frame","description")
official_annotation$gene_id<-gsub("gene_id=","",grep("gene_id=",unlist(strsplit(as.character(official_annotation$description),";")),value=T))

#extract gene name (HUGO ID)
official_annotation$gene_name<-gsub("gene_name=","",grep("gene_name=",unlist(strsplit(as.character(official_annotation$description),";")),value=T))

official_annotation$gene_type<-gsub("gene_type=","",grep("gene_type",unlist(strsplit(as.character(official_annotation$description),";")),value=T))

official_annotation<-official_annotation[,c("gene_name","gene_type","chromosome","start","end","strand")]




setwd(home)

G_genes_with_rep_peaks<-read.delim(paste(home,"G_genes_with_peaks_in_both_rep.tsv",sep=""),header=T)
G_genes_with_rep_peaks_back<-G_genes_with_rep_peaks
G_genes_with_rep_peaks<-unique(as.character(G_genes_with_rep_peaks$gene_name))


H_genes_with_rep_peaks<-read.delim(paste(home,"H_genes_with_peaks_in_both_rep.tsv",sep=""),header=T)
H_genes_with_rep_peaks<-unique(as.character(H_genes_with_rep_peaks$gene_name))

P_genes_with_rep_peaks<-read.delim(paste(home,"P_genes_with_peaks_in_both_rep.tsv",sep=""),header=T)
P_genes_with_rep_peaks<-unique(as.character(P_genes_with_rep_peaks$gene_name))

L_genes_with_rep_peaks<-read.delim(paste(home,"L_genes_with_peaks_in_both_rep.tsv",sep=""),header=T)
L_genes_with_rep_peaks_back<-L_genes_with_rep_peaks
L_genes_with_rep_peaks<-unique(as.character(L_genes_with_rep_peaks$gene_name))

G_genes_with_just_one_rep_peaks<-read.delim(paste(home,"G_genes_with_peaks_in_just_one_rep.tsv",sep=""),header=T)
G_genes_with_just_one_rep_peaks_back<-G_genes_with_just_one_rep_peaks
G_genes_with_just_one_rep_peaks<-unique(as.character(G_genes_with_just_one_rep_peaks$gene_name))


H_genes_with_just_one_rep_peaks<-read.delim(paste(home,"H_genes_with_peaks_in_just_one_rep.tsv",sep=""),header=T)
H_genes_with_just_one_rep_peaks<-unique(as.character(H_genes_with_just_one_rep_peaks$gene_name))

P_genes_with_just_one_rep_peaks<-read.delim(paste(home,"P_genes_with_peaks_in_just_one_rep.tsv",sep=""),header=T)
P_genes_with_just_one_rep_peaks<-unique(as.character(P_genes_with_just_one_rep_peaks$gene_name))

L_genes_with_just_one_rep_peaks<-read.delim(paste(home,"L_genes_with_peaks_in_just_one_rep.tsv",sep=""),header=T)
L_genes_with_just_one_rep_peaks_back<-L_genes_with_just_one_rep_peaks
L_genes_with_just_one_rep_peaks<-unique(as.character(L_genes_with_just_one_rep_peaks$gene_name))


## old code ##
##############
dist<-0

#/home/marcgabriel/Documents/Marina_HOTAIR_lsd1/chip_seq_analysis/

#"/media/marcgabriel/eda138bc-95c2-4778-99df-8915815cb86e/Marina_lsd1_test/"

chip_seq_results_G<-read.delim(paste(home,"G_TSS_with_closest_peak.tsv",sep=""),header=T)

chip_seq_results_L<-read.delim(paste(home,"L_TSS_with_closest_peak.tsv",sep=""),header=T)

chip_seq_IDs_G<-unique(as.character(chip_seq_results_G[which(abs(chip_seq_results_G$dist_peak_to_transcript_loc)>=dist),]$gene_name))

chip_seq_IDs_G<-unique(chip_seq_IDs_G)[unique(chip_seq_IDs_G)!=""]

chip_seq_IDs_L<-unique(as.character(chip_seq_results_L[which(abs(chip_seq_results_L$dist_peak_to_transcript_loc)>=dist),]$gene_name))

chip_seq_IDs_L<-unique(chip_seq_IDs_L)[unique(chip_seq_IDs_L)!=""]


chip_seq_results_H<-read.delim(paste(home,"H_TSS_with_closest_peak.tsv",sep=""),header=T)

chip_seq_results_P<-read.delim(paste(home,"P_TSS_with_closest_peak.tsv",sep=""),header=T)

chip_seq_IDs_H<-unique(as.character(chip_seq_results_H[which(abs(chip_seq_results_H$dist_peak_to_transcript_loc)>=dist),]$gene_name))

chip_seq_IDs_H<-unique(chip_seq_IDs_H)[unique(chip_seq_IDs_H)!=""]

chip_seq_IDs_P<-unique(as.character(chip_seq_results_P[which(abs(chip_seq_results_P$dist_peak_to_transcript_loc)>=dist),]$gene_name))

chip_seq_IDs_P<-unique(chip_seq_IDs_P)[unique(chip_seq_IDs_P)!=""]


#################

result2<-getAllConbinations(list(chip_seq_IDs_G,chip_seq_IDs_L,chip_seq_IDs_H,chip_seq_IDs_P),list("G","L","H","P"))

intersected_all<-intersect(intersect(intersect(chip_seq_IDs_G,chip_seq_IDs_H),chip_seq_IDs_P),chip_seq_IDs_L)

exp_names2<-c(paste("cond1\nG","\n(n = ",format(length(chip_seq_IDs_G),big.mark = " "),")",sep=""),
              paste("cond2\nL","\n(n = ",format(length(chip_seq_IDs_L),big.mark = " "),")",sep=""),
              paste("cond3\nH","\n(n = ",format(length(chip_seq_IDs_H),big.mark = " "),")",sep=""),
              paste("cond4\nP","\n(n = ",format(length(chip_seq_IDs_P),big.mark = " "),")",sep=""))

all_colors<-c("#7f7f7f","#e46c0a","#604a7b","#31859c")

#G,H,P,L
#area1,area3,area4,area2
venn.plot <- draw.quad.venn(
  area1 = length(chip_seq_IDs_G),
  area2 = length(chip_seq_IDs_L),
  area3 =length(chip_seq_IDs_H),
  area4 = length(chip_seq_IDs_P),
  n12 = length(intersect(chip_seq_IDs_G,chip_seq_IDs_L)),
  n13 = length(intersect(chip_seq_IDs_G,chip_seq_IDs_H)),
  n14 = length(intersect(chip_seq_IDs_G,chip_seq_IDs_P)),
  n23 = length(intersect(chip_seq_IDs_L,chip_seq_IDs_H)),
  n24 = length(intersect(chip_seq_IDs_L,chip_seq_IDs_P)),
  n34 = length(intersect(chip_seq_IDs_H,chip_seq_IDs_P)),
  n123 = length(intersect(intersect(chip_seq_IDs_G,chip_seq_IDs_L),chip_seq_IDs_H)),
  n124 = length(intersect(intersect(chip_seq_IDs_G,chip_seq_IDs_L),chip_seq_IDs_P)),
  n134 = length(intersect(intersect(chip_seq_IDs_G,chip_seq_IDs_H),chip_seq_IDs_P)),
  n234 = length(intersect(intersect(chip_seq_IDs_L,chip_seq_IDs_H),chip_seq_IDs_P)),
  n1234 = length(intersected_all),
  category = exp_names2,
  fill = all_colors,
  cex = 1,
  cat.cex = 0.8,
  margin = 0.04,
  alpha=c(0.7,0.7,0.7,0.7),
  cat.dist=c(0.3,0.3,0.2,0.2))


png(filename=paste(home,"intersection_genesWithPeaks_Minus0kb_Plus0kb_G_H_P_L.png",sep=""),res=1700,width=87,height=87,units ="mm")

#grid.draw(venn.plot)
grid.arrange(gTree(children=venn.plot),top="Genes with peaks")

dev.off()

cairo_ps(filename=paste(home,"intersection_genesWithPeaks_Minus0kb_Plus0kb_G_H_P_L.eps",sep=""),width=3.4252,height=3.4252,fallback_resolution = 1700)

grid.arrange(gTree(children=venn.plot),top="Genes with peaks")

dev.off()


venn.plot <- draw.quad.venn(
  area1 = length(chip_seq_IDs_G),
  area2 = length(chip_seq_IDs_L),
  area3 =length(chip_seq_IDs_H),
  area4 = length(chip_seq_IDs_P),
  n12 = length(intersect(chip_seq_IDs_G,chip_seq_IDs_L)),
  n13 = length(intersect(chip_seq_IDs_G,chip_seq_IDs_H)),
  n14 = length(intersect(chip_seq_IDs_G,chip_seq_IDs_P)),
  n23 = length(intersect(chip_seq_IDs_L,chip_seq_IDs_H)),
  n24 = length(intersect(chip_seq_IDs_L,chip_seq_IDs_P)),
  n34 = length(intersect(chip_seq_IDs_H,chip_seq_IDs_P)),
  n123 = length(intersect(intersect(chip_seq_IDs_G,chip_seq_IDs_L),chip_seq_IDs_H)),
  n124 = length(intersect(intersect(chip_seq_IDs_G,chip_seq_IDs_L),chip_seq_IDs_P)),
  n134 = length(intersect(intersect(chip_seq_IDs_G,chip_seq_IDs_H),chip_seq_IDs_P)),
  n234 = length(intersect(intersect(chip_seq_IDs_L,chip_seq_IDs_H),chip_seq_IDs_P)),
  n1234 = length(intersected_all),
  category = exp_names2,
  fill = all_colors,
  cex = 1,
  cat.cex =0.8,
  label.col = rep("transparent", 15),
  margin = 0.04,
  alpha=c(0.7,0.7,0.7,0.7),
  cat.dist=c(0.3,0.3,0.2,0.2))

png(filename=paste(home,"intersection_genesWithPeaks_Minus0kb_Plus0kb_G_H_P_L_NoLabels.png",sep=""),res=1700,width=87,height=87,units ="mm")

#grid.draw(venn.plot)

grid.arrange(gTree(children=venn.plot),top="Genes with peaks")

dev.off()

cairo_ps(filename=paste(home,"intersection_genesWithPeaks_Minus0kb_Plus0kb_G_H_P_L_NoLabels.eps",sep=""),width=3.4252,height=3.4252,fallback_resolution = 1700)

grid.arrange(gTree(children=venn.plot),top="Genes with peaks")

dev.off()

#dev
#restricted overlap

removal_list<-c()
G_clean_all<-result2[,grep("IDs",grep("cond1",names(result2),value=T),value=T)]

tmp<-c()
for(i in 1:ncol(G_clean_all)){
  
  toto<-unlist(strsplit(as.character(unlist(G_clean_all[i]))," "))
  tmp<-c(tmp,toto)
}
G_clean_all<-unique(tmp)
removal_list<-removal_list<-c(removal_list,G_clean_all[!G_clean_all%in%G_genes_with_rep_peaks])
#G_clean_all<-G_clean_all[G_clean_all%in%G_genes_with_rep_peaks]

L_clean_all<-unique(result2[,grep("IDs",grep("cond2",names(result2),value=T),value=T)])

L_clean_all<-result2[,grep("IDs",grep("cond2",names(result2),value=T),value=T)]

tmp<-c()
for(i in 1:ncol(L_clean_all)){
  
  toto<-unlist(strsplit(as.character(unlist(L_clean_all[i]))," "))
  tmp<-c(tmp,toto)
}
L_clean_all<-unique(tmp)
removal_list<-removal_list<-c(removal_list,L_clean_all[!L_clean_all%in%L_genes_with_rep_peaks])



H_clean_all<-unique(result2[,grep("IDs",grep("cond3",names(result2),value=T),value=T)])

H_clean_all<-result2[,grep("IDs",grep("cond3",names(result2),value=T),value=T)]

tmp<-c()
for(i in 1:ncol(H_clean_all)){
  
  toto<-unlist(strsplit(as.character(unlist(H_clean_all[i]))," "))
  tmp<-c(tmp,toto)
}
H_clean_all<-unique(tmp)
removal_list<-removal_list<-c(removal_list,H_clean_all[!H_clean_all%in%H_genes_with_rep_peaks])



P_clean_all<-unique(result2[,grep("IDs",grep("cond4",names(result2),value=T),value=T)])

P_clean_all<-result2[,grep("IDs",grep("cond4",names(result2),value=T),value=T)]

tmp<-c()
for(i in 1:ncol(P_clean_all)){
  
  toto<-unlist(strsplit(as.character(unlist(P_clean_all[i]))," "))
  tmp<-c(tmp,toto)
}
P_clean_all<-unique(tmp)
removal_list<-removal_list<-c(removal_list,P_clean_all[!P_clean_all%in%P_genes_with_rep_peaks])



G_clean_all<-G_clean_all[!G_clean_all%in%removal_list]
L_clean_all<-L_clean_all[!L_clean_all%in%removal_list]
H_clean_all<-H_clean_all[!H_clean_all%in%removal_list]
P_clean_all<-P_clean_all[!P_clean_all%in%removal_list]

G_clean_all<-G_clean_all[!G_clean_all%in%c(H_genes_with_just_one_rep_peaks,L_genes_with_just_one_rep_peaks,P_genes_with_just_one_rep_peaks)]
L_clean_all<-L_clean_all[!L_clean_all%in%c(G_genes_with_just_one_rep_peaks,H_genes_with_just_one_rep_peaks,P_genes_with_just_one_rep_peaks)]
H_clean_all<-H_clean_all[!H_clean_all%in%c(H_genes_with_just_one_rep_peaks,L_genes_with_just_one_rep_peaks,P_genes_with_just_one_rep_peaks)]
P_clean_all<-P_clean_all[!P_clean_all%in%c(G_genes_with_just_one_rep_peaks,H_genes_with_just_one_rep_peaks,L_genes_with_just_one_rep_peaks)]



result3<-getAllConbinations(list(G_clean_all,L_clean_all,H_clean_all,P_clean_all),list("G","L","H","P"))

write.table(result3,file=paste(home,"intersection_genesWithPeaks_Minus0kb_Plus0kb_G_H_P_L",".tsv",sep=""), sep='\t',row.names=F, col.names=T, quote=F)

intersected_all<-intersect(intersect(intersect(G_clean_all,H_clean_all),P_clean_all),L_clean_all)

exp_names2<-c(paste("cond1\nG","\n(n = ",format(length(G_clean_all),big.mark = " "),")",sep=""),
              paste("cond2\nL","\n(n = ",format(length(L_clean_all),big.mark = " "),")",sep=""),
              paste("cond3\nH","\n(n = ",format(length(H_clean_all),big.mark = " "),")",sep=""),
              paste("cond4\nP","\n(n = ",format(length(P_clean_all),big.mark = " "),")",sep=""))

all_colors<-c("#7f7f7f","#e46c0a","#604a7b","#31859c")

#G,H,P,L
#area1,area3,area4,area2
venn.plot <- draw.quad.venn(
  area1 = length(G_clean_all),
  area2 = length(L_clean_all),
  area3 =length(H_clean_all),
  area4 = length(P_clean_all),
  n12 = length(intersect(G_clean_all,L_clean_all)),
  n13 = length(intersect(G_clean_all,H_clean_all)),
  n14 = length(intersect(G_clean_all,P_clean_all)),
  n23 = length(intersect(L_clean_all,H_clean_all)),
  n24 = length(intersect(L_clean_all,P_clean_all)),
  n34 = length(intersect(H_clean_all,P_clean_all)),
  n123 = length(intersect(intersect(G_clean_all,L_clean_all),H_clean_all)),
  n124 = length(intersect(intersect(G_clean_all,L_clean_all),P_clean_all)),
  n134 = length(intersect(intersect(G_clean_all,H_clean_all),P_clean_all)),
  n234 = length(intersect(intersect(L_clean_all,H_clean_all),P_clean_all)),
  n1234 = length(intersected_all),
  category = exp_names2,
  fill = all_colors,
  cex = 1,
  cat.cex = 0.8,
  margin = 0.04,
  alpha=c(0.7,0.7,0.7,0.7),
  cat.dist=c(0.3,0.3,0.2,0.2))

png(filename=paste(home,"intersection_genesWithPeaks_Minus0kb_Plus0kb_G_H_P_L.png",sep=""),res=1700,width=87,height=87,units ="mm")

#grid.draw(venn.plot)
grid.arrange(gTree(children=venn.plot),top="Genes with peaks")

dev.off()

cairo_ps(filename=paste(home,"intersection_genesWithPeaks_Minus0kb_Plus0kb_G_H_P_L.eps",sep=""),width=3.4252,height=3.4252,fallback_resolution = 1700)

grid.arrange(gTree(children=venn.plot),top="Genes with peaks")

dev.off()

venn.plot <- draw.quad.venn(
  area1 = length(G_clean_all),
  area2 = length(L_clean_all),
  area3 =length(H_clean_all),
  area4 = length(P_clean_all),
  n12 = length(intersect(G_clean_all,L_clean_all)),
  n13 = length(intersect(G_clean_all,H_clean_all)),
  n14 = length(intersect(G_clean_all,P_clean_all)),
  n23 = length(intersect(L_clean_all,H_clean_all)),
  n24 = length(intersect(L_clean_all,P_clean_all)),
  n34 = length(intersect(H_clean_all,P_clean_all)),
  n123 = length(intersect(intersect(G_clean_all,L_clean_all),H_clean_all)),
  n124 = length(intersect(intersect(G_clean_all,L_clean_all),P_clean_all)),
  n134 = length(intersect(intersect(G_clean_all,H_clean_all),P_clean_all)),
  n234 = length(intersect(intersect(L_clean_all,H_clean_all),P_clean_all)),
  n1234 = length(intersected_all),
  category = exp_names2,
  fill = all_colors,
  cex = 1,
  cat.cex = 0.8,
  label.col = rep("transparent", 15),
  
  margin = 0.04,
  alpha=c(0.7,0.7,0.7,0.7),
  cat.dist=c(0.3,0.3,0.2,0.2))

png(filename=paste(home,"intersection_genesWithPeaks_Minus0kb_Plus0kb_G_H_P_L_NoLabels.png",sep=""),res=1700,width=87,height=87,units ="mm")

#grid.draw(venn.plot)
grid.arrange(gTree(children=venn.plot),top="Genes with peaks")

dev.off()

cairo_ps(filename=paste(home,"intersection_genesWithPeaks_Minus0kb_Plus0kb_G_H_P_L_NoLabels.eps",sep=""),width=3.4252,height=3.4252,fallback_resolution = 1700)

grid.arrange(gTree(children=venn.plot),top="Genes with peaks")

dev.off()


#end of dev




#########

#LMS_file<-"/home/marcgabriel/Documents/Marina_HOTAIR_lsd1/diff_expression/DESeq2_results/H_P_vs_G_L/sig_diff_upregulated_H_P_vs_G_L.tsv"

#HMS_file<-"/home/marcgabriel/Documents/Marina_HOTAIR_lsd1/diff_expression/DESeq2_results/H_P_vs_G_L/sig_diff_downregulated_H_P_vs_G_L.tsv"

#has been adaptated from the files above : by Marina (filter on foldchange upregulated : > 1.5 ; downregulated < 0.69)
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
genes_withclosest_GL_peaks<-unlist(strsplit(result3$cond1_cond2_IDs," "))

genes_withclosest_GL_peaks_frame<-data.frame(genes_withclosest_GL_peaks)

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

genes_withclosest_GL_peaks_frame<-genes_withclosest_GL_peaks_frame[which(!is.na(genes_withclosest_GL_peaks_frame$dist_peak_to_transcript_loc)),]





colnames(genes_withclosest_GL_peaks_frame)[colnames(genes_withclosest_GL_peaks_frame)=="dist_peak_to_transcript_loc"] <- "G_dist_peak_to_transcript_loc"

colnames(genes_withclosest_GL_peaks_frame)[colnames(genes_withclosest_GL_peaks_frame)=="transcript_start"] <- "G_transcript_start"

genes_withclosest_GL_peaks_frame<-merge(genes_withclosest_GL_peaks_frame
                                        ,L_genes_with_rep_peaks_back[,c("gene_name","transcript_start","dist_peak_to_transcript_loc","chromosome")],
                                        by.x=c("genes_withclosest_GL_peaks","chromosome"),
                                        by.y=c("gene_name","chromosome"),
                                        all.x=TRUE,
                                        all.y=FALSE)

genes_withclosest_GL_peaks_frame<-genes_withclosest_GL_peaks_frame[which(!is.na(genes_withclosest_GL_peaks_frame$dist_peak_to_transcript_loc)),]


colnames(genes_withclosest_GL_peaks_frame)[colnames(genes_withclosest_GL_peaks_frame)=="dist_peak_to_transcript_loc"] <- "L_dist_peak_to_transcript_loc"
colnames(genes_withclosest_GL_peaks_frame)[colnames(genes_withclosest_GL_peaks_frame)=="transcript_start"] <- "L_transcript_start"


genes_withclosest_GL_peaks_frame<-subset(genes_withclosest_GL_peaks_frame,select=-c(gene_type))

genes_withclosest_GL_peaks_frame<-merge(genes_withclosest_GL_peaks_frame,gene_type,
                                        by.x=c("genes_withclosest_GL_peaks"),
                                        by.y=c("ID"),
                                        all.x=TRUE,
                                        all.y=FALSE)

genes_withclosest_GL_peaks_frame<-genes_withclosest_GL_peaks_frame[!duplicated(genes_withclosest_GL_peaks_frame$genes_withclosest_GL_peaks),]

write.table(genes_withclosest_GL_peaks_frame,paste(home,"genes_with_only_GL_peaks.tsv",sep=""),sep="\t",col.names = T,row.names = F,quote=F)







result2<-getAllConbinations(list(genes_withclosest_GL_peaks,LMS_IDs,HMS_IDs),list("G&L","LMS_IDs","HMS_IDs"))

write.table(result2,paste(home,"intersection_genesWithPeaks_G+L","_LMS_HMS.tsv",sep=""),sep="\t",col.names = T,row.names = F,quote=F)


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

png(filename=paste(home,"intersection_genesWithPeaks_G+L","_LMS_HMS.png",sep=""),res=1700,width=87,height=87,units ="mm")

#grid.draw(venn.plot,top="Genes with peaks")

grid.arrange(gTree(children=venn.plot),top="Genes with peaks")

dev.off()


cairo_ps(filename=paste(home,"intersection_genesWithPeaks_G+L","_LMS_HMS.eps",sep=""),width=3.4252,height=3.4252,fallback_resolution = 1700)

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

png(filename=paste(home,"intersection_genesWithPeaks_G+L","_LMS_HMS_NoLabels.png",sep=""),res=1700,width=87,height=87,units ="mm")


grid.arrange(gTree(children=venn.plot),top="Genes with peaks")

dev.off()


cairo_ps(filename=paste(home,"intersection_genesWithPeaks_G+L","_LMS_HMS_NoLabels.eps",sep=""),width=3.4252,height=3.4252,fallback_resolution = 1700)


grid.arrange(gTree(children=venn.plot),top="Genes with peaks")

dev.off()




#HP
########################

genes_withclosest_HP_peaks<-unlist(strsplit(result3$cond3_cond4_IDs," "))

result2<-getAllConbinations(list(genes_withclosest_HP_peaks,LMS_IDs,HMS_IDs),list("H&P","LMS_IDs","HMS_IDs"))

write.table(result2,paste(home,"intersection_genesWithPeaks_H+P","_LMS_HMS.tsv",sep=""),sep="\t",col.names = T,row.names = F,quote=F)

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

png(filename=paste(home,"intersection_genesWithPeaks_H+P","_LMS_HMS.png",sep=""),res=1700,width=87,height=87,units ="mm")

#grid.draw(venn.plot,top="Genes with peaks")

grid.arrange(gTree(children=venn.plot),top="Genes with peaks")

dev.off()


cairo_ps(filename=paste(home,"intersection_genesWithPeaks_H+P","_LMS_HMS.eps",sep=""),width=3.4252,height=3.4252,fallback_resolution = 1700)

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

png(filename=paste(home,"intersection_genesWithPeaks_H+P","_LMS_HMS_NoLabels.png",sep=""),res=1700,width=87,height=87,units ="mm")

#grid.draw(venn.plot,top="Genes with peaks")

grid.arrange(gTree(children=venn.plot),top="Genes with peaks")

dev.off()

cairo_ps(filename=paste(home,"intersection_genesWithPeaks_H+P","_LMS_HMS_NoLabels.eps",sep=""),width=3.4252,height=3.4252,fallback_resolution = 1700)

grid.arrange(gTree(children=venn.plot),top="Genes with peaks")

dev.off()

plot.new()




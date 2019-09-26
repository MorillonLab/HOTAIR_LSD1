#!/usr/bin/env Rscript

args<-commandArgs(TRUE)

library(gridExtra)

library(VennDiagram)


home<-args[1]

conditions<-args[2]

counts<-args[3]

all_total<-args[4]

all_colors<-args[5]

prefix<-args[6]



cat("home is : ",home,"\ncounts are : ",counts,"\nall total are : ",all_total,"\nconditions are : ",conditions,"\nall colors are ",all_colors)

conditions<-unlist(strsplit(conditions,","))

counts<-as.numeric(unlist(strsplit(counts,",")))

all_total<-as.numeric(unlist(strsplit(all_total,",")))

all_colors<-unlist(strsplit(all_colors,","))


exp_names2<-c(paste(conditions[1]," (n=",all_total[1],")",sep=""),
              paste(conditions[2]," (n=",all_total[2],")",sep=""),
              paste(conditions[3]," (n=",all_total[3],")",sep=""),
              paste(conditions[4]," (n=",all_total[4],")",sep=""))






venn.plot <- draw.quad.venn(
  area1 = all_total[1],
  area2 = all_total[2],
  area3 =all_total[3],
  area4 = all_total[4],
  n12 = counts[5]+counts[11]+counts[12]+counts[15],
  n13 = counts[6]+counts[11]+counts[13]+counts[15],
  n14 = counts[7]+counts[12]+counts[13]+counts[15],
  n23 = counts[8]+counts[11]+counts[14]+counts[15],
  n24 = counts[9]+counts[12]+counts[14]+counts[15],
  n34 = counts[10]+counts[13]+counts[14]+counts[15],
  n123 = counts[11]+counts[15],
  n124 = counts[12]+counts[15],
  n134 = counts[13]+counts[15],
  n234 = counts[14]+counts[15],
  n1234 = counts[15],
  category = exp_names2,
  fill = all_colors,
  cex = 1,
  cat.cex = 0.8,
  margin = 0.05,
  alpha=c(0.7,0.7,0.7,0.7),
  cat.dist=c(0.3,0.3,0.2,0.2))

png(filename=paste(home,"venn","_",paste(conditions[1],conditions[2],conditions[3],conditions[4],sep="_"),"_",prefix,".png",sep=""),res=1700,width=87,height=87,units ="mm")

grid.arrange(gTree(children=venn.plot),top="peaks grouped per regions")


dev.off()





venn.plot<-draw.quad.venn(
  area1 = all_total[1],
  area2 = all_total[2],
  area3 =all_total[3],
  area4 = all_total[4],
  n12 = counts[5]+counts[11]+counts[12]+counts[15],
  n13 = counts[6]+counts[11]+counts[13]+counts[15],
  n14 = counts[7]+counts[12]+counts[13]+counts[15],
  n23 = counts[8]+counts[11]+counts[14]+counts[15],
  n24 = counts[9]+counts[12]+counts[14]+counts[15],
  n34 = counts[10]+counts[13]+counts[14]+counts[15],
  n123 = counts[11]+counts[15],
  n124 = counts[12]+counts[15],
  n134 = counts[13]+counts[15],
  n234 = counts[14]+counts[15],
  n1234 = counts[15],
  category = exp_names2,
  fill = all_colors,
  cex = 1,
  cat.cex = 0.8,
  margin = 0.05,
  alpha=c(0.7,0.7,0.7,0.7),
  cat.dist=c(0.3,0.3,0.2,0.2))

cairo_ps(filename=paste(home,"venn","_",paste(conditions[1],conditions[2],conditions[3],conditions[4],sep="_"),"_",prefix,".eps",sep=""),width=3.4252,height=3.4252,fallback_resolution = 1700)

grid.arrange(gTree(children=venn.plot),top="peaks grouped per regions")

dev.off()



venn.plot <- draw.quad.venn(
  area1 = all_total[1],
  area2 = all_total[2],
  area3 =all_total[3],
  area4 = all_total[4],
  n12 = counts[5]+counts[11]+counts[12]+counts[15],
  n13 = counts[6]+counts[11]+counts[13]+counts[15],
  n14 = counts[7]+counts[12]+counts[13]+counts[15],
  n23 = counts[8]+counts[11]+counts[14]+counts[15],
  n24 = counts[9]+counts[12]+counts[14]+counts[15],
  n34 = counts[10]+counts[13]+counts[14]+counts[15],
  n123 = counts[11]+counts[15],
  n124 = counts[12]+counts[15],
  n134 = counts[13]+counts[15],
  n234 = counts[14]+counts[15],
  n1234 = counts[15],
  category = exp_names2,
  fill = all_colors,
  cex = 1,
  label.col = rep("transparent", 15),
  cat.cex = 0.8,
  margin = 0.05,
  alpha=c(0.7,0.7,0.7,0.7),
  cat.dist=c(0.3,0.3,0.2,0.2))

png(filename=paste(home,"venn","_",paste(conditions[1],conditions[2],conditions[3],conditions[4],sep="_"),"_",prefix,"_NoLabels.png",sep=""),res=1700,width=87,height=87,units ="mm")

grid.arrange(gTree(children=venn.plot),top="peaks grouped per regions")

dev.off()


venn.plot<-draw.quad.venn(
  area1 = all_total[1],
  area2 = all_total[2],
  area3 =all_total[3],
  area4 = all_total[4],
  n12 = counts[5]+counts[11]+counts[12]+counts[15],
  n13 = counts[6]+counts[11]+counts[13]+counts[15],
  n14 = counts[7]+counts[12]+counts[13]+counts[15],
  n23 = counts[8]+counts[11]+counts[14]+counts[15],
  n24 = counts[9]+counts[12]+counts[14]+counts[15],
  n34 = counts[10]+counts[13]+counts[14]+counts[15],
  n123 = counts[11]+counts[15],
  n124 = counts[12]+counts[15],
  n134 = counts[13]+counts[15],
  n234 = counts[14]+counts[15],
  n1234 = counts[15],
  category = exp_names2,
  fill = all_colors,
  cex = 1,
  label.col = rep("transparent", 15),
  cat.cex = 0.8,
  margin = 0.05,
  alpha=c(0.7,0.7,0.7,0.7),
  cat.dist=c(0.3,0.3,0.2,0.2))

cairo_ps(filename=paste(home,"venn","_",paste(conditions[1],conditions[2],conditions[3],conditions[4],sep="_"),"_",prefix,"_NoLabels.eps",sep=""),width=3.4252,height=3.4252,fallback_resolution = 1700)

  grid.arrange(gTree(children=venn.plot),top="peaks grouped per regions")

dev.off()


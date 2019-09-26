#!/usr/bin/env Rscript

args<-commandArgs(TRUE)

library(gridExtra)

library(VennDiagram)

library(eulerr)

home<-args[1]

file<-args[2]

prefix<-args[3]

mytitle<-args[4]

#prefix<-"cutoff_50_percent_overlap"


#file<-"/home/marcgabriel/Documents/Marina_P/data_intersections/all_intersections.tsv"

#home<-"/home/marcgabriel/Documents/Marina_P/data_intersections/"

mytable<-read.delim(file=file,sep = "\t",check.names = F,header=T)

mytable$color1<-""

mytable$color2<-""


#set colors for each condition
all_conditions<-unique(c(as.character(mytable$cond1),as.character(mytable$cond2)))

color_table<-data.frame()

all_colors<-rainbow(length(all_conditions))
#all_colors<-c("red","cornflowerblue","gold","green","purple","cyan")

all_colors[1]<-"#7f7f7f"

all_colors[2]<-"#e46c0a"

for(i in 1:length(all_conditions)){
  
  one_frame<-data.frame(condition=all_conditions[i],color=all_colors[i])
  
  color_table<-rbind(color_table,one_frame)
  
}

for(i in 1:nrow(mytable)){
  
  mytable$color1[i]<-as.character(color_table[which(as.character(color_table$condition)==as.character(mytable$cond1[i])),]$color)
  
  mytable$color2[i]<-as.character(color_table[which(as.character(color_table$condition)==as.character(mytable$cond2[i])),]$color)
  
  
}

# venn_list<-list()
for(i in 1:nrow(mytable)){
  
#   
#   #grid.newpage()
#   one_ven<-draw.pairwise.venn(area1 = mytable$cond1_counts[i],
#                               area2 = mytable$cond2_counts[i],
#                               cross.area =mytable$minimal_overlap[i],
#                               category = c(paste(mytable$cond1[i],"\n(overlap with ",mytable$cond2[i]," : ",round((mytable$minimal_overlap[i]/mytable$cond1_counts[i])*100,digits=2),"%)",sep=""),
#                                            paste(mytable$cond2[i],"\n(overlap with ",mytable$cond1[i]," : ",round((mytable$minimal_overlap[i]/mytable$cond2_counts[i])*100,digits=2),"%)",sep="")),
#                               #fill=c(mytable$color1[i],mytable$color2[i]),
#                               fill=c("red","cornflowerblue"),
#                               ind=F,
#                               cat.pos = c(0,0),
#                               cat.dist = rep(0.05, 2),
#                               cat.cex=rep(1.2,2),
#                               cex=2)
#   
#   #venn_list[[length(venn_list)+1]]<-gTree(children=one_ven)
#   
#   png(filename=paste(home,"venn_diagram_annotations_",mytable$cond1[i],"_vs_",mytable$cond2[i],"_",prefix,".png",sep=""),width=1200,height=1000)
#   
#   grid.draw(one_ven)
#   
#   dev.off()
#   
#   

  
    VennDiag <- euler(c("cond1" =mytable$cond1_counts[i]-mytable$max_overlap_cond1[i], "cond2" =mytable$cond2_counts[i]-mytable$max_overlap_cond2[i], "cond1&cond2" = max(c(mytable$max_overlap_cond2[i],mytable$max_overlap_cond1[i]))))
    
    exp_names2<-c(paste(mytable$cond1[i]," (n=",mytable$cond1_counts[i],")",sep=""),paste(mytable$cond2[i]," (n=",mytable$cond2_counts[i],")",sep=""))
    
    
    png(filename=paste(home,"venn_diagram_",mytable$cond1[i],"_vs_",mytable$cond2[i],"_",prefix,".png",sep=""),width=1200,height=800)
    
    plot.new()
    
    print(plot(VennDiag, counts = TRUE, font=1, cex=2, alpha=0.5,
         fill=c(mytable$color1[i],mytable$color2[i]),
         quantities = list(labels=c(format(mytable$cond1_counts[i]-mytable$max_overlap_cond1[i],big.mark=" "),format(mytable$cond2_counts[i]-mytable$max_overlap_cond2[i],big.mark=" "),paste(format(mytable$max_overlap_cond1[i],big.mark=" ")," ","(",mytable$cond1[i],")"," / ",format(mytable$max_overlap_cond2[i],big.mark=" ")," ","(",mytable$cond2[i],")",sep="")),fontsize = 20),
         labels=NULL,
         legend = list(space = "bottom", columns = 2,labels=exp_names2,cex=2)))
         
         title(mytitle)
   
  dev.off()
  
}

# png(filename=paste(home,"venn_diagram_annotations_",prefix,".png",sep=""),width=1200,height=1000)
# 
# grid.arrange(grobs=venn_list,ncol=1)
# 
# dev.off()




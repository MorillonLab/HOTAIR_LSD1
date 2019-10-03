#!/usr/bin/env Rscript

##### remarks :
###############

#all required packages : DESeq2, ggplot2, RColorBrewer, optparse, gplots, pROC, foreach, doParallel, pheatmap, grid, clusterProfiler, DOSE, BiocParallel, data.table

#no need to go inside this script (unless you want to change something, or if Rscript isn't in /usr/bin/)
#Just call it in the terminal like this : ./DESeq2_script_example.R)

#only 3 arguments are required : -f, -a, and -o

#######



##### try to install & load automatically cran and bioconductor packages
#todo : adapt the download of bioconductor packages to the latest version of R (>= version 3.5)
########################################################################

cran_packages<-c("optparse","ggplot2","gplots","RColorBrewer","pROC","foreach","doParallel","grid","BiocParallel","data.table","pheatmap")

for(one_package in cran_packages){
  
  if(!require(one_package,character.only=T)){
    
    install.packages(one_package,repos='https://cloud.r-project.org/')
    
  }
  
}

bioc_packages<-c("DESeq2","clusterProfiler","DOSE")

source("https://bioconductor.org/biocLite.R")

for(one_package in bioc_packages){
  
  if(!require(one_package,character.only=T)){
    
    BiocInstaller::biocLite(one_package)
    
  }
  
}


suppressPackageStartupMessages(suppressMessages(library(DESeq2)))
suppressPackageStartupMessages(suppressMessages(library(ggplot2)))
suppressPackageStartupMessages(suppressMessages(library(RColorBrewer)))
suppressPackageStartupMessages(suppressMessages(library(optparse)))
suppressPackageStartupMessages(suppressMessages(library(gplots)))
suppressPackageStartupMessages(suppressMessages(library(pROC)))
suppressPackageStartupMessages(suppressMessages(library(foreach)))
suppressPackageStartupMessages(suppressMessages(library(doParallel)))
suppressPackageStartupMessages(suppressMessages(library(pheatmap)))
suppressPackageStartupMessages(suppressMessages(library(grid)))
suppressPackageStartupMessages(suppressMessages(library(clusterProfiler)))
suppressPackageStartupMessages(suppressMessages(library(DOSE)))
suppressPackageStartupMessages(suppressMessages(library(BiocParallel)))
suppressPackageStartupMessages(suppressMessages(library(data.table)))


#####

##### list of arguments/options to supply
##########################################

option_list<-list(
  make_option(c("-f","--files_descriptor"),type="character",help = "tab delimited file with 2 columns, indicating count files and their matching conditions (required) like this : 
              \n\t\tpath/to/file\tcondition_1 (replicate1)
              \n\t\tpath/to/file\tcondition_1 (replicate2)
              \n\t\tpath/to/file\tcondition_2 (replicate1)
              \n\t\tpath/to/file\tcondition_2 (replicate2
              \n\t\tpath/to/file\tcondition_3 (replicate1)
              \n\t\tpath/to/file\tcondition_3 (replicate2)
              \n\t\twith --comparison_type=all, we will have this : 1 vs 2, 1 vs 3, 2 vs 3"),
  make_option(c("-a","--gff"),type="character",help = "annotation in GFF3 format ; 9th column should begin with ID= ; IDs should match with those in the count tables (required)"),
  make_option(c("-o","--output_dir"),type="character",default="",help = "directory to store the results (required)"),
  make_option(c("-t","--comparison_type"),type="character",default="all",help = "type of comparison : all (pairwise comparison), or first_vs_others_pairwise (first condition vs all others, in pairwise), or first_vs_others_combined (first condition vs all others combined), or two_firsts_vs_others_combined ; default=all"),
  make_option(c("-s","--padj_threshold"),type="numeric",default=0.05,help = "padj threshold ; default=0.05"),
  make_option(c("-l","--log2FC_threshold"),type="numeric",default=0,help = "log2FC threshold ; default=0"),
  make_option(c("-v","--read_threshold"),type="numeric",default=0,help = "read threshold ; default=0"),
  make_option(c("-m","--mean_threshold"),type="numeric",default=1,help = "Mean threshold ; default=1"),
  make_option(c("-y","--nb_cores"),type="numeric",default=1,help = "nb cores ; default=1"),
  make_option(c("-n","--additional_annotation"),type="character",default="none",help = "additional annotation ; default=none"),
  make_option(c("-c","--compute_AUC"),type="character",default="no",help = "should we compute the AUC ? ; default=no"),
  make_option(c("-u","--use_given_size_factor"),type="character",default="",help = "use these size factors ; default=\"\""),
  make_option(c("-r","--insert_raw_counts"),type="character",default="no",help = "should we put the raw counts in the results ? ; default=no")
)


opt_parser = OptionParser(option_list=option_list)

opt <- parse_args(OptionParser(option_list=option_list))

#if one required parameter is missing, print a warning and stop the script
if(length(opt$gff)==0 | length(opt$output_dir)==0 | length(opt$files_descriptor)==0){
  
  cat("one required parameter is missing !","\n")
  print_help(opt_parser)
  q("no")
}

#####

##### additionnal functions to use
##################################

##### function to rotate names in pheatmap

#https://slowkow.com/notes/heatmap-tutorial/
draw_colnames_45 <- function (coln, gaps, ...) {
  coord <- pheatmap:::find_coordinates(length(coln), gaps)
  x     <- coord$coord - 0.5 * coord$size
  res   <- grid::textGrob(
    coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"),
    vjust = 0.75, hjust = 1, rot = 45, gp = grid::gpar(...)
  )
  return(res)
}
assignInNamespace(
  x = "draw_colnames",
  value = "draw_colnames_45",
  ns = asNamespace("pheatmap")
)

####### function to scale the x axis of the barplots for the GO terms enrichment

getStepAndMax<-function(my_count){
  if(my_count<=10){
    
    my_step<-2
    my_max<-10
  }else if(my_count<=20){
    
    my_step<-4
    my_max<-20
  }else if(my_count<=50){
    
    my_step<-5
    my_max<-50
  }else if(my_count<=100){
    
    my_step<-10
    my_max<-100
  }else if(my_count<=200){
    
    my_step<-20
    my_max<-200
  }else if(my_count<=400){
    
    my_step<-50
    my_max<-400
  }else if(my_count<=500){
    
    my_step<-50
    my_max<-500
  }else if(my_count<=1000){
    
    my_step<-100
    my_max<-1000
  }else if(my_count<=2000){
    
    my_step<-200
    my_max<-2000
  }else if(my_count<=5000){
    
    my_step<-500
    my_max<-5000
    
  }else if(my_count<=8000){
    
    my_step<-1000
    my_max<-8000
    
  }else if(my_count<=10000){
    
    my_step<-1000
    my_max<-10000
  }
  
  return(list(my_step,my_max))
  
}

##########


############## input data
#########################

home<-opt$output_dir

compute_AUC<-opt$compute_AUC

#additional_info<-"/home/marcgabriel/Documents/Julien_Jarroux/hotair_overexp_featureCounts_output/closest_genes/five_three_and_antisense_reduced.txt"
additional_info<-opt$additional_annotation

nb_cores<-opt$nb_cores

given_size_factor<-opt$use_given_size_factor

insert_raw_counts<-opt$insert_raw_counts

if(given_size_factor!=""){
  
  cat("used size factor :\n",given_size_factor,"\n")
  
  split_given_size_factor<-as.numeric(unlist(strsplit(given_size_factor,",")))
  
  
}

#design
design_file<-opt$files_descriptor
#design_file<-"/home/marcgabriel/Documents/Julien_Jarroux/hotair_overexp_featureCounts_output/counting_exons_holdup_design.txt"

#type of comparison
type<-opt$comparison_type
#type<-"all"

#create home if it doesn't exist
#home<-"/home/marcgabriel/Desktop/deseq_test/"
home<-paste(home,"/",sep="")
home<-gsub("//","/",home)
dir.create(home,showWarnings=F,recursive=T)
setwd(home)

### process gff

#adapted for exons !!
#official_annotation<-read.delim(file=pipe(paste("grep -v ^# ",opt$gff," |grep -P \"\\texon\\t\"",sep="")),
#                                sep="\t",check.names=F,header=F)

#for genes
official_annotation<-read.delim(file=pipe(paste("grep -v ^# ",opt$gff," |grep -P \"\\tgene\\t\"",sep="")),
                                sep="\t",check.names=F,header=F)

#for introns
# official_annotation<-read.delim(file=pipe(paste("grep -v ^# ",opt$gff," |grep -P \"\\tintron\\t\"",sep="")),
#                                 sep="\t",check.names=F,header=F)


names(official_annotation) <-c("chromosome","source","feature","start","end","score","strand","frame","description")

official_annotation$description<-gsub("Gene:","",official_annotation$description)

#extract ID from the description
#official_annotation$gene_id<-lapply(strsplit(as.character(official_annotation$description), "\\;"), "[", 1)
#official_annotation$gene_id<-lapply(strsplit(as.character(official_annotation$gene_id), "\\ID="), "[", 2)

#for genes
official_annotation$gene_id<-gsub("gene_id=","",grep("gene_id=",unlist(strsplit(as.character(official_annotation$description),";")),value=T))

#adapted for exons !!
#official_annotation$gene_id<-gsub("exon_id=","",grep("exon_id=",unlist(strsplit(as.character(official_annotation$description),";")),value=T))
#official_annotation<-official_annotation[!duplicated(official_annotation["gene_id"]),]

#for introns
# official_annotation$gene_id<-gsub("intron_id=","",grep("intron_id=",unlist(strsplit(as.character(official_annotation$description),";")),value=T))
# official_annotation<-official_annotation[!duplicated(official_annotation["gene_id"]),]

#extract gene name (HUGO ID)
official_annotation$gene_name<-gsub("gene_name=","",grep("gene_name=",unlist(strsplit(as.character(official_annotation$description),";")),value=T))

official_annotation<-official_annotation[,c("gene_id","gene_name","chromosome","start","end","strand")]

official_annotation_backup<-official_annotation

design_file<-read.table(file=design_file,sep ="\t",header=F)

backup<-design_file

names(design_file)<-c("files","conditions")

UniqueConditions<-unique(as.character(design_file$conditions))

#re-order the data frame, to have grouped conditions
tmp_design<-data.frame()

for(i in 1:length(UniqueConditions)){
  
  tmp_design<-rbind(tmp_design,design_file[which(as.character(design_file$conditions)==UniqueConditions[i]),])
  
  
}

design_file<-tmp_design


all_reps<-list()
for(i in 1:length(UniqueConditions)){
  
  reps_for_OneCondition<-nrow(design_file[which(as.character(design_file$conditions)==as.character(UniqueConditions[i])),])
  all_reps[[length(all_reps)+1]]<-reps_for_OneCondition
}

rep_list<-list()
for(i in 1:length(all_reps)){
  
  for(j in seq(all_reps[[i]])){
    
    rep_list[[length(rep_list)+1]]<-j
    
  }
  
}


design_file<-cbind(design_file,data.frame(rep_number=as.character(rep_list)))

cat("\nthis is the design that has been set :\n")
write.table(design_file,stdout(),sep='\t',row.names=F, col.names=T, quote=F)
cat("\n-------------------\n")

all_frames<-list()
for(i in 1:nrow(design_file)){
  one_frame<-fread(as.character(design_file[i,1]), sep = "\t",check.names = FALSE,header=F,data.table=F)
  names(one_frame)<-c("gene_id",paste(as.character(design_file[i,2]),as.character(design_file[i,3]),sep="_"))
  
  #one_frame<-one_frame[which(one_frame$gene_id!="MAT_1_long" & one_frame$gene_id!="MAT_1_short"),]
 
  all_frames[[length(all_frames)+1]]<-one_frame
}


#one dataframe with conditions side by side
MyCounTable<-Reduce(function(x, y) merge(x, y, all=TRUE), all_frames)

rm(all_frames);gc()



#it will be used for DESeq, all replicates of a condition should have the same name, and order should be kept like in the dataframe above
condition=factor(as.character(design_file$conditions))

##################################################



#the first column is the row name, and removing of this first column
rownames(MyCounTable) <- MyCounTable[,1]
MyCounTable[,1] <- NULL

#with DESeq2, a pseudocount of 1/2 is added when there's full 0 for a comparison, but it's still better to gain time
############## removing IDs with full 0 count #############################
read_threshold<-opt$read_threshold


filter_list<-list()
for(i in 1:ncol(MyCounTable)){
  
  if(i<ncol(MyCounTable)){
    one_filter<-print(paste("MyCounTable","[",i,"]"," > ",read_threshold," |",sep=""))
    
    filter_list[[length(filter_list)+1]]<-one_filter
  }
  
  if(i==ncol(MyCounTable)){
    one_filter<-print(paste("MyCounTable","[",i,"]"," > ",read_threshold,sep=""))
    
    filter_list[[length(filter_list)+1]]<-one_filter
    
    cat(paste("MyCounTable",
              "<-",
              "MyCounTable",
              "[which(",sep=""),
        file="filtering_before_diff.txt",append=F)
    
    lapply(filter_list, write, "filtering_before_diff.txt", append=TRUE) 
    
    cat("),]",file="filtering_before_diff.txt",append=T)
  }
  
}

source("filtering_before_diff.txt")


##############################################################


################################################################################################

######## adding pseudocount +1 to avoid division by 0   #################################

#note :not needed normally anymore, DESeq can automatically add a pseudocount 1/2 to avoid -inf/+inf results (remark, the log2FC will be ridiculously high though), but the padj can be NA if independent filtering=T
MyCounTable<-MyCounTable+1

############ type of comparisons #########

##for all pairwise comparisons
if(type=="all"){
  comparisons_matrix<-NULL
  comparisons<-as.list(UniqueConditions)
  comparisons_matrix<-combn(comparisons,2)
  
}

##for comparisons between the first cond and all others
if(type=="first_vs_others"){
  n=2
  comparisons=as.list(UniqueConditions)
  comparisons_matrix<-NULL
  for(first_cond_pos in 1:(length(comparisons)-1)){  
    
    first_cond<-comparisons[[1]]
    
    second_cond_pos<-first_cond_pos+n
    
    second_cond<-comparisons[[n]]
    
    frame.tmp<-as.matrix(t(data.frame(first_cond,second_cond)))
    row.names(frame.tmp)<-NULL
    comparisons_matrix<-cbind(comparisons_matrix,frame.tmp)
    
    n=n+1
  }
}

if(type=="two_firsts_vs_others_combined"){
  
  
  n=2
  comparisons=as.list(UniqueConditions)
  comparisons_matrix<-NULL
  for(first_cond_pos in 1:(length(comparisons)-1)){  
    
    first_cond<-comparisons[[1]]
    
    second_cond<-comparisons[[n]]
    
    frame.tmp<-as.matrix(t(data.frame(first_cond,second_cond)))
    row.names(frame.tmp)<-NULL
    comparisons_matrix<-cbind(comparisons_matrix,frame.tmp)
    
    n=n+1
  }
  
  group1<-comparisons_matrix[,1]
  group2<-comparisons_matrix[2,2:ncol(comparisons_matrix)]
  
  comparisons_matrix<-matrix(c(paste(group1,collapse="_"),paste(group2,collapse="_")))
  
  cat("group1 : ",group1,"\n")
  cat("group2 : ",group2,"\n")

}

if(type=="first_vs_others_combined"){
  
  
  n=2
  comparisons=as.list(UniqueConditions)
  comparisons_matrix<-NULL
  for(first_cond_pos in 1:(length(comparisons)-1)){  
    
    first_cond<-comparisons[[1]]
    
    second_cond<-comparisons[[n]]
    
    frame.tmp<-as.matrix(t(data.frame(first_cond,second_cond)))
    row.names(frame.tmp)<-NULL
    comparisons_matrix<-cbind(comparisons_matrix,frame.tmp)
    
    n=n+1
  }
  
  group1<-comparisons_matrix[1,1]
  group2<-comparisons_matrix[2,1:ncol(comparisons_matrix)]
  
  
  comparisons_matrix<-matrix(c(paste(group1,collapse="_"),paste(group2,collapse="_")))
  
}
###############################################


########### DESeq2 functions ###################

#samples <-data.frame(row.names=names(MyCounTable),condition=condition,treatment=treatment)

#it's better to split DESq(), because we can control each part of this "multifunction"
samples <- data.frame(row.names=names(MyCounTable),condition=condition)

#dds<-DESeqDataSetFromMatrix(countData=as.matrix(MyCounTable),colData=samples,design=~treatment+condition)

dds<-DESeqDataSetFromMatrix(countData=as.matrix(MyCounTable),colData=samples,design=~condition)

rm(MyCounTable);gc()



if(given_size_factor!=""){
  
  cat("we're going to use our own size factors :\n",split_given_size_factor,"\n")
  
  cat("\n")
  
  sizeFactors(dds)<-split_given_size_factor
  
  
}else{
  
  dds<-estimateSizeFactors(dds)
  
}




#dds<-estimateDispersions(dds)

getDispersions<-function(my_object=""){
  
  dds<-try(estimateDispersions(my_object))
  
  if (class(dds)=="try-error"){
    
    cat("with fitType='parametric', the dispersion trend was not well captured by the function, we will use instead fitType='mean'")
    
    dds<-try(estimateDispersions(my_object,fitType="mean"))
    
  }
  
  if (class(dds)=="try-error"){
    
    cat("all gene-wise dispersion estimates are within 2 orders of magnitude")
    
    dds<-estimateDispersionsGeneEst(my_object)
    
    dispersions(dds)<-mcols(dds)$dispGeneEst
    
  }
  
  return(dds)
}

dds<-getDispersions(dds)

#https://support.bioconductor.org/p/77021/
#betaPrior=FALSE

if(type!="two_firsts_vs_others_combined" & type!= "first_vs_others_combined"){
  
  dds<-nbinomWaldTest(dds)

}else{
  
  dds<-nbinomWaldTest(dds,betaPrior=TRUE)
  
}

cat(colnames(dds),"\n",file="size_factors.tsv")
cat(sizeFactors(dds),"\n",file="size_factors.tsv",append=T)

#checking dist between libraries
pdf("clustering_of_samples.pdf",width=20,height=15)

  
  
  rld<-tryCatch(vst(dds),error=function(e) e, warning=function(w) w)
  
  if(any(class(rld)%in%c("error","warning"))){
    
    cat("we cannot use vst() function, let's try varianceStabilizingTransformation()\n")
    
    rld<-varianceStabilizingTransformation(dds)
    
  }
  
  sampleDists<-dist(t(assay(rld) ) )
  sampleDistMatrix<-as.matrix( sampleDists )
  rownames(sampleDistMatrix)<-colnames(rld)
  colnames(sampleDistMatrix)<-colnames(rld)
  colours=colorRampPalette(rev(brewer.pal(9,"Blues")) )(255)
  annotdf<-data.frame(row.names=colnames(rld),condition=rld$condition)
  # heatmap.2(sampleDistMatrix,trace="none",col=colours,main="clustering of samples",cexRow=1,
  #           cexCol=1,srtRow=45,srtCol=45)
  gg_color_hue <- function(n) {hues = seq(15, 375, length = n + 1);hcl(h = hues, l = 65, c = 100)[1:n]}
  
  #my_conds<-rainbow(length(unique(rld$condition)))
  my_conds<-gg_color_hue(length(unique(rld$condition)))
  names(my_conds)<-as.character(unique(rld$condition))
  
  
  p<-pheatmap(sampleDistMatrix,
              clustering_distance_rows=sampleDists,
              clustering_distance_cols=sampleDists,
              #"ward.D2"
              #clustering_method = "ward.D2",
              col=colours,
              annotation_col = annotdf,
              annotation_row = annotdf,
              annotation_colors=list(condition=my_conds),
              legend=T,
              annotation_names_row =F,
              annotation_names_col =F,
              main="clustering of samples",silent = TRUE,
              #legend_breaks = my_range,
              #legend_labels = c(my_range[2:length(my_range)], "Euclidean dist\n"),
              fontsize =16)
  
  grid.draw(p$gtable)
  
  data<-plotPCA(rld,ntop=nrow(rld),returnData=TRUE)
  
  percentVar<-round(100 * attr(data, "percentVar"))
  
  ggplot(data,aes(PC1,PC2,color=condition))+geom_point()+geom_text(aes(label=name),hjust=0,vjust=0,size=5)+
    xlab(paste0("PC1 : ",percentVar[1],"% variance")) +
    ylab(paste0("PC2 : ",percentVar[2],"% variance"))+
    theme(axis.text.x=element_text(color="black",size=14,hjust=1,face="bold"),axis.text.y = element_text(color="black",size=14,face="bold"),plot.title = element_text(hjust = 0.5,size=16),axis.title=element_text(size=16,face="bold"),legend.text=element_text(size=20),legend.title=element_text(size=16))+scale_color_manual(values = my_conds)


dev.off()

##############################

#initializing data frames/list for further analysis
log2FC_allcomparisons<-list()
padj_allcomparisons<-list()
SigDiffAllComparisons<-data.frame()
score_list<-list()
padj_threshold<-opt$padj_threshold
mean_threshold<-opt$mean_threshold
log2FC_threshold<-opt$log2FC_threshold

white_background<-theme(axis.line = element_line(colour = "black"),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank())

#h vs i
for (j in 1:ncol(comparisons_matrix)){
  
  h<-as.character(comparisons_matrix[1,j])
  i<-as.character(comparisons_matrix[2,j])
  
  official_annotation<-official_annotation_backup
  
  
  one_comparison_path<-paste(home,h,"_vs_",i,"/",sep="")
  dir.create(one_comparison_path, showWarnings = F, recursive = F)
  
  
  print(h)
  print(i)
  print("----")
  
  #1.4.3 http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf
  #with independentFiltering=T (default), IDs with very low counts+high variance (mix of both, or just low counts),  will be discarded, and the padj will be NA
  #with cooksCutoff=F, IDs with high counts and with outliers won't have NA as padj
  
  if(type!="two_firsts_vs_others_combined" & type!="first_vs_others_combined"){
  
  #be carefull ! contrary to DESeq1, here it's 1st arg over 2nd arg, so it's i/h
  result<-results(dds,independentFiltering=F,cooksCutoff=F,contrast=c("condition",i,h),parallel=TRUE,BPPARAM=MulticoreParam(4))
  
  }else{
   
    #independentFiltering=F,cooksCutoff=F 
  result<-results(dds,independentFiltering=F,cooksCutoff=F,contrast=list(paste("condition",group2,sep=""),paste("condition",group1,sep="")), listValues=c(1/length(group2),-1/length(group1)),parallel=TRUE,BPPARAM=MulticoreParam(4))
    
  }
  
  
  #https://www.biostars.org/p/282295/
  #pdf(file=paste(one_comparison_path,"volcanoplot","_",h,"_vs_",i,".pdf",sep=""))

  #    par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)
  #   #
  #    topT <- as.data.frame(result)
  #   #
  #   # #Adjusted P values (FDR Q values)
  #    with(topT, plot(log2FoldChange, -log10(padj), pch=19, main=paste("volcano plot ",h," vs ",i,"\ntotal genes : ",nrow(result),sep=""), cex=1.0, xlab="log2FoldChange", ylab="-log10(padj)"))
  #   #
  #    with(subset(topT, padj<padj_threshold), points(log2FoldChange, -log10(padj), pch=19, col="red", cex=1))
  # 
  # dev.off()
  
all_genes <- result[order(result$padj),]
  
all_genes = as.data.frame(dplyr::mutate(as.data.frame(all_genes), significance=ifelse(all_genes$padj<=padj_threshold, "sig. diff", "not_significant"),id=row.names(all_genes)),row.names=row.names(all_genes))


  #MA plot
  #pdf(file=paste(one_comparison_path,"MAplot","_",h,"_vs_",i,".pdf",sep=""))
  
  png(filename=paste(one_comparison_path,"MAplot","_",h,"_vs_",i,".png",sep=""),width=1200,height=1000)
  
    #mean expression=mean of all normalized reads for one ID in the comparison
    
    plotMA(result,alpha = padj_threshold,ylim=c(-5,5),main=paste("MAplot ",h," vs ",i,"\ntotal genes : ",nrow(result),"\nnb differential : ",nrow(all_genes[which(all_genes$padj<=padj_threshold),]),sep=""),cex=1.2,text.cex = 2)
    #DESeq2::plotMA(result,ylim=c(-5,5),main=paste("MAplot ",h," vs ",i,"\ntotal genes : ",nrow(result),"\nnb differential : ",nrow(all_genes[which(all_genes$padj<=padj_threshold),]),sep=""),cex=1.2,text.cex = 2,col = ifelse(x$padj>=padj_threshold, "black", "red"))
  #text.cex = 1
  
  dev.off()
  
  cairo_ps(filename=paste(one_comparison_path,"MAplot","_",h,"_vs_",i,".eps",sep=""),width=10,height=10)
  
    plotMA(result,alpha = padj_threshold,ylim=c(-5,5),main=paste("MAplot ",h," vs ",i,"\ntotal genes : ",nrow(result),"\nnb differential : ",nrow(all_genes[which(all_genes$padj<=padj_threshold),]),sep=""),cex=1.2,text.cex = 2)
  
  dev.off()
  
  #system(paste("convert ",paste(one_comparison_path,"MAplot","_",h,"_vs_",i,".png",sep="")," eps3:",paste(one_comparison_path,"MAplot","_",h,"_vs_",i,".eps",sep=""),sep=""))
  
  
  ## this part is specific to mouse
 #  y_genes<-as.character(unlist(read.delim("/home/marcgabriel/Documents/mouse_gencode.vM18/Y_genes.tsv",header=F,sep="\t")))
 #  
 #  x_genes<-as.character(unlist(read.delim("/home/marcgabriel/Documents/mouse_gencode.vM18/X_genes.tsv",header=F,sep="\t")))
 #  
 #  
 #  for(one_line in 1:nrow(all_genes)){
 #    
 #    if(as.character(all_genes$id[one_line])%in%y_genes & all_genes$padj[one_line]<=padj_threshold){
 #      
 #      all_genes$significance[one_line]<-"sig. diff & chrY"
 #      
 #    }else if(as.character(all_genes$id[one_line])%in%x_genes & all_genes$padj[one_line]<=padj_threshold){
 #      
 #      all_genes$significance[one_line]<-"sig. diff & chrX"
 #      
 #    }
 #    
 #  }
 #  
 #  #http://www.sthda.com/french/wiki/ggplot2-types-de-points-logiciel-r-et-visualisation-de-donnees
 #  p = ggplot2::ggplot(all_genes, ggplot2::aes(log2FoldChange, -log10(padj))) +
 #    ggplot2::geom_point(ggplot2::aes(col = significance,size=significance)) +
 #    ggplot2::scale_color_manual(values = c("sig. diff"="red", "not_significant"="black","sig. diff & chrY"="cornflowerblue","sig. diff & chrX"="#00A86B")) +
 #    #scale_shape_manual(values=c("female_specific"=0, "male_specific"=4,"unspecific"=19),name="gender specificity")+
 #    scale_size_manual(values=c("sig. diff & chrY"=6, "sig. diff & chrX"=6,"sig. diff"=6,"not_significant"=6),name="significance")+
 #    geom_vline(xintercept=0,color="grey20")+
 #    ggplot2::ggtitle(paste("volcano plot ",h," vs ",i,"\ntotal genes : ",nrow(all_genes),"\nnb differential : ",nrow(all_genes[which(all_genes$padj<=0.05),]),sep=""))
 #  
 #  #ggrepel::geom_text_repel(data=all_genes[1:10, ], ggplot2::aes(label=rownames(all_genes[1:10,])))+
 #  p<-p +theme_bw()+
 #    theme(axis.text.x=element_text(size=34,hjust=1,color="black"),axis.text.y = element_text(color="black",size=34),plot.title = element_text(hjust = 0.5,size=22,face="bold"),axis.title=element_text(size=34,face="bold"),legend.text=element_text(size=34),legend.title=element_text(size=34,face="bold"))+white_background
 #    #scale_y_log10()
 #    #panel.border=element_rect(colour="black",size=1)
 #    #ggrepel::geom_text_repel(data=results[1:10, ], ggplot2::aes(label=rownames(results[1:10, ])))
 #    
 # png(filename=paste(one_comparison_path,"volcanoplot","_",h,"_vs_",i,".png",sep=""),width=1200,height=1000)
 # 
 #    print(p)
 #  
 #  dev.off()
  
  
  
  ## end of specificity to mouse
  
   p = ggplot2::ggplot(all_genes, ggplot2::aes(log2FoldChange, -log10(padj))) +
     ggplot2::geom_point(ggplot2::aes(col = significance,size=significance)) +
     ggplot2::scale_color_manual(values = c("sig. diff"="red", "not_significant"="black")) +
     #scale_shape_manual(values=c("female_specific"=0, "male_specific"=4,"unspecific"=19),name="gender specificity")+
     scale_size_manual(values=c("sig. diff"=6,"not_significant"=6),name="significance")+
     geom_vline(xintercept=0,color="grey20")+
     ggplot2::ggtitle(paste("volcano plot ",h," vs ",i,"\ntotal genes : ",nrow(all_genes),"\nnb differential : ",nrow(all_genes[which(all_genes$padj<=padj_threshold),]),sep=""))

   #ggrepel::geom_text_repel(data=all_genes[1:10, ], ggplot2::aes(label=rownames(all_genes[1:10,])))+
   p<-p +theme_bw()+
     theme(axis.text.x=element_text(size=34,hjust=1,color="black"),axis.text.y = element_text(color="black",size=34),plot.title = element_text(hjust = 0.5,size=22,face="bold"),axis.title=element_text(size=34,face="bold"),legend.text=element_text(size=34),legend.title=element_text(size=34,face="bold"))+white_background
     #scale_y_log10()
     #panel.border=element_rect(colour="black",size=1)
     #ggrepel::geom_text_repel(data=results[1:10, ], ggplot2::aes(label=rownames(results[1:10, ])))

  png(filename=paste(one_comparison_path,"volcanoplot","_",h,"_vs_",i,".png",sep=""),width=1200,height=1000)

     print(p)

   dev.off()
   
   rm(all_genes);gc()
  
  
  #normalized counts 
  NormCount<-as.data.frame(counts(dds,normalized=TRUE ))
  
  raw_counts<-as.data.frame(counts(dds,normalized=FALSE ))
  
  #writing in a file normalized counts
  normalized_counts<-data.frame(id=row.names(NormCount),NormCount,row.names=NULL)
  write.table(normalized_counts,file="normalized_counts.tsv", sep='\t',row.names=F, col.names=T, quote=F)
  
  
  #normalized counts for one compararison with all replicates
  counts_list<-list()
  
  #don't forget : grepl take for example col_t1 and col_t15, we're using underscore (= replicates) just after to distinguish them !!
  #note : be carefull, sometimes there's no underscore between the condition name and the term specific to the replicate  !
  
  separator="_"
  
  h_back<-h;i_back<-i
  my_regex="[0-9]+"
  
  if(type=="two_firsts_vs_others_combined" | type=="first_vs_others_combined"){
    
    h<-paste(paste(paste(group1,'_',sep=""),collapse="[0-9]+|"),my_regex,sep="")
    
    i<-paste(paste(paste(group2,'_',sep=""),collapse="[0-9]+|"),my_regex,sep="")
    
    my_regex=""
    
    separator=""
    
  }
  
  
  l=0
  for (k in 1:length(NormCount)){
    
    
    
    #if condition name is there for one comparison, report the replicates
    if (any(grepl(paste(h,separator,my_regex,sep=""),names(NormCount[k]), perl=TRUE))==TRUE){   
      
      #rep counts cond1 & rep counts cond 1
      counts_list[[length(counts_list)+1]]<-as.data.frame(NormCount[k])
      
      l=l+1
    }
    
  }
  
  for (k in 1:length(NormCount)){
    
    #if condition name is there for the other one comparison, report the replicates
    if(any(grepl(paste(i,separator,my_regex,sep=""),names(NormCount[k]), perl=TRUE)==TRUE)){
      
      #rep counts cond2 & rep counts cond 2
      counts_list[[length(counts_list)+1]]<-as.data.frame(NormCount[k])   
      
    }         
  }
  
  #merge a list of dataframe by col (note : reduce then merge didn't work !)
  #counts rep cond1 and counts rep cond2 are side by side
  new_NormCount <- do.call("cbind",counts_list)
  
  #add Mean of norm. counts, DESeq2 doesn't put them
  mean_1st_cond_name=paste("Mean_",h_back,sep="")
  
  median_1st_cond_name=paste("Median_",h_back,sep="")
  
  mean_1st_cond=data.frame(apply(new_NormCount[1:l],1,mean))
  names(mean_1st_cond)<-mean_1st_cond_name
  
  median_1st_cond=data.frame(apply(new_NormCount[1:l],1,median))
  names(median_1st_cond)<-median_1st_cond_name
  
  mean_2st_cond_name=paste("Mean_",i_back,sep="")
  
  median_2st_cond_name=paste("Median_",i_back,sep="")
  
  mean_2st_cond=data.frame(apply(new_NormCount[(l+1):length(new_NormCount)],1,mean))
  names(mean_2st_cond)<-mean_2st_cond_name
  
  median_2st_cond=data.frame(apply(new_NormCount[(l+1):length(new_NormCount)],1,median))
  names(median_2st_cond)<-median_2st_cond_name
  
  #convert log2FC to FC, DESeq2 has removed them
  result$FoldChange<-2^(result$log2FoldChange)
  
  #construction of a dataframe with gene id, rep cond1, rep cond2, and the rest of DESeq columns
  new_col<-result[,c("FoldChange","log2FoldChange","pvalue","padj")]
  
  new_result <- data.frame(id=row.names(new_col),new_NormCount,mean_1st_cond,mean_2st_cond,median_1st_cond,median_2st_cond,new_col,row.names=NULL)
  
  ################
  
  
  my_regex="[0-9]+"
  
  if(type=="two_firsts_vs_others_combined" | type=="first_vs_others_combined"){
    
    h<-paste(paste(paste(group1,'_',sep=""),collapse="[0-9]+|"),my_regex,sep="")
    
    i<-paste(paste(paste(group2,'_',sep=""),collapse="[0-9]+|"),my_regex,sep="")
    
    my_regex=""
    
    separator=""
    
  }
  
  # counts_list_raw<-list()
  # 
  # l=0
  # for (k in 1:length(raw_counts)){
  #   
  #   
  #   #if condition name is there for one comparison, report the replicates
  #   if (any(grepl(paste(h,separator,my_regex,sep=""),names(raw_counts[k]), perl=TRUE))==TRUE){   
  #     
  #     #rep counts cond1 & rep counts cond 1
  #     counts_list_raw[[length(counts_list_raw)+1]]<-as.data.frame(raw_counts[k])
  #     
  #     l=l+1
  #   }
  #   
  # }
  # 
  # for (k in 1:length(raw_counts)){
  #   
  #   #if condition name is there for the other one comparison, report the replicates
  #   if(any(grepl(paste(i,separator,my_regex,sep=""),names(raw_counts[k]), perl=TRUE)==TRUE)){
  #     
  #     #rep counts cond2 & rep counts cond 2
  #     counts_list_raw[[length(counts_list_raw)+1]]<-as.data.frame(raw_counts[k])   
  #     
  #   }         
  # }
  # 
  # new_raw_counts <- do.call("cbind",counts_list_raw)
  # 
  # raw_result<-data.frame(id=row.names(new_col),new_raw_counts,row.names=NULL)
  
  
  
  print("here 1")
  
  ################
  
  new_result$new_name<-sub("\\.[0-9]+.*","",new_result$id)
  
  print("here 2")
  
  official_annotation$new_name<-sub("\\.[0-9]+.*","",official_annotation$gene_id)
  
  
  official_annotation<-subset(official_annotation, select = -c(gene_id))
  
  
  #if no annotation, deactivate this part and add this instead :  results_annotated<-new_result !!!
  results_annotated<-merge(new_result,official_annotation, 
                           by.x="new_name",
                           by.y="new_name",
                           all.x=TRUE,
                           all.y=FALSE)
  
  print("here 3")
  
  new_result<-subset(new_result, select = -c(new_name))
  
  
  
  
  
  if(additional_info!="none"){
    
    additional_info_list<-as.list((unlist(strsplit(additional_info,","))))
    
    results_annotated$new_name<-sub("\\.[0-9]+.*","",results_annotated$id)
    
    for(one_info in additional_info_list){
      
      additional_info_table<-read.delim(file=one_info,sep="\t",check.names=F,header=T)
      
      additional_info_table$new_name<-sub("\\.[0-9]+.*","",additional_info_table[,1])
      
      
      additional_info_table<-additional_info_table[,c(2:ncol(additional_info_table))]
      
      #names(additional_info_table)<-c("id","type")
      
      #colnames(additional_info_table)[1]
      results_annotated<-merge(results_annotated,additional_info_table, 
                               by.x="new_name",
                               by.y="new_name",
                               all.x=TRUE,
                               all.y=FALSE)
      
      
    }
    
    
    
  }
  
  results_annotated<-subset(results_annotated, select = -c(new_name))
  
  #put as first columns gene_id and gene_name
  results_annotated<-results_annotated[,c(1,ncol(new_result)+1,(ncol(new_result)+2):ncol(results_annotated),2:ncol(new_result))]
  
  if(insert_raw_counts!="no"){
    
  
      results_annotated$new_name<-sub("\\.[0-9]+.*","",results_annotated$id)
      
      names(raw_result)[2:ncol(raw_result)]<-paste(names(raw_result)[2:ncol(raw_result)],"_raw",sep="")
      
      raw_result$new_name<-sub("\\.[0-9]+.*","",raw_result$id)
      
      raw_result<-subset(raw_result, select = -c(id))
      
      results_annotated<-merge(results_annotated,raw_result, 
                               by.x="new_name",
                               by.y="new_name",
                               all.x=TRUE,
                               all.y=FALSE)
      
      results_annotated<-subset(results_annotated, select = -c(new_name))
  
}
  
  #function to get the roc curve values
  getRocCurve<-function(states="",counts=""){
    
    my_function<-function(states1="",counts1=""){
      
      return(roc(states1,counts1,percent=TRUE,ci=TRUE, boot.n=100, ci.alpha=0.9,stratified=FALSE,smooth=T))
      
    }
    
    my_roc_result<-tryCatch(my_function(states,counts),error=function(e) e, warning=function(w) w)
    
    if(any(class(my_roc_result)%in%c("error","warning"))){
      
      cat("we cannot draw the smoothed roc curve, let's try the raw curve\n")
      
      my_roc_result<-roc(states,counts,percent=TRUE,ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE)
      
    }
    
    return(my_roc_result)
  }
  
  if(compute_AUC=="yes"){
    
    print(paste("time to compute AUC for : ",h_back," vs ",i_back," ...",sep=""))
    
    rep_num_cond1<-length(grep(paste(h,'_',sep=""),names(normalized_counts), perl=TRUE,value=T))
    
    rep_num_cond2<-length(grep(paste(i,'_',sep=""),names(normalized_counts), perl=TRUE,value=T))
    
    registerDoParallel(cores=nb_cores)
    
    
    system.time(frame_of_AUC<-foreach(one_row=1:nrow(normalized_counts), .combine=rbind, .inorder=FALSE) %dopar% {
      
      roc1 <-getRocCurve(states=c(rep(h,rep_num_cond1),rep(one_row,rep_num_cond2)),counts=unlist(normalized_counts[one_row,c(grep(paste(h,'_',sep=""),names(normalized_counts), perl=TRUE,value=T),grep(paste(i,'_',sep=""),names(normalized_counts), perl=TRUE,value=T))]))
      
      data.frame(id=as.character((normalized_counts[one_row,1])),
                 AUC_value=round(as.double(roc1$auc),2),
                 AUC_CI_low=round(roc1$ci[1],digits=2),
                 AUC_CI_high=round(roc1$ci[3],digits=2))
      
    })
    
    results_annotated<-merge(results_annotated,frame_of_AUC, 
                             by.x="id",
                             by.y=colnames(frame_of_AUC)[1],
                             all.x=TRUE,
                             all.y=FALSE)
    
  }
  
  
  comparison_name=paste("DESeq_output_",h_back,"_vs_",i_back,".tsv",sep="")
  
  write.table(results_annotated,file=paste(one_comparison_path,comparison_name,sep=""), sep='\t',row.names=F, col.names=T, quote=F)
  
  print("here 6")
  
  ############### extracting log2FC ############################
  
  # log2FC_one_comparison<-results_annotated[,c("id","log2FoldChange")]
  # names(log2FC_one_comparison)[2]<-paste("log2FC_",h_back,"_vs_",i_back,sep="")
  # 
  # padj_one_comparison<-results_annotated[,c("id","padj")]
  # names(padj_one_comparison)[2]<-paste("padj_",h_back,"_vs_",i_back,sep="")
  
  #you can convert this list in data frame, and use it later if needed (merging it with the sigdif dataframe for example)
  #log2FC_allcomparisons[[length(log2FC_allcomparisons)+1]]<-log2FC_one_comparison
  
  #padj_allcomparisons[[length(padj_allcomparisons)+1]]<-padj_one_comparison
  
  print("here 7")
  
  
  ###############################################################
  
  
  ############ score computing ##########################
  # results_annotated_scored<-results_annotated
  # 
  # results_annotated_scored$score<-""
  # 
  # 
  # 
  # #from the original result annotated, give a score for up and down features (1= up ; down = 2, not significant =0)
  # for(l in 1:nrow(results_annotated_scored)){
  #   
  #   
  #   #you can have NA's in DESeq2, if full 0 for a cond, or if independent filtering was T (IDs discarded because of low counts)
  #   if(is.na(results_annotated_scored[l,"padj"])==T){
  #     
  #     results_annotated_scored$score[l]<-sub("$",0,results_annotated_scored$score[l] )
  #     
  #     next
  #   }
  #   
  #   
  #   if(results_annotated_scored[l,"log2FoldChange"] >=log2FC_threshold & is.na(results_annotated_scored[l,"padj"])==F){
  #     
  #     if(results_annotated_scored[l,"padj"] <padj_threshold){
  #       
  #       results_annotated_scored$score[l]<-sub("$",1,results_annotated_scored$score[l] )
  #       
  #       
  #     }
  #     
  #     
  #     
  #     
  #   }
  #   
  #   if(results_annotated_scored[l,"log2FoldChange"] <log2FC_threshold &  is.na(results_annotated_scored[l,"padj"])==F){
  #     
  #     if(results_annotated_scored[l,"padj"] <padj_threshold){
  #       
  #       results_annotated_scored$score[l]<-sub("$",2,results_annotated_scored$score[l] )
  #       
  #     }
  #     
  #     
  #   }
  #   
  #   if(results_annotated_scored[l,"padj"] >padj_threshold){
  #     
  #     results_annotated_scored$score[l]<-sub("$",0,results_annotated_scored$score[l] )
  #   }
  #   
  # }
  # 
  # scores<-results_annotated_scored[,c("id","score")]
  # names(scores)[2]<-paste(h_back,"_vs_",i_back,"_","score",sep="")
  
  #you can convert this list in data frame, and use it later if needed (merging it with the sigdif dataframe for example)
  #score_list[[length(score_list)+1]]<-scores
  
  
  ####################  selecting sig diff results ########################################
  
  sig_diff_results_annotated<- results_annotated[which(results_annotated$padj <= padj_threshold & abs(results_annotated$log2FoldChange)>=log2FC_threshold),]
  
  #dev
  mean_names<-names(head(sig_diff_results_annotated[grep("Mean_",names(sig_diff_results_annotated),value=T)]))
  
  sig_diff_results_annotated<-sig_diff_results_annotated[which(sig_diff_results_annotated[mean_names[1]]>=mean_threshold | sig_diff_results_annotated[mean_names[2]]>=mean_threshold),]
  
  #end of dev
  
  cat(paste("number of diff genes : ",nrow(sig_diff_results_annotated),"\n",sep=""))
  
  
  ######################## sig diff IDs across all comparisons ########################
  #SigDiffAllComparisons<-unique(rbind(SigDiffAllComparisons,data.frame(sig_diff_results_annotated[,c("id","gene_name","chromosome","start","end","strand")])))
  ####################################################################################
  
  print("here 8")
  
  
  diff=paste("sig_diff_",h_back,"_vs_",i_back,".tsv",sep="")
  
  write.table(sig_diff_results_annotated,file=paste(one_comparison_path,diff,sep=""), sep='\t',row.names=F, col.names=T, quote=F)
  
  sig_feature_up_name=paste("sig_diff_upregulated_",h_back,"_vs_",i_back,".tsv",sep="")
  
  sig_upregulated_feature<-sig_diff_results_annotated[which(sig_diff_results_annotated$log2FoldChange >0),] 
  
  write.table(sig_upregulated_feature,file=paste(one_comparison_path,sig_feature_up_name,sep=""), sep='\t',row.names=F, col.names=T, quote=F)
  
  ###dev
  
  #make GO analysis on all diff genes
  if(nrow(sig_diff_results_annotated)!=0){
    
    converted_IDs<-tryCatch(bitr(as.character(sig_diff_results_annotated$gene_name), fromType = "SYMBOL",
                                 toType = c("ENSEMBL", "SYMBOL","ENTREZID"),
                                 OrgDb = "org.Hs.eg.db"),
                            error=function(e) e,
                            warning=function(w) w)
    
    
    if(any(class(converted_IDs)%in%c("error"))){
      
      cat("none of the IDs are recognized (custom IDs ?)\n")
      
      ego <-data.frame()
      
    }else{
      
      converted_IDs<-bitr(as.character(sig_diff_results_annotated$gene_name), fromType = "SYMBOL",
                          toType = c("ENSEMBL", "SYMBOL","ENTREZID"),
                          OrgDb = "org.Hs.eg.db")
      
      ego <- enrichGO(gene          = converted_IDs$ENTREZID,
                      OrgDb         = "org.Hs.eg.db",
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)
      
    }
    
    
    if(length(ego)!=0){
      
      if(nrow(ego)!=0){
        
        
        write.table(ego,file=paste(one_comparison_path,"sig_diff_all_",h_back,"_vs_",i_back,"_GO_terms.tsv",sep=""), sep='\t',row.names=F, col.names=T, quote=F)
        
        png(filename=paste(one_comparison_path,"sig_diff_all_",h_back,"_vs_",i_back,"_GO_terms.png",sep=""),width=1200,height=800)
        
        
        
        
        print(barplot(ego,colorBy = "p.adjust",font.size = 20,title=paste("GO enrichment for ",nrow(sig_diff_results_annotated)," all diff genes  (clusterProfiler package)",sep=""),showCategory=15)+ylab("nb genes")+scale_y_continuous(limits = c(0, getStepAndMax(max(as.data.frame(ego)$Count))[[2]]),breaks = seq(0,getStepAndMax(max(as.data.frame(ego)$Count))[[2]],getStepAndMax(max(as.data.frame(ego)$Count))[[1]]),expand=c(0,0))+
                scale_fill_continuous(name = "adjusted pvalue", guide = guide_colourbar(reverse=T))+theme(legend.text=element_text(size=16),legend.title=element_text(size=16),plot.title = element_text(size=18),axis.text.x=element_text(colour="black",angle = 45, hjust = 1,size=16)))
        
        dev.off()
        
        cairo_ps(filename=paste(one_comparison_path,"sig_diff_all_",h_back,"_vs_",i_back,"_GO_terms.eps",sep=""),width=20,height=10)
        
        print(barplot(ego,colorBy = "p.adjust",font.size = 20,title=paste("GO enrichment for ",nrow(sig_diff_results_annotated)," all diff genes (clusterProfiler package)",sep=""),showCategory=15)+ylab("nb genes")+scale_y_continuous(limits = c(0, getStepAndMax(max(as.data.frame(ego)$Count))[[2]]),breaks = seq(0,getStepAndMax(max(as.data.frame(ego)$Count))[[2]],getStepAndMax(max(as.data.frame(ego)$Count))[[1]]),expand=c(0,0))+
                scale_fill_continuous(name = "adjusted pvalue", guide = guide_colourbar(reverse=T))+theme(legend.text=element_text(size=16),legend.title=element_text(size=16),plot.title = element_text(size=18),axis.text.x=element_text(colour="black",angle = 45, hjust = 1,size=16)))
        
        dev.off()
        
      }
      
      
      
      
    }else{
      
      print(paste("-- no sig. GO terms --"))
    }
    
  }
  
  ##end of dev
  

  
  #make GO analysis on up genes
  if(nrow(sig_upregulated_feature)!=0){
    
    # converted_IDs<- bitr(as.character(sig_upregulated_feature$gene_name), fromType = "SYMBOL",
    #                      toType = c("ENSEMBL", "SYMBOL","ENTREZID"),
    #                      OrgDb = "org.Hs.eg.db")
    
    converted_IDs<-tryCatch(bitr(as.character(sig_upregulated_feature$gene_name), fromType = "SYMBOL",
                                 toType = c("ENSEMBL", "SYMBOL","ENTREZID"),
                                 OrgDb = "org.Hs.eg.db"),
                            error=function(e) e,
                            warning=function(w) w)
    
   
    if(any(class(converted_IDs)%in%c("error"))){
      
      cat("none of the IDs are recognized (custom IDs ?)\n")
      
      ego <-data.frame()
      
    }else{
      
      converted_IDs<-bitr(as.character(sig_upregulated_feature$gene_name), fromType = "SYMBOL",
                                   toType = c("ENSEMBL", "SYMBOL","ENTREZID"),
                                   OrgDb = "org.Hs.eg.db")
    
      ego <- enrichGO(gene          = converted_IDs$ENTREZID,
                      OrgDb         = "org.Hs.eg.db",
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)
      
    }
    
    
    if(length(ego)!=0){
      
      if(nrow(ego)!=0){
      
      
          write.table(ego,file=paste(one_comparison_path,"sig_diff_upregulated_",h_back,"_vs_",i_back,"_GO_terms.tsv",sep=""), sep='\t',row.names=F, col.names=T, quote=F)
          
          png(filename=paste(one_comparison_path,"sig_diff_upregulated_",h_back,"_vs_",i_back,"_GO_terms.png",sep=""),width=1200,height=800)
          
          
          
          
          print(barplot(ego,colorBy = "p.adjust",font.size = 20,title=paste("GO enrichment for ",nrow(sig_upregulated_feature)," upregulated genes (clusterProfiler package)",sep=""),showCategory=15)+ylab("nb genes")+scale_y_continuous(limits = c(0, getStepAndMax(max(as.data.frame(ego)$Count))[[2]]),breaks = seq(0,getStepAndMax(max(as.data.frame(ego)$Count))[[2]],getStepAndMax(max(as.data.frame(ego)$Count))[[1]]),expand=c(0,0))+
            scale_fill_continuous(name = "adjusted pvalue", guide = guide_colourbar(reverse=T))+theme(legend.text=element_text(size=16),legend.title=element_text(size=16),plot.title = element_text(size=18),axis.text.x=element_text(colour="black",angle = 45, hjust = 1,size=16)))
          
          dev.off()
          
          cairo_ps(filename=paste(one_comparison_path,"sig_diff_upregulated_",h_back,"_vs_",i_back,"_GO_terms.eps",sep=""),width=20,height=10)
          
            print(barplot(ego,colorBy = "p.adjust",font.size = 20,title=paste("GO enrichment for ",nrow(sig_upregulated_feature)," upregulated genes  (clusterProfiler package)",sep=""),showCategory=15)+ylab("nb genes")+scale_y_continuous(limits = c(0, getStepAndMax(max(as.data.frame(ego)$Count))[[2]]),breaks = seq(0,getStepAndMax(max(as.data.frame(ego)$Count))[[2]],getStepAndMax(max(as.data.frame(ego)$Count))[[1]]),expand=c(0,0))+
                    scale_fill_continuous(name = "adjusted pvalue", guide = guide_colourbar(reverse=T))+theme(legend.text=element_text(size=16),legend.title=element_text(size=16),plot.title = element_text(size=18),axis.text.x=element_text(colour="black",angle = 45, hjust = 1,size=16)))
            
          dev.off()
      
      }
      
      
      
      
    }else{
      
      print(paste("-- no sig. GO terms --"))
    }
  
  }
  
  
  
  #################################"
  
  
  
  sig_feature_down_name=paste("sig_diff_downregulated_",h_back,"_vs_",i_back,".tsv",sep="")
  
  sig_downregulated_feature<-sig_diff_results_annotated[which(sig_diff_results_annotated$log2FoldChange <0),]   
  
  write.table(sig_downregulated_feature,file=paste(one_comparison_path,sig_feature_down_name,sep=""), sep='\t',row.names=F, col.names=T, quote=F)
  
  #make GO analysis on down genes
  if(nrow(sig_downregulated_feature)!=0){
    
      # converted_IDs<- bitr(as.character(sig_downregulated_feature$gene_name), fromType = "SYMBOL",
      #                      toType = c("ENSEMBL", "SYMBOL","ENTREZID"),
      #                      OrgDb = "org.Hs.eg.db")
    
    
    converted_IDs<-tryCatch(bitr(as.character(sig_downregulated_feature$gene_name), fromType = "SYMBOL",
                                 toType = c("ENSEMBL", "SYMBOL","ENTREZID"),
                                 OrgDb = "org.Hs.eg.db"),
                            error=function(e) e,
                            warning=function(w) w)
    
     if(any(class(converted_IDs)%in%c("error"))){
      
        cat("none of the IDs are recognized (custom IDs ?)\n")
      
        ego <-data.frame()
      
      }else{
        
        converted_IDs<-bitr(as.character(sig_downregulated_feature$gene_name), fromType = "SYMBOL",
                                     toType = c("ENSEMBL", "SYMBOL","ENTREZID"),
                                     OrgDb = "org.Hs.eg.db")
      
        ego <- enrichGO(gene          = converted_IDs$ENTREZID,
                        OrgDb         = "org.Hs.eg.db",
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05,
                        readable      = TRUE)
        
        
      
      }
    
   
      
      if(length(ego)!=0){
        
        if(nrow(ego)!=0){
        
        
          write.table(ego,file=paste(one_comparison_path,"sig_diff_downregulated_",h_back,"_vs_",i_back,"_GO_terms.tsv",sep=""), sep='\t',row.names=F, col.names=T, quote=F)
          
          
          png(filename=paste(one_comparison_path,"sig_diff_downregulated_",h_back,"_vs_",i_back,"_GO_terms.png",sep=""),width=1200,height=800)
          
          print(barplot(ego,colorBy = "p.adjust",font.size = 20,title=paste("GO enrichment for ",nrow(sig_downregulated_feature)," downregulated genes  (clusterProfiler package)",sep=""),showCategory=15)+ylab("nb genes")+scale_y_continuous(limits = c(0,getStepAndMax(max(as.data.frame(ego)$Count))[[2]]),breaks = seq(0,getStepAndMax(max(as.data.frame(ego)$Count))[[2]],getStepAndMax(max(as.data.frame(ego)$Count))[[1]]),expand=c(0,0))+
            scale_fill_continuous(name = "adjusted pvalue", guide = guide_colourbar(reverse=T))+theme(legend.text=element_text(size=16),legend.title=element_text(size=16),plot.title = element_text(size=18),axis.text.x=element_text(colour="black",angle = 45, hjust = 1,size=16)))
          
          dev.off()
          
          cairo_ps(filename=paste(one_comparison_path,"sig_diff_downregulated_",h_back,"_vs_",i_back,"_GO_terms.eps",sep=""),width=20,height=10)
          
            print(barplot(ego,colorBy = "p.adjust",font.size = 20,title=paste("GO enrichment for ",nrow(sig_downregulated_feature)," downregulated genes  (clusterProfiler package)",sep=""),showCategory=15)+ylab("nb genes")+scale_y_continuous(limits = c(0, getStepAndMax(max(as.data.frame(ego)$Count))[[2]]),breaks = seq(0,getStepAndMax(max(as.data.frame(ego)$Count))[[2]],getStepAndMax(max(as.data.frame(ego)$Count))[[1]]),expand=c(0,0))+
                    scale_fill_continuous(name = "adjusted pvalue", guide = guide_colourbar(reverse=T))+theme(legend.text=element_text(size=16),legend.title=element_text(size=16),plot.title = element_text(size=18),axis.text.x=element_text(colour="black",angle = 45, hjust = 1,size=16)))
          
          dev.off()
        
        }
        
        
      }else{
        
        print(paste("-- no sig. GO terms --\n"))
      }
      
  }
  
  
  
  ########################
  
  
  
}

print("here 9")


#write.table(SigDiffAllComparisons,file="SigDiffAllComparisons.tsv", sep='\t',row.names=F, col.names=T, quote=F)
print("all sig. diff. genes across all comparisons are in a table !")


#uncomment this if you need these values
# log2FC_frame<-Reduce(function(x, y) merge(x, y,by=c("id"), all=TRUE),log2FC_allcomparisons)
# 
# write.table(log2FC_frame,file="log2FC_frame.tsv", sep='\t',row.names=F, col.names=T, quote=F)
# 
# print("all log2FC across all comparisons are in a table !")
# 
# padj_frame<-Reduce(function(x,y) merge(x, y,by=c("id"),all=TRUE),padj_allcomparisons)
# 
# write.table(padj_frame,file="padj_frame.tsv", sep='\t',row.names=F, col.names=T, quote=F)
# 
# print("all padj across all comparisons are in a table !")
# 
# score_frame<-Reduce(function(x, y) merge(x, y,by=c("id"),all=TRUE), score_list)
# 
# write.table(score_frame,file="score_frame.tsv", sep='\t',row.names=F, col.names=T, quote=F)
# 
# print("all scores across all comparisons are in a table !")

#dev




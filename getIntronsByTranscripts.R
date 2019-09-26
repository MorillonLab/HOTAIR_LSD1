#!/usr/bin/env Rscript

#check this script for this issue : gene_id & transcrit finishing with "_1" (or globally "_X"), are truncated from the "_" -> gencode27lift37 for example

args <- commandArgs(TRUE)

if (length(args)==0){
  stop("missing arguments !\nUsage : ./getIntronsByTranscripts.R <gff3 annotation> <yes/no> (if gene_name is present)")
}

suppressMessages(library(GenomicFeatures))
suppressMessages(library(data.table))

file<-args[1]

#file<-"/home/marcgabriel/Documents/gencode27/gencode.v27.annotation_sorted.gff3"

response=args[2]

options(ucscChromosomeNames=F)

txdb<-makeTxDbFromGFF(file,format="gff3")

IntronsByTranscripts<-unlist(intronsByTranscript(txdb,use.names=T))

if(length(IntronsByTranscripts)==0){
  
  q("no")
}

if(response=="yes"){
  
  
  gencode<-fread(paste("grep -P \"\\tgene\\t\" ",file,sep=""), header=F, stringsAsFactors=F,data.table=F,sep="\t", showProgress = FALSE)
  
}else{
  
  gencode<-fread(paste("grep -P \"\\tmRNA\\t\" ",file,sep=""), header=F, stringsAsFactors=F,data.table=F,sep="\t", showProgress = FALSE)
  
}

gencode<-gencode[,c(1,4,5,7,9)]

names(gencode)<-c("chr","start","end","strand","name")

if(response=="yes"){
  
  all_gene_IDs<-gsub("gene_id=","",grep("gene_id=",unlist(strsplit(as.character(gencode$name),";")),value = T)) 


}else{
  
  all_gene_IDs<-gsub("ID=","",grep("ID=",unlist(strsplit(as.character(gencode$name),";")),value = T))
  
}

if(response=="yes"){
  
  all_gene_names<-gsub("gene_name=","",grep("gene_name=",unlist(strsplit(as.character(gencode$name),";")),value = T))
  
  
}else{
  
  all_gene_names<-all_gene_IDs
  
}

link_gene_id_name<-data.frame(gene_id=all_gene_IDs,gene_name=all_gene_names)

if(response=="yes"){
  
  all_transcripts<-fread(paste("grep -P \"\\ttranscript\\t|\\tmRNA\\t\" ",file,"|cut -f9",sep=""), header=F, stringsAsFactors=F,data.table=F,sep="\t", showProgress = FALSE)
  
  all_transcripts_IDs<-gsub("transcript_id=","",grep("transcript_id=",unlist(strsplit(unlist(all_transcripts),";")),value = T))
  
  all_transcripts_names<-gsub("transcript_name=","",grep("transcript_name=",unlist(strsplit(unlist(all_transcripts),";")),value = T))
  
  all_gene_types<-gsub("gene_type=","",grep("gene_type=",unlist(strsplit(unlist(all_transcripts),";")),value = T))
  
  
  
}else{
  
  all_transcripts_IDs<-all_gene_IDs
  
  all_transcripts_names<-all_transcripts_IDs
  
}



link_transcript_id_name<-data.frame(transcript_id=all_transcripts_IDs,transcript_name=all_transcripts_names,gene_type=all_gene_types)


if(response=="yes"){
  
  transcriptsByGenes<-unlist(transcriptsBy(txdb,by="gene"))
  
  transcriptsByGenes<-data.frame(transcript_id=transcriptsByGenes$tx_name,gene_id=names(transcriptsByGenes))

}else{
  
  transcriptsByGenes<-data.frame(transcript_id=link_transcript_id_name$transcript_id,gene_id=link_transcript_id_name$transcript_name,gene_type=link_transcript_id_name$all_gene_types)
  
}

link_gene_id_name$fake_gene_id<-link_gene_id_name$gene_id

link_gene_id_name$gene_id<-gsub("_[0-9]+$","",link_gene_id_name$gene_id)

#add gene_name info
transcriptsByGenes<-merge(transcriptsByGenes,link_gene_id_name,
                  by.x="gene_id",
                  by.y="gene_id",
                  all.x=TRUE,
                  all.y=FALSE)

link_transcript_id_name$fake_transcript_id<-link_transcript_id_name$transcript_id

link_transcript_id_name$transcript_id<-gsub("_[0-9]+$","",link_transcript_id_name$transcript_id)


#add transcript_name info
transcriptsByGenes<-merge(transcriptsByGenes,link_transcript_id_name,
                          by.x="transcript_id",
                          by.y="transcript_id",
                          all.x=TRUE,
                          all.y=FALSE)

intronsGFF<-data.frame(transcript_id=names(IntronsByTranscripts),
                       chromosome=seqnames(IntronsByTranscripts),
                       source="annotation",
                       feature="intron",
                       start=start(IntronsByTranscripts),
                       end=end(IntronsByTranscripts),
                       empty1=".",
                       strand=strand(IntronsByTranscripts),
                       empty2=".",
                       ID=paste("ID=",names(IntronsByTranscripts),"_",start(IntronsByTranscripts),"_",end(IntronsByTranscripts),";intron_id=",names(IntronsByTranscripts),"_",start(IntronsByTranscripts),"_",end(IntronsByTranscripts),sep=""),
                       parent=paste(";Parent=",names(IntronsByTranscripts),sep=""),
                       transcript_id2=paste(";transcript_id=",names(IntronsByTranscripts),sep=""))

intronsGFF<-merge(intronsGFF,transcriptsByGenes, 
                  by.x="transcript_id",
                  by.y="transcript_id",
                  all.x=TRUE,
                  all.y=FALSE)

intronsGFF$transcript_id2<-paste(";transcript_id=",intronsGFF$fake_transcript_id,sep="")

intronsGFF$gene_id<-paste(";gene_id=",intronsGFF$fake_gene_id,sep="")

intronsGFF$gene_name<-paste(";gene_name=",intronsGFF$gene_name,sep="")

intronsGFF$transcript_name<-paste(";transcript_name=",intronsGFF$transcript_name,sep="")

intronsGFF$gene_type<-paste(";gene_type=",intronsGFF$gene_type,sep="")

intronsGFF<-intronsGFF[c("chromosome","source","feature","start","end","empty1","strand","empty2","ID","parent","transcript_name","gene_id","gene_name","transcript_id2","gene_type")]

#intronsGFF$transcript_id2
intronsGFF$attributes<-paste(intronsGFF$ID,intronsGFF$transcript_id2,intronsGFF$transcript_name,intronsGFF$parent,intronsGFF$gene_id,intronsGFF$gene_name,intronsGFF$gene_type,sep="")

intronsGFF<-subset(intronsGFF,select=-c(ID,parent,gene_id,transcript_id2,transcript_name,gene_id,gene_name,gene_type))

write.table(intronsGFF,stdout(), sep='\t',row.names=F,col.names=F,quote=F)



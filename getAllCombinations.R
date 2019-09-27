

########### variables ###
#########################

#sets
# set1<-c(1,2,3,4)
# 
# set2<-c(3,4,5,6)
# 
# set3<-c(4,5,6,7)
# 
# set4<-c(4,5,6,8)
# 
# all_sets<-list(set1,set2,set3,set4)
# 
# #name of the sets
# all_names<-list("set1","set2","set3","set4")




getAllConbinations<-function(all_sets,all_names){
  
    if(length(all_sets)!=length(all_names) | length(all_sets)==0 | length(all_names)==0){
      
      stop("number of names should be equal to number of sets !!")
      
    }
  
  
    #unique IDs for each set
    for( i in 1:length(all_sets)){
      
      all_sets[[i]]<-unique(all_sets[[i]])
      
    }
  
    #make a symmetric data frame with the IDs
    all_lengths<-c()
    
    for(i in 1:length(all_sets)){
      
      all_lengths<-c(all_lengths,length(all_sets[[i]]))
      
    }
    
    
    fake_frame<-data.frame(rep(NA,max(all_lengths)))
    
    
    for(i in 1:(length(all_sets)-1)){
      
      fake_frame<-cbind(fake_frame,data.frame(rep(NA,max(all_lengths))))
    }
    
    for(i in 1:(length(all_sets))){
      
      
      for(j in 1:length(all_sets[[i]])){
        
        fake_frame[j,i]<-all_sets[[i]][j]
        
      }
      
      
    }
    
    all_sets_frame<-fake_frame
    
    ######################
    
    
    ############ processing of variables #########
    ##############################################
    
    names(all_sets_frame)<-unlist(all_names)
    
    
    #make all possible combinations, until the length of the sets
    all_comb<-list()
    for(i in 2:length(all_sets)){
      
      
      all_comb[[length(all_comb)+1]]<-combn((1:length(all_sets)),i)
      
      
      
    }
    
    ###############################################
    
    
    pre_dataframe<-data.frame("")
    
    all_unique_counts<-list()
    
    all_unique_IDs<-list()
    
    
    #for each set, compute unique intersections & no intersections
    for(i in 1:(ncol(all_sets_frame))){
      
      unique_IDs<-setdiff(unlist(all_sets_frame[,i])[!is.na(unlist(all_sets_frame[,i]))],
                          unlist(all_sets_frame[,setdiff(1:(length(all_sets)),i)])[!is.na(unlist(all_sets_frame[,setdiff(1:(length(all_sets)),i)]))])
      
      cat("set : ",i,"\n")
      
      cat("\t - name of set : ",names(all_sets_frame[i]),"\n")
      
      #cat("\t - unique IDs : ",unique_IDs,"\n")
      
      cat("\t - length of unique IDs : ",length(unique_IDs),"\n")
      
      pre_dataframe<-cbind(pre_dataframe,data.frame(c(paste("cond",i,sep=""),names(all_sets_frame[i]))))
      
      all_unique_counts[[length(all_unique_counts)+1]]<-data.frame(c(paste("cond",i,"_spe_counts",sep=""),length(unique_IDs)))
      
      all_unique_IDs[[length(all_unique_IDs)+1]]<-data.frame(c(paste("cond",i,"_spe_IDs",sep=""),paste(unique_IDs,collapse="-_-")))
      
    }
    
    all_unique_counts<-as.data.frame(all_unique_counts)
    
    pre_dataframe<-pre_dataframe[,2:ncol(pre_dataframe)]
    
    pre_dataframe<-cbind(pre_dataframe,all_unique_counts)
    
    names(pre_dataframe)<-as.character(unlist(pre_dataframe[1,]))
    
    pre_dataframe<-pre_dataframe[2:nrow(pre_dataframe),]
    
    all_unique_IDs<-as.data.frame(all_unique_IDs)
    
    names(all_unique_IDs)<-as.character(unlist(all_unique_IDs[1,]))
    
    all_unique_IDs<-all_unique_IDs[2:nrow(all_unique_IDs),]
    
    for(i in 1:ncol(all_unique_IDs)){
      
      all_unique_IDs[,i]<-gsub("-_-"," ",all_unique_IDs[,i])
      
    }
    
    
    
    
    
    all_intersections_counts<-list()
    
    all_intersections_IDs<-list()
    
    for(i in 1:length(all_comb)){
      
      
      #one combination of sets (combination of 2, 3, 4...n sets)
      one_comb<-all_comb[[i]]
      
      for(j in 1:ncol(one_comb)){
        
        #one comparison for the combination (comparison of columns 1 & 2, or columns 1 & 2 & 3...)
        one_comp<-one_comb[,j]
        
        wanted_sets<-all_sets_frame[,one_comp]
        
        not_wanted_col<-(1:length(all_sets))[!(1:length(all_sets))%in%one_comp]
        
        not_wanted_sets<-all_sets_frame[,not_wanted_col]
        
        
        #intersection of the sets in the comparison
        for(k in 1:(ncol(wanted_sets)-1)){
          
          if(k==1){
          
            intersected<-intersect(unlist(wanted_sets[k][!is.na(wanted_sets[k])]),unlist(wanted_sets[k+1][!is.na(wanted_sets[k+1])]))
          
          }else{
            
            intersected<-intersect(intersected,unlist(wanted_sets[k+1][!is.na(wanted_sets[k+1])]))
            
          }
          
        }
        
        
        unique_intersected<-intersected[!intersected%in%unlist(not_wanted_sets[!is.na(not_wanted_sets)])]
        
        cat("intersection of sets : ",one_comp,"\n")
        
        cat("\t - name of sets : ",names(all_sets_frame[one_comp]),"\n")
        
        #cat("\t - common unique IDs : ",unique_intersected,"\n")
        
        cat("\t - length of common unique IDs : ",length(unique_intersected),"\n")
        
        all_intersections_counts[[length(all_intersections_counts)+1]]<-data.frame(c(paste(paste("cond",one_comp,sep="",collapse="_"),"_counts",sep=""),length(unique_intersected)))
        
        #if(paste(paste("cond",one_comp,sep="",collapse="_"),"_IDs",sep="")=="cond2_cond3_cond4_IDs"){stop("end")}
        
        all_intersections_IDs[[length(all_intersections_IDs)+1]]<-data.frame(c(paste(paste("cond",one_comp,sep="",collapse="_"),"_IDs",sep=""),paste(unique_intersected,collapse="-_-")))
        
        
      }
      
    }
    
    all_intersections_counts<-as.data.frame(all_intersections_counts)
    
    names(all_intersections_counts)<-as.character(unlist(all_intersections_counts[1,]))
    
    all_intersections_counts<-all_intersections_counts[2:nrow(all_intersections_counts),]
    
    all_intersections_IDs<-as.data.frame(all_intersections_IDs)
    
    names(all_intersections_IDs)<-as.character(unlist(all_intersections_IDs[1,]))
    
    
    #if we just have one column, we lose the type, if we use [li,col], instead of [col] directly...
    if(ncol(all_intersections_IDs)==1){
      
      previous_name<-names(all_intersections_IDs)
      
      all_intersections_IDs<-as.data.frame(all_intersections_IDs[2:nrow(all_intersections_IDs),])
      
      names(all_intersections_IDs)<-previous_name
      
    }else{
      
      all_intersections_IDs<-all_intersections_IDs[2:nrow(all_intersections_IDs),]
      
      
    }
    
    
    for(i in 1:ncol(all_intersections_IDs)){
      
      all_intersections_IDs[,i]<-gsub("-_-"," ",all_intersections_IDs[,i])
      
    }
    
    
    final_results<-cbind(pre_dataframe,all_intersections_counts,all_unique_IDs,all_intersections_IDs)
    
    
    return(final_results)

}




####################




# HOTAIR_LSD1

**HOTAIR lncRNA role in epithelial-mesenchymal transition**

1) Experiments

 ![](https://github.com/MorillonLab/HOTAIR_LSD1/blob/master/images_HOTAIR_LSD1/HOTAIR_LSD1_experiment.png)
 ![](https://github.com/MorillonLab/HOTAIR_LSD1/blob/master/images_HOTAIR_LSD1/analysis.png)

2) ChIP-seq analysis

   2-a) Description of the files
 
       IP BAM files :


               condition G (CTR) :
                           A685C9.bam
                           A685C10.bam

               condition H (HOT) :             
                           A685C11.bam
                           A685C12.bam  

               condition P (HOTdeltaP) :            
                           A685C13.bam
                           A685C14.bam

               condition L (HOTdeltaP) :           
                           A685C15.bam
                           A685C16.bam

       input BAM files :


               condition G (CTR) :
                           A684C1.bam
                           A685C2.bam

               condition H (HOT) :             
                           A684C3.bam
                           A685C4.bam  

               condition P (HOTdeltaP) :            
                           A684C6.bam
                           A684C6.bam

               condition L (HOTdeltaP) :           
                           A684C7.bam
                           A684C8.bam



   2-b) Peak calling script :

      https://github.com/MorillonLab/HOTAIR_LSD1/blob/master/sicer_script.sh

3) RNA-seq analysis

    3-a) DESeq2 script :
    
      https://github.com/MorillonLab/HOTAIR_LSD1/blob/master/DESeq2_script_example.R
  
    3-b) required dependancies :
    
             - Unix-like terminal (debian, ubuntu...) with grep installed.

             - R version >= 3.4.3, with :

               -a) R packages from CRAN

                 ggplot2
                 RColorBrewer
                 optparse
                 gplots
                 pROC
                 foreach
                 doParallel	
                 grid	
                 BiocParallel
                 data.table

                -b) R packages from Bioconductor

                  DESeq2
                  pheatmap
                  DOSE
                  clusterProfiler
                  
    3-c) Usage :
     
       In a Unix-like terminal, type :
     
         chmod 755 DESeq2_script_example.R
     
       To see all the available options, type :
     
          ./DESeq2_script_example.R
     
     
     

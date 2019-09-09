# HOTAIR_LSD1

**HOTAIR lncRNA role in epithelial-mesenchymal transition**

*1) Experiments*

 ![](https://github.com/MorillonLab/HOTAIR_LSD1/blob/master/images_HOTAIR_LSD1/HOTAIR_LSD1_experiment.png)
 ![](https://github.com/MorillonLab/HOTAIR_LSD1/blob/master/images_HOTAIR_LSD1/analysis.png)

 *2) ChIP-seq analysis*

     2-a) Description of the files
 
       IP BAM files :


               condition G (CONTROL) :
                           A685C9.bam
                           A685C10.bam

               condition H (HOTAIR) :             
                           A685C11.bam
                           A685C12.bam  

               condition P (HOTAIR-delta-PRC2) :            
                           A685C13.bam
                           A685C14.bam

               condition L (HOTAIR-delta-LSD1) :           
                           A685C15.bam
                           A685C16.bam

       input BAM files :


               condition G (CONTROL) :
                           A684C1.bam
                           A685C2.bam

               condition H (HOTAIR) :             
                           A684C3.bam
                           A685C4.bam  

               condition P (HOTAIR-delta-PRC2) :            
                           A684C6.bam
                           A684C6.bam

               condition L (HOTAIR-delta-LSD1) :           
                           A684C7.bam
                           A684C8.bam



    2-b) Peak calling script :

   https://github.com/MorillonLab/HOTAIR_LSD1/blob/master/sicer_script.sh

*3) RNA-seq analysis*

    3-a) DESeq2 script :
    
   https://github.com/MorillonLab/HOTAIR_LSD1/blob/master/DESeq2_script_example.R
  
    3-b) required dependancies :*
    
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
     
       - In a Unix-like terminal, type :
     
         chmod 755 DESeq2_script_example.R
     
       - To see all the available options, type :
     
          ./DESeq2_script_example.R
          
       - By default, the outputs are pairwise comparisons ; the minimum required parameters are "-f", "-a", and "-o".
     
       - To compare the conditions GL vs PH, the script has been run with the option :
       
           -t "two_firsts_vs_others_combined"
           
        - The design file supplied to the option "-f" looks like this (2 columns, tab-delimited)  :
        
          /home/marcgabriel/Documents/Marina_HOTAIR_lsd1/diff_expression/G1.counts.tsv	G
          /home/marcgabriel/Documents/Marina_HOTAIR_lsd1/diff_expression/G2.counts.tsv	G
          /home/marcgabriel/Documents/Marina_HOTAIR_lsd1/diff_expression/L2.counts.tsv	L
          /home/marcgabriel/Documents/Marina_HOTAIR_lsd1/diff_expression/L3.counts.tsv	L
          /home/marcgabriel/Documents/Marina_HOTAIR_lsd1/diff_expression/P2.counts.tsv	P
          /home/marcgabriel/Documents/Marina_HOTAIR_lsd1/diff_expression/P3.counts.tsv	P
          /home/marcgabriel/Documents/Marina_HOTAIR_lsd1/diff_expression/H2.counts.tsv	H
          /home/marcgabriel/Documents/Marina_HOTAIR_lsd1/diff_expression/H3.counts.tsv	H
      
      - The results of the comparison GL vs PH are in the directory "G_L_vs_P_H", with the following contents :
      
         - DESeq_output_*           : raw output table from DESeq2 (annotated with the gff file)
         - sig_diff_*               : table of all significant differentially expressed genes according to the filters (default : padj <= 0.05)
         - sig_diff_downregulated_* : table of significant downregulated genes (default : padj <= 0.05 & log2FC < 0)
         - sig_diff_upregulated_*   : table of significant upregulated genes (default : padj <= 0.05 & log2FC > 0)
         - MAplots & volcano-plots for controls

        
       
        
     

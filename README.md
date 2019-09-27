# HOTAIR lncRNA role in epithelial-mesenchymal transition

---

*1) Experiments & diagram of the analysis*

---

 ![](https://github.com/MorillonLab/HOTAIR_LSD1/blob/master/images_HOTAIR_LSD1/HOTAIR_LSD1_experiment.png)
 ![](https://github.com/MorillonLab/HOTAIR_LSD1/blob/master/images_HOTAIR_LSD1/analysis.png)
 
 ---

 *2) ChIP-seq analysis*
 
 ---

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

   -> all the parameters to modify are at the beginning of the script ("input data" part).
   
   -> the variable "blackRegions" can be initialized with this file : https://github.com/MorillonLab/HOTAIR_LSD1/blob/master/data/consensusBlacklist.bed
   
    2-c) Processing of peaks : 
(script : https://github.com/MorillonLab/HOTAIR_LSD1/blob/master/intersection_GPLH_Lsd1_peaks.sh)

       The peaks in gff format can be found at this link : https://www.dropbox.com/sh/m1p915bf8flh7vr/AABiI4GdI2S9H7_BKPZMeU07a?dl=0


       purpose :
       
         -> number of covered nucleotide per peak (images boxplot_distrib_covered_bases_per_peak_* & densityplot_distrib_covered_bases_per_peak_*)
         
         -> peak distance to TSS (images boxplot_distrib_dist_peak_to_TSS_*)
         
         -> peaks according to the genomic features : "exon", "promoter", "intron", "intergenic", "five_prime_utr"  "three_prime_utr" (images *_peak_distrib.png)
         
         -> annotated peaks (tables *_TSS_with_closest_peak.tsv)
         
         -> TSS with peaks around -5/+5 kb (table peaks_around_TSS_Minus5kb_Plus5kb.tsv)
         
         
       
---       
   
*3) RNA-seq analysis*

---

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
          
          -> these files are available at this adress (put the link)
      
      - The results of the comparison GL vs PH are in the directory "G_L_vs_P_H", with the following contents :
      
         - DESeq_output_*           : raw output table from DESeq2 (annotated with the gff file).
         - sig_diff_*               : table of all significant differentially expressed genes according to the filters (default : padj <= 0.05).
         - sig_diff_downregulated_* : table of significant downregulated genes (default : padj <= 0.05 & log2FC < 0).
         - sig_diff_upregulated_*   : table of significant upregulated genes (default : padj <= 0.05 & log2FC > 0).
         - MAplots & volcano-plots for controls.


---       
   
*4) intersection RNA-seq ChIP-seq*

---    
       
        
     

# HOTAIR lncRNA role in epithelial-mesenchymal transition

---

*1) Experiments & diagram of the analysis*

---

 ![](https://github.com/MorillonLab/HOTAIR_LSD1/blob/master/images_HOTAIR_LSD1/HOTAIR_LSD1_experiment.png)
 ![](https://github.com/MorillonLab/HOTAIR_LSD1/blob/master/images_HOTAIR_LSD1/analysis.png)
 
 ---
 
 *2) RNA-seq analysis*

---

    2-a) DESeq2 script :
    
   https://github.com/MorillonLab/HOTAIR_LSD1/blob/master/DESeq2_script_example.R
  
  []()
  
    2-b) required dependancies :
    
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
                  
[]()
    
    2-c) Usage :
     
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
          
          -> these files are available at this adress :
 http://xfer.curie.fr/get/cHZjhx2qJg2/HOTAIR_LSD1_RNAseq_gene_counts.zip
 
      - The results of the comparison GL vs PH are in the directory "G_L_vs_P_H", with the following contents :
      
         - DESeq_output_*           : raw output table from DESeq2 (annotated with the gff file).
         - sig_diff_*               : table of all significant differentially expressed genes according to the filters (default : padj <= 0.05).
         - sig_diff_downregulated_* : table of significant downregulated genes (default : padj <= 0.05 & log2FC < 0).
         - sig_diff_upregulated_*   : table of significant upregulated genes (default : padj <= 0.05 & log2FC > 0).
         - MAplots & volcano-plots for controls.
 
 
      We define as :
          - HMS (High migration set) : genes upregulated in H & P conditions (padj < 0.05, FoldChange > 1.5).
          - LMS (Low migration set)  : genes upregulated in G & L conditions (= downregulated in H & P conditions ; padj < 0.05, FoldChange <  0.69).
    
 
 ---

 *3) ChIP-seq analysis*
 
 ---

     3-a) Description of the files
 
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

[]()

    3-b) Peak calling script :

   https://github.com/MorillonLab/HOTAIR_LSD1/blob/master/sicer_script.sh

   -> all the parameters to modify are at the beginning of the script ("input data" part).
   
   -> the variable "blackRegions" can be initialized with this file : https://github.com/MorillonLab/HOTAIR_LSD1/blob/master/data/consensusBlacklist.bed
   
   -> The peaks in gff format can be found at this link : http://xfer.curie.fr/get/q1GyGKkrEIQ/HOTAIR_LSD1_merged_filteredBR_peaks_gff.zip
    
    3-c) Processing of peaks around TSS (-5/+5kb): 
(script : https://github.com/MorillonLab/HOTAIR_LSD1/blob/master/intersection_GPLH_Lsd1_peaks.sh)
    


The gencode annotation for the variable "annotation" can be found here (uncompress it) : 

ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/GRCh37_mapping/gencode.v26lift37.annotation.gff3.gz

       Results of the script :
       
         -> number of covered nucleotide per peak (images boxplot_distrib_covered_bases_per_peak_* & densityplot_distrib_covered_bases_per_peak_*)
         
         -> peak distance to TSS (images boxplot_distrib_dist_peak_to_TSS_*)
         
         -> genomic distribution of the peaks : "exon", "promoter", "intron", "intergenic", "five_prime_utr"  "three_prime_utr" (images *_peak_distrib.png)
         
         -> annotated peaks (tables *_TSS_with_closest_peak.tsv)
         
         -> TSS with peaks around -5/+5 kb (table peaks_around_TSS_Minus5kb_Plus5kb.tsv)
         
[]()         
          
    
    3-d) Processing of peaks inside genes :    
(script : https://github.com/MorillonLab/HOTAIR_LSD1/blob/master/intersection_GPLH_Lsd1_peaks_inGenes.sh)

      purpose : same as above, but for peaks inside genes
      
[]()

 ---

 *4) Intersection of peaks (-5/+5 around TSS and inside of genes) with RNA-seq (upregulated genes in High & Low migration phenotype)*
 
 ---
 
 -> Remarks : As the final peaks used until there are the ones common between replicates, to make the intersection of the peaks with the differentially expressed genes (HMS & LMS), and avoiding bias, we need to use the peaks of each replicate.
 
 -> The peaks not merged can be found here (go to section 3-a) to see to which condition each file with the prefix A685C* belongs) : http://xfer.curie.fr/get/9JHuo1YBC2w/HOTAIR_LSD1_not_merged_filteredBR_peaks_bed.zip

[]()

    4-a) Run in 1st this script for the peaks around the TSS : 
    
script for TSS (replace the peaks in bed format in "input data" with the ones in the link above) :
https://github.com/MorillonLab/HOTAIR_LSD1/blob/master/intersection_GPLH_Lsd1_peaks_rep_not_merged.sh

    warning : use a new output directory name for the results

[]()

    4-b) Run in 2nd this script for the peaks inside the genes :

script for genebody (replace the peaks in bed format in "input data" with the ones in the link above) :
https://github.com/MorillonLab/HOTAIR_LSD1/blob/master/intersection_GPLH_Lsd1_peaks_rep_not_merged_genebody.sh
     
     warning : use a new output directory name for the results
 
[]()

    4-c) Run in 3rd this script for to combine the results of 4-a) & 4-b):

script to combine results of 4-a) & 4-b) https://github.com/MorillonLab/HOTAIR_LSD1/blob/master/intersection_GPLH_Lsd1_peaks_rep_not_merged_genebody.sh
     
     Usage :

       ./venn_lsd1_diff_TSS_and_gene_body_combined.R home1 home2 home3 official_annotation

       
       home1 : output directory for this script

       home2 : output directory of the script "intersection_GPLH_Lsd1_peaks_rep_not_merged_genebody.sh"

       home3 : output directory of the script "intersection_GPLH_Lsd1_peaks_rep_not_merged.sh"
       
    The "official_annotation" parameter is the file used in the section 3-c) : ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/GRCh37_mapping/gencode.v26lift37.annotation.gff3.gz (uncompress it)

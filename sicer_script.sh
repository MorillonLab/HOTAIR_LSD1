#!/bin/bash

###### purpose of the script #######

#Run SICER on ChIP BAM files (IP & input)
#create for each condition the common peaks between two replicates -> files *merged_filteredBR.peaks.bed
#create the mean  coverage in bigWig format of the filtered islands (enriched regions) -> files *_islands_normalized_mean.bw


####################################

##### input data #######

#path to bedtools
bedtools="bedtools"

#path to SICER
SICER="/home/marcgabriel/Downloads/SICER_V1.1/SICER/SICER.sh"

#location of ChIP bam files
inputDir="/media/marcgabriel/Transcend/Lsd1_chipseq_mapping_files/"

#output dir
outputDir="/media/marcgabriel/Transcend/LSD1_metagenes/GHPL_sicer_outputs/"

#SICER parameters
specie="hg19"
redundancy_threshold=1
window=200
genomeFraction=0.74
fragmentSize=100
Evalue=100
FDR=0.05
gap=1000

#regions to remove from sicer results
blackRegions="/media/marcgabriel/Transcend/LSD1_metagenes/consensusBlacklist.bed"

export wigtobigwig="/home/marcgabriel/Desktop/scripts/wigToBigWig"

#bigwigCompare from deeptools
export bigwigCompare="bigwigCompare"

###########################



####################### design of the ChIP experiment ########

#IP of the first table is related to the input of the second table (order is conserved)

#IP files
rep1_all_IP_bed_files=("A685C9.bam"
                       "A685C11.bam"
                       "A685C13.bam"
                       "A685C15.bam"
                       "A685C10.bam"
                       "A685C12.bam"
                       "A685C14.bam"
                       "A685C16.bam")


#input files                    
rep1_all_input_bed_files=("A684C1.bam"
                          "A684C3.bam"
                          "A684C5.bam"
                          "A684C7.bam"
                          "A684C2.bam"
                          "A684C4.bam"
                          "A684C6.bam"
                          "A684C8.bam")
                          
                          
###############################################################
             



########## run sicer for each pair IP/input ##################

for one_file in $(seq 0 $((${#rep1_all_IP_bed_files[*]}-1)));do

   if [ ! -f "${outputDir}${rep1_all_IP_bed_files[$one_file]}-W200-G1000-FDR0.05-islandfiltered.bed" ];then
   
      echo -e "sicer files missing for ${rep1_all_IP_bed_files[$one_file]}.unique.bam, we're going to re-build them...\n"

      #take only unique mapped reads for IP
      samtools view -h ${inputDir}${rep1_all_IP_bed_files[$one_file]}|grep -v "XS:i"|samtools view -Sbh - >${outputDir}${rep1_all_IP_bed_files[$one_file]}.unique.bam
    
      #take only unique mapped reads for input
      samtools view -h ${inputDir}${rep1_all_input_bed_files[$one_file]}|grep -v "XS:i"|samtools view -Sbh - >${outputDir}${rep1_all_input_bed_files[$one_file]}.unique.bam
    
      #convert both files in bed format (need by SICER)
      bedtools bamtobed -i ${outputDir}${rep1_all_IP_bed_files[$one_file]}.unique.bam >${outputDir}${rep1_all_IP_bed_files[$one_file]}.bed
      bedtools bamtobed -i ${outputDir}${rep1_all_input_bed_files[$one_file]}.unique.bam >${outputDir}${rep1_all_input_bed_files[$one_file]}.bed
    
      #sicer command
      $SICER ${outputDir} ${rep1_all_IP_bed_files[$one_file]}.bed ${rep1_all_input_bed_files[$one_file]}.bed $outputDir $specie $redundancy_threshold $window $fragmentSize $genomeFraction $gap $FDR || { echo "SICER  failure on ${inputDir}${rep1_all_IP_bed_files[$one_file]} & ${inputDir}${rep1_all_input_bed_files[$one_file]} !" 1>&2; exit; }
   
   
   else
   
     echo -e "all files for ${rep1_all_IP_bed_files[$one_file]}.unique.bam are there :)\n"
   
   
   fi
   
   $bedtools intersect -v -a ${outputDir}${rep1_all_IP_bed_files[$one_file]}-W200-G1000-islands-summary-FDR0.05 -b $blackRegions > ${outputDir}${rep1_all_IP_bed_files[$one_file]}.unique-W200-G1000-islands-summary-FDR0.05_filteredBR.bed
   
done

##################################################################


################ take only peaks common in both replicates for a given condition #########

#function to merge the common peaks in the replicates for each condition ######
function mergePeaks {
	
	var1=$1
	var2=$2
	#give peak id to define file source
	awk '{print $0"\t","peakA685C'${var1}'-"NR}' ${outputDir}A685C${var1}.bam.unique-W200-G1000-islands-summary-FDR0.05_filteredBR.bed | cut -f 1,2,3,9 > ${outputDir}A685C${var1}.unique-W200-G1000-islands-summary-FDR0.05_filteredBR_id.bed
	awk '{print $0"\t","peakA685C'${var2}'-"NR}' ${outputDir}A685C${var2}.bam.unique-W200-G1000-islands-summary-FDR0.05_filteredBR.bed | cut -f 1,2,3,9 > ${outputDir}A685C${var2}.unique-W200-G1000-islands-summary-FDR0.05_filteredBR_id.bed
	#sort them by chr & start and then merge with mergeBed
	cat ${outputDir}A685C${var1}.unique-W200-G1000-islands-summary-FDR0.05_filteredBR_id.bed ${outputDir}A685C${var2}.unique-W200-G1000-islands-summary-FDR0.05_filteredBR_id.bed | sort -k1,1 -k2,2n | $bedtools merge -i stdin -o collapse -c 4 > ${outputDir}A685C${var1}.${var2}.unique-W200-G1000-islands-summary-FDR0.05_filteredBR_id.merged.bed
	#retrieve only overlaps
	grep "," ${outputDir}A685C${var1}.${var2}.unique-W200-G1000-islands-summary-FDR0.05_filteredBR_id.merged.bed > ${outputDir}A685C${var1}.${var2}.merged_filteredBR.peaks.bed
	rm ${outputDir}A685C${var1}.${var2}.unique-W200-G1000-islands-summary-FDR0.05_filteredBR_id.merged.bed
	rm ${outputDir}A685C${var1}.unique-W200-G1000-islands-summary-FDR0.05_filteredBR_id.bed
	rm ${outputDir}A685C${var2}.unique-W200-G1000-islands-summary-FDR0.05_filteredBR_id.bed
	
}

#merge replicates by replicates (for a given condition)
#cond G
mergePeaks 9 10
#cond H
mergePeaks 11 12
#cond P
mergePeaks 13 14
#cond L
mergePeaks 15 16

####################################################################################################


##### make the bigwig files ####
# for each condition, make the mean of the normalized islands
################################


cond_list=("G" "H" "P" "L")

            rep1_files=("A685C9.bam"
                       "A685C11.bam"
                       "A685C13.bam"
                       "A685C15.bam")
                       
                       

            rep2_files=("A685C10.bam"
                        "A685C12.bam"
                        "A685C14.bam"
                        "A685C16.bam"  )
                        
subscripts_dir2="${outputDir}subscripts2/"

if [ -d $subscripts_dir2 ];then rm -rf $subscripts_dir2 ;fi

mkdir $subscripts_dir2  
                        
for one_cond in $(seq 0 $((${#cond_list[*]}-1)));do

  echo -e "#!/bin/bash\n" >${subscripts_dir2}${cond_list[$one_cond]}_subscript.sh

  echo -e "processing cond : ${cond_list[$one_cond]}...\n"

  samtools view -H ${inputDir}${rep1_files[$one_cond]}|grep -v -E "user|VN|ID" | awk 'OFS="\t"{print $2,$3}' | sed 's/SN\://g'|sed 's/LN\://g' >${outputDir}${cond_list[$one_cond]}1_chr_size.txt
  samtools view -H ${inputDir}${rep2_files[$one_cond]}|grep -v -E "user|VN|ID" | awk 'OFS="\t"{print $2,$3}' | sed 's/SN\://g'|sed 's/LN\://g' >${outputDir}${cond_list[$one_cond]}2_chr_size.txt
  
  
  if [ ! -f ${outputDir}${cond_list[$one_cond]}1.bw ];then

    echo -e "$wigtobigwig ${outputDir}${rep1_files[$one_cond]}-W200-G1000-FDR0.05-islandfiltered-normalized.wig ${outputDir}${cond_list[$one_cond]}1_chr_size.txt ${outputDir}${cond_list[$one_cond]}1.bw\n\n" >>${subscripts_dir2}${cond_list[$one_cond]}_subscript.sh
    
  fi
  
  if [ ! -f ${outputDir}${cond_list[$one_cond]}2.bw ];then
  
    echo -e "$wigtobigwig ${outputDir}${rep2_files[$one_cond]}-W200-G1000-FDR0.05-islandfiltered-normalized.wig ${outputDir}${cond_list[$one_cond]}2_chr_size.txt ${outputDir}${cond_list[$one_cond]}2.bw\n" >>${subscripts_dir2}${cond_list[$one_cond]}_subscript.sh
  
  fi
  
  
  if [ ! -f ${outputDir}${cond_list[$one_cond]}_islands_normalized_mean.bw ];then
  
    echo -e "$bigwigCompare --bigwig1 ${outputDir}${cond_list[$one_cond]}1.bw --bigwig2 ${outputDir}${cond_list[$one_cond]}2.bw --operation \"mean\" --numberOfProcessors 4 -o ${outputDir}${cond_list[$one_cond]}_islands_normalized_mean.bw --binSize 1 --deepBlueTempDir ${outputDir}\n" >>${subscripts_dir2}${cond_list[$one_cond]}_subscript.sh
    
  fi
  
  chmod 755 ${subscripts_dir2}${cond_list[$one_cond]}_subscript.sh

done  

find ${subscripts_dir2} -name "*_subscript.sh" | xargs -n 1 -P 2 bash || { echo "executing bash subscripts2 failure !" 1>&2; exit; }


##########################################
              





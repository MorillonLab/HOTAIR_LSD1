#!/bin/bash


##### input data #######

export samtools="samtools"

#use bedtools 2.27, with 2.26, groupby is buggy !!
export bedtools="/home/marcgabriel/Desktop/bedtools2/bin/bedtools"

#bedops program
export bedops="/home/marcgabriel/Desktop/bedops-2.4.30/bin/bedops"

#gff
export annotation="/home/marcgabriel/Documents/gencode26lift37/gencode.v26lift37.annotation.gff3"

#create intergenic features
export getIntergenics=$(dirname "$0")/getIntergenics.R

#create intronic features
export getIntronsByTranscripts=$(dirname "$0")/getIntronsByTranscripts.R

#classify the peaks
export classifyPeaks=$(dirname "$0")/classifyPeaks.R

#output dir
output_dir="/media/marcgabriel/eda138bc-95c2-4778-99df-8915815cb86e/Marina_lsd1_test/"

rscript1=$(dirname "$0")/venn_diagrams_with_numerics.R

rscript2=$(dirname "$0")/venn_chipseq_GHPL.R

getBoxplots=$(dirname "$0")/chipeseq_GHPL_boxplots_dist_to_tss.R

getCoveredNuc=$(dirname "$0")/getCoveredNuc.R

#annotation in gff format
ref_annot="/home/marcgabriel/Documents/gencode26lift37/gencode.v26lift37.annotation.gff3"

#enriched regions in gff format
G="/media/marcgabriel/Transcend/LSD1_metagenes/original_peaks_replicates_merged/A685C9.10.merged_filteredBR.peaks.gff"
#G="/media/marcgabriel/Maxtor/Lsd1/Chipseq/peak_calling/sicer/run_1/gap1000/filtered_BR/A685C9.unique-W200-G1000-islands-summary-FDR0.05_filteredBR.bed /media/marcgabriel/Maxtor/Lsd1/Chipseq/peak_calling/sicer/run_1/gap1000/filtered_BR/A685C10.unique-W200-G1000-islands-summary-FDR0.05_filteredBR.bed"

L="/media/marcgabriel/Transcend/LSD1_metagenes/original_peaks_replicates_merged/A685C15.16.merged_filteredBR.peaks.gff"
#L="/media/marcgabriel/Maxtor/Lsd1/Chipseq/peak_calling/sicer/run_1/gap1000/filtered_BR/A685C15.unique-W200-G1000-islands-summary-FDR0.05_filteredBR.bed /media/marcgabriel/Maxtor/Lsd1/Chipseq/peak_calling/sicer/run_1/gap1000/filtered_BR/A685C16.unique-W200-G1000-islands-summary-FDR0.05_filteredBR.bed"

H="/media/marcgabriel/Transcend/LSD1_metagenes/original_peaks_replicates_merged/A685C11.12.merged_filteredBR.peaks.gff"
#H="/media/marcgabriel/Maxtor/Lsd1/Chipseq/peak_calling/sicer/run_1/gap1000/filtered_BR/A685C11.unique-W200-G1000-islands-summary-FDR0.05_filteredBR.bed /media/marcgabriel/Maxtor/Lsd1/Chipseq/peak_calling/sicer/run_1/gap1000/filtered_BR/A685C12.unique-W200-G1000-islands-summary-FDR0.05_filteredBR.bed"

P="/media/marcgabriel/Transcend/LSD1_metagenes/original_peaks_replicates_merged/A685C13.14.merged_filteredBR.peaks.gff"
#P="/media/marcgabriel/Maxtor/Lsd1/Chipseq/peak_calling/sicer/run_1/gap1000/filtered_BR/A685C13.unique-W200-G1000-islands-summary-FDR0.05_filteredBR.bed /media/marcgabriel/Maxtor/Lsd1/Chipseq/peak_calling/sicer/run_1/gap1000/filtered_BR/A685C14.unique-W200-G1000-islands-summary-FDR0.05_filteredBR.bed"


#####################################


output_dir="${output_dir}/"
output_dir=$(echo $output_dir |sed 's/\/\//\//g')

if [ ! -d "$output_dir" ];then


  mkdir $output_dir

fi


grep -v "^#" $annotation |grep -P "\ttranscript\t"|awk 'OFS="\t"{if($7=="+"){$5=$4}else{$4=$5};print}' | sort -k1,1 -k4,4n >${output_dir}all_TSS.gff


original_G=$G
original_L=$L
original_H=$H
original_P=$P

################

####### compute number of covered nucleotide per peak

awk -F'\t' 'OFS="\t"{print "G",($5-$4)+1}' $original_G >${output_dir}covered_nuc.tsv
awk  -F'\t' 'OFS="\t"{print "H",($5-$4)+1}' $original_H >>${output_dir}covered_nuc.tsv
awk  -F'\t' 'OFS="\t"{print "P",($5-$4)+1}' $original_P >>${output_dir}covered_nuc.tsv
awk  -F'\t' 'OFS="\t"{print "L",($5-$4)+1}' $original_L >>${output_dir}covered_nuc.tsv

$getCoveredNuc ${output_dir}


#######################################################


######### processing of bed loc files (dev) ###############

#compute midpeak : |awk 'BEGIN {OFS="\t";printf($1"\t");printf($2"\t");printf($3"\t");printf("%.0f\t", ($4+$5)/2);printf("%.0f\t", ($4+$5)/2);print $6,$7,$8,$9 }'

#concatenate the peaks (from all rep), then merge them
#sort -k1,1 -k2,2n
#cat $G |sort -k1,1 -k2,2n|awk 'OFS="\t"{print $0,"+"}' >${output_dir}G_peaks_all_rep.tmp && $bedtools merge -d 1 -i ${output_dir}G_peaks_all_rep.tmp | awk 'OFS="\t"{print $1,"annotation","peaks",$2+1,$3,".","+",".","G1_peak_number_"NR}'|awk '{OFS="\t";printf($1"\t");printf($2"\t");printf($3"\t");printf("%.0f\t", ($4+$5)/2);printf("%.0f\t", ($4+$5)/2);print $6,$7,$8,$9 }' >${output_dir}G_peaks_all_rep.gff3 && rm ${output_dir}G_peaks_all_rep.tmp
#G="${output_dir}G_peaks_all_rep.gff3"

#cat $L |sort -k1,1 -k2,2n |awk 'OFS="\t"{print $0,"+"}' >${output_dir}L_peaks_all_rep.tmp && $bedtools merge -d 1 -i ${output_dir}L_peaks_all_rep.tmp | awk 'OFS="\t"{print $1,"annotation","peaks",$2+1,$3,".","+",".","L1_peak_number_"NR}'|awk '{OFS="\t";printf($1"\t");printf($2"\t");printf($3"\t");printf("%.0f\t", ($4+$5)/2);printf("%.0f\t", ($4+$5)/2);print $6,$7,$8,$9 }' >${output_dir}L_peaks_all_rep.gff3 && rm ${output_dir}L_peaks_all_rep.tmp
#L="${output_dir}L_peaks_all_rep.gff3"

#cat $H |sort -k1,1 -k2,2n |awk 'OFS="\t"{print $0,"+"}' >${output_dir}H_peaks_all_rep.tmp && $bedtools merge -d 1 -i ${output_dir}H_peaks_all_rep.tmp | awk 'OFS="\t"{print $1,"annotation","peaks",$2+1,$3,".","+",".","H1_peak_number_"NR}'|awk '{OFS="\t";printf($1"\t");printf($2"\t");printf($3"\t");printf("%.0f\t", ($4+$5)/2);printf("%.0f\t", ($4+$5)/2);print $6,$7,$8,$9 }' >${output_dir}H_peaks_all_rep.gff3 && rm ${output_dir}H_peaks_all_rep.tmp
#H="${output_dir}H_peaks_all_rep.gff3"

#cat $P |sort -k1,1 -k2,2n |awk 'OFS="\t"{print $0,"+"}' >${output_dir}P_peaks_all_rep.tmp && $bedtools merge -d 1 -i ${output_dir}P_peaks_all_rep.tmp | awk 'OFS="\t"{print $1,"annotation","peaks",$2+1,$3,".","+",".","P1_peak_number_"NR}'|awk '{OFS="\t";printf($1"\t");printf($2"\t");printf($3"\t");printf("%.0f\t", ($4+$5)/2);printf("%.0f\t", ($4+$5)/2);print $6,$7,$8,$9 }' >${output_dir}P_peaks_all_rep.gff3 && rm ${output_dir}P_peaks_all_rep.tmp

cat $G |sort -k1,1 -k4,4n|awk 'OFS="\t"{print $0,"+"}' |awk 'OFS="\t"{printf($1"\t");printf("annotation""\t");printf("peaks""\t");printf("%.0f\t", $4);printf("%.0f\t", $5);print ".","+",".","G1_peak_number_"NR}' >${output_dir}G_peaks_all_rep.gff3
G="${output_dir}G_peaks_all_rep.gff3"

cat $H |sort -k1,1 -k4,4n|awk 'OFS="\t"{print $0,"+"}' |awk 'OFS="\t"{printf($1"\t");printf("annotation""\t");printf("peaks""\t");printf("%.0f\t", $4);printf("%.0f\t", $5);print ".","+",".","H1_peak_number_"NR}' >${output_dir}H_peaks_all_rep.gff3
H="${output_dir}H_peaks_all_rep.gff3"

cat $L |sort -k1,1 -k4,4n|awk 'OFS="\t"{print $0,"+"}' |awk 'OFS="\t"{printf($1"\t");printf("annotation""\t");printf("peaks""\t");printf("%.0f\t", $4);printf("%.0f\t", $5);print ".","+",".","L1_peak_number_"NR}' >${output_dir}L_peaks_all_rep.gff3
L="${output_dir}L_peaks_all_rep.gff3"

cat $P |sort -k1,1 -k4,4n|awk 'OFS="\t"{print $0,"+"}' |awk 'OFS="\t"{printf($1"\t");printf("annotation""\t");printf("peaks""\t");printf("%.0f\t", $4);printf("%.0f\t", $5);print ".","+",".","P1_peak_number_"NR}' >${output_dir}P_peaks_all_rep.gff3
P="${output_dir}P_peaks_all_rep.gff3"



echo -e "new G is : $G \n new L is $L \n new H is $H\n new P is $P"

#dev
##
cat $original_G |sort -k1,1 -k4,4n|awk -F'\t' 'OFS="\t"{printf($1"\t");printf("annotation""\t");printf("peaks""\t");printf("%.0f\t", $4);printf("%.0f\t", $5);print ".","+",".",$9,$4,$5}'|awk -F'\t' 'OFS="\t"{$9=gensub(" ","+","G",$9);$9=$9"-_-"$10"-_-"$11;$10="";$11="";print}'|cut -f1-9 >${output_dir}G_peaks_all_rep_bis.gff3
G_bis="${output_dir}G_peaks_all_rep_bis.gff3"


cat $original_H |sort -k1,1 -k4,4n|awk -F'\t' 'OFS="\t"{printf($1"\t");printf("annotation""\t");printf("peaks""\t");printf("%.0f\t", $4);printf("%.0f\t", $5);print ".","+",".",$9,$4,$5}'|awk -F'\t' 'OFS="\t"{$9=gensub(" ","+","G",$9);$9=$9"-_-"$10"-_-"$11;$10="";$11="";print}' |cut -f1-9 >${output_dir}H_peaks_all_rep_bis.gff3
H_bis="${output_dir}H_peaks_all_rep_bis.gff3"

cat $original_L |sort -k1,1 -k4,4n|awk -F'\t' 'OFS="\t"{printf($1"\t");printf("annotation""\t");printf("peaks""\t");printf("%.0f\t", $4);printf("%.0f\t", $5);print ".","+",".",$9,$4,$5}' |awk -F'\t' 'OFS="\t"{$9=gensub(" ","+","G",$9);$9=$9"-_-"$10"-_-"$11;$10="";$11="";print}'|cut -f1-9 >${output_dir}L_peaks_all_rep_bis.gff3
L_bis="${output_dir}L_peaks_all_rep_bis.gff3"

cat $original_P |sort -k1,1 -k4,4n|awk -F'\t' 'OFS="\t"{printf($1"\t");printf("annotation""\t");printf("peaks""\t");printf("%.0f\t", $4);printf("%.0f\t", $5);print ".","+",".",$9,$4,$5}' |awk -F'\t' 'OFS="\t"{$9=gensub(" ","+","G",$9);$9=$9"-_-"$10"-_-"$11;$10="";$11="";print}'|cut -f1-9 >${output_dir}P_peaks_all_rep_bis.gff3
P_bis="${output_dir}P_peaks_all_rep_bis.gff3"


#end of dev



###################################### end of dev ##########


#################


#if [ ! -d $output_dir ];then mkdir $output_dir ;fi

#sed 's/ /_/g' $G >${output_dir}file1.gff

#sed 's/ /_/g' $L >${output_dir}file2.gff



#summary=${output_dir}summary.txt

#annot_list=(${output_dir}file1.gff ${output_dir}file2.gff)
            
#prefix=("G" "L")

#echo -e "- annot list : \n\n\t$(echo -e ${annot_list[*]}|sed 's/ /\n\t/g')\n"

#echo -e "- prefixes : \n\n\t$(echo -e ${prefix[*]}|sed 's/ /\n\t/g')\n\n==============\n\n"


#echo -e "cond1\tcond2\tcond1_counts\tcond2_counts\tminimal_overlap\tmax_overlap_cond1\tmax_overlap_cond2\tmax_overlap_cond1_percent\tmax_overlap_cond2_percent" >${output_dir}all_intersections.tsv

##overlap_frac="-f 0.50 -F 0.50"

#overlap_frac=""

##all possible combinations in the annot list (for 3 annotations for example, we have 1 vs 2, 1vs 3, & 2 vs 3)
##in bash, incrementation in lists start from 0, contrary to R

#echo -e "fraction of overlap between annotation 1 & 2 : $overlap_frac\n" >${output_dir}summary.txt

#for i in $(seq 0 $((${#annot_list[*]}-1)));do

  #new_list=($(seq $((i+1)) $((${#annot_list[*]}-1))))


  #for j in ${new_list[*]};do
  
         #type1=""
         
         #type2=""
  
         #if [[ $(echo ${annot_list[$i]} |grep "\.bed" |wc -l) -gt 0 ]];then
         
           #type1="bed"
         
         
         #fi
         
         #if [[ $(echo ${annot_list[$i]}|grep "\.gff" |wc -l) -gt 0 ]];then
         
             #type1="gff"
         
         
         #fi
         
         #if [[ $(echo ${annot_list[j]} |grep "\.bed" |wc -l) -gt 0 ]];then
         
           #type2="bed"
         
         
         #fi
         
         #if [[ $(echo ${annot_list[j]} | grep "\.gff" |wc -l) -gt 0 ]];then
         
             #type2="gff"
         
         
         #fi
         
  
         
         
  
	  #cond1_counts=$(wc -l ${annot_list[$i]} |awk '{print $1}')
	  
	  #cond2_counts=$(wc -l ${annot_list[j]} |awk '{print $1}')
	  
	  #echo -e "- annotations to intersect : ${annot_list[$i]} & ${annot_list[j]} (= ${prefix[$i]} & ${prefix[j]})\n" >>$summary
	  
	  #echo -e "type1 : $type1; type2 : $type2\n" >>$summary
	  
	  #echo -e "\t- number of feat in ${prefix[$i]} : $cond1_counts\n" >>$summary
	  
	  #echo -e "\t- number of feat in ${prefix[j]} : $cond2_counts\n" >>$summary

	  #$bedtools intersect $overlap_frac -e -wb -nobuf -nonamecheck -a ${annot_list[$i]} -b ${annot_list[j]} >${output_dir}intersected.tmp
	 
	  #if [ "$type1" == "gff" ];then
	  
	    #max_overlap_cond1=$(cut -f9 ${output_dir}intersected.tmp |grep -v "^$"|sort -u |wc -l )
	    
	      #if [ "$type2" == "gff" ];then
	    
	         #max_overlap_cond2=$(cut -f18 ${output_dir}intersected.tmp |grep -v "^$"|sort -u|wc -l )
	    
	      #elif [ "$type2" == "bed" ];then
	      
	         #max_overlap_cond2=$(cut -f16 ${output_dir}intersected.tmp |grep -v "^$"|sort -u|wc -l )
	      
	      #fi
	  
	  
	  #fi
	  
	  #if [ "$type1" == "bed" ];then
	  
	  
	      #max_overlap_cond1=$(cut -f4 ${output_dir}intersected.tmp |grep -v "^$"|sort -u |wc -l )
	  
	      #if [ "$type2" == "gff" ];then
	    
	         #max_overlap_cond2=$(cut -f21 ${output_dir}intersected.tmp |grep -v "^$"|sort -u|wc -l )
	    
	      #elif [ "$type2" == "bed" ];then
	      
	         #max_overlap_cond2=$(cut -f16 ${output_dir}intersected.tmp |grep -v "^$"|sort -u|wc -l )
	      
	      #fi
	  
	  
	  #fi
	  
	  #if [[ $max_overlap_cond1 -gt $max_overlap_cond2 ]];then
	  
	      #shared_counts=$max_overlap_cond2
	      
	  #else
	  
	      #shared_counts=$max_overlap_cond1
	  
	  
	  #fi
	  
	  ##shared_counts=$(wc -l ${output_dir}intersected.tmp|awk '{print $1}')
	  
	  #percent_counts=($(awk -v cond1_counts=$cond1_counts -v cond2_counts=$cond2_counts -v shared_counts=$shared_counts 'BEGIN{OFS="\t";print (shared_counts/cond1_counts)*100,(shared_counts/cond2_counts)*100}'))
	  
	  #percent_counts_diff=($(awk -v cond1_counts=$cond1_counts -v cond2_counts=$cond2_counts -v max_overlap_cond1=$max_overlap_cond1 -v max_overlap_cond2=$max_overlap_cond2 'BEGIN{OFS="\t";print (max_overlap_cond1/cond1_counts)*100,(max_overlap_cond2/cond2_counts)*100}'))
	  
	  #percent_in_cond1=${percent_counts[0]}
	  
	  #percent_in_cond2=${percent_counts[1]}
	  
	  #percent_in_cond1_diff=${percent_counts_diff[0]}
	  
	  #percent_in_cond2_diff=${percent_counts_diff[1]}
	  
	  #echo -e "\t\t- number of common features : $shared_counts (= ${percent_in_cond1}% of ${prefix[$i]} & ${percent_in_cond2}% of ${prefix[j]})\n" >>$summary
	  
	  #echo -e "\t\t- number of features for each (diff) : $max_overlap_cond1 = ${percent_in_cond1_diff}% of ${prefix[$i]} ; $max_overlap_cond2 = ${percent_in_cond2_diff}% of ${prefix[j]})\n" >>$summary
	  
	  #echo -e "\n************\n" >>$summary
	  
	  ##if [ "${prefix[$i]}" ==  "class1w" ] && [ "${prefix[$j]}" ==  "lncipedia" ];then exit;fi
	  
	  ##rm ${output_dir}intersected.tmp
	  
	  #echo -e "${prefix[$i]}\t${prefix[j]}\t${cond1_counts}\t${cond2_counts}\t$shared_counts\t$max_overlap_cond1\t$max_overlap_cond2\t${percent_in_cond1_diff}\t${percent_in_cond2_diff}" >${output_dir}${prefix[$i]}_${prefix[j]}_intersections.tsv
	  
	  #cat ${output_dir}${prefix[$i]}_${prefix[j]}_intersections.tsv >>${output_dir}all_intersections.tsv
	  
  #done

#done



#if [ "$overlap_frac" == "" ];then

  #prefix="no_cutoff"

#else

  #prefix=$(echo $overlap_frac |sed 's/ /\n/g' |awk 'NR==2{print "cutoff_"$0*100"_percent_overlap"}')


#fi


#$rscript1 $output_dir ${output_dir}all_intersections.tsv $prefix "ChIP-seq Peaks G vs L"


######### closest peaks to TSS ####

echo -e "condition\tnumber" >${output_dir}peaks_around_TSS_Minus5kb_Plus5kb.tsv


##look for closest G peak to TSS -5kb/2kb
${bedtools} closest -a ${output_dir}all_TSS.gff -b ${G} -t first -D a -sorted -nonamecheck >${output_dir}G_closest_feat.bed || { echo "closestbed  failure 1 !" 1>&2; exit; }

#taking the transcript with the closest TSS to peak (priority to dist 0, then the closest in absolute number ; priority to negative dist if the same absolute dist is present)
awk 'OFS="\t"{split($9,a,";");b="";for(i=1;i<=length(a);i++){if(a[i]~/gene_name/){b=a[i]}};print $10,$1,$4,$7,b,$NF}' ${output_dir}G_closest_feat.bed |sed 's/gene_name=//g'|awk 'OFS="\t"{if($1!="."){if($NF==0){print $2,$3,$4,$5,$6}else{if($NF<0){if($NF>=-5000){print $2,$3,$4,$5,$6}}else{if($NF<=5000){print $2,$3,$4,$5,$6}}}}}'|tee ${output_dir}G_TSS_with_peak_at_0.tsv |awk 'function abs(v){return v < 0 ? -v : v}{OFS="\t";last=abs($NF);if($NF!=0){print $0,last}}' |sort -k4,4 -k6,6n -k5,5n|sort -k4,4 -u --merge|cut -f1-5 >${output_dir}G_TSS_with_peak_at_less_more_0.tsv

grep -P "\t0" ${output_dir}G_TSS_with_peak_at_0.tsv|sort -u -k4,4 >${output_dir}G_TSS_with_peak_at_0.tmp && mv ${output_dir}G_TSS_with_peak_at_0.tmp ${output_dir}G_TSS_with_peak_at_0.tsv

cut -f4 ${output_dir}G_TSS_with_peak_at_0.tsv |awk 'OFS="\t"{print "\t"$1"\t"}' >${output_dir}G_IDs_at_0.tsv

grep -v -f ${output_dir}G_IDs_at_0.tsv ${output_dir}G_TSS_with_peak_at_less_more_0.tsv >${output_dir}G_TSS_with_peak_at_less_more_0.tmp && mv ${output_dir}G_TSS_with_peak_at_less_more_0.tmp ${output_dir}G_TSS_with_peak_at_less_more_0.tsv

echo -e "chromosome\ttranscript_start\ttranscript_strand\tgene_name\tdist_peak_to_transcript_loc" >${output_dir}G_TSS_with_closest_peak.tsv

#concatenate the results
cat ${output_dir}G_TSS_with_peak_at_0.tsv ${output_dir}G_TSS_with_peak_at_less_more_0.tsv >>${output_dir}G_TSS_with_closest_peak.tsv

#retrieve all dist to TSS
cut -f5 ${output_dir}G_TSS_with_closest_peak.tsv >${output_dir}G_dist_to_TSS_around_5kb_2kb.tsv



${bedtools} closest -a ${G_bis} -b ${output_dir}all_TSS.gff -t first -D a -sorted -nonamecheck >${output_dir}G_allPeaks_with_closest_feat.bed || { echo "closestbed  failure 1 !" 1>&2; exit; }

my_num=$(cat ${output_dir}G_allPeaks_with_closest_feat.bed|awk 'OFS="\t"{if($18!="."){if($NF==0){print $18}else{if($NF<0){if($NF>=-5000){print $18}}else{if($NF<=5000){print $18}}}}}' |wc -l|awk '{print $1}')

echo -e "chromosome\tpeak_pos\tpeak_ID\tpeak_start\tpeak_end\tTSS_start\tgene_name\tTSS_dist" >${output_dir}G_annotated_peaks.tsv
cut -f1,4,9,14,18,19 ${output_dir}G_allPeaks_with_closest_feat.bed|sed 's/\+/ /g'|sed 's/-_-/\t/g' >>${output_dir}G_annotated_peaks.tsv
awk -F'\t' 'OFS="\t"{if(NR==1){print}else{b=";"split($7,a,";");for(i=1;i<=length(a);i++){if(a[i]~/gene_name/){b=a[i]}};gene_name=gensub("gene_name=","","G",b);$7=gene_name;split($3,peak_id," ");$3=peak_id[4];print $0}}' ${output_dir}G_annotated_peaks.tsv >${output_dir}G_annotated_peaks.tmp && mv ${output_dir}G_annotated_peaks.tmp ${output_dir}G_annotated_peaks.tsv

echo -e "G\t$my_num" >>${output_dir}peaks_around_TSS_Minus5kb_Plus5kb.tsv




###########################






##look for closest L peak to TSS 5kb/2kb
${bedtools} closest -a ${output_dir}all_TSS.gff -b ${L} -t first -D a -sorted -nonamecheck >${output_dir}L_closest_feat.bed || { echo "closestbed  failure 1 !" 1>&2; exit; }

#taking the transcript with the closest TSS to peak (priority to dist 0, then the closest in absolute number ; priority to negative dist if the same absolute dist is present)
awk 'OFS="\t"{split($9,a,";");b="";for(i=1;i<=length(a);i++){if(a[i]~/gene_name/){b=a[i]}};print $10,$1,$4,$7,b,$NF}' ${output_dir}L_closest_feat.bed |sed 's/gene_name=//g'|awk 'OFS="\t"{if($1!="."){if($NF==0){print $2,$3,$4,$5,$6}else{if($NF<0){if($NF>=-5000){print $2,$3,$4,$5,$6}}else{if($NF<=5000){print $2,$3,$4,$5,$6}}}}}'|tee ${output_dir}L_TSS_with_peak_at_0.tsv |awk 'function abs(v){return v < 0 ? -v : v}{OFS="\t";last=abs($NF);if($NF!=0){print $0,last}}' |sort -k4,4 -k6,6n -k5,5n|sort -k4,4 -u --merge|cut -f1-5 >${output_dir}L_TSS_with_peak_at_less_more_0.tsv

grep -P "\t0" ${output_dir}L_TSS_with_peak_at_0.tsv|sort -u -k4,4 >${output_dir}L_TSS_with_peak_at_0.tmp && mv ${output_dir}L_TSS_with_peak_at_0.tmp ${output_dir}L_TSS_with_peak_at_0.tsv

cut -f4 ${output_dir}L_TSS_with_peak_at_0.tsv |awk 'OFS="\t"{print "\t"$1"\t"}' >${output_dir}L_IDs_at_0.tsv

grep -v -f ${output_dir}L_IDs_at_0.tsv ${output_dir}L_TSS_with_peak_at_less_more_0.tsv >${output_dir}L_TSS_with_peak_at_less_more_0.tmp && mv ${output_dir}L_TSS_with_peak_at_less_more_0.tmp ${output_dir}L_TSS_with_peak_at_less_more_0.tsv

echo -e "chromosome\ttranscript_start\ttranscript_strand\tgene_name\tdist_peak_to_transcript_loc" >${output_dir}L_TSS_with_closest_peak.tsv

#concatenate the results
cat ${output_dir}L_TSS_with_peak_at_0.tsv ${output_dir}L_TSS_with_peak_at_less_more_0.tsv >>${output_dir}L_TSS_with_closest_peak.tsv

#retrieve all dist to TSS
cut -f5 ${output_dir}L_TSS_with_closest_peak.tsv >${output_dir}L_dist_to_TSS_around_5kb_2kb.tsv

${bedtools} closest -a ${L_bis} -b ${output_dir}all_TSS.gff -t first -D a -sorted -nonamecheck >${output_dir}L_allPeaks_with_closest_feat.bed || { echo "closestbed  failure 1 !" 1>&2; exit; }

my_num=$(cat ${output_dir}L_allPeaks_with_closest_feat.bed|awk 'OFS="\t"{if($18!="."){if($NF==0){print $18}else{if($NF<0){if($NF>=-5000){print $18}}else{if($NF<=5000){print $18}}}}}'|wc -l|awk '{print $1}')

echo -e "chromosome\tpeak_pos\tpeak_ID\tpeak_start\tpeak_end\tTSS_start\tgene_name\tTSS_dist" >${output_dir}L_annotated_peaks.tsv
cut -f1,4,9,14,18,19 ${output_dir}L_allPeaks_with_closest_feat.bed|sed 's/\+/ /g'|sed 's/-_-/\t/g' >>${output_dir}L_annotated_peaks.tsv
awk -F'\t' 'OFS="\t"{if(NR==1){print}else{b=";"split($7,a,";");for(i=1;i<=length(a);i++){if(a[i]~/gene_name/){b=a[i]}};gene_name=gensub("gene_name=","","G",b);$7=gene_name;split($3,peak_id," ");$3=peak_id[4];print $0}}' ${output_dir}L_annotated_peaks.tsv >${output_dir}L_annotated_peaks.tmp && mv ${output_dir}L_annotated_peaks.tmp ${output_dir}L_annotated_peaks.tsv

echo -e "L\t$my_num" >>${output_dir}peaks_around_TSS_Minus5kb_Plus5kb.tsv

############################




##look for closest H peak to TSS 5kb/2kb
${bedtools} closest -a ${output_dir}all_TSS.gff -b ${H} -t first -D a -sorted -nonamecheck >${output_dir}H_closest_feat.bed || { echo "closestbed  failure 1 !" 1>&2; exit; }

#taking the transcript with the closest TSS to peak (priority to dist 0, then the closest in absolute number ; priority to negative dist if the same absolute dist is present)
awk 'OFS="\t"{split($9,a,";");b="";for(i=1;i<=length(a);i++){if(a[i]~/gene_name/){b=a[i]}};print $10,$1,$4,$7,b,$NF}' ${output_dir}H_closest_feat.bed |sed 's/gene_name=//g'|awk 'OFS="\t"{if($1!="."){if($NF==0){print $2,$3,$4,$5,$6}else{if($NF<0){if($NF>=-5000){print $2,$3,$4,$5,$6}}else{if($NF<=5000){print $2,$3,$4,$5,$6}}}}}'|tee ${output_dir}H_TSS_with_peak_at_0.tsv |awk 'function abs(v){return v < 0 ? -v : v}{OFS="\t";last=abs($NF);if($NF!=0){print $0,last}}' |sort -k4,4 -k6,6n -k5,5n|sort -k4,4 -u --merge|cut -f1-5 >${output_dir}H_TSS_with_peak_at_less_more_0.tsv

grep -P "\t0" ${output_dir}H_TSS_with_peak_at_0.tsv|sort -u -k4,4 >${output_dir}H_TSS_with_peak_at_0.tmp && mv ${output_dir}H_TSS_with_peak_at_0.tmp ${output_dir}H_TSS_with_peak_at_0.tsv

cut -f4 ${output_dir}H_TSS_with_peak_at_0.tsv |awk 'OFS="\t"{print "\t"$1"\t"}' >${output_dir}H_IDs_at_0.tsv

grep -v -f ${output_dir}H_IDs_at_0.tsv ${output_dir}H_TSS_with_peak_at_less_more_0.tsv >${output_dir}H_TSS_with_peak_at_less_more_0.tmp && mv ${output_dir}H_TSS_with_peak_at_less_more_0.tmp ${output_dir}H_TSS_with_peak_at_less_more_0.tsv

echo -e "chromosome\ttranscript_start\ttranscript_strand\tgene_name\tdist_peak_to_transcript_loc" >${output_dir}H_TSS_with_closest_peak.tsv

#concatenate the results
cat ${output_dir}H_TSS_with_peak_at_0.tsv ${output_dir}H_TSS_with_peak_at_less_more_0.tsv >>${output_dir}H_TSS_with_closest_peak.tsv

#retrieve all dist to TSS
cut -f5 ${output_dir}H_TSS_with_closest_peak.tsv >${output_dir}H_dist_to_TSS_around_5kb_2kb.tsv

${bedtools} closest -a ${H_bis} -b ${output_dir}all_TSS.gff -t first -D a -sorted -nonamecheck >${output_dir}H_allPeaks_with_closest_feat.bed || { echo "closestbed  failure 1 !" 1>&2; exit; }

my_num=$(cat ${output_dir}H_allPeaks_with_closest_feat.bed|awk 'OFS="\t"{if($18!="."){if($NF==0){print $18}else{if($NF<0){if($NF>=-5000){print $18}}else{if($NF<=5000){print $18}}}}}' |wc -l|awk '{print $1}')

echo -e "chromosome\tpeak_pos\tpeak_ID\tpeak_start\tpeak_end\tTSS_start\tgene_name\tTSS_dist" >${output_dir}H_annotated_peaks.tsv
cut -f1,4,9,14,18,19 ${output_dir}H_allPeaks_with_closest_feat.bed|sed 's/\+/ /g'|sed 's/-_-/\t/g' >>${output_dir}H_annotated_peaks.tsv
awk -F'\t' 'OFS="\t"{if(NR==1){print}else{b=";"split($7,a,";");for(i=1;i<=length(a);i++){if(a[i]~/gene_name/){b=a[i]}};gene_name=gensub("gene_name=","","G",b);$7=gene_name;split($3,peak_id," ");$3=peak_id[4];print $0}}' ${output_dir}H_annotated_peaks.tsv >${output_dir}H_annotated_peaks.tmp && mv ${output_dir}H_annotated_peaks.tmp ${output_dir}H_annotated_peaks.tsv

echo -e "H\t$my_num" >>${output_dir}peaks_around_TSS_Minus5kb_Plus5kb.tsv

##############################

##look for closest P peak to TSS 5kb/2kb
${bedtools} closest -a ${output_dir}all_TSS.gff -b ${P} -t first -D a -sorted -nonamecheck >${output_dir}P_closest_feat.bed || { echo "closestbed  failure 1 !" 1>&2; exit; }

#taking the transcript with the closest TSS to peak (priority to dist 0, then the closest in absolute number ; priority to negative dist if the same absolute dist is present)
awk 'OFS="\t"{split($9,a,";");b="";for(i=1;i<=length(a);i++){if(a[i]~/gene_name/){b=a[i]}};print $10,$1,$4,$7,b,$NF}' ${output_dir}P_closest_feat.bed |sed 's/gene_name=//g'|awk 'OFS="\t"{if($1!="."){if($NF==0){print $2,$3,$4,$5,$6}else{if($NF<0){if($NF>=-5000){print $2,$3,$4,$5,$6}}else{if($NF<=5000){print $2,$3,$4,$5,$6}}}}}'|tee ${output_dir}P_TSS_with_peak_at_0.tsv |awk 'function abs(v){return v < 0 ? -v : v}{OFS="\t";last=abs($NF);if($NF!=0){print $0,last}}' |sort -k4,4 -k6,6n -k5,5n|sort -k4,4 -u --merge|cut -f1-5 >${output_dir}P_TSS_with_peak_at_less_more_0.tsv

grep -P "\t0" ${output_dir}P_TSS_with_peak_at_0.tsv|sort -u -k4,4 >${output_dir}P_TSS_with_peak_at_0.tmp && mv ${output_dir}P_TSS_with_peak_at_0.tmp ${output_dir}P_TSS_with_peak_at_0.tsv

cut -f4 ${output_dir}P_TSS_with_peak_at_0.tsv |awk 'OFS="\t"{print "\t"$1"\t"}' >${output_dir}P_IDs_at_0.tsv

grep -v -f ${output_dir}P_IDs_at_0.tsv ${output_dir}P_TSS_with_peak_at_less_more_0.tsv >${output_dir}P_TSS_with_peak_at_less_more_0.tmp && mv ${output_dir}P_TSS_with_peak_at_less_more_0.tmp ${output_dir}P_TSS_with_peak_at_less_more_0.tsv

echo -e "chromosome\ttranscript_start\ttranscript_strand\tgene_name\tdist_peak_to_transcript_loc" >${output_dir}P_TSS_with_closest_peak.tsv

#concatenate the results
cat ${output_dir}P_TSS_with_peak_at_0.tsv ${output_dir}P_TSS_with_peak_at_less_more_0.tsv >>${output_dir}P_TSS_with_closest_peak.tsv

#retrieve all dist to TSS
cut -f5 ${output_dir}P_TSS_with_closest_peak.tsv >${output_dir}P_dist_to_TSS_around_5kb_2kb.tsv


cut -f 4 ${output_dir}G_TSS_with_closest_peak.tsv >${output_dir}G_genes.txt

cut -f 4 ${output_dir}L_TSS_with_closest_peak.tsv >${output_dir}L_genes.txt

grep -F -f ${output_dir}G_genes.txt ${output_dir}L_genes.txt >${output_dir}G_L_genes.txt


${bedtools} closest -a ${P_bis} -b ${output_dir}all_TSS.gff -t first -D a -sorted -nonamecheck >${output_dir}P_allPeaks_with_closest_feat.bed || { echo "closestbed  failure 1 !" 1>&2; exit; }

my_num=$(cat ${output_dir}P_allPeaks_with_closest_feat.bed|awk 'OFS="\t"{if($18!="."){if($NF==0){print $18}else{if($NF<0){if($NF>=-5000){print $18}}else{if($NF<=5000){print $18}}}}}' |wc -l|awk '{print $1}')

echo -e "chromosome\tpeak_pos\tpeak_ID\tpeak_start\tpeak_end\tTSS_start\tgene_name\tTSS_dist" >${output_dir}P_annotated_peaks.tsv
cut -f1,4,9,14,18,19 ${output_dir}P_allPeaks_with_closest_feat.bed|sed 's/\+/ /g'|sed 's/-_-/\t/g' >>${output_dir}P_annotated_peaks.tsv
awk -F'\t' 'OFS="\t"{if(NR==1){print}else{b=";"split($7,a,";");for(i=1;i<=length(a);i++){if(a[i]~/gene_name/){b=a[i]}};gene_name=gensub("gene_name=","","G",b);$7=gene_name;split($3,peak_id," ");$3=peak_id[4];print $0}}' ${output_dir}P_annotated_peaks.tsv >${output_dir}P_annotated_peaks.tmp && mv ${output_dir}P_annotated_peaks.tmp ${output_dir}P_annotated_peaks.tsv

echo -e "P\t$my_num" >>${output_dir}peaks_around_TSS_Minus5kb_Plus5kb.tsv


#######################



######## end of dev ########"


########## make normalize bedgraph #######

#names
all_chip_names=("G" "H" "P" "L")

#tmp_chip="/media/marcgabriel/39f07759-1445-4caf-bf52-06e8bb75d6f3/chip_tmp/"

tmp_chip="/media/marcgabriel/eda138bc-95c2-4778-99df-8915815cb86e/Marina_hotair_lsd1/"

#bam rep 1 IP
rep1_all_IP_bam_files=("/media/marcgabriel/Maxtor/Lsd1/Chipseq/mappings/A685C9.unique.bam"
                       "/media/marcgabriel/Maxtor/Lsd1/Chipseq/mappings/A685C11.unique.bam"
                       "/media/marcgabriel/Maxtor/Lsd1/Chipseq/mappings/A685C13.unique.bam"
                       "/media/marcgabriel/Maxtor/Lsd1/Chipseq/mappings/A685C15.unique.bam")

#bam rep 2 IP                   
rep2_all_IP_bam_files=("/media/marcgabriel/Maxtor/Lsd1/Chipseq/mappings/A685C10.unique.bam"
                       "/media/marcgabriel/Maxtor/Lsd1/Chipseq/mappings/A685C12.unique.bam"
                       "/media/marcgabriel/Maxtor/Lsd1/Chipseq/mappings/A685C14.unique.bam"
                       "/media/marcgabriel/Maxtor/Lsd1/Chipseq/mappings/A685C16.unique.bam")

#bam rep 1 input                     
rep1_all_input_bam_files=("/media/marcgabriel/Maxtor/Lsd1/Chipseq/mappings/A684C1.unique.bam"
                          "/media/marcgabriel/Maxtor/Lsd1/Chipseq/mappings/A684C3.unique.bam"
                          "/media/marcgabriel/Maxtor/Lsd1/Chipseq/mappings/A684C5.unique.bam"
                          "/media/marcgabriel/Maxtor/Lsd1/Chipseq/mappings/A684C7.unique.bam")

#bam rep 2 input                       
rep2_all_input_bam_files=("/media/marcgabriel/Maxtor/Lsd1/Chipseq/mappings/A684C2.unique.bam"
                          "/media/marcgabriel/Maxtor/Lsd1/Chipseq/mappings/A684C4.unique.bam"
                          "/media/marcgabriel/Maxtor/Lsd1/Chipseq/mappings/A684C6.unique.bam"
                          "/media/marcgabriel/Maxtor/Lsd1/Chipseq/mappings/A684C8.unique.bam")


#for i in $(seq 0 $((${#all_chip_names[*]}-1)));do


            #if [ -d "${tmp_chip}${all_chip_names[$i]}_subscripts_dir/" ];then
            
               #rm -rf "${tmp_chip}${all_chip_names[$i]}_subscripts_dir/"
               
            #fi
            
            #mkdir "${tmp_chip}${all_chip_names[$i]}_subscripts_dir/"
            
            ##if [ ! -f ${tmp_chip}${all_chip_names[$i]}_mean_ratio_norm_IP_norm_input.bedgraph.gz ];then
            
				#echo -e "#!/bin/bash" >${tmp_chip}${all_chip_names[$i]}_subscripts_dir/${all_chip_names[$i]}_subscript_1.sh
				#echo -e "#!/bin/bash" >${tmp_chip}${all_chip_names[$i]}_subscripts_dir/${all_chip_names[$i]}_subscript_2.sh
				#echo -e "#!/bin/bash" >${tmp_chip}${all_chip_names[$i]}_subscripts_dir/${all_chip_names[$i]}_subscript_3.sh
				#echo -e "#!/bin/bash" >${tmp_chip}${all_chip_names[$i]}_subscripts_dir/${all_chip_names[$i]}_subscript_4.sh
				
				#echo -e "#!/bin/bash" >${tmp_chip}${all_chip_names[$i]}_subscripts_dir/${all_chip_names[$i]}_subscript_5.sh
				
				#chmod 755 ${tmp_chip}${all_chip_names[$i]}_subscripts_dir/${all_chip_names[$i]}_subscript_1.sh
				#chmod 755 ${tmp_chip}${all_chip_names[$i]}_subscripts_dir/${all_chip_names[$i]}_subscript_2.sh
				#chmod 755 ${tmp_chip}${all_chip_names[$i]}_subscripts_dir/${all_chip_names[$i]}_subscript_3.sh
				#chmod 755 ${tmp_chip}${all_chip_names[$i]}_subscripts_dir/${all_chip_names[$i]}_subscript_4.sh
				
				#chmod 755 ${tmp_chip}${all_chip_names[$i]}_subscripts_dir/${all_chip_names[$i]}_subscript_5.sh

				
				#echo -e "chip name is : ${all_chip_names[$i]}"

				#echo -e "${all_chip_names[$i]}_rep1_IP\t$($samtools view -c ${rep1_all_IP_bam_files[$i]})" >${tmp_chip}${all_chip_names[$i]}_summary.txt
				#read_num=0
				#read_num=$(grep "${all_chip_names[$i]}_rep1_IP" ${tmp_chip}${all_chip_names[$i]}_summary.txt|cut -f2)
				#if [[ $read_num -lt 1 ]];then echo -e "no reads for ${rep1_all_IP_bam_files[$i]} !!";exit 1;fi
				#echo -e "\t- read num. for rep 1 IP  : ${read_num}"
				#echo -e "${bedtools} genomecov -d -split -ibam ${rep1_all_IP_bam_files[$i]}|awk -v read_num=$read_num 'OFS=\"\\\t\"{print \$1,\$2-1,\$2,((\$3)*1000000)/read_num}' | gzip >${tmp_chip}${all_chip_names[$i]}_cov_rep1_IP.bedgraph.gz ||{ echo \"bedGraph 1 failure !\" 1>&2; exit; }" >>${tmp_chip}${all_chip_names[$i]}_subscripts_dir/${all_chip_names[$i]}_subscript_1.sh
				
				#echo -e "${all_chip_names[$i]}_rep2_IP\t$($samtools view -c ${rep2_all_IP_bam_files[$i]})" >>${tmp_chip}${all_chip_names[$i]}_summary.txt
				#read_num=0
				#read_num=$(grep "${all_chip_names[$i]}_rep2_IP" ${tmp_chip}${all_chip_names[$i]}_summary.txt|cut -f2)
				#if [[ $read_num -lt 1 ]];then echo -e "no reads for ${rep2_all_IP_bam_files[$i]} !!";exit 1;fi
				#echo -e "\t- read num. for rep 2 IP  : ${read_num}"
				#echo -e "${bedtools} genomecov -d -split -ibam ${rep2_all_IP_bam_files[$i]}|awk -v read_num=$read_num 'OFS=\"\\\t\"{print \$1,\$2-1,\$2,((\$3)*1000000)/read_num}' | gzip >${tmp_chip}${all_chip_names[$i]}_cov_rep2_IP.bedgraph.gz||{ echo \"bedGraph 2 failure !\" 1>&2; exit; }" >>${tmp_chip}${all_chip_names[$i]}_subscripts_dir/${all_chip_names[$i]}_subscript_2.sh
				
				#echo -e "${all_chip_names[$i]}_rep1_input\t$($samtools view -c ${rep1_all_input_bam_files[$i]})" >>${tmp_chip}${all_chip_names[$i]}_summary.txt
				#read_num=0
				#read_num=$(grep "${all_chip_names[$i]}_rep1_input" ${tmp_chip}${all_chip_names[$i]}_summary.txt|cut -f2)
				#if [[ $read_num -lt 1 ]];then echo -e "no reads for ${rep1_all_input_bam_files[$i]} !!";exit 1;fi
				#echo -e "\t- read num. for rep 1 input  : ${read_num}"
				#echo -e "${bedtools} genomecov -d -split -ibam ${rep1_all_input_bam_files[$i]}|awk -v read_num=$read_num 'OFS=\"\\\t\"{print \$1,\$2-1,\$2,((\$3)*1000000)/read_num}' | gzip >${tmp_chip}${all_chip_names[$i]}_cov_rep1_input.bedgraph.gz||{ echo \"bedGraph 3 failure !\" 1>&2; exit; }" >>${tmp_chip}${all_chip_names[$i]}_subscripts_dir/${all_chip_names[$i]}_subscript_3.sh
				
				#echo -e "${all_chip_names[$i]}_rep2_input\t$($samtools view -c ${rep2_all_input_bam_files[$i]})" >>${tmp_chip}${all_chip_names[$i]}_summary.txt
				#read_num=0
				#read_num=$(grep "${all_chip_names[$i]}_rep2_input" ${tmp_chip}${all_chip_names[$i]}_summary.txt|cut -f2)
				#if [[ $read_num -lt 1 ]];then echo -e "no reads for ${rep2_all_input_bam_files[$i]} !!";exit 1;fi
				#echo -e "\t- read num. for rep 2 input  : ${read_num}"
				#echo -e "${bedtools} genomecov -d -split -ibam ${rep2_all_input_bam_files[$i]}|awk -v read_num=$read_num 'OFS=\"\\\t\"{print \$1,\$2-1,\$2,((\$3)*1000000)/read_num}' | gzip >${tmp_chip}${all_chip_names[$i]}_cov_rep2_input.bedgraph.gz||{ echo \"bedGraph 4 failure !\" 1>&2; exit; }" >>${tmp_chip}${all_chip_names[$i]}_subscripts_dir/${all_chip_names[$i]}_subscript_4.sh
				
				
				#echo -e " find ${tmp_chip}${all_chip_names[$i]}_subscripts_dir/ -name \"*_subscript*.sh\" |grep \"${all_chip_names[$i]}_\"|grep -v \"_5\\\.sh\"| xargs -n 1 -P 2 bash || { echo \"executing bash subscripts failure !\" 1>&2; exit; }" >>${tmp_chip}${all_chip_names[$i]}_subscripts_dir/${all_chip_names[$i]}_subscript_5.sh

				##make the ratio IP/input
				
				##retrieve the smaller value not equal to 0 for input 1 (we will divide this value by 4, and give the result as denominator for division when the input is 0, but not the IP)
				#echo -e " min_val_input_1=\$(zcat ${tmp_chip}${all_chip_names[$i]}_cov_rep1_input.bedgraph.gz|LANG=en_EN sort -k4,4g|head -n2|LANG=en_EN sort -k4,4gr|head -n1)" >>${tmp_chip}${all_chip_names[$i]}_subscripts_dir/${all_chip_names[$i]}_subscript_5.sh
				
				#echo -e " paste <(zcat ${tmp_chip}${all_chip_names[$i]}_cov_rep1_IP.bedgraph.gz) <(zcat ${tmp_chip}${all_chip_names[$i]}_cov_rep1_input.bedgraph.gz) |awk -v min_val_input_1=\$min_val_input_1 'OFS=\"\\\t\"{if(\$1==\$5 && \$2==\$6){if(\$4==0 && \$8==0){\$4=1 ; \$8=1}else if(\$8==0){\$8=min_val_input_1/4};print \$1,\$2,\$3,\$4/\$8}}' | gzip >${tmp_chip}${all_chip_names[$i]}_rep1_ratio_norm_IP_norm_input.gz && rm ${tmp_chip}${all_chip_names[$i]}_cov_rep1_IP.bedgraph.gz ${tmp_chip}${all_chip_names[$i]}_cov_rep1_input.bedgraph.gz" >>${tmp_chip}${all_chip_names[$i]}_subscripts_dir/${all_chip_names[$i]}_subscript_5.sh
				
				##retrieve the smaller value not equal to 0 for input 2 (we will divide this value by 4, and give the result as denominator for division when the input is 0, but not the IP)
				#echo -e " min_val_input_2=\$(zcat ${tmp_chip}${all_chip_names[$i]}_cov_rep2_input.bedgraph.gz|LANG=en_EN sort -k4,4g|head -n2|LANG=en_EN sort -k4,4gr|head -n1)" >>${tmp_chip}${all_chip_names[$i]}_subscripts_dir/${all_chip_names[$i]}_subscript_5.sh
				
				#echo -e " paste <(zcat ${tmp_chip}${all_chip_names[$i]}_cov_rep2_IP.bedgraph.gz) <(zcat ${tmp_chip}${all_chip_names[$i]}_cov_rep2_input.bedgraph.gz) |awk -v min_val_input_2=\$min_val_input_2 'OFS=\"\\\t\"{if(\$1==\$5 && \$2==\$6){if(\$4==0 && \$8==0){\$4=1 ; \$8=1}else if(\$8==0){\$8=min_val_input_2/4};print \$1,\$2,\$3,\$4/\$8}}' | gzip >${tmp_chip}${all_chip_names[$i]}_rep2_ratio_norm_IP_norm_input.gz && rm ${tmp_chip}${all_chip_names[$i]}_cov_rep2_IP.bedgraph.gz ${tmp_chip}${all_chip_names[$i]}_cov_rep2_input.bedgraph.gz" >>${tmp_chip}${all_chip_names[$i]}_subscripts_dir/${all_chip_names[$i]}_subscript_5.sh
				
				#echo -e " paste <(zcat ${tmp_chip}${all_chip_names[$i]}_rep1_ratio_norm_IP_norm_input.gz) <(zcat ${tmp_chip}${all_chip_names[$i]}_rep2_ratio_norm_IP_norm_input.gz) |awk 'OFS=\"\\\t\"{print \$1,\$2,\$3,(\$4+\$8)/2}' | gzip >${tmp_chip}${all_chip_names[$i]}_mean_ratio_norm_IP_norm_input.bedgraph.gz && rm ${tmp_chip}${all_chip_names[$i]}_rep1_ratio_norm_IP_norm_input.gz ${tmp_chip}${all_chip_names[$i]}_rep2_ratio_norm_IP_norm_input.gz" >>${tmp_chip}${all_chip_names[$i]}_subscripts_dir/${all_chip_names[$i]}_subscript_5.sh

            ##fi
            
            #echo -e " if [ ! -f ${tmp_chip}${all_chip_names[$i]}_mean_ratio_norm_IP_norm_input_reduced.bedgraph.gz ];then" >>${tmp_chip}${all_chip_names[$i]}_subscripts_dir/${all_chip_names[$i]}_subscript_5.sh
            
				##https://www.biostars.org/p/267712/
				#echo -e "   zcat ${tmp_chip}${all_chip_names[$i]}_mean_ratio_norm_IP_norm_input.bedgraph.gz|$bedtools groupby -i stdin -g 1,4 -c 2,3 -o min,max|awk 'OFS=\"\\\t\"{print \$1,\$3,\$4,\$2}'|gzip >${tmp_chip}${all_chip_names[$i]}_mean_ratio_norm_IP_norm_input_reduced.bedgraph.gz" >>${tmp_chip}${all_chip_names[$i]}_subscripts_dir/${all_chip_names[$i]}_subscript_5.sh
            
           
             #echo -e " fi" >>${tmp_chip}${all_chip_names[$i]}_subscripts_dir/${all_chip_names[$i]}_subscript_5.sh
            
			
			#echo -e "\n*****************\n"

    
#done

#find ${tmp_chip} -name "*_subscript*.sh" |grep "_5.sh"| xargs -n 1 -P 4 bash || { echo "executing bash subscripts failure !" 1>&2; exit; }


###################################

############


#add the others
cat $original_G $original_L $original_H $original_P |sort -k1,1 -k4,4n  >${output_dir}concatenated_all.gff


$bedtools merge -d 1 -c 9 -o "collapse" -delim " " -i ${output_dir}concatenated_all.gff >${output_dir}merged_all.gff




############## processing for venn for peaks grouped in regions ###############
###############################################################################

#colors for G L H P
colors="#7f7f7f,#e46c0a,#604a7b,#31859c"

not_intersected_G=$(awk -F'\t' '{if($4~/G1_/ && $4!~/L1_/ && $4!~/H1_/ && $4!~/P1_/){print}}' ${output_dir}merged_all.gff|wc -l)

not_intersected_H=$(awk -F'\t' '{if($4~/H1_/ && $4!~/G1_/ && $4!~/L1_/ && $4!~/P1_/){print}}' ${output_dir}merged_all.gff|wc -l)

not_intersected_L=$(awk -F'\t' '{if($4~/L1_/ && $4!~/H1_/ && $4!~/G1_/ && $4!~/P1_/){print}}' ${output_dir}merged_all.gff|wc -l)

not_intersected_P=$(awk -F'\t' '{if($4~/P1_/ && $4!~/G1_/ && $4!~/L1_/ && $4!~/H1_/){print}}' ${output_dir}merged_all.gff|wc -l)

echo -e "- results for not intersected G : $not_intersected_G ; L : $not_intersected_L ; H : $not_intersected_H ; P : $not_intersected_P\n"


intersection_all=$(awk -F'\t' '{if($4~/G1_/ && $4~/H1_/ && $4~/L1_/ && $4~/P1_/){print}}' ${output_dir}merged_all.gff| wc -l)

echo -e "- results for intersection all : $intersection_all\n"

#pairwise intersections :


#G serial
intersection_G_L=$(awk -F'\t' '{if($4~/G1_/ && $4~/L1_/ && $4!~/H1_/ && $4!~/P1_/){print}}' ${output_dir}merged_all.gff| wc -l)

intersection_G_H=$(awk -F'\t' '{if($4~/G1_/ && $4~/H1_/ && $4!~/L1_/ && $4!~/P1_/){print}}' ${output_dir}merged_all.gff|wc -l)

intersection_G_P=$(awk -F'\t' '{if($4~/G1_/ && $4~/P1_/ && $4!~/H1_/ && $4!~/L1_/){print}}' ${output_dir}merged_all.gff|wc -l)


#L serial 

intersection_L_H=$(awk -F'\t' '{if($4~/L1_/ && $4~/H1_/ && $4!~/P1_/ && $4!~/G1_/){print}}' ${output_dir}merged_all.gff|wc -l) 

intersection_L_P=$(awk -F'\t' '{if($4~/L1_/ && $4~/P1_/ && $4!~/G1_/ && $4!~/H1_/){print}}' ${output_dir}merged_all.gff|wc -l) 

#H serial
intersection_H_P=$(awk -F'\t' '{if($4~/H1_/ && $4~/P1_/ && $4!~/G1_/ && $4!~/L1_/){print}}' ${output_dir}merged_all.gff|wc -l) 

echo -e "- results for intersection G vs L : $intersection_G_L ; G vs H : $intersection_G_H ; G vs P : $intersection_G_P ; L vs H : $intersection_L_H ; L vs P : $intersection_L_P ; H vs P : $intersection_H_P\n"

#triple intersections
intersection_G_L_H=$(awk -F'\t' '{if($4~/G1_/ && $4~/L1_/ && $4~/H1_/ && $4!~/P1_/){print}}' ${output_dir}merged_all.gff|wc -l)

intersection_G_L_P=$(awk -F'\t' '{if($4~/G1_/ && $4~/L1_/ && $4~/P1_/ && $4!~/H1_/){print}}' ${output_dir}merged_all.gff|wc -l)
 
intersection_G_H_P=$(awk -F'\t' '{if($4~/G1_/ && $4~/H1_/ && $4~/P1_/ && $4!~/L1_/){print}}' ${output_dir}merged_all.gff|wc -l)

intersection_L_H_P=$(awk -F'\t' '{if($4~/H1_/ && $4~/L1_/ && $4~/P1_/ && $4!~/G1_/){print}}' ${output_dir}merged_all.gff|wc -l)
 
 echo -e "- results for intersections G vs L vs H : $intersection_G_L_H ; G vs L vs P : $intersection_G_L_P ; G vs H vs P : $intersection_G_H_P ; H vs L vs P : $intersection_L_H_P \n"

###########

echo -e "GL is : $intersection_G_L ; GH is : $intersection_G_H ; GP is : $intersection_G_P ; GLH is $intersection_G_L_H ; GHP is : $intersection_G_H_P ; GLP is : $intersection_G_L_P\n==================\n"



all_G=$(awk -F'\t' '{if($4~/G1_/){print}}' ${output_dir}merged_all.gff|wc -l)

all_L=$(awk -F'\t' '{if($4~/L1_/){print}}' ${output_dir}merged_all.gff|wc -l)

all_H=$(awk -F'\t' '{if($4~/H1_/){print}}' ${output_dir}merged_all.gff|wc -l)

all_P=$(awk -F'\t' '{if($4~/P1_/){print}}' ${output_dir}merged_all.gff|wc -l)


$rscript2 ${output_dir} "G,L,H,P" "$not_intersected_G,$not_intersected_L,$not_intersected_H,$not_intersected_P,$intersection_G_L,$intersection_G_H,$intersection_G_P,$intersection_L_H,$intersection_L_P,$intersection_H_P,$intersection_G_L_H,$intersection_G_L_P,$intersection_G_H_P,$intersection_L_H_P,$intersection_all" "$all_G,$all_L,$all_H,$all_P" "$colors" "peaksGroupedByRegions"

#peak dist to TSS
$getBoxplots $output_dir


######### peaks distribution across genomic features #########
##############################################################


getCleanIntervals(){

  gff_to_remove=$1

  gff_to_extract=$2
  
  prefix=$3
  
  file_name=$4
  
  my_output=$5
  
  use_prefix_as_ID=$6
  
  if [ "$gff_to_remove" != "" ];then


	  awk 'OFS="\t"{print $1,$4-1,$5,NR,".",$7}' $gff_to_remove | LANG=en_EN sort -k1,1 -k2,2n >${my_output}feat_to_remove.tmp
	  

	  awk '{if($6=="+" || $6=="."){print}}' ${my_output}feat_to_remove.tmp >${my_output}feat_to_remove_plus.tmp
	  
	  awk '{if($6=="-" || $6=="."){print}}' ${my_output}feat_to_remove.tmp >${my_output}feat_to_remove_minus.tmp
	  

	  #rm ${my_output}feat_to_remove.tmp
	  
	  
  else
  
         >${my_output}feat_to_remove_plus.tmp
         
         >${my_output}feat_to_remove_minus.tmp
  
  
  fi
  
  if [ "$prefix" == "intergenic" ];then
  
    sign_plus="."
    
    sign_minus="."
    
    cat ${my_output}feat_to_remove_minus.tmp ${my_output}feat_to_remove_plus.tmp | LANG=en_EN sort -k1,1 -k2,2n >${my_output}feat_to_remove_plus.tmp2 && mv ${my_output}feat_to_remove_plus.tmp2 ${my_output}feat_to_remove_plus.tmp && >${my_output}feat_to_remove_minus.tmp
  
  else
  
  
    sign_plus="+"
    
    sign_minus="-"
  
  
  fi
  
  
  awk 'OFS="\t"{print $1,$4-1,$5,NR,".",$7}' $gff_to_extract| LANG=en_EN sort -k1,1 -k2,2n >${my_output}feat_to_extract.bed
  
  
  if [ "$prefix" == "intergenic" ];then
  
    cat ${my_output}feat_to_extract.bed >${my_output}feat_to_extract_plus.tmp
    
    >${my_output}feat_to_extract_minus.tmp
    
  
  else
  
  
    awk '{if($6=="+"){print}}' ${my_output}feat_to_extract.bed >${my_output}feat_to_extract_plus.tmp

    awk '{if($6=="-"){print}}' ${my_output}feat_to_extract.bed >${my_output}feat_to_extract_minus.tmp
  
  fi


  
  
   if [ "$gff_to_remove" != "" ];then
   
       if [[ $(wc -l ${my_output}feat_to_extract_plus.tmp|awk '{print $1}') -gt 0 ]];then
       
        if [[ $(wc -l ${my_output}feat_to_remove_plus.tmp|awk '{print $1}') -gt 0 ]];then
        


	  $bedops --difference ${my_output}feat_to_extract_plus.tmp ${my_output}feat_to_remove_plus.tmp |$bedtools merge -i stdin|awk -v sign_plus=$sign_plus 'OFS="\t"{print $0,NR,".",sign_plus}' >${my_output}feat_to_extract_plus.txt
	  
	  
	else
	
	  cat ${my_output}feat_to_extract_plus.tmp >${my_output}feat_to_extract_plus.txt
	  
	
        fi
	
	
	
	else 
	
	  >${my_output}feat_to_extract_plus.txt
	  
	
	fi
	
	
	if [[ $(wc -l ${my_output}feat_to_extract_minus.tmp|awk '{print $1}') -gt 0 ]];then
	
	
	 if [[ $(wc -l ${my_output}feat_to_remove_minus.tmp|awk '{print $1}') -gt 0 ]];then
	 
	 

	  $bedops --difference ${my_output}feat_to_extract_minus.tmp ${my_output}feat_to_remove_minus.tmp|$bedtools merge -i stdin|awk -v sign_minus=$sign_minus 'OFS="\t"{print $0,NR,".",sign_minus}' >${my_output}feat_to_extract_minus.txt
	  
	  
	  else
	  
	   cat ${my_output}feat_to_extract_minus.tmp >${my_output}feat_to_extract_plus.txt
	   
	  
	  fi
	
	
	else 
	
	  >${my_output}feat_to_extract_minus.txt
	  
	
	fi
	
  else
  
    if [[ $(wc -l ${my_output}feat_to_extract_plus.tmp|awk '{print $1}') -gt 0 ]];then
  
      $bedtools merge -i ${my_output}feat_to_extract_plus.tmp|awk -v sign_plus=$sign_plus 'OFS="\t"{print $0,NR,".",sign_plus}' >${my_output}feat_to_extract_plus.txt 
    
    
    else
    
      >${my_output}feat_to_extract_plus.txt
      
    fi
    
    if [[ $(wc -l ${my_output}feat_to_extract_minus.tmp|awk '{print $1}') -gt 0 ]];then
    
      $bedtools merge -i ${my_output}feat_to_extract_minus.tmp|awk -v sign_minus=$sign_minus 'OFS="\t"{print $0,NR,".",sign_minus}' >${my_output}feat_to_extract_minus.txt
    
    
    else
    
      >${my_output}feat_to_extract_minus.txt
    
    
    fi
   
  
  
  fi
  
  if [ "$use_prefix_as_ID" == "T" ];then
  
    cat ${my_output}feat_to_extract_plus.txt ${my_output}feat_to_extract_minus.txt |LANG=en_EN sort -k1,1 -k2,2n |awk -v prefix=$prefix 'OFS="\t"{print $1,"empty",prefix,$2+1,$3,".",$6,".","ID="prefix}' >${my_output}${file_name} && rm ${my_output}feat_to_extract_plus.txt ${my_output}feat_to_extract_minus.txt ${my_output}feat_to_remove_plus.tmp ${my_output}feat_to_remove_minus.tmp ${my_output}feat_to_extract_plus.tmp ${my_output}feat_to_extract_minus.tmp
    
  else
  
  cat ${my_output}feat_to_extract_plus.txt ${my_output}feat_to_extract_minus.txt |LANG=en_EN sort -k1,1 -k2,2n |awk -v prefix=$prefix 'OFS="\t"{print $1,"empty",prefix,$2+1,$3,".",$6,".","ID="NR}' >${my_output}${file_name} && rm ${my_output}feat_to_extract_plus.txt ${my_output}feat_to_extract_minus.txt ${my_output}feat_to_remove_plus.tmp ${my_output}feat_to_remove_minus.tmp ${my_output}feat_to_extract_plus.tmp ${my_output}feat_to_extract_minus.tmp
  
  fi


}

############### end of dev #############

#list of features to search
#exon_number=1
features_list=(exon promoter intron intergenic five_prime_utr three_prime_utr)


start=$(date)

if [ ! -f ${output_dir}ref_annotation_refined.gff ];then

set -x

  #create introns
  LANG=en_EN sort -k1,1 -k4,4n <(grep -v "^#" $annotation) <($getIntronsByTranscripts $annotation "yes") <($getIntergenics $annotation|grep -v "^Created feature") >${output_dir}ref_annotation_refined.gff

fi



#construct the design file
for one_feat in ${features_list[*]};do

  echo -e "\n**** feature is :  $one_feat ****\n"

  if [ "$one_feat" != "promoter" ] ;then
  
         feat=$one_feat
  
         grep -P -i "\t${one_feat}\t" ${output_dir}ref_annotation_refined.gff|awk -v feat=$feat 'OFS="\t"{$3=feat;print}'|LANG=en_EN sort -k1,1 -k4,4n >${output_dir}${one_feat}.gff
         

  elif [ $one_feat == "promoter" ];then
  
       grep -P -i "\tgene\t" ${output_dir}ref_annotation_refined.gff| awk 'OFS="\t"{if($7=="+"){a=$4-1;$4=$4-1001;if($4<=0){$4=1};$5=a}else if($7=="-"){$4=$5+1;if($4<=0){$4=1};$5=$5+1001};$2="promoter";print}' |LANG=en_EN sort -k1,1 -k4,4n >${output_dir}${one_feat}.gff
       
     

  fi
  
 
done


echo -e "\ntry to reduce features...\n"


cat "${output_dir}five_prime_utr.gff" "${output_dir}three_prime_utr.gff" >${output_dir}remove_in_exon.gff

cat ${output_dir}exon.gff ${output_dir}five_prime_utr.gff ${output_dir}three_prime_utr.gff >${output_dir}remove_in_intron.gff

cat ${output_dir}five_prime_utr.gff ${output_dir}three_prime_utr.gff ${output_dir}exon.gff ${output_dir}intron.gff >${output_dir}remove_in_promoter.gff


#create non overlapping features

getCleanIntervals ${output_dir}remove_in_intron.gff ${output_dir}intron.gff "intron" "intron.gff" ${output_dir} "F"

getCleanIntervals ${output_dir}remove_in_exon.gff ${output_dir}exon.gff "exon" "exon.gff" ${output_dir} "F"

getCleanIntervals ${output_dir}remove_in_promoter.gff ${output_dir}promoter.gff "promoter" "promoter.gff" ${output_dir} "F"

getCleanIntervals ${output_dir}promoter.gff ${output_dir}intergenic.gff "intergenic" "intergenic.gff" ${output_dir} "F"

getCleanIntervals "" ${output_dir}five_prime_utr.gff "five_prime_utr" "five_prime_utr.gff" ${output_dir} "F"

getCleanIntervals "" ${output_dir}three_prime_utr.gff "three_prime_utr" "three_prime_utr.gff" ${output_dir} "F"




all_marks=($original_G $original_L $original_P $original_H)

all_names=("G" "L" "P" "H")

echo -e "all marks are : ${all_marks[*]}\n"

features_list=(exon promoter intron intergenic five_prime_utr three_prime_utr)

k=1

for i in ${all_marks[*]};do

        mark_name=${all_names[$((k-1))]}
        
        all_feature_counts="${output_dir}all_${mark_name}_counts.tsv"
        
        >${all_feature_counts}
        
        echo -e "- mark is : $mark_name\n"

	for feature in ${features_list[*]};do
	
	                               #set -x
	
	
	                 echo -e "\t-> feature is : $feature\n"
	                             
				
					 cat ${i} |sed 's/ "/=/g'|sed 's/"; /;/g'|sed 's/"//g' |$bedtools intersect -sorted -a stdin -b ${output_dir}${feature}.gff -loj -nobuf -nonamecheck | LANG=en_EN sort -u -k9,9| awk -F'\t' 'OFS="\t"{if($10=="."){a="F"}else{a="T"};print $9,a}' >${output_dir}${mark_name}_${feature}_counts.txt || { echo "bedtools intersect failure !!" 1>&2; exit; }
			
				
					if [[ $(wc -l ${all_feature_counts}|awk '{print $1}') -gt 1 ]];then
					
					  echo -e "\t-> files to join are : ${output_dir}${mark_name}_${feature}_counts.txt & ${all_feature_counts}\n"
				
					  LANG=en_EN join -t $'\t' -11 -21 ${output_dir}${mark_name}_${feature}_counts.txt ${all_feature_counts} >${all_feature_counts}.tmp && mv ${all_feature_counts}.tmp ${all_feature_counts}
					  
					  nb_peaks=$(wc -l ${all_feature_counts}|awk '{print $1}')
					  
					  echo -e "\t-> nb peaks = $nb_peaks"
					  
					  
				
					else
					
					   echo -e "\t-> this is the first turn ($feature), no files to join, next !\n"
				
					   cat ${output_dir}${mark_name}_${feature}_counts.txt >${all_feature_counts}
					   
                      nb_peaks=$(wc -l ${all_feature_counts}|awk '{print $1}')
					  
					  echo -e "\t-> nb peaks = $nb_peaks"
					   
				
					fi
					
					#set +x
					
					
					if [ "$feature" == "${features_list[$((${#features_list[*]}-1))]}" ];then
										
					
					  cat <(echo -e "ID ${features_list[*]} mark_type"|sed 's/ /\t/g') ${all_feature_counts}|awk -v mark_type=$mark_name 'OFS="\t"{if(NR==1){print}else{print $0,mark_type}}' >${all_feature_counts}.tmp && mv ${all_feature_counts}.tmp ${all_feature_counts}
					
					  nb_peaks=$(wc -l ${all_feature_counts}|awk '{print $1}')
					  
					  echo -e "\t-> nb peaks = $nb_peaks"

					
					
					fi
				
	done
	
	#set -x
	
	$classifyPeaks ${output_dir} ${all_feature_counts}
	
	
	#set +x
	
	echo -e "\n======== next mark ! ============\n"
	
	k=$((k+1))

done


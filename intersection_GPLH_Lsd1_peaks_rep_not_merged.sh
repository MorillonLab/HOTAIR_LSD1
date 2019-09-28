#!/bin/bash


##### input data #######

export samtools="samtools"

#use bedtools 2.27, with 2.26, groupby is buggy !!
export bedtools="/home/marcgabriel/Desktop/bedtools2/bin/bedtools"

#bedops program https://bedops.readthedocs.io/en/latest/
export bedops="/home/marcgabriel/Desktop/bedops-2.4.30/bin/bedops"

#annotation in gff format
export annotation="/home/marcgabriel/Documents/gencode26lift37/gencode.v26lift37.annotation.gff3"

#output dir
output_dir="/media/marcgabriel/Transcend/LSD1_metagenes/Marina_peaks_analysis_rep_not_merged/"

### peaks in bed format (2 per condition) ###
#G="/media/marcgabriel/Transcend/LSD1_metagenes/original_peaks_not_merged/A685C9.10.merged_filteredBR.peaks.gff"
G="/media/marcgabriel/Transcend/LSD1_metagenes/original_peaks_not_merged/A685C9.unique-W200-G1000-islands-summary-FDR0.05_filteredBR.bed /media/marcgabriel/Transcend/LSD1_metagenes/original_peaks_not_merged/A685C10.unique-W200-G1000-islands-summary-FDR0.05_filteredBR.bed"

#L="/media/marcgabriel/Transcend/LSD1_metagenes/original_peaks_not_merged/A685C15.16.merged_filteredBR.peaks.gff"
L="/media/marcgabriel/Transcend/LSD1_metagenes/original_peaks_not_merged/A685C15.unique-W200-G1000-islands-summary-FDR0.05_filteredBR.bed /media/marcgabriel/Transcend/LSD1_metagenes/original_peaks_not_merged/A685C16.unique-W200-G1000-islands-summary-FDR0.05_filteredBR.bed"

#H="/media/marcgabriel/Transcend/LSD1_metagenes/original_peaks_not_merged/A685C11.12.merged_filteredBR.peaks.gff"
H="/media/marcgabriel/Transcend/LSD1_metagenes/original_peaks_not_merged/A685C11.unique-W200-G1000-islands-summary-FDR0.05_filteredBR.bed /media/marcgabriel/Transcend/LSD1_metagenes/original_peaks_not_merged/A685C12.unique-W200-G1000-islands-summary-FDR0.05_filteredBR.bed"

#P="/media/marcgabriel/Transcend/LSD1_metagenes/original_peaks_not_merged/A685C13.14.merged_filteredBR.peaks.gff"

P="/media/marcgabriel/Transcend/LSD1_metagenes/original_peaks_not_merged/A685C13.unique-W200-G1000-islands-summary-FDR0.05_filteredBR.bed /media/marcgabriel/Transcend/LSD1_metagenes/original_peaks_not_merged/A685C14.unique-W200-G1000-islands-summary-FDR0.05_filteredBR.bed"


#create intergenic features
export getIntergenics=$(dirname "$0")/getIntergenics.R

#create intronic features
export getIntronsByTranscripts=$(dirname "$0")/getIntronsByTranscripts.R

#classify the peaks
export classifyPeaks=$(dirname "$0")/classifyPeaks.R


rscript1=$(dirname "$0")/venn_diagrams_with_numerics.R
chmod 755 $rscript1

rscript2=$(dirname "$0")/venn_chipseq_GHPL.R
chmod 755 $rscript2

rscript3=$(dirname "$0")/venn_lsd1_diff.R
chmod 755 $rscript3

getBoxplots=$(dirname "$0")/chipeseq_GHPL_boxplots_dist_to_tss.R
chmod 755 $getBoxplots

getCoveredNuc=$(dirname "$0")/getCoveredNuc.R
chmod 755 $getCoveredNuc

#upregulated genes in high and low migration phenotype
LMS_HMS_file=$(dirname "$0")/HMS_LMS.xlsx

#gene IDs & corresponding gene types (2 columns, tab-delimited)
gene_type=$(dirname "$0")/matching_gene_name_types.tsv

###############



output_dir="${output_dir}/"
output_dir=$(echo $output_dir |sed 's/\/\//\//g')

if [ ! -d "$output_dir" ];then


  mkdir $output_dir

fi


grep -v "^#" $annotation |grep -P "\ttranscript\t"|awk 'OFS="\t"{if($7=="+"){$5=$4}else{$4=$5};print}' | sort -k1,1 -k4,4n >${output_dir}all_TSS.gff


#######################################################


######### processing of bed loc files (dev) ###############

#compute midpeak : |awk 'BEGIN {OFS="\t";printf($1"\t");printf($2"\t");printf($3"\t");printf("%.0f\t", ($4+$5)/2);printf("%.0f\t", ($4+$5)/2);print $6,$7,$8,$9 }'

#concatenate the peaks (from all rep), then merge them

cat <(cat $(echo $G|awk '{print $1}')|awk 'OFS="\t"{$4="G1_peak_number_"NR;print}') <(cat $(echo $G|awk '{print $2}')|awk 'OFS="\t"{$4="G2_peak_number_"NR;print}')|sort -k1,1 -k2,2n >${output_dir}G_combined_peaks.tmp && $bedtools merge -d 1 -c 4 -o "distinct" -i ${output_dir}G_combined_peaks.tmp | awk 'OFS="\t"{print $1,"annotation","peaks",$2+1,$3,".","+",".",$4}'|awk 'OFS="\t"{printf($1"\t");printf("annotation""\t");printf("peaks""\t");printf("%.0f\t", $4);printf("%.0f\t", $5);print ".","+",".",$9}' >${output_dir}G_combined_peaks.tsv
G_merged_rep=${output_dir}G_combined_peaks.tsv

cat $G |sort -k1,1 -k2,2n|awk 'OFS="\t"{print $0,"+"}' | awk 'OFS="\t"{print $1,"annotation","peaks",$2+1,$3,".","+",".","G1_peak_number_"NR}'|awk '{OFS="\t";printf($1"\t");printf($2"\t");printf($3"\t");printf("%.0f\t", $4);printf("%.0f\t", $5);print $6,$7,$8,$9 }' >${output_dir}G_peaks_all_rep.gff3
G="${output_dir}G_peaks_all_rep.gff3"


cat <(cat $(echo $L|awk '{print $1}')|awk 'OFS="\t"{$4="L1_peak_number_"NR;print}') <(cat $(echo $L|awk '{print $2}')|awk 'OFS="\t"{$4="L2_peak_number_"NR;print}')|sort -k1,1 -k2,2n >${output_dir}L_combined_peaks.tmp && $bedtools merge -d 1 -c 4 -o "distinct" -i ${output_dir}L_combined_peaks.tmp | awk 'OFS="\t"{print $1,"annotation","peaks",$2+1,$3,".","+",".",$4}'|awk 'OFS="\t"{printf($1"\t");printf("annotation""\t");printf("peaks""\t");printf("%.0f\t", $4);printf("%.0f\t", $5);print ".","+",".",$9}' >${output_dir}L_combined_peaks.tsv
L_merged_rep=${output_dir}L_combined_peaks.tsv

cat $L |sort -k1,1 -k2,2n |awk 'OFS="\t"{print $0,"+"}' | awk 'OFS="\t"{print $1,"annotation","peaks",$2+1,$3,".","+",".","L1_peak_number_"NR}'|awk '{OFS="\t";printf($1"\t");printf($2"\t");printf($3"\t");printf("%.0f\t", $4);printf("%.0f\t", $5);print $6,$7,$8,$9 }' >${output_dir}L_peaks_all_rep.gff3
L="${output_dir}L_peaks_all_rep.gff3"


cat <(cat $(echo $H|awk '{print $1}')|awk 'OFS="\t"{$4="H1_peak_number_"NR;print}') <(cat $(echo $H|awk '{print $2}')|awk 'OFS="\t"{$4="H2_peak_number_"NR;print}')|sort -k1,1 -k2,2n >${output_dir}H_combined_peaks.tmp && $bedtools merge -d 1 -c 4 -o "distinct" -i ${output_dir}H_combined_peaks.tmp | awk 'OFS="\t"{print $1,"annotation","peaks",$2+1,$3,".","+",".",$4}'|awk 'OFS="\t"{printf($1"\t");printf("annotation""\t");printf("peaks""\t");printf("%.0f\t", $4);printf("%.0f\t", $5);print ".","+",".",$9}' >${output_dir}H_combined_peaks.tsv
H_merged_rep=${output_dir}H_combined_peaks.tsv

cat $H |sort -k1,1 -k2,2n |awk 'OFS="\t"{print $0,"+"}' | awk 'OFS="\t"{print $1,"annotation","peaks",$2+1,$3,".","+",".","H1_peak_number_"NR}'|awk '{OFS="\t";printf($1"\t");printf($2"\t");printf($3"\t");printf("%.0f\t", $4);printf("%.0f\t", $5);print $6,$7,$8,$9 }' >${output_dir}H_peaks_all_rep.gff3
H="${output_dir}H_peaks_all_rep.gff3"


cat <(cat $(echo $P|awk '{print $1}')|awk 'OFS="\t"{$4="P1_peak_number_"NR;print}') <(cat $(echo $P|awk '{print $2}')|awk 'OFS="\t"{$4="P2_peak_number_"NR;print}')|sort -k1,1 -k2,2n >${output_dir}P_combined_peaks.tmp && $bedtools merge -d 1 -c 4 -o "distinct" -i ${output_dir}P_combined_peaks.tmp | awk 'OFS="\t"{print $1,"annotation","peaks",$2+1,$3,".","+",".",$4}'|awk 'OFS="\t"{printf($1"\t");printf("annotation""\t");printf("peaks""\t");printf("%.0f\t", $4);printf("%.0f\t", $5);print ".","+",".",$9}' >${output_dir}P_combined_peaks.tsv
P_merged_rep=${output_dir}P_combined_peaks.tsv

cat $P |sort -k1,1 -k2,2n |awk 'OFS="\t"{print $0,"+"}' | awk 'OFS="\t"{print $1,"annotation","peaks",$2+1,$3,".","+",".","P1_peak_number_"NR}'|awk '{OFS="\t";printf($1"\t");printf($2"\t");printf($3"\t");printf("%.0f\t", $4);printf("%.0f\t", $5);print $6,$7,$8,$9 }' >${output_dir}P_peaks_all_rep.gff3
P="${output_dir}P_peaks_all_rep.gff3"


original_G=$G
original_L=$L
original_H=$H
original_P=$P

####### compute number of covered nucleotide per peak

awk -F'\t' 'OFS="\t"{print "G",($5-$4)+1}' $original_G >${output_dir}covered_nuc.tsv
awk  -F'\t' 'OFS="\t"{print "H",($5-$4)+1}' $original_H >>${output_dir}covered_nuc.tsv
awk  -F'\t' 'OFS="\t"{print "P",($5-$4)+1}' $original_P >>${output_dir}covered_nuc.tsv
awk  -F'\t' 'OFS="\t"{print "L",($5-$4)+1}' $original_L >>${output_dir}covered_nuc.tsv

$getCoveredNuc ${output_dir}


cat $G |sort -k1,1 -k4,4n|awk 'OFS="\t"{print $0,"+"}' |awk 'OFS="\t"{printf($1"\t");printf("annotation""\t");printf("peaks""\t");printf("%.0f\t", $4);printf("%.0f\t", $5);print ".","+",".","G1_peak_number_"NR}' | sort -k1,1 -k4,4n >${output_dir}G_peaks_all_rep.tmp && mv ${output_dir}G_peaks_all_rep.tmp ${output_dir}G_peaks_all_rep.gff3
G="${output_dir}G_peaks_all_rep.gff3"

cat $H |sort -k1,1 -k4,4n|awk 'OFS="\t"{print $0,"+"}' |awk 'OFS="\t"{printf($1"\t");printf("annotation""\t");printf("peaks""\t");printf("%.0f\t", $4);printf("%.0f\t", $5);print ".","+",".","H1_peak_number_"NR}' |sort -k1,1 -k4,4n >${output_dir}H_peaks_all_rep.tmp && mv ${output_dir}H_peaks_all_rep.tmp ${output_dir}H_peaks_all_rep.gff3
H="${output_dir}H_peaks_all_rep.gff3"


cat $L |sort -k1,1 -k4,4n|awk 'OFS="\t"{print $0,"+"}' |awk 'OFS="\t"{printf($1"\t");printf("annotation""\t");printf("peaks""\t");printf("%.0f\t", $4);printf("%.0f\t", $5);print ".","+",".","L1_peak_number_"NR}' |sort -k1,1 -k4,4n >${output_dir}L_peaks_all_rep.tmp && mv ${output_dir}L_peaks_all_rep.tmp ${output_dir}L_peaks_all_rep.gff3

L="${output_dir}L_peaks_all_rep.gff3"

cat $P |sort -k1,1 -k4,4n|awk 'OFS="\t"{print $0,"+"}' |awk 'OFS="\t"{printf($1"\t");printf("annotation""\t");printf("peaks""\t");printf("%.0f\t", $4);printf("%.0f\t", $5);print ".","+",".","P1_peak_number_"NR}' |sort -k1,1 -k4,4n >${output_dir}P_peaks_all_rep.tmp && mv ${output_dir}P_peaks_all_rep.tmp ${output_dir}P_peaks_all_rep.gff3
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



######### closest peaks to TSS ####

echo -e "condition\tnumber" >${output_dir}peaks_around_TSS_Minus5kb_Plus5kb.tsv


##look for closest G peak to TSS -5kb/2kb
${bedtools} closest -a ${output_dir}all_TSS.gff -b ${G_merged_rep} -t first -D a -sorted -nonamecheck >${output_dir}G_closest_feat.bed || { echo "closestbed  failure 1 !" 1>&2; exit; }

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

my_num=$(cat ${output_dir}G_allPeaks_with_closest_feat.bed|awk 'OFS="\t"{if($18!="."){if($NF==0){print $18}else{if($NF<0){if($NF>=-5000){print $18}}else{if($NF<=5000){print $18}}}}}'|wc -l|awk '{print $1}')

echo -e "chromosome\tpeak_pos\tpeak_ID\tpeak_start\tpeak_end\tTSS_start\tgene_name\tTSS_dist" >${output_dir}G_annotated_peaks.tsv
cut -f1,4,9,14,18,19 ${output_dir}G_allPeaks_with_closest_feat.bed|sed 's/\+/ /g'|sed 's/-_-/\t/g' >>${output_dir}G_annotated_peaks.tsv
awk -F'\t' 'OFS="\t"{if(NR==1){print}else{b=";"split($7,a,";");for(i=1;i<=length(a);i++){if(a[i]~/gene_name/){b=a[i]}};gene_name=gensub("gene_name=","","G",b);$7=gene_name;split($3,peak_id," ");$3=peak_id[4];print $0}}' ${output_dir}G_annotated_peaks.tsv >${output_dir}G_annotated_peaks.tmp && mv ${output_dir}G_annotated_peaks.tmp ${output_dir}G_annotated_peaks.tsv

echo -e "G\t$my_num" >>${output_dir}peaks_around_TSS_Minus5kb_Plus5kb.tsv

awk 'OFS="\t"{split($9,a,";");b="";for(i=1;i<=length(a);i++){if(a[i]~/gene_name/){b=a[i]}};print $1,$4,$7,b,$NF,$(NF-1)}' ${output_dir}G_closest_feat.bed |sed 's/gene_name=//g'|awk 'function abs(v){return v < 0 ? -v : v}{OFS="\t";if($5==0){print}else{if(abs($5)<=5000){print}}}' |awk '{if($NF~/G1/ && $NF~/G2/){print}}' |sort -u -k4,4|cat <(echo -e "chromosome\ttranscript_start\ttranscript_strand\tgene_name\tdist_peak_to_transcript_loc\tpeaks_group") - >${output_dir}G_genes_with_peaks_in_both_rep.tsv

awk 'OFS="\t"{split($9,a,";");b="";for(i=1;i<=length(a);i++){if(a[i]~/gene_name/){b=a[i]}};print $1,$4,$7,b,$NF,$(NF-1)}' ${output_dir}G_closest_feat.bed |sed 's/gene_name=//g'|awk 'function abs(v){return v < 0 ? -v : v}{OFS="\t";if($5==0){print}else{if(abs($5)<=5000){print}}}' |awk '{if(($NF~/G1/ && $NF!~/G2/) || ($NF~/G2/ && $NF!~/G1/)){print}}' |sort -u -k4,4|cat <(echo -e "chromosome\ttranscript_start\ttranscript_strand\tgene_name\tdist_peak_to_transcript_loc\tpeaks_group") - >${output_dir}G_genes_with_peaks_in_just_one_rep.tsv 


###########################






##look for closest L peak to TSS 5kb/2kb
${bedtools} closest -a ${output_dir}all_TSS.gff -b ${L_merged_rep} -t first -D a -sorted -nonamecheck >${output_dir}L_closest_feat.bed || { echo "closestbed  failure 1 !" 1>&2; exit; }

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

awk 'OFS="\t"{split($9,a,";");b="";for(i=1;i<=length(a);i++){if(a[i]~/gene_name/){b=a[i]}};print $1,$4,$7,b,$NF,$(NF-1)}' ${output_dir}L_closest_feat.bed |sed 's/gene_name=//g'|awk 'function abs(v){return v < 0 ? -v : v}{OFS="\t";if($5==0){print}else{if(abs($5)<=5000){print}}}' |awk '{if($NF~/L1/ && $NF~/L2/){print}}' |sort -u -k4,4|cat <(echo -e "chromosome\ttranscript_start\ttranscript_strand\tgene_name\tdist_peak_to_transcript_loc\tpeaks_group") - >${output_dir}L_genes_with_peaks_in_both_rep.tsv

awk 'OFS="\t"{split($9,a,";");b="";for(i=1;i<=length(a);i++){if(a[i]~/gene_name/){b=a[i]}};print $1,$4,$7,b,$NF,$(NF-1)}' ${output_dir}L_closest_feat.bed |sed 's/gene_name=//g'|awk 'function abs(v){return v < 0 ? -v : v}{OFS="\t";if($5==0){print}else{if(abs($5)<=5000){print}}}' |awk '{if(($NF~/L1/ && $NF!~/L2/) || ($NF~/L2/ && $NF!~/L1/)){print}}' |sort -u -k4,4|cat <(echo -e "chromosome\ttranscript_start\ttranscript_strand\tgene_name\tdist_peak_to_transcript_loc\tpeaks_group") - >${output_dir}L_genes_with_peaks_in_just_one_rep.tsv

############################


##look for closest H peak to TSS 5kb/2kb
${bedtools} closest -a ${output_dir}all_TSS.gff -b ${H_merged_rep} -t first -D a -sorted -nonamecheck >${output_dir}H_closest_feat.bed || { echo "closestbed  failure 1 !" 1>&2; exit; }

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

my_num=$(cat ${output_dir}H_allPeaks_with_closest_feat.bed|awk 'OFS="\t"{if($18!="."){if($NF==0){print $18}else{if($NF<0){if($NF>=-5000){print $18}}else{if($NF<=5000){print $18}}}}}'|wc -l|awk '{print $1}')

echo -e "chromosome\tpeak_pos\tpeak_ID\tpeak_start\tpeak_end\tTSS_start\tgene_name\tTSS_dist" >${output_dir}H_annotated_peaks.tsv
cut -f1,4,9,14,18,19 ${output_dir}H_allPeaks_with_closest_feat.bed|sed 's/\+/ /g'|sed 's/-_-/\t/g' >>${output_dir}H_annotated_peaks.tsv
awk -F'\t' 'OFS="\t"{if(NR==1){print}else{b=";"split($7,a,";");for(i=1;i<=length(a);i++){if(a[i]~/gene_name/){b=a[i]}};gene_name=gensub("gene_name=","","G",b);$7=gene_name;split($3,peak_id," ");$3=peak_id[4];print $0}}' ${output_dir}H_annotated_peaks.tsv >${output_dir}H_annotated_peaks.tmp && mv ${output_dir}H_annotated_peaks.tmp ${output_dir}H_annotated_peaks.tsv

echo -e "H\t$my_num" >>${output_dir}peaks_around_TSS_Minus5kb_Plus5kb.tsv

awk 'OFS="\t"{split($9,a,";");b="";for(i=1;i<=length(a);i++){if(a[i]~/gene_name/){b=a[i]}};print $1,$4,$7,b,$NF,$(NF-1)}' ${output_dir}H_closest_feat.bed |sed 's/gene_name=//g'|awk 'function abs(v){return v < 0 ? -v : v}{OFS="\t";if($5==0){print}else{if(abs($5)<=5000){print}}}' |awk '{if($NF~/H1/ && $NF~/H2/){print}}' |sort -u -k4,4|cat <(echo -e "chromosome\ttranscript_start\ttranscript_strand\tgene_name\tdist_peak_to_transcript_loc\tpeaks_group") - >${output_dir}H_genes_with_peaks_in_both_rep.tsv

awk 'OFS="\t"{split($9,a,";");b="";for(i=1;i<=length(a);i++){if(a[i]~/gene_name/){b=a[i]}};print $1,$4,$7,b,$NF,$(NF-1)}' ${output_dir}H_closest_feat.bed |sed 's/gene_name=//g'|awk 'function abs(v){return v < 0 ? -v : v}{OFS="\t";if($5==0){print}else{if(abs($5)<=5000){print}}}' |awk '{if(($NF~/H1/ && $NF!~/H2/) || ($NF~/H2/ && $NF!~/H1/)){print}}' |sort -u -k4,4|cat <(echo -e "chromosome\ttranscript_start\ttranscript_strand\tgene_name\tdist_peak_to_transcript_loc\tpeaks_group") - >${output_dir}H_genes_with_peaks_in_just_one_rep.tsv


##############################

##look for closest P peak to TSS 5kb/2kb
${bedtools} closest -a ${output_dir}all_TSS.gff -b ${P_merged_rep} -t first -D a -sorted -nonamecheck >${output_dir}P_closest_feat.bed || { echo "closestbed  failure 1 !" 1>&2; exit; }

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

my_num=$(cat ${output_dir}P_allPeaks_with_closest_feat.bed|awk 'OFS="\t"{if($18!="."){if($NF==0){print $18}else{if($NF<0){if($NF>=-5000){print $18}}else{if($NF<=5000){print $18}}}}}'|wc -l|awk '{print $1}')

echo -e "chromosome\tpeak_pos\tpeak_ID\tpeak_start\tpeak_end\tTSS_start\tgene_name\tTSS_dist" >${output_dir}P_annotated_peaks.tsv
cut -f1,4,9,14,18,19 ${output_dir}P_allPeaks_with_closest_feat.bed|sed 's/\+/ /g'|sed 's/-_-/\t/g' >>${output_dir}P_annotated_peaks.tsv
awk -F'\t' 'OFS="\t"{if(NR==1){print}else{b=";"split($7,a,";");for(i=1;i<=length(a);i++){if(a[i]~/gene_name/){b=a[i]}};gene_name=gensub("gene_name=","","G",b);$7=gene_name;split($3,peak_id," ");$3=peak_id[4];print $0}}' ${output_dir}P_annotated_peaks.tsv >${output_dir}P_annotated_peaks.tmp && mv ${output_dir}P_annotated_peaks.tmp ${output_dir}P_annotated_peaks.tsv

echo -e "P\t$my_num" >>${output_dir}peaks_around_TSS_Minus5kb_Plus5kb.tsv

awk 'OFS="\t"{split($9,a,";");b="";for(i=1;i<=length(a);i++){if(a[i]~/gene_name/){b=a[i]}};print $1,$4,$7,b,$NF,$(NF-1)}' ${output_dir}P_closest_feat.bed |sed 's/gene_name=//g'|awk 'function abs(v){return v < 0 ? -v : v}{OFS="\t";if($5==0){print}else{if(abs($5)<=5000){print}}}' |awk '{if($NF~/P1/ && $NF~/P2/){print}}' |sort -u -k4,4|cat <(echo -e "chromosome\ttranscript_start\ttranscript_strand\tgene_name\tdist_peak_to_transcript_loc\tpeaks_group") - >${output_dir}P_genes_with_peaks_in_both_rep.tsv

awk 'OFS="\t"{split($9,a,";");b="";for(i=1;i<=length(a);i++){if(a[i]~/gene_name/){b=a[i]}};print $1,$4,$7,b,$NF,$(NF-1)}' ${output_dir}P_closest_feat.bed |sed 's/gene_name=//g'|awk 'function abs(v){return v < 0 ? -v : v}{OFS="\t";if($5==0){print}else{if(abs($5)<=5000){print}}}' |awk '{if(($NF~/P1/ && $NF!~/P2/) || ($NF~/P2/ && $NF!~/P1/)){print}}' |sort -u -k4,4|cat <(echo -e "chromosome\ttranscript_start\ttranscript_strand\tgene_name\tdist_peak_to_transcript_loc\tpeaks_group") - >${output_dir}P_genes_with_peaks_in_just_one_rep.tsv



$rscript3 $output_dir $annotation $LMS_HMS_file $gene_type

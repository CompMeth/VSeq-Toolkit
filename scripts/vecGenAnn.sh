outDir=${1}
bedtools=${2}
cd $outDir
for f in *split.txt.cluster.1
do
	sed 's/^ *//g'  $f | awk '{print ($1"\t"$2"\t"$2+1"\t"$3"\t"$4)}' > file2
	$bedtools closest -a file2 -b anno1 -t first > file3
	#RefSeqID,GeneStrand,GeneName,GeneLength,DistTSS,Upstream,Downstream,IntronExon
	awk '{print ($1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$8-$7) }' file3 > file4
	awk '{ if ($2>=$7 && $2<=$8 && $10 ~ /+/) print ($1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$2-$7-1) ; if ($2>=$7 && $2<=$8 && $10 ~ /-/) print ($1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$8-$2); if ($2<$7 || $2>$8) print ($1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$10)}' file4  > file5
	awk '{if($2<$7 && $10 ~ /+/) print ($1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$7-$2+1); if ($2>$8 && $10 ~ /-/) print ($1"\t"$2"\t" $3"\t" $4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$2-$8) ; if ($2>=$7 && $2<=$8 && $10 ~ /+/ || $2>=$7 && $2<=$8 && $10 ~ /-/) print ($1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$10); if ($2>$8 && $10 ~ /+/) print ($1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18,$10); if ($2<$7 && $10~ /-/) print ($1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$10)}' file5  > file6
	awk '{ if ($2>$8 && $10   ~ /+/) print ($1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$2-$8); if ($2<$7 && $10 ~ /-/) print ($1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$7-$2+1); if ($2>=$7 && $2<=$8 && $10 ~ /+/ || $2>=$7 && $2<=$8 && $10 ~ /-/) print ($1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$10) ; if ($2<$7 && $10  ~ /+/) print ($1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$10) ; if ($2>$8 && $10  ~ /-/) print ($1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$10)  } ' file6 > file7
	awk '{ n14 = split($14, t14, ",");n15 = split($15, t15, ","); for (i = 0; ++i < n14;) {print $1,t14[i],t15[i], $10=="+" ? $9 "_Exon" i : $9 "_Exon" n14-i, $10} }' file7 | awk '{print $1"\t"$2"\t"$3"\t"$4"\t""0""\t"$5}' > fileexons.bed
	$bedtools intersect -a fileexons.bed  -b file7  -wb > file8
	awk '{n14 = split($14, t14, ",");n15 = split($15, t15, ",");for (i = 0; ++i < n14-1;) {m14= n14-1; print $1,t15[i],t14[i+1], $10=="+" ? $9 "_Intron" i : $9 "_Intron"  m14-i, $10}}' file7 | awk '{print $1"\t"$2"\t"$3"\t"$4"\t""0""\t"$5}' > fileintrons.bed 
	$bedtools intersect -a fileintrons.bed -b file7 -wb > file9
	sort -k7,7 -k8,8n -k9,9n -k10,10 -k11,11 -k27,27  -u file8 > file10
	sort -k7,7 -k8,8n -k9,9n -k10,10 -k11,11 -k27,27  -u file9 >> file10
	awk '{print $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $13 "\t" $14 "\t" $15 "\t" $16 "\t" $17 "\t" $18 "\t" $19 "\t" $20 "\t" $21 "\t" $22 "\t" $23 "\t" $24 "\t" $25 "\t" $26 "\t" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 }' file10 > file10.1
	echo "NA" >> file10.1
	awk -F' ' 'NR==FNR{c[$1$2]++;next};c[$1$2] == 0' file10.1 file7 > file10.2
	sed -i '$ d' file10.1 
	cat file10.1 file10.2 | sort -k1,1 -k2,2 |cut -f1,2,4,5,9,10,16,17,18,19,20,24 > file11 
	cut -f1,2,3,4,5,6,7,8  file11 > file11a
	cut -f 9,10,11,12 file11 | sed 's/-/ /g'  | sed 's/+/ /g' > file11b
	paste file11a file11b > file11c
	awk '{if ($8!="0") print $0; if ($8=="0") print $1"\t"$2"\t"$3"\t"$4"\t""NA""\t""NA""\t""NA""\t""NA""\t""NA""\t""NA""\t""NA""\t""NA""\t"}' file11c | sed 's/\t/,/g' > $f.anno
done

#This script contains modified version of a part of script from published method Afzal S et al. Mol Ther Nucleic Acid. 2017;6(March):133-139. doi:10.1016/j.omtn.2016.12.001

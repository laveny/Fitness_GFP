#!/bin/bash

Main=$1
RAW=$2
ref=$3
min_length=$4
max_length=$5
call_py=$6
needle_mut_py=$7
cluster_major_py=$8
create_all_sim_bar_illu_py=$9
  
filename=$(basename "$RAW")
filename="${filename%.subreads.*}"

Main=${Main}/${filename}/
mkdir ${Main}
cd ${Main}
Log='pipeline.log'

a=`date`
cat > ${Log} <<EOL
## This file is to mark-down steps for pipeline.sh
## This pipeline runs for pacbio analysis for 202408_OLD_DMS_DATA
Date:${a}
...
EOL


mkdir 0.Logs
mkdir 1.blasr
mkdir 2.ccs
mkdir 3.call_needle_mut


cat >> ${Log} <<EOL
RAWDATA_PATH : ${RAW}
REF_PATH : ${ref}
NAMEBASE : ${filename}
...
EOL

# 1. check md5 + blasr + .sam

nohup md5sum ${RAW} > ./0.Logs/md5check &

a=`date`
cat >> ${Log} <<EOL
1.start running blasr
${a}
...
EOL

blasr --bam --nproc 50 --unaligned ./1.blasr/${filename}.blasr.unaligned.bam \
--out ./1.blasr/${filename}.blasr.t50.bam \
${RAW} ${ref} > ./0.Logs/blasr.log 2>&1

a=`date`
cat >> ${Log} <<EOL
blasr finished.
2.managing the output of blasr
${a}
EOL


variable=./1.blasr/${filename}
samtools view ./1.blasr/${filename}.blasr.t50.bam | awk -v var="$variable" '
  {
	count[$2]++
	if ($2 == 16) {
		print $1 > var".blasr.16.txt"
	} else if ($2 == 0){
		print $1 > var".blasr.0.txt"
	} else {
		print $1 > var".blasr.unuse.txt"
	}
  }
  END {
	for (val in count){	
	    print val, count[val] > var".count.txt"
	}
  }

' -

a=`date`
cat >> ${Log} <<EOL
managing finished

3. start separating negative/positive subreads
${a}
EOL


## 1.2 select raw subreads from rawdata(.sam)
samtools view ${RAW} -H > ${variable}.header.sam

time samtools view ${RAW} |awk 'NR==FNR{a[$1]; next} $1 in a' ${variable}.blasr.16.txt -  > ${variable}.n.subreads.sam

time samtools view ${RAW} |awk 'NR==FNR{a[$1]; next} $1 in a' ${variable}.blasr.0.txt -  > ${variable}.p.subreads.sam

cat ${variable}.header.sam ${variable}.n.subreads.sam | samtools view -bS - > ${variable}.n.subreads.bam

cat ${variable}.header.sam ${variable}.p.subreads.sam | samtools view -bS - > ${variable}.p.subreads.bam

samtools sort -@ 20 -O bam -o ${variable}.p.sorted.subreads.bam ${variable}.p.subreads.bam
samtools sort -@ 20 -O bam -o ${variable}.n.sorted.subreads.bam ${variable}.n.subreads.bam

# 2. run ccs for negative and positive subreads separately

a=`date`
cat >> ${Log} <<EOL
separating finished

4. start CCS
${a}
EOL

ccs --min-length "$min_length" --max-length "$max_length" -j 30 --min-passes 5 ${variable}.p.sorted.subreads.bam ./2.ccs/${filename}.p.ccs.bam

ccs --min-length "$min_length" --max-length "$max_length" -j 30 --min-passes 5 ${variable}.n.sorted.subreads.bam ./2.ccs/${filename}.n.ccs.bam

## 3.3 transfer ccs output into .sam format

samtools view -@ 50 -o ./2.ccs/${filename}.p.ccs.sam ./2.ccs/${filename}.p.ccs.bam

samtools view -@ 50 -o ./2.ccs/${filename}.n.ccs.sam ./2.ccs/${filename}.n.ccs.bam

a=`date`
cat >> ${Log} <<EOL
ccs finished

5.start seq&barcode calling
${a}
EOL


python ${call_py} ${Main}/2.ccs/ ${filename} ${Main}/3.call_needle_mut/

a=`date`
cat >> ${Log} << EOL
calling finished

6.start needle and cal mut
${a}
EOL

python ${needle_mut_py} ${Main}/3.call_needle_mut/${filename}.call.txt ${ref} ${Main}/3.call_needle_mut/${filename}.needle.mut.txt

a=`date`
cat >> ${Log} <<EOL
all is finished

done
${a}
EOL


# Arrange needle output
awk -F"\t" '{print $3}' ${Main}/3.call_needle_mut/${filename}.needle.mut.txt |sort|uniq  > ${Main}/3.call_needle_mut/${filename}.uniq_bar.txt

awk '{print ">"$0"\n"$0}'  ${Main}/3.call_needle_mut/${filename}.uniq_bar.txt > ${Main}/3.call_needle_mut/${filename}.uniq_bar.fa

awk -F"\t" '{print $3"\t"$4}' ${Main}/3.call_needle_mut/${filename}.needle.mut.txt |sort|uniq -c > ${Main}/3.call_needle_mut/${filename}.bar_mut_c.txt

awk -F"\t" '{print $4"\t"$5"\t"$6}' ${Main}/3.call_needle_mut/${filename}.needle.mut.txt |sort|uniq > ${Main}/3.call_needle_mut/${filename}.mut_pro_stop.txt 


# using slidesort to find call barcode clusters
~/SlideSort/SS_v2/slidesort_v2 -d 1 -i ${Main}/3.call_needle_mut/${filename}.uniq_bar.fa > ${Main}/3.call_needle_mut/${filename}.slidesort.ss.txt

python ${cluster_major_py} ${Main}/3.call_needle_mut/ ${filename}


## Final Output

awk -F"\t" 'NR==FNR {a[$1] = $2"\t"$3; next} {if ($2 in a){print $0"\t"a[$2]}}' ${Main}/3.call_needle_mut/${filename}.mut_pro_stop.txt ${Main}/3.call_needle_mut/${filename}.PASS.txt > ${Main}/3.call_needle_mut/${filename}.PASS_pro.txt


## Create all similar barcodes based on PacBio result for illumina clustering

python ${create_all_sim_bar_illu_py} ${filename} ${Main}/3.call_needle_mut/
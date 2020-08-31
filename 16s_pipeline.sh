#!/usr/sh

#############################
# 
# Simple 16s rRNA analysis pipeline
#
# Main tools vsearch
#
# Use:
# chmod 775 16s_pipeline.sh
# sh 16s_pipeline.sh &
#
# Author: Jiabao Cao
############################
set -eu

#需处理的数据在data文件夹并已经去过barcode和引物
#建立临时文件及结果文件
mkdir temp
mkdir result

##工具和数据库位置
tools_dir="/software_users/caojiabao"
database_dir="/software_users/caojiabao/database/RDP_database"

##1. 质控
##合并双端序列
for i in `ls data/*_1.fq | cut -d/ -f 2 | sed 's/_1.fq//'`;
do
vsearch --fastq_mergepairs data/${i}_1.fq --reverse data/${i}_2.fq \
--fastqout temp/${i}.merged.fq
done

##需注意
##质量控制序列长度因测序而异
for i in `ls temp/*.merged.fq`;
do 
vsearch --fastq_filter $i --fastq_maxee 1 --fastq_maxlen 430 --fastaout temp/${i%*.merged.fq*}_qc.fasta;
done

###抽平
minseqsnum=`ls temp/*_qc.fasta | awk '{print "grep \">\" "$1" | wc -l"}' | less | bash | sort -n | head -n1`

for i in `ls temp/*_qc.fasta | cut -d/ -f 2 | sed 's/_qc.fasta//'`;
do 
usearch -fastx_subsample temp/${i}_qc.fasta -sample_size ${minseqsnum} -fastaout temp/${i}.sub${minseqsnum}.fa 
done

###relabel
for i in `ls temp/*.sub*.fa`;
do 
python $tools_dir/python_scripts/relabelFasta_file.py $i temp/${i%*_qc.sub*.fa}_relabel.fa;
done


##2.构建otu
###pool sequences
cat temp/*sub*.fa > temp/all.filtered.fa

####去冗余
vsearch --derep_fulllength temp/all.filtered.fa --output temp/unique.fa --sizeout

####去嵌合体
vsearch --uchime_ref temp/unique.fa --db $database_dir/rdp_16s_v16.fa --nonchimeras temp/unique.nonchimeras.fasta --threads 20
###排序
vsearch --sortbysize temp/unique.nonchimeras.fasta --minsize 2 --output temp/unique.nonchimeras.sorted.fasta

#聚类并写入result文件里
vsearch --cluster_smallmem temp/unique.nonchimeras.sorted.fasta --usersort --id 0.97 --consout result/unique.nonchimeras.sorted.rep_set.fasta --threads 20

#rename
awk '/^>/{print ">OTU_" ++i; next}{print}' < result/unique.nonchimeras.sorted.rep_set.fasta > result/unique.nonchimeras.sorted.rep_set_relabel.fasta

#map reads
vsearch --usearch_global temp/all.filtered.fasta --db result/unique.nonchimeras.sorted.rep_set_relabel.fasta --strand both \
--id 0.97 --uc result/samples.map.uc --threads 20

##创建otu表格
python2 $tools_dir/python_scripts/uc2otutab.py result/samples.map.uc > result/samples.map.otu_table.txt


###注释
java -Xmx8g -jar $tools_dir/RDPTools/classifier.jar classify -c 0.80 -f filterbyconf -o sample.map.tax_table.txt \
temp/unique.nonchimeras.sorted.rep_set_relabel.fasta

#删除临时文件
rm -r temp


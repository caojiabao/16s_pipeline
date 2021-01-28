#1. 生成manifes file
echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" >> manifest
ls data/ |grep _1.fq.gz |sed 's/\_1.fq.gz//'|awk '{print $1"\t$PWD/data/"$1"_1.fq.gz\t$PWD/data/"$1"_2.fq.gz"}' >> manifest

#2. import data
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest \
  --output-path paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33V2

qiime demux summarize --i-data paired-end-demux.qza --o-visualization paired-end-demux.qzv
#qiime tools view paired-end-demux.qzv

#3. 去除引物(如果知道引物序列)
#qiime cutadapt trim-paired \
#--i-demultiplexed-sequences paired-end-demux.qza --p-front-f GTGCCAGCMGCCGCGG \
#--p-front-r CCGTCAATTCMTTTRAGTTT \
#--o-trimmed-sequences paired-end-demux-primer.qza \
#--verbose &> primer_trimming.log

#4.1 dada2(去噪并合并)(220根据实际情况进行筛选)---使用dada2后无需在使用dublur去噪
qiime dada2 denoise-paired \
--i-demultiplexed-seqs paired-end-demux.qza \
--p-trunc-len-f 220 \
--p-trunc-len-r 220 \
--o-table table.qza \
--o-representative-sequences rep-seqs.qza \
--o-denoising-stats denoising-stats.qza

#4.2 deblur去噪##############################
#合并双端序列并进行质量控制
qiime vsearch join-pairs \
  --i-demultiplexed-seqs paired-end-demux.qza \
  --o-joined-sequences demux-joined.qza \
  --p-threads 4
#qiime demux summarize --i-data demux-joined.qza --o-visualization demux-joined.qzv

qiime quality-filter q-score \
  --i-demux demux-joined.qza \
  --o-filtered-sequences demux-joined-filtered.qza \
  --o-filter-stats demux-joined-filter-stats.qza

#qiime demux summarize --i-data demux-joined-filtered.qza --o-visualization demux-joined-filtered.qzv
#qiime metadata tabulate --m-input-file demux-joined-filter-stats.qza --o-visualization demux-joined-filter-stats.qzv

#去噪deblur(包括了去除嵌合体过滤丰度序列)
qiime deblur denoise-16S \
  --i-demultiplexed-seqs demux-joined-filtered.qza \
  --p-trim-length 400 \
  --p-sample-stats \
  --o-representative-sequences rep-seqs.qza \
  --o-table table.qza \
  --o-stats deblur-stats.qza \
  --p-jobs-to-start 30

#分别生成可视化qzv文件：
#qiime deblur visualize-stats --i-deblur-stats deblur-stats.qza --o-visualization deblur-stats.qzv 
#qiime feature-table summarize --i-table table.qza --o-visualization table.qzv
#qiime feature-table tabulate-seqs --i-data rep-seqs.qza --o-visualization rep-seqs.qzv 
###################################

#5. vsearch 对特征表进行聚类
qiime vsearch cluster-features-de-novo \
  --i-table table.qza \
  --i-sequences rep-seqs.qza \
  --p-perc-identity 0.97 \
  --o-clustered-table table-dn-97.qza \
  --o-clustered-sequences rep-seqs-dn-97.qza \
  --p-threads 30

#6.导出特征表
mkdir phyloseq
qiime tools export --input-path table-dn-97.qza --output-path phyloseq

biom convert -i phyloseq/feature-table.biom -o phyloseq/otu_table.tsv --to-tsv
cd phyloseq; sed -i '1d' otu_table.tsv
sed -i 's/#OTU ID/ASV/' otu_table.tsv

#R更改特征表名称
library(tidyverse)
#pacman::p_load(tidyverse,magrittr,stringr)
otu <- "otu_table.tsv" %>%
read.delim(check.names = FALSE,header = T,sep="\t")

rown <- paste0("ASV",seq_len(nrow(otu)))
otu[,1] <- rown
colnames(otu)[1] <- paste0("ASV",colnames(data)[1])
write.table (otu,file ="otu_table.tsv",sep ="\t",row.names = F,quote=FALSE)

# otu_table.tsv --> feature-table.biom ---> otu_table.qza
biom convert -i otu_table.tsv -o feature-table.biom --to-hdf5 --table-type="OTU table"
time qiime tools import \
  --input-path feature-table.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path otu_table-rename.qza

#7. 导出代表性序列
qiime tools export --input-path rep-seqs-dn-97.qza --output-path phyloseq
less dna-sequences.fasta |paste - -|sed '1i ASVID,seq' > rep.fa

library(pacman)
pacman::p_load(tidyverse,magrittr,stringr)
rep <- "rep.fa" %>%
read.delim(check.names = FALSE, row.names = 1) %>%
set_rownames(paste0(">ASV", seq_len(nrow(.))))
write.table (rep,file ="rep.xls", sep ="\t", row.names = T)

less rep.xls|sed '1d'|sed 's/"//g'|sed 's/\r//g'|tr "\t" "\n" > rep-seqs.fasta
#将代表性序列转换成qza格式
time qiime tools import \
--type 'FeatureData[Sequence]' \
--input-path rep-seqs.fasta \
--output-path rep-seqs-rename.qza

rm rep.fa rep.xls dna-sequences.fasta

#代表序列统计
time qiime feature-table tabulate-seqs \
--i-data rep-seqs-rename.qza \
--o-visualization raw.fq.list


#5.注释(参考数据库为silva-138)
time qiime feature-classifier classify-sklearn \
--i-classifier database/silva-138-99-nb-classifier.qza \
--i-reads phyloseq/rep-seqs-rename.qza \
--o-classification taxonomy.qza \
--p-n-jobs 30

qiime tools export \
--input-path taxonomy.qza \
--output-path phyloseq

#分级注释属水平
qiime taxa collapse \
  --i-table otu_table-rename.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table otu_table_genus.qza

qiime tools export --input-path otu_table_genus.qza --output-path phyloseq
biom convert -i phyloseq/feature-table.biom -o phyloseq/otu_genus_table.tsv --to-tsv

#6.构建进化树（比对代表性序列构建系统发育树）
time qiime phylogeny align-to-tree-mafft-fasttree \
 --i-sequences rep-seqs-rename.qza \
 --o-alignment aligned-rep-seqs.qza \
 --o-masked-alignment masked-aligned-rep-seqs.qza \
 --o-tree unrooted-tree.qza \
 --o-rooted-tree rooted-tree.qza \
 --p-n-threads 30

mkdir phyloseq

qiime tools export \
--input-path unrooted-tree.qza \
--output-path phyloseq
cd phyloseq; mv tree.nwk unrooted_tree.nwk; cd ..

qiime tools export \
--input-path rooted-tree.qza \
--output-path phyloseq
cd phyloseq; mv tree.nwk rooted_tree.nwk;cd .. 





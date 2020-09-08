#16s_pipeline2.sh

#原始文件在fastq中
#建立临时文件
#mkdir temp/
cd temp/
#1. Merge reads from one sample #path= /public/chenl/disha_data/16S
#2. filter
#3. remove chimera form quality filtered reads

for i in ../fastq/*_1.fastq;
do
g=`echo $i|cut -d '/' -f 3`
i=${g%*_1.fastq*}
/software_users/caojiabao/16s_pipeline/FLASH-1.2.11-Linux-x86_64/flash -o $i ../fastq/${i}_1.fastq ../fastq/${i}_2.fastq -t 20;
/software_users/caojiabao/16s_pipeline/fastx_toolkit/bin/fastq_quality_filter -i ${i}.extendedFrags.fastq -p 90 -q 25 -Q33|/software_users/caojiabao/16s_pipeline/fastx_toolkit/bin/fastq_to_fasta -o ${i}.fasta -Q33 
vsearch --uchime_ref ${i}.fasta --db /software_users/caojiabao/database/RDP_database/rdp_16s_v16.fa -strand plus --nonchimeras ${i}.good.fasta --threads 20
done

#4. Re-shuffle the sequence of each sample and take 3310 reads out, give a new-id
###抽平
minseqsnum=`ls *.good.fasta | awk '{print "grep \">\" "$1" | wc -l"}' | less | bash | sort -n | head -n1`

for i in `ls *.good.fasta`;
do 
perl /software_users/caojiabao/16s_pipeline/subsample.pl $i ${minseqsnum} > ${i}.sub.fasta
java -jar /software_users/caojiabao/RDPTools/classifier.jar -o ${i}.tax -f fixrank  ${i}.sub.fasta
done

#6. Cat all tax file and produce taxonomical census
cat *tax >ALL.TAX

#7. Make clean taxonomical table
#结果为不同分类单元的table
cd ../
perl /software_users/caojiabao/16s_pipeline/taxa_census_new.pl temp/ALL.TAX



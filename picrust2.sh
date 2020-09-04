#默认预测测序样品中的EC和KO丰度, 并根据预测的EC丰度推断MetaCyc途径丰度.
#文件“data/otu_table.txt”为OTU丰度表格，仅包含丰度组成信息即可，无需添加注释列。
#文件“data/otu.fasta”中包含了OTU丰度表中各OTU的代表序列。
picrust2_pipeline.py -s otu.fasta -i otu_table.txt -o picrust2_result -p 4

#determine KEGG pathway abundances from the predicted KO abundances
#生成KEGG pathway 功能分析
pathway_pipeline.py -i KO_metagenome_out/pred_metagenome_unstrat.tsv -o KEGG_pathways_out 、
--no_regroup --map picrust2/picrust2/default_files/pathway_mapfiles/KEGG_pathways_to_KO.tsv

#添加功能描述，只看KEGG通路可以运行到此处
add_descriptions.py -i KEGG_pathways_out/path_abun_unstrat.tsv.gz \
--custom_map_table picrust2/picrust2/default_files/description_mapfiles/KEGG_pathways_info.tsv.gz \
-o KEGG_pathways_out/path_abun_unstrat_descrip.tsv.gz

#运行整个pipeline(结果包含COG,EC,KO,PFAM,TIGRFAM的丰度信息和MetaCyc途径丰度)
picrust2_pipeline.py -s otu.fasta  -i otu_table.txt -o picrust2_out_pipeline -p 20 \
-r picrust2/picrust2/default_files/prokaryotic/pro_ref/pro_ref --in_traits COG,EC,KO,PFAM,TIGRFAM

#-r 是hmm数据库，默认是default_files/prokaryotic/pro_ref/pro_ref , 用于16S rRNA.
#添加功能描述
add_descriptions.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC \
-o EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz

add_descriptions.py -i KO_metagenome_out/pred_metagenome_unstrat.tsv.gz -m KO \
-o KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz

add_descriptions.py -i COG_metagenome_out/pred_metagenome_unstrat.tsv.gz -m COG \
-o COG_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz

add_descriptions.py -i PFAM_metagenome_out/pred_metagenome_unstrat.tsv.gz -m PFAM \
-o PFAM_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz

add_descriptions.py -i TIGRFAM_metagenome_out/pred_metagenome_unstrat.tsv.gz -m TIGRFAM \
-o TIGRFAM_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz

add_descriptions.py -i pathways_out/path_abun_unstrat.tsv.gz -m METACYC \
-o pathways_out/path_abun_unstrat_descrip.tsv.gz

##文件解释
#KO_metagenome_out/，该结果路径中记录了细菌群落KO（KEGG Orthology）功能的丰度预测结果。
#KO_metagenome_out/seqtab_norm.tsv.gz，对于很多细菌而言，一个个体可能包含多条16S（多拷贝16S），因此在原始OTU 16S rRNA丰度表的基础上，根据物种所含16S rRNA拷贝数对物种丰度进行标准化，得到校正16S rRNA拷贝数后的OTU丰度表。
#KO_metagenome_out/pred_metagenome_unstrat.tsv.gz，该文件中即为预测得到的细菌群落功能丰度表，记录了各样本中所包含KO功能的丰度，丰度计算由上述校正16S rRNA拷贝数标准化后的OTU丰度表推断得到。功能以KO ID为名称，代表了特定的功能基因。
#KO_metagenome_out/weighted_nsti.tsv.gz，各样本预测功能的加权NSTI值，由OTU的NSTI值通过标准化后的丰度加权所得。
#EC_metagenome_out/，该结果路径中记录了细菌群落酶（EC）功能的丰度预测结果。文件结构同上述KO_metagenome_out/，不再展示。
#pathways_out/path_abun_unstrat.tsv.gz，上述为预测得到的以KO ID为名称的KO功能，实则代表了特定的功能基因，将这些功能基因映射到具体的KEGG代谢途径（KEGG pathway）中，并统计各途径在各样本中的丰度，获得该表。
#KO_predicted.tsv.gz和EC_predicted.tsv.gz，两个矩阵文件中记录了OTU对预测功能丰度的贡献，即可以理解为每个OTU所代表的物种个体基因组中，分别有多少数量的基因与对应的KO功能或酶功能有关。如果期望关注哪些OTU是否对群落功能是重要的，这些表格（该表仅代表了单个物种个体基因组的特征，可能还需结合OTU的丰度信息）可以提供参考
#marker_predicted_and_nsti.tsv.gz，记录了OTU代表物种基因组中，16SrRNA拷贝数以及功能预测的NSTI值信息。
#Intermediate/，一些中间文件。


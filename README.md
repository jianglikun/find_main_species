# find_main_species

1、通过基因集与4904个细菌基因组比对结果（blastn）得到taxid

2、利用python ete3包（还需安装six包）通过taxid得到每个细菌对应的种属关系。

3、最后得到的是 cluster  gene  identity taxid  species 门，属，种    的表格

eg：
CAG1554	170627-01_S3_L003_scaffold13_1_gene16	100.000	1519439	Oscillibacter sp. ER4	Firmicutes,Oscillibacter,Oscillibacter sp. ER4

CAG1554	170627-01_S3_L003_scaffold13_1_gene16	93.103	40519	Ruminococcus callidus	Firmicutes,Ruminococcus,Ruminococcus callidus

CAG1554	170627-01_S3_L003_scaffold13_1_gene16	90.541	384638	[Bacteroides] pectinophilus	Firmicutes,[Bacteroides] pectinophilus,NA

CAG1554	170627-01_S3_L003_scaffold13_1_gene16	82.753	1121296	[Clostridium] aminophilum DSM 10710	Firmicutes,Lachnoclostridium,[Clostridium] aminophilum

CAG0050	170627-01_S3_L003_scaffold19_1_gene28	77.519	1229172	Leptolyngbya sp. KIOST-1	Cyanobacteria/Melainabacteria group,Leptolyngbya,Leptolyngbya sp. KIOST-1

CAG0017	170627-01_S3_L003_scaffold21_1_gene32	99.608	821	Bacteroides vulgatus	Bacteroidetes/Chlorobi group,Bacteroidaceae,Bacteroides

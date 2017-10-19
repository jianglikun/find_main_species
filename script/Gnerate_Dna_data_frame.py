#! /usr/bin/python
# -*- coding: UTF-8 -*-
import sys
import time
from ete3 import NCBITaxa

localtime = time.asctime( time.localtime(time.time()) )
print "本地时间为 :", localtime

dic_blastn = {}
dic_blastn_identity = {}
dic_cu={}
dic_genome_bacteria={}
dic_taxid_bacteria={}
dic_taxid_name = {}
with open(sys.argv[2]) as cu_fi:
                for line1 in cu_fi:
                        line1 = line1.strip()
                        blastn_gene_name1 = line1.split("\t")[1]
                        cu_name = line1.split("\t")[0]
			dic_cu[blastn_gene_name1]=cu_name

with open(sys.argv[3],'r') as tax_fi:
                for line2 in tax_fi:
                        line2 = line2.strip()
                        genome_name1 = line2.split("\t")[1]
                        bacteria_name = line2.split("\t")[0]
                    #    genome_name1 = genome_name1.split(".")[0]   
			dic_genome_bacteria[genome_name1]=bacteria_name

with open(sys.argv[4],'r') as tax_id_fi:
                for line3 in tax_id_fi:
                        line3 = line3.strip() 
                        tax_id = line3.split("\t")[2]
                        species_name = line3.split("\t")[3]
                        bacteria_name2 = line3.split("\t")[0]
			dic_taxid_bacteria[bacteria_name2]=tax_id
			dic_taxid_name[tax_id]=species_name

ncbi = NCBITaxa()

with open(sys.argv[1],'r') as blast_file:
                for line in blast_file:
                        line = line.strip()
                        blastn_gene_name =line.split("\t")[0]
			identity = line.split("\t")[2]
                        genome_name = line.split("\t")[1]
			dic_blastn_identity[blastn_gene_name]=identity
                        dic_blastn[blastn_gene_name]=genome_name
			if blastn_gene_name in dic_cu:
				lineage = ncbi.get_lineage(dic_taxid_bacteria[dic_genome_bacteria[dic_blastn[blastn_gene_name]]])
				names = ncbi.get_taxid_translator(lineage)
				tax_seq = [names[taxid] for taxid in lineage]
				if len(tax_seq)>9:
					tax_seq2 = [tax_seq[4],tax_seq[8],tax_seq[9]]
				
				elif len(tax_seq)==9:
					tax_seq2 = [tax_seq[4],tax_seq[8],"NA"]
				elif len(tax_seq)<=9:
					tax_seq2 = [tax_seq[4],"NA","NA"]
				#print tax_seq2
				tax_seq2 = ','.join(tax_seq2)
				tmp= dic_cu[blastn_gene_name]+"\t"+blastn_gene_name+"\t"+dic_blastn_identity[blastn_gene_name]+"\t"+dic_taxid_bacteria[dic_genome_bacteria[dic_blastn[blastn_gene_name]]]+"\t"+dic_taxid_name[dic_taxid_bacteria[dic_genome_bacteria[dic_blastn[blastn_gene_name]]]]+"\t"+tax_seq2+"\n"
				#print tmp
				fi=open(sys.argv[5],'a')
				fi.write(tmp)
				fi.close()

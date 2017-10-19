#! /bin/bash

if [ $# -eq 0 ];then
	echo "Usage: bash Gnerate_Dna_data_frame.sh -d(blastn_out) -i(cu_file) -a(assembly.id.txt) -g(4904MG.info.txt) -o(gene_cu_tax_data_frame)"
	echo "eg:bash Gnerate_Dna_data_frame.sh -d blastn_out -i ASD.40samples_clusters -a assembly.id.txt -g 4904MG.info.txt -o gene_cu_tax_frame"
	exit 1
fi

while getopts ":d:i:a:g:o:h" arg
do
        case $arg in
                d )
                        blastn_out=$OPTARG;;
		i )
			cu_file=$OPTARG;;
		a )
			assembly=$OPTARG;;
		g )
			MG=$OPTARG;;
                o )
                        dna_data=$OPTARG;;
                h )
                        echo "Usage: bash Gnerate_Dna_data_frame.sh -d(blastn_out) -i(cu_file) -a(assembly.id.txt) -g(4904MG.info.txt) -o(gene_cu_tax_data_frame)"
                        exit 0;;
                ? )
	                echo "Usage: bash Gnerate_Dna_data_frame.sh -d(blastn_out) -i(cu_file) -a(assembly.id.txt) -g(4904MG.info.txt) -o(gene_cu_tax_data_frame)"
			exit 1;;
        esac
done

less ${blastn_out}|awk '{if ($4>100) print $0}' > blastn_out
echo "filter mapping length > 100bp"
python Gnerate_Dna_data_frame.py  blastn_out ${cu_file} ${assembly} ${MG} ${dna_data}


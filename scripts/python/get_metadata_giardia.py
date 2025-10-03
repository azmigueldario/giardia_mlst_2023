import pandas as pd 
import numpy as np 

data=pd.read_csv('new_metadata.tsv', sep='\t')
data['location'] =  np.where(data['sample_alias'].str.contains('VANC|Goat river'), 'BC',
                    np.where(data['sample_alias'].str.contains('HAMILTON'), 'ONT',
                    np.where(data['sample_alias'].str.contains(), "NZ" 
                    'other')))

# awk to get accession list
awk 'BEGIN {FS=OFS="\t"} NR>1 {print $4}' metadata/new_metadata.tsv 

# entrez singularity location 
/mnt/cidgoh-object-storage/images/entrez-direct_2.7.70.sif

esearch -db sra -query SAMEA2018796 | efetch -format runinfo | head -n 1 > metadata_edirect.csv

for i in $(cat giardia_accs.txt)
do 
    esearch -db sra -query $i | efetch -format runinfo | tail -1 >> metadata_edirect.csv
    echo "Done with sample: ${i}!"
done


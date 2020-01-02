# Unique Sequence counting for amplicon based sequencing projects. 
Pairs sequencing reads together then counts unique items. 
Outputs csv file per fastq pair. 

## To Run:
nextflow run alemenze/magicuniqueamplicons  --reads '/path_to_reads/*_R{1,2}*.fastq.gz' 

## Requires: 
Docker==19.03.5
Nextflow==19.10.0

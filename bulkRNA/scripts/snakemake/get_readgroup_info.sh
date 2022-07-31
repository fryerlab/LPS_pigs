#!/bin/bash

# change directory
cd ../../rawfastq

# create file with list of R1 samples
ls -1 | grep _R1_ > R1Samples.txt

# change directory 

# loops through list 
touch sampleReadInfo.txt
for sample in `cat R1Samples.txt`; do
    zcat ${sample} | head -1 >> sampleReadInfo.txt	
done;

# mv the files 
mv R1Samples.txt  /LPS_pigs/bulkRNA/scripts/snakemake/R1Samples.txt
mv sampleReadInfo.txt /LPS_pigs/bulkRNA/scripts/snakemake/sampleReadInfo.txt

cd /LPS_pigs/bulkRNA/scripts/snakemake/
paste -d "\t" R1Samples.txt sampleReadInfo.txt > sampleReadGroupInfo.txt
rm R1Samples.txt
rm sampleReadInfo.txt


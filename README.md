# LPS_pigs
Bulk RNAseq of Sus scrofa (pigs) that received either saline or bacterial lipopolysaccharide (LPS).

The goal of this experiment is to identify differentially expressed genes (DEGs) between experimental groups.  Pigs were injected with saline (control) or lipopolysaccharide (LPS) to model sepsis.  Brain, kidney and blood samples were collected and sent for bulk RNA sequencing.

## Explore differentially expressed genes and isoforms and correlation among tissues in our published shiny apps

isoform-level differential analysis for:\
Brain -> https://fryerlab.shinyapps.io/LPS_pigs_brain_isoform/ \
Kidney -> https://fryerlab.shinyapps.io/LPS_pigs_kidney_isoform/ \ 
Blood -> https://fryerlab.shinyapps.io/LPS_pigs_blood_isoform/ \


## Set up conda environment
This workflow uses conda. For information on how to install conda, visit: https://docs.conda.io/projects/conda/en/latest/user-guide/index.html 

To create the pigs environment:

`conda env create -f pigs.yml`

To activate the environment once installed:
`conda activate pigs`

## bulk RNAseq differential expression
We have put together a workflow for inferring differential expression between saline (control) and LPS female pigs using two read aligners STAR and SALMON, and limma/voom for computing differential expression. These tools are publicly available and we ask that if you use this workflow to cite the tools used listed in the table below. 

### 1. set up project directory structure
`cd bulkRNA`\
`mkdir bamstats  featureCounts  kallisto  rObjects  rawQC  rawfastq  results  starAligned  trimmedQC  trimmedReads`\
`cd results`\
`mkdir star kallisto`

### 2. align reads and generate quantification estimates
We ran the data through two pipelines as it has been shown that most of the variation in our data can be explained by the choice of read aligner, see Olney et al. 2020. Therefore, we aligned the reads with star followed by RSEM for quantification and employed the pseudo aligner Kallisto. 
In the snakemake pipeline, the fastq files are trimmed for quality via bbduk package. Information on the specific parameters employed for each step is explained within the snakefile. 

First move to the scripts snakemake folder
`cd LPS_pigs/bulkRNA/scripts/snakemake/`\
Now run the snakefile. 
`snakemake -s Snakefile`\
To run multiple samples in parallel use -j and the number of jobs to run.

### 3. preform differntial expression

move to the results folder. Create sub directories.
 
`cd star`\
`mkidr CPM  DEGs  JSD  MDS  boxplot  density  gprofiler  metascape  volcano  voom`

move into the kallisto results folder.

`cd ../kallisto`\
`mkdir CPM  DEGs  JSD  MDS  boxplot  density  gprofiler  metascape  volcano  voom TPM`

move to the scripts R folder.\
`cd ../../scripts/R`\
read over the gene DE.Rmd file to deteremine which tissue analysis you want to perform (Brain, Kidney, or Blood).

Open up R and run the R script for gene level differential expression. 

`R gene_DE.R`


### publicly available tools used in this analysis


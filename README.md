# LPS_pigs
Bulk RNAseq of Sus scrofa (pigs) that received either saline or bacterial lipopolysaccharide (LPS).

The goal of this experiment is to identify differentially expressed genes (DEGs) between experimental groups.  Pigs were injected with saline (control) or lipopolysaccharide (LPS) to model sepsis.  Brain, kidney and blood samples were collected and sent for bulk RNA sequencing.

## Set up conda environment
This workflow uses conda. For information on how to install conda, visit: https://docs.conda.io/projects/conda/en/latest/user-guide/index.html 

To create the pigs environment:

`conda env create -f pigs.yml`

To activate the environment once installed:
`conda activate pigs`

## bulk RNAseq differential expression
We have put together a workflow for inferring differential expression between saline (control) and LPS female pigs using two read aligners STAR and SALMON, and limma/voom for computing differential expression. These tools are publicly available and we ask that if you use this workflow to cite the tools used listed in the table below. 

To run our suggested workflow for differential expression please see bulkRNA folder. 

### publicly available tools used in this analysis
Tool | usage | citation
--- | --- |  ---
Trimmomatic | Trim RNA-sequences for quality | Bolger AM, Lohse M, Usadel B. Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics. 2014;30: 2114–2120.
SALMON | RNAseq read aligner | <salmon citation>
STAR | RNAseq read aligner | Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, et al. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013;29: 15–21
bamtools | analyzing and processing BAM files | Barnett DW, Garrison EK, Quinlan AR, Strömberg MP, Marth GT. BamTools: a C++ API and toolkit for analyzing and managing BAM files. Bioinformatics. 2011;27: 1691–1692.
FeatureCounts | obtain raw transcriptome counts| Liao Y, Smyth GK, Shi W. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics. 2014;30: 923–930.
Limma/voom | differenital expression analysis | Law CW, Chen Y, Shi W, Smyth GK. voom: Precision weights unlock linear model analysis tools for RNA-seq read counts. Genome Biol. 2014;15: R29.


# Pig reference genome 
Below are instructions for downloading and formatting the pig reference genome and transcriptome for downstream use. 

## 1. Download fasta and GTF from Ensembl

`wget http://ftp.ensembl.org/pub/release-107/fasta/sus_scrofa/dna/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz`

`wget http://ftp.ensembl.org/pub/release-107/gtf/sus_scrofa/Sus_scrofa.Sscrofa11.1.107.gtf.gz`

Unzip the files

`gunzip Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz`

`gunzip Sus_scrofa.Sscrofa11.1.107.gtf.gz`

## 2. Hard mask the Y chromosome
All pigs in this study are XX female. We will hard mask the Y chromosome sequences to improve expression estimates for X-linked genes. See Olney et al. 2020 https://doi.org/10.1186/s13293-020-00312-9 for more information. 

`python hardmaskY.py`

## 3. Build index and dictionary 

Remember that you will need to activate the conda environment first. 

`conda activate pigs`

Make a directory for storing the star reference genome outputs

`mkdir Sus_scrofa.Sscrofa11.1.dna.toplevel_star_Ymask`

`star --runThreadN 12 --runMode genomeGenerate --genomeDir Sus_scrofa.Sscrofa11.1.dna.toplevel_star_Ymask --genomeFastaFiles Sus_scrofa.Sscrofa11.1.dna.toplevel.Ymask.fa --sjdbGTFfile Sus_scrofa.Sscrofa11.1.107.gtf`

KEY:\
--sjdbOverhang ReadLength-1, default value of 100\
--runThreadN # of threads\
--runMode build genome indices (default is aligns reads)\
--genomeFastaFiles /path/to/genome/fasta\
--sjdbGTFfile /path/to/annotations.gtf\
--genomeDir /path/to/genomeDir specifies where the genome indices are stored. Note needs to be made prior to this step\ 






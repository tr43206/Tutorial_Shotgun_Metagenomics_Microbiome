# Microbiome_Tutorial_Shotgun-Metagenomics
Output fastq files of Gut-Lung Axis project from MGI DNBSEQ-G99.



## Step 0 : Preprocessing

1) Quality control (QC) using `fastp` or `faqcs`.
2) Host genome (human or mouse) removal using `bowtie2`.



## Step 1 : Read-based analyses

1) Taxonomic classification using `kraken2` or `metaphlan4` (`metaphlan4` is stricter than `kraken2`).
2) Functional annotation using `humann3`. Could identify contributing species in specific pathways by stratifying.
3) For kraken2 results, relative abundance calculation with bracken is needed.
4) Differential analysis (DA) methods (e.g., LEfSe, ALDEx2, ANCOM-BC2, MaAsLin3). Adjusting `sequence reads` as covariates is recommended.



## Step 2 : Assembly

1) Assemble reads to contigs using `megahit` or `metaspades`.
2) Compare contig qualities of two assembly tools using metaquast.



## Step 3 : Binning and generating MAGs

1) Contigs to bins using `metabat2` or `concoct` or `maxbin2`.
2) Quality check using `checkm2`, and classify quality levels to `High`, `Medium`, and `Low`. Only `High` and `Medium` level MAGs will be used in following steps.



## Step 4 : MAG-based analyses

1) Taxonomic classification of individual MAGs using `gtdbtk2`. Verification of classified results is needed (by comparing ANI identity, etc.). `drep` can be used to identify functional abilities of specific MAGs.
2) Funtional gene annotation using `eggnog-mapper`. Gene classes can be classified in A to Z classes. Results of `eggnog-mapper` can be also coverted to KEGG or EC number level.
3) Carbohydrate classes or candidate carbohydrate substrates can be annotated by `dbcan3`.
4) Gene classes of specific MAGs can be also identified using `bakta`.
5) `Abundance`, `coverage`, and `beadth` calculation of specific MAGs in each samples using `coverm`.

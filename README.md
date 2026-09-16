# Tutorial_Shotgun_Metagenomics_Microbiome
Output fastq files of Gut-Lung Axis project from MGI DNBSEQ-G99.



## Step 0 : Preprocessing

1) Quality control (QC) using `fastp` or `faqcs`.
   ```bash
   fastp \
   -i C01_1.fq.gz \
   -I C01_2.fq.gz \
   -o ../01.fastp.out_q30/C01_1_fastp.fastq.gz \
   -O ../01.fastp.out_q30/C01_2_fastp.fastq.gz \
   --html ../01.fastp.out_q30/C01.html \
   --json ../01.fastp.out_q30/C01.json \
   --thread \
   -q 30
   ```
   ```bash
   FaQCs \
   -1 C01_1.fq.gz \
   -2 C01_2.fq.gz \
   --prefix C01 \
   -d ../02.faqcs_q30/ \
   -t 32 \
   -q 30
   ```
3) Host genome (human or mouse) removal using `bowtie2`.
   ```bash
   bowtie2 \
   --threads 12 \
   -x ../Database/C57BL6J/GRCm39_DB/GRCm39_DB \
   -1 C01_1_fastp.fastq.gz \
   -2 C01_2_fastp.fastq.gz \
   --very-sensitive-local \
   --un-conc-gz ../03.host_removal/C01_host_removed > ../03.host_removal/C01.sam
   ```
   ```bash
   mv C01_host_removed.1 C01_host_removed_R1.fastq.gz && \
   mv C01_host_removed.2 C01_host_removed_R2.fastq.gz
   ```



## Step 1 : Read-based analyses

1) Taxonomic classification using `kraken2` or `metaphlan4` (`metaphlan4` is stricter than `kraken2`).
   ```bash
   kraken2 \
   --db ~/Database/kraken2DB_v260226/ \
   --report ../04.kraken2.out/C01.report \
   --output ../04.kraken2.out/C01.out \
   --paired C01_host_removed_R1.fastq.gz C01_host_removed_R2.fastq.gz \
   --threads 8 \
   --use-names \
   --gzip-compressed
   ```
   ```bash
   metaphlan \
   C01_host_removed_R1.fastq.gz,C01_host_removed_R2.fastq.gz \
   --input_type fastq \
   --nproc 16 \
   --db_dir ~/Database/metaphlan_db_v4.2.4/ \
   --mapout ../04.metaphlan4.out/C01.mapout.bz2 \
   -o ../04.metaphlan4.out/C01_profile.txt
   ```
   ```bash
   merge_metaphlan_tables.py \
   *_profile.txt > merged_abundance_table.txt
   ```
3) Functional annotation using `humann3`. Could identify contributing species in specific pathways by stratifying.
   ```bash
   zcat C01_host_removed_R1.fastq.gz C01_host_removed_R2.fastq.gz > zcat/C01_all_reads.fastq \
   --verbose
   ```
   ```bash
   samples="C01 C02 C03 C09 C10 C11 V01 V02 V03 V04 V05 V06 V07 V08 V10" \
   for s in $samples \
   do \
   humann \
   -i ${s}_all_reads.fastq \
   -o ../../04.humann.out/${s}_humann.out \
   --threads 12 \
   done
   ```
   ```bash
   samples="C01 C02 C03 C09 C10 C11 V01 V02 V03 V04 V05 V06 V07 V08 V10" \
   for s in $samples \
   do \
   humann_renorm_table \
   --input ${s}_humann.out/${s}_all_reads_genefamilies.tsv \
   --output ${s}_humann.out/${s}_all_reads_genefamilies_relab.tsv \
   --units relab \
   -s n \
   done
   ```
   ```bash
   humann_join_tables \
   --input . \
   --output genefamilies_relab.tsv \
   --file_name genefamilies_relab \
   -s
   ```
   ```bash
   humann_split_stratified_table \
   --input genefamilies_relab.tsv \
   --output genefamilies_relab_split
   ```
   ```bash
   humann_regroup_table \
   --input genefamilies_relab.tsv \
   --groups uniref90_ko \
   --output genefamilies_relab_ko.tsv \
   -u N
   ```
   ```bash
   humann_regroup_table \
   --input genefamilies_relab.tsv \
   --groups uniref90_level4ec \
   --output genefamilies_relab_ec.tsv \
   -u N
   ```
5) For kraken2 results, relative abundance calculation with bracken is needed.
   ```bash
   mkdir -p bracken/species bracken/genus
   ```
   ```bash
   samples="C01 C02 C03 C09 C10 C11 V01 V02 V03 V04 V05 V06 V07 V08 V10" \
   for s in $samples \
   do \
   bracken \
   -d ~/Database/kraken2DB_v260226 \
   -i ${s}.report \
   -o bracken/species/${s}.species.bracken \
   -r 150 \
   -l S \
   -t 32 \
   done
   ```
   ```bash
   python /home/sejong/anaconda3/envs/kraken2/bin/combine_bracken_outputs.py \
   --files species/*.species.bracken \
   --names C01,C02,C03,C09,C10,C11,V01,V02,V03,V04,V05,V06,V07,V08,V10 \
   -o species_merged.tsv
   ```
6) Differential analysis (DA) methods (e.g., LEfSe, ALDEx2, ANCOM-BC2, MaAsLin3). Adjusting `sequence reads` as covariates is recommended.



## Step 2 : Assembly

1) Assemble reads to contigs using `megahit` or `metaspades`.
   ```bash
   megahit \
   -t 32 \
   -1 C01_host_removed_R1.fastq.gz \
   -2 C01_host_removed_R2.fastq.gz \
   -o ../04.megahit.out/C01
   ```
   ```bash
   spades.py \
   --meta \
   -t 4 \
   -1 C01_host_removed_R1.fastq.gz \
   -2 C01_host_removed_R2.fastq.gz \
   -o ../04.spades.out/C01
   ```
3) Compare contig qualities of two assembly tools using `quast` or `metaquast`.
   ```bash
   quast.py \
   -t 32 \
   -o ../05.quast.out/C01 \
   -1 C01_host_removed_R1.fastq.gz \
   -2 C01_host_removed_R2.fastq.gz \
   ../04.spades.out/C01/scaffolds.fasta
   ```



## Step 3 : Binning and generating MAGs

1) Contigs to bins using `metabat2` or `concoct` or `maxbin2`.
   ```bash
   jgi_summarize_bam_contig_depths \
   --outputDepth metabat2/C01.depth.txt \
   C01.sorted.bam
   ```
   ```bash
   metabat2 \
   -t 32 \
   -i ../../04.spades.out/C01/scaffolds.fasta \
   -a C01.depth.txt \
   -o C01_metabat2.bins
   ```
   ```bash
   cut_up_fasta.py \
   ../../04.spades.out/C01/scaffolds.fasta \
   -c 10000 \
   -o 0 --merge_last \
   -b C01_contigs_10K.bed > C01_contigs_10K.fa
   ```
   ```bash
   concoct_coverage_table.py \
   C01_contigs_10K.bed \
   ../C01.sorted.bam > C01_coverage_table.tsv
   ```
   ```bash
   samples="C01 C02 C03 C09 C10 C11 V01 V02 V03 V04 V05 V06 V07 V08 V10" \
   for s in $samples \
   do \
   concoct \
   --composition_file ${s}_contigs_10K.fa \
   --coverage_file ${s}_coverage_table.tsv \
   -b ${s}_concoct.bins \
   -t 8 \
   done
   ```
   ```bash
   samples="C01 C02 C03 C09 C10 C11 V01 V02 V03 V04 V05 V06 V07 V08 V10" \
   for s in $samples \
   do \
   merge_cutup_clustering.py \
   ${s}_concoct.bins_clustering_gt1000.csv > ${s}_concoct.bins_clustering_merged.csv \
   done
   ```
   ```bash
   samples="C01 C02 C03 C09 C10 C11 V01 V02 V03 V04 V05 V06 V07 V08 V10" \
   for s in $samples \
   do \
   mkdir ${s}_fasta_bins \
   done
   ```
   ```bash
   samples="C01 C02 C03 C09 C10 C11 V01 V02 V03 V04 V05 V06 V07 V08 V10" \
   for s in $samples \
   do \
   extract_fasta_bins.py \
   --output_path ${s}_fasta_bins/ \
   ../../04.spades.out/${s}/scaffolds.fasta \
   ${s}_concoct.bins_clustering_merged.csv \
   done
   ```
   ```bash
   samples="C01 C02 C03 C09 C10 C11 V01 V02 V03 V04 V05 V06 V07 V08 V10" \
   for s in $samples \
   do \
   jgi_summarize_bam_contig_depths \
   --outputDepth ${s}.depth.txt \
   --noIntraDepthVariance ../${s}.sorted.bam \
   done
   ```
   ```bash
   samples="C01 C02 C03 C09 C10 C11 V01 V02 V03 V04 V05 V06 V07 V08 V10" \
   for s in $samples \
   do \
   sample=$(head -n 1 ${s}.depth.txt | cut -f 4) \
   done
   ```
   ```bash
   samples="C01 C02 C03 C09 C10 C11 V01 V02 V03 V04 V05 V06 V07 V08 V10" \
   for s in $samples \
   do \
   grep -v totalAvgDepth ${s}.depth.txt | cut -f 1,4 > ${s%%.*}.txt \
   done
   ```
   ```bash
   samples="C01 C02 C03 C09 C10 C11 V01 V02 V03 V04 V05 V06 V07 V08 V10" \
   for s in $samples \
   do \
   run_MaxBin.pl \
   -contig ../../04.spades.out/${s}/scaffolds.fasta \
   -abund ${s}.txt \
   -out bins/${s} \
   -thread 32 \
   done
   ```
3) Quality check using `checkm2`, and classify quality levels to `High`, `Medium`, and `Low`. Only `High` and `Medium` level MAGs will be used in following steps.
   ```bash
   samples="C01 C02 C03 C09 C10 C11 V01 V02 V03 V04 V05 V06 V07 V08 V10" \
   for s in $samples \
   do \
   checkm2 predict \
   -t 32 \
   -x fa \
   --input 07.dastool/${s}/${s}_DASTool_bins/ \
   --output-directory 08.checkm/${s}/ \
   done
   ```



## Step 4 : MAG-based analyses

1) Taxonomic classification of individual MAGs using `gtdbtk2`. Verification of classified results is needed (by comparing ANI identity, etc.). `drep` can be used to identify functional abilities of specific MAGs.
   ```bash
   for d in C* V* \
   do \
   sample="$d" \
   bin_dir="$d/${d}_DASTool_bins" \
   for f in "$bin_dir"/concoct* \
   do \
   [ -e "$f" ] || continue \
   base=$(basename "$f") \
   mv "$f" "$bin_dir/${sample}_${base}" \
   done \
   done
   ```
   ```bash
   dRep \
   dereplicate \
   12.drep/drep_all \
   -g 07.dastool/*/*_DASTool_bins/*.fa \
   --genomeInfo 12.drep/quality_report_merged_drep.csv \
   -p 32 \
   -sa 0.95 \
   --S_algorithm ANImf
   ```
   ```bash
   dRep \
   dereplicate \
   drep_Pvulgatus/ \
   -g drep_Pvulgatus/bins/*.fa \
   --genomeInfo quality_report_merged_drep_Pvulgatus.csv \
   -p 32 \
   -sa 0.99 \
   --S_algorithm ANImf
   ```
3) Funtional gene annotation using `eggnog-mapper`. Gene classes can be classified in A to Z classes. Results of `eggnog-mapper` can be also coverted to KEGG or EC number level.
   ```bash
   emapper.py \
   -i 01.bakta/C02_concoct.31/C02_concoct.31.faa \
   --itype proteins \
   -m diamond \
   --cpu 32 \
   --data_dir ~/Database/eggnog-mapper-data/ \
   --output C02_concoct.31 \
   --output_dir 04.eggnog/
   ```
4) Carbohydrate classes or candidate carbohydrate substrates can be annotated by `dbcan3`.
   ```bash
   run_dbcan \
   dereplicated_genomes/C02_concoct.31.fa prok \
   --out_dir dbcan3/C02_concoct.31 \
   --db_dir ../../Database/dbCAN3/db \
   --dia_cpu 32 \
   --hmm_cpu 32 \
   --dbcan_thread 32 \
   --tools all
   ```
   ```bash
   run_dbcan \
   bakta/C02_concoct.31/C02_concoct.31.faa \
   protein \
   --out_dir dbcan3/C02_concoct.31_CGC \
   --db_dir ../../Database/dbCAN3/db \
   --tools hmmer \
   --tf_cpu 32 \
   --stp_cpu 32 \
   --dia_cpu 32 \
   --hmm_cpu 32 \
   -c bakta/C02_concoct.31/C02_concoct.31.gff3
   ```
   ```bash
   run_dbcan \
   bakta/C02_concoct.31/C02_concoct.31.faa \
   protein \
   --out_dir dbcan3/C02_concoct.31_substrate \
   --db_dir ../../Database/dbCAN3/db \
   -c bakta/C02_concoct.31/C02_concoct.31.gff3 \
   --tools hmmer \
   --cgc_substrate \
   --dia_cpu 32 \
   --hmm_cpu 32 \
   --dbcan_thread 32
   ```
   ```bash
   dbcan_plot \
   CGC_plot \
   -i C02_concoct.31_CGC/ \
   --cgcid 'contig_1|CGC3'
   ```
6) Gene classes of specific MAGs can be also identified using `bakta`.
   ```bash
   bakta \
   --db ../../Database/db/ \
   --prefix C02_concoct.31 \
   --output bakta/C02_concoct.31 \
   --threads 32 \
   dereplicated_genomes/C02_concoct.31.fa
   ```
   ```bash
   bakta_plot \
   --output plot_cog_600dpi \
   --prefix C02_concoct.31 \
   --sequences 2 \
   --type cog \
   --label "Phocaeicola vulgatus|representative MAG" \
   --size 8 \
   --dpi 600 \
   C02_concoct.31.json
   ```
7) `Abundance`, `coverage`, and `beadth` calculation of specific MAGs in each samples using `coverm`.
   ```bash
   Control="C01 C02 C03 C09 C10 C11" \
   VNAM="V01 V02 V03 V04 V05 V06 V07 V08 V10" \
   samples="$Control $VNAM" \
   for s in $samples \
   do \
   R1=../../03.host_removal/${s}_host_removed_R1.fastq.gz \
   R2=../../03.host_removal/${s}_host_removed_R2.fastq.gz \
   coverm genome \
   --coupled "$R1" "$R2" \
   --genome-fasta-files dereplicated_genomes/C02_concoct.31.fa \
   --methods relative_abundance trimmed_mean covered_bases length count \
   --min-covered-fraction 0 \
   --min-read-percent-identity-pair 98 \
   --min-read-aligned-percent-pair 90 \
   --proper-pairs-only \
   --exclude-supplementary \
   --output-format sparse \
   -t 32 \
   -o 03.coverm/${s}.tsv \
   done
   ```

##quality control
fastqc -t 10 B73-Chip-RNA-library_L5_A012.R1.clean.fastq.gz B73-Chip-RNA-library_L5_A012.R2.clean.fastq.gz B73-input_L5_A010.R1.clean.fastq.gz B73-input_L5_A010.R2.clean.fastq.gz ## no adapter with low quaility reads
java -jar trimmomatic-0.36.jar PE -phred33 B73-Chip-RNA-library_L5_A012.R1.clean.fastq.gz B73-Chip-RNA-library_L5_A012.R2.clean.fastq.gz B73-RAC_clean1.fq.gz un_RAC1.fq.gz B73-RAC_clean2.fq.gz un_RAC2.fq.gz LEADING:20 TRAILING:20 MINLEN:36
java -jar trimmomatic-0.36.jar PE -phred33 B73-input_L5_A010.R1.clean.fastq.gz B73-input_L5_A010.R2.clean.fastq.gz B73-input_clean1.fq.gz un_input1.fq.gz B73-input_clean2.fq.gz un_input2.fq.gz LEADING:20 TRAILING:20 MINLEN:36

##mapping with hisat2
hisat2 -q -p 12 -x Zea_mays.AGPv3.31 -1 B73-RAC_clean1.fq.gz -2 B73-RAC_clean2.fq.gz -S B73-RAC_clean.sam --known-splicesite-infile genome.ss -t
hisat2 -q -p 12 -x Zea_mays.AGPv3.31 -1 B73-input_clean1.fq.gz -2 B73-input_clean2.fq.gz -S B73-input_clean.sam --known-splicesite-infile genome.ss -t

samtools sort -@ 8 -o B73-RAC_clean_sort.bam B73-RAC_clean.sam
samtools sort -@ 8 -o B73-input_clean_sort.bam B73-input_clean.sam

bamToBed -i B73-RAC_clean_sort.bam | sort -k 1,1 -k 2,2n > B73-RAC_clean_sort.bed
bamToBed -i B73-input_clean_sort.bam | sort -k 1,1 -k 2,2n > B73-input_clean_sort.bed

bedtools genomecov -bg -scale 0.0538776789194207 -i B73-RAC_clean_sort.bed -g Zea_mays.AGPv3.31.fa.fai > B73-Chip-RNA-library_L5.bedGraph ## accroding to the total reads
bedtools genomecov -bg -scale 0.1014564273052524 -i B73-input_clean_sort.bed -g Zea_mays.AGPv3.31.fa.fai > B73-input_L5.bedGraph ## accroding to the total reads

##the bedgrapg files used for visualization in local IGV for Fig 1D, S1E and S1F Fig.
less B73-Chip-RNA-library_L5.bedGraph | awk '{if($1=="2" || $1="5") {print $0}}' | less > B73-Chip-RNA-library_L5_Chr2_Chr5.bedGraph
less B73-input_L5.bedGraph | awk '{if($1=="2" || $1=="5") {print $0}}' | less > B73-input_L5_Chr2_Chr5.bedGraph


## merge reads
SeqPrep -f B73-RAC_clean1.fq -r B73-RAC_clean2.fq -1 B73-RAC_clean_trim1.fq.gz -2 B73-RAC_clean_trim2.fq.gz -q 30 -L 25 -s B73-RAC_clean_trim_merge.fq.gz -E B73-RAC_trim_alignment.txt
SeqPrep -f B73-2-input_clean1.fq -r B73-2-input_clean2.fq -1 B73-2-input_clean_trim1.fq.gz -2 B73-2-input_clean_trim2.fq.gz -q 30 -L 25 -s B73-2-input_clean_trim_merge.fq.gz -E B73-2-input_trim_alignment.txt


## file format from fastq to fasta
gzip -d B73-RAC_clean_trim_merge.fq.gz && fq2fa B73-RAC_clean_trim_merge.fq > B73-RAC_clean_trim_merge.fa
gzip -d B73-2-input_clean_trim_merge.fq.gz && fq2fa B73-2-input_clean_trim_merge.fq > B73-2-input_clean_trim_merge.fa

## the files sequence.fa contain 356bp, CRM1, CRM2, CRM3, CRM4, CentC, all genes DNA sequence.

blastn -query sequence.fa -db B73-RAC_clean_trim_merge.fa -outfmt 6 -evalue 1e-5 -out B73-RAC_clean_trim_merge_maize.blastn
blastn -query sequence.fa -db B73-2-input_clean_trim_merge.fa -outfmt 6 -evalue 1e-5 -out B73-2-input_clean_trim_merge_maize.blastn

## we compare the reletive enrichment between RIP-seq and Input-seq files.

## download illumina reads from GEO SRX1452310, SRR3018834, SRR2000635, SRR2000640, SRR2000648 and GSE124242.
## SeqPrep with these reads
## download Pacbio reads SRX1472849.

## use the blast results to 356bp can idetified the reads containing the back-spliced site in Fig 1E.

blastn -query 356bp.fa -db SRX1452310_merge.fa -outfmt 6 -evalue 1e-5 -out SRX1452310_merge_maize_356.blastn
blastn -query 356bp.fa -db SRR3018834_merge.fa -outfmt 6 -evalue 1e-5 -out SRX1452310_merge_maize_356.blastn
blastn -query 356bp.fa -db SRR2000635_merge.fa -outfmt 6 -evalue 1e-5 -out SRR2000635_merge_maize_356.blastn
blastn -query 356bp.fa -db SRR2000640_merge.fa -outfmt 6 -evalue 1e-5 -out SRR2000640_merge_maize_356.blastn
blastn -query 356bp.fa -db SRR2000648_merge.fa -outfmt 6 -evalue 1e-5 -out SRR2000648_merge_maize_356.blastn
blastn -query 356bp.fa -db GSE124242_merge.fa -outfmt 6 -evalue 1e-5 -out GSE124242_merge_maize_356.blastn
blastn -query 356bp.fa -db SRX1472849_merge.fa -outfmt 6 -evalue 1e-5 -out SRX1472849_merge_maize_356.blastn

## wheat 323 nt circRNA
SeqPrep -f TMU38.R1.clean.fastq.gz -r TMU38.R2.clean.fastq.gz -1 TMU38.R1.clean_trim1.fq.gz -2 TMU38.R2.clean_trim2.fq.gz -q 30 -L 25 -s TMU38.clean_trim_merge.fq.gz -E TMU38.clean.trim_alignment.txt
gzip -d TMU38.clean_trim_merge.fq.gz && fq2fa TMU38.clean_trim_merge.fq > TMU38.clean_trim_merge.fa
## blast the merge fasta file to wheat retrotransposon file

## check the potential back-spliced RNA
blastn -query retrotransposon.fa -db SRX1472849_merge.fa -outfmt 6 -evalue 1e-5 -out SRX1472849_merge_maize_356.blastn

## 323-nt circ-RNA mapping to the chinese spring genome
blastn -query 323-nt.fa -db 161010_Chinese_Spring_v1.0_pseudomolecules  -outfmt 6 -evalue 1e-5 -out 323-nt_161010_Chinese_Spring_v1.0_pseudomolecules.blastn 
format_blastn6_bed.pl -i 323-nt_161010_Chinese_Spring_v1.0_pseudomolecules.blastn -o wheat_323.CS.bed

##the bed file used for visualization in local IGV for S5I Fig.

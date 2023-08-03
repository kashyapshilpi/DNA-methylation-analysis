#fastqc of row data 
fastqc*/*.fastqc.gz
#multiqc 
multiqc.
#trimming of the Row data 
java -jar /home/abcd/shilpi/PRJNA289892/DNA/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 16 /home/abcd//shilpi/PRJNA289892/DNA/$1/$1\_1.fastq.gz /home/abcd//shilpi/PRJNA289892/DNA/$1/$1\_2.fastq.gz /home/abcd//shilpi/PRJNA289892/DNA/$1/output_forward_paired.fq.gz /home/abcd//shilpi/PRJNA289892/DNA/$1/output_forward_unpaired.fq.gz /home/abcd//shilpi/PRJNA289892/DNA/$1/output_reverse_paired.fq.gz /home/abcd//shilpi/PRJNA289892/DNA/$1/output_reverse_unpaired.fq.gz ILLUMINACLIP:/home/abcd//shilpi/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
ls -d SRR21023*|sort|while read i;do sh trim.sh $i;done
#fastqc of the trim data 
fastqc*.fq.gz
#multiqc of the new fastqc files
multiqc. #combine all the html file 
# Running bismark_genome_preparation
time bismark_genome_preparation --path_to_aligner /home/scbb/shilpi/PRJNA289892/DNA/bowtie2-2.5.1-linux-x86_64 --bowtie2 /home/scbb/shilpi/genome/DNA --parallel 6
#running bismark
bismark --parallel 6 /home/scbb/shilpi/genome/DNA/ -1 SRR2102340_output_forward_paired.fq -2 SRR2102340_output_reverse_paired.fq -o sht_1
#Running deduplicate_bismark 
/media/scbb/data3/shilpi_trainee/PRJNA289892/DNA/Bismark-0.24.0/deduplicate_bismark SRR2102339_output_forward_paired_bismark_bt2_pe.bam 
#homeRunning bismark_methylation_extractor
/home/scbb/shilpi/PRJNA289892/DNA/Bismark-0.24.0/bismark_methylation_extractor --gzip --bedGraph -p SRR2102340_output_forward_paired_bismark_bt2_pe.deduplicated.bam 
#Running bismark2report
/home/scbb/shilpi/PRJNA289892/DNA/Bismark-0.24.0/bismark2report SRR2102340_output_forward_paired_bismark_bt2_pe.deduplicated.M-bias.txt
#Running bismark2summary
bismark2summary #This command scans the current working directory for different Bismark alignment, deduplication and methylation extraction (splitting) reports to produce a graphical summary HTML report, as well as a data table

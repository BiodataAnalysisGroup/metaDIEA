
# Bowtie2 aligner
# http://bowtie-bio.sourceforge.net/index.shtml

# tophat2
# https://ccb.jhu.edu/software/tophat/index.shtml

# Cufflinks assembles transcripts
# https://github.com/cole-trapnell-lab/cufflinks

# fasta or fa --> genome file
# fastq --> raw read files to be aligned on the genome file 

cd /app

# 0| Build an index from the FASTA genome file using bowtie

bowtie-build genome/genome.mm10.fa genome.mm10 -p 2

# 1| Map the reads for each sample to the reference genome:

for i in `seq -w 1 10|head -n 6`
do
echo $i
tophat2 -p 16 -G genome/genes.gtf -o sample_${i} genome.mm10 raw/sample_${i}_1.fastq.gz   raw/sample_${i}_2.fastq.gz --no-novel-juncs &> sample_${i}.log.txt 
done


# 2| Assemble transcripts for each sample:

for i in `seq -w 1 10|head -n 6`
do
echo $i
cufflinks -p 16 -o sample_${i}_clout sample_${i}/accepted_hits.bam
done


# 3| Create a file called assemblies.txt that lists the assembly file for each sample. The file should contain the following lines:

for i in `seq -w 1 10|head -n 6`
do
echo "./sample_${i}_clout/transcripts.gtf" 
done > assemblies.txt

# 4| Run Cuffmerge on all your assemblies to create a single merged transcriptome annotation:
cuffmerge -g genome/genes.gtf -s genome/genome.mm10.fa -p 16 assemblies.txt

# 5| Run Cuffdiff by using the merged transcriptome assembly along with the BAM files from TopHat for each replicate:
cuffdiff -o diff_out -b genome/genome.mm10.fa -p 16 -L C1,C2 -u merged_asm/merged.gtf ./sample_01/accepted_hits.bam,./sample_02/accepted_hits.bam,./sample_03/accepted_hits.bam ./sample_04/accepted_hits.bam,./sample_05/accepted_hits.bam,./sample_06/accepted_hits.bam &> cuffdiff.log.txt



#RSEM is a software package for estimating gene and isoform expression levels from RNA-Seq data
#https://github.com/deweylab/RSEM

cd /app

#create index files
mkdir -p genome/RSEM.mm10
genome=./genome
genome_fa=./genome/genome.fa
genes=./genome/genes.gtf
rsem-prepare-reference --gtf $genes --bowtie2 $genome_fa ./genome/RSEM.mm10/mm10

# Sample Analysis
transcripts="./genome/RSEM.mm10/mm10"

# 1| Map the reads for each sample to the reference genome:
for i in `seq -w 1 10|head -n 6`
do
echo $i
F1=sampled_raw/sample_${i}_1.fastq.gz
F2=sampled_raw/sample_${i}_2.fastq.gz
FN=sample_${i}
echo $FN
rsem-calculate-expression -p 20 --paired-end --bowtie2 --estimate-rspd --append-names --output-genome-bam $F1 $F2 $transcripts $FN
done

#data statistics
for i in `seq -w 1 10|head -n 6`
do
FN=sample_${i}
rsem-plot-model ${FN} ${FN}_diagnostic.pdf
done

#Differential Expression Analysis using EBSeq
#differentially expressed isoforms FDR at level 0.05
rsem-generate-ngvector genome/RSEM.mm10/mm10.transcripts.fa human_ref
rsem-generate-data-matrix sample_01.isoforms.results sample_02.isoforms.results sample_03.isoforms.results sample_04.isoforms.results sample_05.isoforms.results sample_06.isoforms.results > IsoMat.txt
rsem-run-ebseq --ngvector human_ref.ngvec  IsoMat.txt 3,3 IsoMat.results
rsem-control-fdr IsoMat.results 0.05 IsoMat.de.txt





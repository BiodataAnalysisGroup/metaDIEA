
# hisat2 aligner
# https://ccb.jhu.edu/software/hisat2/index.shtml

# StringTie Transcript assembly and quantification for RNA-Seq
# http://www.ccb.jhu.edu/software/stringtie/

# DESeq2 R package
# https://bioconductor.org/packages/release/bioc/html/DESeq2.html

cd /app

echo "*** Mapping ***"
# 1| Map the reads for each sample to the reference genome:
for i in `seq -w 1 10|head -n 6`
do
echo $i
hisat2 -p 16 --dta -x genome/genome -1 raw/sample_${i}_1.fastq.gz -2  raw/sample_${i}_2.fastq.gz  -S sample_${i}.sam &> sample_${i}.log.txt 
done

echo "*** Sorting ***"
# 2| Sort and convert the SAM files to BAM:
for i in `seq -w 1 10|head -n 6`
do
echo $i
samtools sort -@ 16 -o sample_${i}.bam sample_${i}.sam
done

rm *sam

echo "*** Estimating transcript abundances ***"

# 3| Estimate transcript abundances and create table counts for Ballgown:
for i in `seq -w 1 10|head -n 6`
do
echo $i
stringtie -e -B -p 16 -G genome/genes.gtf -o ballgown/sample_${i}/sample_${i}.gtf sample_${i}.bam 
done


wget http://www.ccb.jhu.edu/software/stringtie/dl/prepDE.py
chmod +x prepDE.py 
./prepDE.py

echo ids,type >samples.csv
echo sample_01,control >>samples.csv
echo sample_02,control>>samples.csv
echo sample_03,control>>samples.csv
echo sample_04,treated>>samples.csv
echo sample_05,treated>>samples.csv
echo sample_06,treated>>samples.csv

# run DESeq2 script
Rscript Deseq2.hisat.RealQuals.simulated.imb.R .


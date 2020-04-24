# OVERALL PIPELINE
# pipeline was separated into multiple .sh/.slurm files
# as a result, repetitive parts such as module load or cd were excluded
# loaded modules where needed

########################## pre-processing ############################

# wget SRA manually

# quality control
fastqc /scratch/pl9ed/project/wes_data/*.fastq

# fastq > bam
# bowtie2
fileList=($(cat /scratch/pl9ed/project/list.txt))

n=$SLURM_ARRAY_TASK_ID
fileName=${fileList[$n]}

echo "aligning " $fileName
bowtie2 -p 8 -x /scratch/pl9ed/project/refgenome/ref -1 ./"$fileName"_1.fastq -2 ./"$fileName"_2.fastq -S $fileName.sam

########################## variant discovery ############################
# bam > sort > mpileup
# samtools
echo "converting $fileName to bam"
samtools view -bS "$fileName".sam > $fileName.bam
samtools sort -@ 16 "$fileName".bam > "$fileName".sorted.bam
samtools index -@ 16 "$fileName".sorted.bam
samtools mpileup -E -@ 16 -uf /scratch/pl9ed/project/refgenome/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
	"$fileName".sorted.bam > "$fileName".mpileup 

########################## variant prioritization ############################
# VarScan
echo "-------------mpileup2snp-------------"
java -jar VarScan.v2.4.4.jar mpileup2snp $fileName.mpileup --variants > $fileName.v.snp
echo "-------------mpileup2indel-------------"
java -jar VarScan.v2.4.4.jar mpileup2indel $fileName.mpileup --variants > $fileName.v.indel
echo "-------------filter snp-------------"
java -jar VarScan.v2.4.4.jar filter $fileName.v.snp --indel-file $fileName.v.indel --output-file $fileName.v.snp.filter
echo "-------------filter indel-------------"
java -jar VarScan.v2.4.4.jar filter $fileName.v.indel --output-file $fileName.v.indel.filter
echo "-------------readcounts-------------"
java -jar VarScan.v2.4.4.jar readcounts $fileName.mpileup --output-file $fileName.readcounts

# bcftools
echo "Working on $fileName"
echo "n: $n"

bcftools mpileup -f /scratch/pl9ed/project/refgenome/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
	$fileName.sorted.bam --threads 16 -o $fileName.bcf
echo "-------calls-------"
bcftools call $fileName.bcf --threads 16 -mv -Ob -o $fileName.calls.bcf
echo "-------filter-------"
bcftools filter -i'%QUAL>20' $fileName.calls.bcf -Ov -o $fileName.calls.filtered.vcf

########################## annotate ############################
fileList=($(cat /scratch/pl9ed/project/list.txt))
cList=($(cat /scratch/pl9ed/project/clist.txt))
n=$SLURM_ARRAY_TASK_ID

fileName=${fileList[$n]}
cName=${cList[$n]}

echo "fileName: $fileName"
echo "cName: $cName"

cd /scratch/pl9ed/project/annovar

# vcf to avninput - easier to process later
./convert2annovar.pl -format vcf4 ../wes_data/$fileName.calls.filtered.vcf > $fileName.avinput
./convert2annovar.pl -format vcf4 ../control/$cName.calls.filtered.vcf > $cName.avinput

# annotate
./table_annovar.pl $fileName.avinput humandb/ \
	-buildver hg19 \
	-out $fileName \
	-remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a \
	-operation gx,r,f,f,f \
	-nastring . \
	-csvout -polish -xref example/gene_xref.txt
./table_annovar.pl $cName.avinput humandb/ \
	-buildver hg19 \
	-out $cName \
	-remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a \
	-operation gx,r,f,f,f \
	-nastring . \
	-csvout -polish -xref example/gene_xref.txt

# SNV/Indel counts:
# ug.csv and glist.csv from R script
# create CSV of [Gene, Count] for each sample
cd ../annovar

echo "Working on $oName and $cName"

for i in $(cat ../av_scripts/ug.csv)
do

gcount=$(grep $i $oName*.csv | wc -l)
echo "$i, $gcount" >> ../gene_counts/$oName.count.csv

ccount=$(grep $i $cName*.csv | wc -l)
echo "$i, $ccount" >> ../gene_counts/$cName.count.csv
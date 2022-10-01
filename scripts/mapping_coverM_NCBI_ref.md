##2021-10-10
```
#@ptpe2
source /home/PTPE2/Software/miniconda3/bin/activate
conda activate bowtie2

#wkdir
cd /home/PTPE2/User/liupf/Projects_liupf/glacier_virome/mapping_NCBIref208_vgenome

#https://www.ncbi.nlm.nih.gov/refseq/
#1008-2021, 2021
#September 13, 2021, RefSeq Release 208 is available for FTP
#1008-2021, 2021
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.2.1.genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.3.1.genomic.fna.gz
#

#gzip

gzip -d viral.1.1.genomic.fna.gz
gzip -d viral.2.1.genomic.fna.gz
gzip -d viral.3.1.genomic.fna.gz

cat *.fna > viral.123.1.genomicRef208.fasta
#14721

sed -i -e 's/ .*$//g' viral.123.1.genomicRef208.fasta

#file                             format  type  num_seqs      sum_len  min_len   avg_len    max_len
#viral.123.1.genomicRef208.fasta  FASTA   DNA     14,721  459,807,044      136  31,234.8  2,473,870

cat viral.123.1.genomicRef208.fasta |seqkit seq -m 3000 > viral.123.1.genomicRef208m3000.fasta

seqkit stat  viral.123.1.genomicRef208m3000.fasta

bowtie2-build viral.123.1.genomicRef208m3000.fasta viral.123.1.genomicRef208m3000

#raw reads location
cd /home/PTPE2/Project/glacier_metaGV_166/85TPmeta_cleandata

/mapping_NCBIref208_TP85_metadata.txt

./bowtie2_NCBIref208_vgenome.sh mapping_NCBIref208_TP85_metadata.txt 70
#edits


for samfile in *.sam
do
samtools view -@ 60 -bS ${samfile} > "${samfile%.*}".bam
samtools sort -@ 60 "${samfile%.*}".bam > "${samfile%.*}".bam.sorted
done


```


##mapping data of Arctic and Alps to NCBI ref; 45 samples
```
cd /home/PTPE2/User/liupf/Projects_liupf/glacier_virome/mapping_NCBIref208_vgenome


/home/PTPE2/Project/glacier_metaGV_166/45NonTP

./bowtie2_NCBIref208_vgenome.sh mapping_metafile_non-TP45.txt 70

```

#partII + part III another 3 samples
```
./bowtie2_NCBIref208_vgenome.sh mapping_metafile_partII_III.txt 30


```

##mapping data of Antarctic
```
./bowtie2_NCBIref208_vgenome.sh Anartic_TS_mapping_metafile.txt 70


```

###filtering low quality reads and both paired mapped, and properly mapped
##bowtie2 keep singletons that only one pair matched
#check reads mapped only one time
```
#https://www.biostars.org/p/138116/
# samtools flagstat， check mapping stats

#remove singletons -F 0x08； remove low quality， -q 30； #-f 0x2 properly mapped
#samtools view -h -@ 8 -q 30 -F 0x08 -b -f 0x2 SRR1517848.bam > SRR1517848.bam

#missing code here for filtering low quality reads
sortedf


#coverM calculate

```


#2022-05-26 mapping reads of 21 cryoconites samples
```
/home/PTPE2/User/liupf/Projects_liupf/glacier_virome/mapping_NCBIref208_vgenome


PTPE2@user-PowerEdge-R740:~/User/liupf/Projects_liupf/glacier_virome/mapping_NCBIref208_vgenome$ cp /home/PTPE2/User/liupf/Projects_liupf/glacier_virome/mapping_glacier_vOTU/glacier188_bowtie2_mapping/glacie_cryo21_cleanreads_dir.txt ./


#
screen -S bowtie2_NCBIref208_vgenome

source /home/PTPE2/Software/miniconda3/bin/activate

conda activate bowtie2

./bowtie2_NCBIref208_vgenome.sh glacie_cryo21_cleanreads_dir.txt 50

for samfile in *.sam
do
samtools view -@ 60 -bS ${samfile} > "${samfile%.*}".bam

samtools view -h -@ 8 -q 30 -F 0x08 -b -f 0x2 "${samfile%.*}".bam > "${samfile%.*}"f.bam

samtools sort -@ 60 "${samfile%.*}"f.bam > "${samfile%.*}".bam.sortedf

done

~/Software/coverm061/coverm contig --methods count --bam-files *.sortedf --min-covered-fraction 0 --output-file glacier_188_NCBI208_count_filtered.txt -t 60

#get reads number
PTPE2@user-PowerEdge-R740:~/ncbi/public/sra$ wc -l glacier_cryo_21_reads_number.txt
21 glacier_cryo_21_reads_number.txt

```



##coverM counts of mapped reads
```
##

～/Software/coverm061/coverm contig --methods counts --bam-files *.sorted --min-covered-fraction 50 --output-file glacier_167_NCBI208_count50.txt -t 60
#count could use --min-covered fraction
#
#rename 's/glacier_TP85_allvTOU_TPM/glacier_TP85_NCBIRef26_TPM/' *-comverm.txt

# counts based
--bam-files *.sorted  --output-file glacier-t


#cd to TP85 raw reads
seqkit stats -j 50 *.fq* >read.count.txt
seqkit stats -j 50 *.fastq >read.count_fastq.txt
seqkit stats -j 50 1908KQGRS1_R*.gz >read.count_gz.txt

```


##


#phylogeny of Caudovirales
```
#1, VOG 77 HMMS Caudovirales
vog release vog83 (same as ncbi release 83) were used in Low et al., 2019, since the vog id is changing, use the same release as Low et al., 2019 (not as in the paper said 2017 April 1st)

#2, reference sequences, ncbi all Caudovirales sequences
#3, all glacier vOTUs

#criteria
#1, Individual marker alignments were then trimmed to retain positions with less than 50% gaps using trimAl v.1.4 (ref. 90) and concatenated, filling in gaps for missing markers where necessary. Only genomes containing at least three markers and having data at >5% of alignment columns were retained.

#2, The hit with the highest bit score was identified as the best hit and, in cases with multiple hits, the highest scoring sequence was selected. Only genomes with at least three markers were retained for further analysis.
```

#wkdir
```
/home/liupf/Projects/glacier_virome/Cau_tree

#refseqs

cat Caudovirales_taxid10239_su134.nbr|cut -f1 -d$'\t'|tr "," "\n" |tr -d "^ " > Caudovirales_taxid10239_su134_F1.nbr
#4329

conda activate seqkit

seqkit grep -r -f Caudovirales_taxid10239_su134_F1.nbr viral123genomic_ref207.fna -o Caudovirales_taxid10239_su134_F1.fna
#4328 seqs

sed -i -e 's/ .*$//g' Caudovirales_taxid10239_su134_F1.fna

```
#0, get protein sequence for genomes
```
# prodigal
prodigal -p meta -q -i Caudovirales_taxid10239_su134_F1.fna  -a Caudovirales_taxid10239_su134_F1.faa

sed -i -e 's/ .*$//g' Caudovirales_taxid10239_su134_F1.faa

#glacier vOTUs proteins, from vContact2

/home/liupf/Projects/glacier_virome/vContact2/glacier_virome_vOTU_vc2

sed -e 's/ .*$//g' glacier_total_virus_final123_95-85.faa > /home/liupf/Projects/glacier_virome/Cau_tree/glacier_total_virus_final123_95-85fix.faa

cat Caudovirales_taxid10239_su134_F1.faa glacier_total_virus_final123_95-85fix.faa > Cau_hmm_search.faa
sed -i -e 's/ .*$//g'  Cau_hmm_search.faa
```


#1, hmmsearch, use vog83 release hmm profile
```
for file in /home/liupf/Projects/glacier_virome/Cau_tree/vog77hmm/*.hmm
do
hmmsearch --cpu 10 -E 0.001 -Z 2232 --tblout "${file%.*}".txt $file Cau_hmm_search.faa

grep -v '#' "${file%.*}".txt > "${file%.*}"_fix.txt

awk -F ' ' '{print $1,"\t",$0}' "${file%.*}"_fix.txt> "${file%.*}"_fix2.txt

sed -i -e 's/\(.*\)_[0-9].*\t/\1\t/g' "${file%.*}"_fix2.txt
sed -i -e 's/ \+ /\t/g' "${file%.*}"_fix2.txt

#get the top hits
awk '!x[$1]++' "${file%.*}"_fix2.txt > "${file%.*}"_top.txt

rm *fix*.txt
done

```

#test on ref
#compare with Low et al, for e-values, Z, 2232
```
VOG00123.hmm
hmmsearch --cpu 4 --tblout VOG00123.hmm.txt -E 0.001 -Z 2232 ./vog77hmm/VOG00123.hmm  Caudovirales_taxid10239_su134_F1.faa #-E 0.001 -Z 2232
#331
#in Low et al., 2019, 1520 tree, this is 263; Evaluation of a concatenated protein phylogeny for classification of tailed double-stranded DNA viruses belonging to the order Caudovirales
#full dataset of 1804 genomes, this is 338
#

for file in *.txt
do
grep -v '#' $file > "${file%.*}"_fix.txt

awk -F ' ' '{print $1,"\t",$0}' "${file%.*}"_fix.txt> "${file%.*}"_fix2.txt

sed -i -e 's/\(.*\)_[0-9].*\t/\1\t/g' "${file%.*}"_fix2.txt
sed -i -e 's/ \+ /\t/g' "${file%.*}"_fix2.txt

#get the top hits
awk '!x[$1]++' "${file%.*}"_fix2.txt > "${file%.*}"_top.txt

rm *fix*.txt
done
```



#2, top hits
```
#https://www.biostars.org/p/453514/#google_vignette
#http://slhogle.github.io/2015/remove-duplicate-lines/
#https://snakemake-wrappers.readthedocs.io/en/latest/wrappers/hmmer/hmmsearch.html

awk '!x[$3]++' ouput_file.pfam > MYBESTHITS.pfam

```

#3, get sequences and do alignments
```
#get sequences
for file in *_top.txt
do
cat $file |cut -f2 -d$'\t' > "${file%.*}"_tophits_name.txt
sed -i -e 's/ \+//g' "${file%.*}"_tophits_name.txt
sed -i -e 's/-$//g' "${file%.*}"_tophits_name.txt

seqkit grep -n -f "${file%.*}"_tophits_name.txt ../Cau_hmm_search.faa -o "${file%.*}"_tophits.fasta
done

#check
seqkit grep -r -f VOG00022_top_tophits_name.txt ../Cau_hmm_search.faa -o VOG00022_top_tophits2.fasta

VOG00022_top_tophits.fasta.name
VOG00022_top_tophits_name.txt

comm -13 <(sort VOG00022_top_tophits.fasta.name) <(sort VOG00022_top_tophits_name.txt)

#fix
for file in *.fasta
do
sed -i -e 's/\*//g' ${file}
done

# hmm alignment
#hmmalign [-options] <hmmfile> <seqfile>

for file in *.hmm
do
hmmalign -o "${file%.*}".stock $file "${file%.*}"_top_tophits.fasta
esl-reformat -o "${file%.*}".clustal clustal "${file%.*}".stock
done


conda activate emboss6

#reforamt to fasta
for file in *.clustal
do
seqret -sequence "${file%.*}".clustal -outseq "${file%.*}".align.fasta
done



test
hmmalign -o VOG00022.align VOG00022.hmm VOG00022_top_tophits.fasta
esl-reformat -o VOG00022.clustal clustal VOG00022.align
conda activate emboss6
seqret -sequence VOG00022.clustal -outseq VOG00022.align.faa

```

#1, Individual marker alignments were then trimmed to retain positions with less than 50% gaps using trimAl v.1.4 (ref. 90) and concatenated, filling in gaps for missing markers where necessary. Only genomes containing at least three markers and having data at >5% of alignment columns were retained.

#4, trimming
```
conda activate trimal

#keep on genome names for concate
for file in *align.fasta
do
#sed -i -e 's/\(.*\)_[0-9].*$/\1/g' $file
trimal -keepheader -gt 0.5 -in $file -out "${file%.*}"_trim.fasta
done


```


#4, concatnate genomes genomes
```
#PhyloSuite_v1.2.1_Mac
#seqkit concat *_trim.fasta> concat_alignment.faa #not working with fill missing data with gaps

#catfasta2phyml.pl
/home/liupf/scripts/catfasta2phyml.pl -f --concatenate *_trim.fasta  >concat_alignment.faa

# 13870

```


#5, filtering out genomes

#>5% of sequences
```
#on MAC, server is down
/Users/pengfeiliu/Code_scripts_backup/Bioinformatics_demo_workshop/Virus_metaG_pipeline

#filtering by gaps number
python fasta_drop.py concat_alignment.faa concat_alignment_5.faa 0.95
#6459

#inlcuding ref 4112

#ptpe is back
cd /home/liupf/Projects/glacier_virome/Cau_tree/vog77hmm

python /home/liupf/scripts/fasta_drop.py concat_alignment.faa concat_alignment_5.faa 0.95

```

#get a table of marker genes
##keep seq with 3 markers

```
grep '>' concat_alignment_5.faa > concat_alignment_5_seqname.txt

wc -l concat_alignment_5_seqname.txt
6459 concat_alignment_5_seqname.txt

for file in *align_trim.fasta
do
grep -f concat_alignment_5_seqname.txt $file > "${file%.*}"_seqname.txt

sed -i "1i "${file%.*}"" "${file%.*}"_seqname.txt

done

for file in VOG*_seqname.txt
do

#sed -i -e 's/>\(.*\)/>\1\t1/g' ${file}
#sed -i -e 's/VOG\(.*\)/Name\tVOG\1/g' ${file}
sed -i -e 's/Name\t/Name\tVOG/g' ${file}
done

sed -i -e '1i Name' concat_alignment_5_seqname.txt

cp concat_alignment_5_seqname.txt tmp.txt

for file in VOG*_seqname.txt
do
python merge_columns_together_by_id.py tmp.txt ${file} > tmp1.txt

mv tmp1.txt tmp.txt

done

nano glacier_vOTUs_Cau.list

# 6333
```


#build a phyogeny tree
```
#https://github.com/dongzhang0725/PhyloSuite
#
seqkit

seqkit grep -f glacier_vOTUs_Cau.list concat_alignment_5.faa -o glacier_w_ref_concat_alignment_5_gt3.faa

#6333, with 25563 positions

FastTreeMP -gamma -wag -n 1000 glacier_w_ref_concat_alignment_5_gt3.faa > glacier_w_ref_concat_alignment_5_gt3b1000.tree

```



############################################preparation of ref sequences########################
##NCBI ref 207
```
#on MAC
/Users/pengfeiliu/Code_scripts_backup/Bioinformatics_demo_workshop/Virus_metaG_pipeline

cat  taxid10239.nbr |cut -f1,3-5 -d$'\t' |sort| uniq > taxid10239_su.nbr
cat  taxid10239.nbr |cut -f1,3-4 -d$'\t' |sort| uniq > taxid10239_su134.nbr


cat  taxid10239.nbr |cut -f1 -d$'\t' |sort| uniq |wc -l

# 14058 taxid10239_su134.nbr
```


#vog hmm ref vog83
```
#/Users/pengfeiliu/Code_scripts_backup/Bioinformatics_demo_workshop/Virus_metaG_pipeline/vog.hmm

```


##get list of Caudovirales

#Caudovirales
```
##https://en.wikipedia.org/wiki/Caudovirales
#https://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=28883
Caudovirales - 4311 complete genomes	Retrieve sequences:
-- Select data set from the list --
Ackermannviridae  [67]	Autographiviridae  [373]	Chaseviridae  [14]	Demerecviridae  [96]
Drexlerviridae  [114]	Guelinviridae  [9]	Herelleviridae  [128]	Lilyvirus  [1]
Myoviridae  [860]	Podoviridae  [243]	Rountreeviridae  [35]	Salasmaviridae  [34]
Schitoviridae  [84]	Siphoviridae  [2237]	Zobellviridae  [13]	environmental samples  [1]
unclassified Caudovirales [2]

$ touch Caudovirales_NCBI_ref207.txt

Ackermannviridae
Autographiviridae
Chaseviridae
Demerecviridae
Drexlerviridae
Guelinviridae
Herelleviridae
Lilyvirus
Myoviridae
Podoviridae
Rountreeviridae
Salasmaviridae
Schitoviridae
Siphoviridae
Zobellviridae
Caudovirales

##
grep -f Caudovirales_NCBI_ref207.txt taxid10239_su134.nbr > Caudovirales_taxid10239_su134.nbr

#Caudovirales_taxid10239_su134.nbr
#4264 Caudovirales_taxid10239_su134.nbr


#direct download table
taxid28883.acc_lst
#4334 taxid28883.acc_lst

#RefSeq(4,380)
```



#get sequences
```
taxid28883.acc_lst
```




## ref

```
#
Low et al., 2019-Nature Microbiology-Evaluation of a concatenated protein phylogeny for classification of tailed double-stranded DNA viruses belonging to the order Caudovirales

```

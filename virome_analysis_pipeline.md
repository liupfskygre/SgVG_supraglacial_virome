## virome analysis pipeline for the global glacier viruses

#main scripts for the analysis of the GGV catalog

#by Pengfei Liu, Lin Zang, Zhihao Zhang, Tao Ye


#1, assembly of microbial and viral metagenomics by megahit
```
#
#trimming reads 

for samplename in $(cat trim_sample.tsv)
do trimmomatic PE -phred33 ${samplename}.R1.fastq.gz ${samplename}.R2.fastq.gz ${samplename}_trimmed_R1.fastq.gz ${samplename}_trimmed_U_R1.fastq.gz ${samplename}_trimmed_R2.fastq.gz ${samplename}_trimmed_U_R2.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 -threads 20
done

#
#assembly with megahit

megahit -1 "${read1}" -2 "${read2}" --k-min 21 --k-max 141 --k-step 10 --mem-flag 2 -m 0.5 --out-prefix "${line}" -t 40 --out-dir ./"${line}" &> "${line}"_megahit.log

```


#2, viral contigs identification from metagenics derived assemblies
```
#1, virsorter2

virsorter run -w glacier_10k—vs2.out -i glacier_10k.fasta -j 120 --prep-for-dramv --min-length 10000 --include-groups dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae -d /mnt/nfs/software/VirSorter2/db  classify  &> glacier_10k—vs2_reclassify.log

#2, VIBRANT
python3 VIBRANT_run.py -i glacier_metaG10.fasta -t 20 -l 10000


#3, virfinder
R
library(VirFinder)
predResult <- VF.pred("./glacier_metaG10.fasta")
predResult[order(predResult$pvalue),]
predResult$qvalue <- VF.qvalue(predResult$pvalue)
predResult_sorter<-predResult[order(predResult$qvalue),]
write.table(predResult_sorter, file="virfinder-metaG10.txt", quote=F, sep="\t")


#4, DeepVirFinder
python /home/ptpe/software/DeepVirFinder/dvf.py -i glacier_10k.fasta -o glacier_10k_Deepvirfinder -l 10000 -c 200 &>glacier130.deepvirfinder.log


#5, viralrecall
python /home/PTPE2/Software/viralrecall/viralrecall/viralrecall.py -i /home/PTPE2/User/liupf/Projects_liupf/glacier_virome/mapping_glacier_vOTU/glacier166_bowtie2_mapping/reference/glacier_13475vOTUs.fna -p glacier_13475vOTUs_outdir -t 20 -f &> log
```



#3, quality checking of viral contigs by CheckV
```
checkv end_to_end glacier_total_virus_Can.fasta checkv_output_directory -t 30 -d /home/ptpe/biodatabase/checkv_db/checkv-db-v1.0

#discarded list: (virsorter2 SOP: https://www.protocols.io/view/viral-sequence-identification-sop-with-virsorter2-5qpvoyqebg4o/v2)
#virus gene==0 and host gene>=1

awk '$6==0&$7>0' quality_summary.tsv | cut -f1 -d$'\t' > glacier_total_virus_discarded_seq.txt

#kept
seqkit grep -n -v -f glacier_total_virus_discarded_seq.txt  ../virus_tmp.fasta -o glacier_total_virus_final1.fasta

```


#4, clustering of viral contigs to generate vOTUs, and the taxnomic classifcation and annotation of vOTUs
```
#Roux et al., 2020-Nucleic Acids Research-IMG/VR v3: an integrated ecological and evolutionary framework for interrogating genomes of uncultivated viruses
#Clustered into vOTUs following the MIUViG guidelines (14) (95% ANI - Average Nucleotide Identity and 85% AF – Aligned Fraction).

First, create a blast+ database: makeblastdb -in <my_seqs.fna> -dbtype nucl -out <my_db>
makeblastdb -in global_virus_contigs_partI5k.fna -dbtype nucl -out global_virus_contigs_partI5k


Next, use megablast from blast+ package to perform all-vs-all blastn of sequences: 
blastn -query global_virus_contigs_partI5k.fna -db global_virus_contigs_partI5k -outfmt '6 std qlen slen' -max_target_seqs 15000 -o global_virus_contigs_partI5k_blast.tsv -num_threads 32

# blastn -query <my_seqs.fna> -db <my_db> -outfmt '6 std qlen slen' -max_target_seqs 10000 -o <my_blast.tsv> -num_threads 32

Note: using the -perc_identity flag will speed up the search at the cost of sensitivity: blastn -query <my_seqs.fna> -db <my_db> -outfmt '6 std qlen slen' -max_target_seqs 10000 -perc_identity 90 -o <my_blast.tsv> -num_threads 32

Next, calculate pairwise ANI by combining local alignments between sequence pairs: anicalc.py -i <my_blast.tsv> -o <my_ani.tsv>

conda activate vibrant
python ~/User/liupf/Projects_liupf/common_files/anicalc.py -i global_virus_contigs_partI5k_blast.tsv -o global_virus_contigs_partI5k_ani.tsv

global_virus_contigs_partI5k_blast.tsv

Finally, perform UCLUST-like clustering using the MIUVIG recommended-parameters (95% ANI + 85% AF): aniclust.py --fna global_virus_contigs_partI5k.fna --ani global_virus_contigs_partI5k_ani.tsv --out global_virus_contigs_partI5k_clusters.tsv --min_ani 95 --min_tcov 85 --min_qcov 0

python ~/User/liupf/Projects_liupf/common_files/aniclust.py --fna global_virus_contigs_partI5k.fna --ani global_virus_contigs_partI5k_ani.tsv --out global_virus_contigs_partI5k_clusters.tsv --min_ani 95 --min_tcov 85 --min_qcov 0


################taxnomy classifiction of vOTUs
#vContact2
prodigal -i glacier_total_virus_final123_5k_95-85.fna  -a glacier_total_virus_final123_95-85.faa -p meta
vcontact2_gene2genome -p glacier_total_virus_final123_95-85.faa -o glacier_total_virus_final123_95-85.csv -s Prodigal-FAA
vcontact2 --raw-proteins glacier_total_virus_final123_95-85.faa --rel-mode Diamond --proteins-fp glacier_total_virus_final123_95-85.csv --db 'ProkaryoticViralRefSeq201-Merged' --pcs-mode MCL --vcs-mode ClusterONE --c1-bin /home/ptpe/software/cluster_one-1.0.jar  --output-dir glacier_virome_vOTU_vc2  -t 60


#blastn
makeblastdb -in IMGVR_all_nucleotides_replace.fna -dbtype nucl -out IMGVR_all_nucleotides -parse_seqids
blastn -task megablast -query glacier_total_virus_final123_5k_95-85.fna -db /home/ptpe/biodatabase/IMG_VR_2020-10-12_5.1/IMGVR_all_nucleotides -evalue 1e-5 -perc_identity 90 -max_target_seqs 25000  -num_threads 20 -out glacier_total_virus_final123_5k_95-85-IMGvr.txt -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"

#blastp
blastp -query glacier_total_vOTU_prot_fixheader.faa -db /home/PTPE2/User/liupf/Projects_liupf/glacier_virome/NCBI-ref206-virus/viral.123.protein -evalue 1e-3 -max_target_seqs 577642 -num_threads 40 -out glacier_total_vOTU_prot-NCBIvprotein.txt -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"

#VPF-classs
screen -r vpf-class
/home/ptpe/software/vpfclass/vpf-class-x86_64 --data-index /home/ptpe/software/vpfclass/data/data-index.yaml -i ./glacier_total_virus_final123_5k_95-85.fna -o vpfclass-glacier-vOTU --workers 20 --chunk-size 20

#Cau tree (see ** file for details)


#virrecall for NCLDVs
python viralrecall.py -i /home/PTPE2/User/liupf/Projects_liupf/glacier_virome/mapping_glacier_vOTU/glacier166_bowtie2_mapping/reference/glacier_13475vOTUs.fna -p glacier_13475vOTUs_outdir -t 60 -f &> liupf_log

```


#5, in silico prediction of the life style of glacial viruses
```
#checkV, see above output

#VIBRANT, see above output

#bacphlip, 95% cutoff
bacphlip -i glacier_April1_final5k_vOTUs_Cryo21_Anta1.fna --multi_fasta

```


#6, marcordiversity and microdiverstiy analysis of glacial viral communities
```
#####marcordiversity 
##mapping with bowtie2
./bowtie2_2_vgenome.sh mapping_vOTUs_TP85_metadata.txt 24  #see scripts for sh files

#quality filtering

#coverage  ==coverM
conda activate coverm
coverm contig --methods tpm --bam-files *.sorted --min-covered-fraction 10 --output-file
#input into R, microeco package for alpha and beta diversity analysis 

####microdiverstiy
##mapping with bowtie2, same as above for
./bowtie2_2_vgenome.sh mapping_vOTUs_TP85_metadata.txt 60 #see scripts for sh files


#metapop
metapop --input_samples ./bamfile --reference ./reference --norm tp-notp-187-metapop_ctfile.txt --threads 50 --whole_genomes -o all187_metapop_out_whole &> all187_metapop_out_whole.log
```


#7, in silico host prediction of viral vOTUs
```
###1, CRISPR spacer to TG2G genomes databases
#get spacers by CRT

for file in /home/liupf/Projects/Zanglin_virus/MAGs_final/*.fa
do
java -cp /home/ptpe/software/CRT1.2-CLI.jar crt -minNR 4 -maxRL 55 -maxSL 70 $file "${file%.*}".txt
done

####custome provided spacers
##Parsing spacer files

conda activate spacepharer
spacepharer parsespacer *.txt MAGs_1020_queryDB --threads 5

spacepharer createsetdb MAGs_1020_queryDB MAGs_1020_querySetDB tmpFolder --extractorf-spacer 1


###
cd /home/liupf/Projects/glacier_virome/MAGs-TP-2661ana/MAGs_1020_drep099
cd /home/liupf/Projects/glacier_virome/MAGs-TP-2661ana/MAGs_1020_drep099/MAGs_1020_drep099_CRT

spacepharer predictmatch MAGs_1020_querySetDB /home/liupf/Projects/glacier_virome/SpacePharer_MAGs203/glacier_total_virus_final123_5k_95-85.fna.split/targetSetDB /home/liupf/Projects/glacier_virome/SpacePharer_MAGs203/glacier_total_virus_final123_5k_95-85.fna.split/controlSetDB MAGS-1020-outputFileName.tsv tmpFolder


####2, blastn based
blastn -task megablast -query /home/liupf/Projects/glacier_virome/vContact2/glacier_total_virus_final123_5k_95-85.fna -db glacier1020_nucleotides -evalue 1e-5 -perc_identity 90 -max_target_seqs 25000  -num_threads 20 -out glacier_total_vOTUs_MAGs1020.txt -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
#Host predictions were then based on matches of ≥90% nucleotide identity covering ≥2 kb of the virus and (putative) host sequences.


##virus host matcher, VHM
python /home/ptpe/software/VirHostMatcher/vhm.py -v /home/liupf/Projects/glacier_virome/VirHostMatcher/virus -b /home/liupf/Projects/glacier_virome/VirHostMatcher/host -o /home/liupf/Projects/glacier_virome/VirHostMatcher/output-test2 &>v.log2

##tRNA blast match
#tRNA from host
tRNAscan-SE -G -Q --thread 50 -a "${file%.*}"_trna.fasta scaffolds.fna

#tRNA from virus
tRNAscan-SE -G -Q --thread 10 -a vOTUs13839_trna.fasta glacier_total_virus_final123_5k_95-85.fna

makeblastdb -in TP-nonTP-tRNA_replace.fna -dbtype nucl -out TP-nonTP-tRNA_MAGs -parse_seqids
#blastn
blastn -task blastn -query vOTUs13839_trna.fna -db TP-nonTP-tRNA_MAGs -evalue 1e-5 -perc_identity 100 -max_target_seqs 50000  -num_threads 30 -out tRNA_vOTUs_blastn_MAGs1364.txt -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"

#100% perfect match between the tRNA sequences of viruses and host


```

#8, Auxiliary metabolic genes (AMGs) of viral vOTUs
```
#run virsorter with --min-score 0
virsorter run --seqname-suffix-off --viral-gene-enrich-off -w glacier_total_virus_final123_5k_95-85—vs2.out -i glacier_total_virus_final123_5k_95-85.fna --prep-for-dramv --min-length 5000 --min-score 0 -j 40 --include-groups dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae -d /home/ptpe/biodatabase/vs2_db  all  &> glacier_total_virus_final123_5k_95-85—vs2_all.log

#DRAM-v 
DRAM-v.py annotate -i final-viral-combined-for-dramv.part_001.fa -v final-viral-combined-for-dramv.part_001_for-dramv.tab -o final-viral-combined-for-dramv.part_001_annotation --use_uniref --threads 2 --min_contig_size 5000

```

#9, metabolic profile of MAGs
```
#
DRAM.py annotate -i './*.fa' -o Glacier_CRISPR_linked_MAGs --threads 40 --min_contig_size 500 --checkm_quality ~/User/liupf/Projects_liupf/glacier_virome/MAGs-Non-TP-405ana/NonTP-MAGs_344_drep099_checkm_summary.txt --checkm_quality ~/User/liupf/Projects_liupf/glacier_virome/MAGs-TP-2661ana/TP-MAGs_1020_drep099_checkm_summary.txt  --gtdb_taxonomy ~/User/liupf/Projects_liupf/glacier_virome/MAGs-Non-TP-405ana/gtdbtk_out/NonTP-MAGs_344_drep099.bac120.summary.tsv --gtdb_taxonomy ~/User/liupf/Projects_liupf/glacier_virome/MAGs-TP-2661ana/MAGs_1020_drep099/gtdbtk_out/TP-MAGs_1020_drep099.bac120.summary.tsv --gtdb_taxonomy ~/User/liupf/Projects_liupf/glacier_virome/MAGs-TP-2661ana/MAGs_1020_drep099/gtdbtk_out/TP-MAGs_1020_drep099.ar122.summary.tsv --verbose --use_uniref

```

#10, mapping to NCBI ref genome to assess the potential risk to public health of glacier viruses
```
#bowtie2


```

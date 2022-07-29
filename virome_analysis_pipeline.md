## virome analysis pipeline for the global glacier viruses

#by Pengfei Liu, Lin Zang, Zhihao Zhang, Tao Ye

#1, assembly of microbial and viral metagenomics by MEGAhit
```
#
#trimming reads 

for samplename in $(cat trim_sample.tsv)
do trimmomatic PE -phred33 ${samplename}.R1.fastq.gz ${samplename}.R2.fastq.gz ${samplename}_trimmed_R1.fastq.gz ${samplename}_trimmed_U_R1.fastq.gz ${samplename}_trimmed_R2.fastq.gz ${samplename}_trimmed_U_R2.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 -threads 20
done

#
#assembly with megahit

megahit -1 "${read1}" -2 "${read2}" --k-min 41 --k-max 141 --k-step 10 --mem-flag 2 -m 0.5 --out-prefix "${line}" -t 40 --out-dir ./"${line}" &> "${line}"_megahit.log

```


#2, viral contigs identification from metagenics derived assemblies
```
#virsorter2

virsorter run -w glacier_10k—vs2.out -i glacier_10k.fasta -j 120 --prep-for-dramv --min-length 10000 --include-groups dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae -d /mnt/nfs/software/VirSorter2/db  classify  &> glacier_10k—vs2_reclassify.log

#DeepVirFinder
python /home/ptpe/software/DeepVirFinder/dvf.py -i glacier_10k.fasta -o glacier_10k_Deepvirfinder -l 10000 -c 200 &>glacier130.deepvirfinder.log


```



#3, quality checking of viral contigs by CheckV
```


```


#4, clustering of viral contigs to generate vOTUs, and the taxnomic classifcation and annotation of vOTUs
```


```



#5, in silico prediction of the life style of glacial viruses
```


```


#6, marcordiversity and microdiverstiy analysis of glacial viral communities
```


```


#7, in silico host prediction of viral vOTUs
```


```

#8, Auxiliary metabolic genes (AMGs) of viral vOTUs
```


```

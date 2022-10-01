
#!/usr/bin/bash

# by Pengfei Liu on 7/6/21.

#all file with full path

#call ./bbmap_glacier_vOTU.sh mapping_metafile_TP85.txt 24


file1=$1 #file contains the reads and conitgs file info

threads=$2 #cpu used, 36

#check point
echo "${file1}" "${threads}"

IFS=$'\n' #separate by "\n" for cat

##### start of the mapping loop ###

#get the contigs files as ref

for line1 in $(cat "${file1}")
do

sample=$(echo "${line1}"| cut -f1 -d,) #for each sample

contigs=$(echo "${line1}"| cut -f2 -d,) #contigs file with full path

trimmedR1=$(echo "${line1}"| cut -f3 -d,) #trimmed read1 with full path

trimmedR2=$(echo "${line1}"| cut -f4 -d,) #trimmed read2 with full path


#all mapping files here
#one contigs, one dir

#outdir=./"${sample}"_bbmap_binning

if test -f "${sample}".bam
then
echo "file bam exits"

else
echo "bam not present"

#do mapping
/home/ptpe/software/bbmap/bbmap.sh ref="${contigs}" in="${trimmedR1}" in2="${trimmedR2}" xmtag=t ambiguous=random outm="${sample}".bam threads="${threads}" -Xmx100g
fi

done

#end of loop

echo "finish all"

unset IFS



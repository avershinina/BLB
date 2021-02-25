#!/bin/bash

###############################################
# usage: ./DFOIL_runner.sh INT
# Where INT is a number usually {1..n}, accordingly to the number of chromosomes\contigs in the reference genome 
# Dependencies: 
# D-foil
# seqkit
###############################################
# This script was created for horses, thus it can be run as 
# for N in {1..31}; do ./Dfoil_runner.sh $N; done > dfoil_log.txt
# This script runs Dfoil analysis on pseudohaploidized genomes (FASTA files).
# It is also possible to run Dfoil on VCF files, but this pipeline will not work for this purpose.

# Since the analysis is computationally intensive, during the course of this pipeline each genome is splitted into chromosomes and then into windows. This allows faster runs, but creates A LOT of files along the way. Make sure to clean up, remove, or gzip all these numerous files after the analysis is done. 

# If the genome of interest has a numerous number of contigs, and not just ~30 chromosomes, this pipeline will work. However, it will create *unimaginable* amount of files. So it is not ideal. Well, noone is ideal, so...    
# ##############################################

# Before running the pipeline, the following steps are requiered: 

# 1) create a pseudohaploidized FASTA for each of the genomes (with ANGSD -doFasta, ChromCompare pu2fa, or other similar tool);

# 2) split each genome by chromosome
# for f in *.fa.gz; do zcat $f | awk '/^>chr/ {OUT=substr($0,2) "'${f}'.fa"}; OUT {print >OUT}' & done

# 3) Rename each fasta header in a file according to the filename
# for f in chr{1..31}*.fa; do perl -p -e  "s/^>/>${f}-/g" $f > ${f}.renamed & done

# 4) # Simplify each fasta header: split into a string and take chr name and sample name from the string. Change the splitting arguments (_ and .) into whatever works. Rename the fasta header accordingly: 
# for f in chr{1..31}*.renamed; do cat $f | awk ' $1 ~ /^>/ { split($0,a,FS="[_,.]"); print a[1]"_"a[4];next} {print}' > ${f}.header; done

# 5) rename 's/fa.renamed.header/fasta/' *.header

# 6) Use bedtools to create XXXkb windows based on genome size and contig information. I used 200kb, use whatever works, like 100k or 200k.

# 7) Split the file with window coordinates into one file per chromosome. Use awk to create change '-' into ':' for the window coordinates.

# 8) Seqtk can't into 0-based coordinates, so add one for the start coordinate in each file:
# for f in chr*.windows.bed; do sed '0,/0:/{s/0:/1:/}' $f > $f.1fix & done

# The final bed file should look like this:
# chr10 1:200000
# chr10 200000:400000
# chr10 400000:600000
# chr10 600000:800000
# chr10 800000:1000000

# 9) Look at hard-coded variables in this script and modify accordinly to the needs (for example names of the genome files and such).

# ##############################################
# FINALLY
# This script is not ideal in any way, there are ways to optimize it. 
# Feel free to modify and improve. 
# ##############################################


#### Sort out the files by chromosome

mkdir chr${1}_outdir
CHR_OUTDIR=chr${1}_outdir

mv chr${1}_*dfoil_rmTrans_*.fasta ${CHR_OUTDIR}

echo "Files moved to "${CHR_OUTDIR}

CHR_COORDS=/200kb_coords/chr${1}.bed.1fix

echo "Processing "${CHR_COORDS}

# 
# Combine fastas into a form of Dfoil pops: P1,P2,P3,P4,O
# In this case P1-P3, O are fixed, P4 is variable. 

P1=genome1
P2=genome2
P3=genome3
P41=genome4
P42=genome5
P43=genome6
P44=genome7
P5=genome8

for POP in $P41 $P42 $P43 $P44; do
  cat ${CHR_OUTDIR}/chr${1}_dfoil_rmTrans_${P1}.fasta ${CHR_OUTDIR}/chr${1}_dfoil_rmTrans_${P2}.fasta ${CHR_OUTDIR}/chr${1}_dfoil_rmTrans_${P3}.fasta ${CHR_OUTDIR}/chr${1}_dfoil_rmTrans_${POP}.fasta ${CHR_OUTDIR}/chr${1}_dfoil_rmTrans_${P5}.fasta > ${CHR_OUTDIR}/chr${1}_P4_${POP}.fasta ;
done

wait
echo "Generating chromosome-wide alignments"
wait

for POP in $P41 $P42 $P43 $P44; do mkdir ${CHR_OUTDIR}/chr${1}_P4_${POP}_dir ; done

for POP in $P41 $P42 $P43 $P44; do mv ${CHR_OUTDIR}/chr${1}_P4_${POP}.fasta ${CHR_OUTDIR}/chr${1}_P4_${POP}_dir ; done 

for POP in $P41 $P42 $P43 $P44; do echo "Alignments generated and moved to ${CHR_OUTDIR}/chr${1}_P4_${POP}_dir" ; done

echo "----------------------------------"
echo "Beware, the nested loop is coming!"
echo "----------------------------------"


# Split each chromosome-wide alignment into windows  

for POP in $P41 $P42 $P43 $P44; do
  for W in $(cat ${CHR_COORDS} | awk '{print $2}'); do
    seqkit subseq -r $W ${CHR_OUTDIR}/chr${1}_P4_${POP}_dir/chr${1}_P4_${POP}.fasta > ${CHR_OUTDIR}/chr${1}_P4_${POP}_dir/chr${1}_P4_${POP}_${W}.fasta ; 
  done;
done

echo "Starting D-foil estimates"
# c=0; for f in $(grep ">" test.cat.fa); do TEST[$c]=$f; c=$c+1; done
# c=0; for f in $(grep ">" test.cat.fa | sed "s/^>//g"); do TEST[$c]=$f; c=$c+1; done

# To check - do echo ${TEST[0]}

# Edit here to change where Dfoil scripts are stored

mkdir Dfoil_results 
DFOIL_OUT=/Dfoil_results
DFOIL_SCRIPTS=/dfoil

# Loop through concatenated files to get fasta headers, that we then supply into the dfoil sample name flag (-n).

for POP in $P41 $P42 $P43 $P44; do
  c=0; # we are making an array, so need a counter. Each array element will have an assigned ID accordingly to this counter. 
  for f in $(grep ">" ${CHR_OUTDIR}/chr${1}_P4_${POP}_dir/chr${1}_P4_${POP}.fasta | sed "s/^>//g"); do 
    SAMPLE[$c]=$f; 
  c=$c+1; # Reset the counter to assign the next value to the $f variable (each $f variable is a hasta header without '>')
  done
  echo "Processing P4 population $POP sample P1 ${SAMPLE[0]}"
  echo "Processing P4 population $POP sample P2 ${SAMPLE[1]}"
  echo "Processing P4 population $POP sample P3 ${SAMPLE[2]}"
  echo "Processing P4 population $POP sample P4 ${SAMPLE[3]}"
  echo "Processing P4 population $POP sample P5 ${SAMPLE[4]}"

  wait

  python3.8 ${DFOIL_SCRIPTS}/fasta2dfoil.py ${CHR_OUTDIR}/chr${1}_P4_${POP}_dir/chr${1}_P4_${POP}_*.fasta -o ${DFOIL_OUT}/chr${1}_dfoil_rmTrans_P4_${POP}.Dcounts -n ${SAMPLE[0]},${SAMPLE[1]},${SAMPLE[2]},${SAMPLE[3]},${SAMPLE[4]}
  wait

  python3.8 ${DFOIL_SCRIPTS}/dfoil.py --infile ${DFOIL_OUT}/chr${1}_dfoil_rmTrans_P4_${POP}.Dcounts --out ${DFOIL_OUT}/chr${1}_dfoil_rmTrans_P4_${POP}.Dcounts.stats
  wait

  python3.8 ${DFOIL_SCRIPTS}/dfoil_analyze.py ${DFOIL_OUT}/chr${1}_dfoil_rmTrans_P4_${POP}.Dcounts.stats > ${DFOIL_OUT}/chr${1}_dfoil_rmTrans_P4_${POP}.Dcounts.stats.result
  wait
done

echo "All done"

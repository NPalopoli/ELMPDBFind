#!/bin/bash

if [ ! -d ./ELMPDBFind_files ]
then
  mkdir ./ELMPDBFind_files
fi

if [ ! -f ELMInstances.lst ]
then
  ls /home/npalopoli/SLiMBench/ELMmap/getELMInstancesSeqs_output >ELMInstances.lst
  sed -i 's/.fasta//g' ELMInstances.lst
fi

while read line
do
#  ./ELMPDBFind.py /home/npalopoli/SLiMBench/ELMmap/getELMInstancesSeqs_output/"$line".fasta /home/npalopoli/DBs/SEQATOMs/SEQATOMs_split_all.fasta >./ELMPDBFind_files/ELMPDBFind_"$line".tsv
  ./ELMPDBFind.py /home/npalopoli/SLiMBench/ELMmap/getELMInstancesSeqs_output/"$line".fasta /home/npalopoli/PDBss/PDB_ss_dis_SEQATOM_all.txt >./ELMPDBFind_files/ELMPDBFind_"$line".tsv
done<ELMInstances.lst 1>ELMPDBFind.out 2>ELMPDBFind.err

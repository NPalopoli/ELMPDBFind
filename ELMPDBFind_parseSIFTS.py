#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
  Parse pdb_chain_uniprot.tsv from SIFTS to map ELM entries to PDB ss and disorder information
  File name: parseSIFTS_PDBssSEQATOM.py
  Author: Nicolas Palopoli
  Date created: 2015/10/30
  Date last modified: 2015/11/02
  Python Version: 2.7
'''

import sys
import csv
from collections import OrderedDict
from Bio import SeqIO

# Read input files

try:
  inELMinstances = open(sys.argv[1])
  inSIFTSpdbuniprot = open(sys.argv[2])
  inELMPDBFind = open(sys.argv[3])
except IndexError:
  print("Input file not specified. Format: ./parseSIFTS_PDBssSEQATOM.py <elm_instances[.date].tsv> <pdb_chain_uniprot.tsv> <ELMPDBFind_<ELMInstance>.tsv>")
  exit()
except IOError:
  print("Input file not found. Format: ./parseSIFTS_PDBssSEQATOM.py <elm_instances[.date].tsv> <pdb_chain_uniprot.tsv> <ELMPDBFind_<ELMInstance>.tsv>")
  exit()

# Functions

def readELMPDBFind(infile):
  '''Store list of PDB where ELM instance is found'''
  listpdbfind = csv.DictReader(filter(lambda row: row[0]!='#', infile), delimiter='\t', quotechar='"', fieldnames=("Accession","PDB","CHAIN","SEQRESELM","SEQRESFLANK","SEQATOMELM","SEQATOMFLANK"))
  return listpdbfind
#  listfile = csv.DictReader(filter(lambda row: row[0]!='#', infile), delimiter='\t', quotechar='"')
#  return listfile

'''
npalopoli@antares:~/20150924_ELM-Struct/ELMPDBFind$ ls ELMPDBFind_files/*ELMI000004*
ELMPDBFind_files/ELMPDBFind_ELMI000004.tsv
npalopoli@antares:~/20150924_ELM-Struct/ELMPDBFind$ head ELMPDBFind_files/ELMPDBFind_ELMI000004.tsv
#NAM:sp|P97764|WBP1_MOUSE|ELMI000004|LIG|LIG_WW_1
#SEQ:MARASSRNSSEEAWGSLQAPQQQQSPAASSLEGAIWRRAGTQTRALDTILYHPQQSHLLRELCPGVNTQPYLCETGHCCGETGCCTYYYELWWFWLLWTVLILFSCCCAFRHRRAKLRLQQQQRQREINLLAYHGACHGAGPVPTGSLLDLRLLSAFKPPAYEDVVHHPGTPPPPYTVGPGYPWTTSSECTRCSSESSCSAHLEGTNVEGVSSQQSALPHQEGEPRAGLSPVHIPPSCRYRRLTGDSGIELCPCPDSSEGEPLKEARASASQPDLEDHSPCALPPDSVSQVPPMGLASSCGDIP
#ELM:PPPY
#FLK:VVHHPGTPPPPYTVGPGYPW
ELMI000004	3FBI	C	17	-	17	-
ELMI000004	3FBI	A	17	-	17	-
ELMI000004	1V55	A	499	-	499	-
ELMI000004	1V55	N	499	-	499	-
ELMI000004	4V4B	BM	108	-	108	-
ELMI000004	1I8L	D	101	-	101	-
'''

def parseELMPDBFind(listpdbfind):
  '''Make dict with ELMPDBfind information'''
  pdbfindpos = {}
  for row in listpdbfind:
    key = '%s:%s:%s' % (row['Accession'],row['PDB'],row['CHAIN'].upper())  # index by ELMAccession:PDB:CHAIN
    pdbfindpos[key] = [row['Accession'],row['PDB'],row['CHAIN'],row["SEQRESELM"],row["SEQRESFLANK"],row["SEQATOMELM"],row["SEQATOMFLANK"]]
  return pdbfindpos

def mapELMPDBFind(pdbfindpos, ELMpos):
  '''Get Primary_Acc for PDBs where ELM instance is found'''
  for key in pdbfindpos:
    if pdbfindpos[key][0] in ELMpos:
#    if key in ELMpos:  # match ELMAccession
      selectELMpos = [ELMpos[pdbfindpos[key][0]][0], ELMpos[pdbfindpos[key][0]][1], ELMpos[pdbfindpos[key][0]][2]]
      pdbfindpos[key].extend(selectELMpos)
#      pdbfindpos[key].append(ELMpos[pdbfindpos[key][0]][0])
  return pdbfindpos

def printELMPDBFind(pdbfindpos):
  '''Format ELMPDBFind with Primary_Acc for printing'''
  for key in pdbfindpos:
    print '\t'.join(map(str, pdbfindpos[key]))

def readELMinstances(infile):
  '''Store ELM instances information as list of dicts'''
  listfile = csv.DictReader(filter(lambda row: row[0]!='#', infile),delimiter='\t', quotechar='"')
  return listfile

def mapELMpositions(listELM):
  '''Make dict with {ELMAccession:[Primary_acc,Start,End,[PDB1,...,PDBN]]}'''
  ELMpos = {}
  for row in listELM:
    ELMpos[row['Accession']] = [row['Primary_Acc'],row['Start'],row['End'],row['PDB'].split()]
  return ELMpos

def readSIFTSpdbuniprot(infile):
  '''Store SIFTS information as list of dicts'''
  listfile = csv.DictReader(filter(lambda row: row[0]!='#', infile),delimiter='\t', quotechar='"')
  return listfile

def mapSIFTSpdbpositions(listSIFTS):
  '''Make dict with {'SP_PRIMARY:PDB':[PDB,CHAIN,RES_BEG,RES_END,PDB_BEG,PDB_END,SP_PRIMARY,SP_BEG,SP_END]}'''
  SIFTSpos = {}
  for row in listSIFTS:
    key = '%s:%s' % (row['SP_PRIMARY'],row['PDB'].upper())  # index by SP_PRIMARY:PDB
    SIFTSpos[key] = [row['PDB'].upper(),row['CHAIN'].upper(),row['RES_BEG'],row['RES_END'],row['PDB_BEG'],row['PDB_END'],row['SP_PRIMARY'],row['SP_BEG'],row['SP_END']]
  return SIFTSpos

def mapELMPDBFindposSIFTS(pdbfindpos,SIFTSpdbpos):
  '''Map ELMPDBFind start to Uniprot'''
  for key in pdbfindpos:
    target = '{}:{}'.format(pdbfindpos[key][7],pdbfindpos[key][1])  # index by Accession:PDB
    if target in SIFTSpdbpos:
      continue  # program this!!
   
'''
ELMPDBFind_files/ELMPDBFind_ELMI000004.tsv
npalopoli@antares:~/20150924_ELM-Struct/ELMPDBFind$ head ELMPDBFind_files/ELMPDBFind_ELMI000004.tsv
#NAM:sp|P97764|WBP1_MOUSE|ELMI000004|LIG|LIG_WW_1
#SEQ:MARASSRNSSEEAWGSLQAPQQQQSPAASSLEGAIWRRAGTQTRALDTILYHPQQSHLLRELCPGVNTQPYLCETGHCCGETGCCTYYYELWWFWLLWTVLILFSCCCAFRHRRAKLRLQQQQRQREINLLAYHGACHGAGPVPTGSLLDLRLLSAFKPPAYEDVVHHPGTPPPPYTVGPGYPWTTSSECTRCSSESSCSAHLEGTNVEGVSSQQSALPHQEGEPRAGLSPVHIPPSCRYRRLTGDSGIELCPCPDSSEGEPLKEARASASQPDLEDHSPCALPPDSVSQVPPMGLASSCGDIP
#ELM:PPPY
#FLK:VVHHPGTPPPPYTVGPGYPW
ELMI000004      3FBI    C       17      -       17      -
ELMI000004      3FBI    A       17      -       17      -
ELMI000004      1V55    A       499     -       499     -
'''

'''
PDB	CHAIN	SP_PRIMARY	RES_BEG	RES_END	PDB_BEG	PDB_END	SP_BEG	SP_END
101m	A	P02185	1	154	0	153	1	154
102l	A	P00720	1	40	1	40	1	40
102l	A	P00720	42	165	41	164	41	164
102m	A	P02185	1	154	0	153	1	154
103l	A	P00720	1	40	1	40	1	40
103l	A	P00720	44	167	41	164	41	164
103m	A	P02185	1	154	0	153	1	154
104l	A	P00720	1	44	1	44	1	44
'''


def getELMPDBFindpos(pdbfindpos,SIFTSpdbpos):
  '''Add Uniprot accession to ELMPDBFind list '''
  for key in pdbfindpos:
    target = '{}:{}'.format(pdbfindpos[key][7],pdbfindpos[key][1])  # index by Accession:PDB
    if target in SIFTSpdbpos:
      if SIFTSpdbpos[target][0] == pdbfindpos[key][1] and SIFTSpdbpos[target][1] == pdbfindpos[key][2] and SIFTSpdbpos[target][6] == pdbfindpos[key][7] and SIFTSpdbpos[target][7] <= pdbfindpos[key][4] and SIFTSpdbpos[target][8] >= pdbfindpos[key][5]:  # PDB, Chain and Primary Accession match, within limits
#      if SIFTSpdbpos[target][0] == pdbfindpos[key][1] and SIFTSpdbpos[target][1] == pdbfindpos[key][2]:  # PDB, Chain and Primary Accession match, within limits
        pdbfindpos[key].append(SIFTSpdbpos[target][6])
#    pdbfindpos[key].append(SIFTSpdbpos)
  return pdbfindpos

def getPDBpos(ELMpos,SIFTSpdbpos):
  ''' '''
  for keyELM,valELM in ELMpos.items():
    for PDB in valELM[3]:
      target = '%s:%s' % (valELM[0],PDB)  # index by ELMAccession:PDB
      if target in SIFTSpdbpos:  # subset of ELMs where PDB sequence is included in UniprotID -- need to extend this for cases where PDB is not in UniprotID
        if SIFTSpdbpos[target][0] == PDB and SIFTSpdbpos[target][6] == valELM[0] and SIFTSpdbpos[target][7] <= valELM[1] and SIFTSpdbpos[target][8] >= valELM[2]:  # PDB and Primary Accession match, within limits
          print parsePDBss('../PDBss/PDB_ss_dis_SEQATOM_all.txt',SIFTSpdbpos[target][0],SIFTSpdbpos[target][1],SIFTSpdbpos[target][2],SIFTSpdbpos[target][6],SIFTSpdbpos[target][7],keyELM)
#  return PDBss

def parsePDBss(infile,pdb,chain,resbeg,uniprot,spbeg,keyELM):
  '''Read 2nd struct from PDB_ss_dis_SEQATOM_all.txt (unified from PDB files ss.txt and ss_dis.txt and with SEQRES and SEQATOM from SEQATOMs db) in FASTA format'''
  resbeg = int(resbeg)
  seqgaps = '-' * (int(spbeg)-1)
#  seqgaps = '-' * (int(spbeg)-1+resbeg-1)
#  seqgaps = '-' * (int(spbeg)+resbeg-1)
  seq = '-'
  secstr = '-'
  disorder = '-'
  seqres = '-'
  seqatom = '-'
  fastaseqs = SeqIO.parse(open(infile),'fasta')
  for fasta in fastaseqs:
    if fasta.id[0:4] == pdb and fasta.id[5] == chain:
      if 'sequence' in fasta.id:
        seq = "%s%s" % (seqgaps,fasta.seq.tostring()[resbeg-1:])
      elif 'secstr' in fasta.id:
        secstr = "%s%s" % (seqgaps,fasta.seq.tostring()[resbeg-1:])
      elif 'disorder' in fasta.id:
        disorder = "%s%s" % (seqgaps,fasta.seq.tostring()[resbeg-1:])
      elif 'SEQRES' in fasta.id:
        seqres = "%s%s" % (seqgaps,fasta.seq.tostring()[resbeg-1:])
      elif 'SEQATOM' in fasta.id:
        seqatom = "%s%s" % (seqgaps,fasta.seq.tostring()[resbeg-1:])
  name = "%s:%s:%s:%s" % (keyELM,uniprot,pdb,chain)
  res = ">%s:sequence\n%s\n>%s:secstr\n%s\n>%s:disorder\n%s\n>%s:SEQRES\n%s\n>%s:SEQATOM\n%s\n" % (name,seq,name,secstr,name,disorder,name,seqres,name,seqatom)
#  res = ">%s:%s\n%s\n>%s:%s\n%s\n>%s:%s\n%s\n" % (uniprot,pdb,seq,uniprot,pdb,secstr,uniprot,pdb,disorder)
  return res

# START

# Make list of one dict per ELM instance in input file
listELM = readELMinstances(inELMinstances)
inELMinstances.close()
# Get ELM start and end positions on UniProt Primary Accession, and get all PDBs annotated for the ELM
ELMpos = mapELMpositions(listELM)

listpdbfind = readELMPDBFind(inELMPDBFind)
pdbfindpos = parseELMPDBFind(listpdbfind)
pdbfindpos = mapELMPDBFind(pdbfindpos, ELMpos)

# Make list of one dict per SIFTS entry in input file
listSIFTS = readSIFTSpdbuniprot(inSIFTSpdbuniprot)
inSIFTSpdbuniprot.close()
# Make list of dicts with positions from UniProt Primary Accession to PDB Chain SEQRES
# e.g.: { (...), '3DNL:P04578': ['3DNL', 'E', '1', '67', '330', '396', 'P04578', '330', '396'], (...) }
SIFTSpdbpos = mapSIFTSpdbpositions(listSIFTS)

pdbfindpos = getELMPDBFindpos(pdbfindpos,SIFTSpdbpos)
printELMPDBFind(pdbfindpos)

exit()

# Make list of one dict per SIFTS entry in input file
parsedSIFTS = readSIFTSpdbuniprot(inSIFTSpdbuniprot)
inSIFTSpdbuniprot.close()

# Make list of dicts with positions from UniProt Primary Accession to PDB Chain SEQRES
# e.g.: { (...), '3DNL:P04578': ['3DNL', 'E', '1', '67', '330', '396', 'P04578', '330', '396'], (...) }
SIFTSpdbpos = mapSIFTSpdbpositions(parsedSIFTS)

# Get PDB, chain, start and end positions for entries in ELM
PDBss = getPDBpos(ELMpos,SIFTSpdbpos)
#print PDBss

exit()

#PDBss = parsePDBss('../PDBss/PDB_ss_dis.txt',SIFTSpdbpos['Q60795'][0].upper(),SIFTSpdbpos['Q60795'][1].upper())
#PDBss = parsePDBss('../PDBss/PDB_ss_dis_all.txt',SIFTSpdbpos[primaryAcc][0].upper(),SIFTSpdbpos[primaryAcc][1].upper())
PDBss = {}
for key in SIFTSpdbpos:
  print parsePDBss('../PDBss/PDB_ss_dis_all.txt',SIFTSpdbpos[key][0],SIFTSpdbpos[key][1])
  PDBss[key] = parsePDBss('../PDBss/PDB_ss_dis_all.txt',SIFTSpdbpos[key][0],SIFTSpdbpos[key][1])
print PDBss
# {'disorder': '-----------------------------------', 'seq': 'MDLIDILWRQDIDLGVSREVFDFSQRQKDYELEKQ', 'pdb': '3WN7', 'chain': 'M', 'secstr': '--HHHHHHTTTGGGT--GGGG--------------'}



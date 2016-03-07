#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
  Find ELM Instances in PDB
  File name: ELMPDBFind.py
  Author: Nicolas Palopoli
  Date created: 2016/01/24
  Date last modified: 2016/02/26
  Python Version: 2.7
'''

import sys
from collections import OrderedDict
import csv
from Bio import SeqIO

# Read input files
try:
  in_elm_fasta = open(sys.argv[1])
  in_seqatom_fasta = open(sys.argv[2])
except IndexError:
  print("Input file(s) not specified. Format: ./ELMPDBFind.py <in_ELM-instance-flank.fasta> <in_PDB-SEQATOM.fasta>")
  # ./ELMPDBFind.py /home/npalopoli/SLiMBench/ELMmap/getELMInstancesSeqs_output/"$line".fasta /home/npalopoli/DBs/SEQATOMs/SEQATOMs_split_all.fasta
  exit()
except IOError:
  print("Input file(s) not found. Format: ./ELMPDBFind.py <in_ELM-instance-flank.fasta> <in_PDB_SEQATOM.fasta>")
  exit()

# Functions

def read_motif(in_elm_fasta):
  '''Store sequences of ELM instance and flanking residues'''
  records = SeqIO.parse(in_elm_fasta, "fasta")
  for record in records:
    record.seq.strip("-")
  return records
  
def read_fasta(infasta):
  '''Store fasta sequences from file.'''
  seqs = OrderedDict()
#  seqs={}  # dict for raw seqs
  readFirstSeq = False
  for line in infasta:
    if line.strip() or line not in ['\n', '\r\n']:  # avoid empty or only whitespace lines
      line=line.rstrip()  # discard newline at the end (if any)
      if line[0]=='>':  # or line.startswith('>'); distinguish header
        if readFirstSeq:  # exit if more than 2 sequences
          break
        readFirstSeq = True
        words=line.split()
        name=words[0][1:]
        seqs['res']=''
      else :  # sequence, not header, possibly multi-line
        seqs['res'] = seqs['res'] + line
  seqs['res'] = list(seqs['res'])
  seqs['position'] = range(1,len(seqs['res'])+1)
  seqs['name'] = [name] * len(seqs['res'])
  return seqs

def readELMinstances(infile):
  '''Store ELM instances information as list of dicts'''
  elm = csv.DictReader(filter(lambda row: row[0]!='#', infile),delimiter='\t', quotechar='"')
  return elm

def mapELMpositions(parsedELM,primaryAcc):
  '''Make dict with [ELMAccession:[Start,End,ELMType,ELMIdentifier,Primary_Acc]]'''
  ELMpos = {}
  for row in parsedELM:
    if primaryAcc == row['Primary_Acc']:
      ELMpos[row['Accession']] = [row['Start'],row['End'],row['ELMType'],row['ELMIdentifier'],row['Primary_Acc']]
  return ELMpos

def placeELM(seq,ELMpos):
  '''Map ELM to fasta sequence'''
  seqlen = len(seq['res'])
  seq['ELMpos'] = list('-' * seqlen)
  seq['ELMacc'] = list('-' * seqlen)
  seq['ELMType'] = list('-' * seqlen)
  seq['ELMIdentifier'] = list('-' * seqlen)
  seq['ELMflank'] = list('-' * seqlen)
#  for accession, limits in ELMpos.iteritems():
#    for pos in range(int(limits[0])-1,int(limits[1])):
  for accession, vals in ELMpos.iteritems():
    for pos in range(int(vals[0])-1,int(vals[1])):
      seq['ELMpos'][pos] = seq['res'][pos]
      if '-' in seq['ELMacc'][pos]: 
        seq['ELMacc'][pos] = accession
        seq['ELMType'][pos] = vals[2]
        seq['ELMIdentifier'][pos] = vals[3]
      else:
        seq['ELMacc'][pos] = seq['ELMacc'][pos] + accession
        seq['ELMType'][pos] = seq['ELMType'][pos] + vals[2]
        seq['ELMIdentifier'][pos] = seq['ELMIdentifier'][pos] + vals[3]
#    flanksize = ( 20 + int(vals[0]) - int(vals[1])) / 2  # Compute flanking sequence
    flanksize = 20 + int(vals[0]) - 1 - int(vals[1])  # Compute flanking sequence
    flanksizeodd = 0 if flanksize % 2 == 0 else 1
    flanksize = flanksize / 2
#    flankfirst = max(0, int(vals[0]) - flanksize)
    flankfirst = max(1, int(vals[0]) - flanksize)
    flanklast = min(int(vals[1]) + flanksize + flanksizeodd, seqlen)
#    if (flankfirst == 0 and flanklast == seqlen):
    if (flankfirst == 1 and flanklast == seqlen):
      continue
#    elif (flankfirst == 0):
    elif (flankfirst == 1):
#      flanklast = flanklast + int(vals[0]) - 1
#      flanklast = flanklast + flanksize + flanksizeodd + int(vals[0]) - 1
      flanklast = flanklast + flanksize - ( int(vals[0]) - 1 )
    elif (flanklast == seqlen):
#      flankfirst = flankfirst - (flanksize - (flanklast - len(seq['res'])))
#      flankfirst = flankfirst - (flanklast - seqlen))
#      flankfirst = flankfirst - flanksize - (flanklast - seqlen)
      flankfirst = flankfirst - (int(vals[1]) + flanksize + flanksizeodd - seqlen)
    for pos in range(flankfirst-1,flanklast):
      seq['ELMflank'][pos] = seq['res'][pos]
  return seq

def OLD_find_elm_pdb(pdb_seqs, elm_seqs):
  '''Find ELM Instance and flanking sequences in PDB sequences'''
#  len_adaptor = len(adaptor) #cache this for later
  pdb_match = {}
  for record in pdb_seqs:
    # if 'SEQRES' in pdb_seqs.id:
    recordname = str(record.id.rsplit('|', 1)[0])
    index = record.id.find('SEQRES')
    if index != -1:  # reading SEQRES sequence
      pdb_match[recordname] = {}
      instance_index = record.seq.find(elm_seqs['instance'])
      if instance_index == -1:  # if instance not found in SEQRES
       seqres_match = 0
       pdb_match[recordname]['instance-seqres'] = '-'
       pdb_match[recordname]['flanks-seqres'] = '-'
      else:  # if instance found in SEQRES
        seqres_match = 1
        pdb_match[recordname]['instance-seqres'] = instance_index + 1
        flanks_index = record.seq.find(elm_seqs['flanks'])
        if flanks_index == -1:  # if flanks not found in SEQRES
          pdb_match[recordname]['flanks-seqres'] = '-'
        else:  # # if flanks found in SEQRES
          pdb_match[recordname]['flanks-seqres'] = flanks_index + 1
    else:  # reading SEQATOM sequence
      if seqres_match != 1:  # if instance not found in corresponding SEQRES
        pdb_match[recordname]['instance-seqatom'] = '-'
        pdb_match[recordname]['flanks-seqatom'] = '-'
      else:  # if instance found in corresponding SEQRES
        instance_index = record.seq.find(elm_seqs['instance'])
        if instance_index == -1:  # if instance not found in SEQATOM
          pdb_match[recordname]['instance-seqatom'] = '-'
          pdb_match[recordname]['flanks-seqatom'] = '-'
        else:  # if instance found in SEQATOM
          pdb_match[recordname]['instance-seqatom'] = instance_index + 1
          flanks_index = record.seq.find(elm_seqs['flanks'])
          if flanks_index == -1:  # if flanks not found in SEQATOM
            pdb_match[recordname]['flanks-seqatom'] = '-'
          else:  # if flanks found in SEQATOM
            pdb_match[recordname]['flanks-seqatom'] = flanks_index + 1
      seqres_match = 0
  return pdb_match  

def find_elm_pdb(pdb_seqs, elm_seqs, elm_seqs_key):
  '''Find ELM Instance and flanking sequences in PDB sequences'''
#  len_adaptor = len(adaptor) #cache this for later
  pdb_match = {}
  elm_seqs_key_sequence = '-'.join([elm_seqs_key,'sequence'])
  elm_seqs_key_seqres = '-'.join([elm_seqs_key,'seqres'])
  elm_seqs_key_seqatom = '-'.join([elm_seqs_key,'seqatom'])
  sequence_match = 0
  seqres_match = 0

  for record in pdb_seqs:
    # if 'SEQRES' in pdb_seqs.id:
    recordname = str(record.id.rsplit('|', 1)[0])
    indexseqres = record.id.find('SEQRES')
    indexseqatom = record.id.find('SEQATOM')
#    indexsequence = record.id.find('sequence')
#    if indexsequence != -1:  # reading PDB_ss.txt sequence; activate this and indexsequence definition above if line below fails
    if record.id.find('sequence') != -1:  # reading PDB_ss.txt sequence (identical to Uniprot?)
      pdb_match[recordname] = {}

      instance_index = record.seq.find(elm_seqs['instance'])
      if instance_index == -1:  # if instance not found in PDB_ss sequence
        sequence_match = 0
        pdb_match[recordname]['instance-sequence'] = '-'
        pdb_match[recordname]['flanks-sequence'] = '-'
        pdb_match[recordname][elm_seqs_key_sequence] = '-'
      else: # if instance found in PDB_ss sequence
        sequence_match = 1
        pdb_match[recordname]['instance-sequence'] = instance_index + 1

        flanks_index = record.seq.find(elm_seqs['flanks'])
        if flanks_index == -1:  # if flanks not found in PDB_ss sequence
          pdb_match[recordname]['flanks-sequence'] = '-'
          pdb_match[recordname][elm_seqs_key_sequence] = '-'
        else:  # if flanks found in PDB_ss sequence
          pdb_match[recordname]['flanks-sequence'] = flanks_index + 1

        flankn_index = record.seq.find(elm_seqs[elm_seqs_key])
        if flankn_index == -1:  # if flankN not found in PDB_ss sequence
          pdb_match[recordname][elm_seqs_key_sequence] = '-'
        else:  # if flankN found in PDB_ss sequence
          pdb_match[recordname][elm_seqs_key_sequence] = flankn_index + 1

    elif indexseqres != -1:  # reading SEQRES sequence
      if sequence_match != 1:  # if instance not found in corresponding PDB_ss sequence
        pdb_match[recordname]['instance-seqres'] = '-'
        pdb_match[recordname]['flanks-seqres'] = '-'
        pdb_match[recordname][elm_seqs_key_seqres] = '-'

      else:  # if instance found in corresponding PDB_ss sequence
        instance_index = record.seq.find(elm_seqs['instance'])
        if instance_index == -1:  # if instance not found in SEQRES
         seqres_match = 0
         pdb_match[recordname]['instance-seqres'] = '-'
         pdb_match[recordname]['flanks-seqres'] = '-'
         pdb_match[recordname][elm_seqs_key_seqres] = '-'
        else:  # if instance found in SEQRES
          seqres_match = 1
          pdb_match[recordname]['instance-seqres'] = instance_index + 1

          flanks_index = record.seq.find(elm_seqs['flanks'])
          if flanks_index == -1:  # if flanks not found in SEQRES
            pdb_match[recordname]['flanks-seqres'] = '-'
            pdb_match[recordname][elm_seqs_key_seqres] = '-'
          else:  # if flanks found in SEQRES
            pdb_match[recordname]['flanks-seqres'] = flanks_index + 1
        
          flankn_index = record.seq.find(elm_seqs[elm_seqs_key])
          if flankn_index == -1:  # if flankN not found in SEQRES
            pdb_match[recordname][elm_seqs_key_seqres] = '-'
          else:  # if flankN found in SEQRES
            pdb_match[recordname][elm_seqs_key_seqres] = flankn_index + 1

    elif indexseqatom != -1:  # reading SEQATOM sequence
      if seqres_match != 1:  # if instance not found in corresponding SEQRES
        pdb_match[recordname]['instance-seqatom'] = '-'
        pdb_match[recordname]['flanks-seqatom'] = '-'
        pdb_match[recordname][elm_seqs_key_seqatom] = '-'

      else:  # if instance found in corresponding SEQRES
        instance_index = record.seq.find(elm_seqs['instance'])
        if instance_index == -1:  # if instance not found in SEQATOM
          pdb_match[recordname]['instance-seqatom'] = '-'
          pdb_match[recordname]['flanks-seqatom'] = '-'
          pdb_match[recordname][elm_seqs_key_seqatom] = '-'
        else:  # if instance found in SEQATOM
          pdb_match[recordname]['instance-seqatom'] = instance_index + 1

          flanks_index = record.seq.find(elm_seqs['flanks'])
          if flanks_index == -1:  # if flanks not found in SEQATOM
            pdb_match[recordname]['flanks-seqatom'] = '-'
            pdb_match[recordname][elm_seqs_key_seqatom] = '-'
          else:  # if flanks found in SEQATOM
            pdb_match[recordname]['flanks-seqatom'] = flanks_index + 1

          flankn_index = record.seq.find(elm_seqs[elm_seqs_key])
          if flankn_index == -1:  # if flankN not found in SEQATOM
            pdb_match[recordname][elm_seqs_key_seqatom] = '-'
          else:  # if flankN found in SEQATOM
            pdb_match[recordname][elm_seqs_key_seqatom] = flankn_index + 1

      sequence_match = 0
      seqres_match = 0

  return pdb_match


# Start      

# Make dict of sequences in input ELM file
elm_seqs = {}
for index, record in enumerate(SeqIO.parse(in_elm_fasta, "fasta")):
  if 'sequence' in record.id:
    elm_seqs['name'] = record.id.rsplit('|',1)[0]
    elm_seqs['sequence'] = record.seq
    print '#NAM:{}\n#SEQ:{}'.format(elm_seqs['name'],elm_seqs['sequence'])
  elif 'instance' in record.id:
    elm_seqs['instance'] = record.seq.strip("-")
    print '#ELM:{}'.format(elm_seqs['instance'])
  elif 'flanks' in record.id:
    elm_seqs['flanks'] = record.seq.strip("-")
    print '#FLK:{}'.format(elm_seqs['flanks'])
  else:  # parse flankn
    elm_seqs_key = record.id.rsplit('|',1)[1]  # e.g. flank10
    elm_seqs[elm_seqs_key] = record.seq.strip("-")
    elm_seqs_flanklen = len(elm_seqs[elm_seqs_key]) - len(elm_seqs['instance'])
    print '#F{}:{}'.format(elm_seqs_flanklen,elm_seqs[elm_seqs_key])
    print '{}\t{}\t{}\t{}\t{}\t{}{}\t{}\t{}\t{}{}\t{}\t{}\t{}{}'.format('ELMAccession','PDBID','PDBChain','ELM-SeqPDBss','flanks-SeqPDBss',elm_seqs_key,'-SeqPDBss','ELM-SEQRES','flanks-SEQRES',elm_seqs_key,'-SEQRES','ELM-SEQATOM','flanks-SEQATOM',elm_seqs_key,'-SEQATOM')
#    print '#FLN:{}'.format(elm_seqs[elm_seqs_key])
#print '#NAM:{}\n#SEQ:{}\n#ELM:{}\n#FLK:{}\n#FLN:{}'.format(elm_seqs['name'],elm_seqs['sequence'],elm_seqs['instance'],elm_seqs['flanks'],elm_seqs[elm_seqs_key])

# Find ELM sequences in PDB sequences
pdb_seqs = SeqIO.parse(in_seqatom_fasta, "fasta")
#pdb_match = find_elm_pdb(pdb_seqs, elm_seqs)
elm_seqs_key_sequence = '-'.join([elm_seqs_key,'sequence'])
elm_seqs_key_seqres = '-'.join([elm_seqs_key,'seqres'])
elm_seqs_key_seqatom = '-'.join([elm_seqs_key,'seqatom'])
pdb_match = find_elm_pdb(pdb_seqs, elm_seqs, elm_seqs_key)
#print 'pdb_match:\n{}'.format(pdb_match)
#for pdbchain, match in pdb_match:
for pdbchain, match in pdb_match.iteritems():
#  for matchkey, matchval in match.iteritems():
#    if '-' not in matchkey.values():
#  for matchkey in match.iterkeys():
#    if matchkey.['instance-seqres'] != '-':

#  if match['instance-seqres'] != '-':
#    print '{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(elm_seqs['name'].split('|')[3],pdbchain.split('|')[0],pdbchain.split('|')[1],match['instance-seqres'],match['flanks-seqres'],match['instance-seqatom'],match['flanks-seqatom'])
  if match['instance-sequence'] != '-':
    print '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(elm_seqs['name'].split('|')[3],pdbchain.split('|')[0],pdbchain.split('|')[1],match['instance-sequence'],match['flanks-sequence'],match[elm_seqs_key_sequence],match['instance-seqres'],match['flanks-seqres'],match[elm_seqs_key_seqres],match['instance-seqatom'],match['flanks-seqatom'],match[elm_seqs_key_seqatom])
exit()



count = SeqIO.write(pdb_match, "pdb_match.fasta", "fasta")
print("Saved %i reads" % count)

exit()

parsedELM = readELMinstances(inELMinstances)
inELMinstances.close()
ELMpos = mapELMpositions(parsedELM,primaryAcc)
'''
print ELMpos
seq = placeELM(seq,ELMpos)

for instance in ELMpos:
  header = '|'.join([seq['name'][0],instance,ELMpos[instance][2],ELMpos[instance][3]])
print '{}{}|sequence'.format('>',''.join(header))
print ''.join(seq['res'])
#print '{}{}|{}|{}|{}|ELMInstance'.format('>',''.join(seq['name'][0]),ELMpos.keys(),ELMpos['ELMType'],ELMpos['ELMIdentifier'])
print ''.join(seq['ELMpos'])
'''
for instance in ELMpos:
  singleELMpos = {}
  singleELMpos[instance] = ELMpos[instance]
  seq = placeELM(seq,singleELMpos)
  header = '|'.join([seq['name'][0],instance,ELMpos[instance][2],ELMpos[instance][3]])
  outname = instance + '.fasta'
  outfile=open(outname, 'w+')
  print >>outfile, '{}{}|sequence'.format('>',''.join(header))
  print >>outfile, ''.join(seq['res'])
  print >>outfile, '{}{}|instance'.format('>',''.join(header))
  print >>outfile, ''.join(seq['ELMpos'])
  print >>outfile, '{}{}|flanks'.format('>',''.join(header))
  print >>outfile, ''.join(seq['ELMflank'])
  outfile.close()

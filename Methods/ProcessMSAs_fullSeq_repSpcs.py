#!/usr/bin/python3
#
# This program pre-check the input MSA before futher coevolution calculations.
# The first record in MSA must be query sequence.
# Choose the seq that has Max(identy_query) as the representative homolog from one species.
# Max(identy_query) must be >0.1 to avoid mismatched domains.
# Select each MSA row that only corresponds to the query sequence and save as output.  
# ----------------------------------------------------------------------
# Written by Haiqing Zhao
# Honig Group at Columbia University
# Last Update: 02/01/2023
# ----------------------------------------------------------------------

import sys,os,subprocess
from Bio import AlignIO,SeqIO,Seq,SeqRecord
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import difflib

if len(sys.argv)< 3:
    print ("\nUsage: python ProcessMSA_repSpcs.py Output.msa Query_&_MSA.msa (Query_fasta2.msa) \n")
    print ("Default MSA format: fasta or stockholm")
    #print ("Supported MSA format: clustal, emboss, fasta, fasta-m10, ig, maf, mauve, nexus, phylip, phylip-sequential, phylip-relaxed, stockholm")
    exit()

if len(sys.argv) > 3:
    multi_MSA = True
else:
    multi_MSA = False

try:
    alignment1 = AlignIO.read(sys.argv[2], 'fasta')
    format='fasta'
except ValueError:
    try:
        alignment1 = AlignIO.read(sys.argv[2], 'stockholm')
        format='stockholm'
    except ValueError:
        print('Input hmmer.msa not in fasta or stockholm format.')
        exit()

if multi_MSA:
    for i in range(0,len(sys.argv)-3):
        alignment_i = AlignIO.read(sys.argv[i+3], 'fasta')
        alignment1.extend(alignment_i[1:])

def pick_subseq (long_seq,index_list):
    new_seq=""
    for i in index_list:
        new_seq = new_seq + long_seq[i]
    return Seq(new_seq)

query_inMSA = alignment1[0] #First record in MSA must be the query.

if format =='stockholm':
    if len(query_inMSA.id.split('|')) > 1:
        uni = query_inMSA.id.split('|')[1]
    elif len(query_inMSA.id.split('|')) == 1:
        uni = query_inMSA.id.split('|')[0]

    gapindex_list = [i for i, x in enumerate(query_inMSA) if x == "-"]
    nongap_origindex_list = set(range(0,len(query_inMSA))) - set(gapindex_list)
    nongap_origindex_list = sorted(list(nongap_origindex_list))
elif format =='fasta':
    uni = query_inMSA.description.split()[1]

if len(query_inMSA.description.split("OX="))>1:
    query_Tax = query_inMSA.description.split("OX=")[1].split(' ')[0]
else:
    query_Tax = 'unknown'

print("Total # of sequences/Uni/Query_length: ",len(alignment1),uni, len(query_inMSA.seq.ungap('-')))

"""
#--- Read aln; Pick subseq if stockholm (from jackhmmer) ---#
"""

species_fullseq, species_genenames = {},{}
for record in alignment1[1:]:
    if format =='stockholm':
        seq_full = pick_subseq(record.seq.upper(),nongap_origindex_list)
        [species,genename] = record.id.split('.')[0:2]
    elif format =='fasta':
        species = record.id; genename = record.description
        seq_full= record.seq

    if species in species_fullseq:
        species_fullseq[species].append(seq_full)
        species_genenames[species].append(genename)
    else:
        species_fullseq.update({species:[seq_full]})
        species_genenames.update({species:[genename]})

"""
#--- Choose 1 representative seq for 1 species ---#
"""
species_fullseq_RepSeq,species_genenames_RepSeq = {},{}
for species in species_fullseq:
    id_list=[]; 
    for seq in species_fullseq[species]:
        if format =='fasta':
            idt = difflib.SequenceMatcher(None, query_inMSA.seq.ungap('-'),seq,autojunk=False).ratio()
        elif format =='stockholm':
            idt = difflib.SequenceMatcher(None, query_inMSA.seq.ungap('-'),seq.ungap('-'),autojunk=False).ratio()
        id_list.append(idt)
    max_index = id_list.index(max(id_list))
    repSeq_species = species_fullseq[species][max_index] 
    repSeq_genename = species_genenames[species][max_index]
    if max(id_list) >= 0.1:
        species_fullseq_RepSeq.update({species:repSeq_species})
        species_genenames_RepSeq.update({species:repSeq_genename})

print("Total # of species: ",len(species_fullseq_RepSeq))
sorted_species = sorted([i for i in species_fullseq_RepSeq.keys()])

"""
#--- Write 1rep_seq of 1species into Fasta ---#
"""
outfile1 = sys.argv[1]
print("Writing to: ",outfile1)

query_record = SeqRecord(query_inMSA.seq.ungap("-"),id=query_Tax,description = uni+" ")

with open(outfile1, "w") as out:
    SeqIO.write(query_record, out, "fasta")
    for species in sorted_species:
        species = str(species)
        record1 = SeqRecord(species_fullseq_RepSeq[species],id=species,description = species_genenames_RepSeq[species])
        SeqIO.write(record1, out, "fasta")

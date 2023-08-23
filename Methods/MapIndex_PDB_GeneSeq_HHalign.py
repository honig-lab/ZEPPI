#!/usr/bin/python3
#
# Name: MapIndex_PDB_GeneSeq_HHalign.py
# This program reads a PDB file, output the seqmapping between PDB residue indices
# and indices from fasta file (gene sequence).
# Insertions in PDB are not considered.
# The hhalign (run from c2b2) from HHBlits aligns much better than pairwise2.local module.
# Output saves in /ifs/data/c2b2/bh_lab/hz2592/Bacterial_PDB_MI/PDB_IFC_MI/HH_output/
# ----------------------------------------------------------------------
# Author: Haiqing Zhao
# Honig Group at Columbia University
# Created: 06/20/2020
# Last Update: 02/01/2023
# Supports modified amino acids.
# Takes care of PDB insertion & missing cases
# Added namecheck_4bash for space issue in filename
# Added important feature to check index match or not
# Skip case like E65A, R65B 
# ----------------------------------------------------------------------

import os,sys,tempfile,subprocess
import Bio.PDB.PDBParser
from Bio import Seq,SeqIO,AlignIO
from Bio.PDB import Polypeptide as pp
from Bio.Data import SCOPData

#Command: python MapIndex_PDB_GeneSeq_HHalign.py PDB chainID Fastafile SeqMap_output

def is_regular_res(res):
    #res_id = res.get_id()[1]
    if res.get_id()[2] ==' ':
        res_id = res.get_id()[1]
    #else:
        #res_id = str(res.get_id()[1])+res.get_id()[2]
        #print(res_id)
    #if res_id.isnumeric():
        if pp.is_aa(res) and res.has_id('CA') and res.has_id('O') and res_id>0: # Modified amino acids are supported
            return True
        else:
            return False
    #return False

def res_fullid(res):
    res_One = SCOPData.protein_letters_3to1[res.resname] #pp.three_to_one(res.resname)
    if res.get_id()[2] ==' ':
        res_id = res.get_id()[1]
    else:
        res_id = str(res.get_id()[1])+res.get_id()[2]
        #print(res_id)
    res_fullid = res_One + str(res_id)
    return res_fullid

def sequence_obtain(chain):
    sequence=[]
    for res in chain:
        if is_regular_res(res):
            aminos = SCOPData.protein_letters_3to1[res.resname] #included modified aa.pp.three_to_one(res.resname)
            sequence.append(aminos)
    sequence = Seq.Seq("".join(sequence)) # Join with nothing, return Seq object.
    if len(sequence)==0:
        print("Warning: no amino acid sequence obtained in ",chain)
    return sequence #Seq object/class, not SeqRecord

def fullid_obtain(chain):
    fullid_list = []
    for res in chain:
        if is_regular_res(res):
            resi_fullid = res_fullid(res)
            fullid_list.append(resi_fullid)
    if len(fullid_list)==0:
        print("Warning: no amino acid fullid obtained in ",chain)
    return fullid_list

"""
"""
hh_dir = "/ifs/data/c2b2/bh_lab/hz2592/Bacterial_PDB_MI/PDB_IFC_MI/HH_output/" #tempfile.TemporaryDirectory()
hh_cmd = "/ifs/data/c2b2/bh_lab/shares/hhpred/bin/hhalign -i "
#hh_dir = "/Users/hz2592/OneDrive - cumc.columbia.edu/WORK/Ecoli_MI/Decoys/Dockground_Unbound/Unbound_Decoy_Set2/2grx/"
#hh_cmd = "hhalign -i "

def namecheck_4bash(file_name):
    return file_name.replace(' ', '\ ')

def call_hhalign(PDBseq_object1,GeneseqRecord_object2,filename): #Seq object vs. SeqRecord object
    if '/' not in filename:
        filename = hh_dir+filename
    tmp_seqfile1 = tempfile.NamedTemporaryFile(delete=False)
    tmp_seqfile2 = tempfile.NamedTemporaryFile(delete=False)

    with open(tmp_seqfile1.name, "w") as handle1, open(tmp_seqfile2.name, "w") as handle2:
        handle1.write(">PDB\n"+str(PDBseq_object1)+"\n")
        SeqIO.write(GeneseqRecord_object2,handle2,"fasta")
    
    filename_bash =namecheck_4bash(filename)
    cmd_hhalign = [hh_cmd+ tmp_seqfile2.name + " -t " +tmp_seqfile1.name +" -nocons -ofas " + filename_bash ] #for c2b2
#   cmd_hhalign = [hh_cmd+ tmp_seqfile2.name + " -t " +tmp_seqfile1.name +" -hide_cons -Ofas " + filename_bash ] #for iMac
    
    #/ifs/data/c2b2/bh_lab/shares/hhpred/bin/hhalign  -i OrthGOPHER/Uni_FASTA/Q6P503.fasta  -t 1tmp.fa -ofas test
    #cd "+ temp_dir.name +";
    
    process = subprocess.Popen(cmd_hhalign, shell=True)
    output, error = process.communicate()
    p_status = process.wait()
  
    # temp_dir.cleanup()
    # os.remove(tmp_seqfile1.name)
    # os.remove(tmp_seqfile2.name)
    if os.path.exists(filename) and os.path.getsize(filename) > 0:
        return True
    else:
        return False

### Main Function Starts Here ##
def main(PDBid,chainname,fastafile,output):
    if os.path.exists(output) and os.path.getsize(output) > 0:
        print("Seqmap already exists; redo ",output)
    else:
        pass
    out = open(output, "w")

### Read Gene seq and PDB seq   
    for record in SeqIO.parse(fastafile, "fasta"):
        geneseqRecord_full = record
        geneseq_full = record.seq

    if PDBid[-4:].lower()==".pdb":
        pdb_file = PDBid
    else:
        pdb_file = PDBid + ".pdb"
    p = Bio.PDB.PDBParser(QUIET = True)
    s = p.get_structure(pdb_file, pdb_file)
    for chain in s[0]:
        if chain.id == chainname:
            pdbseq_full = sequence_obtain(chain)
            pdb_fullid_list = fullid_obtain(chain)

#    print(geneseq_full)
#    for id in pdb_fullid_list:
        #print (id)
#        if id[0] == geneseq_full[int(id[1:])+7]:
#            print (id,id[0])

    if not pdbseq_full or not pdb_fullid_list:
        print("Cannot obtain regular amino acids. Exit function.")
        return
    if len(pdbseq_full) != len(pdb_fullid_list):
        print("Error: sequence read error",len(pdbseq_full),len(pdb_fullid_list) )
        return

#### Align Gene seq and PDB seq   
    hh_filename = hh_dir+"hh_"+pdb_file[-8:-4]+"_"+chainname
    hh_success = call_hhalign(pdbseq_full,geneseqRecord_full,hh_filename)

    if hh_success:
        alignment = AlignIO.read(open(hh_filename), "fasta")
    else:
        print("Cannot find the hhalign output file. Exit function.",hh_filename)
        return

#### Read Alignment   
    if 'PDB' in alignment[0].id: 
        pdbseq_aligned = alignment[0].seq
        geneseq_alignedpart = alignment[1].seq
    else:   
        geneseq_alignedpart = alignment[0].seq
        pdbseq_aligned = alignment[1].seq

    geneseq_alignedpart_nogap = geneseq_alignedpart.ungap("-")
    pdbseq_aligned_nogap = pdbseq_aligned.ungap("-")
#    print("geneseq_full",geneseq_full)
#    print("geneseq_alignedpart_nogap",geneseq_alignedpart_nogap)
    index_genepart_start = str(geneseq_full).find(str(geneseq_alignedpart_nogap))
    if index_genepart_start >=0:
        index_genepart_start = index_genepart_start+1 #real count from 1
    else:
        print("Error: aligned seq_substring not found in geneSeq. Exit.")
        return

    geneseq_alignedpart_nogap_fullgeneindex=[];n_gap_insert=0
    for idx, i in enumerate(geneseq_alignedpart): #idx from 0
        if i != "-":
            i_fullid = i + str(idx+index_genepart_start - n_gap_insert)
            geneseq_alignedpart_nogap_fullgeneindex.append(i_fullid)
        else:
            n_gap_insert = n_gap_insert + 1

    index_pdbpart_start = str(pdbseq_full).find(str(pdbseq_aligned_nogap)) #pdb starting position
    if index_pdbpart_start < 0:
        print("Error: aligned seq_substring not found in PdbSeq. Exit.")
        return

    n_gap_insert, n_gap_miss = 0,0
    chain_orig_index = index_pdbpart_start
    for hhalned_idx, i in enumerate(pdbseq_aligned): #start from 0
        if i != "-": # Skip case of PDB missing residues 
            chain_orig_index = chain_orig_index + 1 # real count, after 1, become 1.
            if geneseq_alignedpart[hhalned_idx] == "-": # case of PDB insertion
                n_gap_insert = n_gap_insert + 1   #del pdb_fullid_list[chain_orig_index-1]
                continue
            else:
                try:
                    out.write('%s %s\n' %(pdb_fullid_list[chain_orig_index-1],geneseq_alignedpart_nogap_fullgeneindex[hhalned_idx-n_gap_insert]))
                except IndexError:
                    print('IndexError for PDB_fullid or GeneSeq_aligned:',chain_orig_index,hhalned_idx)
        else:
            n_gap_miss = n_gap_miss + 1

if __name__ == '__main__': 
#When executed by importing a function. e.g. MapIndex_PDB_GeneSeq.main(PDBid,chainname1,seq1file)
    if len(sys.argv) < 5:
        print ("Command: python MapIndex_PDB_GeneSeq_HHalign.py PDB chainID Fastafile Map_output \n")
        print ("Check your inputs.")
        print ("python MapIndex_PDB_GeneSeq_HHalign.py",sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
        exit()
    pairs = main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])

#----------------------------------------------------
"""
"""

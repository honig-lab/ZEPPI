#!/usr/bin/python3
#
# This program compares PDB complexes like AB.
# output the interfacial contact in the format K12-L31
# K12 is from A, L31 from B.
# Added interface HBonds option; needs more tests.
# Different distances like r, r_CA, r_CB. Here uses r_nonH. 04/01/2020
# ----------------------------------------------------------------------
# Author: Haiqing Zhao, Honig Group at Columbia University
# Created: 11/1/2019
# Last Update: 06/22/2020
# ----------------------------------------------------------------------

import sys
import numpy
import Bio.PDB
from Bio.PDB import Polypeptide as pp
from Bio import Seq

def is_regular_res(res):
    res_id = res.get_id()[1]
    if pp.is_aa(res,standard=True) and res.has_id('CA') and res.has_id('O') and res_id > 0:
        return True
    else:
        return False
    
def calc_residue_distCB(residue_one, residue_two):
    if residue_one.resname !="GLY" and residue_two.resname !="GLY":
        dist = numpy.linalg.norm(residue_one["CB"].get_coord() - residue_two["CB"].get_coord())
    elif residue_one.resname !="GLY" and residue_two.resname =="GLY":
        dist = numpy.linalg.norm(residue_one["CB"].get_coord() - residue_two["O"].get_coord())
    elif residue_one.resname =="GLY" and residue_two.resname !="GLY":
        dist = numpy.linalg.norm(residue_one["O"].get_coord() - residue_two["CB"].get_coord())
    else:
        dist = numpy.linalg.norm(residue_one["O"].get_coord() - residue_two["O"].get_coord())
    return dist
    
def calc_residue_dist(residue_one, residue_two):
    pairwise_dists=[]
    for atom1 in residue_one.get_atoms():
        for atom2 in residue_two.get_atoms():
            if atom1.get_name()[0] !="H" and atom2.get_name()[0] !="H" :
                dist = numpy.linalg.norm(atom1.get_coord() - atom2.get_coord())
                pairwise_dists.append(dist)
    return min(pairwise_dists)

def sequence_obtain(chain):
    sequence=[]
    for res in chain:
        if is_regular_res(res):
            aminos = pp.three_to_one(res.resname)
            sequence.append(aminos)
    sequence = Seq.Seq("".join(sequence)) # Join all amino acids with nothing, return Seq object.
    return sequence

def res_fullid(res):
    res_One = pp.three_to_one(res.resname)
    if res.get_id()[2] ==' ':
        res_id = res.get_id()[1]
    else:
        res_id = str(res.get_id()[1])+res.get_id()[2]
    res_fullid = res_One + str(res_id)
    return res_fullid

def computeIntContact(chainA,chainB,cut=6.05):
    r_cutoff = cut # The closest inter-chain atom distance cut-off
    r_cutoff_CB = 12.0 # 8 in some papers as well; but not good for long residues which often seen on interfaces; example R231-E18, 4PKN_AO
    out_list=[]
    for resA in chainA:
        for resB in chainB:
            resA_id = resA.get_id()[1]
            resB_id = resB.get_id()[1]
            if is_regular_res(resA) and is_regular_res(resB) and resA_id > 0 and resB_id > 0:
                r = calc_residue_dist(resA,resB)
                #r_CB = calc_residue_distCB(resA,resB)
                if r <= float(r_cutoff):# and r_CB <= float(r_cutoff_CB):
                #if r_CB <= float(r_cutoff_CB) and r_CA <= float(r_cutoff_CA):
                    resA_fullid = res_fullid(resA)
                    resB_fullid = res_fullid(resB)
                    out_list.append([resA_fullid,resB_fullid])
    return out_list

def main(PDBid,chainname1,chainname2):
    if PDBid[-4:].lower()==".pdb":
        pdb_file = PDBid
    else:
        pdb_file = PDBid + ".pdb"
    p = Bio.PDB.PDBParser(QUIET = True)
    s = p.get_structure(pdb_file, pdb_file)
    for i in s[0]:
        for j in s[0]:
            if i.id == chainname1 and j.id == chainname2:
                ifc_pairs = computeIntContact(i,j)
    return ifc_pairs

if __name__ == '__main__':
    if len(sys.argv) < 5:
        print ("\nCalcPdbContact.py PDB chainname1 chainname2 Output\n")
        print ("Check your inputs;")
        exit()
    else:
        print ("Program executed. Will write to:", sys.argv[4])
    
    pairs = main(sys.argv[1],sys.argv[2],sys.argv[3])
    
    with open(sys.argv[4], 'w') as file:
        for pair in pairs:
            file.write("%s-%s\n" %(pair[0],pair[1]))

#----------------------------------------------------
"""
"""

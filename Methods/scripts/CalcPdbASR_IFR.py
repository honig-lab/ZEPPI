#!/usr/bin/python3
#
# This program read *.ASA files for protein1 and protein2 (and protein12), output surface residues(SR), ASA etc.
# ASA files can be obtained by running Surfv (https://honig.c2b2.columbia.edu/surface-algorithms).
# ----------------------------------------------------------------------
# Author: Haiqing Zhao
# Honig Group at Columbia University
# Created: 11/1/2019
# Last Update: 09/09/2023
# ----------------------------------------------------------------------

import sys
import Bio.PDB
from Bio.PDB import Polypeptide as pp

r_cutoff = 6.05 # The closed inter-chain atom distance cut-off

def is_regular_res(res):
    res_id = res.get_id()[1]
    if pp.is_aa(res,standard=True) and res.get_id()[0].isspace() and res.has_id('CA') and res.has_id('O') and res_id > 0:
        return True
    else:
        return False
    
def res_fullid(res):
    res_One = pp.three_to_one(res.resname)
    if res.get_id()[2] ==' ':
        res_id = res.get_id()[1]
    else:
        res_id = str(res.get_id()[1])+res.get_id()[2]
    res_fullid = res_One + str(res_id)
    return res_fullid
    
#ALA    1    54.92
#GLY   -2   150.59 (ignore)
def readASA_getSurfaceResi_fullid(ASAfile):
    infile = open(ASAfile, 'r')
    lines = infile.readlines()
    surfaceResi = []
    for line in lines:
        sline = " ".join(line.split())
        try:
            resname = sline.split(' ')[0]
            res_One = pp.three_to_one(resname)
            res_id = sline.split(' ')[1]
            asa = sline.split(' ')[2]
            if float(asa) > 0.0 and res_id[0] != '-' and res_id.isdigit(): #int(res_id) >0: ##could be 45B, (2bum_A); 81A, 81B (2bum_B)
                res_fullid = res_One + str(res_id)
                surfaceResi.append(res_fullid)
        except KeyError:
            continue
    return surfaceResi

def readASA_getSurfaceResi(ASAfile):
    infile = open(ASAfile, 'r')
    lines = infile.readlines()
    surfaceResi = []
    for line in lines:
        sline = " ".join(line.split())
        try:
            resname = sline.split(' ')[0]
            res_One = pp.three_to_one(resname)
            res_id = sline.split(' ')[1]
            asa = sline.split(' ')[2]
            if float(asa) > 0.0 and res_id.isdigit(): # exclude 45B, (2bum_A); 81A, 81B (2bum_B)
                res_fullid = res_One + str(res_id)
                surfaceResi.append(res_fullid)
        except KeyError:
            continue
    if len(surfaceResi) < 1:
        print('Warning: no surface residues are found!')
    return surfaceResi

def readASA_getResi(ASAfile):
    infile = open(ASAfile, 'r')
    lines = infile.readlines()
    allResi = []#{}#collections.OrderedDict()
    for line in lines:
        sline = " ".join(line.split())
        try:
            resname = sline.split(' ')[0]
            res_One = pp.three_to_one(resname)
            res_id = sline.split(' ')[1]
            if res_id.isdigit() and res_id[0] != '-':
                asa = sline.split(' ')[2]
                res_fullid = res_One + str(res_id)
                allResi.append([res_fullid,float(asa)])
        except KeyError:
            continue
    return allResi

def readBASA_getIntfResi(ASAfile1,ASAfile2,ASAfile12):
    list1 = readASA_getResi(ASAfile1); N1 = len(list1)
    list2 = readASA_getResi(ASAfile2); N2 = len(list2)
    list12_1 = readASA_getResi(ASAfile12)[:N1]
    list12_2 = readASA_getResi(ASAfile12)[N1:]
    buried_interface1,buried_interface2={},{}
    for i in range(N1):
        if list1[i][0] == list12_1[i][0] and list12_1[i][1] < list1[i][1]:
            basa = list12_1[i][1] - list1[i][1]
            buried_interface1.update({list1[i][0]:basa})
    for i in range(N2):
        if list2[i][0] == list12_2[i][0] and list12_2[i][1] < list2[i][1]:
            basa = list12_2[i][1] - list2[i][1]
            buried_interface2.update({list2[i][0]:basa})
    print("# of resi/surf_resi in 1stASAfile_protein:",N1,len(buried_interface1))
    print("# of resi/surf_resi in 2ndASAfile_protein:",N2,len(buried_interface2))

    return buried_interface1,buried_interface2

def main_ASA12(ASAfile1,ASAfile2,ASAfile12):
    ASR1 = readASA_getSurfaceResi(ASAfile1)
    ASR2 = readASA_getSurfaceResi(ASAfile2)
    IFR1,IFR2 = readBASA_getIntfResi(ASAfile1,ASAfile2,ASAfile12)
    return ASR1,ASR2,IFR1,IFR2 #,ifc_pairs

def main_ASA(ASAfile1,ASAfile2):
    ASR1 = readASA_getSurfaceResi(ASAfile1)
    ASR2 = readASA_getSurfaceResi(ASAfile2)
    return ASR1,ASR2

if __name__ == '__main__':
    if len(sys.argv) not in [4,5]:
        print ("\nCalcPdbASR_IFR.py ASAfile1 ASAfile2 (ASA12) Output_ASR\n") 
        # potentially can consider (ASA12) for BASA 
        print ("Check your inputs.")
        exit()
    else:
        print ("Program executed. Will write to:", sys.argv[-1])
        n_input = len(sys.argv)

    if n_input ==5:
        oASR1,oASR2,oIFR1,oIFR2 = main_ASA12(sys.argv[1],sys.argv[2],sys.argv[3])
        with open(sys.argv[-1], 'w') as file:
            file.write(','.join(oASR1) + '\n')
            file.write(','.join(oASR2) + '\n')
            file.write(','.join(oIFR1) + '\n')
            file.write(','.join(oIFR2) + '\n')
    elif n_input ==4:
        oASR1,oASR2 = main_ASA(sys.argv[1],sys.argv[2])
        with open(sys.argv[-1], 'w') as file:
            file.write(','.join(oASR1) + '\n')
            file.write(','.join(oASR2) + '\n')


#!/usr/bin/python3
#
# Name: MapIFCtoIFCseq.py 
# A simple matching program. 
# Match *.asr and *.seq_map file, => *.asr_seq.
# Add Donald_style res_map for *.seq_map
# Add PrePPI IFC 
# ----------------------------------------------------------------------
# Author: Haiqing Zhao
# Honig Group at Columbia University
# Created: 05/27/2020
# Last Update: 12/10/2022
# ----------------------------------------------------------------------
import os,sys
if len(sys.argv) < 5:
    print ("COMMAND: python MapPDBIndex_toSeqIndex.py IFCfile_or_ASRfile SeqMapf1 SeqMapf2 IFCseq_file\n")
    print ("Check your inputs")
    exit()

#c2b2pdb_db="/ifs/data/c2b2/bh_lab/shares/databases/pdb/"
#Bacfasta_db="/ifs/data/c2b2/bh_lab/hz2592/Bacterial_PDB_MI/OrthGOPHER/Uni_FASTA/"
#Seqmap_dir="/ifs/data/c2b2/bh_lab/hz2592/Bacterial_PDB_MI/PDB_IFC_MI/Seqmap/"
#ASA_dir="/ifs/data/c2b2/bh_lab/hz2592/Bacterial_PDB_MI/PDB_IFC_MI/ASA/"
#IFC_dir="/ifs/data/c2b2/bh_lab/hz2592/Bacterial_PDB_MI/PDB_IFC_MI/IFC_BASA/"
IFCfile=False;ASRfile=False

#== Check IFC or ASR file; if file exists; read

if sys.argv[1][-4:] == ".ifc":
    IFCfile=sys.argv[1]; onlyindex=False
    if os.path.exists(IFCfile) and os.path.getsize(IFCfile) >0:
        contactpairs_list = []
        with open(IFCfile,'r') as file:
            for line in file:
                rp1,rp2 = line.rstrip().split('-')
                contactpairs_list.append([rp1,rp2])
        if contactpairs_list[0][0].isnumeric(): onlyindex=True
    else:
        print(IFCfile, 'not exist or empty. Exit.')
        exit()

elif sys.argv[1][-4:] == ".asr":
    ASRfile=sys.argv[1]
    if os.path.exists(ASRfile) and os.path.getsize(ASRfile) >0:
        with open(ASRfile, "r") as srfile:
            lines = srfile.readlines()
            if len(lines) == 4:
                p1_sr_indexlist, p2_sr_indexlist, p1_ifr_indexlist, p2_ifr_indexlist = lines[0:4]
                p1_sr_indexlist = p1_sr_indexlist.strip('\n').split(',')
                p2_sr_indexlist = p2_sr_indexlist.strip('\n').split(',')
                p1_ifr_indexlist = p1_ifr_indexlist.strip('\n').split(',')
                p2_ifr_indexlist = p2_ifr_indexlist.strip('\n').split(',')
            elif len(lines) == 2:
                p1_sr_indexlist, p2_sr_indexlist = lines[0:2]
                p1_sr_indexlist = p1_sr_indexlist.strip('\n').split(',')
                p2_sr_indexlist = p2_sr_indexlist.strip('\n').split(',')
                p1_ifr_indexlist, p2_ifr_indexlist=[],[]
    else:
        print(ASRfile, 'not exist or empty. Exit.')
        exit()
else:
    print('Input file must be *.ifc or *.asr')
    exit()

#== Check and read Seq_Mapfiles
mapfile_list = [sys.argv[2],sys.argv[3]]
mapdict_list = [{},{}]
for i in [0,1]:
    if os.path.exists(mapfile_list[i]):
        pass
    else:
        print("No seqmap file found; Check path or run MapIndex_PDB_GeneSeq_HHalign: ")
        exit()

    with open(mapfile_list[i], 'r') as file:
        first_line = file.readline()
        if ',' in first_line: maptype=1;
        else: maptype=2
    with open(mapfile_list[i], 'r') as file:
        if maptype==1:            
            for line in file:
                (key,val) = line.rstrip().split(',')
                mapdict_list[i][key] = val
        elif maptype==2:
            for line in file:
                (key,val) = line.rstrip().split()
                mapdict_list[i][key] = val

#== Map IFC contacts to seq_maps, save IFC_seq
if IFCfile:
    if onlyindex: 
        for i in list(mapdict_list[0]): mapdict_list[0][i[1:]] = mapdict_list[0][i][1:]
        for i in list(mapdict_list[1]): mapdict_list[1][i[1:]] = mapdict_list[1][i][1:]

    contactpairs_list_new = []
    if mapdict_list[0] and mapdict_list[1]:
        for pair in contactpairs_list:    
            if pair[0] in mapdict_list[0] and pair[1] in mapdict_list[1]:
                res1 = mapdict_list[0][pair[0]]
                res2 = mapdict_list[1][pair[1]]
                contactpairs_list_new.append([res1,res2])
            elif pair[0] not in mapdict_list[0]:
                print("Warning: ifc_pair[0] not in seqmap; skip ",pair[0],IFCfile)
            elif pair[1] not in mapdict_list[1]:
                print("Warning: ifc_pair[1] not in seqmap; skip ",pair[1],IFCfile)
    else:
        print("Warning: mapdict_list empty.")
   #== Write IFC_seq file
    if len(contactpairs_list_new) != 0:
        with open(sys.argv[4], 'w') as file:
            for pair in contactpairs_list_new:
                file.write("%s-%s\n" %(pair[0],pair[1]))

elif ASRfile:
    newlists=[]
    for idx, oldlist in enumerate([p1_sr_indexlist, p2_sr_indexlist, p1_ifr_indexlist,p2_ifr_indexlist]):
        newlist=[]
        for iID in oldlist:
            if (idx % 2) == 0 and iID in mapdict_list[0]:
                iID_new = mapdict_list[0][iID]
                newlist.append(iID_new)
            if (idx % 2) != 0 and iID in mapdict_list[1]:
                iID_new = mapdict_list[1][iID]
                newlist.append(iID_new)
        newlists.append(newlist)
    #== Write ASR_seq file
    with open(sys.argv[4], 'w') as f:
        for newlist in newlists:
            f.write(str(','.join(newlist)))
            f.write('\n')

 
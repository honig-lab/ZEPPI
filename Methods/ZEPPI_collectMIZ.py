#!/usr/bin/python3
#
# This program reads the MIZ files for the PPIs in the given list
# Outputs detailed properties of PPI, and the calculated MI/Con-related metrics.
#
# ----------------------------------------------------------------------
# Written by Haiqing Zhao
# Honig Group at Columbia University
# Last Update: 08/10/2023
# ----------------------------------------------------------------------

import numpy as np
import sys

if len(sys.argv) < 4:
    print ("\n python ZEPPI_collectMIZ.py PPI_csvfile output_csv metricfile_path \n")
    exit()

tiny = 1e-6
def read_MIZfile (MIZfile):
    with open(MIZfile) as f:
        lines = f.readlines()
        if len(lines) > 4:
            N_msa = int(float(lines[1].split(" ")[4]))
            N_eff = int(float(lines[1].split(" ")[8]))
            N_ifr_eff1 = int(lines[2].split(",")[1])
            N_ifr_eff2 = int(lines[2].split(",")[4])
            N_ifr = int(lines[1].split(" ")[3].split(",")[0].split("(")[1]) + int(lines[1].split(" ")[3].split(",")[1].split(")")[0])
            #N_ifr_eff = N_ifr #_eff1 + N_ifr_eff2
            N_ifr_eff = N_ifr_eff1 + N_ifr_eff2
            Max_NM_array = np.array(lines[4].split()[4:8]).astype(np.float32)
            Mean_NM_array = np.array(lines[5].split()[4:8]).astype(np.float32)
            Max_IFC_array = np.array(lines[6].split()[4:8]).astype(np.float32)
            Mean_IFC_array = np.array(lines[7].split()[4:8]).astype(np.float32)
            Z_Max_NM_array = np.array(lines[8].split()[4:8]).astype(np.float32)
            Z_Mean_NM_array = np.array(lines[9].split()[4:8]).astype(np.float32)
            Z_Max_IFC_array = np.array(lines[10].split()[4:8]).astype(np.float32)
            Z_Mean_IFC_array = np.array(lines[11].split()[4:8]).astype(np.float32)
            Mean_Max_IFC_array = Mean_IFC_array/(Max_IFC_array+tiny)
            Mean_Max_NM_array = Mean_NM_array/(Max_NM_array+tiny)
            return Max_NM_array, Mean_NM_array, Max_IFC_array, Mean_IFC_array,N_msa, N_eff,N_ifr, N_ifr_eff, Mean_Max_IFC_array,Mean_Max_NM_array,Z_Max_NM_array, Z_Mean_NM_array, Z_Max_IFC_array, Z_Mean_IFC_array

Zscore_allpairs_Max_NM, Zscore_allpairs_Mean_NM, Zscore_allpairs_Max_IFC, Zscore_allpairs_Mean_IFC =[],[],[],[]
Neff_allpairs, Nifr_allpairs = [],[]; pair_avail = [];LR_allpairs=[]
real_allpairs_Max_IFC,real_allpairs_Mean_IFC=[],[]
Zscore_allpairs_MeanMax_IFC, Zscore_allpairs_MeanMax_NM =[],[]

headerlist1 = ['PPI','Nmsa','Neff','Nifr', 'Nifr_eff']

metriclist = ['Max_NM_ar','Mean_NM_ar','Max_IFC_ar','Mean_IFC_ar','MeanMax_IFC_ar','MeanMax_NM_ar','Z_Max_NM_ar','Z_Mean_NM_ar','Z_Max_IFC_ar', 'Z_Mean_IFC_ar']
headerlist = headerlist1
for metric in metriclist:
    for i in range(1,5,1):
        headerlist.append(metric+str(i))
headerlist=','.join(headerlist)

pairfile=sys.argv[1]
filetowrite=sys.argv[2]

with open(filetowrite, 'w') as f:
    f.write(headerlist+'\n') # 
towrite = open(filetowrite,'ab')

ppipair_list = [];dic_LR={}
pairfile=open(pairfile,'r')
next(pairfile)
for line in pairfile:
    #[uni1,uni2] = line.strip().split(",")[0:2];
    #pairname=uni1+'_'+uni2
    #LR=line.strip().split(",")[4]
    #dic_LR.update({uni1+'_'+uni2:LR})
    pdbid = line.split(",")[0]; [chain1,chain2] = line.split(",")[1:3]
    pairname=pdbid+"_"+chain1+chain2; 
    
    ppipair_list.append(pairname)

mifile_path=sys.argv[3]
if mifile_path[-1] != "/": mifile_path=mifile_path+"/"
for pair in ppipair_list[:]:
    mifilename = mifile_path+pair+".miZ"
    try:
        Max_NM_ar, Mean_NM_ar, Max_IFC_ar, Mean_IFC_ar,Nmsa, Neff, Nifr, Nifr_eff, Mean_Max_IFC_ar, Mean_Max_NM_ar, Z_Max_NM_ar, Z_Mean_NM_ar, Z_Max_IFC_ar, Z_Mean_IFC_ar = read_MIZfile(mifilename);
    except (FileNotFoundError,TypeError) as fe:
        print("No",mifilename)
        continue

    Neff_allpairs.append(Neff); 
    data = np.concatenate((np.array([pair]),np.array([Nmsa]),np.array([Neff]),np.array([Nifr]),np.array([Nifr_eff]),Max_NM_ar,Mean_NM_ar, Max_IFC_ar, Mean_IFC_ar,Mean_Max_IFC_ar,Mean_Max_NM_ar,Z_Max_NM_ar, Z_Mean_NM_ar, Z_Max_IFC_ar, Z_Mean_IFC_ar))
    dt_Name = ('','U32'); dt_score = ('','float16'); 
    dt_all = list((dt_Name, ) * len(headerlist1)) + list((dt_score, ) * (len(data)-len(headerlist1)))
    data = [tuple(data.tolist())];
    data = np.array(data, dtype=dt_all)
    fmt_all = '%s' + ', %f' * (len(data)-1)
    np.savetxt(towrite, np.array(data), fmt= fmt_all,delimiter=',') 

print("Number of input PPIs and that of effective calculations:",len(ppipair_list),len(Neff_allpairs))
towrite.close()

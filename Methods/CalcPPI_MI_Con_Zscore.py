#!/usr/bin/python3
#
# This program computes the the mutual information of (a,b) based on multiple sequence alignments.
# Inputs be the MSA files(from hmmalign) whose 1st alignment is the query itself.
# Added checking rowwise coverage, 50% for SR of both proteins.
# Added checking columnwise coverage, 50% for SR of both proteins.
# A new change is to use numpy instead of character in MI calculation 
# Added APC for p1-all to p2-all.
# Added numba acceleration 
# Added Z-score for the real IFC among 100 fake interfaces {from nonInf Surf residues}
# Added PrePPI domain recognization from filename
# Added homo option, save 25% cpu time
# ----------------------------------------------------------------------
# Copyright (C) Haiqing Zhao
# Honig Group at Columbia University
# Last Update: 01/12/2023
# ----------------------------------------------------------------------

import sys,time
import numpy as np
from Bio import AlignIO
import numba as nb
import itertools,random

if len(sys.argv) != 6:
    print ("\n python MSA_MI_Gremlin.py MSA_file1 MSA_file2 IFC_file SurfaceRes_file Output_file\n")
    print ("Check your inputs; IFC files format: D12-R45 or 12-45 ")
    exit()

def mk_msa(seqs): #'''one hot encode msa'''
    alphabet = list("ARNDCQEGHILKMFPSTWYV-")
    states = len(alphabet)
    alpha = np.array(alphabet, dtype='|S1').view(np.uint8)
    msa = np.array([list(s) for s in seqs], dtype='|S1').view(np.uint8)  
    for n in range(states):
        msa[msa == alpha[n]] = n  
    msa[msa > states] = states-1
    return np.eye(states)[msa]

if sys.argv[1] == sys.argv[2]: homo=True
else: homo=False

#--- READ the MSA information from MSA of protein1 and protein2  ---------------------------#

alignment_p1 = AlignIO.read(sys.argv[1], 'fasta')
Num_rows_p1 = len(alignment_p1)
Num_columns_p1 = alignment_p1.get_alignment_length()

alignment_p2 = AlignIO.read(sys.argv[2], 'fasta')
Num_rows_p2 = len(alignment_p2)
Num_columns_p2 = alignment_p2.get_alignment_length()

if len(sys.argv[1].split('/')[-1].split('.'))==3: firstindex_p1=int(sys.argv[1].split('/')[-1].split('.')[1].split('_')[0])
else: firstindex_p1 = 1
if len(sys.argv[2].split('/')[-1].split('.'))==3: firstindex_p2=int(sys.argv[2].split('/')[-1].split('.')[1].split('_')[0])
else: firstindex_p2 = 1

#print("First Indexes as:",firstindex_p1,firstindex_p2)
#--- Pair common species for p1-orthologs and p2-orthologs -------------------------#

species_proteinseq1,species_proteinseq2 = {},{}
species_all1,species_all2 = [],[]
for record in alignment_p1:
    species = record.id.split(' ')[0] #.split(":")[0]
    species_proteinseq1.update({species:record.seq.upper()})

for record in alignment_p2:
    species = record.id.split(' ')[0] #.split(":")[0]
    species_proteinseq2.update({species:record.seq.upper()})

print("Number of species in two MSAs:", len(species_proteinseq1),len(species_proteinseq2))
common_species=[]
def intersection(lst1, lst2):
    return list(set(lst1) & set(lst2)) 

common_species = intersection(species_proteinseq1,species_proteinseq2)
common_species = sorted(common_species) #important for co-existing

#--- READ the interface contact/residues information from IFC file ---------------------------#

resipairs = open(sys.argv[3], "r") # Default input numbering from 1
rp_fullname_list,rp_seqindex_list,rp_resname_list = [],[],[]
for rp in resipairs:
    rp = rp.strip(' \t\n\r').split('-')
    if len(rp) !=2:
        print("Warning: empty line in IFC file", rp)  
        continue
    else: rp_fullname = [x for x in rp]
    if rp[0].isdigit():
        rp_seqindex = [int(rp[0])-firstindex_p1,int(rp[1])-firstindex_p2]
        #rp_seqindex = [int(x)-1 for x in rp] # ifc_res_index starts from 1;for ifc format: 12-45
    elif rp[0][1:].isdigit():
        rp_seqindex = [int(rp[0][1:])-firstindex_p1,int(rp[1][1:])-firstindex_p2] #for ifc format: D12-R45
    rp_resname = [str(x[0]) for x in rp_fullname]
    rp_fullname_list.append(rp_fullname)
    rp_seqindex_list.append(rp_seqindex)
    rp_resname_list.append(rp_resname)
resipairs.close()

if len(rp_fullname_list) == 0: print("Emply IFC file; exit.") ;exit()

p1_ifcres_indexlist = sorted(list(set([i[0] for i in rp_seqindex_list])))
p2_ifcres_indexlist = sorted(list(set([i[1] for i in rp_seqindex_list])))
Nifc=[len(p1_ifcres_indexlist),len(p2_ifcres_indexlist)]
print('Number of IFC residues in protein1 and protein2:',len(p1_ifcres_indexlist),len(p2_ifcres_indexlist))

#--- READ SurfaceRes information; Generate fake_IFCs---------------------------#

with open(sys.argv[4], "r") as srfile:
    p1_sr_indexlist, p2_sr_indexlist =srfile.readlines()[0:2]
    p1_sr_indexlist = p1_sr_indexlist.strip('\n').split(',')
    p2_sr_indexlist = p2_sr_indexlist.strip('\n').split(',')
    if p1_sr_indexlist[0][0].isdigit():
        p1_allsr_indexlist = [int(x[0])-firstindex_p1 for x in p1_sr_indexlist]
        p2_allsr_indexlist = [int(x[0])-firstindex_p2 for x in p2_sr_indexlist]
    else:
        p1_allsr_indexlist = [int(x[1:])-firstindex_p1 for x in p1_sr_indexlist]
        p2_allsr_indexlist = [int(x[1:])-firstindex_p2 for x in p2_sr_indexlist]
p1_allsr_indexlist=list(set(p1_allsr_indexlist)); p2_allsr_indexlist=list(set(p2_allsr_indexlist))
print('Number of Surface Residues:',len(p1_sr_indexlist),len(p2_sr_indexlist))

p1_sr_indexlist = [item for item in p1_allsr_indexlist if item not in p1_ifcres_indexlist]
p2_sr_indexlist = [item for item in p2_allsr_indexlist if item not in p2_ifcres_indexlist]

print('Number of noninterfacial Surface Residues:',len(p1_sr_indexlist),len(p2_sr_indexlist))

if len(p1_ifcres_indexlist) < len(p1_sr_indexlist): pass
else:
    print("The # IFR_P1 is more than nonIFR_SurfaceRes:",len(p1_ifcres_indexlist),len(p1_sr_indexlist))
    p1_sr_indexlist = list(range(0,Num_columns_p1))

if len(p2_ifcres_indexlist) < len(p2_sr_indexlist): pass
else:
    print("The # IFR_P2 is more than nonIFR_SurfaceRes:",len(p2_ifcres_indexlist),len(p2_sr_indexlist))
    p2_sr_indexlist = list(range(0,Num_columns_p2))

def generatefake_rp_index(real_rp_index_list,p1_fake_ifr,p2_fake_ifr):
    def numberlink(indexin1,list1,list2):
        dict_oldtonew = dict(zip(list1, list2))
        if indexin1 in list1:
            return dict_oldtonew[indexin1]
    fake_rp_index_list=[]
    p1_ifr_indexlist = sorted(list(set([i[0] for i in real_rp_index_list])))
    p2_ifr_indexlist = sorted(list(set([i[1] for i in real_rp_index_list])))
    for rpindex in real_rp_index_list:
        index_fa = numberlink(rpindex[0],p1_ifr_indexlist,p1_fake_ifr)
        index_fb = numberlink(rpindex[1],p2_ifr_indexlist,p2_fake_ifr)
        if index_fa != None and index_fb != None:
            fake_rp_index_list.append([index_fa,index_fb])      
    return fake_rp_index_list

if homo:
    p12_fifr_lists=[]
    p12_sr_indexlist = list(set(p1_sr_indexlist + p2_sr_indexlist))
    p12_ifcres_indexlist = set(p1_ifcres_indexlist + p2_ifcres_indexlist)
    for i in range(0,100):
        p12_fifr_lists.append(random.sample(p12_sr_indexlist,len(p12_ifcres_indexlist)))
else:
    p1_fifr_lists=[];p2_fifr_lists=[]
    for i in range(0,100):
        p1_fifr_lists.append(random.sample(p1_sr_indexlist,len(p1_ifcres_indexlist)))
        p2_fifr_lists.append(random.sample(p2_sr_indexlist,len(p2_ifcres_indexlist)))

def generatefake_rp_index_homo(real_rp_index_list,p12_fake_ifr):
    fake_rp_index_list=[]
    p12_ifr_indexlist = sorted(set([i[0] for i in real_rp_index_list] + [i[1] for i in real_rp_index_list]))
    p12_fifr_indexlist = p12_fake_ifr
    dict_oldtonew = dict(zip(p12_ifr_indexlist, p12_fifr_indexlist))
    for rpindex in real_rp_index_list:
        index_fa = dict_oldtonew[rpindex[0]]
        index_fb = dict_oldtonew[rpindex[1]] 
        if index_fa != None and index_fb != None:
            fake_rp_index_list.append([index_fa,index_fb])      
    return fake_rp_index_list

#--- Select row-wise those cover >50% of SR on p1 and p2 ----#

def check_row_cov (long_seq,index_list,j=0):
    for i in index_list:
        if long_seq[i] == "-":
            j=j+1
    return 1-j/len(index_list)

droprow_species=[]
for k in common_species:
    portion1 = check_row_cov(species_proteinseq1[k],p1_allsr_indexlist)
    portion2 = check_row_cov(species_proteinseq2[k],p2_allsr_indexlist)
    if portion1 > 0.5 and portion2 > 0.5: pass
    else: droprow_species.append(k)

common_species_cov = list(set(common_species) - set(droprow_species))
common_species_cov = sorted(common_species_cov)

if len(common_species_cov) == 0 or len(common_species_cov) == 1:
    print("WARNING: no common_species_cov_nRd. Exit.",len(common_species),len(common_species_cov),len(common_species_cov))
    with open(sys.argv[5], 'w') as out:
        out.write("# Protein1_MSA1 and Protein2_MSA2: (%i, %i), (%i, %i).\n" % (Num_rows_p1,Num_columns_p1,Num_rows_p2,Num_columns_p2))
        out.write("# %i IFCs; %i common_species, %i common_cov\n" % (len(rp_fullname_list),len(common_species),len(common_species_cov)))
    exit()

#--- Select column-wise: MSA columns cover >50% (nongap) ----#
def dropcol_belowcov(position_list,MSA,use_rows,cutoff=0.5):
    todrop_list=[]
    for i in position_list:
        j=0;
        for s in use_rows:
            if MSA[s][i] == "-": j=j+1
        gap_portion = j/len(common_species_cov)
        if gap_portion >= cutoff: todrop_list.append(i)
    return [item for item in position_list if item not in todrop_list]

p1_sr_indexlist =dropcol_belowcov(p1_sr_indexlist,species_proteinseq1,common_species_cov)
p2_sr_indexlist =dropcol_belowcov(p2_sr_indexlist,species_proteinseq2,common_species_cov)
p1_ifcres_indexlist =dropcol_belowcov(p1_ifcres_indexlist,species_proteinseq1,common_species_cov)
p2_ifcres_indexlist =dropcol_belowcov(p2_ifcres_indexlist,species_proteinseq2,common_species_cov)

print('Number of IFC residues (after column-check):',len(p1_ifcres_indexlist),len(p2_ifcres_indexlist))

#--- concatenate MSA for P1-P2 ----#
def contanate_seq12_allMSA12(MSA1,MSA2,use_row):
    out_seq12={}
    for k in use_row:
        out_seq12.update({k:MSA1[k]+MSA2[k]})
    return list(out_seq12.values())

sequences_p12 = contanate_seq12_allMSA12(species_proteinseq1,species_proteinseq2,common_species_cov)

"""
"""

@nb.njit()
def entropy_col(msa_np,column_i):
    x=msa_np[:,column_i].sum(axis=0)
    x=x[x!=0]
    p=x/x.sum(axis=0)
    entr = (-p*np.log(p)).sum(axis=0)
    return entr

@nb.njit()
def entropy(arr2d_p,entr=0):
    for i in arr2d_p:
        for j in i:
            if j != 0:
                entr += -j*np.log(j)
    return entr

@nb.njit()
def jointentropy(msa_np,column_i,column_j):
    colmsa_i=msa_np[:,column_i]
    colmsa_j=msa_np[:,column_j]
    joint=np.zeros((21, 21));consv=np.zeros((21, 21))
    for n in range(0,len(colmsa_i)):
        loc1=colmsa_i[n].argmax() # faster than np.nonzero, loc1=np.nonzero(colmsa_i[n])[0][0]
        loc2=colmsa_j[n].argmax()
        joint[loc1][loc2]=joint[loc1][loc2]+ colmsa_i[n][loc1]*colmsa_j[n][loc2]
        #consv[loc1][loc2]=consv[loc1][loc2]+ 1
    p = joint/joint.sum();#p=p[p!=0]#entr = (-p*np.log(p)).sum(axis=0)
    entr = entropy(p) #.sum()
    #p = consv/consv.sum();#p=p[p!=0]
    conserv = 1 - entr/(2*np.log(21))
    return entr,conserv

@nb.njit()
def MI(msa_np,column_i,column_j):
    return entropy_col(msa_np,column_i) + entropy_col(msa_np,column_j) - jointentropy(msa_np,column_i,column_j)

@nb.njit()
def compute_MI_p1top2(paired_msa,N1,N2):
    N1toN2 = np.zeros((N1,N2), dtype='float');N1toN2_consv= np.zeros((N1,N2), dtype='float')
    entr_N1=[];entr_N2=[]
    for i in range(N1):
        entr_N1.append(entropy_col(paired_msa,i))
    for j in range(N2):
        entr_N2.append(entropy_col(paired_msa,N1+j))
    for i in range(N1):
        for j in range(N2):
            jointentr,consv = jointentropy(paired_msa,i,N1+j)
            mi = entr_N1[i]+entr_N2[j]-jointentr
            N1toN2[i][j] = mi; N1toN2_consv[i][j] = consv
    if N1toN2.min()<0:
        print("Warning: negative MI value. May cause very large APC.")
    return N1toN2,N1toN2_consv

@nb.njit()
def compute_MI_p1top1(paired_msa,N1):
    N2=N1
    N1toN2 = np.zeros((N1,N2), dtype='float');N1toN2_consv= np.zeros((N1,N2), dtype='float')
    entr_N1=[];
    for i in range(N1):
        entr_N1.append(entropy_col(paired_msa,i))
    for i in range(N1):
        for j in range(i+1,N1):
            jointentr,consv = jointentropy(paired_msa,i,N1+j)
            mi = entr_N1[i]+entr_N1[j]-jointentr
            N1toN2[i][j] = mi; N1toN2_consv[i][j] = consv
    N1toN2=N1toN2 + N1toN2.T
    N1toN2_consv=N1toN2_consv+N1toN2_consv.T
    if N1toN2.min()<0:
        print("Warning: negative MI value. May cause very large APC.")
    return N1toN2,N1toN2_consv

def get_asc(raw):
    asc = np.sum(raw,axis=1)/raw.shape[1] + np.sum(raw,axis=0)/raw.shape[0] - np.sum(raw)/(raw.shape[0]*raw.shape[1])
    return (raw-asc)

def get_apc(raw):
    a = np.sum(raw,axis=1,keepdims=True)/raw.shape[1]
    b = np.sum(raw,axis=0,keepdims=True)/raw.shape[0]
    c = np.sum(raw)/(raw.shape[0]*raw.shape[1]) + 1e-10
    apced = raw - a*b/c
    #if c !=0: #or + tiny
    #    apced = raw - a*b/c
    #else:
    #    print("Warning: apc=0; changed to ASC instead.")
    #    apced = get_asc(raw) #np.fill_diagonal(apced,0)
    return apced

### --- Computations and Output --- ###
N1,N2 =len(p1_ifcres_indexlist),len(p2_ifcres_indexlist)
end0 = time.time()
msa_numpy_p12 = mk_msa(sequences_p12)
if homo:
    MI_mat_P12,Consv_mat_P12 = compute_MI_p1top1(msa_numpy_p12,Num_columns_p1)
else:
    MI_mat_P12,Consv_mat_P12 = compute_MI_p1top2(msa_numpy_p12,Num_columns_p1,Num_columns_p2)
    
MI_APC_mat=get_apc(MI_mat_P12)
Consv_APC_mat=get_apc(Consv_mat_P12)
end1 = time.time()
print("Time used for MI & APC: %.4f" %(end1-end0))

msa_weights = np.ones(len(msa_numpy_p12))

## Compute NM_ | IFC_ max and mean, Z-score
def max_mean_submatrix(matrix,indexpair):
    sub=[]
    for p in indexpair:
        sub.append(matrix[p[0]][p[1]])
    return np.max(sub),np.mean(sub)

def max_mean_NM_submatrix(matrix,indexpair):
    sub=[]
    p1_indexlist = sorted(list(set([i[0] for i in indexpair])))
    p2_indexlist = sorted(list(set([i[1] for i in indexpair])))
    NM_indexpair = list(itertools.product(p1_indexlist, p2_indexlist))
    for p in NM_indexpair:
        sub.append(matrix[p[0]][p[1]])
    try:
        return [np.max(sub),np.mean(sub)]
    except ValueError:
        return [None,None]

def Col_Zscore (inputmatrix,reference):
    inputmatrix = np.array(inputmatrix)
    reference = np.array(reference)
    #print("Z-score function checking shapes of input:",inputmatrix.shape,reference.shape) #(10, 6, 2) (6, 2)
    if reference.shape[0] != inputmatrix.shape[1] :
        print ('Error of unequal dimensions:', reference.shape[0],inputmatrix.shape[1])
        exit()
    real_fake = np.append([reference],inputmatrix,axis= 0)
    column_mean = np.nanmean(real_fake, axis=0);
    column_std = np.nanstd(real_fake, axis=0) + 1e-10;
    try:
        column_zscore = ((reference)-(column_mean))/(column_std)
    except (ValueError, ZeroDivisionError):
        print('ValueError or ZeroDivisionError')
        pass
    else:
        return column_zscore

MI_mat = MI_mat_P12;Consv_mat=Consv_mat_P12
NM_maxmean_list=[];IFC_maxmean_list=[]
matr_toStat=[MI_mat,MI_APC_mat,Consv_mat,Consv_APC_mat]

for matr in matr_toStat:
    NM_maxmean_list.append(max_mean_NM_submatrix(matr,rp_seqindex_list))
    IFC_maxmean_list.append(max_mean_submatrix(matr,rp_seqindex_list))

# Generate fake_rp_seqindex_list based on rp_seqindex_list and get Z-score
NM_maxmean_f_list=[];IFC_maxmean_f_list=[]
for i in range(0,100):
    if homo:
        f_rp_seqindex_list = generatefake_rp_index_homo(rp_seqindex_list,p12_fifr_lists[i])
    else:
        f_rp_seqindex_list = generatefake_rp_index(rp_seqindex_list,p1_fifr_lists[i],p2_fifr_lists[i])
    NM_maxmean_f=[];IFC_maxmean_f=[]
    for matr in matr_toStat:
        NM_maxmean_f.append(max_mean_NM_submatrix(matr,f_rp_seqindex_list))
        IFC_maxmean_f.append(max_mean_submatrix(matr,f_rp_seqindex_list))
    NM_maxmean_f_list.append(NM_maxmean_f)
    IFC_maxmean_f_list.append(IFC_maxmean_f)

Zscore_NM_list=Col_Zscore(NM_maxmean_f_list, NM_maxmean_list)
Zscore_IFC_list=Col_Zscore(IFC_maxmean_f_list, IFC_maxmean_list)

end2 = time.time()
print ("Time used for Z-score: %.4f" %(end2-end1))

with open(sys.argv[5], 'w') as out:
    out.write("# Protein1_MSA and Protein2_MSA: (%i, %i), (%i, %i).\n" % (Num_rows_p1,Num_columns_p1,Num_rows_p2,Num_columns_p2))
    out.write("# %i IFCs (%i,%i); %i common_species, %i common_cov, %.2f Neff.\n" % (len(rp_fullname_list),Nifc[0],Nifc[1],len(common_species),len(common_species_cov),msa_weights.sum()))
    #out.write("# P1_IFC_MSA and P2_IFC_MSA: (%i, %i), (%i, %i).\n" % (len(common_species_cov),len(p1_ifcres_indexlist),len(common_species_cov),len(p2_ifcres_indexlist)))
    out.write("# P1_MSA_IFC_SR and P2_MSA_IFC_SR: (%i, %i, %i), (%i, %i, %i).\n" % (len(common_species),len(p1_ifcres_indexlist),len(p1_sr_indexlist),len(common_species),len(p2_ifcres_indexlist),len(p2_sr_indexlist)))

    out.write("# Header in order of: MI, MI_APC, Consv, Consv_APC\n")
    out.write("# Max in NMmatrix: %.4f\t%.4f\t%.4f\t%.4f\n" %(NM_maxmean_list[0][0],NM_maxmean_list[1][0],NM_maxmean_list[2][0],NM_maxmean_list[3][0]))
    out.write("# Mean in NMmatrx: %.4f\t%.4f\t%.4f\t%.4f\n" %(NM_maxmean_list[0][1],NM_maxmean_list[1][1],NM_maxmean_list[2][1],NM_maxmean_list[3][1]))
    out.write("# Max in IFCpairs: %.4f\t%.4f\t%.4f\t%.4f\n" %(IFC_maxmean_list[0][0],IFC_maxmean_list[1][0],IFC_maxmean_list[2][0],IFC_maxmean_list[3][0]))
    out.write("# Mean in IFCpair: %.4f\t%.4f\t%.4f\t%.4f\n" %(IFC_maxmean_list[0][1],IFC_maxmean_list[1][1],IFC_maxmean_list[2][1],IFC_maxmean_list[3][1]))

    out.write("# Z of max_in__NM: %.4f\t%.4f\t%.4f\t%.4f\n" %(Zscore_NM_list[0][0],Zscore_NM_list[1][0],Zscore_NM_list[2][0],Zscore_NM_list[3][0]))
    out.write("# Z of mean_in_NM: %.4f\t%.4f\t%.4f\t%.4f\n" %(Zscore_NM_list[0][1],Zscore_NM_list[1][1],Zscore_NM_list[2][1],Zscore_NM_list[3][1]))
    out.write("# Z of max_in_IFC: %.4f\t%.4f\t%.4f\t%.4f\n" %(Zscore_IFC_list[0][0],Zscore_IFC_list[1][0],Zscore_IFC_list[2][0],Zscore_IFC_list[3][0]))
    out.write("# Z of meanin_IFC: %.4f\t%.4f\t%.4f\t%.4f\n" %(Zscore_IFC_list[0][1],Zscore_IFC_list[1][1],Zscore_IFC_list[2][1],Zscore_IFC_list[3][1]))
    for rp in rp_fullname_list:
        a_fullname = rp[0]
        b_fullname = rp[1]
        if rp[0].isdigit():
            index_a = int(rp[0]) - firstindex_p1
            index_b = int(rp[1]) - firstindex_p2
        else:
            index_a = int(rp[0][1:]) - firstindex_p1
            index_b = int(rp[1][1:]) - firstindex_p2
        if (index_a in p1_ifcres_indexlist) and (index_b in p2_ifcres_indexlist):
            mi,mi_apc = MI_mat[index_a][index_b], MI_APC_mat[index_a][index_b]
            consv,consv_apc = Consv_mat[index_a][index_b], Consv_APC_mat[index_a][index_b]
            out.write("%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\n" %(a_fullname,b_fullname,mi,mi_apc,consv,consv_apc))

"""
Example output:
# Protein1_MSA and Protein2_MSA: (286, 161), (220, 161).
# 132 IFCs (46,40); 208 common_species, 206 common_cov, 206.00 Neff.
# P1_MSA_IFC_SR and P2_MSA_IFC_SR: (208, 46, 100), (208, 40, 108).
# Header in order of: MI, MI_APC, Consv, Consv_APC
# Max in NMmatrix: 1.5238   0.6424  0.8993  0.1087
# Mean in NMmatrx: 0.5730   0.0274  0.7304  0.0034
# Max in IFCpairs: 0.9896   0.2774  0.8481  0.0468
# Mean in IFCpair: 0.5449   0.0049  0.7469  -0.0004
# Z of max_in__NM: 1.1138   -0.7570 -0.0971 -0.4966
# Z of mean_in_NM: -0.3665  6.5892  7.2810  6.0533
# Z of max_in_IFC: -1.1972  -0.6636 0.3517  -0.7079
# Z of meanin_IFC: -0.9531  -0.0059 6.6532  -0.8316
"""
"""
"""

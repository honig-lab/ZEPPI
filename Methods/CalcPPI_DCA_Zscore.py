#!/usr/bin/python3
#
# This program computes direct coupling analysis and direct information for PPIs.
# It uses mean-field approximation to get the maximum entropy model solution of the partition function. 
# Details refer to Morcos, Faruck, et al. PNAS, 108.49 (2011): E1293-E1301.
# ----------------------------------------------------------------------
# Code by Haiqing Zhao
# Honig Group at Columbia University
# Last Update: 01/12/2023
# ----------------------------------------------------------------------

import numpy as np
from Bio import AlignIO
import itertools,random
from scipy.spatial.distance import *
import numba as nb
import time,sys

if len(sys.argv) != 6:
    print ("\n python MSA_DCA_p1p2.py MSA_file1 MSA_file2 IFC_file SurfaceRes_file Output_file\n")
    print ("Check your inputs; IFC files format: D12-R45 or 12-45 ")
    exit()

msafile1=sys.argv[1]
msafile2=sys.argv[2]
ifcfile=sys.argv[3]
asrfile=sys.argv[4]
outfile=sys.argv[5]

if sys.argv[1] == sys.argv[2]: homo=True
else: homo=False

if len(sys.argv[1].split('/')[-1].split('.'))==3: firstindex_p1=int(sys.argv[1].split('/')[-1].split('.')[1].split('_')[0])
else: firstindex_p1 = 1
if len(sys.argv[2].split('/')[-1].split('.'))==3: firstindex_p2=int(sys.argv[2].split('/')[-1].split('.')[1].split('_')[0])
else: firstindex_p2 = 1


# %%
alignment_p1 = AlignIO.read(msafile1, 'fasta')
species_proteinseq1={}
for record in alignment_p1:
    species = record.id.split(' ')[0] #.split(":")[0]
    species_proteinseq1.update({species:record.seq.upper()})

# %%
species_proteinseq2={}
alignment_p2 = AlignIO.read(msafile2, 'fasta')
def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3
for record in alignment_p2:
    species = record.id.split(' ')[0] #.split(":")[0]
    species_proteinseq2.update({species:record.seq.upper()})


# %%
common_species = intersection(species_proteinseq1,species_proteinseq2)
Num_rows_p1 = len(alignment_p1)
Num_columns_p1 = alignment_p1.get_alignment_length()
Num_rows_p2 = len(alignment_p2)
Num_columns_p2 = alignment_p2.get_alignment_length()


# %%
resipairs = open(ifcfile, "r") # Default input numbering from 1
rp_fullname_list,rp_seqindex_list,rp_resname_list = [],[],[]
for rp in resipairs:
    rp = rp.strip(' \t\n\r').split('-')
    if len(rp) !=2:
        print("Warning: empty line in IFC file", rp)  
        continue
    else: rp_fullname = [x for x in rp]
    if rp[0].isdigit():
        rp_seqindex = [int(rp[0])-firstindex_p1,int(rp[1])-firstindex_p2]
    elif rp[0][1:].isdigit():
        rp_seqindex = [int(rp[0][1:])-firstindex_p1,int(rp[1][1:])-firstindex_p2] #IFC format: D12-R45
    rp_resname = [str(x[0]) for x in rp_fullname]
    rp_fullname_list.append(rp_fullname)
    rp_seqindex_list.append(rp_seqindex)
    rp_resname_list.append(rp_resname)
resipairs.close()

if len(rp_fullname_list) == 0: print("Emply IFC file; exit."); exit()

p1_ifcres_indexlist = sorted(list(set([i[0] for i in rp_seqindex_list])))
p2_ifcres_indexlist = sorted(list(set([i[1] for i in rp_seqindex_list])))
Nifc=[len(p1_ifcres_indexlist),len(p2_ifcres_indexlist)]
print('Number of IFC residues in protein1 and protein2:',len(p1_ifcres_indexlist),len(p2_ifcres_indexlist))

#--- READ SurfaceRes information; Generate fake_IFCs---------------------------#
with open(asrfile, "r") as srfile:
    p1_sr_indexlist, p2_sr_indexlist =srfile.readlines()[0:2]
    p1_sr_indexlist = p1_sr_indexlist.strip('\n').split(',')
    p2_sr_indexlist = p2_sr_indexlist.strip('\n').split(',')
    if p1_sr_indexlist[0][0].isdigit():
        p1_allsr_indexlist = [int(x[0])-firstindex_p1 for x in p1_sr_indexlist]#Format:12,45
        p2_allsr_indexlist = [int(x[0])-firstindex_p2 for x in p2_sr_indexlist]
    else:
        p1_allsr_indexlist = [int(x[1:])-firstindex_p1 for x in p1_sr_indexlist]#Format:D12,R45
        p2_allsr_indexlist = [int(x[1:])-firstindex_p2 for x in p2_sr_indexlist]
p1_allsr_indexlist=list(set(p1_allsr_indexlist)); p2_allsr_indexlist=list(set(p2_allsr_indexlist))
print('Number of Surface Residues:',len(p1_sr_indexlist),len(p2_sr_indexlist))
p1_sr_indexlist = [item for item in p1_allsr_indexlist if item not in p1_ifcres_indexlist]
p2_sr_indexlist = [item for item in p2_allsr_indexlist if item not in p2_ifcres_indexlist]
print('Number of non-interfacial Surface Residues:',len(p1_sr_indexlist),len(p2_sr_indexlist))

if len(p1_ifcres_indexlist) < len(p1_sr_indexlist): pass
else:
    print("The # IFR_P1 is more than nonIFR_SurfaceRes:",len(p1_ifcres_indexlist),len(p1_sr_indexlist))
    p1_sr_indexlist = list(range(0,Num_columns_p1))

if len(p2_ifcres_indexlist) < len(p2_sr_indexlist): pass
else:
    print("The # IFR_P2 is more than nonIFR_SurfaceRes:",len(p2_ifcres_indexlist),len(p2_sr_indexlist))
    p2_sr_indexlist = list(range(0,Num_columns_p2))

if homo:
    p12_fifr_lists=[]
    p12_sr_indexlist = list(set(p1_sr_indexlist + p2_sr_indexlist))
    p12_ifcres_indexlist = list(set(p1_ifcres_indexlist + p2_ifcres_indexlist))
    for i in range(0,100):
        p12_fifr_lists.append(random.sample(p12_sr_indexlist,len(p12_ifcres_indexlist)))
else:
    p1_fifr_lists=[];p2_fifr_lists=[]
    for i in range(0,100):
        p1_fifr_lists.append(random.sample(p1_sr_indexlist,len(p1_ifcres_indexlist)))
        p2_fifr_lists.append(random.sample(p2_sr_indexlist,len(p2_ifcres_indexlist)))

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
    p12_ifcres_indexlist = list(set(p1_ifcres_indexlist + p2_ifcres_indexlist))
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


# %%
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

#--- concatenate MSA for P1-P2 ----#
def contanate_seq12_allMSA12(MSA1,MSA2,use_row):
    out_seq12={}
    for k in use_row:
        out_seq12.update({k:str(MSA1[k])+str(MSA2[k])})
    return list(out_seq12.values())

if homo:
    sequences_p1=[]
    for k in common_species_cov:
        sequences_p1.append(species_proteinseq1[k])
else: sequences_p12 = contanate_seq12_allMSA12(species_proteinseq1,species_proteinseq2,common_species_cov)


# %%
# convert seqs to array
#msa_numpy_p12 = gremlin_tf_v2.mk_msa(sequences_p12)
alphabet ='ACDEFGHIKLMNPQRSTVWY-'
a2n = {}
for a,n in zip(alphabet,range(21)):
  a2n[a] = n
def aa2num(aa):
  if aa in a2n: return a2n[aa]
  else: return a2n['-']
def get_eff(msa,eff_cutoff=0.8):
  '''compute effective weight for each sequence'''
  msa_sm = 1.0 - squareform(pdist(msa,"hamming"))
  msa_w = (msa_sm >= eff_cutoff).astype(np.float32)
  msa_w = 1/np.sum(msa_w,-1)
  return msa_w

def mk_msa(seqs):
  msa_ori = []
  for seq in seqs:
    msa_ori.append([aa2num(aa) for aa in seq])
  msa_ori = np.array(msa_ori)
  msa = msa_ori
  msa_weights = get_eff(msa,0.8)
  #msa_weights=np.ones(len(msa))
  ncol = msa.shape[1] # length of sequence
  #w_idx=[]
  #for i in range(0,ncol1):
  #  for j in range(ncol1,ncol1+ncol2):
  #    w_idx.append((i,j))
  #w_idx = np.array(w_idx)
  return {"msa_ori":msa_ori,
          "msa":msa,
          "weights":msa_weights,
          "neff":np.sum(msa_weights),
  #        "w_idx":w_idx,
          "nrow":msa.shape[0],
          "ncol":ncol,
          "ncol_ori":msa_ori.shape[1]}

if homo: msa_numpy = mk_msa(sequences_p1)
else: msa_numpy = mk_msa(sequences_p12)


# %%

@nb.njit()
def get_Pi(input_msa,msa_w):
    nrow=input_msa.shape[0]
    ncol=input_msa.shape[1]
    Pi=np.zeros((ncol,21),dtype=np.float32)
    for m in range(nrow):
        for i in range(ncol):
            Pi[i,input_msa[m,i]]=Pi[i,input_msa[m,i]] + msa_w[m]
    Pi=Pi/np.sum(msa_w)
    return Pi

@nb.njit()
def get_Pij(input_msa,msa_w):
    nrow=input_msa.shape[0]
    ncol=input_msa.shape[1]
    Pij=np.zeros((ncol,ncol,21,21),dtype=np.float32) #row m, col i & j
    for m in range(nrow):
        for i in range(ncol):
            for j in range(i+1,ncol):
                Pij[i,j,input_msa[m,i],input_msa[m,j]] = Pij[i,j,input_msa[m,i],input_msa[m,j]] + msa_w[m]
    Pij = Pij/np.sum(msa_w)
    return Pij

def Compute_Weighted_Frequencies(msa_num_dict):
    ncol,weights=msa_num_dict['ncol'],msa_num_dict['weights']
    msa=msa_num_dict['msa']
    Pi_w = get_Pi(msa,weights)
    Pij_w = get_Pij(msa,weights)
    Pij_w = Pij_w + np.transpose(Pij_w,axes=[1,0,3,2])
    for i in range(ncol):
        for a in range(21):
            Pij_w[i,i,a,a]=Pi_w[i,a]
    return Pi_w,Pij_w

#Fi_w,Fij_w = Compute_Weighted_Frequencies(msa_numpy_p11)


# %%
# Refer to PNAS 2011 paper by Morcos et al.
end0 = time.time()

def compute_mu(i,j,W_mf,Pi,epsilon=1e-4):
    diff = 1.0;
    mui, muj = np.ones(21)/21, np.ones(21)/21
    pi, pj = Pi[i,:], Pi[j,:]
    while (diff > epsilon):
        newi = pi/np.dot(muj, W_mf[i,j,:,:].T)
        newi = newi/sum(newi)
        newj = pj/np.dot(mui, W_mf[i,j,:,:])
        newj = newj/sum(newj)
        diffi=np.max(np.abs(newi-mui))
        diffj=np.max(np.abs(newj-muj))
        diff=max(diffi,diffj)
        mui = np.copy(newi)
        muj = np.copy(newj)
    return mui,muj

def mfDCA_p1p2(msa_dict,ncol1,ncol2,alpha=0.5):
    #alpha: psedo_weight = lamda/(lamda + Meff), Eq.1
    if homo: ncol=ncol1; print("Calculating DCA for Homodimer with length:",ncol1)
    else: ncol = ncol1+ncol2

    # computes weighted frequencies, Eq.1
    Freq_i,Freq_ij = Compute_Weighted_Frequencies(msa_dict)

    # adds pseudocount, Eq.1
    Pi = (1.0-alpha)*Freq_i + (alpha/21)*np.ones((ncol,21),dtype=np.float32)  
    Pij = (1.0-alpha)*Freq_ij + (alpha/21/21)*np.ones((ncol,ncol,21,21),dtype=np.float32) 
    for i in range(ncol):
        Pij[i,i,:,:] = (1.0 - alpha)*Freq_ij[i,i,:,:] + (alpha/21)*np.eye(21)
    # computes correlation/covariance matrix and its inverse (coupling matrix), Eq.9-10
    C = np.float32(Pij-Pi[:,None,:,None]*Pi[None,:,None,:])
    C = C[:,:,0:20,0:20]
    C = C.transpose([0,2,1,3])
    C = np.reshape(C,(ncol*(20),ncol*(20)))
    invC = np.linalg.inv(C)

    # computes (weighted) entropy and mutual information
    tiny = 1e-6
    Si = -np.sum(Freq_i * np.log(Freq_i+tiny),1)
    Sij = -np.sum(Freq_ij * np.log(Freq_ij+tiny),(2,3))
    MI=Si[None,:]+ Si[:,None] - Sij
    if not homo:
        MI=MI[:ncol1,ncol1:]

    # computes mean field W, Eq.8
    W_mf = np.ones((ncol,21,ncol,21),dtype=np.float32)
    W_mf[:,0:20,:,0:20] = np.exp(-invC.reshape(ncol,20,ncol,20))
    W_mf = np.transpose(W_mf,axes=[0,2,1,3])

    # computes direct information between p1 and p2
    if not homo:
        DI = np.zeros((ncol1,ncol2),dtype=np.float32)
        for i in range(ncol1):
            for j in range(ncol1,ncol1+ncol2):
            # computes hi & hj, Eq.12
                mui, muj = compute_mu(i,j,W_mf,Pi)
            # Eq.11
                Pdir = W_mf[i,j,:,:] * (np.tensordot(mui,muj,0))    
                Pdir = Pdir/np.sum(Pdir)
            # Eq.13
                Pfac = np.tensordot(Pi[i,:].T,Pi[j,:],0)
                DI[i,j-ncol1] = np.trace(np.dot(Pdir.T, np.log((Pdir+tiny)/(Pfac+tiny))))
    else:
        DI = np.zeros((ncol1,ncol1),dtype=np.float32)
        for i in range(ncol1):
            for j in range(i+1,ncol1):
            # computes hi & hj, Eq.12
                mui, muj = compute_mu(i,j,W_mf,Pi)
            # Eq.11
                Pdir = W_mf[i,j,:,:] * (np.tensordot(mui,muj,0))    
                Pdir = Pdir/np.sum(Pdir)
            # Eq.13
                Pfac = np.tensordot(Pi[i,:].T,Pi[j,:],0)
                DI[i,j] = np.trace(np.dot(Pdir.T, np.log((Pdir+tiny)/(Pfac+tiny))))
        DI = DI + DI.T
    end3 = time.time(); #print ("Time till mfDCA: ", end3-end0)
    out = {"di": DI,}
#        "mi":MI,} 
    return out
   
dc_out = mfDCA_p1p2(msa_numpy,Num_columns_p1,Num_columns_p2)

# %%
def checkTruth_triu_tril(mat):
    return np.allclose(mat, np.triu(mat)), np.allclose(mat, np.tril(mat))

# %%
def get_apc(raw_p12):
    raw_sq=raw_p12
    t=checkTruth_triu_tril(raw_sq)
    if t[0] and not t[1]:
        raw_sq = raw_sq + raw_sq.T
    elif not t[0] and not t[1]:
        raw_sq = raw_sq
    ap_sq = np.sum(raw_sq,0,keepdims=True)*np.sum(raw_sq,1,keepdims=True)/np.sum(raw_sq)
    apc = raw_sq - ap_sq
    out = {
        "raw": raw_sq,
        "apc": apc,}
    return out
    
di_apc_out=get_apc(dc_out['di'])
di_raw_p1p2 = di_apc_out['raw']
di_apc_p1p2 = di_apc_out['apc']

end1 = time.time()
print("Time used for DI & APC: %.4f" %(end1-end0))

# %%
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
    real_fake = np.append([reference],inputmatrix,axis=0)
    column_mean = np.nanmean(real_fake, axis=0);
    column_std = np.nanstd(real_fake, axis=0) + 1e-10;
    try:
        column_zscore = ((reference)-(column_mean))/(column_std)
    except (ValueError, ZeroDivisionError):
        print('ValueError or ZeroDivisionError')
        pass
    else:
        return column_zscore

# %%
matr_toStat=[di_raw_p1p2,di_apc_p1p2]

NM_maxmean_list=[];IFC_maxmean_list=[]
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

end2 = time.time(); 
print("Time used for Z-score: %.4f" %(end2-end1))

# %%

with open(outfile, 'w') as out:
    out.write("# Protein1_MSA and Protein2_MSA: (%i, %i), (%i, %i).\n" % (Num_rows_p1,Num_columns_p1,Num_rows_p2,Num_columns_p2))
    out.write("# %i IFCs (%i,%i); %i common_species, %i common_cov, %.2f Neff.\n" % (len(rp_fullname_list),Nifc[0],Nifc[1],len(common_species),len(common_species_cov),msa_numpy['neff']))
    out.write("# P1_MSA_IFC_SR and P2_MSA_IFC_SR: (%i, %i, %i), (%i, %i, %i).\n" % (len(common_species),len(p1_ifcres_indexlist),len(p1_sr_indexlist),len(common_species),len(p2_ifcres_indexlist),len(p2_sr_indexlist)))
    out.write("# Header in order of: DI, DI_APC\n")
    out.write("# Max in NMmatrix: %.4f\t%.4f\n" %(NM_maxmean_list[0][0],NM_maxmean_list[1][0]))
    out.write("# Mean in NMmatrx: %.4f\t%.4f\n" %(NM_maxmean_list[0][1],NM_maxmean_list[1][1]))
    out.write("# Max in IFCpairs: %.4f\t%.4f\n" %(IFC_maxmean_list[0][0],IFC_maxmean_list[1][0]))
    out.write("# Mean in IFCpair: %.4f\t%.4f\n" %(IFC_maxmean_list[0][1],IFC_maxmean_list[1][1]))

    out.write("# Z of max_in__NM: %.4f\t%.4f\n" %(Zscore_NM_list[0][0],Zscore_NM_list[1][0]))
    out.write("# Z of mean_in_NM: %.4f\t%.4f\n" %(Zscore_NM_list[0][1],Zscore_NM_list[1][1]))
    out.write("# Z of max_in_IFC: %.4f\t%.4f\n" %(Zscore_IFC_list[0][0],Zscore_IFC_list[1][0]))
    out.write("# Z of meanin_IFC: %.4f\t%.4f\n" %(Zscore_IFC_list[0][1],Zscore_IFC_list[1][1]))
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
            di,di_apc = di_raw_p1p2[index_a][index_b], di_apc_p1p2[index_a][index_b]
            out.write("%s\t%s\t%.4f\t%.4f\n" %(a_fullname,b_fullname,di,di_apc))
"""
Example output:
# Protein1_MSA and Protein2_MSA: (286, 161), (220, 161).
# 132 IFCs (46,40); 208 common_species, 206 common_cov, 66.55 Neff.
# P1_MSA_IFC_SR and P2_MSA_IFC_SR: (208, 46, 100), (208, 40, 108).
# Header in order of: DI, DI_APC
# Max in NMmatrix: 0.4139   0.3963
# Mean in NMmatrx: 0.0058   0.0027
# Max in IFCpairs: 0.0070   0.0067
# Mean in IFCpair: 0.0020   0.0001
# Z of max_in__NM: 0.2831   0.2607
# Z of mean_in_NM: -3.1235  5.6344
# Z of max_in_IFC: -1.5256  -1.4694
# Z of meanin_IFC: -2.8269  -0.4286
"""


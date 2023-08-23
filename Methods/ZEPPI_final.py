import pandas as pd
import sys,os

dca_csvfile=False
if len(sys.argv) < 3:
    print ("\n python ZEPPI_final.py miz_csvfile (dca_csvfile) output \n")
    exit()
elif len(sys.argv) == 3 and os.path.isfile(sys.argv[1]):
    miz_csvfile=sys.argv[1]
    #print ("\nFinalizing ZEPPI based on Z-score of MI and Con.")
elif len(sys.argv) == 4 and os.path.isfile(sys.argv[1]) and os.path.isfile(sys.argv[2]):
    #print ("\nFinalizing ZEPPI based on Z-score of MI, Con and DCA.")
    miz_csvfile=sys.argv[1]
    dca_csvfile=sys.argv[2]

outputfile=sys.argv[-1]

def read_MI_csvfile(MIZcsv):
    df1 = pd.read_csv(MIZcsv,sep=',',skiprows=0).drop_duplicates()
    df1 = df1[['PPI','Neff','Nifr','Nifr_eff','Z_Mean_IFC_ar1','Z_Mean_IFC_ar2','Z_Mean_IFC_ar3','Z_Mean_IFC_ar4',\
    'Z_Max_IFC_ar1','Z_Max_IFC_ar2','Z_Max_IFC_ar3','Z_Max_IFC_ar4',]]
    df3= df1.dropna()
    print('Number of PPIs in MIZ file:',len(df1.index))
    df3['Ratio'] = (df3['Neff']/df3['Nifr'])
    df3['Zmean_MI']=df3[['Z_Mean_IFC_ar1','Z_Mean_IFC_ar2',]].max(axis=1)
    df3['Zmean_Con']=df3[['Z_Mean_IFC_ar3','Z_Mean_IFC_ar4']].max(axis=1)
    df3['Ztop_MI']=df3[['Z_Max_IFC_ar1','Z_Max_IFC_ar2',]].max(axis=1)
    df3['Ztop_Con']=df3[['Z_Max_IFC_ar3','Z_Max_IFC_ar4']].max(axis=1)
    df3['Zmaxm'] =  df3[['Zmean_MI','Zmean_Con']].max(axis=1); df3['Zmaxt'] = df3[['Ztop_MI','Ztop_Con']].max(axis=1);
    df3['ZEPPI'] = df3[['Zmaxm','Zmaxt']].max(axis=1)
    return df3

def read_MI_DCA_csvfiles(MIcsv,DCAcsv):
    df1 = pd.read_csv(MIcsv,sep=',',skiprows=0)
    df2 = pd.read_csv(DCAcsv,sep=',',skiprows=0)
    df1 = df1[['PPI','Nmsa','Neff','Nifr','Nifr_eff','Mean_IFC_ar1','Mean_IFC_ar2','Mean_IFC_ar3','Mean_IFC_ar4',\
    'Z_Mean_IFC_ar1','Z_Mean_IFC_ar2','Z_Mean_IFC_ar3','Z_Mean_IFC_ar4','Z_Max_IFC_ar1','Z_Max_IFC_ar2','Z_Max_IFC_ar3','Z_Max_IFC_ar4',]]
    df2 =df2.rename(columns={"Z_Mean_IFC_ar1": "Z_Mean_IFC_ar5","Z_Mean_IFC_ar2": "Z_Mean_IFC_ar6"})
    df2 =df2.rename(columns={"Z_Max_IFC_ar1": "Z_Max_IFC_ar5","Z_Max_IFC_ar2": "Z_Max_IFC_ar6"})
    df2=df2[["PPI","Z_Mean_IFC_ar5","Z_Mean_IFC_ar6","Z_Max_IFC_ar5","Z_Max_IFC_ar6"]]
    df3 = pd.merge(df1, df2, on="PPI", how="left")
    df3= df3.dropna()
    print('Number of PPIs in MIZ and DCA file:',len(df1.index),len(df2.index))
    df3['Zmean_MI']=df3[['Z_Mean_IFC_ar1','Z_Mean_IFC_ar2',]].max(axis=1)
    df3['Zmean_Con']=df3[['Z_Mean_IFC_ar3','Z_Mean_IFC_ar4']].max(axis=1)
    df3['Zmean_DCA']=df3[['Z_Mean_IFC_ar5','Z_Mean_IFC_ar6']].max(axis=1)
    df3['Ztop_MI']=df3[['Z_Max_IFC_ar1','Z_Max_IFC_ar2',]].max(axis=1)
    df3['Ztop_Con']=df3[['Z_Max_IFC_ar3','Z_Max_IFC_ar4']].max(axis=1)
    df3['Ztop_DCA']=df3[['Z_Max_IFC_ar5','Z_Max_IFC_ar6']].max(axis=1)
    df3['Ratio'] = (df3['Neff']/df3['Nifr']) 
    df3['ZEPPI'] = df3[['Zmaxm','Zmaxt']].max(axis=1)
    return df3


if not dca_csvfile:
    df1 = read_MI_csvfile(miz_csvfile)
else:
    df1 = read_MI_DCA_csvfiles(miz_csvfile,dca_csvfile)

df1=df1[['PPI','Neff','Nifr','Zmean_MI','Zmean_Con','Ztop_MI','Ztop_Con','ZEPPI',]]
    
df1.to_csv(outputfile, float_format='%.4f',index=False)


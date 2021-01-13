import numpy as np
from sklearn.neighbors.kde import KernelDensity
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numpy import array, linspace
import pandas as pd
from scipy.signal import argrelextrema

cov = pd.read_csv(input_csv,index_col = 0)

colors = ['r','g','b','c','m','y']

thres_dic = {}
for i in range(len(cov)):
    gene = cov.index[i]
    a = array(cov.iloc[i,:]).reshape(-1,1)
    kde = KernelDensity(kernel='gaussian', bandwidth=0.075).fit(a)
    s = linspace(min(a)-0.05,max(a)+0.05)
    e = kde.score_samples(s.reshape(-1,1))
    mi, ma = argrelextrema(e,np.less)[0], argrelextrema(e,np.greater)[0]
    n_thress = len(mi)
    n_bin = len(ma)
    thres_dic[gene] = s[mi]
    if n_thress > 0:
        plt.figure()
        plt.plot(s[:mi[0]+1],e[:mi[0]+1], colors[0])
        for j in range(1,n_bin-1):
            plt.plot(s[mi[j-1]:mi[j]+1],e[mi[j-1]:mi[j]+1], colors[j])
        plt.plot(s[mi[n_bin-2]:],e[mi[n_bin-2]:],colors[n_bin-1])
        plt.savefig('process20200828/gene_copy/fig/combined.'+gene+'.png')
    else:
        print(gene+' had zero threshold')

genelist = ['KIR3DS1','KIR3DL1','KIR2DS4','KIR2DS3;KIR2DS5','KIR2DS2','KIR2DS1','KIR2DP1','KIR2DL5A;KIR2DL5B','KIR2DL3','KIR2DL2','KIR2DL1','KIR3DL3','KIR3DL2','KIR2DL4']

tmp = np.zeros((len(genelist),len(cov.columns)),dtype = "int")
copy = pd.DataFrame(tmp)
copy.columns = cov.columns
copy.index = genelist

for i in range(len(genelist)-1):
    gene = genelist[i]
    cut_thres = np.hstack(([0],np.ravel(thres_dic[gene]),[4]))
    ratio = cov.loc[gene,:]
    copy.iloc[i,:] = np.array(pd.cut(ratio, cut_thres, labels=False))

# for KIR2DL4
i=len(genelist)-1
gene = 'KIR2DL4'
cut_thres = np.hstack(([0],np.ravel(thres_dic[gene]),[4]))
labels = np.array([1,2,3,4])
ratio = cov.loc[gene,:]
copy.iloc[i,:] = np.array(pd.cut(ratio, cut_thres, labels=labels))

copy.to_csv(output_ploidy_csv)

OUT = open(output_ploidy_file, "w")
for i in range(len(copy.columns)):
    sample = copy.columns[i]
    for j in range(len(copy.index)):
        gene = copy.index[j]
        n_copy = copy.iloc[j,i]
        if n_copy > 0 :
            out = sample + "\t" + gene + "\t" + str(n_copy)
            print >> OUT, out
OUT.close()
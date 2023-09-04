import numpy as np
import pandas as pd
import ADlasso
from ADlasso.AD_utils import *
import time
import os
import scipy
outPath = '/home/yincheng23/ADlasso_manuscript/Datasets/data/ADlass_results/'
dataSetsPath = '/home/yincheng23/ADlasso_manuscript/Datasets/Studies'
os.chdir(dataSetsPath)
dataSets = sorted(os.listdir(dataSetsPath))
bigdatasets = ['ArcticTransects','GWMC_ASIA_NA','GWMC_HOT_COLD','ob_goodrich']
for dataSet in dataSets:
    if dataSet not in bigdatasets:
        continue
    os.chdir(dataSetsPath)
    os.chdir(dataSet)
    
    Data = pd.read_csv("ASV_vst.txt",sep = "\t")
    Data = Data.T
    Data = Data.sort_index(axis = 0)
    Data_std = scipy.stats.zscore(Data, axis=0,ddof=0)
    RawData = pd.read_csv("ASV_table.txt",sep = "\t")
    RawData = RawData.T
    RawData = RawData.loc[Data.index]
    Cohort = pd.read_csv("metadata.txt",sep = "\t")
    Cohort = Cohort.sort_index(axis = 0)
    
    label = Cohort['Class'].tolist()
    pvl = get_prevalence(RawData, np.arange(RawData.shape[0]))
    
    start = time.time()
    lmbd_range = auto_scale(Data_std, RawData, label, step=50, device='cuda', max_iter=100000)
    lmbd_rng = np.exp(lmbd_range)
    outpath_dataset = outPath + dataSet
    result_dict =lambda_tuning(Data_std, RawData, label, lmbd_rng, 5, outpath_dataset,device='cuda', max_iter=100000)
    end = time.time()
    print(dataSet, ", total cost：%f min" % round((end - start)/60,3))

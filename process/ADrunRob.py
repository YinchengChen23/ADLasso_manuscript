import numpy as np
import pandas as pd
import ADlasso as AD
from ADlasso.AD_utils import *
import time
import os
import scipy

dataSetsPath = '/home/yincheng23/ADlasso_manuscript/Datasets/Testing_Bias_robustness/Diarrhea'
outPath = '/home/yincheng23/ADlasso_manuscript/Datasets/data/ADlass_results_robustness/Diarrhea_'

os.chdir(dataSetsPath)
dataSets = os.listdir(dataSetsPath)

for dataSet in dataSets:
    print(dataSet)
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

    k_fold = 5
    start = time.time()
    lmbd_range = auto_scale(Data_std, RawData, label, step=50, training_echo=True, max_iter=100000)
    lmbd_rng = np.exp(lmbd_range)
    outPath_dataSet = outPath + dataSet
    result_dict =lambda_tuning(Data_std, RawData, label, lmbd_rng, k_fold, outPath_dataSet, max_iter=100000)
    end = time.time()
    print(dataSet, ", total cost：%f min" % round((end - start)/60,3))


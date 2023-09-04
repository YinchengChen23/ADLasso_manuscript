import numpy as np
import pandas as pd
import ADlasso as AD
from ADlasso.AD_utils import *
import time
import os
import scipy

pathlist = '/home/yincheng23/ADlasso_manuscript/Datasets/LIBSVM'

outpathlist = '/home/yincheng23/ADlasso_manuscript/Datasets/data/ADlass_results/LIBSVM_rs'

X, y = get_data("/home/yincheng23/ADlasso_manuscript/Datasets/LIBSVM/real-sim")
pvl = get_prevalence(X, np.arange(X.shape[0]))

print('start to ADlasso cpu')
start = time.time()
res = ADlasso(lmbd = 6.951927961775601e-08,device='cpu', echo=True)
res.fit(X, y, pvl)
end = time.time()
print('total selected feature :',np.sum(res.feature_set))
print("Total cost：%f sec" % (end - start))

print('start to ADlasso gpu')
start = time.time()
res = ADlasso(lmbd = 6.951927961775601e-08,device='cuda', echo=True)
res.fit(X, y, pvl)
end = time.time()
print('total selected feature :',np.sum(res.feature_set))
print("Total cost：%f sec" % (end - start))

#k_fold = 5
#lmbd_range = auto_scale(X, X, y, step=20, training_echo=True, device='cuda')
#lmbd_rng = np.exp(lmbd_range)
#result_dict =lambda_tuning(X, X, y, lmbd_rng, k_fold, outpathlist, device='cuda')

#lmbd_range = auto_scale(X, X, y, step=50, training_echo=True)
#lmbd_rng = np.exp(lmbd_range)
#result_dict =lambda_tuning(X, X, y, lmbd_rng, 10, outpathlist, training_echo=True)

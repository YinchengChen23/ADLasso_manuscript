import os
import sys
import time
import scipy
import itertools
import numpy as np
import pandas as pd
import random as rnd
from random import randrange
from sklearn import metrics
from sklearn.svm import LinearSVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.calibration import CalibratedClassifierCV
from sklearn.linear_model import LogisticRegression, ElasticNet, Lasso
from sklearn.feature_selection import mutual_info_classif
from sklearn.model_selection import train_test_split
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from sklearn.metrics import accuracy_score
from sklearn.metrics import matthews_corrcoef
from sklearn.metrics import precision_recall_curve
from sklearn.metrics.pairwise import pairwise_distances
from xgboost import XGBClassifier
import mrmr
from mrmr import mrmr_classif
from warnings import simplefilter
from sklearn.exceptions import ConvergenceWarning
simplefilter("ignore", category=ConvergenceWarning)

def reliefF(X, y, **kwargs):
    if "k" not in kwargs.keys():
        k = 5
    else:
        k = kwargs["k"]
    n_samples, n_features = X.shape
    distance = pairwise_distances(X, metric='manhattan')
    score = np.zeros(n_features)
    for iter in range(n_samples):
        idx = randrange(0, n_samples, 1)
        near_hit = []
        near_miss = dict()
        self_fea = X[idx, :]
        c = np.unique(y).tolist()
        stop_dict = dict()
        for label in c:
            stop_dict[label] = 0
        del c[c.index(y[idx])]
        p_dict = dict()
        p_label_idx = float(len(y[y == y[idx]]))/float(n_samples)
        for label in c:
            p_label_c = float(len(y[y == label]))/float(n_samples)
            p_dict[label] = p_label_c/(1-p_label_idx)
            near_miss[label] = []
        distance_sort = []
        distance[idx, idx] = np.max(distance[idx, :])
        for i in range(n_samples):
            distance_sort.append([distance[idx, i], int(i), y[i]])
        distance_sort.sort(key=lambda x: x[0])
        for i in range(n_samples):
            if distance_sort[i][2] == y[idx]:
                if len(near_hit) < k:
                    near_hit.append(distance_sort[i][1])
                elif len(near_hit) == k:
                    stop_dict[y[idx]] = 1
            else:
                if len(near_miss[distance_sort[i][2]]) < k:
                    near_miss[distance_sort[i][2]].append(distance_sort[i][1])
                else:
                    if len(near_miss[distance_sort[i][2]]) == k:
                        stop_dict[distance_sort[i][2]] = 1
            stop = True
            for (key, value) in stop_dict.items():
                    if value != 1:
                        stop = False
            if stop:
                break
        near_hit_term = np.zeros(n_features)
        for ele in near_hit:
            near_hit_term = np.array(abs(self_fea-X[ele, :]))+np.array(near_hit_term)
        near_miss_term = dict()
        for (label, miss_list) in near_miss.items():
            near_miss_term[label] = np.zeros(n_features)
            for ele in miss_list:
                near_miss_term[label] = np.array(abs(self_fea-X[ele, :]))+np.array(near_miss_term[label])
            score += near_miss_term[label]/(k*p_dict[label])
        score -= near_hit_term/k
    return score

def fisherScore(X, y, k):
    u1 = np.mean(X[np.where(y == 0),:], axis=1)
    u2 = np.mean(X[np.where(y == 1),:], axis=1)
    v1 = np.var(X[np.where(y == 0),:], axis=1)
    v2 = np.var(X[np.where(y == 1),:], axis=1)
    fisherS = ((u2 - u1)**2) / (v2 - v1)
    fisherS = fisherS.reshape(-1)
    selected = np.argsort(fisherS)[0:k]
    return selected

def FDC(X, k):
    n_sample, n_feature = X.shape
    FD_ = []
    for i in range(n_feature):
        FD_.append(np.log(np.sum(np.exp(X[:,i]))) - np.mean(X[:,i]))
    selected = np.argsort(FD_)[::-1][0:k]
    return selected

def create_grid(params):
    comb = []
    for t in itertools.product(*params):
        comb.append(t)
    return comb

def writeList(featureSet, featureName, weight, method, path):
    outpath = path + '/' + method + '.txt'
    w = open(outpath,'w')
    for i in featureSet:
        w.writelines(featureName[i] + "\t" + str(weight[i]) + '\n')
    w.close()

def writeList_No(featureSet, featureName, method, path):
    outpath = path + '/' + method + '.txt'
    w = open(outpath,'w')
    for i in featureSet:
        w.writelines(featureName[i] + '\n')
    w.close()

def getParaSet():
    MLmethods = {'LA':'LASSO', 'SV':'SVM-LASSO', 'EN':'EN','RF':'RF', 'XG':'XGBoost', 'MI':'MI','re':'reliefF'}
    ParaDist = {}
    parmfile = open("parm_info.txt", "r")
    parmlines = parmfile.readlines()
    for parmline in parmlines:
        parmline = parmline.split('\n')[0]
        if parmline[0:4] == 'Best' :
            method = parmline[25:27]
            condition_ = parmline.split('is: (')[1]
            condition = condition_.split(') with prevalence')[0]
            condition = condition.replace(' ','')
            if condition[-1] == ',':
                condition = condition[:-1]
            ParaDist[MLmethods[method]] = condition.split(',')
    parmfile.close()
    return ParaDist

def assigntParaSet(bestPara):
    parameters = {'LASSO'     : ['C'],
              'SVM-LASSO' : ['C'],
              'EN'        : ['alpha', 'l1_ratio'],
              'RF'        : ['n_estimators', 'max_depth'],
              'XGBoost'   : ['learning_rate', 'n_estimators', 'max_depth', 'reg_alpha', 'reg_lambda'],
              'MI'        : ['n_neighbors'],
              'reliefF'   : ['k']}
    final_parameters = {}
    for ml in bestPara:
        subPara = {}
        for ix, para in enumerate(bestPara[ml]):
            if parameters[ml][ix] in ['n_estimators', 'max_depth', 'n_neighbors', 'k']:
                subPara[parameters[ml][ix]] = int(para)
            else:
                subPara[parameters[ml][ix]] = float(para)
        final_parameters[ml] = subPara
    return final_parameters

class FeatureSelector:
    def __init__(self, model=None, name=None, params=None, max_num_feat=None):
        self.name = name
        self.model = model
        self.params = params
        self.weight = None
        self.selectedset = None
        self.sortedImp = None
        self.max_num_feat = max_num_feat

    def fit(self, Data, X, y, outdir):
        if self.name == 'LASSO':
            classifier = LogisticRegression(C=self.params['C'], penalty="l1", solver='liblinear').fit(X, y)
            self.model = classifier
            weight = classifier.coef_[0]
            fullset = [i for i, w in enumerate(weight) if w != 0]
            weight_abs = np.abs(weight)
            limitedset = np.argsort(-weight_abs)[0:self.max_num_feat]
            writeList(fullset, Data.columns.tolist(), weight, 'LASSO', outdir)
            writeList(limitedset, Data.columns.tolist(), weight, 'LASSO_limited', outdir)
        if self.name == 'SVM-LASSO':
            classifier = LinearSVC(C=self.params['C'], penalty="l1", dual=False).fit(X, y)
            weight = classifier.coef_[0]
            fullset = [i for i, w in enumerate(weight) if w != 0]
            weight_abs = np.abs(weight)
            limitedset = np.argsort(-weight_abs)[0:self.max_num_feat]
            writeList(fullset, Data.columns.tolist(), weight, 'SVMLASSO', outdir)
            writeList(limitedset, Data.columns.tolist(), weight, 'SVMLASSO_limited', outdir)
        if self.name == 'EN':
            classifier = ElasticNet(alpha=self.params['alpha'], l1_ratio=self.params['l1_ratio']).fit(X, y)
            self.model = classifier
            weight = classifier.coef_
            fullset = [i for i, w in enumerate(weight) if w != 0]
            weight_abs = np.abs(weight)
            limitedset = np.argsort(-weight_abs)[0:self.max_num_feat]
            writeList(fullset, Data.columns.tolist(), weight, 'EN', outdir)
            writeList(limitedset, Data.columns.tolist(), weight, 'EN_limited', outdir)
        if self.name == 'RF':
            classifier = RandomForestClassifier(n_estimators=self.params['n_estimators'], max_depth=self.params['max_depth']).fit(X, y)
            self.model = classifier
            weight = classifier.feature_importances_
            fullset = [i for i, w in enumerate(weight) if w != 0]
            weight_abs = np.abs(weight)
            limitedset = np.argsort(-weight_abs)[0:self.max_num_feat]
            writeList(fullset, Data.columns.tolist(), weight, 'RF', outdir)
            writeList(limitedset, Data.columns.tolist(), weight, 'RF_limited', outdir)
        if self.name == 'XGBoost':
            classifier = XGBClassifier(learning_rate=self.params['learning_rate'], booster = 'gbtree', use_label_encoder = False, n_jobs = 20,
                                       n_estimators=self.params['n_estimators'], max_depth=self.params['max_depth'],
                                       reg_alpha=self.params['reg_alpha'], reg_lambda=self.params['reg_lambda']).fit(X, y, eval_metric = 'logloss')
            self.model = classifier
            weight = classifier.feature_importances_
            fullset = [i for i, w in enumerate(weight) if w != 0]
            weight_abs = np.abs(weight)
            limitedset = np.argsort(-weight_abs)[0:self.max_num_feat]
            writeList(fullset, Data.columns.tolist(), weight, 'XGBoost', outdir)
            writeList(limitedset, Data.columns.tolist(), weight, 'XGBoost_limited', outdir)
        if self.name == 'MI':
            weight = mutual_info_classif(X, y, n_neighbors=self.params['n_neighbors'])
            sordedidx = np.argsort(weight)[::-1]
            featureset = sordedidx[0:self.max_num_feat]
            writeList(featureset, Data.columns.tolist(), weight, 'MI', outdir)
        if self.name == 'reliefF':
            weight = reliefF(X, y, k=self.params['k'])
            sordedidx = np.argsort(weight)[::-1]
            featureset = sordedidx[0:self.max_num_feat]
            writeList(featureset, Data.columns.tolist(), weight, 'ReliefF', outdir)
        if self.name == 'mRMR':
            X_ = pd.DataFrame(X)
            y_ = pd.Series(y)
            featureset = mrmr_classif(X=X_, y=y_, K=self.max_num_feat)
            writeList_No(featureset, Data.columns.tolist(), 'mRMR', outdir)
        if self.name == 'fisher':
            featureset = fisherScore(X, y, k=self.max_num_feat)
            writeList_No(featureset, Data.columns.tolist(), 'fisher', outdir)
        if self.name == 'FD':
            featureset = FDC(X, k=self.max_num_feat)
            writeList_No(featureset, Data.columns.tolist(), 'FDC', outdir)


MLmethods = {'LA':'LASSO', 'SV':'SVM-LASSO', 'EN':'EN','RF':'RF', 'XG':'XGBoost', 'MI':'MI','re':'reliefF'}
dataSetsPath = '/home/yincheng23/ADlasso_manuscript/Datasets/Testing_Bias_robustness/'
outPath = '/home/yincheng23/ADlasso_manuscript/Datasets/data/ML_parm_info_rob/'
refPath = '/home/yincheng23/ADlasso_manuscript/Datasets/data/selected_results_robustness/'
paramsPath = '/home/yincheng23/ADlasso_manuscript/Datasets/data/ML_parm_info_rob/'
os.chdir(refPath)
dataSets = sorted(os.listdir())


for dataSet in dataSets:
    #======================== Data loading=============================================
    ach1 = dataSet.split('_')[0]; ach2 = '_'.join(dataSet.split('_')[1:])
    os.chdir(dataSetsPath); os.chdir(ach1); os.chdir(ach2)
    print(dataSet)

    Data = pd.read_csv("ASV_vst.txt",sep = "\t")
    Data = Data.T
    Data = Data.sort_index(axis = 0)
    Data_std = scipy.stats.zscore(Data, axis=0,ddof=0)

    RawData = pd.read_csv("ASV_table.txt",sep = "\t")
    RawData = RawData.T
    RawData = RawData.loc[Data.index]
    RawData[RawData > 0] = 1
    prevalence = [len(RawData.index[RawData.iloc[:,i] > 0])*100/RawData.shape[0] for i in range(RawData.shape[1])]

    Cohort = pd.read_csv("metadata.txt",sep = "\t")
    Cohort = Cohort.sort_index(axis = 0)
    Label = Cohort['Class'].tolist()

    X = Data_std.to_numpy()
    n_samples, n_featunre = X.shape
    anchor = Label[0]; y = [0 if yi == anchor else 1 for yi in Label]
    y = np.array(y); y = y.reshape(-1)

    refile = refPath + dataSet + '/ADlasso.txt'
    if os.path.isfile(refile):
        refile = pd.read_csv(refile,sep = "\t", header = None)
    else:
        continue
    max_num = refile.shape[0]
    print(dataSet, max_num)
    #=================================================================================
    os.chdir(paramsPath)
    os.chdir(dataSet)
    BestPara = getParaSet()
    BestPara = assigntParaSet(BestPara)
    models_FS = {'LASSO'   : FeatureSelector(name='LASSO', params=BestPara['LASSO'], max_num_feat=max_num),
                 'SVM-LASSO' : FeatureSelector(name='SVM-LASSO', params=BestPara['SVM-LASSO'], max_num_feat=max_num),
                 'EN'      : FeatureSelector(name='EN', params=BestPara['EN'], max_num_feat=max_num),
                 'RF'      : FeatureSelector(name='RF', params=BestPara['RF'], max_num_feat=max_num),
                 'XGBoost' : FeatureSelector(name='XGBoost', params=BestPara['XGBoost'], max_num_feat=max_num),
                 'MI'      : FeatureSelector(name='MI', params=BestPara['MI'], max_num_feat=max_num),
                 'reliefF' : FeatureSelector(name='reliefF', params=BestPara['reliefF'], max_num_feat=max_num),
                 'mRMR' : FeatureSelector(name='mRMR', max_num_feat=max_num),
                 'fisher' : FeatureSelector(name='fisher', max_num_feat=max_num),
                 'FD' : FeatureSelector(name='FD', max_num_feat=max_num)}
    outpath_ = paramsPath + dataSet
    for fs_name, fs_model in models_FS.items():
        print(fs_name)
        fs_model.fit(Data, X, y, outpath_)

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

class FeatureSelector:
    def __init__(self, model=None, name=None, params=None):
        self.name = name
        self.model = model
        self.params = params
        self.weight = None
        self.selectedset = None
        self.sortedImp = None
        self.max_num_feat = None
    def setParams(self, comb_par, params_name, params):
        for par_name, par in zip(params_name, comb_par):
            params[par_name] = par
        self.params = params
    def fit(self, X, y):
        if self.name == 'LASSO':
            classifier = LogisticRegression(C=self.params['C'], penalty="l1", solver='liblinear').fit(X, y)
            self.model = classifier
            self.weight = classifier.coef_[0]
            self.selectedset = [i for i, w in enumerate(self.weight) if w != 0]
        if self.name == 'SVM-LASSO':
            classifier = LinearSVC(C=self.params['C'], penalty="l1", dual=False).fit(X, y)
            self.weight = classifier.coef_[0]
            self.selectedset = [i for i, w in enumerate(self.weight) if w != 0]
        if self.name == 'EN':
            classifier = ElasticNet(alpha=self.params['alpha'], l1_ratio=self.params['l1_ratio']).fit(X, y)
            self.model = classifier
            self.weight = classifier.coef_
            self.selectedset = [i for i, w in enumerate(self.weight) if w != 0]
        if self.name == 'RF':
            classifier = RandomForestClassifier(n_estimators=self.params['n_estimators'], max_depth=self.params['max_depth']).fit(X, y)
            self.model = classifier
            self.weight = classifier.feature_importances_
            self.selectedset = [i for i, w in enumerate(self.weight) if w != 0]
        if self.name == 'XGBoost':
            classifier = XGBClassifier(learning_rate=self.params['learning_rate'], booster = 'gbtree', use_label_encoder = False,
                                       n_estimators=self.params['n_estimators'], max_depth=self.params['max_depth'],n_jobs=1,
                                       reg_alpha=self.params['reg_alpha'], reg_lambda=self.params['reg_lambda']).fit(X, y, eval_metric = 'logloss')
            self.model = classifier
            self.weight = classifier.feature_importances_
            self.selectedset = [i for i, w in enumerate(self.weight) if w != 0]
        if self.name == 'MI':
            self.weight = mutual_info_classif(X, y, n_neighbors=self.params['n_neighbors'])
            sordedidx = np.argsort(self.weight)[::-1]
            self.selectedset = sordedidx[0:self.max_num_feat]
        if self.name == 'reliefF':
            self.weight = reliefF(X, y, k=self.params['k'])
            sordedidx = np.argsort(self.weight)[::-1]
            self.selectedset = sordedidx[0:self.max_num_feat]
        if self.name == 'mRMR':
            X_ = pd.DataFrame(X)
            y_ = pd.Series(y)
            self.selectedset = mrmr_classif(X=X_, y=y_, K=self.max_num_feat)
        if self.name == 'fisher':
            self.selectedset = fisherScore(X, y, k=self.max_num_feat)
        if self.name == 'FD':
            self.selectedset = FDC(X, k=self.max_num_feat)

    def get_pvl(self, pvl):
        feature_idx = self.selectedset
        if len(feature_idx) == 0:
            return 0
        else:
            return np.median([pvl[i] for i in feature_idx])
    def evaluate(self, X_train, X_test, y_train, y_test):
        feature_idx = self.selectedset
        if len(feature_idx) == 0:
            return 0, 0, 0
        X_train, X_test = X_train[:, feature_idx], X_test[:,feature_idx]
        classifier = LogisticRegression(penalty = 'none').fit(X_train, y_train)
        y_pred_proba = classifier.predict_proba(X_test)
        auc_roc = roc_auc_score(y_test, y_pred_proba[:, 1])
        precision, recall, _ = precision_recall_curve(y_test, y_pred_proba[:, 1])
        auc_prc = metrics.auc(recall, precision)
        acc = accuracy_score(y_test, classifier.predict(X_test))
        return auc_roc, auc_prc, acc

def voting(curr_fs):
    combs = curr_fs.keys()
    combs = sorted(combs)
    voting_matrix = {}
    for comb in combs:
        voting_matrix[comb] = 0
    
    _competitors = {}
    for comb in combs:
        _competitors[comb] = curr_fs[comb]['AUC_mean']
    winners = [comb for m in [max(_competitors.values())] for comb, val in _competitors.items() if val == m]
    for winner in winners:
        voting_matrix[winner] += 1
    _competitors = {}
    for comb in combs:
        _competitors[comb] = curr_fs[comb]['AUPR_mean']
    winners = [comb for m in [max(_competitors.values())] for comb, val in _competitors.items() if val == m]
    for winner in winners:
        voting_matrix[winner] += 1
    _competitors = {}
    for comb in combs:
        _competitors[comb] = curr_fs[comb]['ACC_mean']
    winners = [comb for m in [max(_competitors.values())] for comb, val in _competitors.items() if val == m]
    for winner in winners:
        voting_matrix[winner] += 1
        
    _competitors = {}
    for comb in combs:
        if curr_fs[comb]['AUC_mean'] == 0.5 or curr_fs[comb]['AUC_mean'] == 0:
            _competitors[comb] = 99999
        else:
            _competitors[comb] = curr_fs[comb]['AUC_std']
    winners = [comb for m in [min(_competitors.values())] for comb, val in _competitors.items() if val == m]
    for winner in winners:
        voting_matrix[winner] += 1
    _competitors = {}
    for comb in combs:
        if curr_fs[comb]['AUPR_mean'] == 0.5 or curr_fs[comb]['AUPR_mean'] == 0:
            _competitors[comb] = 99999
        else:
            _competitors[comb] = curr_fs[comb]['AUPR_std']
    winners = [comb for m in [min(_competitors.values())] for comb, val in _competitors.items() if val == m]
    for winner in winners:
        voting_matrix[winner] += 1
    _competitors = {}
    for comb in combs:
        if curr_fs[comb]['ACC_mean'] == 0.5 or curr_fs[comb]['ACC_mean'] == 0:
            _competitors[comb] = 99999
        else:
            _competitors[comb] = curr_fs[comb]['ACC_std']
    winners = [comb for m in [min(_competitors.values())] for comb, val in _competitors.items() if val == m]
    for winner in winners:
        voting_matrix[winner] += 1
    return voting_matrix

FS = {}
params = {'LASSO'   : {'C': 1},
          'SVM-LASSO' : {'C': 1},
          'EN'      : {'alpha': 1},
          'RF'      : {'n_estimators': 100, 'max_depth': 6},
          'XGBoost' : {'learning_rate': 0.05, 'n_estimators': 300, 'max_depth': 6,
                       'reg_alpha':1, 'reg_lambda': 1},
          'MI'      : {'n_neighbors' : 3},
          'reliefF' : {'k' : 3}}
models_FS = {'LASSO'   : FeatureSelector(name='LASSO', params=params['LASSO']),
             'SVM-LASSO' : FeatureSelector(name='SVM-LASSO', params=params['SVM-LASSO']),
             'EN'      : FeatureSelector(name='EN', params=params['EN']),
             'RF'      : FeatureSelector(name='RF', params=params['RF']),
             'XGBoost' : FeatureSelector(name='XGBoost', params=params['XGBoost']),
             'MI'      : FeatureSelector(name='MI', params=params['MI']),
             'reliefF' : FeatureSelector(name='reliefF', params=params['reliefF'])}
tuned_parameters = {'LASSO'   : {'C': [1e-15, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1, 10, 50, 100, 500, 1000, 5000, 1000]},
                    'SVM-LASSO' : {'C': [1e-15, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1, 10, 50, 100, 500, 1000, 5000, 1000]},
                    'EN'      : {'alpha': [1e-15, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1, 10, 100, 1000], 'l1_ratio' : [0.25,0.5,0.75,1]},
                    'RF'      : {'n_estimators': [50,100,200,300,400], 'max_depth' : [2,4,6]},
                    'XGBoost' : {'learning_rate': [0.05, 0.001], 'n_estimators': [50,100,300], 'max_depth': [2,6],
                                 'reg_alpha': [0.1,1,10], 'reg_lambda': [0.1,1,10]},
                    'MI'        : {'n_neighbors': [1,2,3,5,7]},
                    'reliefF'   : {'k':[1,2,3,5,7]}}

dataSetsPath = '/home/yincheng23/ADlasso_manuscript/Datasets/Studies'
outPath = '/home/yincheng23/ADlasso_manuscript/Datasets/data/ML_parm_info/'
refPath = '/home/yincheng23/ADlasso_manuscript/Datasets/data/selected_results/'
os.chdir(dataSetsPath)
dataSets = sorted(os.listdir(dataSetsPath))
for dataSet in dataSets:
    #------------------ chick is working or not------------------
    refile = refPath + dataSet + '/ADlasso.txt'
    if os.path.isfile(refile):
        refile = pd.read_csv(refile,sep = "\t", header = None)
    else:
        print(dataSet, 'did not have reference')
        continue
    tagfolder = outPath + dataSet
    if os.path.exists(tagfolder) == False :
        os.mkdir(tagfolder)
    else:
        continue
    #-------------------------------------------------------------
    os.chdir(dataSetsPath)
    os.chdir(dataSet)
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
    y = [0 if yi == Label[0] else 1 for yi in Label]
    y = np.array(y); y = y.reshape(-1)

    max_num_feat = refile.shape[0]
    k_fold = 5
    kf = StratifiedKFold(n_splits=k_fold, shuffle=True)
    for fs_name, fs_model in models_FS.items():
        print(fs_name.__str__())
        FS.update({fs_name: {}})
        comb = []
        params_name = []
        for name, tun_par in tuned_parameters[fs_name].items():
            comb.append(tun_par)
            params_name.append(name)
        combs = create_grid(comb)
        for comb in combs:
            FS[fs_name].update({comb: {}})
            CV = np.ones([k_fold])
            examineAUC = []; examinePRC = []; examineACC = []
            pev_distribution = []
            fs_model.setParams(comb,params_name,params[fs_name])
            for train_ix, test_ix in kf.split(X, y):
                X_train, X_test = X[train_ix,:], X[train_ix,:]
                y_train, y_test = y[train_ix], y[train_ix]
                fs_model.max_num_feat = max_num_feat
                fs_model.fit(X_train, y_train)
                roc, prc, acc = fs_model.evaluate(X_train, X_test, y_train, y_test)
                examineAUC.append(roc); examinePRC.append(prc); examineACC.append(acc)
                pev_distribution.append(fs_model.get_pvl(prevalence))
            FS[fs_name][comb]['AUC_mean'] = np.mean(examineAUC)
            FS[fs_name][comb]['AUC_std'] = np.std(examineAUC)
            FS[fs_name][comb]['AUPR_mean'] = np.mean(examinePRC)
            FS[fs_name][comb]['AUPR_std'] = np.std(examinePRC)
            FS[fs_name][comb]['ACC_mean'] = np.mean(examineACC)
            FS[fs_name][comb]['ACC_std'] = np.std(examineACC)
            FS[fs_name][comb]['prevalence_mean'] = np.mean(pev_distribution)
            FS[fs_name][comb]['prevalence_std'] = np.std(pev_distribution)
    
    foldpath = tagfolder + '/parm_info.txt'
    with open(foldpath, "a") as f:
        for curr_fs_name,curr_FS_res in FS.items():
            voted_res = voting(curr_FS_res)
            combs = curr_FS_res.keys()
            combs = sorted(combs)
            for comb in combs:
                print('Parameters set: '+ comb.__str__() +' got votes: ' + voted_res[comb].__str__(), file=f)
            performance_winners = [comb for m in [max(voted_res.values())] for comb, val in voted_res.items() if val == m]
            pvl_dict = {win:curr_FS_res[win]['prevalence_mean'] for win in performance_winners}
            prevalence_winners = [comb for m in [max(pvl_dict.values())] for comb, val in pvl_dict.items() if val == m]
            prevalence_winners = prevalence_winners[0]
            print('Best parameters set for: ' + curr_fs_name + ' is: ' +  str(prevalence_winners) + ' with prevalence : ' + 
                str(round(FS[curr_fs_name][prevalence_winners]['prevalence_mean'],3)), '+-',
                str(round(FS[curr_fs_name][prevalence_winners]['prevalence_std'],3)), file=f)
            print('-------------', file=f)

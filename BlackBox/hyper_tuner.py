# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 00:33:42 2021

@author: vardg
"""
import pandas  as pd
import numpy   as np
import sklearn.ensemble as sklens 
import sklearn.model_selection as sklms

excel_CD28_low = pd.read_excel('ERK_times_CD28_low.xlsx', header=0)
CD28_raw = np.array(excel_CD28_low)
param_input = CD28_raw[:,0:-1]
time_output = CD28_raw[:,-1]
param_set ={'learning_rate':    [0.005, 0.01,   0.05,  0.1,  0.5,     1],  \
            'max_features':     [    5,   10,     15,   30,   40,    48],  \
            'subsample':        [  0.2,  0.3,    0.4,  0.5,  0.6,   0.7],  \
            'min_samples_leaf': [    1,    5,     10,   50,  100,   500] }

gradtree = sklens.GradientBoostingRegressor(n_estimators=5000, random_state=19951212)
gradtree_tuned = sklms.GridSearchCV(estimator=gradtree, param_grid=param_set, verbose=1, cv=5, n_jobs=23, pre_dispatch=50)
print("Best Score: ", gradtree_tuned.best_score_, "Best Params: ", gradtree_tuned.best_params_)

# Best Score was found to be XXX, for parameter set XXX
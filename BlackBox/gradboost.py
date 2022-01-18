# -*- coding: utf-8 -*-
"""
Created on Wed Mar  3 21:43:20 2021

@author: vardg
"""
import pandas  as pd
import numpy   as np
import sklearn.ensemble as sklens 
import sklearn.model_selection as sklms
import sklearn.metrics as sklmet
import sklearn.inspection as skli

def gradient_boosted_analysis(filename):
    excel_dframe = pd.read_excel(filename, header=0, engine='openpyxl')
    data_raw = np.array(excel_dframe)
    param_input = data_raw[:,0:-1]
    time_output = data_raw[:,-1]
    
    x_train = param_input
    y_train = time_output
    gradtree_regress = sklens.GradientBoostingRegressor(n_estimators=7500,  max_features="log2", subsample=0.5, n_jobs=24, random_state=19951212)
    scores = sklms.cross_validate(gradtree_regress, x_train, y_train, scoring=['r2', 'neg_mean_absolute_error'])

    print("Source File: ",  filename)
    print("Rsq  Mean: ",  np.mean(scores['test_r2']), "Rsq SEM: ", np.sqrt(1/5)*np.std(scores['test_r2']))
    print("MAPE Mean: ",  np.mean(scores['test_neg_mean_absolute_error']), "MAPE SEM: ", np.sqrt(1/5)*np.std(scores['test_neg_mean_absolute_error']))
    
    param_names = ['ant_Kd', 'ktrans', 'KmA2', 'KmB1', 'KmB2', 'KmC1', 'Xi', 'Kcat_LCKPU_CD3z', 'cat_UU_394', 'cat_PU_505', \
                   'CSKon', 'CSKoff_PU', 'CSKoff_UU', 'cat_CSK_PU', 'on', 'off_PU_PP', 'ZAPu_on', 'GADS_SLP_on', 'Kcat_LCKPU_ZAP315', \
                   'at_ZAP', 'LAT_Grb2_on', 'on_PLC_LAT', 'Kcat_PLCg', 'on_RasGRP', 'Kcat_RasGRP', 'Kcat1', 'Km1', 'Vmax2', 'Km2', \
                   'Kcat3', 'Kcat4', 'Km4', 'Vmax5', 'Km5', 'Vmax6', 'Km6', 'Ki1', 'Kcat7', 'Kcat8', 'Vmax9', 'Ki2', 'km', \
                   'Kcat_CD45_LCK505', 'Kcat_CD45_A1', 'SHP1_on', 'kcat_LCKpu_SHP1', 'kcat_SHP1', 'kcat_SHP1_SHP1']
    
    per_imp = skli.permutation_importance(gradtree_regress, x_train, y_train, n_repeats=5, n_job=24, random_state=19951212)
    print("Permutation Importance Means: ", per_imp.importances_mean)
    print("Permutation Importance SEMs:  ", np.sqrt(1/5)*per_imp.importances_std)
#    mplt.bar(np.arange(48), per_imp.importances_mean, yerr=per_imp.importances_std)
#    mplt.set_title('Permutation Importances for Parameters, data set ' + filename)
#    mplt.set_xticklabels(param_names)
#    mplt.xticks(rotation=75) 
    
    print('____________________________________________________')
    


gradient_boosted_analysis('ERK_times_CD28_low.xlsx')
#gradient_boosted_analysis('ERK_times_CD3z_low.xlsx')
#gradient_boosted_analysis('ERK_times_CD28_high.xlsx')
#gradient_boosted_analysis('ERK_times_CD3z_high.xlsx')

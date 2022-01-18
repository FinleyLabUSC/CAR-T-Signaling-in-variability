# -*- coding: utf-8 -*-
"""
Created on Sat Mar  6 17:12:51 2021

@author: vardg
"""

import numpy as np
import scipy.stats as scistat
import matplotlib.pyplot as mplt


def perimp_analysis(filename, per_means, per_sems):
    
    t_scores = per_means/(per_sems)
    p_val = scistat.t.sf(t_scores, df=4)
    
    param_names = ['ant_Kd', 'ktrans', 'KmA2', 'KmB1', 'KmB2', 'KmC1', 'Xi', 'Kcat_LCKPU_CD3z', 'cat_UU_394', 'cat_PU_505', \
               'CSKon', 'CSKoff_PU', 'CSKoff_UU', 'cat_CSK_PU', 'on', 'off_PU_PP', 'ZAPu_on', 'GADS_SLP_on', 'Kcat_LCKPU_ZAP315', \
               'kcat_ZAP', 'LAT_Grb2_on', 'on_PLC_LAT', 'Kcat_PLCg', 'on_RasGRP', 'Kcat_RasGRP', 'Kcat1', 'Km1', 'Vmax2', 'Km2', \
               'Kcat3', 'Kcat4', 'Km4', 'Vmax5', 'Km5', 'Vmax6', 'Km6', 'Ki1', 'Kcat7', 'Kcat8', 'Vmax9', 'Ki2', 'km', \
               'Kcat_CD45_LCK505', 'Kcat_CD45_A1', 'SHP1_on', 'kcat_LCKpu_SHP1', 'kcat_SHP1', 'kcat_SHP1_SHP1']
    
    mplt.figure()
    mplt.bar(np.arange(1,49), per_means, yerr=per_sems, capsize=5)
    #mplt.title('Permutation Importances for Parameters, Data Set: '+filename, fontdict={'fontsize':25})
    mplt.xticks(ticks=np.arange(1,49), labels=param_names, fontsize=35, rotation=90, fontweight='bold')
    
    mplt.ylim(0,0.7)
    mplt.yticks(ticks=np.arange(0,0.7,0.1), fontsize=35, fontweight='bold') 
    
    mplt.xlabel("Parameter Name", fontsize=75, fontweight='bold')
    mplt.ylabel("Importance Score", fontsize=75, fontweight='bold')
    
    #mplt.plot(np.arange(1,49), (per_means+0.02)*(p_val<0.01)+0.01, marker='*', ms=17, color='k', linestyle="None")
    #mplt.plot(np.arange(1,49), (per_means+0.03)*(p_val<0.01)+0.01, marker='*', ms=17, color='k', linestyle="None")
    return p_val
  
    
def perimp_analysis_o(filename, per_means, per_sems):
    
    t_scores = per_means/(per_sems)
    p_val = scistat.t.sf(t_scores, df=4)
    
    param_names = ['ant_Kd', 'ktrans', 'KmA2', 'KmB1', 'KmB2', 'KmC1', 'Xi', 'Kcat_LCKPU_CD3z', 'cat_UU_394', 'cat_PU_505', \
               'CSKon', 'CSKoff_PU', 'CSKoff_UU', 'cat_CSK_PU', 'on', 'off_PU_PP', 'ZAPu_on', 'GADS_SLP_on', 'Kcat_LCKPU_ZAP315', \
               'kcat_ZAP', 'LAT_Grb2_on', 'on_PLC_LAT', 'Kcat_PLCg', 'on_RasGRP', 'Kcat_RasGRP', 'Kcat1', 'Km1', 'Vmax2', 'Km2', \
               'Kcat3', 'Kcat4', 'Km4', 'Vmax5', 'Km5', 'Vmax6', 'Km6', 'Ki1', 'Kcat7', 'Kcat8', 'Vmax9', 'Ki2', 'km', \
               'Kcat_CD45_LCK505', 'Kcat_CD45_A1', 'SHP1_on', 'kcat_LCKpu_SHP1', 'kcat_SHP1', 'kcat_SHP1_SHP1']
    
    mplt.figure()
    mplt.bar(np.arange(1,49), per_means, yerr=per_sems, capsize=5, color=[1,0.4,0,1])
    #mplt.title('Permutation Importances for Parameters, Data Set: '+filename, fontdict={'fontsize':25})
    mplt.xticks(ticks=np.arange(1,49), labels=param_names, fontsize=35, rotation=90, fontweight='bold')
    
    
    mplt.ylim(0,0.7)
    mplt.yticks(ticks=np.arange(0,0.7,0.1), fontsize=35, fontweight='bold') 
    
    mplt.xlabel("Parameter Name", fontsize=75, fontweight='bold')
    mplt.ylabel("Importance Score", fontsize=75, fontweight='bold')
    
    #mplt.plot(np.arange(1,49), (per_means+0.02)*(p_val<0.01)-0.001, marker='*', ms=17, color='k', linestyle="None")
    #mplt.plot(np.arange(1,49), (per_means+0.03)*(p_val<0.01)-0.001, marker='*', ms=17, color='k', linestyle="None")
    return p_val
###############################################################################

filename = 'ERK_times_CD28_low.xlsx'
per_means = np.array([ 1.90527528e-01,  1.18208803e-02,  3.35660948e-03,  4.90536541e-03, \
                       1.19282822e-02,  7.51768270e-03,  2.39021713e-03,  5.95192953e-01, \
                       2.41342752e-02,  4.12524062e-02,  6.90514128e-02,  5.71039377e-02, \
                       3.90876977e-04,  4.31438880e-02,  8.43017065e-03,  1.47107432e-03, \
                       1.29628240e-01,  5.36060142e-03,  9.87273388e-02,  3.46166219e-01, \
                       1.03499567e-04,  3.65075949e-04,  1.02582179e-01, -7.43320517e-05, \
                       1.03942423e-01,  9.74912223e-02,  2.52705107e-03,  5.90703896e-02, \
                       2.32390248e-02,  7.46248818e-03,  1.45497600e-02,  6.20984863e-03, \
                       1.35416717e-02,  8.25593992e-03,  6.11128962e-03,  5.44230417e-03, \
                       1.70972007e-02,  5.91862855e-03,  3.54272002e-03,  2.72763299e-03, \
                       3.62873029e-04,  5.35787442e-02,  1.73525724e-01,  8.45720223e-02, \
                       3.78829265e-03,  2.50768625e-04,  4.66679503e-03,  7.95457963e-06])

per_sems = np.array([1.85416390e-03, 1.92416081e-04, 1.09523743e-04, 1.55411004e-04, \
                     2.53712833e-04, 8.12694146e-05, 4.81414774e-05, 3.84005021e-03, \
                     1.66473437e-04, 3.69106813e-04, 7.89429849e-04, 3.81332341e-04, \
                     3.49498500e-05, 4.82906169e-04, 2.26263060e-04, 6.10680614e-05, \
                     6.33763407e-04, 1.85031716e-04, 6.64536234e-04, 1.44556961e-03, \
                     8.61292218e-05, 3.92099765e-05, 3.92256460e-04, 4.32401573e-05, \
                     6.35598533e-04, 4.97375604e-04, 7.47383230e-05, 1.04742366e-04, \
                     3.26250668e-04, 1.37609601e-04, 2.29825067e-04, 1.69392828e-04, \
                     1.89333629e-04, 1.24812756e-04, 1.60724620e-04, 1.53040989e-04, \
                     2.08212930e-04, 1.27911234e-04, 1.04710788e-04, 9.58724602e-05, \
                     3.12602980e-05, 2.67183681e-04, 1.67156293e-03, 4.69072032e-04, \
                     7.17202892e-05, 8.63105513e-05, 6.19880785e-05, 4.28273827e-05])

p_val = perimp_analysis_o(filename, per_means, per_sems)

print("Source File: ERK_times_CD28_low.xlsx")
print("Rsq  Mean:  0.8996310730752842, Rsq SEM:  0.00040363963530977263")
print("EV   Mean:  0.8996389479536469, EV  SEM:  0.0004028415656165327")
print("______________________________________________________________________")

###############################################################################

filename = 'ERK_times_CD3z_low.xlsx'
per_means = np.array([5.43921387e-01, 4.56800913e-02, 4.79615398e-03, 1.14298364e-02, \
                      1.37918541e-02, 2.52025127e-03, 3.15947827e-03, 6.65119003e-01, \
                      2.41494768e-02, 7.10471423e-02, 1.29647495e-01, 8.18686300e-02, \
                      9.87281891e-04, 8.71663834e-02, 2.40104688e-02, 2.39460977e-03, \
                      1.15715435e-01, 9.44497206e-03, 8.35056105e-02, 3.22832623e-01, \
                      5.50712065e-05, 4.24316550e-04, 9.81551472e-02, 1.26713087e-04, \
                      9.92777728e-02, 9.20950583e-02, 2.97430861e-03, 8.18602623e-02, \
                      2.70672832e-02, 9.04336624e-03, 1.75873773e-02, 8.38145371e-03, \
                      1.63792374e-02, 9.49132834e-03, 7.49791072e-03, 6.61448693e-03, \
                      2.27524535e-02, 6.08027245e-03, 4.87685516e-03, 2.87989659e-03, \
                      4.04570356e-04, 6.96147178e-02, 2.36964807e-01, 1.69988995e-01, \
                      7.47549379e-03, 1.01513164e-04, 1.07280788e-02, 6.04291085e-04])

per_sems = np.array([1.79198557e-03, 4.80745481e-04, 9.94697293e-05, 1.78134525e-04, \
                     1.71113072e-04, 8.54071327e-05, 5.55776472e-05, 2.50968340e-03, \
                     2.83519823e-04, 2.90666736e-04, 9.45243895e-04, 3.41159640e-04, \
                     7.12185750e-05, 5.68593397e-04, 3.30006446e-04, 7.32454885e-05, \
                     9.12893707e-04, 1.76728493e-04, 9.31314194e-04, 6.03777948e-04, \
                     3.12972842e-05, 8.16889244e-05, 4.08086741e-04, 3.23531243e-05, \
                     3.27257487e-04, 8.59715219e-04, 9.90764040e-05, 3.14170512e-04, \
                     2.10279747e-04, 1.94847800e-04, 1.31882877e-04, 1.02050000e-04, \
                     1.49434328e-04, 1.28672103e-04, 2.06014655e-04, 1.27841174e-04, \
                     1.06115189e-04, 4.59458712e-05, 1.42950082e-04, 7.86396415e-05, \
                     2.19556236e-05, 3.67915075e-04, 7.88304246e-04, 1.32074810e-03, \
                     7.57582413e-05, 6.46042932e-05, 1.61699554e-04, 3.71951541e-05])

p_val = perimp_analysis(filename, per_means, per_sems)
print("Source File: ERK_times_CD3z_low.xlsx")
print("Rsq  Mean:  0.9047267669058787, Rsq SEM:  0.0008211004013158351")
print("EV   Mean:  0.9047337157001543, EV  SEM:  0.0008225519441557797")
print("______________________________________________________________________")

###############################################################################

filename = 'ERK_times_CD28_high.xlsx'
per_means = np.array([ 3.49684174e-02,  4.78550447e-02,  4.84849534e-03,  1.11764711e-02, \
                       2.35956452e-02,  8.53781697e-03,  1.78969468e-02,  2.80762061e-01, \
                       6.16966106e-02,  1.66624141e-02,  1.56069470e-02,  1.72390295e-02, \
                       9.37156671e-04,  9.71969490e-03,  4.74645737e-03,  1.90681318e-03, \
                       1.18704696e-01,  2.50511784e-03,  2.11917427e-01,  2.86753169e-01, \
                       1.42271904e-03,  5.27648902e-04,  1.35949746e-01, -2.96852831e-04, \
                       1.12729068e-01,  1.21640922e-01,  2.48345330e-03,  3.39308216e-02, \
                       1.44128170e-02,  1.12894636e-02,  1.81818841e-02,  6.24526882e-03, \
                       1.06563468e-02,  7.15074419e-03,  5.87941827e-03,  6.54116662e-03, \
                       1.00430003e-02,  9.06602829e-03,  3.27767026e-03,  1.83885378e-03, \
                      -2.72606066e-04,  5.28205505e-02,  9.31068948e-03,  1.68170202e-02, \
                       3.15601235e-02,  1.14796463e-06,  4.35378495e-02,  1.71931460e-03])

per_sems = np.array([0.00054132, 0.00144089, 0.00043655, 0.00032674, 0.0008891 , 0.00022062, \
                     0.00029354, 0.00178685, 0.00129409, 0.00018291, 0.00039157, 0.00052727, \
                     0.00018657, 0.00030246, 0.00055032, 0.00019516, 0.00209835, 0.00024125, \
                     0.00463505, 0.00139834, 0.00015441, 0.00019015, 0.00127446, 0.0001492 , \
                     0.00094109, 0.00087522, 0.0001138 , 0.00072516, 0.0002981 , 0.00045707, \
                     0.00063264, 0.00015976, 0.00017576, 0.00019391, 0.00033846, 0.00027743, \
                     0.00047181, 0.00057033, 0.00014624, 0.00026869, 0.00024326, 0.00107567, \
                     0.00030093, 0.00069276, 0.00035555, 0.0002349 , 0.00070474, 0.00019427])


p_val = perimp_analysis_o(filename, per_means, per_sems)

print("Source File: ERK_times_CD28_high.xlsx")
print("Rsq  Mean:  0.8072100271894289, Rsq SEM:  0.003520363803811731")
print("EV   Mean:  0.8072467886307004, EV  SEM:  0.003535817899467444")
print("______________________________________________________________________")

###############################################################################

filename = 'ERK_times_CD3z_high.xlsx'
per_means = np.array([6.51832578e-02,  1.25297281e-01,  5.45188021e-03,  9.19604030e-03, \
                      1.97086250e-02,  4.81252266e-03,  2.38940665e-02,  3.06360026e-01, \
                      4.32408852e-02,  1.79799963e-02,  2.57305405e-02,  3.04953757e-02, \
                     -1.10824245e-04,  2.00650485e-02,  1.07981349e-02,  1.52252765e-03, \
                      9.26557614e-02,  1.56433275e-03,  2.06926490e-01,  2.84864443e-01, \
                      2.83868750e-04,  3.94733748e-04,  1.02157100e-01, -2.90529638e-04, \
                      9.43199052e-02,  9.65919100e-02,  1.50154333e-03,  3.21664085e-02, \
                      1.82899558e-02,  7.83103154e-03,  1.31345134e-02,  3.72399577e-03, \
                      8.40746270e-03,  5.93278369e-03,  5.72567833e-03,  6.85747505e-03, \
                      1.19568336e-02,  6.97265604e-03,  2.37990528e-03,  2.15541238e-03, \
                     -4.92859713e-05,  3.42320747e-02,  1.08374794e-02,  3.01827134e-02, \
                      5.93633566e-02,  2.53331660e-03,  9.66606286e-02,  2.13212664e-03])

per_sems = np.array([6.93058212e-04, 6.26443994e-04, 1.55555496e-04, 3.95442766e-04, \
                     3.47334708e-04, 1.85298131e-04, 6.09101065e-04, 3.10920997e-03, \
                     9.22107030e-04, 5.09200422e-04, 6.47273480e-04, 5.61010057e-04, \
                     1.78460068e-04, 3.12241859e-04, 1.94711951e-04, 9.86323852e-05, \
                     1.15284629e-03, 1.91179824e-04, 1.50692870e-03, 2.50246974e-03, \
                     8.06928052e-05, 1.43633783e-04, 5.17188993e-04, 9.31263080e-05, \
                     9.70974774e-04, 7.86581786e-04, 1.80159447e-04, 5.11514043e-04, \
                     4.52855622e-04, 2.23120328e-04, 2.83454431e-04, 1.23098910e-04, \
                     1.99235176e-04, 3.35280821e-04, 2.85991798e-04, 2.82856474e-04, \
                     1.73700473e-04, 1.24968058e-04, 1.84178762e-04, 7.84806903e-05, \
                     1.42670204e-04, 2.78545840e-04, 2.32408148e-04, 6.20595125e-04, \
                     7.57439482e-04, 1.06849907e-04, 7.20967966e-04, 2.73341003e-04])

p_val = perimp_analysis(filename, per_means, per_sems)

print("Source File: ERK_times_CD3z_high.xlsx")
print("Rsq  Mean:  0.8157799277034739, Rsq SEM:  0.0023970894901433333")
print("EV   Mean:  0.8158011376630616, EV  SEM:  0.0024110760177379962")
print("______________________________________________________________________")
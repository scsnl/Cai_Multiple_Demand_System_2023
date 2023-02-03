##########################
#### Weidong Cai, 2021

import os
import sys
import numpy as np
import pandas as pd
import scipy.stats as spss
import sklearn
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn import preprocessing
from sklearn.cross_decomposition import CCA

### read BSDS state OR data
state_OR_fname = ''
state_OR_df = pd.read_csv(state_OR_fname)

### read behav data
behav_fname = ''
behav_df = pd.read_csv(behav_fname)

# select variables for CCA
X_varnames = [] # example: ['SH_OR', 'S2_OR', 'S3_OR', 'S4_OR']
y_varnames = [] # example: ['Con_ACC', 'Incon_ACC', 'Con_RT', 'Incon_RT']

X = state_OR_df[X_varnames]
y = behav_df[y_varnames]

X_scaled =  preprocessing.scale(X)
y_scaled =  preprocessing.scale(y)

if (y.shape[1] < X.shape[1]):
  cca = CCA(n_components=y.shape[1])
else:
  cca = CCA(n_components=X.shape[1])

cca.fit(X_scaled, y_scaled)

key_component_id = 0

# return canoncial coefficient
X_c, y_c = cca.transform(X_scaled, y_scaled)
X_c = [i[key_component_id] for i in X_c]
y_c = [i[key_component_id] for i in y_c]

# print canonical correlation for the 1st component
print(spss.pearsonr(X_c,y_c))
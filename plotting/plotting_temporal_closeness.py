import os
import sys
import numpy as np
import pandas as pd
import scipy.io as spio
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

data_fname = '' # BSDS model statematch results path
ofname = '' # figure output file name

data = spio.loadmat(data_fname)

corr_mtx = data['matched_state_cor_mtx']
train_ids = data['train_state_ids']
test_ids = data['matched_test_state_ids']

train_state_labels = ['S' + str(ii) for ii in train_ids[0]]
test_state_labels = ['S' + str(ii) for ii in test_ids[0]]

sns.set(style="white")
sns.set(font_scale=1.1)
annot_size = 10

# use matlab default colormap
cmap_list = []
cmap_fname = 'matlab_default_cmap.txt'
cmap_data = np.loadtxt(cmap_fname)
for i in range(len(cmap_data)):
  icmap = cmap_data[i]
  icmap_new = np.append(icmap, 1.0)
  cmap_list.append(icmap_new)

cmap_tuple = tuple(tuple(i) for i in cmap_list)
myColors = cmap_tuple

cmap = LinearSegmentedColormap.from_list('Custom', myColors, 20)

f, ax = plt.subplots(figsize=(5, 4))

# Draw the heatmap with the mask and correct aspect ratio
sns.heatmap(corr_mtx, annot=True, yticklabels='', xticklabels='', cmap=cmap, center=0, square=True, linewidths=.5, vmax=0.8, vmin=-0.8)
ax.xaxis.set_ticks_position('top')

plt.savefig(ofname, dpi=1000)
plt.show()
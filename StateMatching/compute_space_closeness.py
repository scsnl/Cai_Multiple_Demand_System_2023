import numpy as np
import h5py
import scipy.io as spio
import sys

from symmetrized_kld import SymmetrizedKullbackLeiblerDivergenceGaussians

def compute_space_closeness(train_model_fname, test_model_fname):
  '''
  This script is used to compute latent space distance between states from BSDS models
  
  Inputs:
  train_model_fname: BSDS model from train data
  test_model_fname: BSDS model from test data

  Output:
  closeness_mtx : an n*m matrix for space closeness
  n is the number of states in the train model, m is the number of states
  in the test model
  
  wdcai@stanford.edu (2022)
  '''
  print('test')
  model_train = h5py.File(train_model_fname, 'r')
  model_test = h5py.File(test_model_fname, 'r')
  
  modelinTrain = model_train['model']
  modelinTest = model_test['model']
  
  state_idx_train = np.where(np.array(modelinTrain['fractional_occupancy_group_wise'][:]) > 0)[0]
  state_idx_test = np.where(np.array(modelinTest['fractional_occupancy_group_wise'][:]) > 0)[0]
  
  skld_mtx = np.zeros((len(state_idx_train), len(state_idx_test)))
  
  for i in range(len(state_idx_train)):
    istate_idx_train = state_idx_train[i]
    istate_mean_train = model_train[modelinTrain['estimated_mean'][istate_idx_train,0]].value[0]
    istate_cov_train = model_train[modelinTrain['estimated_covariance'][istate_idx_train,0]].value
    for j in range(len(state_idx_test)):
      jstate_idx_test = state_idx_test[j]
      jstate_mean_test = model_test[modelinTest['estimated_mean'][jstate_idx_test,0]].value[0]
      jstate_cov_test = model_test[modelinTest['estimated_covariance'][jstate_idx_test,0]].value
      skld = SymmetrizedKullbackLeiblerDivergenceGaussians(cov1=istate_cov_train, cov2=jstate_cov_test, mean1=istate_mean_train, mean2=jstate_mean_test)
      skld_mtx[i, j] = skld.evaluate()

  closeness_mtx = 1./skld_mtx

  return closeness_mtx

closeness_mtx = compute_space_closeness('model1.mat', 'model2.mat')
print(closeness_mtx)
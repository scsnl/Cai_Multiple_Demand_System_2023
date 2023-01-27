function corr_mtx_r = compute_temporal_closeness(train_model_fname, test_model_fname, test_data_fname)

%% -----------------------------------------------------------------
% This script is used for matching states from different BSDS models.
% A trained model was applied on the test data to compute posterior probability of each state
% from the trained model in the test data, from which temporal correlation between the
% posterior probability of each state from the train model in the test data and the posterior
% probability of each state from test model in the test data (saved in the test model).
%
% Inputs:
% train_model_fname: BSDS model from train data
% test_model_fname: BSDS model from test data
% test_data_fname: test data
%
% Output:
% corr_mtx_r: an n*m matrix for temporal closeness (correlation)
% n is the number of states in the train model, m is the number of states
% in the test model
%
% wdcai@stanford.edu (2022)
%% -----------------------------------------------------------------------

% load trained model
model_train = load(train_model_fname);

% load test data
data_test = load(test_data_fname);

% load test model
model_test = load(test_model_fname);

% compute posterior probability of each state from the train model in test data
predQns = computeQnsFromGivenNetForNewData(data_test.data, model_train.model.net);

% initialize temporal correlation matrix
corr_mtx_r = zeros(size(predQns,2), size(model_test.model.net.hidden.Qns,2));
corr_mtx_p = corr_mtx_r;

for i = 1:size(predQns,2)
    for j = 1:size(model_test.model.net.hidden.Qns,2)
        [corr_mtx_r(i,j), corr_mtx_p(i,j)] = corr(predQns(:,i),model_test.model.net.hidden.Qns(:,j));
    end
end

import os
import numpy as np
import pandas as pd
import lightgbm as lgb
from sklearn.metrics import mean_squared_error

#load affinity data
df_affinity_train = pd.read_csv('F:/ComData/df_affinity_train.csv')
df_affinity_toBePredicted = pd.read_csv('F:/ComData/df_affinity_test_toBePredicted.csv')
df_affinity_toBePredicted['Ki'] = 6
data = pd.concat([df_affinity_train,df_affinity_toBePredicted])

os.chdir('F:/output/')
# load molecule data
df_molecule = pd.read_csv('molecule.csv')
data = data.merge(df_molecule, on='Molecule_ID', how='left')
#load ai coefficient
feature = pd.read_csv('aliphatic_index_result.csv')
data = data.merge(feature, on='Protein_ID', how='left')
del feature
#load GRAVY
feature = pd.read_csv('GRAVY_result.csv')
data = data.merge(feature, on='Protein_ID', how='left')
del feature
#load word2vec
feature = pd.read_csv('w2v_result.csv')
data = data.merge(feature, on='Protein_ID', how='left')
del feature
#load Ext
feature = pd.read_csv('Ext_result.csv')
data = data.merge(feature, on='Protein_ID', how='left')
del feature
#load protein size & pI
feature = pd.read_csv('Protein_size_pi.csv')
data = data.merge(feature, on='Protein_ID', how='left')
del feature

#LGB setting
train_feature = data[data['Ki'] > 6]
test_feature = data[data['Ki'] <= 6]

offline_choice_list = [1, 2]
idx = train_feature['Protein_ID'].isin(offline_choice_list)
offline_test = train_feature[idx]
train_xy = train_feature[(~idx)]
label_train = train_xy['Ki']
label_test = offline_test['Ki']

submission = test_feature[['Protein_ID', 'Molecule_ID']]
train_xy = train_xy.drop('Ki', axis=1)
train_xy = train_xy.drop('Protein_ID', axis=1)
train_xy = train_xy.drop('Molecule_ID', axis=1)

test_feature = test_feature.drop('Ki', axis=1)
test_feature = test_feature.drop('Protein_ID', axis=1)
test_feature = test_feature.drop('Molecule_ID', axis=1)

offline_test = offline_test.drop('Ki', axis=1)
offline_test = offline_test.drop('Protein_ID', axis=1)
offline_test = offline_test.drop('Molecule_ID', axis=1)

train = lgb.Dataset(train_xy, label=label_train)
test = lgb.Dataset(offline_test, label=label_test, reference=train)

#parameters
params = {
    'task': 'predict',
    'boosting_type': 'gbdt',
    'objective': 'regression_l2',
    'metric': ['l2'],
    'min_child_weight': 3,
    'max_depth': 4,
    'num_leaves': 78,
    'is_unbalance': True,
    'lambda_l2': 5,
    'learning_rate': 0.045,
    'seed': 3000,
    'nthread': 4,
}

# Train
num_round = 5000
gbm = lgb.train(params,
                train,
                num_round,
                verbose_eval=50,
                valid_sets=[train, test]
                )

# Predicts
preds_offline = gbm.predict(offline_test)
print('MSE:', mean_squared_error(label_test, preds_offline))
preds_sub = gbm.predict(test_feature)
submission['Ki'] = preds_sub
submission.to_csv('F:/output/lgb.csv', index=False)

import pandas as pd
import numpy as np
from os.path import join
from random import shuffle

data_dir = '/Users/daiwei89/storage/cs_proj/dap/data'
os_core = 'GBM_OS_core.txt'
num_folds = 100

# Development set.
dev_file = 'GBM_dev_sample_list.txt'
non_dev_file = 'GBM_non_dev_sample_list.txt'

# 100-fold split training set within dev
dev_train_file = 'GBM_dev_train_sample_list.txt'

# 100-fold split test set within dev
dev_test_file = 'GBM_dev_test_sample_list.txt'

df = pd.read_csv(join(data_dir, os_core), sep='\t')
sample_mat = df.feature.as_matrix()
num_samples = sample_mat.shape[0]
num_samples_dev = int(num_samples * 0.8)
print 'num_dev:', num_samples_dev, 'num_test:', num_samples - num_samples_dev
np.random.shuffle(sample_mat)
dev_mat = sample_mat[:num_samples_dev]
non_dev_mat = sample_mat[num_samples_dev:]

np.savetxt(join(data_dir, dev_file), dev_mat, delimiter='\t', fmt='%s')
np.savetxt(join(data_dir, non_dev_file), non_dev_mat, delimiter='\t', fmt='%s')

multi = np.tile(dev_mat, (num_folds, 1)).transpose()
# shuffle only shuffles along 0-th index (exchanging rows)
for i in range(multi.shape[1]):
  np.random.shuffle(multi[:,i])

num_samples_dev = multi.shape[0]
num_samples_dev_train = int(num_samples_dev * 0.8)
dev_train = multi[:num_samples_dev_train, :]
dev_test = multi[num_samples_dev_train:, :]
print 'num_samples_dev:', num_samples_dev, 'num_samples_dev_train:', \
  num_samples_dev_train, 'num_samples_dev_test:', \
  num_samples_dev - num_samples_dev_train
np.savetxt(join(data_dir, dev_train_file), dev_train, delimiter='\t', fmt='%s')
np.savetxt(join(data_dir, dev_test_file), dev_test, delimiter='\t', fmt='%s')
print 'saved to', join(data_dir, dev_train_file)
print 'saved to', join(data_dir, dev_test_file)

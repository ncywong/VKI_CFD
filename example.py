import os
import numpy as np
import sys
from shutil import copy
import time
from util_funcs import *
from auto_run import *

#%% Generate FFD box and designs, the sampling method is provided in t_gendesign_doe.py
# alternatively, provide a DoE in the root directory 
# os.system('python t_gendesign_doe.py')
# X = np.loadtxt('....')

#%% make the directory structure with deformed geometries according to the designs
data_dir_name = 'random_designs'

makedirs(X, data_dir_name)

#%% run CFD
run_locally(data_dir_name)

#%% For example, get loss values
num_designs = 1000
loss_all = np.zeros(num_designs)

for i in range(num_designs):
    loss_all[i] = get_loss_from_file(data_dir_name)

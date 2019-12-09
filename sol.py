import os
import numpy as np
import sys
from shutil import copy
import time

num_designs = 51

for n in range(num_designs):
        os.chdir('design_%04d' % n)
        os.system('SU2_SOL sol_cfd.cfg')
        os.chdir('..')

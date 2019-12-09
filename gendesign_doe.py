import numpy as np
from pyDOE import lhs
import os
'''
This script generates a latin-hypercube DOE within a bounding box using pyDOE
and also generates the mesh with FFD box
'''
# This is also the number of FFD control points. This should be divisible by 2.
num_vars = 20 

bounds = np.array([0.01 for _ in range(num_vars)])

num_designs = 2000

x_design = lhs(num_vars, num_designs)
x_design = 2.0*x_design - 1.0 # so that they are now btn [-1,1]
x_design *= bounds

np.savetxt('VP_data2_designs.txt', x_design, delimiter=',')

cfg_file = open('gen_ffd_box.cfg','r')
configs = cfg_file.readlines()
cfg_file.close()

line_389 = "FFD_DEGREE= (1, " + str((num_vars//2) - 1) + ", 0)\n"

configs[388] = line_389

cfg_file = open('gen_ffd_box.cfg','w')
cfg_file.writelines(configs)
cfg_file.close()

os.system('SU2_DEF gen_ffd_box.cfg')

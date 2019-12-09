import os
import numpy as np
import sys
from shutil import copy
import time

def makedirs(x, dir_name, bound = 0.01, vis=False):

    title = dir_name

    # Import design params and locations
    all_designs = x * bound

    num_exps = all_designs.shape[0]
    num_vars = all_designs.shape[1]

    os.chdir('/home/nick/Documents/vki/data')
    try:
        os.mkdir(title)
    except:
        pass

    start_time = time.time()
    for n in range(num_exps):

        current_design = all_designs[n,:]

        # Create cfg file for SU2_DEF
        cfg_file = open('../deform_ffd.cfg','r')
        configs = cfg_file.readlines()
        cfg_file.close()

        line_342 = "DV_KIND= "
        line_367 = "DV_PARAM= "
        line_370 = "DV_VALUE= "

        for i in range(num_vars):
            line_342 += "FFD_CONTROL_POINT_2D, "
            line_367 += "(MAIN_BOX," + str(i // (num_vars // 2)) + ", " + str(i % (num_vars // 2)) + ", 0.0, 1.0); "
            line_370 += str(current_design[i]) + ", "

        line_342 = line_342[:-2] + '\n'
        line_367 = line_367[:-2] + '\n'
        line_370 = line_370[:-2] + '\n'

        configs[341] = line_342
        configs[366] = line_367
        configs[369] = line_370
        if vis:
            configs[372] = 'VISUALIZE_SURFACE_DEF= YES\n'
            configs[375] = 'VISUALIZE_VOLUME_DEF= YES\n'

        os.chdir(title)
        copy('../../submit_jobarray', '.')
        copy('../../sol.py', '.')
        dir_name = 'design_%04d' % (n)
        try:
            os.mkdir(dir_name)
        except:
            # Directory already exists?
            pass
        copy('../../mesh_ffd_box.su2',dir_name)
        copy('../../run_cfd.cfg', dir_name)
        copy('../../sol_cfd.cfg', dir_name)
        copy('../../baseline_flow.dat', dir_name)

        os.chdir(dir_name)

        new_cfg_file = open('deform_ffd.cfg','w')
        new_cfg_file.writelines(configs)
        new_cfg_file.close()

        # Run SU2_DEF on this
        os.system('SU2_DEF deform_ffd.cfg > /dev/null')

        # Return to parent dir
        os.chdir('../..')

    print('done def')

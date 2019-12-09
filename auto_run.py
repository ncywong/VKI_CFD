import numpy as np
import os
import sys


def run_locally(dir_name):
    os.chdir('data/' + dir_name)
    num_designs = len([i for i in os.listdir('.') if 'design' in i])
    with open('sol.py', 'r') as f:
        sol_configs = f.readlines()
        sol_configs[6] = 'num_designs = %d\n' % (num_designs)
    with open('sol.py', 'w') as f:
        f.writelines(sol_configs)

    for i in range(num_designs):
        os.chdir('design_%04d' % i)
        os.system('mpirun.mpich -n 4 SU2_CFD run_cfd.cfg > /dev/null')
        os.chdir('..')
    print('done cfd')
    os.system('python sol.py > /dev/null')
    output_title = dir_name + '_fin'
    os.chdir('..')
    try:
        os.system('rm -r ' + output_title)
    except:
        pass
    try:
        os.mkdir(output_title)
    except:
        pass
    os.chdir(dir_name)
    for folder in [i for i in os.listdir('.') if 'design' in i]:
        print(folder)
        os.system('cp -r ' + folder + ' ../' + output_title + '/' + folder)
    print('all done')

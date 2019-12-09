import sys
import numpy as np
import os
from shutil import copy
import pyvista
import pandas as pd
from makedirs import makedirs
import matplotlib
matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
from vtk.numpy_interface import algorithms as algs
from vtk.numpy_interface import dataset_adapter as dsap

def massav(slice, array):

    y = slice.points[:, 1]
    ro = slice['Density']
    U = slice['Umag']
    array = slice[array]

    index = np.argsort(y)
    y = y[index]
    array = array[index]
    ro = ro[index]
    U = U[index]

    avarray = 0.0
    totmass = 0.0
    for j in range(len(y) - 1):
        dy = y[j + 1] - y[j]
        mass = 0.5 * (ro[j + 1] + ro[j]) * 0.5 * (U[j + 1] + U[j]) * dy
        totmass += mass
        avarray += mass * 0.5 * (array[j + 1] + array[j])

    avarray = avarray / totmass

    return avarray

def get_mdot_from_file(path_from_data):
    os.chdir('/home/nick/Documents/vki/data')
    os.chdir(path_from_data)
    # Extract loss value
    try:
        grid = pyvista.read('flow.vtk')
    except:
        print('diverged')
        return np.inf

    df = pd.read_csv('history.csv')
    final_res = df['Res_Flow[0]'].values[-1]
    if final_res > -8.0:
        print('diverged!')
        return np.inf

    inlet_slice = grid.slice(origin=(-0.054999, 0, 0))
    ro = inlet_slice['Density']
    Ux = inlet_slice['Momentum'][:,0]/inlet_slice['Density']
    y = inlet_slice.points[:, 1]
    totmass = 0.0
    for j in range(len(y) - 1):
        dy = y[j + 1] - y[j]
        mass = 0.5 * (ro[j + 1] + ro[j]) * 0.5 * (Ux[j + 1] + Ux[j]) * dy
        totmass += mass
    print('mdot is %f' % totmass)
    return totmass

def get_loss_from_file(path_from_data):
    os.chdir('/home/nick/Documents/vki/data')
    os.chdir(path_from_data)
    # Extract loss value
    try:
        grid = pyvista.read('flow.vtk')
    except:
        print('diverged')
        return np.inf

    df = pd.read_csv('history.csv')
    final_res = df['Res_Flow[0]'].values[-1]
    if final_res > -8.0:
        print('diverged!')
        return np.inf

    cp = 1005
    gamma = 1.4
    R = 287
    gm1 = gamma - 1.0

    T = grid['Temperature']
    p = grid['Pressure']
    ro = grid['Density']
    Mach = grid['Mach']
    dyn = Mach ** 2.0 * gm1 / 2.0

    grid['T0'] = T * (1.0 + dyn)  # Stag temp
    grid['p0'] = p * (1.0 + dyn) ** (gamma / gm1)  # Stag pres
    grid['a'] = np.sqrt(gamma * p / ro)  # speed of sound?
    grid['Umag'] = Mach * grid['a']  # absolute Speed

    inlet_slice = grid.slice(origin=(-0.054999, 0, 0))
    outlet_slice = grid.slice(origin=(.07999, 0, 0))

    T01av = massav(inlet_slice, 'T0')
    p01av = massav(inlet_slice, 'p0')
    p1av = massav(inlet_slice, 'Pressure')

    T02av = massav(outlet_slice, 'T0')
    p02av = massav(outlet_slice, 'p0')
    p2av = massav(outlet_slice, 'Pressure')

    grid['T0s'] = T01av * (grid['p0'] / p01av) ** (gm1 / gamma)

    # Outlet Re

    loss = (p01av - p02av) / (p02av - p2av)
    print('loss is %f' % loss)
    return loss

def get_Misen_from_file(path_from_data):
    os.chdir('/home/nick/Documents/vki/data')
    os.chdir(path_from_data)
    # Extract loss value
    try:
        grid = pyvista.read('flow.vtk')
    except:
        print('diverged')
        return np.inf, None

    df = pd.read_csv('history.csv')
    final_res = df['Res_Flow[0]'].values[-1]
    if final_res > -8.0:
        print('diverged!')
        return np.inf, None

    gamma = 1.4
    gm1 = gamma - 1.0

    T = grid['Temperature']
    p = grid['Pressure']
    ro = grid['Density']
    Mach = grid['Mach']
    dyn = Mach ** 2.0 * gm1 / 2.0

    inlet_slice = grid.slice(origin=(-0.054999, 0, 0))

    p01av = massav(inlet_slice, 'p0')
    surface_data = pyvista.read('surface_flow.vtk')
    pressure_data = surface_data.point_arrays['Pressure']

    x_data = surface_data.points[:,0]
    Misen = np.sqrt((2.0 / gm1) * ((p01av / pressure_data) ** (gm1 / gamma) - 1.0))
    print('done extracting Misen')
    return Misen, x_data

def get_pres_from_file(path_from_data):
    os.chdir('/home/nick/Documents/vki/data')
    os.chdir(path_from_data)
    # Extract loss value
    try:
        grid = pyvista.read('flow.vtk')
    except:
        print('diverged')
        return np.inf, np.inf

    df = pd.read_csv('history.csv')
    final_res = df['Res_Flow[0]'].values[-1]
    if final_res > -8.0:
        print('diverged!')
        return np.inf, np.inf

    cp = 1005
    gamma = 1.4
    R = 287
    gm1 = gamma - 1.0

    T = grid['Temperature']
    p = grid['Pressure']
    ro = grid['Density']
    Mach = grid['Mach']
    dyn = Mach ** 2.0 * gm1 / 2.0

    grid['T0'] = T * (1.0 + dyn)  # Stag temp
    grid['p0'] = p * (1.0 + dyn) ** (gamma / gm1)  # Stag pres
    grid['a'] = np.sqrt(gamma * p / ro)  # speed of sound?
    grid['Umag'] = Mach * grid['a']  # absolute Speed

    inlet_slice = grid.slice(origin=(-0.054999, 0, 0))
    outlet_slice = grid.slice(origin=(.07999, 0, 0))

    p02av = massav(outlet_slice, 'p0')
    p2av = massav(outlet_slice, 'Pressure')

    return p02av, p2av
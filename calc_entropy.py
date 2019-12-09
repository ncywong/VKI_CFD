from equadratures import *
import os
import numpy as np

# files = ...


def integrate(xslice, var):
    var = xslice[var]
    y = xslice.points[:, 1]
    index = np.argsort(y)
    y = y[index]
    var = var[index]
    intvar = 0.0
    for j in range(len(y) - 1):
        dy = y[j + 1] - y[j]
        if (dy > 0.0005 and xslice.points[0, 0] < 0.03695): dy = 0
        intvar += dy * 0.5 * (var[j + 1] + var[j])
    return intvar


def massav(xslice, var):
    var = xslice[var]
    y = xslice.points[:, 1]
    index = np.argsort(y)
    y = y[index]
    var = var[index]
    ro = xslice['Density'][index]
    U = xslice['Umag'][index]
    av = 0.0
    totmass = 0.0
    for j in range(len(y) - 1):
        dy = y[j + 1] - y[j]
        mass = 0.5 * (ro[j + 1] + ro[j]) * 0.5 * (U[j + 1] + U[j]) * dy
        totmass += mass
        av += mass * 0.5 * (var[j + 1] + var[j])
    av = av / totmass
    return av


def redim(vtk):
    muref = 2667.71
    roref = 6.46971
    uref = 412.338
    Tref = 592.295
    pref = 1.1e+06

    vtk['Density'] *= roref
    vtk['Momentum'] *= roref * uref
    vtk['Temperature'] *= Tref
    vtk['Pressure'] *= pref
    vtk['Laminar_Viscosity'] *= muref
    vtk['Eddy_Viscosity'] *= muref

    return vtk


Tf = 372  # bulk static temp (take as average of inflow and outflow static temp)
Prlam = 0.72
Prturb = 0.9
gamma = 1.4
R = 287.058
gm1 = gamma - 1.0
cp = (gamma / gm1) * R

# fig1, ax1 = plt.subplots()
# fig2, ax2 = plt.subplots()
# fig3, ax3 = plt.subplots()
# fig4, ax4 = plt.subplots()
# ax1.set_xlabel('x')
# ax2.set_xlabel('x')
# ax3.set_xlabel('x')
# ax4.set_xlabel('x')
# ax1.set_ylabel('$S_{ht}$')
# ax2.set_ylabel('$S_v$')
# ax3.set_ylabel('$S_t$')
# ax4.set_ylabel('$S_{tot}$')

# for file in files:

def calc_entropy(dir_name):
    os.chdir(base_dir)
    os.chdir('data')
    N = 1000
    print('\n', dir_name)
    try:
        vtk = pyvista.read(dir_name + '/flow.vtk')
    except:
        return np.inf * np.ones(N)

    # Re-dim variables on file
    vtk = redim(vtk)
    dsa = dsap.WrapDataObject(vtk)

    ro = dsa.PointData['Density']
    U = (dsa.PointData['Momentum'] / ro)
    T = dsa.PointData['Temperature']
    p = dsa.PointData['Pressure']
    mu = dsa.PointData['Laminar_Viscosity']
    mut = dsa.PointData['Eddy_Viscosity']
    conduc = mu / (ro * Prlam) + mut / (ro * Prturb)
    vtk['MomentumX'] = vtk['Momentum'][:, 0]

    # Take gradients of U and T
    Jt = algs.gradient(U)  # Jt is this one as algs uses j,i ordering
    dUdx = algs.apply_dfunc(np.transpose, Jt, (0, 2, 1))
    dTdx = algs.gradient(T)

    # Heat transfer entropy gen
    Sht = (conduc / T ** 2) * (dTdx[:, 0] ** 2 + dTdx[:, 1] ** 2 + dTdx[:, 2] ** 2)

    # Friction entropy gen
    phi = 2 * ((dUdx[:, 0, 0] ** 2 + dUdx[:, 1, 1] ** 2)) + (dUdx[:, 0, 1] + dUdx[:, 1, 0]) ** 2
    Sv = mu * phi / T
    St = mut * phi / T
    #    Sv = (mu/Tf)*( 2*(dUdx[:,0,0]**2 + dUdx[:,1,1]**2 + dUdx[:,2,2]**2) + (dUdx[:,0,1]+dUdx[:,1,0])**2 +
    #            (dUdx[:,0,2]+dUdx[:,2,0])**2 + (dUdx[:,1,2]+dUdx[:,2,1])**2)
    #    St = (mut/Tf)*( 2*(dUdx[:,0,0]**2 + dUdx[:,1,1]**2 + dUdx[:,2,2]**2) + (dUdx[:,0,1]+dUdx[:,1,0])**2 +
    #            (dUdx[:,0,2]+dUdx[:,2,0])**2 + (dUdx[:,1,2]+dUdx[:,2,1])**2)

    # Save entropy gen rates to vtk
    vtk['Sht'] = Sht
    vtk['Sv'] = Sv
    vtk['St'] = St
    vtk['Stot'] = Sht + Sv + St

    vtk.save(dir_name + '/flow_withS.vtk')

    # Integrate in y-dir and plot at each x
    points = dsa.GetPoints()
    xmin = np.min(points[:, 0]) + 0.001
    xmax = np.max(points[:, 0]) - 0.001


    xslices = np.linspace(xmin, xmax, N)
    Sht_int = np.zeros(N)
    Sv_int = np.zeros(N)
    St_int = np.zeros(N)
    Sht_xint = 0
    Sv_xint = 0
    St_xint = 0

    #    i = 341
    #    xslice = vtk.slice(origin=(xslices[i], 0, 0),normal='x')
    #    temp = integrate(xslice,'Sv')
    #    print(xslices[i],temp)
    #    quit()

    for i, x in enumerate(xslices):
        xslice = vtk.slice(origin=(x, 0, 0), normal='x')
        Sht_int[i] = integrate(xslice, 'Sht')
        Sv_int[i] = integrate(xslice, 'Sv')
        St_int[i] = integrate(xslice, 'St')

        if (i > 0):
            dx = xslices[i] - xslices[i - 1]
            Sht_xint += dx * 0.5 * (Sht_int[i - 1] + Sht_int[i])
            Sv_xint += dx * 0.5 * (Sv_int[i - 1] + Sv_int[i])
            St_xint += dx * 0.5 * (St_int[i - 1] + St_int[i])
        #            ax2.plot(xslices[i],Sv_xint,'o')

    # ax1.plot(xslices, Sht_int, label=file)
    # ax2.plot(xslices, Sv_int, label=file)
    # ax3.plot(xslices, St_int, label=file)

    # 0.03698267
    upper_x = .04
    in_chord = [i for i in range(len(xslices)) if 0 < xslices[i] < upper_x]
    Stot_int = Sht_int + Sv_int + St_int
    plot_x = np.linspace(0, upper_x / 0.03698267,len(in_chord))

    # ax4.plot(plot_x, Stot_int[in_chord], label=file)
    # ax4.set_yscale('log')
    # ax4.set_ylim([1,2500])
    # ax4.set_xlim([0,upper_x / 0.03698267])

    # ax1.legend()
    # ax2.legend()
    # ax3.legend()
    # ax4.legend()

    # print('Sht = ', Sht_xint)
    # print('Sv  = ', Sv_xint)
    # print('St  = ', St_xint)
    dStot = Sht_xint + St_xint + Sv_xint
    # print('Stot = ', dStot)

    # Get mass-averaged properties at xmin and xmax
    Mach = vtk['Mach']
    dyn = Mach ** 2.0 * gm1 / 2.0
    vtk['T0'] = T * (1.0 + dyn)  # Stag temp
    vtk['p0'] = p * (1.0 + dyn) ** (gamma / gm1)  # Stag pres
    vtk['Umag'] = np.sqrt(U[:, 0] ** 2 + U[:, 1] ** 2 + U[:, 2] ** 2)

    inlet_slice = vtk.slice(origin=(xmin, 0, 0))
    outlet_slice = vtk.slice(origin=(xmax, 0, 0))

    T01av = massav(inlet_slice, 'T0')
    p01av = massav(inlet_slice, 'p0')
    p1av = massav(inlet_slice, 'Pressure')
    T1av = massav(inlet_slice, 'Temperature')
    T02av = massav(outlet_slice, 'T0')
    p02av = massav(outlet_slice, 'p0')
    p2av = massav(outlet_slice, 'Pressure')
    T2av = massav(outlet_slice, 'Temperature')
    U2mag = massav(outlet_slice, 'Umag')

    # Pressure loss coefficient
    loss = (p01av - p02av) / (p02av - p2av)
    # print('Loss coeff is %f' % loss)

    # Entropy (gain) "loss" coefficient
    ds = -R * np.log(p02av / p01av)
    entropy_coeff = T2av * ds / (0.5 * U2mag ** 2)
    # print('entropy coeff is %f' % entropy_coeff)
    # print('ds from p02/p01 %f:' % ds)

    mdot = integrate(outlet_slice, 'MomentumX')
    ds = dStot / mdot
    # print('ds from integrating entropy gen rate %f:' % ds)

    # print('p02-p01 %f:' % (p02av - p01av))
    # print('ds*mdot %f:' % dStot)
    return Stot_int
    # print('mdot %f:' % mdot)

# Module contains classes with launch vehicle data
# Inclusion of launch vehicle performance model
# Data collected from publically available sources

# TODO ==== add additional US launchers: Delta-IV Heavy, Atlav V 551 etc.

import numpy as np
import time
import matplotlib.pyplot as plt

# Launch Vehicle Data

# Data: Europe

class SOYUZ():
    # data source: http://www.arianespace.com/wp-content/uploads/2015/09/Soyuz-Users-Manual-March-2012.pdf
    name = "Soyuz - ArianeSpace"
    C3 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]
    performance = [2180, 2118.7, 2060, 2005, 1950, 1900, 1851.7, 1805, 1755, 1703.5, 1659.9, 1616.2, 1577.6, 1534,
                   1490.3, 1445.7, 1405, 1368.5, 1328.9, 1289.3, 1253.8, 1213.2, 1176.6, 1138.1, 1098.4, 1061.9]
    params = [-1.38369063e-02, 9.49540276e-01, -5.97332462e+01, 2.17749685e+03]

class A64():
    # data source:
    name = "Ariane 64 - ArianeSpace"
    C3 = []
    performance = []
    params = []

class PROTON():
    # data source: http://www.ilslaunch.com/sites/default/files/pdf/Proton%20Mission%20Planner%27s%20Guide%20Revision%207%20%28LKEB-9812-1990%29.pdf
    name = "Proton - ILS"
    C3 = [0, 1, 4, 9, 16, 25, 36]
    performance = [6475, 6355, 6002,  5454,  4745,  3971,  3111]
    params = [-8.53221832e-03, 1.15478692e+00, -1.24005320e+02, 6.47763737e+03]

# Data: USA

class SLS_1B():
    # data source: https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20170005323.pdf
    name = "SLS Block 1 B"
    C3 = [0, 10, 20, 30, 40, 50, 60, 70, 80, 83, 90, 100, 110, 120]
    performance = [36739, 31450, 26787, 22605, 18876, 15562, 12618, 10004, 7678, 7031, 5605, 3753, 2096, 607]
    params = [-6.59136287e-03, 2.86530512e+00, -5.49937736e+02, 3.67067181e+04]

class SLS_2():
    # data source: https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20170005323.pdf
    name = "SLS Block 2"
    C3 = [0, 10, 20, 30, 40, 50, 60, 70, 80, 83, 90, 100, 110, 120, 130]
    performance = [44306, 37626, 31739, 26597, 22125, 18237, 14852, 11893, 9299, 8583, 7014, 4994, 3200, 1602, 171]
    params = [-1.03389681e-02, 4.12421189e+00, -7.01074973e+02, 4.42499654e+04]

class F9():
    # data source: https://www.spaceflightnow.com/falcon9/001/f9guide.pdf
    # TODO ==== update with block 5 data
    name = "Falcon 9 (block 2) - SpaceX"
    C3 = [0, 4, 7, 11, 14, 19, 23, 28, 33, 39, 45]
    performance = [2473, 2248, 2023, 1798, 1573, 1348, 1123, 898, 673, 448, 223]
    params = [-9.54151620e-04, 4.65690789e-01, -6.92797732e+01, 2.48933308e+03]

class FH_R():
    # data source: https://elvperf.ksc.nasa.gov/Pages/Query.aspx
    name = "Falcon Heavy (recovered) - SpaceX"
    C3 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 40]
    performance = [5830, 5675, 5520, 5365, 5205, 5050, 4895, 4740, 4585, 4430, 4275, 3935, 3595, 3280, 2970, 2680,
                2395, 2130, 1870, 930]
    params = [-5.13090812e-03, 1.35199881e+00, -1.68756191e+02, 5.84881584e+03]

class FH_E():
    # data source: https://elvperf.ksc.nasa.gov/Pages/Query.aspx
    name = "Falcon Heavy (expended) - SpaceX"
    C3 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 40, 50, 60]
    performance = [12430, 12170, 11915, 11655, 11400, 11140, 10880, 10625, 10365, 10110, 9850, 9295, 8735, 8220, 7710,
                7235, 6760, 6320, 5885, 4700, 3010, 1850]
    params = [-1.86786978e-02, 3.01282040e+00, -2.91177291e+02, 1.24902721e+04]

class ATLAS_V401():
    # data source: https://elvperf.ksc.nasa.gov/Pages/Query.aspx
    name = "Atlas V (401) - ULA"
    C3 = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45]
    performance = [3035, 2720, 2425, 2145, 1880, 1630, 1400, 1185, 985, 800]
    params = [-1.32090132e-04, 3.32400932e-01, -6.43393551e+01, 3.03464336e+03]

class ATLAS_V551():
    # data source: https://elvperf.ksc.nasa.gov/Pages/Query.aspx
    name = "Atlas V (551) - ULA"
    C3 = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60]
    performance = [6105, 5570, 5060, 4585, 4140, 3730, 3345, 2995, 2670, 2380, 2120, 1910, 1695]
    params = [-9.28983210e-17, 6.20579421e-01, -1.10652348e+02, 6.10576923e+03]

class DELTA_IVH():
    # data source: https://elvperf.ksc.nasa.gov/Pages/Query.aspx
    name = "Delta IV (Heavy) - ULA"
    C3 = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 80, 100]
    performance = [10185, 9285, 8460, 7695, 6995, 6350, 5755, 5205, 4700, 4225, 3790, 3395, 3000, 1710, 705]
    params = [-4.36559445e-03, 1.31450361e+00, -1.82528619e+02, 1.01665722e+04]

# Data: JAPAN

class HIIA202():
    # https://www.mhi.com/jp/products/pdf/manual.pdf
    name = "H-IIA 202 - JAXA"
    C3 = [0, 10, 20, 40]
    performance = [2592, 1940, 1372, 464]
    params = [-1.000e-03, 4.500e-01, -6.960e+01, 2.592e+03]

# LAUNCH VEHICLE PERFORMANCE EVALUATION

# launch vehicle performance
def lv_performance(launch_vehicle = F9, C3 = 2.0):
    # evaluates launch vehicle performance for a given C3
    params = launch_vehicle.params
    performance = sum([A * C3 ** i for A, i in zip(reversed(params), range(len(params)))])
    return performance

# launch vehicle performance - curve parameters
def get_lv_params(launch_vehicle = F9, C3 = 2.0, poly = 3):
    params = np.polyfit(launch_vehicle.C3, launch_vehicle.performance, poly)
    return params

# PLOT LAUNCHER MODEL

def plot_lv_performance(lvs):

    # Array of launch vehicles to plot
    lvs = lvs

    # create figure
    fig, ax = plt.subplots(figsize = (10,10))



    # color map range
    n = len(lvs)
    colors = plt.cm.jet(np.linspace(0,1,n))

    # data generation and plot creation
    for lv, n in zip(lvs,range(n)):
        x = []
        y = []
        for C3 in np.arange(-5,36,1):
            x.append(C3)
            y.append(lv_performance(lv,C3))
        ax.plot(x,y, label = lv.name, color = colors[n])

    # add legend
    ax.legend(loc='upper right')

    # adjust plot limits
    plt.xlim(-5,35)
    plt.ylim(0, 12000)

    # settings for plot grid lines amd ticks
    plt.grid(which='major',axis='both',color='gainsboro', linestyle='-')
    plt.grid(b = True, which='minor', axis='both', color='gainsboro', linestyle='-')
    ax.spines['top'].set_visible(True)
    ax.spines['right'].set_visible(True)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_visible(True)
    ax.spines['bottom'].set_color('dimgray')
    ax.spines['left'].set_color('dimgray')
    ax.tick_params(axis='x', colors='dimgray')
    ax.tick_params(axis='y', colors='dimgray')
    ax.spines['bottom'].set_color('dimgray')
    ax.spines['top'].set_color('dimgray')
    ax.spines['right'].set_color('dimgray')
    ax.spines['left'].set_color('dimgray')

    # add plot title
    fig.suptitle('MERLOT launcher performance model', fontsize=14, x = 0.5, y = 0.925)

    # add axis title
    plt.xlabel('C3 (km2/s2)', fontsize=12)
    plt.ylabel('Injected mass (kg)', fontsize=12)

    plt.show()

plot_lv_performance([SOYUZ,F9,HIIA202,ATLAS_V401,FH_R,ATLAS_V551,PROTON,DELTA_IVH,FH_E,SLS_1B,SLS_2])

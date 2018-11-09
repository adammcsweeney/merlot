# Module contains classes with launch vehicle data (data collected from publicly available sources)
# Note: C3 v performance relationship considered linear for extrapolated data points beyond stated in source data

# TODO ==== add launcher data from Robic Biesbroek - Lunar and Interplanetary Trajectories
# TODO ==== compare data with existing in file, update for launchers
# TODO ==== produce better model than linear extrapolation for modelling performance curve at higher isp

import numpy as np
import matplotlib.pyplot as plt

# Launch Vehicle Data

# Data: Europe

class VEGA():
    # data source "Lunar and Interplanetary Trajectories - Robin Biesbroek"
    name = "Vega - ArianeSpace"
    C3 = []
    performance = []
    params = []

class SOYUZ():
    # data source: http://www.arianespace.com/wp-content/uploads/2015/09/Soyuz-Users-Manual-March-2012.pdf
    name = "Soyuz - ArianeSpace"
    shortname = 'soyuz'
    C3 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 30, 35, 40,
          50, 54]
    performance = [2180.0, 2118.7, 2060.0, 2005.0, 1950.0, 1900.0, 1851.7, 1805.0, 1755.0, 1703.5, 1659.9, 1616.2,
                   1577.6, 1534.0, 1490.3, 1445.7, 1405.0, 1368.5, 1328.9, 1289.3, 1253.8, 1213.2, 1176.6, 1138.1,
                   1098.4, 1061.9, 879.4, 696.9, 514.4, 149.4, 3.4]
    params = [-5.53513526e-03, 5.79861593e-01, -5.54029246e+01, 2.16773719e+03]

class PROTON():
    # data source: http://www.ilslaunch.com/sites/default/files/pdf/Proton%20Mission%20Planner%27s%20Guide%20Revision%207%20%28LKEB-9812-1990%29.pdf
    name = "Proton - ILS"
    shortname = 'proton'
    C3 = [0.0, 0.3, 1.0, 2.3, 4.0, 6.3, 9.0, 12.3, 16.0, 20.3, 25.0, 30.3, 36.0, 40.0, 50.0, 60.0, 70.0, 75.0, 76.0]
    performance = [6475, 6445, 6355, 6205, 6002, 5748, 5454, 5124, 4745, 4361, 3971, 3556, 3111, 2801.46, 2027.55,
                   1253.64, 479.73, 92.775, 15.384]
    params = [-6.17046459e-03, 9.13323799e-01, -1.18922266e+02, 6.46782776e+03]

# Data: USA

class SLS_1B():
    # data source: https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20170005323.pdf
    name = "SLS Block 1 B"
    shortname = 'SLS1B'
    C3 = [0, 10, 20, 30, 40, 50, 60, 70, 80, 83, 90, 100, 110, 120, 124]
    performance = [36739.00, 31450.00, 26787.00, 22605.00, 18876.00, 15562.00, 12618.00, 10004.00, 7678.00, 7031.00,
                   5605.00, 3753.00, 2096.00, 607.00, 11.4]
    params = [-6.57563021e-03, 2.86290854e+00, -5.49846816e+02, 3.67061828e+04]

class SLS_2():
    # data source: https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20170005323.pdf
    name = "SLS Block 2"
    shortname = 'SLS2'
    C3 = [0, 10, 20, 30, 40, 50, 60, 70, 80, 83, 90, 100, 110, 120, 130]
    performance = [44306, 37626, 31739, 26597, 22125, 18237, 14852, 11893, 9299, 8583, 7014, 4994, 3200, 1602, 171]
    params = [-1.03389681e-02, 4.12421189e+00, -7.01074973e+02, 4.42499654e+04]

class F9():
    # data source: https://www.spaceflightnow.com/falcon9/001/f9guide.pdf
    # TODO ==== update with block 5 data
    name = "Falcon 9 (block 2) - SpaceX"
    shortname = 'falcon9'
    C3 = [0, 4, 7, 11, 14, 19, 23, 28, 33, 39, 45]
    performance = [2473, 2248, 2023, 1798, 1573, 1348, 1123, 898, 673, 448, 223]
    params = [-9.54151620e-04, 4.65690789e-01, -6.92797732e+01, 2.48933308e+03]

class FH_R():
    # data source: https://elvperf.ksc.nasa.gov/Pages/Query.aspx
    name = "Falcon Heavy (recovered) - SpaceX"
    shortname = 'falconheavy-rec'
    C3 = [0, 10, 20, 30, 40, 50, 60, 64]
    performance = [6690, 5130, 3845, 2740, 1805, 1005, 320, 46]
    params = [-6.01040586e-03, 1.38697478e+00, -1.67836559e+02, 6.68648485e+03]

class FH_E():
    # data source: https://elvperf.ksc.nasa.gov/Pages/Query.aspx
    name = "Falcon Heavy (expended) - SpaceX"
    shortname = 'falconheavy-exp'
    C3 = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    performance = [15010, 12345, 10115, 8230, 6640, 5280, 4100, 3080, 2060, 1040, 20]
    params = [-9.83877234e-03, 2.36468531e+00, -2.88051088e+02, 1.50063636e+04]

class ATLAS_V401():
    # data source: https://elvperf.ksc.nasa.gov/Pages/Query.aspx
    name = "Atlas V (401) - ULA"
    shortname = 'atlasv401'
    C3 = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65]
    performance = [3035, 2720, 2425, 2145, 1880, 1630, 1400, 1185, 985, 800, 615, 430, 245, 60]
    params = [-2.75126834e-03, 5.09869542e-01, -6.73110443e+01, 3.04213235e+03]

class ATLAS_V551():
    # data source: https://elvperf.ksc.nasa.gov/Pages/Query.aspx
    name = "Atlas V (551) - ULA"
    shortname = 'atlasv551'
    C3 = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 99]
    performance = [6105, 5570, 5060, 4585, 4140, 3730, 3345, 2995, 2670, 2380, 2120, 1910, 1695, 1265, 835, 405, 18]
    params = [-4.22395608e-03, 9.97234264e-01, -1.19098875e+02, 6.13690297e+03]

class DELTA_IVH():
    # data source: https://elvperf.ksc.nasa.gov/Pages/Query.aspx
    name = "Delta IV (Heavy) - ULA"
    shortname = 'deltaiv-heavy'
    C3 = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 80, 100, 114]
    performance = [10185, 9285, 8460, 7695, 6995, 6350, 5755, 5205, 4700, 4225, 3790, 3395, 3000, 1710, 705, 1.5]
    params = [-4.26489744e-03, 1.30121464e+00, -1.82097720e+02, 1.01641693e+04]

# Data: JAPAN

class HIIA202():
    # https://www.mhi.com/jp/products/pdf/manual.pdf
    name = "H-IIA 202 - JAXA"
    shortname = 'HIIA202'
    C3 = [0, 10, 20, 40, 50]
    performance = [2592.40, 1939.47, 1371.54, 463.62, 9.2]
    params = [-6.56818627e-03, 7.61948950e-01, -7.33579398e+01, 2.59494370e+03]

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

# print all params
# for lv in [SOYUZ, PROTON, SLS_1B, SLS_2, F9, FH_R, FH_E, ATLAS_V401, ATLAS_V551,DELTA_IVH, HIIA202]:
#     print(f"{lv.name}\t{get_lv_params(launch_vehicle=lv)}")

# PLOT LAUNCHER MODEL
def plot_lv_performance(lvs,C3_range=[0,75]):

    # Array of launch vehicles to plot
    lvs = lvs

    # create figure
    fig, ax = plt.subplots(figsize = (8,8))

    # color map range
    n = len(lvs)
    colors = plt.cm.jet(np.linspace(0,1,n))

    # data generation and plot creation
    for lv, n in zip(lvs,range(n)):
        x = []
        y = []
        for C3 in np.arange(min(C3_range),max(C3_range)+1,1):
            x.append(C3)
            y.append(lv_performance(lv,C3))
        ax.plot(x,y, label = lv.name, color = colors[n])

    # add legend
    ax.legend(loc='upper right')

    # # adjust plot limits
    plt.xlim(min(C3_range),max(C3_range))
    plt.ylim(0)

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

plot_lv_performance(lvs=[SOYUZ,F9,HIIA202,ATLAS_V401,FH_R,ATLAS_V551,PROTON,DELTA_IVH],C3_range=[0,75])
plot_lv_performance(lvs=[SOYUZ,F9,HIIA202,ATLAS_V401,FH_R,ATLAS_V551,PROTON,DELTA_IVH,FH_E,SLS_1B,SLS_2],C3_range=[0,40])








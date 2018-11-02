# Module contains classes including EP thruster data
# Note: data currently based o estimates from publically available sources (see comments for source data)

# TODO ==== update thrusters included
# TODO ==== figure out why higher power engines do not work
# TODO ==== check acceptability of polyfit solutions (i.e. variation in isp from source data)
# TODO ==== experiment, set power input profile also to power cubed (i.e. as in MARGO, and not JPL model)


import numpy as np
import matplotlib.pyplot as plt


# THRUSTER DATA


# EUROPEAN THRUSTERS


class T6:
    # source data
    # http://arc.aiaa.org | DOI: 10.2514/1.B34173

    # data
    name = "T6"
    power = [2430, 3160, 3920, 4500]
    thrust = [0.074, 0.099, 0.123, 0.143]
    isp = [3710, 3940, 4080, 4120]
    eff = [0.55, 0.60, 0.62, 0.64]
    throughput = "unknown"

    # modelling parameters
    poly = 3
    T_c = [1.91177665e-12, -1.99713503e-08, 1.00827186e-04, -8.05131412e-02]
    Isp_c =[8.79420628e-10, -9.61874290e-05, 8.32028893e-01, 2.24352822e+03]


class HT20k:
    # source data:
    # http://epic-src.eu/wp-content/uploads/24_EPICWorkshop2017_SITAEL_EPIC_WS_1b.pdf

    # data
    name = 'HT20K'
    power = [10000, 20000]
    thrust = [0.3, 1.1]
    isp = [2000, 3800]
    throughput = 'unknown'
    eff = [.7, .7]

    # modelling parameters
    poly = 2
    T_c = [8.e-05, -5.e-01]
    Isp_c = [1.8e-01, 2.0e+02]


# AMERICAN THRUSTERS


class XIPS:
    # source data
    # https://arc.aiaa.org/doi/abs/10.2514/6.2008-4914

    # data
    name = "XIPS"
    power = [439,1995,3002,4009,5117]
    thrust = [0.015,0.072,0.101,0.140,0.174]
    isp = [1665,3291,3361,3416,3487]
    throughput = "unknown"

    # modelling parameters
    poly = 3
    T_c = [ 7.76150193e-14, -9.85236348e-10,  3.72957509e-05, -9.20265575e-04]
    Isp_c = [ 6.13263037e-08, -6.64082847e-04,  2.32488758e+00,  7.72555989e+02]


class NSTAR:
    # source data
    # https://ieeexplore.ieee.org/document/878371

    # data
    name = "NSTAR"
    power = [490, 860, 1240, 1590, 1950, 2310]
    thrust = [0.021, 0.031, 0.047, 0.062, 0.077, 0.094]
    isp = [1939, 2753, 3039, 3132, 3177, 3195]
    eff = [0.39, 0.48, 0.56, 0.60, 0.61, 0.64]
    throughput = "unknown"

    # modelling parameters
    poly = 3
    T_c = [-3.80038872e-12, 2.08717752e-08, 7.16087245e-06, 1.26546815e-02]
    Isp_c = [ 5.47940922e-07, -2.96015772e-03, 5.30269375e+00, -1.69255737e+00]


class BPT4000:
    # source data
    # https://www.researchgate.net/publication/272830179_Robotic_Mars_Exploration_Trajectories_Using_Hall_Thrusters

    # data (estimated from source)
    name = "BPT-4000"
    power = [302,350,400,450,500,600,1000,1500,2000,2500,3000,3500,4000,4500,4839]
    thrust = [0.015,0.018,0.022,0.026,0.029,0.037,0.066,0.100,0.131,0.160,0.187,0.213,0.238,0.263,0.281]
    isp = [650,724,786,838,883,957,1156,1322,1452,1560,1649,1721,1780,1831,1865]
    eff = []
    throughput = 580
    # throughput reference:
    # http://richard.hofer.com/pdf/iepc-2009-085.pdf

    # modelling parameters
    # TODO ==== fix parameters
    poly = 3
    # T_c = [4.71936824e-23, -7.41236203e-19, 4.59526182e-15, -1.35754920e-11, 1.44640126e-08, 6.78685611e-05, -6.70562752e-03]
    # Isp_c = [-2.96518852e-18, 5.11251303e-14, -3.48693175e-10, 1.19806054e-06, -2.21379650e-03, 2.36467728e+00, 1.19556548e+02]

    T_c = [6.33727204e-13, -8.04029989e-09, 8.42169901e-05, -1.04841379e-02]
    Isp_c = [1.92395714e-08, -2.05440557e-04, 8.33121283e-01, 4.81872636e+02]


class NEXT_HighThrust:
    # source data
    # https://arc.aiaa.org/doi/abs/10.2514/6.2008-4914

    # data (estimated from source)
    name = "NEXT (high-thrust)"
    power = [485,2005,4018,5995,7542]
    thrust = [0.023,0.076,0.156,0.212,0.239]
    isp = [1374,3086,3173,3723,4186]
    eff = []
    throughput = "unknown"

    # modelling paramters
    poly = 3
    T_c = [-4.04741615e-13, 3.02178551e-09, 3.09517010e-05, 6.76930978e-03]
    Isp_c = [2.83542264e-08, -3.87892632e-04, 1.78869817e+00, 6.52196584e+02]


class NEXT_HighIsp:
    # source data
    # https://arc.aiaa.org/doi/abs/10.2514/6.2008-4914

    # data (estimated from source)
    name = "NEXT (high-isp)"
    power = [485, 2005, 4018, 5995, 7542]
    thrust = [0.023,0.069,0.125,0.196,0.242]
    isp = [1076,3440,4257,4131,4257]
    eff = []
    throughput = "unknown"

    # modelling paramters
    poly = 3
    T_c = [-9.51248864e-14, 1.52765068e-09, 2.45492623e-05, 1.16399821e-02]
    Isp_c = [3.31864848e-08, -5.29086442e-04, 2.68081192e+00, -9.66900768e+01]

plot_compare_thrust([NSTAR,T6,BPT4000, XIPS, NEXT_HighIsp, NEXT_HighThrust])
plot_compare_isp([NSTAR,T6,BPT4000, XIPS, NEXT_HighIsp, NEXT_HighThrust])


# EVALUATING THRUSTER PERFORMANCE


# script for polynomial coefficient calculation
def get_thruster_params(engine):
    # thrust coefficients
    Thr_coeffs = np.polyfit(engine.power,engine.thrust,engine.poly)
    print(Thr_coeffs)
    # isp coefficients
    Isp_coeffs = np.polyfit(engine.power, engine.isp, engine.poly)
    print(Isp_coeffs)

def eval_thrust(thruster, power):
    params = thruster.T_c
    thrust = sum([A * power ** i for A, i in zip(reversed(params), range(len(params)))])
    return thrust

def eval_isp(thruster, power):
    params = thruster.Isp_c
    isp = sum([A * power ** i for A, i in zip(reversed(params), range(len(params)))])
    return isp


# PLOTTING THRUSTER PERFORMANCE


def plot_compare_thrust(thrusters):
    # Array of launch vehicles to plot
    thrusters = thrusters

    # create figure
    fig, ax = plt.subplots(figsize=(8, 8))

    # color map range
    n = len(thrusters)
    colors = plt.cm.jet(np.linspace(0, 1, n))

    # data generation and plot creation
    for thruster, n in zip(thrusters, range(n)):
        x = []
        y = []
        for power in np.arange(min(thruster.power), max(thruster.power)+1, 1):
            x.append(power)
            y.append(eval_thrust(thruster,power))
        ax.plot(x, y, label=thruster.name, color=colors[n])

    # add legend
    ax.legend(loc='upper left')

    # # adjust plot limits
    # plt.xlim(0,7500)
    # plt.ylim(0, 0.3)

    # settings for plot grid lines amd ticks
    plt.grid(which='major', axis='both', color='gainsboro', linestyle='-')
    plt.grid(b=True, which='minor', axis='both', color='gainsboro', linestyle='-')
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
    fig.suptitle('MERLOT engine thrust comparison', fontsize=14, x=0.5, y=0.925)

    # add axis title
    plt.xlabel('Power (W)', fontsize=12)
    plt.ylabel('Thrust (N)', fontsize=12)

    plt.show()

def plot_compare_isp(thrusters):
    # Array of launch vehicles to plot
    thrusters = thrusters

    # create figure
    fig, ax = plt.subplots(figsize=(8, 8))

    # color map range
    n = len(thrusters)
    colors = plt.cm.jet(np.linspace(0, 1, n))

    # data generation and plot creation
    for thruster, n in zip(thrusters, range(n)):
        x = []
        y = []
        for power in np.arange(min(thruster.power), max(thruster.power)+1, 1):
            x.append(power)
            y.append(eval_isp(thruster,power))
        ax.plot(x, y, label=thruster.name, color=colors[n])

    # add legend
    ax.legend(loc='upper left')

    # # adjust plot limits
    # plt.xlim(0,7500)
    # plt.ylim(1000, 4000)

    # settings for plot grid lines amd ticks
    plt.grid(which='major', axis='both', color='gainsboro', linestyle='-')
    plt.grid(b=True, which='minor', axis='both', color='gainsboro', linestyle='-')
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
    fig.suptitle('MERLOT engine isp comparison', fontsize=14, x=0.5, y=0.925)

    # add axis title
    plt.xlabel('Power (W)', fontsize=12)
    plt.ylabel('isp (s)', fontsize=12)

    plt.show()













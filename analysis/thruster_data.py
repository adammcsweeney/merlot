# Module contains classes including EP thruster data
# Note: data currently based o estimates from publically available sources (see comments for source data)
# Note: linear relationships currently assumed between power and thrust/isp

# TODO ==== update thrusters included
# TODO ==== curves for thruster data
# TODO ==== figure out why higher power engines do not work

import numpy as np
import matplotlib.pyplot as plt

# EUROPEAN THRUSTERS

class HT20k():
    # source data:
    # http://epic-src.eu/wp-content/uploads/24_EPICWorkshop2017_SITAEL_EPIC_WS_1b.pdf
    name = 'HT20K'
    power = [10000, 20000]
    thrust = [0.3, 1.1]
    isp = [2000, 3800]
    throughput = 'unknown'
    eff = [.7, .7]
    poly = 1
    T_c = [8.e-05, -5.e-01]
    Isp_c = [1.8e-01, 2.0e+02]

# AMERICAN THRUSTERS

class NSTAR():
    # source data: https://pdssbn.astro.umd.edu/holdings/ds1-c-micas-3-rdr-visccd-borrelly-v1.0/document/doc_Apr04/int_reports/IPS_Integrated_Report.pdf
    name = 'NSTAR'
    power = [423, 2288]
    thrust = [0.019, 0.093]
    isp = [1814, 3127]
    throughput = 125
    poly = 1
    T_c = [0.00004, 0.002]
    Isp_c = [0.07, 1516.199]

class BPT4000():
    # source data:
    # https://www.researchgate.net/profile/Jon_Sims/publication/228909022_Benefits_of_Using_Hall_Thrusters_for_a_Mars_Sample_Return_Mission/links/55f1985308aedecb6900e3f2.pdf
    # http://erps.spacegrant.org/uploads/images/images/iepc_articledownload_1988-2007/2009index/IEPC-2009-144.pdf
    name = 'BPT-4000'
    power = [300, 4500]
    thrust = [0.017, 0.252]
    isp = [2060, 2060]
    throughput = 250
    poly = 1
    T_c = [0.00006, 0.0002]
    Isp_c = [0, 2060]




# script for polynomial coefficient calculation
def get_thruster_params(engine = SPT140):
    # thrust coefficients
    Thr_coeffs = np.polyfit(engine.power,engine.thrust,engine.poly)
    print(Thr_coeffs)
    # isp coefficients
    Isp_coeffs = np.polyfit(engine.power, engine.isp, engine.poly)
    print(Isp_coeffs)

def eval_thrust(thruster = SPT140, power = 2500):
    params = thruster.T_c
    thrust = sum([A * power ** i for A, i in zip(reversed(params), range(len(params)))])
    return thrust

def eval_isp(thruster = SPT140, power = 2500):
    params = thruster.Isp_c
    isp = sum([A * power ** i for A, i in zip(reversed(params), range(len(params)))])
    return isp


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

plot_compare_isp([HT20k])
plot_compare_thrust([HT20k])

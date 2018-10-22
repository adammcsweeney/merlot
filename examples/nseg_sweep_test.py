# example script for parametric sweep
# evaluates optimisation run for similar input parameters, varying n_seg



from merlot.interplanetary import merlot_sep_e2pl
import pygmo as pg
from matplotlib import pyplot as plt
from pykep.examples import add_gradient
import pygmo_plugins_nonfree as ppn
import merlot.physical_data as pd
import merlot.thruster_data as td
import merlot.launcher_data as lv
from merlot.orbit import spiral_in, spiral_out
import pykep as pk
import numpy as np

import time

# 0 - set-up

# computation start time
start = time.time()

# initial run number
ID = 1

# USER INPUTS
planet = pd.Mars
_power = 9000
_tof = 450
_n_seg = 10
engine = td.BPT4000
no_engine = 1
LV = lv.F9
target_alt = 320

# 1 - Algorithm

# WORHP as user-defined algorithm
uda = ppn.worhp(screen_output=True,library=r'C:\Users\andrew mcsweeney\Documents\Project_Orbiter_Tool\worhp.dll')

# WORHP setup
uda.set_numeric_option('TolFeas', 1e-3)
uda.set_numeric_option('AcceptTolFeas',1e-4)
uda.set_numeric_option('AcceptTolOpti',1e-4)
uda.set_integer_option('MaxIter', 500)

for _n_seg in np.arange(3,101,1):

    try:
        # algorithm definition
        algo = pg.algorithm(uda)

        # 2 - Problem

        udp = add_gradient(merlot_sep_e2pl(planet,
                                           n_seg=int(_n_seg),
                                           year=2026,
                                           tof=int(_tof),
                                           Psa=int(_power),
                                           engine=engine,
                                           no_engine=no_engine,
                                           LV=LV,
                                           vinf_arr=1e-5,
                                           sep=True),
                           with_grad=False)

        prob = pg.problem(udp)

        prob.c_tol = [1e-5] * prob.get_nc()

        # 3 - Population
        pop = pg.population(prob, 1, 200)
        print(pop.get_seed())

        # 4 - Solve the problem (evolve)
        pop = algo.evolve(pop)

        # 5 - Inspect the solution
        print("Feasible?:", prob.feasibility_x(pop.champion_x))
        if prob.feasibility_x(pop.champion_x):

            # 6 - Spiraling to target orbit

            spirals = spiral_in(pl=planet,
                                orbit_alt=target_alt,
                                Psa=int(_power),
                                i_init=0,
                                i_fin=0,
                                engine=engine,
                                no_engine=no_engine,
                                epoch=pop.champion_x[0] + int(_tof),
                                mi=pop.champion_x[1])

            # Data collection on spirals
            sp_T, sp_isp = spirals.sep_model()
            sp_t, sp_mp = spirals.spiral_prop_time()
            sp_dv = spirals.spiral_dv()
            spiral_data = sp_T, sp_isp, sp_t, sp_dv, pk.epoch(pop.champion_x[0] + int(_tof) + sp_t), sp_mp

            end = time.time()

            # 7 - Write data to csv file
            udp.udp_inner.data_reporting(x=pop.champion_x, runtime=(end - start), sp_data=spiral_data, ID=ID)

        #
        # # 8 - plot trajectory
        # axis = udp.udp_inner.plot_traj(pop.champion_x, plot_segments=True)
        # plt.title("The trajectory in the heliocentric frame")
        #
        # # 9 - plot control
        # udp.udp_inner.plot_dists_thrust(pop.champion_x)
        #
        # # 10 - pretty
        # udp.udp_inner.pretty(pop.champion_x)
        #
        # # 11 - pretty - f
        # udp.udp_inner.pretty_f(pop.champion_f)
        #
        # plt.show()


    except:
        pass

    ID += 1


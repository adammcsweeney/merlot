# example script for parametric sweep
# evaluates optimisation run for similar input parameters, varying power and time of flight

from merlot.interplanetary import merlot_sep_e2pl
import pygmo as pg
from pykep.examples import add_gradient
import pygmo_plugins_nonfree as ppn
import merlot.physical_data as pd
import merlot.thruster_data as td
import merlot.launcher_data as lv
from merlot.orbit import spiral_in, spiral_out
import pykep as pk
import numpy as np
from merlot.worhp_setup import worhp_setup
import examples.data_store_writer as dsw

import time

# 1 - set-up

# USER INPUTS

# number of leg segments
# _n_seg = 10
# target planet
_planet = pd.Mars
# departure year
_year = [2020,2022,2024]
# EP thruster
_engine = td.RIT_A
# number of active EP thrusters
_no_engine = 1
# power at 1 AU
_power = 15000
# time of flight (heliocentric)
_tof = 450
# launch vehicle
_LV = lv.F9
# target altitude at target (km)
target_alt = 320

# create new file to store run data
dsw.outbound_datafile_make(_planet.name, _year, _engine.name,_no_engine,_LV.shortname)

# initial run number
# (can be set to continue data run)
ID = 1

# 2 - Algorithm

# WORHP as user-defined algorithm
uda = ppn.worhp(screen_output=False,
                library=r'C:\Users\andrew mcsweeney\Documents\Project_Orbiter_Tool\worhp.dll')

# add WORHP user settings
worhp_setup(uda)

# algorithm definition
algo = pg.algorithm(uda)

for _n_seg in np.arange(3,101,1):

    # computation start time
    end = time.time()
    start = time.time()



    try:

        # 3 - Optimisation problem definition

        # user defined problem (taken from inputs)
        udp = add_gradient(merlot_sep_e2pl(target=_planet,
                                           n_seg=int(_n_seg),
                                           year=_year,
                                           tof=int(_tof),
                                           Psa=int(_power),
                                           engine=_engine,
                                           no_engine=_no_engine,
                                           LV=_LV,
                                           vinf_arr=1e-5,
                                           sep=True),
                           with_grad=False)

        prob = pg.problem(udp)

        # tolerance for examining solution feasibility
        prob.c_tol = [1e-5] * prob.get_nc()

        # 3 - Population definition
        # set SEED = int for manually configuring the initial guess
        pop = pg.population(prob, 1, 2412496532)
        print(pop.get_seed())

        # 4 - Solve the problem (evolve)
        pop = algo.evolve(pop)

        # 5 - Inspect the solution
        print("Feasible?:", prob.feasibility_x(pop.champion_x))

        # continue with orbit spiralling and data writing if solution is feasible
        if prob.feasibility_x(pop.champion_x):

            # 6 - Spiraling to target orbit

            spirals = spiral_in(pl=_planet,
                                orbit_alt=target_alt,
                                Psa=int(_power),
                                i_init=0,
                                i_fin=0,
                                engine=_engine,
                                no_engine=_no_engine,
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


    except:
        pass

    ID += 1


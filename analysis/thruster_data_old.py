



# Module contains classes including EP thruster data
# Note: data currently based o estimates from publically available sources (see comments for source data)
# Note: linear relationships currently assumed between power and thrust/isp

# TODO ==== update thrusters included
# TODO ==== curves for thruster data

import numpy as np

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

# # plt.show()
# engine = BPT4000
#
# # thrust coefficients
# Thr_coeffs = np.polyfit(engine.power,engine.thrust,engine.poly)
# print(Thr_coeffs)
# # isp coefficients
# Isp_coeffs = np.polyfit(engine.power, engine.isp, engine.poly)
# print(Isp_coeffs)















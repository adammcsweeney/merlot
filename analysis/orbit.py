# Module for evaluating orbital low-thrust spiralling between departure/arrival altitude and target orbit altitude
# Arrival/departure altitude set as hill sphere
# Spiral model is analytical. Currently based on the Edelbaum approximation. Assumptions seen in source data.
# Source Data: https://arc.aiaa.org/doi/pdf/10.2514/1.51024

# TODO === implement thrust-off for eclipses
# TODO ==== implement thrust-off for solar conjunction
# TODO ==== include target J2 effect (is this necessary)
# TODO ==== check validity of arrival/departure orbit altitude
# TODO ==== check validity of analytical method.
# TODO ==== implement spiral in/out as a single class with bool for in/out
# TODO ==== add spiral plots
# TODO ==== add class for data writing

import math
import pykep as pk
import merlot.thruster_data as td
import merlot.physical_data as pd

class spiral_in():

    def __init__(self,
                 pl = pd.Mars,
                 orbit_alt = 320,
                 i_init = 0,
                 i_fin = 0,
                 Psa = 10000,
                 engine = td.BPT4000,
                 no_engine = 1,
                 epoch = 1000,
                 mi = 4000):

        # 1 - Class Data Members

        # planet and orbit
        self.pl = pl
        self.target_orbit = orbit_alt

        # initial and final inclination
        self.i_init = i_init
        self.i_fin = i_fin

        # spacecraft initial mass
        self.mi = mi

        # solar arrays
        self.Psa = Psa

        # propulsion subsystem
        # thrusters
        self.engine = engine
        self.no_engine = no_engine
        # max performance
        Thr_max = max(engine.thrust)
        Isp_max = max(engine.isp)
        Tmax = Thr_max * no_engine

        # planet
        planet = pk.planet.jpl_lp(self.pl.name)
        # helio-radius and velocity at target epoch
        r, v = planet.eph(pk.epoch(epoch))
        # solar distance
        self.r = math.sqrt(r[0]**2+r[1]**2+r[2]**2)

    def sep_model(self):

        # 0) inputs

        # adjust distance for AU
        r = self.r / pk.AU

        # spacecraft bus power supply (W)
        Psupply = 700.0

        # 1) solar array performance model

        # performance parameters
        SA = [1.32077, -0.10848, -0.11665, 0.10843, -0.01279]

        # power generated at r
        P = ((self.Psa / r ** 2) * ((SA[0] + (SA[1] / r) + (SA[2] / r ** 2)) / (1 + SA[3] * r + SA[4] * r ** 2)))

        # power available for sep
        P -= Psupply

        # 2) thruster performance model

        # engine parameters

        # operational power range
        P_min = min(self.engine.power)
        P_max = max(self.engine.power)

        # efficiency TODO ====> check mix with PPU efficiency and engine efficiency
        eff = 1.0

        # power available to thruster cluster
        Pin = eff * P

        # power available per thruster
        Pin /= self.no_engine

        # Thruster operational power range
        if Pin > P_max:
            Pin = P_max
        if Pin < P_min:
            Pin = 0

        # thrust coefficients
        Thr_coeffs = self.engine.T_c

        # isp coefficients
        Isp_coeffs = self.engine.Isp_c

        # thrust produced per thruster
        Tmax = sum([K * Pin ** i for K, i in zip(reversed(Thr_coeffs), range(len(Thr_coeffs)))])

        # thrust produced from thuster cluster
        Tmax *= self.no_engine

        if Tmax < 0:
            Tmax = 0

        # isp at r
        # TODO ==== isp currently assumed constant. r changes during spirals, implement fix?
        Isp = sum([J * Pin ** i for J, i in zip(reversed(Isp_coeffs), range(len(Isp_coeffs)))])

        return Tmax, Isp

    def hill_sphere(self):
        rH = self.pl.a * ((self.pl.m / (3 * pd.Sun.m)) ** (1 / 3))
        return rH

    def spiral_dv(self):
        # initial and final radius
        Ri = self.hill_sphere()
        Rf = self.pl.r + self.target_orbit
        # initial and final velocity
        Vi = math.sqrt((pd.G * self.pl.m)/Ri)
        Vf = math.sqrt((pd.G * self.pl.m)/Rf)
        # inclination change
        d_inc = abs(self.i_fin - self.i_init)

        a = Vi ** 2
        b = Vf ** 2
        c = (2 * Vi * Vf) * math.cos((math.pi/2) * d_inc)

        # delta v for orbit change
        d_V = math.sqrt(a + b - c)

        return d_V

    def spiral_prop_time(self):

        # inputs
        mi = self.mi
        T, isp = self.sep_model()
        dv = self.spiral_dv()

        # spiral time of flight (seconds)
        tof = dv / (T / mi)
        #spiral time of flight (days)
        tof = tof / (60 * 60 * 24)

        mf = mi * math.exp(-dv / (isp * 9.81))
        mp = mi - mf

        return tof, mp

class spiral_out():

    def __init__(self,
                 pl = pd.Mars,
                 orbit_alt = 320,
                 i_init = 0,
                 i_fin = 0,
                 Psa = 10000,
                 engine = td.BPT4000,
                 no_engine = 1,
                 epoch = 1000,
                 mi = 4000):

        # 1 - Class Data Members

        # planet and orbit
        self.pl = pl
        self.target_orbit = orbit_alt

        # initial and final inclination
        self.i_init = i_init
        self.i_fin = i_fin

        # spacecraft initial mass
        self.mi = mi

        # solar arrays
        self.Psa = Psa

        # propulsion subsystem
        # thrusters
        self.engine = engine
        self.no_engine = no_engine
        # max performance
        Thr_max = max(engine.thrust)
        Isp_max = max(engine.isp)
        Tmax = Thr_max * no_engine

        # planet
        planet = pk.planet.jpl_lp(self.pl.name)
        # helio-radius and velocity at target epoch
        r, v = planet.eph(pk.epoch(epoch))
        # solar distance
        self.r = math.sqrt(r[0] ** 2 + r[1] ** 2 + r[2] ** 2)

    def sep_model(self):

        # 0) inputs

        # adjust distance for AU
        r = self.r / pk.AU

        # spacecraft bus power supply (W)
        Psupply = 700.0

        # 1) solar array performance model

        # performance parameters
        SA = [1.32077, -0.10848, -0.11665, 0.10843, -0.01279]

        # power generated at r
        P = ((self.Psa / r ** 2) * ((SA[0] + (SA[1] / r) + (SA[2] / r ** 2)) / (1 + SA[3] * r + SA[4] * r ** 2)))

        # power available for sep
        P -= Psupply

        # 2) thruster performance model

        # engine parameters

        # operational power range
        P_min = min(self.engine.power)
        P_max = max(self.engine.power)

        # efficiency TODO ====> check mix with PPU efficiency and engine efficiency
        eff = 1.0

        # power available to thruster cluster
        Pin = eff * P

        # power available per thruster
        Pin /= self.no_engine

        # Thruster operational power range
        if Pin > P_max:
            Pin = P_max
        if Pin < P_min:
            Pin = 0

        # thrust coefficients
        Thr_coeffs = self.engine.T_c

        # isp coefficients
        Isp_coeffs = self.engine.Isp_c

        # thrust produced per thruster
        Tmax = sum([K * Pin ** i for K, i in zip(reversed(Thr_coeffs), range(len(Thr_coeffs)))])

        # thrust produced from thuster cluster
        Tmax *= self.no_engine

        if Tmax < 0:
            Tmax = 0

        # isp at r
        # TODO ==== isp currently assumed constant. r changes during spirals, implement fix?
        Isp = sum([J * Pin ** i for J, i in zip(reversed(Isp_coeffs), range(len(Isp_coeffs)))])

        return Tmax, Isp

    def hill_sphere(self):
        rH = self.pl.a * ((self.pl.m / (3 * pd.Sun.m)) ** (1 / 3))
        return rH

    def spiral_dv(self):
        # initial and final radius
        Ri = self.pl.r + self.target_orbit
        Rf = self.hill_sphere()
        # initial and final velocity
        Vi = math.sqrt((pd.G * self.pl.m)/Ri)
        Vf = math.sqrt((pd.G * self.pl.m)/Rf)
        # inclination change
        d_inc = abs(self.i_fin - self.i_init)

        a = Vi ** 2
        b = Vf ** 2
        c = (2 * Vi * Vf) * math.cos((math.pi/2) * d_inc)

        # delta v for orbit change
        d_V = math.sqrt(a + b - c)

        return d_V

    def spiral_prop_time(self):

        # inputs
        mi = self.mi
        T, isp = self.sep_model()
        dv = self.spiral_dv()

        # spiral time of flight (seconds)
        tof = dv / (T / mi)
        #spiral time of flight (days)
        tof = tof / (60 * 60 * 24)

        mf = mi * math.exp(-dv / (isp * 9.81))
        mp = mi - mf

        return tof, mp





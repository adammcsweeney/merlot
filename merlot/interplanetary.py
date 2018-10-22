# Module containing classes for modelling interplanetary transfers
# Currently considers outbound and return transfers from Earth to Solar System planets

# TODO ==== merge into single class, specifying outbound or inbound leg
# TODO ==== make thruster model compatible with thruster_data module
# TODO ==== include link to solar system planetary data module
# TODO ==== include link to launch vehicle model
# TODO ==== solve --> RuntimeError: Maximum number of iteration reached
# TODO ==== make plots 2D
# TODO ==== include asteroids as targets
# TODO ==== include planet data in writer

import pykep as pk
import numpy as np
import math
import csv
import merlot.thruster_data as td
import merlot.launcher_data as lv
import merlot.physical_data as pd

class merlot_sep_e2pl:
    """""
    This class can be used as a User Defined Problem (UDP) in the pygmo2 software and if successfully solved,
    represents a low-thrust interplanetary trajectory from the Earth to a target planet.
    
    The problem takes as input the year of departure from Earth (year), the heliocentric time of flight (tof), 
    the power of the solar arrays at 1 AU, the EP thruster (engine), the number of thrusters on the spacecraft 
    (no_engine), and the launch vehicle (LV).
    
    The arrival hyperbolic excess velocity can also be set by the user. 
    (Anticipated that this will be used in a future version for modelling hybrid mission concepts).

    The trajectory is modeled using the Sims-Flanagan model.

    The propulsion model can be both nuclear (NEP) or solar (SEP). In the case of SEP, the thrust profile is a function 
    of Solar Array power provided at 1 AU, and distance from the Sun.  
    
    The launch vehicle model includes the launcher performance (mass delivered for a given C3) within the optimisation.

    The decision vector (chromosome) is::
      [t0] + [mf] + [vxi, vyi, vzi] + [vxf, vyf, vzf] + [throttles1] + [throttles2] + ...
      [launch epoch] + [mass delivered to Mars] + [departure velocity] + [arrival velocity] + [engine throttles] + ...
    """""

    def __init__(self,
                 target=pd.Mars,
                 n_seg=20,
                 year=2026,
                 tof=500,
                 Psa=10000,
                 engine=td.NSTAR,
                 no_engine=1,
                 LV=lv.F9,
                 vinf_arr=0.01,
                 sep=True):
        """
        :param target: (``pykep.planet``): target planet
        :param n_seg: (``int``): number of segments to use in the problem transcription (time grid)
        :param t0: (``int``): departure year from earth
        :param tof: (``int``): amount of time allocated for the time of flight (days)
        :param m0: (``float``): initial mass of the spacecraft (kg)
        :param Psa: (``float``): power provided by the solar array at 1 AU and BOL (W)
        :param engine: (``str``): type of ion engine on board the spacecraft
        :param no_engine: (``int``): number of active ion engines during the transfer
        :param LV: (``str``): launch vehicle used for inserting spacecraft in earth escape trajectory
        :param vinf_arr: (``float``): upper bound on target arrival excess velocity (km/s)
        :param sep: (``bool``): activates a Solar Electric Propulsion model for the thrust - distance dependency
        """

        # 1) Checks of input data

        # if target not in ["Mercury", "Venus", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"]:
        #     raise ValueError("Departure planet must be either 'Mercury', 'Venus', 'Mars', 'Jupiter', 'Saturn', "
        #                      "'Uranus', 'Neptune', 'Pluto'")
        # TODO ==== update launcher references and thruster references for data checks
        # if LV not in ["SOYUZ", "FALCON9", "PROTON", "A64", "H-IIA(202)", "DELTA_IV(H)", "ATLAS_V(551)", "ATLAS_V(401)"]:
        #     raise ValueError("Launch vehicle must be either 'SOYUZ', 'FALCON9', 'PROTON', 'A64', 'H-IIA(202)', "
        #                      "'DELTA_IV(H)', 'ATLAS_V(551)', 'ATLAS_V(401)")

        # 2) Class data members

        # planets
        self.target_name = target.name
        self.target = pk.planet.jpl_lp(self.target_name)
        self.earth = pk.planet.jpl_lp('earth')

        # segments
        self.__n_seg = n_seg

        # departure date range
        t_0 = pk.epoch((float(year - 2001) * 365.25) + 366)
        t_f = pk.epoch((float((year + 1) - 2001) * 365.25) + 365.9999)
        date_range = [t_0, t_f]

        # time of flight
        self.tof = tof

        # launch vehicle
        self.LV = LV
        # self.C3_max = max(LV.C3)
        # (below)test to check extended C3 ranges
        self.C3_max = 50
        # TODO ==== update and delete string lv references
        # # launcher parameter selection
        # if self.LV == 'SOYUZ':
        #     self.C3_max = 25.0
        # elif self.LV == 'PROTON':
        #     self.C3_max = 36.0
        # elif self.LV == 'FALCON9':
        #     self.C3_max = 45.0
        # elif self.LV == 'A64':
        #     self.C3_max = 35.0
        # elif self.LV == 'H-IIA(202)':
        #     self.C3_max = 40.0
        # elif self.LV == 'DELTA_IV(H)':
        #     self.C3_max = 100.0
        # elif self.LV == 'ATLAS_V(551)':
        #     self.C3_max = 60.0
        # elif self.LV == 'ATLAS_V(401)':
        #     self.C3_max = 45.0

        # propulsion subsystem
        self.engine = engine
        self.no_engine = no_engine
        Thr_max = max(engine.thrust)
        Isp_max = max(engine.isp)
        Tmax = Thr_max * no_engine

        # power subsystem
        self.Psa = Psa

        # spacecraft
        mass_max = self._lv_model(0)
        # TODO ==== check if this can be moved to fitness method
        self.__sc = pk.sims_flanagan._sims_flanagan.spacecraft(mass_max, Tmax, Isp_max)

        # solar electric propulsion active
        self.sep = sep

        # grid construction
        grid = np.array([i / n_seg for i in range(n_seg + 1)])

        # index corresponding to the middle of the transfer
        fwd_seg = int(np.searchsorted(grid, 0.5, side='right'))
        bwd_seg = n_seg - fwd_seg
        fwd_grid = grid[:fwd_seg + 1]
        bwd_grid = grid[fwd_seg:]
        self.__fwd_seg = fwd_seg
        self.__fwd_grid = fwd_grid
        self.__fwd_dt = np.array([(fwd_grid[i + 1] - fwd_grid[i])
                                  for i in range(fwd_seg)]) * pk.DAY2SEC
        self.__bwd_seg = bwd_seg
        self.__bwd_grid = bwd_grid
        self.__bwd_dt = np.array([(bwd_grid[i + 1] - bwd_grid[i])
                                  for i in range(bwd_seg)]) * pk.DAY2SEC

        # 3) Bounds definition

        # vinf max values
        self.vinf_dep = ((self.C3_max * 1000000) ** 0.5)  # (in m)
        self.vinf_arr = vinf_arr * 1000  # (in m/s)

        # lower and upper bounds
        lb = [date_range[0].mjd2000] + [self.__sc.mass * 0.2] + [-self.vinf_dep] * 3 + [-self.vinf_arr] * 3 \
             + [-1, -1, -1] * n_seg
        ub = [date_range[1].mjd2000] + [self.__sc.mass] + [self.vinf_dep] * 3 + [self.vinf_arr] * 3 + [1, 1, 1] * n_seg
        self.__lb = lb
        self.__ub = ub


    # Fitness function
    def fitness(self, x):
        # fitness function objective - maximise arrival mass (mf)
        obj = -x[1]

        ceq = list()
        cineq = list()
        throttles_constraints = []
        mismatch_constraints = []

        # equality constraints

        rfwd, rbwd, vfwd, vbwd, mfwd, mbwd, _, _, _, _, dfwd, dbwd = self._propagate(x)

        # mismatch constraints
        mismatch_constraints.extend([a - b for a, b in zip(rfwd[-1], rbwd[0])])
        mismatch_constraints.extend([a - b for a, b in zip(vfwd[-1], vbwd[0])])
        mismatch_constraints.append(mfwd[-1] - mbwd[0])
        ceq.extend(mismatch_constraints)

        # Making the mismatches non dimensional
        ceq[0] /= pk.AU
        ceq[1] /= pk.AU
        ceq[2] /= pk.AU
        ceq[3] /= pk.EARTH_VELOCITY
        ceq[4] /= pk.EARTH_VELOCITY
        ceq[5] /= pk.EARTH_VELOCITY
        # TODO ==== more appropriate value to use?
        ceq[6] /= self._lv_model(0.1)

        # inequality constraints

        # throttle constraints
        throttles = [x[8 + 3 * i: 11 + 3 * i] for i in range(self.__n_seg)]
        for t in throttles:
            throttles_constraints.append(t[0] ** 2 + t[1] ** 2 + t[2] ** 2 - 1.)
        cineq.extend(throttles_constraints)

        # vinf constraints
        # departure
        # v_dep_con_1 = (0.0 - (x[2] ** 2 + x[3] ** 2 + x[4] ** 2))
        v_dep_con_2 = ((x[2] ** 2 + x[3] ** 2 + x[4] ** 2) - (self.vinf_dep ** 2))
        # arrival
        v_arr_con = ((x[5] ** 2 + x[6] ** 2 + x[7] ** 2) - (self.vinf_arr ** 2))

        # nondimensionalize vinf inequality constraints
        # v_dep_con_1 /= pk.EARTH_VELOCITY ** 2
        v_dep_con_2 /= pk.EARTH_VELOCITY ** 2
        v_arr_con /= pk.EARTH_VELOCITY ** 2

        # return np.hstack(([obj], ceq, cineq, v_dep_con_1, v_dep_con_2, v_arr_con))

        return np.hstack(([obj], ceq, cineq, v_dep_con_2, v_arr_con))

    # Get bounds
    def get_bounds(self):
        return (self.__lb, self.__ub)

    # Get number of inequality contraints
    def get_nic(self):
        return self.__n_seg + 2

    # Get number of equality contraints
    def get_nec(self):
        return 7

    def _lv_model(self, C3):
        # TODO ==== are the units coming in for C3 correct?
        # TODO ==== find faster way to evaluate performance?
        # launcher parameter selection

        A = self.LV.params
        # # launcher performance model
        # # TODO ==== check if there is a better approach for this
        if C3 > self.C3_max:
            m_i = 0.
        elif C3 < 0.0:
            m_i = 0.
        else:
            m_i = lv.lv_performance(self.LV,C3)


        return m_i

    # SEP model
    def _sep_model(self, r):

        # spacecraft bus power supply (W)
        Psupply = 700.0

        # 1) solar array performance model
        SA = [1.32077, -0.10848, -0.11665, 0.10843, -0.01279]

        P = ((self.Psa / r ** 2) * ((SA[0] + (SA[1] / r) + (SA[2] / r ** 2)) / (1 + SA[3] * r + SA[4] * r ** 2)))

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

        # Thruster power limits
        if Pin > P_max:
            Pin = P_max
        if Pin < P_min:
            Pin = 0

        # thrust coefficients
        Thr_coeffs = self.engine.T_c

        # isp coefficients
        Isp_coeffs = self.engine.Isp_c


        # thrust (N)

        Tmax = sum([K * Pin ** i for K, i in zip(reversed(Thr_coeffs), range(len(Thr_coeffs)))])

        # TODO ==== check if this is appropriate means to modle multiple engines firing
        Tmax *= self.no_engine

        if Tmax < 0:
            Tmax = 0

        # Isp (s)

        Isp = sum([J * Pin ** i for J, i in zip(reversed(Isp_coeffs), range(len(Isp_coeffs)))])

        return Tmax, Isp

    # Propagates the trajectory
    def _propagate(self, x):
        # 1 - Decode the chromosome
        t0 = x[0]
        T = self.tof
        m_f = x[1]

        # Extract the number of segments for forward and backward
        # propagation
        n_seg = self.__n_seg
        fwd_seg = self.__fwd_seg
        bwd_seg = self.__bwd_seg
        # We extract information on the spacecraft
        m_i = self._lv_model(((x[2] / 1000) ** 2) + ((x[3] / 1000) ** 2) + ((x[4] / 1000) ** 2))
        max_thrust = self.__sc.thrust
        isp = self.__sc.isp
        veff = isp * pk.G0
        # And on the leg
        throttles = [x[8 + 3 * i: 11 + 3 * i] for i in range(n_seg)]
        # Return lists
        n_points_fwd = fwd_seg + 1
        n_points_bwd = bwd_seg + 1
        rfwd = [[0.0] * 3] * (n_points_fwd)
        vfwd = [[0.0] * 3] * (n_points_fwd)
        mfwd = [0.0] * (n_points_fwd)
        ufwd = [[0.0] * 3] * (fwd_seg)
        dfwd = [[0.0] * 3] * (fwd_seg)  # distances E/S
        rbwd = [[0.0] * 3] * (n_points_bwd)
        vbwd = [[0.0] * 3] * (n_points_bwd)
        mbwd = [0.0] * (n_points_bwd)
        ubwd = [[0.0] * 3] * (bwd_seg)
        dbwd = [[0.0] * 3] * (bwd_seg)

        # 2 - We compute the initial and final epochs and ephemerides
        t_i = pk.epoch(t0)
        r_i, v_i = self.earth.eph(t_i)
        t_f = pk.epoch(t0 + T)
        r_f, v_f = self.target.eph(t_f)

        # 3 - Forward propagation
        fwd_grid = t0 + T * self.__fwd_grid  # days
        fwd_dt = T * self.__fwd_dt  # seconds
        # Initial conditions
        rfwd[0] = r_i
        # vfwd[0] = v_i
        # TODO ==== addition of launcher vinf to v_i, is this correct?
        vfwd[0] = tuple(sum(i) for i in zip((x[2], x[3], x[4]), v_i))
        mfwd[0] = m_i
        # Propagate
        # TODO ====> add in forced coast constraints t0+30 T_max = 0?
        for i, t in enumerate(throttles[:fwd_seg]):
            if self.sep:
                r = math.sqrt(rfwd[i][0] ** 2 + rfwd[i][1]
                              ** 2 + rfwd[i][2] ** 2) / pk.AU
                max_thrust, isp = self._sep_model(r)
                veff = isp * pk.G0
            ufwd[i] = [max_thrust * thr for thr in t]
            rfwd[i + 1], vfwd[i + 1], mfwd[i + 1] = pk.propagate_taylor(
                rfwd[i], vfwd[i], mfwd[i], ufwd[i], fwd_dt[i], pk.MU_SUN, veff, -10, -10)

        # 4 - Backward propagation
        bwd_grid = t0 + T * self.__bwd_grid  # days
        bwd_dt = T * self.__bwd_dt  # seconds
        # Final conditions
        rbwd[-1] = r_f
        # TODO === is this correct? should arrival vinf be subtracted?
        # vbwd[-1] = v_f
        vbwd[-1] = tuple(sum(i) for i in zip((x[5], x[6], x[7]), v_f))
        mbwd[-1] = m_f
        # Propagate
        for i, t in enumerate(throttles[-1:-bwd_seg - 1:-1]):
            if self.sep:
                r = math.sqrt(rbwd[-i - 1][0] ** 2 + rbwd[-i - 1]
                [1] ** 2 + rbwd[-i - 1][2] ** 2) / pk.AU
                max_thrust, isp = self._sep_model(r)
                veff = isp * pk.G0
            ubwd[-i - 1] = [max_thrust * thr for thr in t]
            rbwd[-i - 2], vbwd[-i - 2], mbwd[-i - 2] = pk.propagate_taylor(
                rbwd[-i - 1], vbwd[-i - 1], mbwd[-i - 1], ubwd[-i - 1], -bwd_dt[-i - 1], pk.MU_SUN, veff, -10, -10)

        return rfwd, rbwd, vfwd, vbwd, mfwd, mbwd, ufwd, ubwd, fwd_dt, bwd_dt, dfwd, dbwd

    # Visualizes the trajectory
    def plot_traj(self, x, units=pk.AU, plot_segments=True, plot_thrusts=False, axes=None):
        """
        ax = prob.plot_traj(self, x, units=AU, plot_segments=True, plot_thrusts=False, axes=None)
        Args:
            - x (``list``, ``tuple``, ``numpy.ndarray``): Decision chromosome, e.g. (``pygmo.population.champion_x``).
            - units (``float``): the length unit to be used in the plot
            - plot_segments (``bool``): when True plots also the segments boundaries
            - plot_thrusts (``bool``): when True plots also the thrust vectors
        Returns:
            matplotlib.axes: axes where to plot
        Visualizes the the trajectory in a 3D plot
        """

        if not len(x) == len(self.get_bounds()[0]):
            raise ValueError("Invalid length of the decision vector x")

        import matplotlib as mpl
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D

        # Creating the axes if necessary
        if axes is None:
            mpl.rcParams['legend.fontsize'] = 10
            fig = plt.figure()
            axes = fig.gca(projection='3d')

        n_seg = self.__n_seg
        fwd_seg = self.__fwd_seg
        bwd_seg = self.__bwd_seg
        t0 = x[0]
        # T = x[1]
        T = self.tof
        isp = self.__sc.isp
        veff = isp * pk.G0
        fwd_grid = t0 + T * self.__fwd_grid  # days
        bwd_grid = t0 + T * self.__bwd_grid  # days

        throttles = [x[8 + 3 * i: 11 + 3 * i] for i in range(n_seg)]
        alphas = [min(1., np.linalg.norm(t)) for t in throttles]

        times = np.concatenate((fwd_grid, bwd_grid))

        rfwd, rbwd, vfwd, vbwd, mfwd, mbwd, ufwd, ubwd, fwd_dt, bwd_dt, dfwd, dbwd = self._propagate(x)

        # Plotting the Sun, the Earth and the target
        axes.scatter([0], [0], [0], color='gold', s=100)
        pk.orbit_plots.plot_planet(self.earth, pk.epoch(
            t0), units=units, legend=True, color=(0.7, 0.7, 1), ax=axes)
        pk.orbit_plots.plot_planet(self.target, pk.epoch(
            t0 + T), units=units, legend=True, color=(0.7, 0.7, 1), ax=axes)

        # set axes range
        # units in AU
        axes.set_xlim(-2, 2)
        axes.set_ylim(-2, 2)
        axes.set_zlim(-2, 2)

        # Forward propagation
        xfwd = [0.0] * (fwd_seg + 1)
        yfwd = [0.0] * (fwd_seg + 1)
        zfwd = [0.0] * (fwd_seg + 1)
        xfwd[0] = rfwd[0][0] / units
        yfwd[0] = rfwd[0][1] / units
        zfwd[0] = rfwd[0][2] / units

        for i in range(fwd_seg):
            if self.sep:
                r = math.sqrt(rfwd[i][0] ** 2 + rfwd[i][1]
                              ** 2 + rfwd[i][2] ** 2) / pk.AU
                _, isp = self._sep_model(r)
                veff = isp * pk.G0
            pk.orbit_plots.plot_taylor(rfwd[i], vfwd[i], mfwd[i], ufwd[i], fwd_dt[
                i], pk.MU_SUN, veff, N=10, units=units, color=(alphas[i], 0, 1 - alphas[i]), ax=axes)
            xfwd[i + 1] = rfwd[i + 1][0] / units
            yfwd[i + 1] = rfwd[i + 1][1] / units
            zfwd[i + 1] = rfwd[i + 1][2] / units
        if plot_segments:
            axes.scatter(xfwd[:-1], yfwd[:-1], zfwd[:-1],
                         label='nodes', marker='o', s=5, c='k')

        # Backward propagation
        xbwd = [0.0] * (bwd_seg + 1)
        ybwd = [0.0] * (bwd_seg + 1)
        zbwd = [0.0] * (bwd_seg + 1)
        xbwd[-1] = rbwd[-1][0] / units
        ybwd[-1] = rbwd[-1][1] / units
        zbwd[-1] = rbwd[-1][2] / units

        for i in range(bwd_seg):
            if self.sep:
                r = math.sqrt(rbwd[-i - 1][0] ** 2 + rbwd[-i - 1]
                [1] ** 2 + rbwd[-i - 1][2] ** 2) / pk.AU
                _, isp = self._sep_model(r)
                veff = isp * pk.G0
            pk.orbit_plots.plot_taylor(rbwd[-i - 1], vbwd[-i - 1], mbwd[-i - 1], ubwd[-i - 1], -bwd_dt[-i - 1],
                                       pk.MU_SUN, veff, N=10, units=units,
                                       color=(alphas[-i - 1], 0, 1 - alphas[-i - 1]), ax=axes)
            xbwd[-i - 2] = rbwd[-i - 2][0] / units
            ybwd[-i - 2] = rbwd[-i - 2][1] / units
            zbwd[-i - 2] = rbwd[-i - 2][2] / units
        if plot_segments:
            axes.scatter(xbwd[1:], ybwd[1:], zbwd[1:], marker='o', s=5, c='k')

        # Plotting the thrust vectors
        if plot_thrusts:
            throttles = np.array(throttles)
            xlim = axes.get_xlim()
            xrange = xlim[1] - xlim[0]
            ylim = axes.get_ylim()
            yrange = ylim[1] - ylim[0]
            zlim = axes.get_zlim()
            zrange = zlim[1] - zlim[0]

            scale = 0.005

            throttles[:, 0] *= xrange
            throttles[:, 1] *= yrange
            throttles[:, 2] *= zrange

            throttles *= scale

            for (x, y, z, t) in zip(xfwd[:-1] + xbwd[:-1], yfwd[:-1] + ybwd[:-1], zfwd[:-1] + zbwd[:-1], throttles):
                axes.plot([x, x + t[0]], [y, y + t[1]], [z, z + t[2]], c='r')

        return axes

    def plot_dists_thrust(self, x, axes=None):
        """
        axes = prob.plot_dists_thrust(self, x, axes=None)
        Args:
            - x (``list``, ``tuple``, ``numpy.ndarray``): Decision chromosome, e.g. (``pygmo.population.champion_x``).
        Returns:
            matplotlib.axes: axes where to plot
        Plots the distance of the spacecraft from the Earth/Sun and the thrust profile
        """

        if not len(x) == len(self.get_bounds()[0]):
            raise ValueError("Invalid length of the decision vector x")

        import matplotlib as mpl
        import matplotlib.pyplot as plt

        # Creating the axes if necessary
        if axes is None:
            mpl.rcParams['legend.fontsize'] = 10
            fig = plt.figure()
            axes = fig.add_subplot(111)

        n_seg = self.__n_seg
        fwd_seg = self.__fwd_seg
        bwd_seg = self.__bwd_seg
        t0 = x[0]
        T = self.tof
        fwd_grid = t0 + T * self.__fwd_grid  # days
        bwd_grid = t0 + T * self.__bwd_grid  # days

        throttles = [np.linalg.norm(x[8 + 3 * i: 11 + 3 * i])
                     for i in range(n_seg)]

        dist_earth = [0.0] * (n_seg + 2)  # distances spacecraft - Earth
        dist_sun = [0.0] * (n_seg + 2)  # distances spacecraft - Sun
        times = np.concatenate((fwd_grid, bwd_grid))

        rfwd, rbwd, vfwd, vbwd, mfwd, mbwd, ufwd, ubwd, fwd_dt, bwd_dt, _, _ = self._propagate(x)

        # Forward propagation
        xfwd = [0.0] * (fwd_seg + 1)
        yfwd = [0.0] * (fwd_seg + 1)
        zfwd = [0.0] * (fwd_seg + 1)
        xfwd[0] = rfwd[0][0] / pk.AU
        yfwd[0] = rfwd[0][1] / pk.AU
        zfwd[0] = rfwd[0][2] / pk.AU
        r_E = [ri / pk.AU for ri in self.earth.eph(pk.epoch(fwd_grid[0]))[0]]
        dist_earth[0] = np.linalg.norm(
            [r_E[0] - xfwd[0], r_E[1] - yfwd[0], r_E[2] - zfwd[0]])
        dist_sun[0] = np.linalg.norm([xfwd[0], yfwd[0], zfwd[0]])

        for i in range(fwd_seg):
            xfwd[i + 1] = rfwd[i + 1][0] / pk.AU
            yfwd[i + 1] = rfwd[i + 1][1] / pk.AU
            zfwd[i + 1] = rfwd[i + 1][2] / pk.AU
            r_E = [
                ri / pk.AU for ri in self.earth.eph(pk.epoch(fwd_grid[i + 1]))[0]]
            dist_earth[
                i + 1] = np.linalg.norm([r_E[0] - xfwd[i + 1], r_E[1] - yfwd[i + 1], r_E[2] - zfwd[i + 1]])
            dist_sun[
                i + 1] = np.linalg.norm([xfwd[i + 1], yfwd[i + 1], zfwd[i + 1]])

        # Backward propagation
        xbwd = [0.0] * (bwd_seg + 1)
        ybwd = [0.0] * (bwd_seg + 1)
        zbwd = [0.0] * (bwd_seg + 1)
        xbwd[-1] = rbwd[-1][0] / pk.AU
        ybwd[-1] = rbwd[-1][1] / pk.AU
        zbwd[-1] = rbwd[-1][2] / pk.AU
        r_E = [
            ri / pk.AU for ri in self.earth.eph(pk.epoch(bwd_grid[-1]))[0]]
        dist_earth[-1] = np.linalg.norm([r_E[0] - xbwd[-1],
                                         r_E[1] - ybwd[-1], r_E[2] - zbwd[-1]])
        dist_sun[-1] = np.linalg.norm([xbwd[-1], ybwd[-1], zbwd[-1]])

        for i in range(bwd_seg):
            xbwd[-i - 2] = rbwd[-i - 2][0] / pk.AU
            ybwd[-i - 2] = rbwd[-i - 2][1] / pk.AU
            zbwd[-i - 2] = rbwd[-i - 2][2] / pk.AU
            r_E = [
                ri / pk.AU for ri in self.earth.eph(pk.epoch(bwd_grid[-i - 2]))[0]]
            dist_earth[-i - 2] = np.linalg.norm(
                [r_E[0] - xbwd[-i - 2], r_E[1] - ybwd[-i - 2], r_E[2] - zbwd[-i - 2]])
            dist_sun[-i -
                     2] = np.linalg.norm([xbwd[-i - 2], ybwd[-i - 2], zbwd[-i - 2]])

        axes.set_xlabel("t [mjd2000]")
        # Plot Earth distance
        axes.plot(times, dist_earth, c='b', label="S/C-Earth")
        # Plot Sun distance
        axes.plot(times, dist_sun, c='gold', label="S/C-Sun")
        axes.set_ylabel("distance [AU]", color='k')
        axes.set_ylim(bottom=0.)
        axes.tick_params('y', colors='k')
        axes.legend(loc=2)

        # Plot thrust profile
        axes = axes.twinx()
        if self.sep:
            max_thrust = self.__sc.thrust
            thrusts = np.linalg.norm(
                np.array(ufwd + ubwd), axis=1) / max_thrust
            # plot maximum thrust achievable at that distance from the Sun
            distsSun = dist_sun[:fwd_seg] + \
                       dist_sun[-bwd_seg:] + [dist_sun[-1]]
            Tmaxs = [self._sep_model(d)[0] / max_thrust for d in distsSun]
            axes.step(np.concatenate(
                (fwd_grid, bwd_grid[1:])), Tmaxs, where="post", c='lightgray', linestyle=':')
        else:
            thrusts = throttles.copy()
        # duplicate the last for plotting
        thrusts = np.append(thrusts, thrusts[-1])
        axes.step(np.concatenate(
            (fwd_grid, bwd_grid[1:])), thrusts, where="post", c='k', linestyle='--')
        axes.set_ylabel("T/Tmax$_{1AU}$", color='k')
        axes.tick_params('y', colors='k')
        axes.set_xlim([times[0], times[-1]])
        axes.set_ylim([0, max(thrusts) + 0.2])

        print(f"{'max dist'}\t{max(dist_sun)}")
        print(f"{'min dist'}\t{min(dist_sun)}")

        return axes

    def pretty(self, x):
        """
        prob.pretty(x)
        Args:
            - x (``list``, ``tuple``, ``numpy.ndarray``): Decision chromosome, e.g. (``pygmo.population.champion_x``).
        Prints human readable information on the trajectory represented by the decision vector x
        """

        if not len(x) == len(self.get_bounds()[0]):
            raise ValueError("Invalid length of the decision vector x")

        n_seg = self.__n_seg
        self.m_i = self._lv_model((x[2] / 1000) ** 2 + (x[3] / 1000) ** 2 + (x[4] / 1000) ** 2)
        t0 = x[0]
        T = self.tof
        self.m_f = x[1]
        thrusts = [np.linalg.norm(x[8 + 3 * i: 11 + 3 * i])
                   for i in range(n_seg)]

        tf = t0 + T
        self.mP = self.m_i - self.m_f
        self.deltaV = self.__sc.isp * pk.G0 * np.log(self.m_i / self.m_f)

        dt = np.append(self.__fwd_dt, self.__bwd_dt) * T / pk.DAY2SEC
        time_thrusts_on = sum(dt[i] for i in range(
            len(thrusts)) if thrusts[i] > 0.1)

        print(f"{'Departure: '}\t{pk.epoch(t0)}")
        print(f"{'Arrival: '}\t{pk.epoch(tf)}")
        print(f"{'Time of flight: '}\t{T}")
        print(f"{'Delta-v: '}\t{self.deltaV}")
        print(f"{'Spacecraft Initial Mass: '}\t{self.m_i}")
        print(f"{'Spacecraft Final Mass: '}\t{self.m_f}")
        print(f"{'Propellant Consumed: '}\t{self.mP}")
        print(f"{'Engine: '}\t{self.engine.name}")
        print(f"{'Number of active engines: '}\t{self.no_engine}")

        print(f"{'Solar array power at 1 AU BOL (W): '}\t{self.Psa}")
        print(f"{'Thrust-on time: '}\t{time_thrusts_on}")

        print(f"{'Launch Vehicle: '}\t{self.LV.name}")
        print(f"{'Launcher C3 (km^2/s^2): '}\t{(x[2] / 1000) ** 2 + (x[3] / 1000) ** 2 + (x[4] / 1000) ** 2}")
        print(
            f"{'Departure excess velocity (km/s): '}\t{((x[2] / 1000) ** 2 + (x[3] / 1000) ** 2 + (x[4] / 1000) ** 2) ** 0.5}")

        print(f"{'Arrival velocity: '}\t{((x[5]) ** 2 + (x[6]) ** 2 + (x[7]) ** 2) ** 0.5}\t{x[5]}\t{x[6]}\t{x[7]}")

    def pretty_f(self, f):
        print("=====Tolerance Check=====")
        print("=====mismatch constraints=====")
        print("mf (kg): ", -f[0])
        print("x mismatch (m): ", f[1] * pk.AU)
        print("x mismatch: ", f[1])
        print("y mismatch (m): ", f[2] * pk.AU)
        print("y mismatch: ", f[2])
        print("z mismatch (m): ", f[3] * pk.AU)
        print("z mismatch: ", f[3])
        print("vx mismatch (m/s): ", f[4] * pk.EARTH_VELOCITY)
        print("vx mismatch: ", f[4])
        print("vy mismatch (m/s): ", f[5] * pk.EARTH_VELOCITY)
        print("vy mismatch: ", f[5])
        print("vz mismatch (m/s): ", f[6] * pk.EARTH_VELOCITY)
        print("vz mismatch: ", f[6])
        print("mass mismatch (kg): ", f[7] * self._lv_model(0.1))
        print("mass mismatch: ", f[7])

    def data_reporting(self, x, runtime, sp_data, ID):
        # current time
        import datetime
        now = datetime.datetime.now()
        # number of segments
        n_seg = self.__n_seg
        # launch mass
        self.m_i = self._lv_model((x[2] / 1000) ** 2 + (x[3] / 1000) ** 2 + (x[4] / 1000) ** 2)
        # launch epoch
        t0 = x[0]
        # helio time of flight
        T = self.tof
        # mass at target arrival
        self.m_f = x[1]
        # target arrival epoch
        tf = t0 + T
        # helio propellant mass
        self.mP = self.m_i - self.m_f
        # effective helio delta-v
        self.deltaV = self.__sc.isp * pk.G0 * np.log(self.m_i / self.m_f)
        # mass delivered to target orbit
        delivered_mass = x[1] - sp_data[5]

        # data collection for max/min helio distance
        fwd_seg = self.__fwd_seg
        bwd_seg = self.__bwd_seg
        dist_sun = [0.0] * (n_seg + 2)  # distances spacecraft - Sun
        rfwd, rbwd, vfwd, vbwd, mfwd, mbwd, ufwd, ubwd, fwd_dt, bwd_dt, _, _ = self._propagate(x)
        # Forward propagation
        xfwd = [0.0] * (fwd_seg + 1)
        yfwd = [0.0] * (fwd_seg + 1)
        zfwd = [0.0] * (fwd_seg + 1)
        xfwd[0] = rfwd[0][0] / pk.AU
        yfwd[0] = rfwd[0][1] / pk.AU
        zfwd[0] = rfwd[0][2] / pk.AU
        dist_sun[0] = np.linalg.norm([xfwd[0], yfwd[0], zfwd[0]])
        for i in range(fwd_seg):
            xfwd[i + 1] = rfwd[i + 1][0] / pk.AU
            yfwd[i + 1] = rfwd[i + 1][1] / pk.AU
            zfwd[i + 1] = rfwd[i + 1][2] / pk.AU
            dist_sun[
                i + 1] = np.linalg.norm([xfwd[i + 1], yfwd[i + 1], zfwd[i + 1]])
        # Backward propagation
        xbwd = [0.0] * (bwd_seg + 1)
        ybwd = [0.0] * (bwd_seg + 1)
        zbwd = [0.0] * (bwd_seg + 1)
        xbwd[-1] = rbwd[-1][0] / pk.AU
        ybwd[-1] = rbwd[-1][1] / pk.AU
        zbwd[-1] = rbwd[-1][2] / pk.AU
        dist_sun[-1] = np.linalg.norm([xbwd[-1], ybwd[-1], zbwd[-1]])
        for i in range(bwd_seg):
            xbwd[-i - 2] = rbwd[-i - 2][0] / pk.AU
            ybwd[-i - 2] = rbwd[-i - 2][1] / pk.AU
            zbwd[-i - 2] = rbwd[-i - 2][2] / pk.AU
            dist_sun[-i -
                     2] = np.linalg.norm([xbwd[-i - 2], ybwd[-i - 2], zbwd[-i - 2]])

        with open('merlot_outbound_data.csv', 'a+', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter=',')
            writer.writerow(  # RUN DATA
                [ID, str(now), n_seg, runtime,
                 # INPUTS
                 'Earth', self.target_name, self.Psa, self.engine.name, self.no_engine,self.LV.name,
                 # TIMELINE DATA
                 pk.epoch(t0), T, pk.epoch(t0+T), sp_data[2], sp_data[4],
                 # DEPARTURE AND ARRIVAL EXCESS VELOCITY
                  (x[2] / 1000) ** 2 + (x[3] / 1000) ** 2 + (x[4] / 1000) ** 2,
                 ((x[2] / 1000) ** 2 + (x[3] / 1000) ** 2 + (x[4] / 1000) ** 2) ** 0.5,
                 ((x[5]) ** 2 + (x[6]) ** 2 + (x[7]) ** 2) ** 0.5,
                 # MASS DATA
                 self.m_i, self.mP, self.m_f, sp_data[5], delivered_mass,
                 # DELTA-V
                 self.deltaV, sp_data[3],
                 # SPIRAL THRUST PARAMETERS
                 sp_data[0], sp_data[1],
                 # HELIO DISTANCE
                 max(dist_sun), min(dist_sun)])

class merlot_sep_pl2e:
    """""
    This class can be used as a User Defined Problem (UDP) in the pygmo2 software and if successfully solved,
    represents a low-thrust interplanetary trajectory from a target planet to Earth. 
    
    The problem takes as input the year of departure from target planet (year), the heliocentric time of flight (tof), 
    the power of the solar arrays at 1 AU, the EP thruster (engine), the number of thrusters on the spacecraft 
    (no_engine), and the mass returned to Earth.
    
    The departure and arrival hyperbolic excess velocity can also be set by the user. 
    (Anticipated that these will be used in a future version for modelling hybrid mission concepts).

    The trajectory is modeled using the Sims-Flanagan model.

    The propulsion model can be both nuclear (NEP) or solar (SEP). In the case of SEP, the thrust profile is a function 
    of Solar Array power provided at 1 AU, and distance from the Sun.
    
    The 

    The decision vector (chromosome) is::
      [t0] + [mp] + [vxi, vyi, vzi] + [vxf, vyf, vzf] + [throttles1] + [throttles2] + ...
      [launch epoch] + [propellant mass] + [departure velocity] + [arrival velocity] + [engine throttles] + ...
    """""


    def __init__(self,
                 target=pd.Mars,
                 n_seg=20,
                 year=2028,
                 tof=500,
                 Psa=10000,
                 engine=td.NSTAR,
                 no_engine=1,
                 vinf_dep = 1e-5,
                 vinf_arr=6.8,
                 return_mass = 1000.0,
                 sep=True):
        """
        :param target: (``pykep.planet``): target planet
        :param n_seg: (``int``): number of segments to use in the problem transcription (time grid)
        :param t0: (``int``): departure year from earth
        :param tof: (``int``): amount of time allocated for the time of flight (days)
        :param m0: (``float``): initial mass of the spacecraft (kg)
        :param Psa: (``float``): power provided by the solar array at 1 AU and BOL (W)
        :param engine: (``str``): type of ion engine on board the spacecraft
        :param no_engine: (``int``): number of active ion engines during the transfer
        :param LV: (``str``): launch vehicle used for inserting spacecraft in earth escape trajectory
        :param vinf_arr: (``float``): upper bound on target arrival excess velocity (km/s)
        :param sep: (``bool``): activates a Solar Electric Propulsion model for the thrust - distance dependency
        """

        # 1) Checks of input data

        # if target not in ["Mercury", "Venus", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"]:
        #     raise ValueError("Departure planet must be either 'Mercury', 'Venus', 'Mars', 'Jupiter', 'Saturn', "
        #                      "'Uranus', 'Neptune', 'Pluto'")
        # TODO ==== DELETE LV REFERENCES
        # if LV not in ["SOYUZ", "FALCON9", "PROTON", "A64", "H-IIA(202)", "DELTA_IV(H)", "ATLAS_V(551)", "ATLAS_V(401)"]:
        #     raise ValueError("Launch vehicle must be either 'SOYUZ', 'FALCON9', 'PROTON', 'A64', 'H-IIA(202)', "
        #                      "'DELTA_IV(H)', 'ATLAS_V(551)', 'ATLAS_V(401)")
        # TODO ==== fix engine checks with thruster_data module
        # if engine not in ["T6", "RIT2X", "HT20k"]:
        #     raise ValueError("Engine must be either 'T6', 'RIT2X', 'HT20k'")

        # 2) Class data members

        # planets
        self.target_name = target
        self.target = pk.planet.jpl_lp(self.target_name)
        self.earth = pk.planet.jpl_lp('earth')

        # segments
        self.__n_seg = n_seg

        # departure date range
        t_0 = pk.epoch((float(year - 2001) * 365.25) + 366)
        t_f = pk.epoch((float((year + 1) - 2001) * 365.25) + 365.9999)
        date_range = [t_0, t_f]

        # time of flight
        self.tof = tof

        # TODO ==== DELETE LV REFERENCES
        # # launch vehicle

        # propulsion subsystem
        self.engine = engine
        self.no_engine = no_engine
        Thr_max = max(engine.thrust)
        Isp_max = max(engine.isp)
        # if engine == td.T6:
        #     Thr_max = max(engine.thrust)
        #     Isp_max = max(engine.isp)
        # elif engine == td.NSTAR:
        #     Thr_max = max(engine.thrust)
        #     Isp_max = max(engine.isp)
        # elif engine == td.HT20k:
        #     Thr_max = max(engine.thrust)
        #     Isp_max = max(engine.isp)
        Tmax = Thr_max * no_engine

        # power subsystem
        self.Psa = Psa

        # spacecraft
        # mass_max = self._lv_model(0)
        # TODO ==== check if changing the number below impacts results. Check validity of this value (currently, 1.3)
        mass_max = return_mass*1.3
        # TODO ==== check if this can be moved to fitness method
        self.__sc = pk.sims_flanagan._sims_flanagan.spacecraft(mass_max, Tmax, Isp_max)

        # solar electric propulsion active
        self.sep = sep

        # grid construction
        grid = np.array([i / n_seg for i in range(n_seg + 1)])

        # index corresponding to the middle of the transfer
        fwd_seg = int(np.searchsorted(grid, 0.5, side='right'))
        bwd_seg = n_seg - fwd_seg
        fwd_grid = grid[:fwd_seg + 1]
        bwd_grid = grid[fwd_seg:]
        self.__fwd_seg = fwd_seg
        self.__fwd_grid = fwd_grid
        self.__fwd_dt = np.array([(fwd_grid[i + 1] - fwd_grid[i])
                                  for i in range(fwd_seg)]) * pk.DAY2SEC
        self.__bwd_seg = bwd_seg
        self.__bwd_grid = bwd_grid
        self.__bwd_dt = np.array([(bwd_grid[i + 1] - bwd_grid[i])
                                  for i in range(bwd_seg)]) * pk.DAY2SEC

        # MASS RETURNED TO EARTH

        self.return_mass = return_mass

        # 3) Bounds definition

        # vinf max values
        # self.vinf_dep = ((self.C3_max * 1000000) ** 0.5)  # (in m)
        self.vinf_dep = vinf_dep * 1000  # (in m/s)
        self.vinf_arr = vinf_arr * 1000  # (in m/s)

        # lower and upper bounds
        lb = [date_range[0].mjd2000] + [0] + [-self.vinf_dep] * 3 + [-self.vinf_arr] * 3 + [-1, -1, -1] * n_seg
        ub = [date_range[1].mjd2000] + [self.__sc.mass-self.return_mass] + [self.vinf_dep] * 3 + [self.vinf_arr] * 3 \
             + [1, 1, 1] * n_seg
        self.__lb = lb
        self.__ub = ub

    # Fitness function
    def fitness(self, x):

        # fitness function objective - MINIMISE PROPELLANT MASS (mp)
        obj = x[1]

        ceq = list()
        cineq = list()
        throttles_constraints = []
        mismatch_constraints = []

        # equality constraints

        rfwd, rbwd, vfwd, vbwd, mfwd, mbwd, _, _, _, _, dfwd, dbwd = self._propagate(x)

        # mismatch constraints
        mismatch_constraints.extend([a - b for a, b in zip(rfwd[-1], rbwd[0])])
        mismatch_constraints.extend([a - b for a, b in zip(vfwd[-1], vbwd[0])])
        mismatch_constraints.append(mfwd[-1] - mbwd[0])
        ceq.extend(mismatch_constraints)

        # Making the mismatches non dimensional
        ceq[0] /= pk.AU
        ceq[1] /= pk.AU
        ceq[2] /= pk.AU
        ceq[3] /= pk.EARTH_VELOCITY
        ceq[4] /= pk.EARTH_VELOCITY
        ceq[5] /= pk.EARTH_VELOCITY
        # TODO ==== more appropriate value to use?
        # ceq[6] /= self._lv_model(0.1)
        ceq[6] /= x[1]+self.return_mass

        # inequality constraints

        # throttle constraints
        throttles = [x[8 + 3 * i: 11 + 3 * i] for i in range(self.__n_seg)]
        for t in throttles:
            throttles_constraints.append(t[0] ** 2 + t[1] ** 2 + t[2] ** 2 - 1.)
        cineq.extend(throttles_constraints)

        # vinf constraints
        # departure
        # v_dep_con_1 = (0.0 - (x[2] ** 2 + x[3] ** 2 + x[4] ** 2))
        v_dep_con_2 = ((x[2] ** 2 + x[3] ** 2 + x[4] ** 2) - (self.vinf_dep ** 2))
        # arrival
        v_arr_con = ((x[5] ** 2 + x[6] ** 2 + x[7] ** 2) - (self.vinf_arr ** 2))

        # nondimensionalize vinf inequality constraints
        # v_dep_con_1 /= pk.EARTH_VELOCITY ** 2
        v_dep_con_2 /= pk.EARTH_VELOCITY ** 2
        v_arr_con /= pk.EARTH_VELOCITY ** 2

        # return np.hstack(([obj], ceq, cineq, v_dep_con_1, v_dep_con_2, v_arr_con))

        return np.hstack(([obj], ceq, cineq, v_dep_con_2, v_arr_con))

        # Get bounds

    def get_bounds(self):
        return (self.__lb, self.__ub)
        # Get number of inequality contraints

    def get_nic(self):
        return self.__n_seg + 2
        # Get number of equality contraints

    def get_nec(self):
        return 7

    # SEP model
    def _sep_model(self, r):

        # spacecraft bus power supply (W)
        Psupply = 700.0

        # 1) solar array performance model
        SA = [1.32077, -0.10848, -0.11665, 0.10843, -0.01279]

        P = ((self.Psa / r ** 2) * ((SA[0] + (SA[1] / r) + (SA[2] / r ** 2)) / (1 + SA[3] * r + SA[4] * r ** 2)))

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

        # Thruster power limits
        # TODO introduce engine switching loop. i.e. if insufficient power, go to two engines?
        if Pin > P_max:
            Pin = P_max
        if Pin < P_min:
            Pin = 0

        # thrust coefficients
        Thr_coeffs = self.engine.T_c

        # isp coefficients
        Isp_coeffs = self.engine.Isp_c

        # thrust (N)

        Tmax = sum([K * Pin ** i for K, i in zip(reversed(Thr_coeffs), range(len(Thr_coeffs)))])

        # TODO ==== check if this is appropriate means to modle multiple engines firing
        Tmax *= self.no_engine

        if Tmax < 0:
            Tmax = 0

        # Isp (s)

        Isp = sum([J * Pin ** i for J, i in zip(reversed(Isp_coeffs), range(len(Isp_coeffs)))])

        return Tmax, Isp

    # Propagates the trajectory
    def _propagate(self, x):

        # 1 - Decode the chromosome
        t0 = x[0]
        T = self.tof
        # TODO ==== Change to fix arrival mass
        # m_f = x[1]
        m_f = self.return_mass
        # TODO ==== inclusion of propellant mass as objective
        m_p = x[1]


        # Extract the number of segments for forward and backward
        # propagation
        n_seg = self.__n_seg
        fwd_seg = self.__fwd_seg
        bwd_seg = self.__bwd_seg
        # We extract information on the spacecraft
        # TODO ==== initial mass found from final mass plus objective (propellant mass)
        # m_i = self._lv_model(((x[2] / 1000) ** 2) + ((x[3] / 1000) ** 2) + ((x[4] / 1000) ** 2))
        m_i = m_f + m_p
        max_thrust = self.__sc.thrust
        isp = self.__sc.isp
        veff = isp * pk.G0
        # And on the leg
        throttles = [x[8 + 3 * i: 11 + 3 * i] for i in range(n_seg)]
        # Return lists
        n_points_fwd = fwd_seg + 1
        n_points_bwd = bwd_seg + 1
        rfwd = [[0.0] * 3] * (n_points_fwd)
        vfwd = [[0.0] * 3] * (n_points_fwd)
        mfwd = [0.0] * (n_points_fwd)
        ufwd = [[0.0] * 3] * (fwd_seg)
        dfwd = [[0.0] * 3] * (fwd_seg)  # distances E/S
        rbwd = [[0.0] * 3] * (n_points_bwd)
        vbwd = [[0.0] * 3] * (n_points_bwd)
        mbwd = [0.0] * (n_points_bwd)
        ubwd = [[0.0] * 3] * (bwd_seg)
        dbwd = [[0.0] * 3] * (bwd_seg)

        # 2 - We compute the initial and final epochs and ephemerides
        t_i = pk.epoch(t0)
        # r_i, v_i = self.earth.eph(t_i)
        r_i, v_i = self.target.eph(t_i)
        t_f = pk.epoch(t0 + T)
        # r_f, v_f = self.target.eph(t_f)
        r_f, v_f = self.earth.eph(t_f)

        # 3 - Forward propagation
        fwd_grid = t0 + T * self.__fwd_grid  # days
        fwd_dt = T * self.__fwd_dt  # seconds

        # Initial conditions
        rfwd[0] = r_i
        # vfwd[0] = v_i
        # TODO ==== addition of launcher vinf to v_i, is this correct?
        vfwd[0] = tuple(sum(i) for i in zip((x[2], x[3], x[4]), v_i))
        mfwd[0] = m_i

        # Propagate
        # TODO ====> add in forced coast constraints t0+30 T_max = 0?
        for i, t in enumerate(throttles[:fwd_seg]):
            if self.sep:
                r = math.sqrt(rfwd[i][0] ** 2 + rfwd[i][1]
                              ** 2 + rfwd[i][2] ** 2) / pk.AU
                max_thrust, isp = self._sep_model(r)
                veff = isp * pk.G0
            ufwd[i] = [max_thrust * thr for thr in t]
            rfwd[i + 1], vfwd[i + 1], mfwd[i + 1] = pk.propagate_taylor(
                rfwd[i], vfwd[i], mfwd[i], ufwd[i], fwd_dt[i], pk.MU_SUN, veff, -10, -10)

        # 4 - Backward propagation
        bwd_grid = t0 + T * self.__bwd_grid  # days
        bwd_dt = T * self.__bwd_dt  # seconds
        # Final conditions
        rbwd[-1] = r_f
        # TODO === is this correct? should arrival vinf be subtracted?
        # vbwd[-1] = v_f
        vbwd[-1] = tuple(sum(i) for i in zip((x[5], x[6], x[7]), v_f))
        mbwd[-1] = self.return_mass
        # Propagate
        for i, t in enumerate(throttles[-1:-bwd_seg - 1:-1]):
            if self.sep:
                r = math.sqrt(rbwd[-i - 1][0] ** 2 + rbwd[-i - 1]
                [1] ** 2 + rbwd[-i - 1][2] ** 2) / pk.AU
                max_thrust, isp = self._sep_model(r)
                veff = isp * pk.G0
            ubwd[-i - 1] = [max_thrust * thr for thr in t]
            rbwd[-i - 2], vbwd[-i - 2], mbwd[-i - 2] = pk.propagate_taylor(
                rbwd[-i - 1], vbwd[-i - 1], mbwd[-i - 1], ubwd[-i - 1], -bwd_dt[-i - 1], pk.MU_SUN, veff, -10, -10)

        return rfwd, rbwd, vfwd, vbwd, mfwd, mbwd, ufwd, ubwd, fwd_dt, bwd_dt, dfwd, dbwd

# Visualizes the trajectory
    def plot_traj(self, x, units=pk.AU, plot_segments=True, plot_thrusts=False, axes=None):
        """
        ax = prob.plot_traj(self, x, units=AU, plot_segments=True, plot_thrusts=False, axes=None)
        Args:
            - x (``list``, ``tuple``, ``numpy.ndarray``): Decision chromosome, e.g. (``pygmo.population.champion_x``).
            - units (``float``): the length unit to be used in the plot
            - plot_segments (``bool``): when True plots also the segments boundaries
            - plot_thrusts (``bool``): when True plots also the thrust vectors
        Returns:
            matplotlib.axes: axes where to plot
        Visualizes the the trajectory in a 3D plot
        """

        if not len(x) == len(self.get_bounds()[0]):
            raise ValueError("Invalid length of the decision vector x")

        import matplotlib as mpl
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D

        # Creating the axes if necessary
        if axes is None:
            mpl.rcParams['legend.fontsize'] = 10
            fig = plt.figure()
            axes = fig.gca(projection='3d')

        n_seg = self.__n_seg
        fwd_seg = self.__fwd_seg
        bwd_seg = self.__bwd_seg
        t0 = x[0]
        # T = x[1]
        T = self.tof
        isp = self.__sc.isp
        veff = isp * pk.G0
        fwd_grid = t0 + T * self.__fwd_grid  # days
        bwd_grid = t0 + T * self.__bwd_grid  # days

        throttles = [x[8 + 3 * i: 11 + 3 * i] for i in range(n_seg)]
        alphas = [min(1., np.linalg.norm(t)) for t in throttles]

        times = np.concatenate((fwd_grid, bwd_grid))

        rfwd, rbwd, vfwd, vbwd, mfwd, mbwd, ufwd, ubwd, fwd_dt, bwd_dt, dfwd, dbwd = self._propagate(x)

        # Plotting the Sun, the Earth and the target
        axes.scatter([0], [0], [0], color='gold', s=100)
        # pk.orbit_plots.plot_planet(self.earth, pk.epoch(
        #     t0), units=units, legend=True, color=(0.7, 0.7, 1), ax=axes)
        # pk.orbit_plots.plot_planet(self.target, pk.epoch(
        #     t0 + T), units=units, legend=True, color=(0.7, 0.7, 1), ax=axes)

        pk.orbit_plots.plot_planet(self.target, pk.epoch(
            t0), units=units, legend=True, color=(0.7, 0.7, 1), ax=axes)
        pk.orbit_plots.plot_planet(self.earth, pk.epoch(
            t0 + T), units=units, legend=True, color=(0.7, 0.7, 1), ax=axes)

        # set axes range
        # units in AU
        axes.set_xlim(-2, 2)
        axes.set_ylim(-2, 2)
        axes.set_zlim(-2, 2)

        # Forward propagation
        xfwd = [0.0] * (fwd_seg + 1)
        yfwd = [0.0] * (fwd_seg + 1)
        zfwd = [0.0] * (fwd_seg + 1)
        xfwd[0] = rfwd[0][0] / units
        yfwd[0] = rfwd[0][1] / units
        zfwd[0] = rfwd[0][2] / units

        for i in range(fwd_seg):
            if self.sep:
                r = math.sqrt(rfwd[i][0] ** 2 + rfwd[i][1]
                              ** 2 + rfwd[i][2] ** 2) / pk.AU
                _, isp = self._sep_model(r)
                veff = isp * pk.G0
            pk.orbit_plots.plot_taylor(rfwd[i], vfwd[i], mfwd[i], ufwd[i], fwd_dt[
                i], pk.MU_SUN, veff, N=10, units=units, color=(alphas[i], 0, 1 - alphas[i]), ax=axes)
            xfwd[i + 1] = rfwd[i + 1][0] / units
            yfwd[i + 1] = rfwd[i + 1][1] / units
            zfwd[i + 1] = rfwd[i + 1][2] / units
        if plot_segments:
            axes.scatter(xfwd[:-1], yfwd[:-1], zfwd[:-1],
                         label='nodes', marker='o', s=5, c='k')

        # Backward propagation
        xbwd = [0.0] * (bwd_seg + 1)
        ybwd = [0.0] * (bwd_seg + 1)
        zbwd = [0.0] * (bwd_seg + 1)
        xbwd[-1] = rbwd[-1][0] / units
        ybwd[-1] = rbwd[-1][1] / units
        zbwd[-1] = rbwd[-1][2] / units

        for i in range(bwd_seg):
            if self.sep:
                r = math.sqrt(rbwd[-i - 1][0] ** 2 + rbwd[-i - 1]
                [1] ** 2 + rbwd[-i - 1][2] ** 2) / pk.AU
                _, isp = self._sep_model(r)
                veff = isp * pk.G0
            pk.orbit_plots.plot_taylor(rbwd[-i - 1], vbwd[-i - 1], mbwd[-i - 1], ubwd[-i - 1], -bwd_dt[-i - 1],
                                       pk.MU_SUN, veff, N=10, units=units,
                                       color=(alphas[-i - 1], 0, 1 - alphas[-i - 1]), ax=axes)
            xbwd[-i - 2] = rbwd[-i - 2][0] / units
            ybwd[-i - 2] = rbwd[-i - 2][1] / units
            zbwd[-i - 2] = rbwd[-i - 2][2] / units
        if plot_segments:
            axes.scatter(xbwd[1:], ybwd[1:], zbwd[1:], marker='o', s=5, c='k')

        # Plotting the thrust vectors
        if plot_thrusts:
            throttles = np.array(throttles)
            xlim = axes.get_xlim()
            xrange = xlim[1] - xlim[0]
            ylim = axes.get_ylim()
            yrange = ylim[1] - ylim[0]
            zlim = axes.get_zlim()
            zrange = zlim[1] - zlim[0]

            scale = 0.005

            throttles[:, 0] *= xrange
            throttles[:, 1] *= yrange
            throttles[:, 2] *= zrange

            throttles *= scale

            for (x, y, z, t) in zip(xfwd[:-1] + xbwd[:-1], yfwd[:-1] + ybwd[:-1], zfwd[:-1] + zbwd[:-1], throttles):
                axes.plot([x, x + t[0]], [y, y + t[1]], [z, z + t[2]], c='r')

        return axes
# Plots distance thrust profile
    def plot_dists_thrust(self, x, axes=None):
        """
        axes = prob.plot_dists_thrust(self, x, axes=None)
        Args:
            - x (``list``, ``tuple``, ``numpy.ndarray``): Decision chromosome, e.g. (``pygmo.population.champion_x``).
        Returns:
            matplotlib.axes: axes where to plot
        Plots the distance of the spacecraft from the Earth/Sun and the thrust profile
        """

        if not len(x) == len(self.get_bounds()[0]):
            raise ValueError("Invalid length of the decision vector x")

        import matplotlib as mpl
        import matplotlib.pyplot as plt

        # Creating the axes if necessary
        if axes is None:
            mpl.rcParams['legend.fontsize'] = 10
            fig = plt.figure()
            axes = fig.add_subplot(111)

        n_seg = self.__n_seg
        fwd_seg = self.__fwd_seg
        bwd_seg = self.__bwd_seg
        t0 = x[0]
        T = self.tof
        fwd_grid = t0 + T * self.__fwd_grid  # days
        bwd_grid = t0 + T * self.__bwd_grid  # days

        throttles = [np.linalg.norm(x[8 + 3 * i: 11 + 3 * i])
                     for i in range(n_seg)]

        dist_earth = [0.0] * (n_seg + 2)  # distances spacecraft - Earth
        dist_sun = [0.0] * (n_seg + 2)  # distances spacecraft - Sun
        times = np.concatenate((fwd_grid, bwd_grid))

        rfwd, rbwd, vfwd, vbwd, mfwd, mbwd, ufwd, ubwd, fwd_dt, bwd_dt, _, _ = self._propagate(x)

        # Forward propagation
        xfwd = [0.0] * (fwd_seg + 1)
        yfwd = [0.0] * (fwd_seg + 1)
        zfwd = [0.0] * (fwd_seg + 1)
        xfwd[0] = rfwd[0][0] / pk.AU
        yfwd[0] = rfwd[0][1] / pk.AU
        zfwd[0] = rfwd[0][2] / pk.AU
        # TODO ==== rename r_E
        r_E = [ri / pk.AU for ri in self.earth.eph(pk.epoch(fwd_grid[0]))[0]]
        # r_E = [ri / pk.AU for ri in self.target.eph(pk.epoch(fwd_grid[0]))[0]]
        dist_earth[0] = np.linalg.norm(
            [r_E[0] - xfwd[0], r_E[1] - yfwd[0], r_E[2] - zfwd[0]])
        dist_sun[0] = np.linalg.norm([xfwd[0], yfwd[0], zfwd[0]])

        for i in range(fwd_seg):
            xfwd[i + 1] = rfwd[i + 1][0] / pk.AU
            yfwd[i + 1] = rfwd[i + 1][1] / pk.AU
            zfwd[i + 1] = rfwd[i + 1][2] / pk.AU
            r_E = [
                ri / pk.AU for ri in self.earth.eph(pk.epoch(fwd_grid[i + 1]))[0]]
            # TODO ==== rename r_E
            # r_E = [
            #     ri / pk.AU for ri in self.target.eph(pk.epoch(fwd_grid[i + 1]))[0]]
            dist_earth[
                i + 1] = np.linalg.norm([r_E[0] - xfwd[i + 1], r_E[1] - yfwd[i + 1], r_E[2] - zfwd[i + 1]])
            dist_sun[
                i + 1] = np.linalg.norm([xfwd[i + 1], yfwd[i + 1], zfwd[i + 1]])

        # Backward propagation
        xbwd = [0.0] * (bwd_seg + 1)
        ybwd = [0.0] * (bwd_seg + 1)
        zbwd = [0.0] * (bwd_seg + 1)
        xbwd[-1] = rbwd[-1][0] / pk.AU
        ybwd[-1] = rbwd[-1][1] / pk.AU
        zbwd[-1] = rbwd[-1][2] / pk.AU
        r_E = [
            ri / pk.AU for ri in self.earth.eph(pk.epoch(bwd_grid[-1]))[0]]
        # TODO ==== rename r_E
        # r_E = [
        #     ri / pk.AU for ri in self.target.eph(pk.epoch(bwd_grid[-1]))[0]]
        dist_earth[-1] = np.linalg.norm([r_E[0] - xbwd[-1],
                                         r_E[1] - ybwd[-1], r_E[2] - zbwd[-1]])
        dist_sun[-1] = np.linalg.norm([xbwd[-1], ybwd[-1], zbwd[-1]])

        for i in range(bwd_seg):
            xbwd[-i - 2] = rbwd[-i - 2][0] / pk.AU
            ybwd[-i - 2] = rbwd[-i - 2][1] / pk.AU
            zbwd[-i - 2] = rbwd[-i - 2][2] / pk.AU
            r_E = [
                ri / pk.AU for ri in self.earth.eph(pk.epoch(bwd_grid[-i - 2]))[0]]
            # TODO ==== rename r_E
            # r_E = [
            #     ri / pk.AU for ri in self.target.eph(pk.epoch(bwd_grid[-i - 2]))[0]]
            dist_earth[-i - 2] = np.linalg.norm(
                [r_E[0] - xbwd[-i - 2], r_E[1] - ybwd[-i - 2], r_E[2] - zbwd[-i - 2]])
            dist_sun[-i -
                     2] = np.linalg.norm([xbwd[-i - 2], ybwd[-i - 2], zbwd[-i - 2]])

        axes.set_xlabel("t [mjd2000]")
        # Plot Earth distance
        axes.plot(times, dist_earth, c='b', label="S/C-Earth")
        # Plot Sun distance
        axes.plot(times, dist_sun, c='gold', label="S/C-Sun")
        axes.set_ylabel("distance [AU]", color='k')
        axes.set_ylim(bottom=0.)
        axes.tick_params('y', colors='k')
        axes.legend(loc=2)

        # Plot thrust profile
        axes = axes.twinx()
        if self.sep:
            max_thrust = self.__sc.thrust
            thrusts = np.linalg.norm(
                np.array(ufwd + ubwd), axis=1) / max_thrust
            # plot maximum thrust achievable at that distance from the Sun
            distsSun = dist_sun[:fwd_seg] + \
                       dist_sun[-bwd_seg:] + [dist_sun[-1]]
            Tmaxs = [self._sep_model(d)[0] / max_thrust for d in distsSun]
            axes.step(np.concatenate(
                (fwd_grid, bwd_grid[1:])), Tmaxs, where="post", c='lightgray', linestyle=':')
        else:
            thrusts = throttles.copy()
        # duplicate the last for plotting
        thrusts = np.append(thrusts, thrusts[-1])
        axes.step(np.concatenate(
            (fwd_grid, bwd_grid[1:])), thrusts, where="post", c='k', linestyle='--')
        axes.set_ylabel("T/Tmax$_{1AU}$", color='k')
        axes.tick_params('y', colors='k')
        axes.set_xlim([times[0], times[-1]])
        axes.set_ylim([0, max(thrusts) + 0.2])

        # print(f"{'max dist'}\t{max(dist_sun)}")
        # print(f"{'min dist'}\t{min(dist_sun)}")

        return axes

    def pretty(self, x):
        # TODO === fix parameter readout for return data
        """
        prob.pretty(x)
        Args:
            - x (``list``, ``tuple``, ``numpy.ndarray``): Decision chromosome, e.g. (``pygmo.population.champion_x``).
        Prints human readable information on the trajectory represented by the decision vector x
        """
        m_p = x[1]
        m_f = self.return_mass

        if not len(x) == len(self.get_bounds()[0]):
            raise ValueError("Invalid length of the decision vector x")

        n_seg = self.__n_seg
        # self.m_i = self._lv_model((x[2] / 1000) ** 2 + (x[3] / 1000) ** 2 + (x[4] / 1000) ** 2)
        self.m_i = m_f + m_p
        t0 = x[0]
        T = self.tof
        thrusts = [np.linalg.norm(x[8 + 3 * i: 11 + 3 * i])
                   for i in range(n_seg)]

        tf = t0 + T
        # self.mP = self.m_i - self.return_mass
        self.deltaV = self.__sc.isp * pk.G0 * np.log(self.m_i / self.return_mass)

        dt = np.append(self.__fwd_dt, self.__bwd_dt) * T / pk.DAY2SEC
        time_thrusts_on = sum(dt[i] for i in range(
            len(thrusts)) if thrusts[i] > 0.1)

        print(f"{'Departure: '}\t{pk.epoch(t0)}")
        print(f"{'Arrival: '}\t{pk.epoch(tf)}")
        print(f"{'Time of flight: '}\t{T}")
        print(f"{'Delta-v: '}\t{self.deltaV}")
        print(f"{'Spacecraft Initial Mass: '}\t{self.m_i}")
        print(f"{'Spacecraft Final Mass: '}\t{self.return_mass}")
        print(f"{'Propellant Consumed: '}\t{m_p}")
        print(f"{'Engine: '}\t{self.engine.name}")
        print(f"{'Number of active engines: '}\t{self.no_engine}")
        print(f"{'Solar array power at 1 AU BOL (W): '}\t{self.Psa}")
        print(f"{'Thrust-on time: '}\t{time_thrusts_on}")
        print(f"{'Launcher C3 (km^2/s^2): '}\t{(x[2] / 1000) ** 2 + (x[3] / 1000) ** 2 + (x[4] / 1000) ** 2}")
        print(f"{'Departure excess velocity (km/s): '}\t{((x[2] / 1000) ** 2 + (x[3] / 1000) ** 2 + (x[4] / 1000) ** 2) ** 0.5}")
        print(f"{'Arrival velocity: '}\t{((x[5]) ** 2 + (x[6]) ** 2 + (x[7]) ** 2) ** 0.5}\t{x[5]}\t{x[6]}\t{x[7]}")

    def pretty_f(self, f):
        print("=====Tolerance Check=====")
        print("=====mismatch constraints=====")
        print("mf (kg): ", -f[0])
        print("x mismatch (m): ", f[1] * pk.AU)
        print("x mismatch: ", f[1])
        print("y mismatch (m): ", f[2] * pk.AU)
        print("y mismatch: ", f[2])
        print("z mismatch (m): ", f[3] * pk.AU)
        print("z mismatch: ", f[3])
        print("vx mismatch (m/s): ", f[4] * pk.EARTH_VELOCITY)
        print("vx mismatch: ", f[4])
        print("vy mismatch (m/s): ", f[5] * pk.EARTH_VELOCITY)
        print("vy mismatch: ", f[5])
        print("vz mismatch (m/s): ", f[6] * pk.EARTH_VELOCITY)
        print("vz mismatch: ", f[6])
        # TODO ==== make equal value set in constraints in fitness
        print("mass mismatch (kg): ", f[7] * (x[1]+self.return_mass))
        # print("mass mismatch (kg): ", f[7] * self.mass_dim)
        print("mass mismatch: ", f[7])
        # print("====total====")
        # for i in np.arange(0, f.size, 1):
        #     print(f"{i:.0f}\t{f[i]}")

    def data_reporting(self,x,runtime,sp_data, ID):

        import datetime

        # current time
        now = datetime.datetime.now()
        # number of segments
        n_seg = self.__n_seg
        # target departure date
        t0 = x[0]
        # helio time of flight
        T = self.tof
        # mass of helio propellant
        mp = x[1]
        # mass returned to Earth
        mf = self.return_mass
        # mass on target departure
        mi = mp  + mf
        # Earth return epoch
        tf = t0 + T
        # equivalent delta-v
        self.deltaV = self.__sc.isp * pk.G0 * np.log(mi / mf)
        # mass in target orbit
        mass_in_orbit = mi + sp_data[5]

        # data collection for max/min helio distance
        fwd_seg = self.__fwd_seg
        bwd_seg = self.__bwd_seg
        dist_sun = [0.0] * (n_seg + 2)  # distances spacecraft - Sun
        rfwd, rbwd, vfwd, vbwd, mfwd, mbwd, ufwd, ubwd, fwd_dt, bwd_dt, _, _ = self._propagate(x)
        # Forward propagation
        xfwd = [0.0] * (fwd_seg + 1)
        yfwd = [0.0] * (fwd_seg + 1)
        zfwd = [0.0] * (fwd_seg + 1)
        xfwd[0] = rfwd[0][0] / pk.AU
        yfwd[0] = rfwd[0][1] / pk.AU
        zfwd[0] = rfwd[0][2] / pk.AU
        dist_sun[0] = np.linalg.norm([xfwd[0], yfwd[0], zfwd[0]])
        for i in range(fwd_seg):
            xfwd[i + 1] = rfwd[i + 1][0] / pk.AU
            yfwd[i + 1] = rfwd[i + 1][1] / pk.AU
            zfwd[i + 1] = rfwd[i + 1][2] / pk.AU
            dist_sun[
                i + 1] = np.linalg.norm([xfwd[i + 1], yfwd[i + 1], zfwd[i + 1]])
        # Backward propagation
        xbwd = [0.0] * (bwd_seg + 1)
        ybwd = [0.0] * (bwd_seg + 1)
        zbwd = [0.0] * (bwd_seg + 1)
        xbwd[-1] = rbwd[-1][0] / pk.AU
        ybwd[-1] = rbwd[-1][1] / pk.AU
        zbwd[-1] = rbwd[-1][2] / pk.AU
        dist_sun[-1] = np.linalg.norm([xbwd[-1], ybwd[-1], zbwd[-1]])
        for i in range(bwd_seg):
            xbwd[-i - 2] = rbwd[-i - 2][0] / pk.AU
            ybwd[-i - 2] = rbwd[-i - 2][1] / pk.AU
            zbwd[-i - 2] = rbwd[-i - 2][2] / pk.AU
            dist_sun[-i -
                     2] = np.linalg.norm([xbwd[-i - 2], ybwd[-i - 2], zbwd[-i - 2]])


        with open('merlot_return_data.csv', 'a+',newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter=',')
            writer.writerow(  # RUN DATA
                            [ID, str(now), n_seg, runtime,
                             # INPUTS
                             self.target_name, 'Earth', self.Psa, self.engine.name, self.no_engine, mf,
                             # TIMELINE DATA
                             sp_data[4], sp_data[2], pk.epoch(t0), T, pk.epoch(t0+T),
                             # DEPARTURE AND ARRIVAL EXCESS VELOCITY
                             ((x[2]) ** 2 + (x[3]) ** 2 + (x[4]) ** 2) ** 0.5,
                             ((x[5]) ** 2 + (x[6]) ** 2 + (x[7]) ** 2) ** 0.5,
                             # MASS DATA
                             mass_in_orbit, sp_data[5], mi, mp, mf,
                             # DELTA-V
                             sp_data[3], self.deltaV,
                             # SPIRAL THRUST PARAMETERS
                             sp_data[0], sp_data[1],
                             # HELIO DISTANCE
                             max(dist_sun), min(dist_sun)])







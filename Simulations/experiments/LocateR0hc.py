# Author: Nathaniel Ruhl

# This class locates the value of r0_hc, the expected position at crossing that depends on starECI and
# the model orbit during the observation (h_unit particularly... in observation_dict).
# The current algorithm only finds r0_hc for a circular orbit (r_model default argument)

# import standard libraries
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d
import datetime as datetime

# import local modules
import tools


class LocateR0hc:
    t_step_size = 0.01   # second step size between definitions of model orbit position
    n_step_size = 0.1   # km step size along the line-of-sight (LOS)
    max_los_dist = 4000   # km, max distance that we look for graze point along the LOS

    def __init__(self, observation_dict, earth_shape_string, r_model_type='circle', interp_aster_kind="cubic"):
        if r_model_type != 'circle' and r_model_type != 'aster':
            raise RuntimeWarning('This class is only defined for r_model_type="circle" and "aster"')
        if interp_aster_kind != "cubic" and interp_aster_kind != "linear":
            raise RuntimeWarning('This class is only defined for interp_aster_kind="cubic" and "linear"')
        self.interp_aster_kind = interp_aster_kind
        self.obs_dict = observation_dict
        self.earth_shape_string = earth_shape_string
        self.r_model_type = r_model_type
        # The following will be defined in self.get_orbit_vec_circle()
        self.A_2d = None   # Theoretical 2d distance along LOS to graze point
        self.A_3d = None   # Actual 3d distance along LOS to graze point
        self.r0_2d = None  # Expected r0 for a 2d horizon crossing

        if self.r_model_type == "circle":
            self.t_orbit = np.arange(0, int(self.obs_dict['T']), LocateR0hc.t_step_size)
            self.orbit_model = self.get_orbit_vec_circle()    # array of ECI positions for an entire orbit
            self.h_unit_approx = self.obs_dict['h_unit']
            self.R_orbit_approx = self.obs_dict['R_orbit']
        elif self.r_model_type == "aster":
            orbit_model_10sec, v_orbit, utc_aster = self.get_orbit_vec_aster()
            t_orbit_10sec = self.define_time_array_aster(utc_aster)
            self.t_orbit = np.arange(min(t_orbit_10sec), max(t_orbit_10sec)+1, 1)   # longer version
            self.orbit_model = self.interpolate_aster_position_scipy(r_vec=orbit_model_10sec, time_seconds=t_orbit_10sec,
                                                                     new_time=self.t_orbit, kind=self.interp_aster_kind)
            # TODO: We need to be able to define the fields below based on the orbit model.
            # Now, we kind of know where to look
            self.h_unit_approx = self.obs_dict['h_unit']
            self.R_orbit_approx = self.obs_dict['R_orbit']
        self.t0_guess_list, self.r0_guess_list = self.get_t0_guess_indices()   # positions within 5% of r0_2d
        self.t0_model = None  # will be defined in the method below
        self.r0_hc = self.locate_r0_numerical()

    def get_orbit_vec_circle(self):
        # Define the orbital circle in the perifocal frame, then transform to ECI

        # Make use of python broadcasting to fill the perifocal array
        t_orbit_column = self.t_orbit.reshape((len(self.t_orbit), 1))
        x_per = self.obs_dict['R_orbit'] * np.cos(self.obs_dict['OMEGA_ORB'] * t_orbit_column)
        y_per = self.obs_dict['R_orbit'] * np.sin(self.obs_dict['OMEGA_ORB'] * t_orbit_column)
        z_per = np.zeros(t_orbit_column.shape)
        orbit_vec_per = np.hstack((x_per, y_per, z_per))

        orbit_vec_eci = np.zeros(orbit_vec_per.shape)

        # I don't know how to do this operation with broadcasting
        for ti, time in enumerate(self.t_orbit):
            orbit_vec_eci[ti] = np.dot(np.matrix.transpose(self.obs_dict['Q']), orbit_vec_per[ti])

        return orbit_vec_eci

    def get_orbit_vec_aster(self):
        df = pd.read_csv("~/Desktop/ISS_03Feb2020_Orbit/HorizonCrossingOrbit_GMATOutput2020.txt", delimiter="\s+")
        utc = df['DefaultSC.UTCGregorian'].to_numpy()
        rx = df['DefaultSC.EarthMJ2000Eq.X'].to_numpy()
        ry = df['DefaultSC.EarthMJ2000Eq.Y'].to_numpy()
        rz = df['DefaultSC.EarthMJ2000Eq.Z'].to_numpy()
        vx = df['DefaultSC.EarthMJ2000Eq.VX'].to_numpy()
        vy = df['DefaultSC.EarthMJ2000Eq.VY'].to_numpy()
        vz = df['DefaultSC.EarthMJ2000Eq.VZ'].to_numpy()
        r_vec = np.vstack((rx, ry, rz)).T  # km
        v_vec = np.vstack((vx, vy, vz)).T  # km/s
        return r_vec, v_vec, utc

    # this function is a helper for reading the aster labs model
    def define_time_array_aster(self, utc_array):
        time_array = []  # sec since the beginning of the day
        for i in range(len(utc_array)):
            time_list = utc_array[i].split(":")
            time_seconds = float(time_list[0]) * 3600 + float(time_list[1]) * 60 + float(time_list[2])
            time_array.append(time_seconds)
        return np.array(time_array)

    def interpolate_aster_position_scipy(self, r_vec, time_seconds, new_time, kind):
        x_interp = interp1d(time_seconds, r_vec[:, 0], kind=kind)
        y_interp = interp1d(time_seconds, r_vec[:, 1], kind=kind)
        z_interp = interp1d(time_seconds, r_vec[:, 2], kind=kind)
        r_interpolated = np.vstack((x_interp(new_time), y_interp(new_time), z_interp(new_time))).T
        return r_interpolated

    def get_t0_guess_indices(self):
        # Use the 2d formulas to guess where r0 may be
        g_unit_proj = np.cross(self.obs_dict['starECI_proj'], self.h_unit_approx)
        self.A_2d = np.sqrt(self.R_orbit_approx ** 2 - tools.R_EARTH ** 2)
        self.r0_2d = tools.R_EARTH * g_unit_proj - self.A_2d * self.obs_dict['starECI_proj']
        r0_guess_indices = np.isclose(self.orbit_model, self.r0_2d, 0.05)

        # 1.5% corresponds to ~50km (or more depending on x,y,z component)
        # depends on |r0_3d-r0_2d|

        t0_guess_list = []  # indices in orbit_vec

        for index, value in enumerate(r0_guess_indices):
            if all(value) == True:
                t0_guess_list.append(index)

        # get the positions that corresponds to the t0 list
        # t0 indices are for orbit_vec
        r0_guess_list = self.orbit_model[min(t0_guess_list):max(t0_guess_list)]

        # Plot the algorithm
        # self.plot_locate_r02d(t0_guess_list)

        return t0_guess_list, r0_guess_list

    # Line of sight from the predicted satellite position r(t)
    def los_line(self, time_index, n_list):
        if isinstance(n_list, int) or isinstance(n_list, float):
            # n_list is not a list, but a single number
            n = n_list
            return self.r0_guess_list[time_index] + n * self.obs_dict['starECI']
        else:
            n_column_vec = n_list.reshape((len(n_list), 1))
            starArray = np.ones((len(n_list), 3)) * self.obs_dict['starECI']
            return self.r0_guess_list[time_index] + n_column_vec * starArray

    # Locate r0 via aligning the LOS to be tangent to earth
    def locate_r0_numerical(self):
        # Loop through different times, different lines of sight during the crossing
        for time_index, time in enumerate(self.t0_guess_list):
            # Lists to check radial altitude at different points along the LOS
            n_list = np.arange(0, LocateR0hc.max_los_dist, LocateR0hc.n_step_size)
            los_points = self.los_line(time_index, n_list)  # all points along the LOS

            # Lists to check radial altitude at different points along the LOS
            # Variable below is distance from point to origin, not los length
            los_mag_list = np.sqrt(los_points[:, 0] ** 2 + los_points[:, 1] ** 2 + los_points[:, 2] ** 2)
            if self.earth_shape_string == 'sphere':
                phi_list = np.arccos(
                    los_points[:, 2] / los_mag_list)  # polar angle at every point along the line of sight
                earth_radius_list = np.ones_like(los_mag_list) * tools.R_EARTH
            elif self.earth_shape_string == 'ellipsoid':
                phi_list = np.arccos(
                    los_points[:, 2] / los_mag_list)  # polar angle at every point along the line of sight
                # Find the radius of earth with the same polar angle as points along the line of sight
                earth_points = tools.point_on_earth(np.zeros_like(phi_list), phi_list)
                earth_radius_list = np.sqrt(earth_points[:, 0] ** 2 + earth_points[:, 1] ** 2 + earth_points[:, 2] ** 2)
            else:
                raise RuntimeError("Please specify if earth_shape_string = 'sphere' or 'ellipsoid'")

            # Check if we reached the tangent point, all entries greater than R_EARTH
            # This will be achieved for multiple lines of sight
            if all(los_mag_list >= earth_radius_list):
                # Find the point of closest approach, the tangent point
                n_graze_index = np.argmin(los_mag_list)
                self.A_3d = n_list[n_graze_index]
                # The 2 below definitions are insightful, but not currently being used
                graze_point = los_points[n_graze_index]
                graze_phi = phi_list[n_graze_index]   # polar angle at graze_point
                self.t0_model = time_index / LocateR0hc.t_step_size
                return self.r0_guess_list[time_index]
            else:
                # keep moving through time until the whole LOS is above earth
                continue
        print('Tangent point not located in specified time range')
        return 0, 0, 0

    @classmethod
    def set_t_step_size(cls, t_step):
        cls.t_step_size = t_step

    def plot_locate_r02d(self):
        plt.figure()
        plt.title(r'Locating $r_{0,2d}$ in an orbit model')
        plt.plot(self.t_orbit, self.orbit_model[:, 0], 'r', label=r'$r_x$')
        plt.plot(self.t_orbit, self.orbit_model[:, 1], 'g', label=r'$r_y$')
        plt.plot(self.t_orbit, self.orbit_model[:, 2], 'b', label=r'$r_z$')
        plt.hlines(self.r0_2d[0], 0, self.t_orbit[-1], 'r', linestyle='--', label='$r_{0,2d}(x)$')
        plt.hlines(self.r0_2d[1], 0, self.t_orbit[-1], 'g', linestyle='--', label='$r_{0,2d}(y)$')
        plt.hlines(self.r0_2d[2], 0, self.t_orbit[-1], 'b', linestyle='--', label='$r_{0,2d}(z)$')
        plt.vlines(self.t_orbit[int(min(self.t0_guess_list))], -self.obs_dict['R_orbit'], self.obs_dict['R_orbit'], 'k', linestyle='--',
                   label='Intersection')
        plt.vlines(self.t_orbit[int(max(self.t0_guess_list))], -self.obs_dict['R_orbit'], self.obs_dict['R_orbit'], 'k', linestyle='--')

        plt.xlabel('Seconds From Ascending Node (circle)')
        plt.ylabel('km')
        plt.legend()
        plt.show()
        return


def main():
    from Observations.Crab.obs4522010103 import crab03_dict as od
    # from Observations.V4641_Feb3 import feb3_dict as od
    # ellipsoid earth
    r0hc_ellipsoid_obj = LocateR0hc(observation_dict=od.crab03_dict, earth_shape_string='ellipsoid',
                                    r_model_type='circle')
    print(f'r0_hc_ellipsoid = {r0hc_ellipsoid_obj.r0_hc}')
    print(r0hc_ellipsoid_obj.t0_model)
    r0hc_ellipsoid_obj.plot_locate_r02d()

    # spherical earth
    # r0hc_sphere_obj = LocateR0hc(observation_dict=od.crab03_dict, earth_shape_string='sphere', r_model_type='circle')
    # print(f'r0_hc_sphere = {r0hc_sphere_obj.r0_hc}')

    # print(f"r0_2d: = {r0hc_sphere_obj.r0_2d}")
    return 0


if __name__ == '__main__':
    import time
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))
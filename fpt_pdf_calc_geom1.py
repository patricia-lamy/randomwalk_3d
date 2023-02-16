import numpy as np
import pandas as pd
import time
from fpt_calc_methods import reflect

# We set the random seed
# np.random.seed(42)

# Parameters for outputs wanted
n_particles = 10000  # Number of simulations

# Setting the conditions of the problem.
sphere_radius = 1.02       # Sphere radius, also domain size (cm)
diff_coef = 3.852e-3       # Diffusion coefficient (cm2/h)

# Angle (rad) defining the lens radius (from the top of the geometry):
phi_lens = np.pi - 2.6

# Angle (rad) defining the intersection of the hyaloid membrane with the ILM
# (from the top of the geometry):
phi_ilm = np.pi - 2.2

# Permeability parameters:
perm_param_hya = 0.06876  # cm/h, permeability parameter for hyaloid membrane
perm_param_ilm = 6.516e-4  # cm/h, permeability parameter for ILM

# Initialisation for simulations
time_max = 1500  # Time (in h) when the simulation stops
delta_t = 0.001  # Time step length (h)

# Corresponding variables for the simulation (do not change)
max_time_steps = int(time_max/delta_t)   # Nb of max time steps in simulation.
step_length = np.sqrt(6*diff_coef*delta_t)  # Step length (cm) of random walk
print('Step length is defined as:', step_length, ' cm')
# Based on Erban and Chapman 2007:
prob_absorb_hya = (np.sqrt(delta_t) *
                   (perm_param_hya*np.sqrt(np.pi))/(np.sqrt(diff_coef)))
prob_absorb_ilm = (np.sqrt(delta_t) *
                   (perm_param_ilm*np.sqrt(np.pi))/(np.sqrt(diff_coef)))

# Initial position: center of the sphere
rho0 = 0
theta0 = 0
phi0 = 0

# Initialisation of array collecting the first passage times
fpt_list = np.empty(n_particles)


# Initialisation of array collecting membrane where particle exited
# If exited through hya: 0
# If exited through ilm: 1
# If not exited: nan
exit_point = np.empty(n_particles)

start = time.time()
for i in range(n_particles):
    # Initial position of particle
    rho = rho0
    theta = theta0
    phi = phi0

    # Setting time step counter to first time step
    t = delta_t
    t_index = 1

    while t < time_max:
        # Sampling random angles for variation of theta and phi
        d_theta = (np.random.uniform() * (2 * np.pi))
        cos_d_phi = np.random.uniform(-1, 1)
        sin_d_phi = np.sqrt(1 - cos_d_phi**2)

        # Calculating current cartesian position:
        x1 = rho*np.sin(phi)*np.cos(theta)
        y1 = rho*np.sin(phi)*np.sin(theta)
        z1 = rho*np.cos(phi)
        # Calculating next cartesian position:
        xt = x1 + step_length*sin_d_phi*np.cos(d_theta)
        yt = y1 + step_length*sin_d_phi*np.sin(d_theta)
        zt = z1 + step_length*cos_d_phi

        # Calculating the corresponding spherical position particle at time t
        rhot = np.sqrt(xt**2 + yt**2 + zt**2 + 1e-14)
        thetat = np.arctan2(yt, xt)
        phit = np.arccos(zt/rhot)

        if rhot > sphere_radius:  # If the particle hits the boundary at time t
            if phit < phi_lens:  # Lens region
                [rho_r, theta_r, phi_r] = reflect((rho, theta, phi),
                                                  (rhot, thetat, phit),
                                                  sphere_radius)[0]
                [rhot, thetat, phit] = [rho_r, theta_r, phi_r]

            elif phi_lens <= phit <= phi_ilm:  # Region of hyaloid membrane
                random_num = np.random.uniform()
                if random_num >= prob_absorb_hya:  # Reflected
                    [rho_r, theta_r, phi_r] = reflect((rho, theta, phi),
                                                      (rhot, thetat, phit),
                                                      sphere_radius)[0]
                    [rhot, thetat, phit] = [rho_r, theta_r, phi_r]
                else:   # Absorbed
                    fpt_list[i] = t  # Record fpt in the list for particle i
                    exit_point[i] = 0  # Records particle exited through hya
                    break

            elif phit > phi_ilm:  # Region of ILM
                random_num = np.random.uniform()
                if random_num >= prob_absorb_ilm:  # Reflected
                    [rho_r, theta_r, phi_r] = reflect((rho, theta, phi),
                                                      (rhot, thetat, phit),
                                                      sphere_radius)[0]
                    [rhot, thetat, phit] = [rho_r, theta_r, phi_r]
                else:  # Absorbed
                    fpt_list[i] = t  # Record fpt in the list for particle i
                    exit_point[i] = 1  # Records that particle exited via ilm
                    break

            else:
                print('shouldnt come here')

        # Resetting indexes
        t = round(t+delta_t, 10)
        t_index = t_index + 1
        rho = rhot
        theta = thetat
        phi = phit

    else:  # If we reach max time step before particle exited
        fpt_list[i] = np.nan
        exit_point[i] = np.nan


end = time.time()
print('Time elapsed (s):', end - start)


# Saving the information
df = pd.DataFrame({'#Fpt': fpt_list,
                   '#Exit point (0:hya and 1:ilm)': exit_point})
df.to_csv("data/fpt_array_geomD_14-02-2023_deltat0.001.csv", index=False)

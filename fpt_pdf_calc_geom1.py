import numpy as np
import pandas as pd
import time
from reflective_bc_functions import reflected_point_method_sph

# We set the random seed
# np.random.seed(42)

# Setting the conditions of the problem.
sphere_radius = 1.02       # Sphere radius, also domain size (cm)
diff_coef = 3.852e-3       # Diffusion coefficient (cm2/h)

# Angle (rad) defining the lens radius (from the top of the geometry):
phi_lens = np.pi - 2.6
#phi_lens = 0   #Temporary for Geom0

# Angle (rad) defining the intersection of the hyaloid membrane with the ILM
# (from the top of the geometry):
phi_ilm = np.pi - 2.2

# Permeability parameters:
perm_param_hya = 0.04284   # cm/h, permeability parameter for hyaloid membrane
perm_param_ilm = 6.516e-4  # cm/h, permeability parameter for ILM
# Supposing that semi-permeable are fully absorbing for now:
perm_param_ilm = 1
perm_param_hya = 1

# Initialisation for simulations
n_particles = 10000  # Number of simulations
time_max = 1500  # Time (in h) when the simulation stops
delta_t = 0.001  # Time step length (h)

# Corresponding variables for the simulation (do not change)
max_time_steps = int(time_max/delta_t)   # Nb of max time steps in simulation.
step_length = np.sqrt(6*diff_coef*delta_t)  # Step length (cm) of random walk
# Based on Erban and Chapman 2007:
prob_absorb_hya = (np.sqrt(delta_t) *
                   (perm_param_hya*np.sqrt(np.pi))/(2*np.sqrt(diff_coef)))
prob_absorb_hya = 1

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

#Temporary thing
#t_list = []

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
        d_theta = (np.random.uniform() * (2 * np.pi)) + 1e-14
        d_phi = (np.random.uniform() * (np.pi)) + 1e-14

        # Calculating current cartesian position:
        x1 = rho*np.sin(phi)*np.cos(theta)
        y1 = rho*np.sin(phi)*np.sin(theta)
        z1 = rho*np.cos(phi)
        # Calculating next cartesian position:
        xt = x1 + step_length*np.sin(d_phi)*np.cos(d_theta)
        yt = y1 + step_length*np.sin(d_phi)*np.sin(d_theta)
        zt = z1 + step_length*np.cos(d_phi)

        # Calculating the corresponding spherical position particle at time t
        rhot = np.sqrt(xt**2 + yt**2 + zt**2 + 1e-14)
        thetat = np.arctan2(yt, xt)
        phit = np.arccos(zt/rhot)

        if rhot > sphere_radius:  # If the particle hits the boundary at time t
            if phit < phi_lens:  # Lens region
                reflected_point = reflected_point_method_sph((rho, theta, phi),
                                                             (rhot, thetat,
                                                              phit),
                                                             sphere_radius)
                rhot = reflected_point[0]
                thetat = reflected_point[1]
                phit = reflected_point[2]
                if rhot >= sphere_radius:
                    print('New point not in sphere, rho =', rhot)
                    rhot = sphere_radius
            elif phi_lens <= phit <= phi_ilm:  # Region of hyaloid membrane
                random_num = np.random.uniform()
                if random_num >= prob_absorb_hya:
                    print('shouldnt be here1')
                    # We suppose the particle did not move (BBC)
                    rhot = rho
                    thetat = theta
                    phit = phi
                else:  # If absorbed
                    fpt_list[i] = t  # Record fpt in the list for particle i
                    exit_point[i] = 0  # Records particle exited through hya
                    break
            elif phit > phi_ilm:  # Region of ILM
                random_num = np.random.uniform()
                if random_num >= perm_param_ilm:
                    print('shouldnt be here2')
                    rhot = sphere_radius  # Effectively a reflective boundary
                else:  # If absorbed
                    fpt_list[i] = t  # Record fpt in the list for particle i
                    exit_point[i] = 1  # Records that the particle exited through ilm
                    break
            else:
                print('shouldnt come here')

        # Resetting indexes
        t = round(t+delta_t, 10)
        t_index += 1
        #t_list.append(t)
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
df.to_csv("fpt_array_geom0.1_23-01-2023_deltat0.001.csv", index=False)

#df_time = pd.DataFrame({'Time points': t_list})
#df_time.to_csv("time_array_geom0_20-01-2023_tests.csv", index=False)


#np.savetxt("fpt_array_geom1_0.01deltat_02-12-22.csv", fpt_list, delimiter=",",
#header='First passage time (in hours) for off-lattice random walk simulation
#on Geometry 1 (reflectice boundary condition on lens surface),
#with delta_t = 0.01 hours.')

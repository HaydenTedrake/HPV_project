import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import gala.integrate as gi
import gala.dynamics as gd
import gala.potential as gp
from gala.potential import NFWPotential
from gala.units import galactic

# define potential
# from tutorial
# pot = gp.NFWPotential.from_circular_velocity(v_c=200*u.km/u.s, r_s=10.*u.kpc, units=galactic)
# other method?
pot = NFWPotential.from_M200_c(M200=1e12 * u.M_sun, c=10.0, units=galactic)

# set initial conditions and specify time-stepping
ics = gd.PhaseSpacePosition(pos=[10,0,0.] * u.kpc, vel=[0,175,0] * u.km/u.s)
orbit = pot.integrate_orbit(ics, dt=1*u.Myr, n_steps=100)

# integrating orbits in parallel
norbits = 128
new_pos = np.random.normal(ics.pos.xyz.to(u.pc).value, 100., size=(norbits,3)).T * u.pc
new_vel = np.random.normal(ics.vel.d_xyz.to(u.km/u.s).value, 1., size=(norbits,3)).T * u.km/u.s
new_ics = gd.PhaseSpacePosition(pos=new_pos, vel=new_vel)
orbits = pot.integrate_orbit(new_ics, dt=1., n_steps=100)

# # plot
# grid = np.linspace(-15,15,64)
# fig,ax = plt.subplots(1, 1, figsize=(5,5))
# fig = pot.plot_contours(grid=(grid,grid,0), cmap='Greys', ax=ax)
# fig = orbits[-1].plot(['x', 'y'], color='#9ecae1', s=1., alpha=0.5, axes=[ax], auto_aspect=False)
# plt.show()


# # plot
# grid = np.linspace(-15,15,64)
# fig,ax = plt.subplots(1, 1, figsize=(5,5))
# fig = pot.plot_contours(grid=(grid,grid,0), cmap='Greys', ax=ax)
# fig = orbits[:,0].plot(['x', 'y'],)
# plt.show()

# print(orbits.shape)

# Plot the first particle
times = np.arange(orbits.ntimes)
pos_output = orbits[:,0].pos.xyz.value

plt.figure()
plt.scatter(pos_output[0,:], pos_output[1,:],
            c=times, cmap='Spectral',s=1)
plt.colorbar()
plt.xlabel('X Value')
plt.ylabel('Y Value')
plt.title('Y vs X')
# plt.show()

print((1*u.kpc/(500*u.km/u.s)).to(u.Myr))
print((1*u.kpc/(1*u.Myr)).to(u.km/u.s))
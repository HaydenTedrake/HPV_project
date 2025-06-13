import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import gala.integrate as gi
import gala.dynamics as gd
import gala.potential as gp
from gala.potential import NFWPotential
from gala.units import galactic

# define potential
pot = NFWPotential.from_M200_c(M200=1e12 * u.M_sun, c=10.0, units=galactic)

# set initial conditions and specify time-stepping
ics = gd.PhaseSpacePosition(pos=[10,0,0.] * u.kpc, vel=[0,175,0] * u.km/u.s)
orbit = pot.integrate_orbit(ics, dt=1*u.Myr, n_steps=100)

# integrating orbits in parallel
norbits = 128
new_pos = np.random.normal(ics.pos.xyz.to(u.pc).value, 100., size=(norbits,3)).T * u.pc
new_vel = np.random.normal(ics.vel.d_xyz.to(u.km/u.s).value, 1., size=(norbits,3)).T * u.km/u.s
new_ics = gd.PhaseSpacePosition(pos=new_pos, vel=new_vel)
orbits = pot.integrate_orbit(new_ics, dt=1*u.Myr, n_steps=100)

# plot the first particle
times = np.arange(orbits.ntimes)
pos_output = orbits[:,0].pos.xyz.value

plt.figure()
plt.scatter(pos_output[0,:], pos_output[1,:],
            c=times, cmap='Spectral',s=1)
plt.colorbar()
plt.xlabel('X Value')
plt.ylabel('Y Value')
plt.title('Y vs X')
plt.show()
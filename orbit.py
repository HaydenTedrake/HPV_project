import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import gala.dynamics as gd
import gala.potential as gp
from gala.potential import NFWPotential
from gala.units import galactic

# Define potential
pot = NFWPotential.from_M200_c(M200=1e12 * u.M_sun, c=10.0, units=galactic)

# Initial conditions
ics = gd.PhaseSpacePosition(pos=[10, 0, 0] * u.kpc, vel=[0, 175, 0] * u.km/u.s)

# Define time step values to compare
dt_values = [1, 10, 50]  # in Myr
colors = ['blue', 'orange', 'red']
labels = ['dt = 1 Myr', 'dt = 10 Myr', 'dt = 50 Myr']

plt.figure()

for dt_val, color, label in zip(dt_values, colors, labels):
    dt = dt_val * u.Myr
    n_steps = int(500 * u.Myr / dt)
    orbit = pot.integrate_orbit(ics, dt=dt, n_steps=n_steps)

    pos = orbit.pos.xyz.value
    plt.scatter(pos[0], pos[1], color=color, label=label, s=5)

plt.xlabel('X [kpc]')
plt.ylabel('Y [kpc]')
plt.title('Comparison of Orbits with Varying Time Steps')
plt.legend()
plt.axis('equal')
plt.grid(True)
plt.show()
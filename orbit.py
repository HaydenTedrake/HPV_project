import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import gala.integrate as gi
import gala.dynamics as gd
import gala.potential as gp
from gala.units import galactic

pot = gp.NFWPotential(mvir = 1e12 * u.M_sun, conc=15.0, units=galactic)
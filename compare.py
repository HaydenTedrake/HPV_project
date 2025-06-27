# Imports
import astropy.units as u
import astropy.coordinates as coord
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import numpy as np
import gala.integrate as gi
import gala.dynamics as gd
import gala.potential as gp
from gala.units import galactic
import scipy.optimize as opt

# Define potential
pot = gp.MilkyWayPotential()

# Stars
hvs5 = {
    'ra': 139.4978105603, #degree
    'e_ra': 0.0001, #degree
    'dec': 67.3773227272, #degree
    'e_dec': 0.0001, #degree
    'dist': 44.20, #kpc
    'e_dist': 5.09, #kpc
    'pmra': 0.00, #mas/yr
    'e_pmra': 0.08, #mas/yr
    'pmdec': -0.99, #mas/yr
    'e_pmdec': 0.11, #mas/yr
    'rv': 545.50, #km/s
    'e_rv': 4.30, #km/s
    }
hvs5['SkyCoord'] = coord.SkyCoord(
    ra=hvs5['ra'] * u.degree,
    dec=hvs5['dec'] * u.degree,
    distance=hvs5['dist'] * u.kpc,
    pm_ra_cosdec=hvs5['pmra'] * (u.mas / u.yr),
    pm_dec=hvs5['pmdec'] * (u.mas / u.yr),
    radial_velocity=hvs5['rv'] * (u.km / u.s),
    frame="icrs"
)
hvs5['GC'] = hvs5['SkyCoord'].transform_to(coord.Galactocentric)
hvs5['ics'] = gd.PhaseSpacePosition(pos=hvs5['GC'].cartesian.xyz, vel=hvs5['GC'].velocity.d_xyz)
hvs5['orbit1'] = pot.integrate_orbit(hvs5['ics'], dt=-0.1*u.Myr, n_steps=10000) # 1 Gyr
hvs5['orbit2'] = pot.integrate_orbit(hvs5['ics'], dt=-0.1*u.Myr, n_steps=5000) # 500 Myr
hvs5['orbit3'] = pot.integrate_orbit(hvs5['ics'], dt=-0.1*u.Myr, n_steps=1000) # 100 Myr
hvs8 = {
    'ra': 145.5584849227, #degree
    'e_ra': 0.0001, #degree
    'dec': 20.0561234065, #degree
    'e_dec': 0.0001, #degree
    'dist': 53.19, #kpc
    'e_dist': 9.80, #kpc
    'pmra': -0.88, #mas/yr
    'e_pmra': 0.16, #mas/yr
    'pmdec': -0.28, #mas/yr
    'e_pmdec': 0.14, #mas/yr
    'rv': 499.30, #km/s
    'e_rv': 2.90, #km/s
    }
hvs8['SkyCoord'] = coord.SkyCoord(
    ra=hvs8['ra'] * u.degree,
    dec=hvs8['dec'] * u.degree,
    distance=hvs8['dist'] * u.kpc,
    pm_ra_cosdec=hvs8['pmra'] * (u.mas / u.yr),
    pm_dec=hvs8['pmdec'] * (u.mas / u.yr),
    radial_velocity=hvs8['rv'] * (u.km / u.s),
    frame="icrs"
)
hvs8['GC'] = hvs8['SkyCoord'].transform_to(coord.Galactocentric)
hvs8['ics'] = gd.PhaseSpacePosition(pos=hvs8['GC'].cartesian.xyz, vel=hvs8['GC'].velocity.d_xyz)
hvs8['orbit1'] = pot.integrate_orbit(hvs8['ics'], dt=-0.1*u.Myr, n_steps=10000) # 1 Gyr
hvs8['orbit2'] = pot.integrate_orbit(hvs8['ics'], dt=-0.1*u.Myr, n_steps=5000) # 500 Myr
hvs8['orbit3'] = pot.integrate_orbit(hvs8['ics'], dt=-0.1*u.Myr, n_steps=1000) # 100 Myr
hvs14 = {
    'ra': 161.0072805184, #degree
    'e_ra': 0.0001, #degree
    'dec': 6.1941762509, #degree
    'e_dec': 0.0001, #degree
    'dist': 102.66, #kpc
    'e_dist': 16.55, #kpc
    'pmra': -2.17, #mas/yr
    'e_pmra': 1.38, #mas/yr
    'pmdec': 2.28, #mas/yr
    'e_pmdec': 1.68, #mas/yr
    'rv': 537.30, #km/s
    'e_rv': 7.20, #km/s
    }
hvs14['SkyCoord'] = coord.SkyCoord(
    ra=hvs14['ra'] * u.degree,
    dec=hvs14['dec'] * u.degree,
    distance=hvs14['dist'] * u.kpc,
    pm_ra_cosdec=hvs14['pmra'] * (u.mas / u.yr),
    pm_dec=hvs14['pmdec'] * (u.mas / u.yr),
    radial_velocity=hvs14['rv'] * (u.km / u.s),
    frame="icrs"
)
hvs14['GC'] = hvs14['SkyCoord'].transform_to(coord.Galactocentric)
hvs14['ics'] = gd.PhaseSpacePosition(pos=hvs14['GC'].cartesian.xyz, vel=hvs14['GC'].velocity.d_xyz)
hvs14['orbit1'] = pot.integrate_orbit(hvs14['ics'], dt=-0.1*u.Myr, n_steps=10000) # 1 Gyr
hvs14['orbit2'] = pot.integrate_orbit(hvs14['ics'], dt=-0.1*u.Myr, n_steps=5000) # 500 Myr
hvs14['orbit3'] = pot.integrate_orbit(hvs14['ics'], dt=-0.1*u.Myr, n_steps=1000) # 100 Myr
hvs17 = {
    'ra': 250.4849449351, #degree
    'e_ra': 0.0001, #degree
    'dec': 47.3961264077, #degree
    'e_dec': 0.0001, #degree
    'dist': 49.82, #kpc
    'e_dist': 3.90, #kpc
    'pmra': -1.13, #mas/yr
    'e_pmra': 0.09, #mas/yr
    'pmdec': -0.93, #mas/yr
    'e_pmdec': 0.10, #mas/yr
    'rv': 250.20, #km/s
    'e_rv': 2.90, #km/s
    }
hvs17['SkyCoord'] = coord.SkyCoord(
    ra=hvs17['ra'] * u.degree,
    dec=hvs17['dec'] * u.degree,
    distance=hvs17['dist'] * u.kpc,
    pm_ra_cosdec=hvs17['pmra'] * (u.mas / u.yr),
    pm_dec=hvs17['pmdec'] * (u.mas / u.yr),
    radial_velocity=hvs17['rv'] * (u.km / u.s),
    frame="icrs"
)
hvs17['GC'] = hvs17['SkyCoord'].transform_to(coord.Galactocentric)
hvs17['ics'] = gd.PhaseSpacePosition(pos=hvs17['GC'].cartesian.xyz, vel=hvs17['GC'].velocity.d_xyz)
hvs17['orbit1'] = pot.integrate_orbit(hvs17['ics'], dt=-0.1*u.Myr, n_steps=10000) # 1 Gyr
hvs17['orbit2'] = pot.integrate_orbit(hvs17['ics'], dt=-0.1*u.Myr, n_steps=5000) # 500 Myr
hvs17['orbit3'] = pot.integrate_orbit(hvs17['ics'], dt=-0.1*u.Myr, n_steps=1000) # 100 Myr
hvs23 = {
    'ra': 329.1209186025, #degree
    'e_ra': 0.0001, #degree
    'dec': 0.9122807743, #degree
    'e_dec': 0.0001, #degree
    'dist': 114.87, #kpc
    'e_dist': 20.10, #kpc
    'pmra': -1.21, #mas/yr
    'e_pmra': 1.29, #mas/yr
    'pmdec': -2.46, #mas/yr
    'e_pmdec': 1.50, #mas/yr
    'rv': 259.30, #km/s
    'e_rv': 9.80, #km/s
    }
hvs23['SkyCoord'] = coord.SkyCoord(
    ra=hvs23['ra'] * u.degree,
    dec=hvs23['dec'] * u.degree,
    distance=hvs23['dist'] * u.kpc,
    pm_ra_cosdec=hvs23['pmra'] * (u.mas / u.yr),
    pm_dec=hvs23['pmdec'] * (u.mas / u.yr),
    radial_velocity=hvs23['rv'] * (u.km / u.s),
    frame="icrs"
)
hvs23['GC'] = hvs23['SkyCoord'].transform_to(coord.Galactocentric)
hvs23['ics'] = gd.PhaseSpacePosition(pos=hvs23['GC'].cartesian.xyz, vel=hvs23['GC'].velocity.d_xyz)
hvs23['orbit1'] = pot.integrate_orbit(hvs23['ics'], dt=-0.1*u.Myr, n_steps=10000) # 1 Gyr
hvs23['orbit2'] = pot.integrate_orbit(hvs23['ics'], dt=-0.1*u.Myr, n_steps=5000) # 500 Myr
hvs23['orbit3'] = pot.integrate_orbit(hvs23['ics'], dt=-0.1*u.Myr, n_steps=1000) # 100 Myr

# Dwarf galaxy
dict = {
    'ra': 152.117175, #degree
    'e_ra': 0.0001, #degree
    'dec': 12.3065, #degree
    'e_dec': 0.0001, #degree
    'dist': 254.0, #kpc
    'e_dist': 15.5, #kpc
    'pmra': -0.007, #mas/yr
    'e_pmra': 0.035, #mas/yr
    'pmdec': -0.119, #mas/yr not updated
    'e_pmdec': 0.026, #mas/yr not updated
    'rv': 282.9, #km/s
    'e_rv': 0.5, #km/s
    }
dwarf = coord.SkyCoord(
    ra=dict['ra'] * u.degree,
    dec=dict['dec'] * u.degree,
    distance=dict['dist'] * u.kpc,
    pm_ra_cosdec=dict['pmra'] * (u.mas / u.yr),
    pm_dec=dict['pmdec'] * (u.mas / u.yr),
    radial_velocity=dict['rv'] * (u.km / u.s),
    frame="icrs"
)
dwarf_GC = dwarf.transform_to(coord.Galactocentric)

dwarf_ics = gd.PhaseSpacePosition(pos=dwarf_GC.cartesian.xyz, vel=dwarf_GC.velocity.d_xyz)
dwarf_orbit1 = pot.integrate_orbit(dwarf_ics, dt=-0.1*u.Myr, n_steps=10000) # 1 Gyr
dwarf_orbit2 = pot.integrate_orbit(dwarf_ics, dt=-0.1*u.Myr, n_steps=5000) # 500 Myr
dwarf_orbit3 = pot.integrate_orbit(dwarf_ics, dt=-0.1*u.Myr, n_steps=1000) # 100 Myr

n_samples = 100
ra_samples = np.random.normal(dict['ra'], dict['e_ra'], n_samples)
dec_samples = np.random.normal(dict['dec'], dict['e_dec'], n_samples)
dist_samples = np.random.normal(dict['dist'], dict['e_dist'], n_samples)
pmra_samples = np.random.normal(dict['pmra'], dict['e_pmra'], n_samples)
pmdec_samples = np.random.normal(dict['pmdec'], dict['e_pmdec'], n_samples)
rv_samples = np.random.normal(dict['rv'], dict['e_rv'], n_samples)

ics_list = []
for ra, dec, dist, pmra, pmdec, rv in zip(ra_samples, dec_samples, dist_samples, pmra_samples, pmdec_samples, rv_samples):
    sc = coord.SkyCoord(ra=ra * u.degree,
                  dec=dec * u.degree,
                  distance=dist * u.kpc,
                  pm_ra_cosdec=pmra * (u.mas / u.yr),
                  pm_dec=pmdec * (u.mas / u.yr),
                  radial_velocity=rv * (u.km / u.s),
                  frame="icrs")
    gc = sc.transform_to(coord.Galactocentric)
    ics = gd.PhaseSpacePosition(pos=gc.cartesian.xyz, vel=gc.velocity.d_xyz)
    ics_list.append(ics)

orbits1 = []
orbits2 = []
orbits3 = []
orbits1.append(dwarf_orbit1)
orbits2.append(dwarf_orbit2)
orbits3.append(dwarf_orbit3)
for ics in ics_list:
    orbit1 = pot.integrate_orbit(ics, dt=-0.1*u.Myr, n_steps=10000)
    orbit2 = pot.integrate_orbit(ics, dt=-0.1*u.Myr, n_steps=5000)
    orbit3 = pot.integrate_orbit(ics, dt=-0.1*u.Myr, n_steps=1000)
    orbits1.append(orbit1)
    orbits2.append(orbit2)
    orbits3.append(orbit3)

# Define comparing technique
def compare(x_vals1, y_vals1, x_vals2, y_vals2):
    pos1 = np.column_stack((x_vals1, y_vals1))
    mean1 = np.mean(pos1, axis=0)
    cov1 = np.cov(pos1, rowvar=False)

    pos2 = np.column_stack((x_vals2, y_vals2))
    mean2 = np.mean(pos2, axis=0)
    cov2 = np.cov(pos2, rowvar=False)

    dx = mean1[0] - mean2[0]
    dy = mean1[1] - mean2[1]

    combined_sigma_x = np.sqrt(cov1[0, 0] + cov2[0, 0])
    combined_sigma_y = np.sqrt(cov1[1, 1] + cov2[1, 1])

    Dx = dx / combined_sigma_x
    Dy = dy / combined_sigma_y

    return Dx, Dy

# 1 Gyr in the past
red_points = []

for i, orbit in enumerate(orbits1):
    final_pos = orbit[-1].pos.xyz.to(u.kpc).value[:3]
    if i == 0:
        pass
    else:
        red_points.append(final_pos)

red_points = np.array(red_points)

x = red_points[:, 0]
y = red_points[:, 1]
z = red_points[:, 2]

for star in [hvs5, hvs8, hvs14, hvs17, hvs23]:
    pass
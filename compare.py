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
dwarf_ra_samples = np.random.normal(dict['ra'], dict['e_ra'], n_samples)
dwarf_dec_samples = np.random.normal(dict['dec'], dict['e_dec'], n_samples)
dwarf_dist_samples = np.random.normal(dict['dist'], dict['e_dist'], n_samples)
dwarf_pmra_samples = np.random.normal(dict['pmra'], dict['e_pmra'], n_samples)
dwarf_pmdec_samples = np.random.normal(dict['pmdec'], dict['e_pmdec'], n_samples)
dwarf_rv_samples = np.random.normal(dict['rv'], dict['e_rv'], n_samples)

dwarf_ics_list = []
for ra, dec, dist, pmra, pmdec, rv in zip(dwarf_ra_samples, dwarf_dec_samples, dwarf_dist_samples, dwarf_pmra_samples, dwarf_pmdec_samples, dwarf_rv_samples):
    sc = coord.SkyCoord(ra=ra * u.degree,
                  dec=dec * u.degree,
                  distance=dist * u.kpc,
                  pm_ra_cosdec=pmra * (u.mas / u.yr),
                  pm_dec=pmdec * (u.mas / u.yr),
                  radial_velocity=rv * (u.km / u.s),
                  frame="icrs")
    gc = sc.transform_to(coord.Galactocentric)
    ics = gd.PhaseSpacePosition(pos=gc.cartesian.xyz, vel=gc.velocity.d_xyz)
    dwarf_ics_list.append(ics)

dwarf_orbits1 = []
dwarf_orbits2 = []
dwarf_orbits3 = []
dwarf_orbits1.append(dwarf_orbit1)
dwarf_orbits2.append(dwarf_orbit2)
dwarf_orbits3.append(dwarf_orbit3)
for ics in dwarf_ics_list:
    orbit1 = pot.integrate_orbit(ics, dt=-0.1*u.Myr, n_steps=10000)
    orbit2 = pot.integrate_orbit(ics, dt=-0.1*u.Myr, n_steps=5000)
    orbit3 = pot.integrate_orbit(ics, dt=-0.1*u.Myr, n_steps=1000)
    dwarf_orbits1.append(orbit1)
    dwarf_orbits2.append(orbit2)
    dwarf_orbits3.append(orbit3)

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

    return (Dx, Dy)

map = {
    0: 'hvs5',
    1: 'hvs8',
    2: 'hvs14',
    3: 'hvs17',
    4: 'hvs23',
}

# Calculate
for num, star in enumerate([hvs5, hvs8, hvs14, hvs17, hvs23]):
    SkyCoord = coord.SkyCoord(
    ra=star['ra'] * u.degree,
    dec=star['dec'] * u.degree,
    distance=star['dist'] * u.kpc,
    pm_ra_cosdec=star['pmra'] * (u.mas / u.yr),
    pm_dec=star['pmdec'] * (u.mas / u.yr),
    radial_velocity=star['rv'] * (u.km / u.s),
    frame="icrs"
    )
    GC = SkyCoord.transform_to(coord.Galactocentric)
    ics = gd.PhaseSpacePosition(pos=GC.cartesian.xyz, vel=GC.velocity.d_xyz)
    orbit1 = pot.integrate_orbit(ics, dt=-0.1*u.Myr, n_steps=10000) # 1 Gyr
    orbit2 = pot.integrate_orbit(ics, dt=-0.1*u.Myr, n_steps=5000) # 500 Myr
    orbit3 = pot.integrate_orbit(ics, dt=-0.1*u.Myr, n_steps=1000) # 100 Myr
    ra_samples = np.random.normal(star['ra'], star['e_ra'], n_samples)
    dec_samples = np.random.normal(star['dec'], star['e_dec'], n_samples)
    dist_samples = np.random.normal(star['dist'], star['e_dist'], n_samples)
    pmra_samples = np.random.normal(star['pmra'], star['e_pmra'], n_samples)
    pmdec_samples = np.random.normal(star['pmdec'], star['e_pmdec'], n_samples)
    rv_samples = np.random.normal(star['rv'], star['e_rv'], n_samples)

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
    orbits1.append(orbit1)
    orbits2.append(orbit2)
    orbits3.append(orbit3)
    for ics in ics_list:
        orbit1 = pot.integrate_orbit(ics, dt=-0.1*u.Myr, n_steps=10000)
        orbit2 = pot.integrate_orbit(ics, dt=-0.1*u.Myr, n_steps=5000)
        orbit3 = pot.integrate_orbit(ics, dt=-0.1*u.Myr, n_steps=1000)
        orbits1.append(orbit1)
        orbits2.append(orbit2)
        orbits3.append(orbit3)
        
    points = []
    for i, orbit in enumerate(orbits1):
        final_pos = orbit[-1].pos.xyz.to(u.kpc).value[:3]
        if i == 0:
            pass
        else:
            points.append(final_pos)

    points = np.array(points)

    x = points[:, 0]
    y = points[:, 1]
    z = points[:, 2]

    # 1 Gyr
    dwarf_points = []

    for i, orbit in enumerate(dwarf_orbits1):
        final_pos = orbit[-1].pos.xyz.to(u.kpc).value[:3]
        dwarf_points.append(final_pos)

    dwarf_points = np.array(dwarf_points)

    dwarf_x = dwarf_points[:, 0]
    dwarf_y = dwarf_points[:, 1]
    dwarf_z = dwarf_points[:, 2]

    XY = compare(dwarf_x, dwarf_y, x, y)
    distance_xy = np.sqrt(XY[0]**2 + XY[1]**2)
    XZ = compare(dwarf_x, dwarf_z, x, z)
    distance_xz = np.sqrt(XZ[0]**2 + XZ[1]**2)
    YZ = compare(dwarf_y, dwarf_z, y, z)
    distance_yz = np.sqrt(YZ[0]**2 + YZ[1]**2)

    print(f"1 Gyr, {map[num]}, XY: {distance_xy}, XZ: {distance_xz}, YZ: {distance_yz}")

    # 500 Myr
    dwarf_points = []

    for i, orbit in enumerate(dwarf_orbits2):
        final_pos = orbit[-1].pos.xyz.to(u.kpc).value[:3]
        dwarf_points.append(final_pos)

    dwarf_points = np.array(dwarf_points)

    dwarf_x = dwarf_points[:, 0]
    dwarf_y = dwarf_points[:, 1]
    dwarf_z = dwarf_points[:, 2]

    XY = compare(dwarf_x, dwarf_y, x, y)
    distance_xy = np.sqrt(XY[0]**2 + XY[1]**2)
    XZ = compare(dwarf_x, dwarf_z, x, z)
    distance_xz = np.sqrt(XZ[0]**2 + XZ[1]**2)
    YZ = compare(dwarf_y, dwarf_z, y, z)
    distance_yz = np.sqrt(YZ[0]**2 + YZ[1]**2)

    print(f"500 Myr, {map[num]}, XY: {distance_xy}, XZ: {distance_xz}, YZ: {distance_yz}")

    # 100 Myr
    dwarf_points = []

    for i, orbit in enumerate(dwarf_orbits3):
        final_pos = orbit[-1].pos.xyz.to(u.kpc).value[:3]
        dwarf_points.append(final_pos)

    dwarf_points = np.array(dwarf_points)

    dwarf_x = dwarf_points[:, 0]
    dwarf_y = dwarf_points[:, 1]
    dwarf_z = dwarf_points[:, 2]

    XY = compare(dwarf_x, dwarf_y, x, y)
    distance_xy = np.sqrt(XY[0]**2 + XY[1]**2)
    XZ = compare(dwarf_x, dwarf_z, x, z)
    distance_xz = np.sqrt(XZ[0]**2 + XZ[1]**2)
    YZ = compare(dwarf_y, dwarf_z, y, z)
    distance_yz = np.sqrt(YZ[0]**2 + YZ[1]**2)

    print(f"100 Myr, {map[num]}, XY: {distance_xy}, XZ: {distance_xz}, YZ: {distance_yz}")
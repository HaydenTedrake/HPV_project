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

def compare(x_vals1, y_vals1, x_vals2, y_vals2):
    pos1 = np.column_stack((x_vals1, y_vals1))
    mean1 = np.mean(pos1, axis=0)
    cov1 = np.cov(pos1, rowvar=False)

    pos2 = np.column_stack((x_vals2, y_vals2))
    mean2 = np.mean(pos2, axis=0)
    cov2 = np.cov(pos2, rowvar=False)

    dx = mean1[0] - mean2[0]
    dy = mean1[1] - mean2[1]

    combined_sigma_x = np.sqrt(np.max([cov1[0, 0],cov2[0, 0]]))
    combined_sigma_y = np.sqrt(np.max([cov1[1, 1],cov2[1, 1]]))

    Dx = dx / combined_sigma_x
    Dy = dy / combined_sigma_y

    return (Dx, Dy)

def final_points(dict, ics, num_steps):
    orbit = pot.integrate_orbit(ics, dt=-0.1*u.Myr, n_steps=num_steps)

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

    orbits = []
    orbits.append(orbit)
    for ics in ics_list:
        orbit = pot.integrate_orbit(ics, dt=-0.1*u.Myr, n_steps=num_steps)
        orbits.append(orbit)

    points = []
    for i, orbit in enumerate(orbits):
        final_pos = orbit[-1].pos.xyz.to(u.kpc).value[:3]
        points.append(final_pos)
    points = np.array(points)

    x = points[:, 0]
    y = points[:, 1]
    z = points[:, 2]

    return x, y, z

pot = gp.MilkyWayPotential()

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

map = {
    0: 'hvs5',
    1: 'hvs8',
    2: 'hvs14',
    3: 'hvs17',
    4: 'hvs23',
    5: hvs5,
    6: hvs8,
    7: hvs14,
    8: hvs17,
    9: hvs23,
}

# Leo I
dwarf_dict = {
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
# # Sextans (I)
# dwarf_dict = {
#     'ra': 153.262584, #degree
#     'e_ra': 0.0001, #degree
#     'dec': -1.6147, #degree
#     'e_dec': 0.0001, #degree
#     'dist': 95.0, #kpc
#     'e_dist': 3.0, #kpc
#     'pmra': -0.41, #mas/yr
#     'e_pmra': 0.01, #mas/yr
#     'pmdec': 0.04, #mas/yr
#     'e_pmdec': 0.01, #mas/yr
#     'rv': 224.3, #km/s
#     'e_rv': 0.1, #km/s
#     }
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

# 1 Gyr
dx1, dy1, dz1 = final_points(dwarf_dict, dwarf_ics, 10000)
# 500 Myr
dx2, dy2, dz2 = final_points(dwarf_dict, dwarf_ics, 5000)
# 100 Myr
dx3, dy3, dz3 = final_points(dwarf_dict, dwarf_ics, 1000)

result5 = []
result8 = []
result14 = []
result17 = []
result23 = []

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

    # 1 Gyr
    x1, y1, z1 = final_points(star, ics, 10000)
    # 500 Myr
    x2, y2, z2 = final_points(star, ics, 5000)
    # 100 Myr
    x3, y3, z3 = final_points(star, ics, 1000)

    # 1 Gyr
    XY = compare(dx1, dy1, x1, y1)
    distance_xy = np.sqrt(XY[0]**2 + XY[1]**2)
    XZ = compare(dx1, dz1, x1, z1)
    distance_xz = np.sqrt(XZ[0]**2 + XZ[1]**2)
    YZ = compare(dy1, dz1, y1, z1)
    distance_yz = np.sqrt(YZ[0]**2 + YZ[1]**2)
    total = np.sqrt(distance_xy**2 + distance_xz**2 + distance_yz**2)

    if num == 0:
        result5.append(total)
    if num == 1:
        result8.append(total)
    if num == 2:
        result14.append(total)
    if num == 3:
        result17.append(total)
    if num == 4:
        result23.append(total)
    print(f"1 Gyr, {map[num]}, {total}")

    # 500 Myr
    XY = compare(dx2, dy2, x2, y2)
    distance_xy = np.sqrt(XY[0]**2 + XY[1]**2)
    XZ = compare(dx2, dz2, x2, z2)
    distance_xz = np.sqrt(XZ[0]**2 + XZ[1]**2)
    YZ = compare(dy2, dz2, y2, z2)
    distance_yz = np.sqrt(YZ[0]**2 + YZ[1]**2)
    total = np.sqrt(distance_xy**2 + distance_xz**2 + distance_yz**2)

    print(f"500 Myr, {map[num]}, {total}")

    # 100 Myr
    XY = compare(dx3, dy3, x3, y3)
    distance_xy = np.sqrt(XY[0]**2 + XY[1]**2)
    XZ = compare(dx3, dz3, x3, z3)
    distance_xz = np.sqrt(XZ[0]**2 + XZ[1]**2)
    YZ = compare(dy3, dz3, y3, z3)
    distance_yz = np.sqrt(YZ[0]**2 + YZ[1]**2)
    total = np.sqrt(distance_xy**2 + distance_xz**2 + distance_yz**2)

    if num == 0:
        result5.append(total)
    if num == 1:
        result8.append(total)
    if num == 2:
        result14.append(total)
    if num == 3:
        result17.append(total)
    if num == 4:
        result23.append(total)
    print(f"100 Myr, {map[num]}, {total}")

print(result5)
print(result8)
print(result14)
print(result17)
print(result23)

for num, result in enumerate([result5, result8, result14, result17, result23]):
    num += 5
    star = map[num]
    if result[0]<result[1]:
        min = result[0]
        print(min)
        num_steps = 10000
        val = 0
        # need to go farther back in time
        while val < min:
            min = val
            num_steps += 1
            print(num_steps)
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

            dx, dy, dz = final_points(dwarf_dict, dwarf_ics, num_steps)
            x, y, z = final_points(star, ics, num_steps)

            XY = compare(dx, dy, x, y)
            distance_xy = np.sqrt(XY[0]**2 + XY[1]**2)
            XZ = compare(dx, dz, x, z)
            distance_xz = np.sqrt(XZ[0]**2 + XZ[1]**2)
            YZ = compare(dy, dz, y, z)
            distance_yz = np.sqrt(YZ[0]**2 + YZ[1]**2)
            val = np.sqrt(distance_xy**2 + distance_xz**2 + distance_yz**2)
            print(val)
        num -= 5
        print(f"{map[num]} is closest at {(num_steps-1)} steps with a distance of {val}")
    if result[0]>result[1]:
        min = result[1]
        num_steps = 1000
        val = 0
        # need to go farther back in time
        while val < min:
            min = val
            num_steps -= 1
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

            dx, dy, dz = final_points(dwarf_dict, dwarf_ics, num_steps)
            x, y, z = final_points(star, ics, num_steps)

            XY = compare(dx, dy, x, y)
            distance_xy = np.sqrt(XY[0]**2 + XY[1]**2)
            XZ = compare(dx, dz, x, z)
            distance_xz = np.sqrt(XZ[0]**2 + XZ[1]**2)
            YZ = compare(dy, dz, y, z)
            distance_yz = np.sqrt(YZ[0]**2 + YZ[1]**2)
            val = np.sqrt(distance_xy**2 + distance_xz**2 + distance_yz**2)
        num -= 5
        print(f"{map[num]} is closest at {(num_steps-1)} steps with a distance of {val}")
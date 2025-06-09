from astropy import units as u
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord
import math

### Simple Coordinate Transformations for Draco and Sagittarius ###

draco_icrs = SkyCoord('17h20m12s', '+57d54m55s', frame='icrs', distance=76 * u.kpc)
draco_galactic = draco_icrs.transform_to('galactic')
draco_galactocentric = draco_icrs.transform_to(coord.Galactocentric)

sagittarius_icrs = SkyCoord('18h55m20s', '-30d32m43s', frame='icrs', distance=20 * u.kpc)
sagittarius_galactic = sagittarius_icrs.transform_to('galactic')
sagittarius_galactocentric = sagittarius_icrs.transform_to(coord.Galactocentric)

# print("DRACO")
# print(draco_icrs)
# print(draco_galactic)
# print(draco_galactocentric)
# print("\n\n")
# print("SAGITTARIUS")
# print(sagittarius_icrs)
# print(sagittarius_galactic)
# print(sagittarius_galactocentric)

### Incorporate Velocities ###

hvs1 = coord.SkyCoord(
    ra=136.93746493225 * u.degree,
    dec=2.75190851109 * u.degree,
    distance=102.24 * u.kpc,
    pm_ra_cosdec=-0.60 * (u.mas / u.yr),
    pm_dec=-0.47 * (u.mas / u.yr),
    radial_velocity=0.37 * (u.km / u.s),
    frame="icrs"
)

### 1 mas/yr, 50 kpc. find velocity in km/s ###

radians = 1/(1000*3600)*math.pi/180
radianspersecond = radians * 1/(365*24*3600)

kilometers = 50 * 3.086e16

kilometerspersecond = radianspersecond * kilometers
print(kilometerspersecond)
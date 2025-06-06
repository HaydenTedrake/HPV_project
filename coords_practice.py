from astropy import units as u
from astropy.coordinates import SkyCoord

draco_icrs = SkyCoord('17h20m12s', '+57d54m55s', frame='icrs', distance=76 * u.kpc)
draco_galactic = draco_icrs.transform_to('galactic')
draco_galactocentric = draco_icrs.transform_to('galactocentric')

sagittarius_icrs = SkyCoord('18h55m20s', '-30d32m43s', frame='icrs', distance=20 * u.kpc)
sagittarius_galactic = sagittarius_icrs.transform_to('galactic')
sagittarius_galactocentric = sagittarius_icrs.transform_to('galactocentric')

print("DRACO")
print(draco_icrs)
print(draco_galactic)
print(draco_galactocentric)
print("\n\n")
print("SAGITTARIUS")
print(sagittarius_icrs)
print(sagittarius_galactic)
print(sagittarius_galactocentric)
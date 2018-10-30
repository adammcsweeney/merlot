# Module contains Astronomical Constants, and data on solar Solar System planets
# Classes currently contain following data: name, radius - r (m), mass - m (kg), semi-major axis of orbit - a (m)

# TODO ==== add sidereal rotation period, inclination of equator to orbit plane, orbit eccentricity, inclination of orbit to the ecliptic plane, orbit sidereal period

# 1 - PHYSICAL CONSTANTS
# universal gravitational constant (m3 kg-1 s-2)
G = 6.674e-11


# 2 - SOLAR SYSTEM DATA

class Sun():
    name = 'Sun'
    r = 696000e3
    m = 1.989e30

class Mercury():
    name = 'Mercury'
    r = 2440e3
    m = 330.2e21
    a = 57.91e9

class Venus():
    name = 'Venus'
    r = 6052e3
    m = 4.869e24
    a = 108.2e9

class Earth():
    name = 'Earth'
    r = 6878e3
    m = 5.974e24
    a = 149.6e9

class Moon():
    name = 'Moon'
    r = 1737e3
    m = 73.48e21
    a = 384.4e6

class Mars():
    name = 'Mars'
    r = 3396e3
    m = 641.9e21
    a = 227.9e9

class Jupiter():
    name = 'Jupiter'
    r = 71490e3
    m = 1.899e27
    a = 778.6e9

class Saturn():
    name = 'Saturn'
    r = 60270e3
    m = 568.5e24
    a = 1.433e12

class Uranus():
    name = 'Uranus'
    r = 25560e3
    m = 86.83e24
    a = 2.872e12

class Neptune():
    name = 'Neptune'
    r = 24760e3
    m = 102.4e24
    a = 4.495e12

class Pluto():
    name = 'Pluto'
    r = 1195e3
    m = 12.5e21
    a = 5.87e12










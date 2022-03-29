import numpy as np


## Physical parameters & constants in cgs unit
EARTH_MASS = np.float64(5.972e27)       # Earth mass in g
SOLAR_MASS = np.float64(1.989e33)       # Solar mass in g
JUPITER_MASS = np.float64(1.898e30)     # Jupiter mass in g

## Unit converting factor
AU_TO_CM = np.float64(1.496e13)         # from AU to cm
CM_TO_AU = np.float64(6.684e-14)        # from cm to AU
YEAR_TO_SEC = np.float64(3.1536e7)      # from year to sec
SEC_TO_YEAR = np.float64(3.171e-8)      # from sec to year

DISTANCE_TO_EARTH = AU_TO_CM            # Distance to earth in cm
DISTANCE_TO_JUPITER = 5.2*AU_TO_CM      # Distance to jupiter in cm 

G = np.float64(6.67e-8)                 # Gravitational constant in cgs
GG = G*SOLAR_MASS*YEAR_TO_SEC**2./AU_TO_CM**3 # Normalized G

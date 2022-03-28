import numpy as np


## Physical parameters & constants in cgs unit
EARTH_MASS = np.float64(5.972e27)       # Earth mass in cgs
SOLAR_MASS = np.float64(1.989e33)       # Solar mass in cgs
JUPITER_MASS = np.float64(1.898e30)     # Jupiter mass in cgs

G = np.float64(6.67e-8)                 # Gravitational constant in cgs

## Unit converting factor
AU_TO_CM = np.float64(1.496e13)         # from AU to cm
CM_TO_AU = np.float64(6.684e-14)        # from cm to AU
YEAR_TO_SEC = np.float64(3.1536)        # from year to sec
SEC_TO_YEAR = np.float64(3.171e-8)      # from sec to year

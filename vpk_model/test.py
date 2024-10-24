import vpk_model as vpk
import numpy as np
# Example usage
del2 = -0.5
xb = 0.2
Q2 = 2.0
phi_g = np.pi / 4
E = 10.6
heli = 1
mesonmass = 0.1349766

result = vpk.dvmpx_vectorized(del2, xb, Q2, phi_g, E, heli, mesonmass)
print(result)

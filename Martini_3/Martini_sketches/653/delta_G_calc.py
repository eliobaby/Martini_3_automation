import numpy as np

# your data
water   = np.array([ 1.38, 28.29, 37.89])
octanol = np.array([15.59, 28.40, 23.88])

# constants
R = 8.314          # J/(mol·K)
T = 298.15         # K

# 1) partition coefficients
P = octanol / water

# 2) ΔG in J/mol
deltaG_J = -R * T * np.log(P)

# 3) ΔG in kJ/mol
deltaG_kJ = deltaG_J / 1000

print(deltaG_kJ)
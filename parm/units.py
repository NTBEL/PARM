# Avogadro's Number from scipy
from scipy.constants import N_A

# Conversion factors for concentration units.
# microMolar to molecule/pL
microM_to_molec_per_pL = 1e-6 * N_A * 1e-12
# molecule/pL to microM
molec_per_pL_to_microM = 1 / microM_to_molec_per_pL

# nanoMolar to molecule/pL
nM_to_molec_per_pL = 1e-9 * N_A * 1e-12
# molecule/pL to nanoMolar
molec_per_pL_to_nM = 1 / nM_to_molec_per_pL

# Cubic micron to picoliter:
#               cubic-micron to mL * mL to L * L to pL
cubicmicron_to_pL = (1e-4) ** 3 * 1e-3 * 1e12

# milliliter to picoliter
mL_to_pL = 1e-3 * 1e12
# picoliter to milliliter
pL_to_mL = 1 / mL_to_pL

# liter to picoliter
L_to_pL = 1e12

# picoliter to liter
pL_to_L = 1e-12

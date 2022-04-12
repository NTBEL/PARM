from . import units

# Default forward, reverse, and catalytic rates:
# Default forward binding rate (for cell volume of 1 pL) from
# Aldridge et al. https://doi.org/10.1038/ncb1497
# Note that KF_BIND is equivalent to kf/(1e-12) / N_A for kf in 1/(M*s).
KF_BIND = 1e-6  # 1/(molecule/pL*s)
# Default dissociation rate from Albeck et al.
# https://doi.org/10.1371/journal.pbio.0060299
# is 1e-3 1/s.
# Note that the range of plausible dissociation rates suggested by Aldridgle et
# al. is 0.001 - 0.01 1/s
KR_BIND = 1e-3  # 1/s for concentrations in number/pL
# kcat for the "average enzyme" from
# Bar-Even et al. https://doi.org/10.1021/bi2002289 is 10 1/s.
# Note the broader range of plausible catalytic rates suggested by Aldridgle et al.
# is 1 - 10 1/s.
KCAT = 10  # 1/s

# Default signaling protein concentration range
# of 1 nM to 1 microM range assumed by
# Albeck et al. https://doi.org/10.1371/journal.pbio.0060299
# for which they reference
# Wu and Pollard https://doi.org/10.1126/science.1113230
# Also see Aldridge et al. https://doi.org/10.1038/ncb1497
# Set to 100 nM or 0.1 microM
SPC = 0.1 * units.microM_to_molec_per_pL

# Ion channel transport rate: up to 1e8 ions/s https://www.ncbi.nlm.nih.gov/books/NBK26910/
K_ION_CHANNEL = 1e8

# Default molecule degradation rate.
K_DEGRADE = 1  # 1/s

# Default unidirection conversion rate.
K_CONVERT = 1  # 1/s

# Default cellular volume of 1 pL (or 10^-12 L) from
# Aldridge et al. https://doi.org/10.1038/ncb1497
# Also assumed by Albeck et al. https://doi.org/10.1371/journal.pbio.0060299
# for their apoptosis model.
VOL_CELL = 1  # pL

# Default membrane thickness for effective volume of membrane compartments.
# Effective cell-membrane thickness
# Assume 10 nm (0.01 micron) as in
#  https://github.com/RuleWorld/BNGTutorial/blob/master/CBNGL/LR_comp.bngl
MEMBRANE_THICKNESS = 0.01

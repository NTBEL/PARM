# PySB components
from pysb import (
    Model,
)

from parm import compartments, receptor_modules


Model()

# Define the comparments and their volumes.
compartments.hek293_cell()
# Receptor monomer.
receptor_modules.receptor_monomer_only()
# Initials for the receptor.
receptor_modules.receptor_initials()
# 1st-order degradation of the resting inactive receptor.
receptor_modules.inactive_par2_degradation()
# zero-order synthesis of PAR2.
receptor_modules.par2_synthesis()

# observables.
receptor_modules.observables()

# PySB components
from pysb import (
    Model,
)

from parm import compartments, receptor_modules, gprotein_modules, calcium_modules


Model()

# Define the comparments and their volumes.
compartments.hek293_cell()

# PAR2 activation using single-state model of activation:
#   2AT + PAR2_I <---> TAT:PAR2_A
receptor_modules.single_state_par2_activation()

# Include the observables.
# Receptor occupation, etc.
receptor_modules.observables()

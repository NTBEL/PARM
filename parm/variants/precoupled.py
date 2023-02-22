# PySB components
from pysb import (
    Model,
)

from parm import receptor_modules, gprotein_modules

# Import the base parm model.
from parm import parm

Model(base=parm.model)
# Add the new mechanistic parts for precoupling:
# PAR_I + Gaq:GDP:Gbg <---> PAR2_I:Gaq:GTP:Gbg
gprotein_modules.heterotrimer_precouples_free_inactive_par2()
# 2AT + PAR_I:Gaq:GDP:Gbg ---> 2AT:PAR2_A + Gaq:GTP + Gbg
receptor_modules.addon_single_state_precoupled_par2_activation_catalytic()
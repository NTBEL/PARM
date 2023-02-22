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
# Receptor degradation:
#   TAT:PAR2_A ---> None
receptor_modules.occupied_par2_degradation()

# Gprotein activation and regulation using the heterotrimeric G protein cycle
# mechanism:
#  i) G protein activation:
#   2AT:PAR2_A + Gaq:GDP:Gbg ---> 2AT:PAR2_A + Gaq:GTP + Gbg
#  ii) GTP hydrolosis:
#       a) (slow) hydrolosis by Gaq alone:
#           Gaq:GTP ---> Gaq:GDP
#       b) RGS enhanced hydrolosis:
#           Gaq:GTP --RGS--> Gaq:GDP
#  iii) heterotrimer reassociation:
#   Gaq:GDP + Gbg ---> Gaq:GDP:Gbq
gprotein_modules.yi2003_heterotrimeric_gprotein_cycle()
# PAR_I + Gaq:GDP:Gbg <---> PAR2_I:Gaq:GTP:Gbg
gprotein_modules.heterotrimer_precouples_free_inactive_par2()
# 2AT + PAR_I:Gaq:GDP:Gbg ---> 2AT:PAR2_A + Gaq:GTP + Gbg
receptor_modules.addon_single_state_precoupled_par2_activation_catalytic()

# Calcium signaling via PLC and IP3 that leads to release of the ER calcium
# store.
calcium_modules.gaq_activated_calcium_signaling_simplified()
# PLC enhances Gaq's hydrolosis of GTP.
gprotein_modules.addon_plc_enhances_gaq_hydrolosis_of_gtp_to_gdp()

# Include the observables.
# Receptor occupation, etc.
receptor_modules.observables()
# Active Gaq and ratio.
gprotein_modules.observables()
# Cytosolic calcium levels, etc.
calcium_modules.observables()
# Parameters and Expressions for the TN-XXL FRET sensor.
calcium_modules.fret_calcium_indicator_tnxxl()
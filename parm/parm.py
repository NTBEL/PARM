# PySB components
from pysb import (
    Model,
)

from . import compartments, receptor_modules, gprotein_modules, calcium_modules


Model()

# Define the comparments and their volumes.
compartments.hek293_cell()
# compartments.default_cell()
# PAR2 activation.
# Minimal two-state model of activation.
# 2AT + PAR2_I <---> TAT:PAR2_I <---> TAT:PAR2_A
receptor_modules.minimal_two_state_par2_activation()
# Receptor degradation:
#  i) bound and inactive PAR2:
#     TAT:PAR2_I ---> None
# ii) bound and active PAR2:
#     TAT:PAR2_A ---> None
receptor_modules.occupied_par2_degradation()
# Gprotein activation:
# classic, no precoupling, mechanism.
# i) G protein heterotrimer binds activated PAR2:
#  PAR2_A + Gaq:GDP:Gbg <---> PAR2_A:Gaq:GDP:Gbg
# ii) GDP preferentially unbinds from Gaq:
#  PAR2_A:Gaq:GDP:Gbc <---> PAR2_A:Gaq:Gbc + GDP
# iii) GTP preferentially binds Gaq:
#  PAR2_A:Gaq:Gbc + GTP <---> PAR2_A:Gaq_A:GTP:Gbg
# iv) Gbg dissociates from Gaq (i.e., heterotrimer dissociation):
#  PAR2_A:Gaq:GTP:Gbc <---> PAR2_A:Gaq:GTP + Gbc
# v) Gaq:GTP dissociates from PAR2, Gaq is now active (G protein dissociation from the receptor):
#  PAR2_A:Gaq:GTP <---> PAR2_A + Gaq:GTP
gprotein_modules.classic_activation_mechanism()
# Grpotein signal regulation
# a) Slow hydrolosis by Gaq alone
#    Gaq:GTP ---> Gaq:GDP
# b) RGS enhanced hydrolosis
#    Gaq:GTP + RGS <---> Gaq:GTP:RGS ---> Gaq:GDP + RGS
gprotein_modules.gaq_hydrolyzes_gtp_to_gdp()
gprotein_modules.rgs_enhances_gaq_hydrolosis_of_gtp_to_gdp()
# Calcium signaling via PLC and IP3 that leads to release of the ER calcium
# store.
# calcium_modules.gaq_activated_calcium_signaling()
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

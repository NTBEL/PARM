"""PARM: PAR2 Activation-driven calcium Release Model, classic, noncompetitive

Variant of parm that incorporates noncompetitive antagonism.

"""

from . import antagonist_modules

# Import the base parm model.
from parm.parm import model

# Add the noncompetitive antognist.
antagonist_modules.noncompetitive_par2_antagonist()

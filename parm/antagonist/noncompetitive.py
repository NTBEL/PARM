"""PARM: PAR2 Activation-driven calcium Release Model, classic, noncompetitive

Variant of parm that incorporates noncompetitive antagonism.

"""
# PySB components
from pysb import (
    Model,
)
from . import antagonist_modules

# Import the base parm model.
from parm import parm

Model(base=parm.model)
# Add the noncompetitive antogonist.
antagonist_modules.noncompetitive_par2_antagonist()

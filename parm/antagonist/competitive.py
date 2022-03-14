"""PARM: PAR2 Activation-driven calcium Release Model, classic, competitive

Variant of parm that incorporates competitive antagonism.

"""
# PySB components
from pysb import (
    Model,
)

from . import antagonist_modules

# Import the base parm model.
from parm import parm

Model(base=parm.model)
# Add the competitive antognist.
antagonist_modules.competitive_par2_antagonist()

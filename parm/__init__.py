import numpy as np

# Import the main model.
from .parm import model

# Make a vector of the default model parameters.
default_param_values = [param.value for param in model.parameters]

# Create a dictionary of pamameter masks for fancy boolean indexing of the
# the parameter vector.
parameter_masks = dict()
for param in model.parameters:
    parameter_masks[param.name] = [par.name == param.name for par in model.parameters]

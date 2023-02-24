# NumPy, Matplotlib, and seaborn
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
# Seaborn settings
sns.set(rc={"figure.dpi": 100, "savefig.dpi": 300})
sns.set_palette("bright")
sns.set_context("notebook")
sns.set_style("ticks")
# PySB simulator
from pysb.simulator import ScipyOdeSimulator  # http://pysb.org/download
# The model
from parm.variants.LR import (
    model,
)  # path_to/PARM needs to be added to your PYTHONPATH
# The units module
from parm import units


tspan = np.arange(0, 100)

# PySB solver options.
integrator_options = {"rtol": 1e-6, "atol": 1e-6}
solver_args = {"integrator": "lsoda", "integrator_options": integrator_options}
solver = ScipyOdeSimulator(model, tspan=tspan, **solver_args)

# Get the nominal parameter values.
param_values = np.array([param.value for param in model.parameters])


# Mask for initial concentration of 2AT
twoat_mask = [param.name == "TAT_0" for param in model.parameters]
# Mask for the resting cytosolic Ca2+ amount for the FRET estimation.
kd_mask = [param.name == "Kd_PAR2_bind_TAT" for param in model.parameters]
# Mask for the initial amount of PAR2.
par2_mask = [param.name == "PAR2_0" for param in model.parameters]
# Mask for the resting cytosolic Ca2+ amount for the FRET estimation.
kf_mask = [param.name == "kf_PAR2_bind_TAT" for param in model.parameters]

tat_concs = np.array([10, 31.6, 100, 110, 316, 1000, 3160])
Kd = 110  # nM

v_extra = model.parameters["Vextra"].value
param_values[kd_mask] = Kd * units.nM_to_molec_per_pL * v_extra
#param_values[kf_mask] *= 10

for tat_conc in tat_concs:
    tat_num = tat_conc * units.nM_to_molec_per_pL * v_extra
    param_values[twoat_mask] = tat_num
    sim = solver.run(param_values=param_values, initials=None).all
    # Get the simulated FRET response
    y_sim = sim["par2_occupancy"]
    plt.plot(tspan, y_sim, label=str(tat_conc))
plt.hlines(0.5, 0, 100, linestyle="--")
plt.legend()
plt.show()
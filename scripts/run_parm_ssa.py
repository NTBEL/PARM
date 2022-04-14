import numpy as np
from pysb.simulator.bng import BngSimulator
import matplotlib.pyplot as plt
from parm.parm import model
from parm import units

tspan = np.linspace(0, 300, num=600, endpoint=True)
simulator = BngSimulator(model)
param_values = np.array([parm.value for parm in model.parameters])
# Mask for initial concentration of 2AT
twoat_mask = [parm.name == "TAT_0" for parm in model.parameters]
param_values[twoat_mask] = (
    316 * units.nM_to_molec_per_pL * model.parameters["Vextra"].value
)
print("running PARM...")
sim_run = simulator.run(tspan=tspan, param_values=param_values, n_runs=5, method="ssa")
tout = sim_run.tout
yout = np.array(sim_run.observables)
print("Done running...")

plt.plot(tout.T, yout["FRET"].T, label="FRET")
plt.ylabel("FRET Ratio")
plt.xlabel("Time (s)")
plt.show()

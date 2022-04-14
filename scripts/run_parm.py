import numpy as np
from pysb.simulator import ScipyOdeSimulator
import matplotlib.pyplot as plt
from parm.parm import model
from parm import units

tspan = np.linspace(0, 300, num=600, endpoint=True)
solver = ScipyOdeSimulator(model, tspan=tspan, integrator="lsoda")
param_values = np.array([parm.value for parm in model.parameters])
# Mask for initial concentration of 2AT
twoat_mask = [parm.name == "TAT_0" for parm in model.parameters]
param_values[twoat_mask] = (
    316 * units.nM_to_molec_per_pL * model.parameters["Vextra"].value
)
print("running PARM...")
m_run = solver.run(param_values=param_values)
yout = m_run.all
print("Done running...")

plt.plot(tspan, yout["FRET"], label="FRET")
plt.ylabel("FRET Ratio")
plt.xlabel("Time (s)")
plt.show()

import numpy as np
import matplotlib.pyplot as plt

from parm import model, default_param_values
from parm import util

tspan = np.linspace(0, 300, num=600, endpoint=True)

# Input the inital 2AT concentration as 316.0 nM.
param_values = util.set_tat_initial_nM(default_param_values, 316.0)

print("running PARM...")
traj_out = util.run_model(model, tspan, param_values=param_values)
print("Done running...")

plt.plot(tspan, traj_out["FRET"], label="FRET")
plt.ylabel("FRET Ratio")
plt.xlabel("Time (s)")
plt.show()

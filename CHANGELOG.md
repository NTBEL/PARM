## [branch: setup-script-and-run-function] - 2022-10-25

### Changed
* Updated the `parm.util.run_model` function to accept an option `param_values` input for parameter vectors to be passed on to the `ScipyOdeSimulator`. Also added type hints and a docstring to the the function.

### Added
* A `setup.py` setuptools install script.
* An `environment.yml` conda environment file that includes pysb and cython.

### Fixed
* A variable misspelling error in the `parm.util.run_model` function for the return variable `yout`.


## [branch: reorganize_modularize] - 2022-03-04

### Changed
* Major reorganization of the model code into modular functions (PySB modules) and definition of models using modular functions.
* Reduced the models down to a single core parm.py model (i.e., the main model) and any encoded variations of the model are now defined in the `variants` module.
* antogonist models (in `parm.antagonist`) were updated to reflect changes in the core model and to also use the modular design.

### Added
* Module functions for additional mechanistic elements including calcium feedback in the cytosol (IP3R inhibition, PLC enhancement), calcium buffering, and additional reactions to control calcium homeostasis. These mechanistic elements were added to the core parm model.
* With restructuring of model the following new modules were created: `comparments`, `receptor_modules`, `gprotein_modules`, `calcium_modules`, `defaults`, `units`, and `parm`. Additionally, in `parm.antagonist` the new `antagonist_modules` module was added.
* Addition of the `util` module with functions for running the model, pre-equlibration, and expanding the time points in a previously generated (or experimental) time series.


## [branch: new_reactions] - 2022-02-17

### Added
  * classic_2m model: New reactions for degradation of PAR2 (3 degrade macro calls with new parameters `k_PAR2_bound_degradation` and `k_PAR2_denatured_degradation`).
  * classic_2m model: Explicit setting of the Initial for cytosolic calcium (Initial Ca_C_0).
  * classic_2m model: Define and set extracellular concentration of calcium (Initial Ca_extra_0).
  * classic_2m model: New reaction for influx of extracellular calcium into the cytosol (Rule Ca_extra_to_cyt with Parameter k_Ca_extra_to_cyt).
  * classic_2m model: expression named `occupancy`.

### Changed
  * classic_2m model: Changed the forward rate constant for 2AT binding to PAR2 from the default KF_BIND to be 3*KF_BIND which matches the ligand-receptor binding forward rate constant in the  G-protein cycle model of  Yi et al. 2003 PNAS https://doi.org/10.1073/pnas.1834247100.
  * classic_2m model: the Rule for degradation of cytosolic calcium was replaced with a first order conversion of cytosolic calcium to extracellular calcium.
  * classic_2m model: occupancy_ratio expression was changed to inactive_bound_ratio.

### Fixed
 * A mistake in the equation given in the comment about the default KF_BIND value. Updated to `KF_BIND is equivalent to kf/(Vcell*1e-12) / N_A for kf in 1/(M*s)`

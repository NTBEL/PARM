## v0.3.0 - 2023-02-09

### Changed
* `util.run_model` has an input parameter `nprocs` that is passed to the solver's `run` function as `num_processors` to allow for parallel execution of the inputs if multiple parameter sets or initial conditions are provided. 

### Added
* `util.load_pydream_chains` function that can load the sampled parameter vectors output from PyDREAM sampling chains.
* `visualize.display_expcomp_single` function that takes in a single parameter vector, simulates the model at various 2AT concentrations, and plots a comparison with the experimental FRET ratio data.
* `visualize.display_expcomp_multi_grid_mean_ci` function that takes in a list/array of parameter vectors, simulates the model at various 2AT concentrations for each parameter vector, and plots a comparison with the experimental FRET ratio data with a 2x3 grid of plots for each of the 6 2AT concentrations. The function plots the mean and standard deviation over all the input parameter vectors for a given 2AT concentration.

## v0.2.0 - 2023-02-05

Note: merging in from branch `model-variant/extrusion-inlflux-single`

### Changed
* Moved the `parm.parm` model to `parm.variants.classic`. This model also uses `gaq_activated_calcium_signaling_simplified` now instead of `gaq_activated_calcium_signaling`.
* Moved the `parm.variants.heterogprot_cycle` model to be `parm.parm`, so that this model is now the main model. This model also uses `gaq_activated_calcium_signaling_simplified` now instead of `gaq_activated_calcium_signaling`.
* The `calcium_modules.gaq_activated_calcium_signaling_simplified` module-function was updated to use the `calcium_extrusion_and_influx_single` along with the `cytosolic_calcium_positive_feedback`.
* The `util.expand_times` function no longer uses `np.linspace` but instead just uses `np.arange` to make time points at 1 second intervals. This works with the experimental time points because they are all whole number times.


### Added
* Created new `calcium_modules.calcium_extrusion_and_influx_single` module-function that defines a single 1st order reaction to pump excess Ca2+ from the cytosol.
* Created new `calcium_modules.cytosolic_calcium_positive_feedback` module-function that just defines the positive feedback of Ca2+ on PLC catalysis of PIP2 to IP3.  
* New `parm.variants.LR` model variant that just defines ligand-receptor binding using the `receptor_modules.single_state_par2_activation` module-function. This model can be used to test the binding/dissociation.

## v0.1.0 - 2022-10-27

### Changed
* Starting semantic versioning of the model code at v0.1.0.
* Made the `parm.util.pre_equilibration` function private by renaming to `parm.util._pre_equilibration`.
* The main PARM model `parm.parm` is now imported in the package `__init__.py` so that it can be imported like `from parm import model`.
* Moved the pysb dependency in setup.py to the `install_requires` list.

### Added
* Additional sections (Install, License, etc.) and information in the README, as well as shields/badges towards the top of the README.
* New wrapper function for pre-equilibration in the `util` module, `parm.util.pre_equilibrate` that calls the `_pre_equilibration` function but also does addtional steps to automate the process of setting 2AT to zero concentration and adjusting the baseline cytosolic calcium concentration used for the FRET ratio. This function also has a simpler output that contains an updated parameter vector and the initial concentrations that can be directly passed to the `run_model` function.
* New utility function `parm.util.calcium_homeostasis_reverse_rate_coupling` that adjusts a parameter vector to account for coupling of the reverse rates of the calcium homeostasis reactions to the forward rates. This function can be used adjust the parameters when sampling values of the reverse rate parameters for those reactions: for example, during model training.
* New default parameter vector, `default_param_values` created in the package `__init__.py` for the main `parm.parm` model which can be imported like `from parm import default_param_values`.
* New dictionary of parameter vector masks for fancy indexing, `parameter_masks`, created for the main `parm.parm` model which can imported like `from parm import parameter_masks`.
* Utility function `parm.util.get_parameter_vector` that takes a PySB model and returns a vector of the parameter values.
* Utiltiy function `parm.util.get_parameter_mask` that takes a PySB model and parameter name and returns a boolean mask to fancy index the corresponding paramter value vector for that particular parameter.
* Utility function `parm.util.set_tat_initial_nM` that takes a parameter value vector and an initial concentration of 2AT in nM and updates the parameter vector with the desired 2AT concentration.

### Fixed
* Converted the `extras_require` argument in the setup.py to be dictionary instead of a list, now with just the optional Cython dependency. Note that the pysb dependency was moved into the `install_requires` list.
* Commented out the `long_description` in setup.py because it was causing error when trying to read the README file.  


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


## [branch: new_reactions] - 2022-02-17

### Added
  * classic_2m model: New reactions for synthesis and degradation of PAR2 (1 synthesize and 4 degrade macro calls with new Parameter k_PAR2_synthesis, k_PAR2_I_degradation, k_PAR2_bound_degradation, and k_PAR2_denatured_degradation).
  * classic_2m model: Explicit setting of the Initial for cytosolic calcium (Initial Ca_C_0).
  * classic_2m model: Define and set extracellular concentration of calcium (Initial Ca_extra_0).
  * classic_2m model: New reaction for influx of extracellular calcium into the cytosol (Rule Ca_extra_to_cyt with Parameter k_Ca_extra_to_cyt).

### Changed
  * classic_2m model: Changed the forward rate constant for 2AT binding to PAR2 from the default KF_BIND to be 3*KF_BIND which matches the ligand-receptor binding forward rate constant in the  G-protein cycle model of  Yi et al. 2003 PNAS https://doi.org/10.1073/pnas.1834247100.
  * classic_2m model: the Rule for degradation of cytosolic calcium was replaced with a first order conversion of cytosolic calcium to extracellular calcium.

### Fixed
 * A mistake in the equation given in the comment about the default KF_BIND value. Updated to `KF_BIND is equivalent to kf/(Vcell*1e-12) / N_A for kf in 1/(M*s)`

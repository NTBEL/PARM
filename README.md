# PARM


![Python version badge](https://img.shields.io/badge/python-3.8-blue.svg)
[![PySB version badge](https://img.shields.io/badge/PySB->%3D1.13.2-9cf.svg)](https://pysb.org/)
[![license](https://img.shields.io/github/license/NTBEL/PARM.svg)](LICENSE)
![version](https://img.shields.io/badge/version-0.1.0-orange.svg)
[![release](https://img.shields.io/github/release-pre/NTBEL/PARM.svg)](https://github.com/NTBEL/PARM/releases/tag/v0.1.0)

**P**AR2 **A**ctivation and calcium signaling **R**eaction **M**odel (PARM)

PARM contains rules-based models of PAR2 (proteinase-activated receptor isoform 2) activation, G-protein activation, and calcium signaling via the phospholipase C and IP3 (inositol triphosphate) pathway. Models are encoded as Python program modules using the [PySB](http://pysb.org/) modeling framework. The models are designed to mathematically model the underlying signaling dynamics relevant to the inactivation of PAR2 in HEK-293 cells by Molecular Hyperthermia as described in:

  Kang et al.,  Transient Photoinactivation of Cell Membrane Protein Activity without Genetic Modification by Molecular Hyperthermia, ACS Nano 2019, 13, 11, 12487–12499 [https://doi.org/10.1021/acsnano.9b01993](https://doi.org/10.1021/acsnano.9b01993)

## Table of Contents

  1. [Install](#install)
    1. [Dependencies](#dependencies)
    2. [pip install](#pip-install)
    3. [Manual install](#manual-install)
    4. [Recommended additional software](#recommended-additional-software)
  2. [Documentation and Usage](#documentation-and-usage)
    1. [Models in PARM](#models-in-parm)
    2. [Example usage](#example-usage)    
  3. [License](#license)
  4. [Change Log](#change-log)
  5. [Contact](#contact)
  6. [Citing](#citing)  

  ------

# Install

**PARM** installs as the `parm` Python package. It is tested with Python 3.8.

### Dependencies

Note that `parm` has the following core dependency:
   * [PySB](https://pysb.org/) >= 1.13.2

### pip install

First, install [PySB](https://pysb.org/download).

You can then install `parm` version 0.1.0 with `pip` sourced from the GitHub repo:

##### with git installed:

Fresh install:
```
pip install git+https://github.com/NTBEL/PARM@v0.1.0
```
Or to upgrade from an older version:
```
pip install --upgrade git+https://github.com/NTBEL/PARM@v0.1.0
```

##### without git installed:

Fresh install:
```
pip install https://github.com/NTBEL/diffusion-fit/archive/refs/tags/v0.1.0.zip
```
Or to upgrade from an older version:
```
pip install --upgrade https://github.com/NTBEL/diffusion-fit/archive/refs/tags/v0.1.0.zip
```

### Manual install

First, install [PySB](https://pysb.org/download).
Then, download the repository and from the `PARM` folder/directory run
```
pip install .
```

### Recommended additional software

The following software is not required for the basic operation of **parm**, but provide extra capabilities and features when installed.

#### Cython

[Cython](https://cython.org/) is used by PySB to compile the ODE reactions on-the-fly, which can greatly improve model performance when running with the `ScipyOdeSimulator`.

pip:
```
pip install Cython
```
conda:
```
conda install cython
```
------    

# Documentation and Usage

## Models in PARM

The core model of PARM is defined in `parm.parm` and can be imported at the package level like `from parm import model`.

 Additionally, PARM contains 2 extensions of the `parm.parm` model which incorporate an antagonist:

  * `parm.antagonist.competitive` - Adds a competitive antagonist.
  * `parm.antagonist.noncompetitive` - Adds a noncompetitive antagonist which operates via negative allosteric modulation of the agonist binding affinity. (Note: The factor
    which controls the allosteric modulation could also be set such that the antagonist induces positive allosteric modulation, increasing agonist binding affinity.)

There are also 3 models with mechanistic variations defined in `parm.variants`:
  * `parm.variants.precoupled` - Adds pre-coupling between PAR2 and the G-protein heterotrimer such that some PAR2 can bind to the heterotrimer under resting conditions (without any agonist).
  * `parm.variants.heterogprot_cycle` - The receptor binding and G-protein interaction mechanism is adapted from the yeast G-protein cycle model of [Yi et al. 2003](https://doi.org/10.1073/pnas.1834247100).
  * `parm.variants.par2_synthesis_degradation` - This model only contains PAR2 with reactions for its resting synthesis and degradation.

### Example usage

```
from parm import model
from parm.util import run_model

tspan = list(range(0, 10, 1))

traj_out = run_model(model, tspan)
```

#### Note on pre-equlibration

The main parm model includes some calcium homeostasis reactions that may require pre-eqilibration before running the actual simulation. This affects the calcium concentrations in different compartments and can affect the estimate of the FRET ratio. The model can pre-equilibrated using the `parm.util.pre_equilibrate` function. Here is an example:

```
from parm import model
from parm import util
import numpy as np
from pysb.simulator import ScipyOdeSimulator

# set the time span for pre-eqilibration.
tspan_pre = list(range(0, 3000, 1))
# Run the pre-equlibration.
param_values_eq, initials_eq = util.pre_equilibrate(model, tspan_pre)

# set the time span for the simulation.
tspan = list(range(0, 300, 1))                                                      
# Setup the PySB solver/simulator.
solver = ScipyOdeSimulator(model, tspan=tspan, integrator='lsoda')
# Run the simulation.
sim = solver.run(param_values=param_values_eq, initials=initial_eq)

```

------

# License

This project is licensed under the MIT License - see the [LICENSE](LICENSE.md) file for details

------

# Change Log

See: [CHANGELOG](CHANGELOG.md)

------

# Contact

Please open a [GitHub Issue](https://github.com/NTBEL/PARM/issues) to
report any problems/bugs or make any comments, suggestions, or feature requests for PARM.

If you need assistance with PySB-specific issues then you can also try the pysb gitter channel: https://gitter.im/pysb/pysb

------

# Citing

If this model or other package features are useful in your research and you wish to cite it, you can use the following software citation:
> B. A. Wilson, “PARM: PAR2 Activation and calcium signaling Reaction Model” (v0.1.0), https://github.com/NTBEL/PARM, 2022.

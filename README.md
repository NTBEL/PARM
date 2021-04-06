# PARM
**P**AR2 **A**ctivation and calcium signaling **R**eaction **M**odel

PARM contains rules-based models of PAR2 (proteinase-activated receptor isoform 2) activation and calcium signaling via the phospholipase C and IP3 (inositol triphosphate) pathway. Models are encoded as Python program modules using the [PySB](http://pysb.org/) modeling framework. The models are designed to mathematically model the underlying signaling dynamics relevant to the inactivation of PAR2 in HEK-293 cells by Molecular Hyperthermia as described in:

  Kang et al.,  Transient Photoinactivation of Cell Membrane Protein Activity without Genetic Modification by Molecular Hyperthermia, ACS Nano 2019, 13, 11, 12487–12499 [https://doi.org/10.1021/acsnano.9b01993](https://doi.org/10.1021/acsnano.9b01993)


## Models in PARM
PARM currently contains 2 model variants which encode alternative hypotheses for the GPCR/G-protein activation mechanism:

  * `parm.classic` - Classic GPCR/G-protein activation mechanism: G-protein heterotrimers only interact with the receptor after receptor-activation.
  * `parm.precoupled` - Pre-coupled GPCR/G-protein activation mechanism: A pre-defined fraction of receptors are initialized with G-protein heterotrimers pre-coupled to the receptor.

We refer you to Fig 2 of Oliveira et al. [https://doi.org/10.3389/fnagi.2019.00089](https://doi.org/10.3389/fnagi.2019.00089) for a shematic depicting the main elements of these two GPCR/G-protein activation mechanisms.

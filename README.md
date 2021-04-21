# PARM
**P**AR2 **A**ctivation and calcium signaling **R**eaction **M**odel

PARM contains rules-based models of PAR2 (proteinase-activated receptor isoform 2) activation and calcium signaling via the phospholipase C and IP3 (inositol triphosphate) pathway. Models are encoded as Python program modules using the [PySB](http://pysb.org/) modeling framework. The models are designed to mathematically model the underlying signaling dynamics relevant to the inactivation of PAR2 in HEK-293 cells by Molecular Hyperthermia as described in:

  Kang et al.,  Transient Photoinactivation of Cell Membrane Protein Activity without Genetic Modification by Molecular Hyperthermia, ACS Nano 2019, 13, 11, 12487–12499 [https://doi.org/10.1021/acsnano.9b01993](https://doi.org/10.1021/acsnano.9b01993)


## Models in PARM
PARM currently contains 6 model variants which encode alternative hypotheses for PAR2 activation and the GPCR/G-protein activation mechanisms:

  * `parm.classic_1` - Single-state PAR2 activation (L + R <---> LR*) and classic GPCR/G-protein activation mechanism: G-protein heterotrimers only interact with the receptor after receptor-activation.
  * `parm.classic_2m` - Minimal two-state PAR2 activation (L + R <---> LR <---> LR*) and classic GPCR/G-protein activation mechanism: G-protein heterotrimers only interact with the receptor after receptor-activation.
  * `parm.classic_2f` - Full two-state PAR2 activation (R <---> R*, L + R <---> LR, L + R* <---> LR*, LR <---> LR*) and classic GPCR/G-protein activation mechanism: G-protein heterotrimers only interact with the receptor after receptor-activation.  
  * `parm.precoupled_1` - Single-state PAR2 activation and pre-coupled GPCR/G-protein activation mechanism: A pre-defined fraction of receptors are initialized with G-protein heterotrimers pre-coupled to the receptor.
  * `parm.precoupled_2m` - Minimal two-state PAR2 activation and pre-coupled GPCR/G-protein activation mechanism: A pre-defined fraction of receptors are initialized with G-protein heterotrimers pre-coupled to the receptor.
  * `parm.precoupled_2f` - Full two-state PAR2 activation and pre-coupled GPCR/G-protein activation mechanism: A pre-defined fraction of receptors are initialized with G-protein heterotrimers pre-coupled to the receptor.

We refer you to Fig 2 of Oliveira et al. [https://doi.org/10.3389/fnagi.2019.00089](https://doi.org/10.3389/fnagi.2019.00089) for a schematic depicting the main elements of the classic and pre-coupled GPCR/G-protein activation mechanisms.

Additionally, PARM contains 2 extensions of the `parm.classic_2m` model which incorporate an antagonist:

  * `parm.antagonist.competitive` - Classic GPCR/G-protein activation mechanism with a competitive antagonist.
  * `parm.antagonist.noncompetitive` - Classic GPCR/G-protein activation mechanism with a noncompetitive antagonist which operates via negative allosteric modulation of the agonist binding affinity. (Note: The factor which controls the allosteric modulation could also be set such that the antagonist induces positive allosteric modulation, increasing agonist binding affinity.)

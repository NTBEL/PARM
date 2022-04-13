# PARM
**P**AR2 **A**ctivation and calcium signaling **R**eaction **M**odel

PARM contains rules-based models of PAR2 (proteinase-activated receptor isoform 2) activation, G-protein activation, and calcium signaling via the phospholipase C and IP3 (inositol triphosphate) pathway. Models are encoded as Python program modules using the [PySB](http://pysb.org/) modeling framework. The models are designed to mathematically model the underlying signaling dynamics relevant to the inactivation of PAR2 in HEK-293 cells by Molecular Hyperthermia as described in:

  Kang et al.,  Transient Photoinactivation of Cell Membrane Protein Activity without Genetic Modification by Molecular Hyperthermia, ACS Nano 2019, 13, 11, 12487–12499 [https://doi.org/10.1021/acsnano.9b01993](https://doi.org/10.1021/acsnano.9b01993)


## Models in PARM
The core model of PARM is defined in `parm.parm`.

 Additionally, PARM contains 2 extensions of the `parm.parm` model which incorporate an antagonist:

  * `parm.antagonist.competitive` - Adds a competitive antagonist.
  * `parm.antagonist.noncompetitive` - Adds a noncompetitive antagonist which operates via negative allosteric modulation of the agonist binding affinity. (Note: The factor
    which controls the allosteric modulation could also be set such that the antagonist induces positive allosteric modulation, increasing agonist binding affinity.)

There are also 3 models with mechanistic variations defined in `parm.variants`:
  * `parm.variants.precoupled` - Adds pre-coupling between PAR2 and the G-protein heterotrimer such that some PAR2 can bind to the heterotrimer under resting conditions (without any agonist).
  * `parm.variants.heterogprot_cycle` - The receptor binding and G-protein interaction mechanism is adapted from the yeast G-protein cycle model of [Yi et al. 2003](https://doi.org/10.1073/pnas.1834247100).
  * `parm.variants.par2_synthesis_degradation` - This model only contains PAR2 with reactions for its resting synthesis and degradation. 

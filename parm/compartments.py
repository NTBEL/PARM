"""Functions to define the model compartments and their parameters.

This module defines two functions to setup the model compartments:
    * default_cell - defines a cell with the default cytosolic volume of 1 pL.
    * hek293_cell - defines a cell with size based on values more specific
        to HEK-293 cells.

Note that the two options each define the same 5 compartments:
    * EXTRACELLULAR - the extracellular space surrounding the cell.
        dimension = 3.
    * CELL_MEMB - the cell membrane.
        dimension = 2. parent = EXTRACELLULAR.
    * CYTOSOL - the cell's cytosol. dimension = 3. parent = CELL_MEMB.
    * ER_MEMB - the ER membrane.
        dimension = 2. parent = CYTOSOL.
    * ER_LUMEN - the ER lumen.
        dimension = 3. parent = ER_MEMB.

The volumes in pL of each compartment are set by parameters and the size of
each compartment is set as the ratio between that compartment's volume and the
volume of the cytosol (i.e., compartment sizes are fractional volumes taken
relative to the CYTOSOL compartment.) Under this convention all forward binding
rate constants can be defined relative to the cytosol compartment volume (and they
will be properly rescaled by the fractional compartment size). Concentrations
should also be defined in molecules per cell.

When comparment volume scaling is applied it should yield:

    kf_bind / (Vcompartment/Vcyto) = kf_bind * (Vcyto/Vcomparment),

such that kf_bind*Vcyto is the binding rate in pL/(s * molecule/cell) and then
division by Vcompartment returns the compartment-specific scaled binding
rate in 1/(s * molecule/cell).

Note:
    "For  elementary  bimolecular  reactions  with  the  rate  constant  given
    in  units  of volume/time,  the  rate  constant  is  divided  by the  volume
    of  the  reactant  compartment  with  the  highest  dimension" --
    https://www.informs-sim.org/wsc09papers/087.pdf Although the above says
    "volume/time" for units of the forward binding rate, it looks like it is
    more explicitly volume/(time*molecule/cell) so that dividing by the
    compartment volume would give you 1/(time * molecule/cell) (e.g.,
    1/s*molecule/cell), corresponding to concentrations given in molecule per
    cell. Also see BNGL example:
    https://github.com/RuleWorld/BNGTutorial/blob/master/CBNGL/LRR_comp.bngl

"""


# PySB components
from pysb import (
    Parameter,
    Expression,
    Annotation,
    Compartment,
    ANY,
)

from pysb.util import alias_model_components

# NumPy
import numpy as np

from . import units, defaults


def _define_compartments():
    """Defines the model compartments and the expressions defining their size."""
    alias_model_components()
    Expression("V_EXTRA", Vextra / Vcyto)
    Expression("V_CM", Vcm / Vcyto)
    Expression("V_C", Vcyto / Vcyto)
    Expression("V_ERM", Verm / Vcyto)
    Expression("V_ERL", Ver / Vcyto)
    alias_model_components()
    # Compartments
    # ============
    Compartment("EXTRACELLULAR", dimension=3, size=V_EXTRA)
    alias_model_components()
    # Cell Membrane
    Compartment("CELL_MEMB", dimension=2, parent=EXTRACELLULAR, size=V_CM)
    alias_model_components()
    # Cytosol
    Compartment("CYTOSOL", dimension=3, parent=CELL_MEMB, size=V_C)
    alias_model_components()
    #  ER membrane
    Compartment("ER_MEMB", dimension=2, parent=CYTOSOL, size=V_ERM)
    alias_model_components()
    # ER lumen volume
    Compartment("ER_LUMEN", dimension=3, parent=ER_MEMB, size=V_ERL)


def default_cell():
    """Defines the model compartments for a default cell.

    The default cell is defined with a default 1 pL cytosolic volume.

    Note that all volume parameters are given in pL, surface areas are given in
    micron^2, and thicknesses are given in micron.

    Adds 9 parameters:
        * Vcyto - the cytosolic volume. Taken as 1 pL.
        * SAcell - the surface area of the cell. Estimated from the
            cytosolic volume using a spherical cell assumption.
        * CMthickness - the cell membrane thickness. Taken as 0.01 micron.
        * Vcm - the effective volume of the cell membrane. Estimated from
            product of SAcell and CMthickness.
        * Vextra - the volume of the cytosolic space. Taken as 10*Vcyto (or
            10 pL).
        * Ver - the volume of the ER lumen. Taken as 0.185*Vcyto.
        * SAer - the surface area of the ER membrane. Estimated from Ver
            using a spherical compartment assumption.
        * ERMthickness - the thickness of ER membrane. Also taken as 0.01 micron.
        * Verm - the effective volume of the ER membrane. Estimated from
            product of SAer and ERMthickness.

    Adds 5 expressions:
        * V_EXTRA - fractional volume of the extracellular space relative to the
            cytosol volume: Vextra / Vcyto.
        * V_CM - fractional volume of the cell membrane relative to the
            cytosol volume: Vcm / Vcyto.
        * V_C - fractional volume of the cytosol relative itself:
            Vcyto / Vcyto = 1.
        * V_ERM - fractional volume of the ER membrane relative to the
            cytosol volume: Verm / Vcyto.
        * V_ERL - fractional volume of the ER lumen relative to the
            cytosol volume: Ver / Vcyto.

    Adds 5 compartments:

        * EXTRACELLULAR - the extracellular space surrounding the cell.
            dimension = 3. size = V_EXTRA.
        * CELL_MEMB - the cell membrane.
            dimension = 2. parent = EXTRACELLULAR. size = V_CM.
        * CYTOSOL - the cell's cytosol.
            dimension = 3. parent = CELL_MEMB. size = V_C.
        * ER_MEMB - the ER membrane.
            dimension = 2. parent = CYTOSOL. size = V_ERM.
        * ER_LUMEN - the ER lumen.
            dimension = 3. parent = ER_MEMB. size = V_ERL.

    """
    Parameter("Vcyto", defaults.VOL_CELL)  # 1 pL
    # Define the effective surface area to be a function of the spherical radius
    #  corresponding to the cytosolic volume.
    Rcell = (3 / 4 * np.pi * defaults.VOL_CELL / units.cubicmicron_to_pL) ** (1 / 3)
    # Cell-membrane surface area.
    Parameter("SAcell", 4 * np.pi * Rcell ** 2)
    # Effective cell-membrane thickness.
    # Assume 10 nm (0.01 micron) is reasonable.
    Parameter("CMthickness", defaults.MEMBRANE_THICKNESS)
    alias_model_components()
    # Effective volume of the cell-membrane
    Parameter("Vcm", SAcell.value * CMthickness.value * units.cubicmicron_to_pL)
    # Volume of the extracellular space
    # The following BNGL examples use 1000x the cell volume:
    #   https://github.com/RuleWorld/BNGTutorial/blob/master/CBNGL/LRR_comp.bngl
    #   https://github.com/RuleWorld/BNGTutorial/blob/master/CBNGL/LR_comp.bngl
    # We'll use 10 time the cytosolic volume
    Parameter("Vextra", Vcyto.value * 10)
    # Volume of the ER lumen/cisternal space.
    # It is often >10% of cell volume according to
    # Alberts et al. https://www.ncbi.nlm.nih.gov/books/NBK26841/ .
    # Lemon et al. 2003 use ratio of 0.185 ER-lumen/cytosol
    # Politi et al. 2006 also use ratio of 0.185 ER-lumen/cytosol
    # We can also use 0.185 ratio.
    Parameter("Ver", Vcyto.value * 0.185)
    alias_model_components()
    # Spherical radius corresponding to the ER volume.
    Rer = (3 / 4 * np.pi * Ver.value / units.cubicmicron_to_pL) ** (1 / 3)
    # Effective spherical surface area of the ER membrane.
    Parameter("SAer", 4 * np.pi * Rer ** 2)
    # Effective thickness of ER membrane.
    # Assume 10 nm is still reasonable.
    Parameter("ERMthickness", 0.01)
    alias_model_components()
    # Effective Volume of the ER membrane.
    Parameter("Verm", ERMthickness.value * SAer.value * units.cubicmicron_to_pL)
    alias_model_components()
    _define_compartments()
    return


def hek293_cell():
    """Defines the model compartments for a HEK-293 cell.

    The HEK-293 cell is defined with a size based on an experimentally measured
    cell diameter value for spherical HEK-293 cells in
    suspension. The cell diameter measurement of 13.9 +- 0.1 microns is from
    Mateus et al. Mol. Pharmaceutics 2013 10:6 2467-2478, https://doi.org/10.1021/mp4000822.

    Note that all volume parameters are given in pL, surface areas are given in
    micron^2, and thicknesses are given in micron.

    Adds 9 parameters:
        * Vcyto - the cytosolic volume. Estimated using a spherical cell
            diameter of 13.9 microns to define the total cell volume and an
            assumption that the cytosol is only 1/3 of the total cell volume
            (due to relatively large size of the nucleus):
                Vcyto = {(pi / 6) * (13.9 micron)^3} /  3 ~ 3.8 pL .
        * SAcell - the surface area of the cell. Estimated from the
            spherical cell diameter of 13.9 microns: ~ 607 micron^2.
        * CMthickness - the cell membrane thickness. Taken as 0.01 micron.
        * Vcm - the effective volume of the cell membrane. Estimated from
            product of SAcell and CMthickness.
        * Vextra - the volume of the cytosolic space. Taken as 10*Vcyto.
        * Ver - the volume of the ER lumen. Taken as 0.185*Vcyto.
        * SAer - the surface area of the ER membrane. Estimated from Ver
            using a spherical compartment assumption.
        * ERMthickness - the thickness of ER membrane. Also taken as 0.01 micron.
        * Verm - the effective volume of the ER membrane. Estimated from
            product of SAer and ERMthickness.

    Adds 5 expressions:
        * V_EXTRA - fractional volume of the extracellular space relative to the
            cytosol volume: Vextra / Vcyto.
        * V_CM - fractional volume of the cell membrane relative to the
            cytosol volume: Vcm / Vcyto.
        * V_C - fractional volume of the cytosol relative itself:
            Vcyto / Vcyto = 1.
        * V_ERM - fractional volume of the ER membrane relative to the
            cytosol volume: Verm / Vcyto.
        * V_ERL - fractional volume of the ER lumen relative to the
            cytosol volume: Ver / Vcyto.

    Adds 5 compartments:

        * EXTRACELLULAR - the extracellular space surrounding the cell.
            dimension = 3. size = V_EXTRA.
        * CELL_MEMB - the cell membrane.
            dimension = 2. parent = EXTRACELLULAR. size = V_CM.
        * CYTOSOL - the cell's cytosol.
            dimension = 3. parent = CELL_MEMB. size = V_C.
        * ER_MEMB - the ER membrane.
            dimension = 2. parent = CYTOSOL. size = V_ERM.
        * ER_LUMEN - the ER lumen.
            dimension = 3. parent = ER_MEMB. size = V_ERL.

    """
    # HEK-293 cells are roughly spherical with average diameter of
    # 13.9 +- 0.1 micron in suspension as per Mateus et al. 2013 https://doi.org/10.1021/mp4000822
    diameter_cell = 13.9  # micron
    radius_cell = 0.5 * diameter_cell  # micron
    V_cell = (4 / 3) * np.pi * diameter_cell ** 3 * units.cubicmicron_to_pL  # pL
    # The nucleus of HEK293 cells is quite large. I don't have an exact
    # volume fraction for the cytosol, but using a rough estimation from the
    # IHC images the cytosol is only about 30-50% of the cell volume.
    Parameter("Vcyto", V_cell / 3)
    # Estimate cell-membrane surface area assuming a spherical cell with
    # radius corresponding to the cell diameter of 13.9 microns.
    # This yields a surface area of ~ 607 micron^2.
    Parameter("SAcell", 4 * np.pi * radius_cell ** 2)
    # Note that Blumlein et al. Scientific Reports 2017 7:7346,
    #   https://doi.org/10.1038/s41598-017-07813-5
    # report a surface area ~ 706.5 micron^2 for Hek cells.

    # Effective cell-membrane thickness
    # Assume 10 nm (0.01 micron) is reasonable.
    Parameter("CMthickness", defaults.MEMBRANE_THICKNESS)
    alias_model_components()
    # Effective volume of the cell-membrane
    Parameter("Vcm", SAcell.value * CMthickness.value * units.cubicmicron_to_pL)
    # Volume of the extracellular space
    # The following BNGL examples use 1000x the cell volume:
    #   https://github.com/RuleWorld/BNGTutorial/blob/master/CBNGL/LRR_comp.bngl
    #   https://github.com/RuleWorld/BNGTutorial/blob/master/CBNGL/LR_comp.bngl
    # We'll use 10 time the cytosolic volume
    Parameter("Vextra", Vcyto.value * 10)
    # Volume of the ER lumen/cisternal space.
    # It is often >10% of cell volume according to
    # Alberts et al. https://www.ncbi.nlm.nih.gov/books/NBK26841/ .
    # Lemon et al. 2003 use ratio of 0.185 ER-lumen/cytosol
    # Politi et al. 2006 also use ratio of 0.185 ER-lumen/cytosol
    # We can also use 0.185 ratio.
    Parameter("Ver", Vcyto.value * 0.185)
    alias_model_components()
    # Spherical radius corresponding to the ER volume.
    Rer = (3 / 4 * np.pi * Ver.value / units.cubicmicron_to_pL) ** (1 / 3)
    # Effective spherical surface area of the ER membrane.
    Parameter("SAer", 4 * np.pi * Rer ** 2)
    # Effective thickness of ER membrane.
    # Assume 10 nm (0.01 micron) is still reasonable.
    Parameter("ERMthickness", defaults.MEMBRANE_THICKNESS)
    alias_model_components()
    # Effective Volume of the ER membrane.
    Parameter("Verm", ERMthickness.value * SAer.value * units.cubicmicron_to_pL)
    alias_model_components()
    _define_compartments()
    return

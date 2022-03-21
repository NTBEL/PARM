"""Defines functions for add PAR2 antagonism.

It includes competitive and noncompetive antagonsim:
    * competitive_par2_antagonist
    * noncompetitive_par2_antagonist

"""


# PySB components
from pysb import (
    Monomer,
    Parameter,
    Initial,
    Rule,
    Observable,
    Expression,
    Annotation,
    ANY,
    WILD,
)

# PySB macros
from pysb.macros import (
    bind,
    bind_complex,
    catalyze,
    catalyze_complex,
    catalyze_state,
    degrade,
    synthesize,
)


from pysb.util import alias_model_components
from .. import defaults, units


def antagonist_monomer_initial():
    """Defines the antagonist and it's initial condition.

    Adds 1 monomer and 1 parameter. Defines the Initial for the antagonist.

    Monomers:
        * Ant - antagonist molecule.

    Parameters:
        Ant_0 - the initial amount of antagonist in the extracellurar space.
    """
    # Antagonist
    Monomer("Ant", ["b"])
    alias_model_components()
    # Just set a nominal value of 100 nM.
    Parameter("Ant_0", 100 * units.nM_to_molec_per_pL * Vextra.value)
    alias_model_components()
    Initial(Ant(b=None) ** EXTRACELLULAR, Ant_0)
    return


def competitive_par2_antagonist():
    """Defines a competitive PAR2 antagonist that binds the same location as 2AT.

    This function adds an antagonist moleculer and reactions for competitive
    antagonism of the orthosteric 2AT binding location on PAR2. Adds the
    reaction:
        PAR2_I + Ant <---> PAR2_I:Ant

    Calls:
        * antagonist_monomer_initial

    Adds 2 additional parameters and 1 additioanl expression.

    Parameters:
        * kf_PAR2_bind_Ant - forward rate constant for antagonist binding to
            PAR2 (orthosteric site).
        * Kd_PAR2_bind_Ant - dissociation constant for antagonist binding to
            PAR2 (orthosteric site).

    Expressions:
        * kr_PAR2_bind_Ant - reverse rate consant for antagonist binding to
            PAR2 (orthosteric site).
    """
    antagonist_monomer_initial()
    alias_model_components()
    # Antagonist binding to orthosteric site of PAR2
    # 'Diffusion limited' forward rate of 1000 1/s*microM used in CORM https://github.com/LoLab-VU/CORM
    Parameter("kf_PAR2_bind_Ant", (1000 / units.microM_to_molec_per_pL) / Vextra.value)
    # For nominal value just assume Kd is 100 nM
    Parameter("Kd_PAR2_bind_Ant", 100 * units.nM_to_molec_per_pL * Vextra.value)
    alias_model_components()
    # Division by V_EXTRA is take into account the compartment scaling for kf.
    Expression("kr_PAR2_bind_Ant", Kd_PAR2_bind_Ant * kf_PAR2_bind_Ant / V_EXTRA)
    alias_model_components()
    # Competitive binding of the Antoginst
    #   Ant + PAR2_I <---> Ant:PAR2_I
    bind(
        Ant() ** EXTRACELLULAR,
        "b",
        PAR2(state="I") ** CELL_MEMB,
        "bortho",
        [kf_PAR2_bind_Ant, kr_PAR2_bind_Ant],
    )

    return


def noncompetitive_par2_antagonist():
    """Defines a noncompetitive PAR2 antagonist that binds an allosteric location relative to 2AT.

    This function adds an antagonist molecule and reactions for noncompetitive
    antagonism through allosteric modulation of 2AT binding/dissociation on
    PAR2. Adds the reactions:
        i. PAR2_I + Ant <---> PAR2_I:Ant
       ii. Ant + PAR2_A <---> Ant:PAR2_A
      iii. 2AT + PAR2_I:Ant <---> TAT:PAR2_I:Ant
       iv. 2AT + PAR2_A:Ant <---> TAT:PAR2_A:Ant

    Calls:
        * antagonist_monomer_initial

    Adds 3 additional parameters and 3 additional expressions.

    Parameters:
        * kf_PAR2_bind_Ant - forward rate constant for antagonist binding to
            PAR2 (orthosteric site).
        * Kd_PAR2_bind_Ant - dissociation constant for antagonist binding to
            PAR2 (orthosteric site).
        * allosteric_modulation_factor - allosteric modulation factor that
            alters the dissociation of 2AT from PAR2 when the antagonist is
            bound in allosteric site.

    Expressions:
        * kr_PAR2_bind_Ant - reverse rate consant for antagonist binding to
            PAR2 (orthosteric site).
        * Kd_PAR2_bind_TAT_allo - allosterically modulated dissociation constant
            for 2AT binding to PAR2. Given by:
                Kd_PAR2_bind_TAT_allo = Kd_PAR2_bind_TAT * allosteric_modulation_factor
        * kr_PAR2_bind_TAT_allo - reverse rate consant for 2AT binding PAR2
            when antagonist is bound in the allosteric site of PAR2.
    """
    antagonist_monomer_initial()
    alias_model_components()
    # Antagonist binding to orthosteric site of PAR2
    # 'Diffusion limited' forward rate of 1000 1/s*microM used in CORM https://github.com/LoLab-VU/CORM
    Parameter("kf_PAR2_bind_Ant", (1000 / units.microM_to_molec_per_pL) / Vextra.value)
    # For nominal value just assume Kd is 100 nM
    Parameter("Kd_PAR2_bind_Ant", 100 * units.nM_to_molec_per_pL * Vextra.value)
    alias_model_components()
    # Division by V_EXTRA is take into account the compartment scaling for kf.
    Expression("kr_PAR2_bind_Ant", Kd_PAR2_bind_Ant * kf_PAR2_bind_Ant / V_EXTRA)
    # Assume the allosteric modulation affects the Kd for agonist binding, but
    # that the forward binding rate is unaffected.
    # For nominal value assume that the Kd of 2AT is increased by a factor of 2.
    Parameter("allosteric_modulation_factor", 2)
    alias_model_components()
    Expression("Kd_PAR2_bind_TAT_allo", Kd_PAR2_bind_TAT * allosteric_modulation_factor)
    alias_model_components()
    Expression(
        "kr_PAR2_bind_TAT_allo", Kd_PAR2_bind_TAT_allo * kf_PAR2_bind_TAT / V_EXTRA
    )
    alias_model_components()
    # nonompetitive binding of the Antoginst to an allosteric site
    #   Ant + PAR2_I <---> Ant:PAR2_I
    bind(
        Ant() ** EXTRACELLULAR,
        "b",
        PAR2(state="I", bortho=WILD, bgaq=WILD) ** CELL_MEMB,
        "ballo",
        [kf_PAR2_bind_Ant, kr_PAR2_bind_Ant],
    )
    #   Ant + PAR2_A <---> Ant:PAR2_A
    bind(
        Ant() ** EXTRACELLULAR,
        "b",
        PAR2(state="A", bortho=WILD, bgaq=WILD) ** CELL_MEMB,
        "ballo",
        [kf_PAR2_bind_Ant, kr_PAR2_bind_Ant],
    )
    #    2AT + PAR2_I:Ant <---> TAT:PAR2_I:Ant
    bind(
        TAT() ** EXTRACELLULAR,
        "b",
        PAR2(state="I", ballo=ANY, bgaq=WILD) ** CELL_MEMB,
        "bortho",
        [kf_PAR2_bind_TAT, kr_PAR2_bind_TAT_allo],
    )
    #    2AT + PAR2_A:Ant <---> TAT:PAR2_A:Ant
    bind(
        TAT() ** EXTRACELLULAR,
        "b",
        PAR2(state="A", ballo=ANY, bgaq=WILD) ** CELL_MEMB,
        "bortho",
        [kf_PAR2_bind_TAT, kr_PAR2_bind_TAT_allo],
    )

    return

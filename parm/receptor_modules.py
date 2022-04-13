"""Defines functions for receptor activation, degradation, synthesis, and receptor observables.

It includes three alternate mechanisms for receptor activation:
    * minimal_two_state_par2_activation
    * single_state_par2_activation
    * full_two_state_par2_activation
As well as functions to add receptor degradation and synthesis:
    * occupied_par2_degradation
    * unoccupied_active_par2_degradation
    * denatured_par2_degradation
    * par2_synthesis
And addon functions for precoupled receptor activation:
    * addon_single_state_precoupled_par2_activation
    * addon_minimal_twostate_precoupled_par2_activation
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
    equilibrate,
)

from pysb.util import alias_model_components

# NumPy
import numpy as np

# Avogadro's Number from scipy
from scipy.constants import N_A
from sympy.functions.elementary.miscellaneous import Max

from . import defaults, units


def receptor_agonist_monomers():
    """Declares the monomers for 2AT (TAT) and PAR2.
    Adds two monomers and 1 annotation.

    Monomers:
        * TAT - PAR2-agonist 2AT.
        * PAR2 - the receptor.
    """
    # PAR2 agonist 2AT as in Kang et al. https://doi.org/10.1021/acsnano.9b01993
    Monomer("TAT", ["b"])
    # PAR2:
    #  binding site:
    #    bortho = orthosteric 2AT binding site
    #    ballo = allosteric binding site
    #    bgaq = G protein binding site
    #  states: I = inactive, A = active, D = denatured
    Monomer("PAR2", ["bortho", "ballo", "bgaq", "state"], {"state": ["I", "A", "D"]})
    alias_model_components()
    # Annotations
    # ===========
    Annotation(PAR2, "https://identifiers.org/uniprot:P55085")
    return


def receptor_monomer_only():
    """Declares the monomer PAR2.
    Adds one monomers and 1 annotation.

    Monomers:
        * PAR2 - the receptor.
    """
    # PAR2:
    #  binding site:
    #    bortho = orthosteric 2AT binding site
    #    ballo = allosteric binding site
    #    bgaq = G protein binding site
    #  states: I = inactive, A = active, D = denatured
    Monomer("PAR2", ["bortho", "ballo", "bgaq", "state"], {"state": ["I", "A", "D"]})
    alias_model_components()
    # Annotations
    # ===========
    Annotation(PAR2, "https://identifiers.org/uniprot:P55085")
    return


def agonist_initial():
    """Declares the initial condition for 2AT (TAT).

    Adds 1 parameter and 1 annotation.

    Parameters:
        * TAT_0 - initial concentration of 2AT in the extracellular space.
    """
    alias_model_components()
    # PAR2 agonist 2AT,
    # 330 nM for Fig 2D data from Kang et al. https://doi.org/10.1021/acsnano.9b01993
    C_2AT = 330  # nM
    Parameter("TAT_0", C_2AT * units.nM_to_molec_per_pL * Vextra.value)
    alias_model_components()
    Initial(TAT(b=None) ** EXTRACELLULAR, TAT_0)
    Annotation(
        TAT_0, "https://doi.org/10.1021/acsnano.9b01993", predicate="isDescribedBy"
    )
    return


def receptor_initials():
    """Declares initial conditions for different PAR2 states.
    Adds 2 parameters and 1 expression.

    Parameters:
        * PAR2_0 - the initial concentration of all PAR2 in the cell membrane.
        * f_denature - fraction of PAR2 denatured by photoinactivation.

    Expressions:
        * PAR2_0_D - amount of PAR2 that has been denatured through photoinactivation.
        * PAR2_0_I -  amount of initial inactive PAR2.
    """
    # total PAR2
    # From Falkenburger et al. 2010 https://dx.doi.org/10.1085%2Fjgp.200910344
    # # tsA201 cells
    #     Endogenous receptor density: 1/micrometer^2
    #     Overexpressed receptor density: 3,000/micrometer^2
    # From Brinkerhoff, Choi, and Linderman 2008 https://doi.org/10.1016/j.jtbi.2008.01.002
    #    Receptor concentration 2e3 to 2e4 /cell
    # Our current estimate of PAR2/HEK293-cell from the AuNP/cell binding
    # saturation curve is ~160.
    Parameter("PAR2_0", 160)
    Parameter("f_denature", 0.0)  # By default no PAR2 has been denatured.
    alias_model_components()
    Expression("PAR2_0_D", f_denature * PAR2_0)
    alias_model_components()
    Expression("PAR2_0_I", Max(PAR2_0 - PAR2_0_D, 0))
    alias_model_components()
    Initial(PAR2(state="I", bortho=None, ballo=None, bgaq=None) ** CELL_MEMB, PAR2_0_I)
    Initial(PAR2(state="D", bortho=None, ballo=None, bgaq=None) ** CELL_MEMB, PAR2_0_D)
    return


def minimal_two_state_par2_activation():
    """Defines a minimal two-state activation mechansim of PAR2 by 2AT.

    The miminal two-state activation mechanism is given by:
        2AT + PAR2_I <---> TAT:PAR2_I <---> TAT:PAR2_A

    Calls functions:
        * receptor_agonist_monomers
        * agonist_initial
        * receptor_initials

    Declares initial conditions for PAR2 in state I and state D.

    Adds 4 additional parameters and 2 additional expressions.

    Parameters:
        * kf_PAR2_bind_TAT - forward rate constant for 2AT binding to PAR2.
        * Kd_PAR2_bind_TAT - equilibrium dissociation constant for 2AT binding to PAR2.
        * k_activate_PAR2 - forward rate constant for conversion of occupied PAR2 to active PAR2.
        * k_inactivate_PAR2 - reverse rate constant for conversion of occupied PAR2 to active PAR2.
    Expressions:
        * PAR2_0_I - initial concentration of free inactive PAR2.
        * kr_PAR2_bind_TAT - reverse rate constant for 2AT binding to PAR2.

    """

    receptor_agonist_monomers()
    agonist_initial()
    receptor_initials()
    # Needed to recognize the monomer and parameter names in the present scope
    alias_model_components()
    # Binding forward rate constant for ligand-receptor binding in the yeast
    # G-protein cycle model of  Yi et al. 2003 PNAS https://doi.org/10.1073/pnas.1834247100
    # is 2x10^6 1/(M*s). Converting for our model that would be:
    #    2x10^6 1/(M*s) / (10^-12 L) / N_A = 3x10^-6 1/(molecule/pL*s) = 3*KF_BIND
    # We'll use that value as our initial estimate.
    # Note we divide KF_BIND by Vcyto to convert into 1/(molecule*s) which will
    # ultimately be rescaled by Vextra/Vcyto.
    Parameter("kf_PAR2_bind_TAT", defaults.KF_BIND * 3 / Vcyto.value)
    # PAR2 activation by 2AT
    # Note: Ca2+ signal Max. FRET Dose-Response for 2AT activation of PAR2
    # has EC50 = 101.7 +- 28.7 nM, Kang et al. https://doi.org/10.1021/acsnano.9b01993
    # PAR2 agonists in HEK 293T cells - LeSarge et al. https://doi.org/10.1021/acsmedchemlett.9b00094
    #   2f-LIGRLO(Sulfo-Cy5)-NH2 has Kd = 430 nM with EC50 = 296 nM
    #   Isox-Cha-Chg-ARK(Sulfo-Cy5)-NH2 has Kd = 38 nM with EC50 = 16 nM
    # Since 2AT has EC50 = 101.8 nM in Hek 293 cells its probably safe to
    # assume that the Kd for 2AT is somewhere between those two compounds.
    # 142 nM = (430-38)/(296-16) *101.5
    Parameter("Kd_PAR2_bind_TAT", 142 * units.nM_to_molec_per_pL * Vextra.value)
    Parameter("k_activate_PAR2", defaults.K_CONVERT * 10)
    Parameter("k_inactivate_PAR2", defaults.K_CONVERT / 10)
    alias_model_components()
    # Division by V_EXTRA is take into account the compartment scaling for kf.
    Expression("kr_PAR2_bind_TAT", Kd_PAR2_bind_TAT * kf_PAR2_bind_TAT / V_EXTRA)
    alias_model_components()
    # Alias the TAT:PAR2 complexes
    tat_PAR2_i = (
        TAT(b=1) ** EXTRACELLULAR
        % PAR2(state="I", bortho=1, ballo=None, bgaq=None) ** CELL_MEMB
    )
    tat_PAR2_a = (
        TAT(b=1) ** EXTRACELLULAR
        % PAR2(state="A", bortho=1, ballo=None, bgaq=None) ** CELL_MEMB
    )
    # 2-step activation of PAR2 by 2AT agonist:
    # i.   2AT + PAR2_I <---> TAT:PAR2_I
    Rule(
        "tat_bind_PAR2",
        TAT(b=None) ** EXTRACELLULAR
        + PAR2(state="I", bortho=None, ballo=None, bgaq=None) ** CELL_MEMB
        | tat_PAR2_i,
        kf_PAR2_bind_TAT,
        kr_PAR2_bind_TAT,
    )
    # ii.  TAT:PAR2_I <---> TAT:PAR2_A
    Rule(
        "tat_activate_PAR2", tat_PAR2_i | tat_PAR2_a, k_activate_PAR2, k_inactivate_PAR2
    )
    return


def single_state_par2_activation():
    """Defines a single-state activation mechansim of PAR2 by 2AT.
    The single-state activation mechanism is given by:
        2AT + PAR2_I <---> TAT:PAR2_A

    Calls functions:
        * receptor_agonist_monomers
        * agonist_initial
        * receptor_initials

    Declares initial conditions for PAR2 in state I and state D.

    Adds 2 additional parameters and 2 additional expressions.

    Parameters:
        * kf_PAR2_bind_TAT - forward rate constant for 2AT binding to PAR2.
        * Kd_PAR2_bind_TAT - equilibrium dissociation constant for 2AT binding to PAR2.

    Expressions:
        * PAR2_0_I - initial concentration of free inactive PAR2.
        * kr_PAR2_bind_TAT - reverse rate constant for 2AT binding to PAR2.
    """
    receptor_agonist_monomers()
    agonist_initial()
    receptor_initials()
    alias_model_components()
    # Binding forward rate constant for ligand-receptor binding in the yeast
    # G-protein cycle model of  Yi et al. 2003 PNAS https://doi.org/10.1073/pnas.1834247100
    # is 2x10^6 1/(M*s). Converting for our model that would be:
    #    2x10^6 1/(M*s) / (10^-12 L) / N_A = 3x10^-6 1/(molecule*s) = 3*KF_BIND
    # We'll use that value as our initial estimate.
    Parameter("kf_PAR2_bind_TAT", defaults.KF_BIND * 3 / Vcyto.value)
    # PAR2 activation by 2AT
    # Note: Ca2+ signal Max. FRET Dose-Response for 2AT activation of PAR2
    # has EC50 = 101.7 +- 28.7 nM, Kang et al. https://doi.org/10.1021/acsnano.9b01993
    # PAR2 agonists in HEK 293T cells - LeSarge et al. https://doi.org/10.1021/acsmedchemlett.9b00094
    #   2f-LIGRLO(Sulfo-Cy5)-NH2 has Kd = 430 nM with EC50 = 296 nM
    #   Isox-Cha-Chg-ARK(Sulfo-Cy5)-NH2 has Kd = 38 nM with EC50 = 16 nM
    # Since 2AT has EC50 = 101.8 nM in Hek 293 cells its probably safe to
    # assume that the Kd for 2AT is somewhere between those two compounds.
    # 142 nM = (430-38)/(296-16) *101.5
    Parameter("Kd_PAR2_bind_TAT", 142 * units.nM_to_molec_per_pL * Vcyto.value)
    alias_model_components()
    # Division by V_EXTRA is take into account the compartment scaling for kf.
    Expression("kr_PAR2_bind_TAT", Kd_PAR2_bind_TAT * kf_PAR2_bind_TAT / V_EXTRA)
    alias_model_components()
    # Alias the TAT:PAR2 complexes
    tat_PAR2_i = (
        TAT(b=1) ** EXTRACELLULAR
        % PAR2(state="I", bortho=1, ballo=None, bgaq=None) ** CELL_MEMB
    )
    tat_PAR2_a = (
        TAT(b=1) ** EXTRACELLULAR
        % PAR2(state="A", bortho=1, ballo=None, bgaq=None) ** CELL_MEMB
    )
    # Single state receptor activation (induced fit) by 2AT.
    #    2AT + PAR2_I <---> TAT:PAR2_A
    Rule(
        "tat_bind_PAR2",
        TAT(b=None) ** EXTRACELLULAR
        + PAR2(state="I", bortho=None, ballo=None, bgaq=None) ** CELL_MEMB
        | tat_PAR2_a,
        kf_PAR2_bind_TAT,
        kr_PAR2_bind_TAT,
    )
    return


def constitutive_par2_activity():
    """Defines constitutive activity for PAR2.
    Adds isomerization of PAR2 between its inactive and active forms so that
    PAR2 can become active without agonist:
        PAR2_I <---> PAR2_A

    Also includes enhanced binding of 2AT to the active form of PAR2:
        2AT + PAR2_A <---> TAT:PAR2_A

    Declares a binding_amplification_factor parametr that controls how much
    stronger the 2AT binding is to the active form of PAR2.

    Adds 3 parameters and 2 expressions.

    Parameters:
        * kf_PAR2_constitutively_active - forward rate constant for
            isomerization of PAR2 between inactive and active states:
            R ---> R* , kf
        * kr_PAR2_constitutively_active - reverse rate constant for
            isomerization of PAR2 between inactive and active states:
            R* ---> R , kr
        * binding_amplification_factor - the factor by which 2AT binding to
                the active form of PAR2 is amplified. Should be > 0. The
                dissociation constant for 2AT binding to the inactive form
                of PAR2 is divided by this factor to define the
                corresponding dissociation constant for 2AT binding to the
                active form of PAR2:
                    > 1 = stronger binding than to inactive PAR2.
                    1 = same binding strength as to inactive PAR2.
                    < 1 = weaker bindng than to inactive PAR2.

    Expressions:
        * Kd_PAR2_A_bind_TAT - equilibrium dissociation constant for 2AT
            and the active form of PAR2.
        * kr_PAR2_A_bind_TAT - reverse rate constant for 2AT binding to
            active PAR2.
    """
    Parameter("kf_PAR2_constitutively_active", 1e-3)  # 1/s
    Parameter("kr_PAR2_constitutively_active", 1e-1)  # 1/s
    alias_model_components()
    # Equilibrium between the inactive and active form of the receptor.
    # R <---> R*
    # PAR2_I <---> PAR2_A
    equilibrate(
        PAR2(state="I", bortho=None, ballo=None, bgaq=None) ** CELL_MEMB,
        PAR2(state="A", bortho=None, ballo=None, bgaq=None) ** CELL_MEMB,
        [kf_PAR2_constitutively_active, kr_PAR2_constitutively_active],
    )
    #    2AT + PAR2_A <---> TAT:PAR2_A
    # We'll assume that agonist binds to active receptor with the same forward
    # rate constant (assumed to be diffusion limited) as for the inactive
    # receptor but that the binding affinity can be amplified so the active receptor
    # has a smaller dissociation constant and thus a slower unbinding rate
    # rate constant.
    Parameter("binding_amplification_factor", 1.5)
    alias_model_components()
    Expression("Kd_PAR2_A_bind_TAT", Kd_PAR2_bind_TAT / binding_amplification_factor)
    alias_model_components()
    # Division by V_EXTRA is take into account the compartment scaling for kf.
    Expression("kr_PAR2_A_bind_TAT", Kd_PAR2_A_bind_TAT * kf_PAR2_bind_TAT / V_EXTRA)
    alias_model_components()
    tat_PAR2_a = (
        TAT(b=1) ** EXTRACELLULAR
        % PAR2(state="A", bortho=1, ballo=None, bgaq=None) ** CELL_MEMB
    )
    Rule(
        "tat_bind_PAR2_A",
        TAT(b=None) ** EXTRACELLULAR
        + PAR2(state="A", bortho=None, ballo=None, bgaq=None) ** CELL_MEMB
        | tat_PAR2_a,
        kf_PAR2_bind_TAT,
        kr_PAR2_A_bind_TAT,
    )
    return


def full_two_state_par2_activation():
    """Defines full two-state PAR2 activation by 2AT.
    Combines minimal two-state activation with constitutive PAR2 activity to
    encode a full two-state activation mechanism:
        L + R  <---> LR
          |           |
        L + R* <---> LR*
    Or reactions:
        i. 2AT + PAR2_I <---> TAT:PAR2_I <---> TAT:PAR2_A
       ii. PAR2_I <---> PAR2_A
      iii. 2AT + PAR2_A <---> TAT:PAR2_A

    Calls:
        * minimal_two_state_par2_activation
        * constitutive_par2_activity
    """
    minimal_two_state_par2_activation()
    constitutive_par2_activity()
    return


def occupied_par2_degradation():
    """Defines degradation of 2AT bound PAR2.
    1st order degradation of occupied PAR2 receptors:
        PAR2:2AT ---> None

    Adds 1 parameter.

    Parameters:
        * k_PAR2_bound_degradation - rate constant for 1st order degradation of
            2AT-occupied PAR2.
    """
    # PAR2 degradation

    # Occupied PAR2:
    # rate constant for ligand-bound receptor degradation from
    # Yi et al. 2003 PNAS https://doi.org/10.1073/pnas.1834247100
    # used in yeast G-protein cycle model is: 4x10^-3 1/s
    Parameter("k_PAR2_bound_degradation", 4e-3)
    alias_model_components()
    # Alias the TAT:PAR2 complex

    tat_PAR2 = (
        TAT(b=1) ** EXTRACELLULAR % PAR2(bortho=1, ballo=WILD, bgaq=None) ** CELL_MEMB
    )
    # PAR2 degradation
    # Occupied any state:
    degrade(tat_PAR2, k_PAR2_bound_degradation)
    return


def denatured_par2_degradation():
    """Defines degradation of denature PAR2.
    1st order denatured PAR2 degradation:
        PAR2_D:2AT ---> None

    Adds one parameter.

    Parameters:
        * k_PAR2_denatured_degradation - rate constant for 1st order degradation
            of denatured (state='D') PAR2.
    """
    # Denatured PAR2:
    # As a default we can set the denatured PAR2 degradation to have the same
    # rate constant as bound PAR2 degradation, but we'll assume it could be
    # different.
    Parameter("k_PAR2_denatured_degradation", 4e-3)
    alias_model_components()
    # Denatured:
    degrade(
        PAR2(state="D", bortho=None, ballo=None, bgaq=None) ** CELL_MEMB,
        k_PAR2_denatured_degradation,
    )
    return


def unoccupied_active_par2_degradation():
    """Defines degradation of free but active PAR2.
    1st order active PAR2 degradation:
        PAR2_A ---> None

    Adds 1 parameter.

    Parameters:
        * k_PAR2_free_degradation - rate constant for 1st order degradation of
            free but active PAR2.
    """
    # Free unbound but active PAR2:
    # rate constant for ligand-bound receptor degradation from
    # Yi et al. 2003 PNAS https://doi.org/10.1073/pnas.1834247100
    # used in yeast G-protein cycle model is: 4x10^-3 1/s
    Parameter("k_PAR2_free_degradation", 4e-3)
    alias_model_components()
    degrade(
        PAR2(state="A", bortho=None, ballo=None, bgaq=None) ** CELL_MEMB,
        k_PAR2_free_degradation,
    )
    return


def inactive_par2_degradation():
    """Defines degradation of free inactive PAR2.
    1st order PAR2 degradation:
        PAR2_I ---> None

    Adds 1 parameter.

    Parameters:
        * k_PAR2_I_degradation - rate constant for 1st order degradation of
            free inactive PAR2.
    """
    # Free unbound and inactive PAR2:
    # rate constant for receptor degradation from
    # Yi et al. 2003 PNAS https://doi.org/10.1073/pnas.1834247100
    # used in yeast G-protein cycle model is: 4x10^-4 1/s
    Parameter("k_PAR2_I_degradation", 4e-4)
    alias_model_components()
    degrade(
        PAR2(state="I", bortho=None, ballo=None, bgaq=None) ** CELL_MEMB,
        k_PAR2_I_degradation,
    )
    return


def par2_synthesis():
    """Defines synthesis of PAR2.
    Zero order synthesis of PAR2:
        None ---> PAR2_I

    Adds 1 parameter.

    Parameters:
        * k_PAR2_synthesis - rate constant for zero-order synthesis of PAR2.
    """
    alias_model_components()
    # Zero-order synthesis of PAR2.
    # The rate constant for receptor synthesis from
    # Yi et al. 2003 PNAS https://doi.org/10.1073/pnas.1834247100
    # used in yeast G-protein cycle model is: 4 molecules/s
    # Rate gets multiplied by the compartment size so we divide by the
    # fractional volume of the cell membrane to keep right value in molecules/cell*s.
    Parameter("k_PAR2_synthesis", 4 / V_CM.get_value())
    alias_model_components()
    synthesize(
        PAR2(state="I", bortho=None, ballo=None, bgaq=None) ** CELL_MEMB,
        k_PAR2_synthesis,
    )
    # Rule("synthesize_PAR2", None >> PAR2(state="I", bortho=None, ballo=None, bgaq=None) ** CELL_MEMB, k_PAR2_synthesis)
    return


def addon_minimal_twostate_precoupled_par2_activation():
    """Adds activation of PAR2 precoupled to G-protein to the mimimal two-state receptor activation mechanism.
    This addon function uses the G-protein heterotrimer so must be called after
    the G-protein mechanism modules. It adds two additional rules for 2AT binding
    to the preceoupled receptor and reversible activation of receptor once
    2AT is bound:
        TAT + PAR2_I:Gaq:GDP:Gbg <---> TAT:PAR2_I:Gaq:GDP:Gbg
        TAT:PAR2_I:Gaq:GDB:Gbg <---> TAT:PAR2_A:Gaq:GDP:Gbg

    Adds no additional parameters or expressions.
    """
    alias_model_components()
    # Alias the pre-coupled PAR2-Gprotein complex
    PAR2_i_Gaq_gdp_Gbg = (
        PAR2(state="I", bortho=None, ballo=None, bgaq=2) ** CELL_MEMB
        % Gaq(bpar=2, bgdp=3, bgbg=4) ** CELL_MEMB
        % GDP(b=3) ** CELL_MEMB
        % Gbg(b=4) ** CELL_MEMB
    )
    # Alias the complex 2AT:PAR2_I:Gaq:GDP:Gbg
    tat_PAR2_i_Gaq_gdp_Gbg = (
        TAT(b=1) ** EXTRACELLULAR
        % PAR2(state="I", bortho=1, ballo=None, bgaq=2) ** CELL_MEMB
        % Gaq(bpar=2, bgdp=3, bgbg=4) ** CELL_MEMB
        % GDP(b=3) ** CELL_MEMB
        % Gbg(b=4) ** CELL_MEMB
    )
    # Alias the complex 2AT:PAR2_A:Gaq:GDP:Gbg
    tat_PAR2_a_Gaq_gdp_Gbg = (
        TAT(b=1) ** EXTRACELLULAR
        % PAR2(state="A", bortho=1, ballo=None, bgaq=2) ** CELL_MEMB
        % Gaq(bpar=2, bgdp=3, bgbg=4) ** CELL_MEMB
        % GDP(b=3) ** CELL_MEMB
        % Gbg(b=4) ** CELL_MEMB
    )
    #    2AT + PAR2_I:Gaq:GDP:Gbg <---> TAT:PAR2_I:Gaq:GDP:Gbg
    Rule(
        "tat_bind_PAR2_pre",
        TAT(b=None) ** EXTRACELLULAR + PAR2_i_Gaq_gdp_Gbg | tat_PAR2_i_Gaq_gdp_Gbg,
        kf_PAR2_bind_TAT,
        kr_PAR2_bind_TAT,
    )
    #    TAT:PAR2_I:Gaq:GDB:Gbg <---> TAT:PAR2_A:Gaq:GDP:Gbg
    Rule(
        "tat_activate_PAR2_pre",
        tat_PAR2_i_Gaq_gdp_Gbg | tat_PAR2_a_Gaq_gdp_Gbg,
        k_activate_PAR2,
        k_inactivate_PAR2,
    )
    return


def addon_single_state_precoupled_par2_activation():
    """Adds activation of PAR2 precoupled to G-protein to the single state receptor activation mechanism.
    This addon function uses the G-protein heterotrimer so must be called after
    the G-protein mechanism modules. It adds one additional rule for concerted
    2AT binding to the preceoupled receptor and receptor activation:
        TAT + PAR2_I:Gaq:GDP:Gbg <---> TAT:PAR2_A:Gaq:GDP:Gbg

    Adds no additional parameters or expressions.
    """
    alias_model_components()
    # Alias the pre-coupled PAR2-Gprotein complex
    PAR2_i_Gaq_gdp_Gbg = (
        PAR2(state="I", bortho=None, ballo=None, bgaq=2) ** CELL_MEMB
        % Gaq(bpar=2, bgdp=3, bgbg=4) ** CELL_MEMB
        % GDP(b=3) ** CELL_MEMB
        % Gbg(b=4) ** CELL_MEMB
    )
    # Alias the complex 2AT:PAR2_A:Gaq:GDP:Gbg
    tat_PAR2_a_Gaq_gdp_Gbg = (
        TAT(b=1) ** EXTRACELLULAR
        % PAR2(state="A", bortho=1, ballo=None, bgaq=2) ** CELL_MEMB
        % Gaq(bpar=2, bgdp=3, bgbg=4) ** CELL_MEMB
        % GDP(b=3) ** CELL_MEMB
        % Gbg(b=4) ** CELL_MEMB
    )
    #    2AT + PAR2_I:Gaq:GDP:Gbg <---> TAT:PAR2_A:Gaq:GDP:Gbg
    Rule(
        "tat_bind_PAR2_pre",
        TAT(b=None) ** EXTRACELLULAR + PAR2_i_Gaq_gdp_Gbg | tat_PAR2_a_Gaq_gdp_Gbg,
        kf_PAR2_bind_TAT,
        kr_PAR2_bind_TAT,
    )

    return


def observables():
    """Defines observables for the receptor states.

    Defines 4 observables and 3 expressions.

    Observables:
        iPAR2 = inactive PAR2
        aPAR2 = active PAR2
        dPAR2 = denatured PAR2
        occupied_PAR2 = PAR2 with 2AT bound

    Expressions:
        totPAR2 = total amount of inactive and active PAR2
        par2_active_ratio = ratio of aPAR2 to totPAR2
        par2_occupancy = ratio of occupied_PAR2 to totPAR2 (receptor occupancy)
    """
    alias_model_components()
    # Inactive PAR2
    Observable("iPAR2", PAR2(state="I"))
    # Active PAR2
    Observable("aPAR2", PAR2(state="A"))
    # Denatured PAR2
    Observable("dPAR2", PAR2(state="D"))
    Observable("allPAR2", PAR2())
    alias_model_components()
    # Ro
    Expression("totPAR2", iPAR2 + aPAR2)
    try:
        tat_PAR2 = (
            TAT(b=1) ** EXTRACELLULAR
            % PAR2(bortho=1, ballo=WILD, bgaq=WILD) ** CELL_MEMB
        )
        Observable("occupied_PAR2", tat_PAR2)
        Expression("par2_occupancy", occupied_PAR2 / totPAR2)
    except:
        pass
    alias_model_components()
    Expression("par2_active_ratio", aPAR2 / totPAR2)

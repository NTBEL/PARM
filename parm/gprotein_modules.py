"""Defines functions for G-protein activation, regulation, and G-protein observables.

It includes four alternate mechanisms for G-protein activation:
    * classic_activation_mechanism
    * classic_activation_mechanism_with_constitutive_activity
    * precoupled_activation_mechanism
    * precoupled_activation_mechanism_with_constitutive_activity

As well as functions to add G-protein regulation:
    * gaq_hydrolyzes_gtp_to_gdp
    * rgs_enhances_gaq_hydrolosis_of_gtp_to_gdp

And an addon function for additional regulation of G-protein by PLC:
    * addon_plc_enhances_gaq_hydrolosis_of_gtp_to_gdp

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

# NumPy
import numpy as np

# Avogadro's Number from scipy
from scipy.constants import N_A
from sympy.functions.elementary.miscellaneous import Max

from . import defaults, units


def gprotein_monomers():
    """Declares the G-protein monomers.
    Adds 4 monomers and 3 annotations.

    Monomers:
        * Gaq - G-protein alpha subunit
        * Gbg - The beta-gamma G-protein subunit
        * GDP - guanosine diphosphate
        * GTP - guanosine triphosphate

    """
    # G-alpha_q G-protein unit, states: I = inactive, A = active
    Monomer("Gaq", ["bpar", "bgbg", "bgdp"])
    # G-beta-gamma G-protein units
    Monomer("Gbg", ["b"])
    # GDP
    Monomer("GDP", ["b"])
    # GTP
    Monomer("GTP", ["b"])
    alias_model_components()
    # Annotations
    # ===========
    Annotation(Gaq, "https://identifiers.org/uniprot:P50148")
    Annotation(GDP, "https://identifiers.org/CHEBI:17552")
    Annotation(GTP, "https://identifiers.org/CHEBI:15996")


def gprotein_initials():
    """Declares initial conditions for G-protein monomers.

    G-proteins are initialized as the heterotrimer consisting of the alpha
    subunit Gaq and beta-gamma subunit Gbg with GDP bound to Gaq.

    Adds 3 parameters.

    Parameters:
        * Gprotein_0 - initial concentration of G-protein heterotrimer in
            the cell membrane. The Gaq subunit is initially bound to GDP.
        * GTP_0 - initial concentration of free GTP in the cytosol.
        * GDP_0 - initial concentration of free GDP in the cytosol.
    """
    alias_model_components()
    # From Falkenburger et al. 2010 https://dx.doi.org/10.1085%2Fjgp.200910344
    #     - Endogenous G-protein density: 40/micrometer^2
    #     - Overexpressed G-protein density: 3,000/micrometer^2
    # Brinkerhoff et al. 2008: G-protein concentration 1e4 /cell
    Parameter("Gprotein_0", 1e4)
    # GTP
    # Physiolocal concentration of GTP in mammalian cells is generally 468 +/- 224 microM
    # and for human cells it is 305 microM
    # as per Traut https://doi.org/10.1007/bf00928361
    Parameter("GTP_0", 305 * units.microM_to_molec_per_pL * Vcyto.value)
    # GDP
    # Physiolocal concentration of GDP in human cells is 36 microM
    # as per Traut https://doi.org/10.1007/BF00928361
    Parameter("GDP_0", 36 * units.microM_to_molec_per_pL * Vcyto.value)
    alias_model_components()
    # inactive G-protein heterotrimer Gaq-GDP:Gbg
    # Note that the beta and gamma portions are modeled as a single unit.
    # Alias the free Gprotein heterotrimer
    Gaq_gdp_Gbg = (
        Gaq(bpar=None, bgbg=3, bgdp=4) ** CELL_MEMB
        % GDP(b=3) ** CELL_MEMB
        % Gbg(b=4) ** CELL_MEMB
    )
    Initial(Gaq_gdp_Gbg, Gprotein_0)
    # GTP
    Initial(GTP(b=None) ** CYTOSOL, GTP_0)
    # GDP
    Initial(GDP(b=None) ** CYTOSOL, GDP_0)

    return


def gprotein_activation_by_2at_bound_active_par2():
    """Defines activation of G-protein by 2AT-bound active PAR2.

    This function encodes activation of Gaq by 2AT-occupied active PAR2 with
    reactions from PAR2 and G-protein binding to release of the active GTP-bound
    Gaq subunit. In all, there are 5 steps:
        i) G protein heterotrimer binds activated PAR2:
            TAT:PAR2_A + Gaq:GDP:Gbg <---> TAT:PAR2_A:Gaq:GDP:Gbg
        ii) GDP unbinds from Gaq:
            TAT:PAR2_A:Gaq:GDP:Gbc <---> TAT:PAR2_A:Gaq:Gbc + GDP
        iii) GTP binds Gaq:
            TAT:PAR2_A:Gaq:Gbc + GTP <---> TAT:PAR2_A:Gaq_A:GTP:Gbg
        iv) Gbg dissociates from Gaq (i.e., heterotrimer dissociation):
            TAT:PAR2_A:Gaq:GTP:Gbc <---> TAT:PAR2_A:Gaq:GTP + Gbc
        v) Gaq:GTP dissociates from PAR2:
            TAT:PAR2_A:Gaq:GTP <---> TAT:PAR2_A + Gaq:GTP

    Adds 10 parameters but no expressions.

    Parameters:
        * kf_TAT_PAR2_A_bind_Gaq - forward rate constant for PAR2 binding
            to Gaq in the heterotrimer.
                TAT:PAR2_A + Gaq:GDP:Gbg ---> TAT:PAR2_A:Gaq:GDP:Gbg, kf
        * kr_TAT_PAR2_A_bind_Gaq - reverse rate constant for PAR2 binding to
            Gaq in the heterotrimer.
                TAT:PAR2_A:Gaq:GDP:Gbg ---> TAT:PAR2_A + Gaq:GDP:Gbg, kr
        * k_gdp_release - rate constant for unbinding of GDP from Gaq after
            PAR2 binding.
                TAT:PAR2_A:Gaq:GDP:Gbg ---> TAT:PAR2_A:Gaq:Gbg + GDP, k_release
        * k_gdp_bind - rate constant for binding of GDP to Gaq after PAR2
            binding.
                TAT:PAR2_A:Gaq:Gbg + GDP ---> TAT:PAR2_A:Gaq:GDP:Gbg, k_bind
        * k_gtp_bind - rate constant for GTP binding to Gaq after PAR2 binding
            and the release of GDP.
                TAT:PAR2_A:Gaq:Gbg + GTP ---> TAT:PAR2_A:Gaq:GTP:Gbg, k_bind
        * k_gtp_release - rate constant for GTP unbinding from Gaq post PAR2
            binding and release of GDP.
                TAT:PAR2_A:Gaq:GTP:Gbg ---> TAT:PAR2_A:Gaq:Gbg + GTP, k_release
        * kf_heterotrimer_dissociation - forward rate constant for dissociation
            of the G-protein heterotrimer. i.e., Gbg unbinds from PAR2 bound
            Gaq-GTP.
        * kr_heterotrimer_dissociation - reverse rate constant for dissociation
            of the G-protein heterotrimer. i.e., Gbg re-binds to the PAR2 bound
            Gaq-GTP.
        * kf_gaq_dissociation - forward rate constant for dissociation of Gaq
            from PAR2. i.e., GTP bound Gaq unbinds from PAR2.
        * kr_gaq_dissociation - reverse rate constant for dissociation of Gaq
            from PAR2. i.e., GTP bound Gaq re-binds to PAR2.
    """
    alias_model_components()

    # Gaq binding activated-PAR2
    Parameter("kf_TAT_PAR2_A_bind_Gaq", defaults.KF_BIND / Vcyto.value)
    Parameter("kr_TAT_PAR2_A_bind_Gaq", defaults.KR_BIND)
    # Gaq release GDP
    Parameter("k_gdp_release", defaults.KR_BIND * 100)
    Parameter("k_gdp_bind", defaults.KF_BIND / Vcyto.value)
    # Gaq bind GTP
    Parameter("k_gtp_bind", defaults.KF_BIND / Vcyto.value)
    Parameter("k_gtp_release", defaults.KR_BIND / 10)
    # Gbg dissociates from Gaq
    Parameter("kf_heterotrimer_dissociation", defaults.KR_BIND * 100)
    Parameter("kr_heterotrimer_dissociation", defaults.KF_BIND / Vcyto.value)
    # Gaq:GTP dissociates from PAR2
    Parameter("kf_gaq_dissociation", defaults.KR_BIND * 100)
    Parameter("kr_gaq_dissociation", defaults.KF_BIND / Vcyto.value)
    alias_model_components()

    # Gaq activation by activated-PAR2:
    #    PAR2_A + Gaq_I <---> PAR2_A:Gaq_I <---> PAR2_A + Gaq_A
    tat_PAR2_a = (
        TAT(b=1) ** EXTRACELLULAR
        % PAR2(state="A", bortho=1, ballo=None, bgaq=None) ** CELL_MEMB
    )
    # Alias the free Gprotein heterotrimer
    Gaq_gdp_Gbg = (
        Gaq(bpar=None, bgbg=3, bgdp=4) ** CELL_MEMB
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
    # PAR2 bindings the G protein heterotrimer
    Rule(
        "par2_bind_gaq",
        tat_PAR2_a + Gaq_gdp_Gbg | tat_PAR2_a_Gaq_gdp_Gbg,
        kf_TAT_PAR2_A_bind_Gaq,
        kr_TAT_PAR2_A_bind_Gaq,
    )
    # Alias the complex  2AT:PAR2_A:Gaq:Gbg
    tat_PAR2_a_Gaq_Gbg = (
        TAT(b=1) ** EXTRACELLULAR
        % PAR2(state="A", bortho=1, ballo=None, bgaq=2) ** CELL_MEMB
        % Gaq(bpar=2, bgdp=None, bgbg=4) ** CELL_MEMB
        % Gbg(b=4) ** CELL_MEMB
    )
    # GDP unbinds from Gaq
    Rule(
        "gaq_releases_gdp",
        tat_PAR2_a_Gaq_gdp_Gbg | tat_PAR2_a_Gaq_Gbg + GDP(b=None) ** CYTOSOL,
        k_gdp_release,
        k_gdp_bind,
    )
    # Alias the complex 2AT:PAR2_A:Gaq:GTP:Gbg
    tat_PAR2_a_Gaq_gtp_Gbg = (
        TAT(b=1) ** EXTRACELLULAR
        % PAR2(state="A", bortho=1, ballo=None, bgaq=2) ** CELL_MEMB
        % Gaq(bpar=2, bgdp=3, bgbg=4) ** CELL_MEMB
        % GTP(b=3) ** CELL_MEMB
        % Gbg(b=4) ** CELL_MEMB
    )
    # GTP binds to Gaq
    Rule(
        "gaq_binds_gtp",
        tat_PAR2_a_Gaq_Gbg + GTP(b=None) ** CYTOSOL | tat_PAR2_a_Gaq_gtp_Gbg,
        k_gtp_bind,
        k_gtp_release,
    )
    # Alias the complex 2AT:PAR2_A:Gaq:GTP
    tat_PAR2_a_Gaq_gtp = (
        TAT(b=1) ** EXTRACELLULAR
        % PAR2(state="A", bortho=1, ballo=None, bgaq=2) ** CELL_MEMB
        % Gaq(bpar=2, bgdp=3, bgbg=None) ** CELL_MEMB
        % GTP(b=3) ** CELL_MEMB
    )
    # The Beta-Gamma G protein units unbind from Gaq
    Rule(
        "release_gbg",
        tat_PAR2_a_Gaq_gtp_Gbg | tat_PAR2_a_Gaq_gtp + Gbg(b=None) ** CELL_MEMB,
        kf_heterotrimer_dissociation,
        kr_heterotrimer_dissociation,
    )
    # Alias the complex Gaq:GTP
    Gaq_gtp = Gaq(bpar=None, bgdp=3, bgbg=None) ** CELL_MEMB % GTP(b=3) ** CELL_MEMB
    # Gaq unbinds from PAR2
    Rule(
        "release_gaq",
        tat_PAR2_a_Gaq_gtp | Gaq_gtp + tat_PAR2_a,
        kf_gaq_dissociation,
        kr_gaq_dissociation,
    )

    return


def heterotrimer_precouples_free_inactive_par2():
    """Defines precoupling of the G-protein heterotrimer to free inactive PAR2.

    Allows precouling of G-protein to PAR2 by reversible binding:
        PAR2_i + Gaq:GDP:Gbg <---> PAR2_i:Gaq:GDP:Gbg

    Adds 2 parameters and 1 expression.

    Parameters:
        Kd_gprotein_precoupling_i - dissociation constant for PAR2 and G-protein
            binding.
        kf_gprotein_precoupling_i - forward rate constant for PAR2 and G-protein
            binding.

    Expressions:
        kr_gprotein_precoupling_i - reverse rate constant for PAR2 and G-protein
            binding.
    """
    alias_model_components()
    Parameter("Kd_gprotein_precoupling_i", 100 * units.nM_to_molec_per_pL * Vcm.value)
    Parameter("kf_gprotein_precoupling_i", defaults.KF_BIND / Vcyto.value)
    alias_model_components()
    # Division by V_CM is take into account the compartment scaling for kf.
    Expression(
        "kr_gprotein_precoupling_i",
        Kd_gprotein_precoupling_i * kf_gprotein_precoupling_i / V_CM,
    )
    alias_model_components()
    # Alias the inactive PAR2
    PAR2_i = PAR2(state="I", bortho=None, ballo=None, bgaq=None) ** CELL_MEMB
    # Alias the free Gprotein heterotrimer
    Gaq_gdp_Gbg = (
        Gaq(bpar=None, bgbg=3, bgdp=4) ** CELL_MEMB
        % GDP(b=3) ** CELL_MEMB
        % Gbg(b=4) ** CELL_MEMB
    )
    # Alias the complex 2AT:PAR2_I:Gaq:GDP:Gbg
    PAR2_i_Gaq_gdp_Gbg = (
        PAR2(state="I", bortho=None, ballo=None, bgaq=2) ** CELL_MEMB
        % Gaq(bpar=2, bgdp=3, bgbg=4) ** CELL_MEMB
        % GDP(b=3) ** CELL_MEMB
        % Gbg(b=4) ** CELL_MEMB
    )
    Rule(
        "rule_grotein_precoupling",
        PAR2_i + Gaq_gdp_Gbg | PAR2_i_Gaq_gdp_Gbg,
        kf_gprotein_precoupling_i,
        kr_gprotein_precoupling_i,
    )

    return


def heterotrimer_binds_free_active_par2():
    """Defines binding of the G-protein heterotrimer to free active PAR2.

    Allows coupling of G-protein to active PAR2 by reversible binding:
        PAR2_A + Gaq:GDP:Gbg <---> PAR2_A:Gaq:GDP:Gbg

    Adds 2 parameters and 1 expression.

    Parameters:
        Kd_gprotein_binds_PAR2_A - dissociation constant for PAR2 and G-protein
            binding.
        kf_gprotein_binds_PAR2_A - forward rate constant for PAR2 and G-protein
            binding.

    Expressions:
        kr_gprotein_binds_PAR2_A - reverse rate constant for PAR2 and G-protein
            binding.
    """
    alias_model_components()
    Parameter("Kd_gprotein_binds_PAR2_A", 50 * units.nM_to_molec_per_pL * Vcm.value)
    Parameter("kf_gprotein_binds_PAR2_A", defaults.KF_BIND / Vcyto.value)
    alias_model_components()
    # Division by V_CM is take into account the compartment scaling for kf.
    Expression(
        "kr_gprotein_binds_PAR2_A",
        Kd_gprotein_precoupling_a * kf_gprotein_precoupling_a / V_CM,
    )
    alias_model_components()
    # Alias the inactive PAR2
    PAR2_a = PAR2(state="A", bortho=None, ballo=None, bgaq=None) ** CELL_MEMB
    # Alias the free Gprotein heterotrimer
    Gaq_gdp_Gbg = (
        Gaq(bpar=None, bgbg=3, bgdp=4) ** CELL_MEMB
        % GDP(b=3) ** CELL_MEMB
        % Gbg(b=4) ** CELL_MEMB
    )
    # Alias the complex 2AT:PAR2_A:Gaq:GDP:Gbg
    PAR2_a_Gaq_gdp_Gbg = (
        PAR2(state="A", bortho=None, ballo=None, bgaq=2) ** CELL_MEMB
        % Gaq(bpar=2, bgdp=3, bgbg=4) ** CELL_MEMB
        % GDP(b=3) ** CELL_MEMB
        % Gbg(b=4) ** CELL_MEMB
    )
    Rule(
        "rule_grotein_binds_PAR2_A",
        PAR2_i + Gaq_gdp_Gbg | PAR2_i_Gaq_gdp_Gbg,
        kf_gprotein_binds_PAR2_A,
        kr_gprotein_binds_PAR2_A,
    )

    return


def gprotein_activation_by_free_active_par2():
    """Defines activation of G-protein by 2AT-bound active PAR2.

    This function encodes activation of Gaq by 2AT-occupied active PAR2 with
    reactions from PAR2 and G-protein binding to release of the active GTP-bound
    Gaq subunit.

    Calls:
        * heterotrimer_binds_free_active_par2

    Adds 4 additional steps:
        i) GDP unbinds from Gaq:
            PAR2_A:Gaq:GDP:Gbc <---> PAR2_A:Gaq:Gbc + GDP
        ii) GTP binds Gaq:
            PAR2_A:Gaq:Gbc + GTP <---> PAR2_A:Gaq_A:GTP:Gbg
        iii) Gbg dissociates from Gaq (i.e., heterotrimer dissociation):
            PAR2_A:Gaq:GTP:Gbc <---> PAR2_A:Gaq:GTP + Gbc
        iv) Gaq:GTP dissociates from PAR2:
            PAR2_A:Gaq:GTP <---> PAR2_A + Gaq:GTP

    Adds 8 additional parameters but no expressions.

    Parameters:
        * k_gdp_release_free_PAR2_A - rate constant for unbinding of GDP from
            Gaq after PAR2 binding.
                PAR2_A:Gaq:GDP:Gbg ---> PAR2_A:Gaq:Gbg + GDP, k_release
        * k_gdp_bind_free_PAR2_A - rate constant for binding of GDP to Gaq after
            PAR2 binding.
                PAR2_A:Gaq:Gbg + GDP ---> PAR2_A:Gaq:GDP:Gbg, k_bind
        * k_gtp_bind_free_PAR2_A - rate constant for GTP binding to Gaq after
            PAR2 binding and the release of GDP.
                PAR2_A:Gaq:Gbg + GTP ---> PAR2_A:Gaq:GTP:Gbg, k_bind
        * k_gtp_release_free_PAR2_A - rate constant for GTP unbinding from Gaq
            post PAR2 binding and release of GDP.
                PAR2_A:Gaq:GTP:Gbg ---> PAR2_A:Gaq:Gbg + GTP, k_release
        * kf_heterotrimer_dissociation_free_PAR2_A - forward rate constant for
            dissociation of the G-protein heterotrimer. i.e., Gbg unbinds from
            PAR2 bound Gaq-GTP.
        * kr_heterotrimer_dissociation_free_PAR2_A - reverse rate constant for
            dissociation of the G-protein heterotrimer. i.e., Gbg re-binds to
            the PAR2 bound Gaq-GTP.
        * kf_gaq_dissociation_free_PAR2_A - forward rate constant for
            dissociation of Gaq from PAR2. i.e., GTP bound Gaq unbinds from PAR2.
        * kr_gaq_dissociation_free_PAR2_A - reverse rate constant for
            dissociation of Gaq from PAR2. i.e., GTP bound Gaq re-binds to PAR2.
    """
    heterotrimer_binds_free_active_par2()
    alias_model_components()
    # Gaq release GDP
    Parameter("k_gdp_release_free_PAR2_A", defaults.KR_BIND * 100)
    Parameter("k_gdp_bind_free_PAR2_A", defaults.KF_BIND)
    # Gaq bind GTP
    Parameter("k_gtp_bind_free_PAR2_A", defaults.KF_BIND / Vcyto.value)
    Parameter("k_gtp_release_free_PAR2_A", defaults.KR_BIND / 10)
    # Gbg dissociates from Gaq
    Parameter("kf_heterotrimer_dissociation_free_PAR2_A", defaults.KR_BIND * 100)
    Parameter(
        "kr_heterotrimer_dissociation_free_PAR2_A", defaults.KF_BIND / Vcyto.value
    )
    # Gaq:GTP dissociates from PAR2
    Parameter("kf_gaq_dissociation_free_PAR2_A", defaults.KR_BIND * 100)
    Parameter("kr_gaq_dissociation_free_PAR2_A", defaults.KF_BIND / Vcyto.value)
    alias_model_components()
    # Gaq activation by activated-PAR2:
    #    PAR2_A + Gaq_I <---> PAR2_A:Gaq_I <---> PAR2_A + Gaq_A
    # Alias the complex PAR2_A:Gaq:GDP:Gbg
    PAR2_a_Gaq_gdp_Gbg = (
        PAR2(state="A", bortho=None, ballo=None, bgaq=2) ** CELL_MEMB
        % Gaq(bpar=2, bgdp=3, bgbg=4) ** CELL_MEMB
        % GDP(b=3) ** CELL_MEMB
        % Gbg(b=4) ** CELL_MEMB
    )

    # Alias the complex  PAR2_A:Gaq:Gbg
    PAR2_a_Gaq_Gbg = (
        PAR2(state="A", bortho=None, ballo=None, bgaq=2) ** CELL_MEMB
        % Gaq(bpar=2, bgdp=None, bgbg=4) ** CELL_MEMB
        % Gbg(b=4) ** CELL_MEMB
    )
    # GDP unbinds from Gaq
    Rule(
        "gaq_releases_gdp_free_par2_a",
        PAR2_a_Gaq_gdp_Gbg | PAR2_a_Gaq_Gbg + GDP(b=None) ** CYTOSOL,
        k_gdp_release_free_PAR2_A,
        k_gdp_bind_free_PAR2_A,
    )
    # Alias the complex PAR2_A:Gaq:GTP:Gbg
    PAR2_a_Gaq_gtp_Gbg = (
        PAR2(state="A", bortho=None, ballo=None, bgaq=2) ** CELL_MEMB
        % Gaq(bpar=2, bgdp=3, bgbg=4) ** CELL_MEMB
        % GTP(b=3) ** CELL_MEMB
        % Gbg(b=4) ** CELL_MEMB
    )
    # GTP binds to Gaq
    Rule(
        "gaq_binds_gtp_free_par2_a",
        PAR2_a_Gaq_Gbg + GTP(b=None) ** CYTOSOL | PAR2_a_Gaq_gtp_Gbg,
        k_gtp_bind_free_PAR2_A,
        k_gtp_release_free_PAR2_A,
    )
    # Alias the complex PAR2_A:Gaq:GTP
    PAR2_a_Gaq_gtp = (
        PAR2(state="A", bortho=None, ballo=None, bgaq=2) ** CELL_MEMB
        % Gaq(bpar=2, bgdp=3, bgbg=None) ** CELL_MEMB
        % GTP(b=3) ** CELL_MEMB
    )
    # The Beta-Gamma G protein units unbind from Gaq
    Rule(
        "release_gbg_free_par2_a",
        PAR2_a_Gaq_gtp_Gbg | PAR2_a_Gaq_gtp + Gbg(b=None) ** CELL_MEMB,
        kf_heterotrimer_dissociation_free_PAR2_A,
        kr_heterotrimer_dissociation_free_PAR2_A,
    )
    # Alias the complex Gaq:GTP
    Gaq_gtp = Gaq(bpar=None, bgdp=3, bgbg=None) ** CELL_MEMB % GTP(b=3) ** CELL_MEMB
    # Alias the free active PAR2
    PAR2_a = PAR2(state="A", bortho=None, ballo=None, bgaq=None) ** CELL_MEMB
    # Gaq unbinds from PAR2
    Rule(
        "release_gaq_free_par2_a",
        PAR2_a_Gaq_gtp | Gaq_gtp + PAR2_a,
        kf_gaq_dissociation_free_PAR2_A,
        kr_gaq_dissociation_free_PAR2_A,
    )
    return

def single_step_catalytic_gprotein_activation_by_par2():
    """Defines catalytic activation of the G-protein by agonist-bound active PAR2.

    This function encodes a catalytic G-protein activation mechanism like that
    used by Yi et al. PNAS 2003 https://doi.org/10.1073/pnas.1834247100 in their
    yeast G-protein cycle model whereby the G-protein heterotrimer binds and Gaq
    gets catalytically activated by the receptor (PAR2, in this case) in a
    single 2nd-order irreversible binding and conversion reaction:
        2AT:PAR2_A + Gaq:GDP:Gbg ---> 2AT:PAR2_A + Gaq:GTP + Gbg

    Adds 1 parameter and 1 annotation.

    Parameters:
        kf_PAR2_activate_Gprotein - binding rate constant controlling the
            catalytic activation of Gaq by PAR2.

    """

    alias_model_components()
    # The corresponding parameter in Yi et al. 2003 is:
    #     k_Ga = 1x10^-5  1/(molecule/cell * s)
    # which was obtained from from fitting their model to FRET data.
    # We can use that value as our nominal value here.
    # Since the binding rate is already molecules per cell we can multiply by
    # the V_CM value here to negate the compartment scaling applied later to the
    # reaction.
    Parameter('kf_PAR2_activate_Gprotein', 1e-5 * V_CM.get_value())
    alias_model_components()
    tat_PAR2_a = PAR2(state='A', bortho=ANY, bgaq=None)**CELL_MEMB
    # Alias the free Gprotein heterotrimer
    Gaq_gdp_Gbg = Gaq(bpar=None, bgbg=3, bgdp=4)**CELL_MEMB % GDP(b=3)**CELL_MEMB % Gbg(b=4)**CELL_MEMB
    # Alias the complex Gaq:GTP
    Gaq_gtp = (Gaq(bpar=None, bgdp=3, bgbg=None)**CELL_MEMB % GTP(b=3)**CELL_MEMB)
    # Define the reaction rule.
    Rule('par2_single_step_activate_Gprotein', tat_PAR2_a + Gaq_gdp_Gbg >>
                              tat_PAR2_a + Gaq_gtp + Gbg(b=None)**CELL_MEMB,
                              kf_PAR2_activate_Gprotein)

    Annotation(par2_activate_Gprotein,
               'https://identifiers.org/doi:10.1073/pnas.1834247100',
                predicate='isDerivedFrom')
    return

def gaq_hydrolyzes_gtp_to_gdp():
    """Defines hydrolosis of GTP to GDP when bound to free Gaq.

    Encodes the 1st order hydrolosis of GTP to GDP by Gaq:
        Gaq:GTP ---> Gaq:GDP

    Adds 1 parameter.

    Parameters:
            * k_gtp_to_gdp_auto - 1st order rate constant for the conversion.
    """
    # Hydrolosis of GTP bound to Gaq
    # 1. Autocatalysis rate for Gaq is ~0.8 1/min = 0.0133 1/s
    # Bernstein et al. https://doi.org/10.1016/0092-8674(92)90165-9
    # Also see Sprang https://dx.doi.org/10.1002%2Fbip.22836
    # Note that the corresponding value from Yi et al. PNAS 2003
    # https://doi.org/10.1073/pnas.1834247100 is:
    #    k_Gd0 = 0.004 1/s
    Parameter("k_gtp_to_gdp_auto", 1.33e-2)
    alias_model_components()
    # Alias the complex Gaq:GTP
    Gaq_gtp = Gaq(bpar=None, bgdp=3, bgbg=None) ** CELL_MEMB % GTP(b=3) ** CELL_MEMB
    # Alias the complex Gaq:GDP
    Gaq_gdp = Gaq(bpar=None, bgdp=3, bgbg=None) ** CELL_MEMB % GDP(b=3) ** CELL_MEMB
    # Gaq can (slowly) hydolyze GTP to GDP
    Rule("gtp_hydrolosis_auto", Gaq_gtp >> Gaq_gdp, k_gtp_to_gdp_auto)
    return


def rgs_enhances_gaq_hydrolosis_of_gtp_to_gdp():
    """Defines a regulator of g-protein signaling that promotes hydrolosis of GTP.

    Adds a an enahanced catlalytic conversion of GTP to GDP by Gaq when RGS
    protein is bound to Gaq:
        i. Gaq:GTP + RGS <---> Gaq:GTP:RGS
        ii. Gaq:GTP:RGS ---> Gaq:GDP + RGS

    Adds 1 monomer, 4 parameters, 1 annotation.
    Also defines the initial condition for the RGS protein.

    Monomers:
        * RGS - regulator of G-protein signaling (e.g., RGS4)

    Parameters:
        * RGS_0 - initial concentration of RGS protein in the cytosol.
        * kf_rgs_bind_gaq - forward rate constant for binding of RGS to Gaq.
        * kr_rgs_bind_gaq - reverse rate constant for binding of RGS to Gaq.
        * k_gtp_to_gdp_rgs - 1st order catalytic rate for the conversion of
            GTP to GDP by RGS-bound Gaq.
    """
    # Regulator of G protein Signaling (RGS)
    Monomer("RGS", ["b"])
    alias_model_components()
    # RGS
    # For RAW 264.7 Cell model between 0.008 and 0.012 microM as
    # per Maurya and Subramaniam 2007 https://doi.org/10.1529/biophysj.106.097469
    # Assume 0.010 microM is a reasonable starting point.
    Parameter("RGS_0", 0.010 * units.microM_to_molec_per_pL * Vcyto.value)
    # 2. RGS binding and enhanced conversion of GTP to GDP
    Parameter("kf_rgs_bind_gaq", defaults.KF_BIND / Vcyto.value)
    Parameter("kr_rgs_bind_gaq", defaults.KR_BIND)
    # Enhanced catalysis rate is for Gaq catalysis by RGS4
    # is 100x higher than autocatalsis as per
    # Chidiac and Ross https://doi.org/10.1074/jbc.274.28.19639
    # Also see Sprang https://dx.doi.org/10.1002%2Fbip.22836
    Parameter("k_gtp_to_gdp_rgs", 1.33e-2 * 100)
    alias_model_components()

    Initial(RGS(b=None) ** CYTOSOL, RGS_0)

    # Alias the complex Gaq:GTP
    Gaq_gtp = Gaq(bpar=None, bgdp=3, bgbg=None) ** CELL_MEMB % GTP(b=3) ** CELL_MEMB
    # Alias the complex Gaq:GDP
    Gaq_gdp = Gaq(bpar=None, bgdp=3, bgbg=None) ** CELL_MEMB % GDP(b=3) ** CELL_MEMB
    # Alias the complex Gaq:GTP:RGS
    Gaq_gtp_RGS = (
        Gaq(bpar=None, bgdp=3, bgbg=1) ** CELL_MEMB
        % GTP(b=3) ** CELL_MEMB
        % RGS(b=1) ** CYTOSOL
    )
    # Gaq:GTP binds to RGS protein
    Rule(
        "gaq_gtp_binds_rgs",
        Gaq_gtp + RGS(b=None) ** CYTOSOL | Gaq_gtp_RGS,
        kf_rgs_bind_gaq,
        kr_rgs_bind_gaq,
    )
    # Binding of RGS protein promotes (faster) hydrolysis of GTP to GDP
    Rule(
        "gtp_hydrolosis_rgs",
        Gaq_gtp_RGS >> Gaq_gdp + RGS(b=None) ** CYTOSOL,
        k_gtp_to_gdp_rgs,
    )

    Annotation(RGS, "https://identifiers.org/uniprot:P49798", predicate="hasVersion")
    return

def rgs_enhances_gaq_hydrolosis_of_gtp_to_gdp_implicit_single_step():
    """Defines an implicit regulator of g-protein signaling that promotes hydrolosis of GTP.

    Adds an enahanced catlalytic conversion of GTP to GDP by Gaq promoted by an
    implicit RGS protein using a single-step 1st-order irreversible reaction
    like that in the heterotrimeric G protein cycle of Yi et al. PNAS 2003
    https://doi.org/10.1073/pnas.1834247100:
        Gaq:GTP --RGS--> Gaq:GDP

    Adds 1 parameter.

    Parameters:
        * k_gtp_to_gdp_rgs - 1st-order rate constant for the conversion.

    """

    alias_model_components()
    # The corresponding parameter in Yi et al. 2003 is:
    #     k_Gd1 = 0.11 1/s
    # We can use that value as our nominal value here.
    Parameter("k_gtp_to_gdp_rgs", 0.11)
    alias_model_components()
    # Alias the complex Gaq:GTP
    Gaq_gtp = Gaq(bpar=None, bgdp=3, bgbg=None) ** CELL_MEMB % GTP(b=3) ** CELL_MEMB
    # Alias the complex Gaq:GDP
    Gaq_gdp = Gaq(bpar=None, bgdp=3, bgbg=None) ** CELL_MEMB % GDP(b=3) ** CELL_MEMB
    # Gaq can hydolyze GTP to GDP faster due to effect of RGS
    Rule("gtp_hydrolosis_rgs", Gaq_gtp >> Gaq_gdp, k_gtp_to_gdp_rgs)
    return

def reversible_heterotrimer_reassociation():
    """Defines reassociation of the G-protein heterotrimer by an reversible binding reaction.

    This function encodes a reversible G-protein reassociation mechanism whereby
    the inactive GDP-bound alhpa subunit (Gaq) reforms the G-protein
    heterotrimer with the beta-gamma subunits (Gbg) via reversible binding:
        Gaq:GDP + Gbg <---> Gaq:GDP:Gbq

    Adds 2 parameters.

    Parameters:
        kf_reversible_heterotrimer_reassociation - forward binding rate constant
            for reformation of the G-protein heterotrimer.
        kr_reversible_heterotrimer_reassociation - reverse binding rate constant
            for reformation of the G-protein heterotrimer.

    """

    alias_model_components()
    # The corresponding parameter in Yi et al. 2003 is:
    #     k_G1 = 1 / (molecule/cell * s)
    # For kf we can use set the nominal value to the irreversible binding
    # rate constant from Yi et al. PNAS 2003
    # https://doi.org/10.1073/pnas.1834247100 (yeast G-protein cycle for SSTR),
    # which is:
    #     k_G1 = 1 / (molecule/cell * s)
    # Since the binding rate is already molecules per cell we can multiply by
    # the V_CM value here to negate the compartment scaling applied later to the
    # reaction.
    Parameter("kf_reversible_heterotrimer_reassociation", 1 * V_CM.get_value())
    Parameter("kr_reversible_heterotrimer_reassociation", defaults.KR_BIND)
    alias_model_components()
    Gaq_gdp = Gaq(bpar=None, bgdp=3) ** CELL_MEMB % GDP(b=3) ** CELL_MEMB
    bind_complex(
        Gaq_gdp, 'bgbg',
        Gbg()**CELL_MEMB, 'b',
        [kf_reversible_heterotrimer_reassociation,
         kr_reversible_heterotrimer_reassociation],
    )
    return

def irreversible_heterotrimer_reassociation():
    """Defines reassociation of the G-protein heterotrimer by an irreversible binding reaction.

    This function encodes an irreversible G-protein reassociation mechanism like
    that used by Yi et al. PNAS 2003 https://doi.org/10.1073/pnas.1834247100 in
    their yeast G-protein cycle model (SST receptor) whereby the inactive GDP-bound alhpa
    subunit (Gaq) reforms the G-protein heterotrimer with the beta-gamma
    subunits (Gbg) via a 2nd-order irreversible binding reaction:
        Gaq:GDP + Gbg ---> Gaq:GDP:Gbq

    Adds 1 parameter.

    Parameters:
        k_irrev_heterotrimer_reassociation - binding rate constant describing
            the irreversible reformation of the G-protein heterotrimer.

    """

    alias_model_components()
    # The corresponding parameter in Yi et al. 2003 is:
    #     k_G1 = 1 / (molecule/cell * s)
    # We can use that value as our nominal value here.
    # Since the binding rate is already molecules per cell we can multiply by
    # the V_CM value here to negate the compartment scaling applied later to the
    # reaction.
    Parameter("kf_irrev_heterotrimer_reassociation", 1 * V_CM.get_value())
    alias_model_components()
    # Alias the complex Gaq:GDP
    Gaq_gdp = Gaq(bpar=None, bgdp=3, bgbg=None) ** CELL_MEMB % GDP(b=3) ** CELL_MEMB
    # Alias the free Gprotein heterotrimer
    Gaq_gdp_Gbg = (
        Gaq(bpar=None, bgbg=3, bgdp=4) ** CELL_MEMB
        % GDP(b=3) ** CELL_MEMB
        % Gbg(b=4) ** CELL_MEMB
    )

    Rule(
        'irrev_heterotrimer_reassociation',
         Gaq_gdp + Gbg(b=None)**CELL_MEMB >> Gaq_gdp_Gbg,
         kf_irrev_heterotrimer_reassociation,
    )
    return

def classic_activation_mechanism():
    """Defines the classic mechanism of G-protein activation.

    In the classic activation mechanism G-proteins only interact with and get
    actived by active PAR2 without any precoupling to the inactive receptor.

    Calls:
        * gprotein_monomers
        * gprotein_initials
        * gprotein_activation_by_2at_bound_active_par2
        * reversible_heterotrimer_reassociation
    """
    gprotein_monomers()
    gprotein_initials()
    gprotein_activation_by_2at_bound_active_par2()
    reversible_heterotrimer_reassociation()
    return


def classic_activation_mechanism_with_constitutive_activity():
    """Defines the classic mechanism of G-protein activation but adds constitutive PAR2 activity.

    In the classic activation mechanism G-proteins only interact with and get
    actived by active PAR2 without any precoupling to the inactive receptor.

    Calls:
        * classic_activation_mechanism
        * gprotein_activation_by_free_active_par2
    """
    classic_activation_mechanism()
    gprotein_activation_by_free_active_par2()
    return


def precoupled_activation_mechanism():
    """Defines a mechanism of G-protein activation that includes receptor precoupling.

    In the precoupled activation mechanism some portion of G-proteins can bind
    to inactive PAR2. This precoupling of G-protein to inactive PAR2 can
    facilitate faster G-protein activation once a receptor activating stimulus
    is added.

    Calls:
        * gprotein_monomers
        * gprotein_initials
        * heterotrimer_precouples_free_inactive_par2
        * gprotein_activation_by_2at_bound_active_par2
        * reversible_heterotrimer_reassociation
    """
    gprotein_monomers()
    gprotein_initials()
    heterotrimer_precouples_free_inactive_par2()
    gprotein_activation_by_2at_bound_active_par2()
    reversible_heterotrimer_reassociation()
    return


def precoupled_activation_mechanism_with_constitutive_activity():
    """Defines a mechanism of G-protein activation that includes receptor precoupling and constitutive PAR2 activity.

    In the precoupled activation mechanism some portion of G-proteins can bind
    to inactive PAR2. This precoupling of G-protein to inactive PAR2 can
    facilitate faster G-protein activation once a receptor activating stimulus
    is added. This mechanism also incorporates constitutive activity of PAR2.

    Calls:
        * gprotein_monomers
        * gprotein_initials
        * heterotrimer_precouples_free_inactive_par2
        * gprotein_activation_by_2at_bound_active_par2
        * gprotein_activation_by_free_active_par2
        * reversible_heterotrimer_reassociation
    """
    gprotein_monomers()
    gprotein_initials()
    heterotrimer_precouples_free_inactive_par2()
    gprotein_activation_by_2at_bound_active_par2()
    gprotein_activation_by_free_active_par2()
    reversible_heterotrimer_reassociation()
    return


def yi2003_heterotrimeric_gprotein_cycle():
    """Defines a mechanism of G-protein activation consistent with that of Yi et al. 2003.

    This function encodes the G protein activation portion of the heterotrimeric
    G protein cycle described by Yi et al. PNAS 2003
    https://doi.org/10.1073/pnas.1834247100 adpated for PAR2 as the receptor.
    The nominal parameters also match those reported in their paper, which are
    based on yeast.

    Note that the yeast heterotrimer G protein cycle would be a type of classic
    activation mechanism where only the ligand-bound receptor interacts with and
    activates G proteins. Also note that this mechanism doesn't use free GTP or
    GDP in cytosol (although they initialized with the gprotein_initials
    function).

    Calls:
        * gprotein_monomers
        * gprotein_initials
        * single_step_catalytic_gprotein_activation_by_par2
        * gaq_hydrolyzes_gtp_to_gdp
        * irreversible_heterotrimer_reassociation
        * rgs_enhances_gaq_hydrolosis_of_gtp_to_gdp_implicit_single_step
    """
    
    gprotein_monomers()
    gprotein_initials()
    single_step_catalytic_gprotein_activation_by_par2()
    gaq_hydrolyzes_gtp_to_gdp()
    alias_model_components()
    # Update the nominal value of the autocatalytic conversion of GTP to GDP
    # to match the value in Yi et al.:
    #    k_Gd0 = 0.004 1/s
    k_gtp_to_gdp_auto.value = 0.004 # 1/s
    irreversible_heterotrimer_reassociation()
    rgs_enhances_gaq_hydrolosis_of_gtp_to_gdp_implicit_single_step()
    return

def addon_plc_enhances_gaq_hydrolosis_of_gtp_to_gdp():
    """Adds enhanced hydrolosis of GTP to GDP by Gaq when bound to PLC.
    This addon function uses PLC so must be called after
    calcium signaling modules that define the PLC monomer.
    It adds one additional rule for 1st order catalytic hydrolosis of
    GTP to GDP by Gaq when Gaq is bound to PLC (binding to PLC is part
    of the calcium signaling modules):
        Gaq:GTP:PLC ---> Gaq:GDP + PLC

    Adds 1 additional parameter.

    Parameters:
        * k_gtp_to_gdp_plc -
    """
    # 3. PLC binding enhanced conversion of GTP to GDP
    Parameter("k_gtp_to_gdp_plc", 1.33e-2 * 10)
    alias_model_components()
    Gaq_gtp_PLC = (
        Gaq(bpar=None, bgdp=3, bgbg=1) ** CELL_MEMB
        % GTP(b=3) ** CELL_MEMB
        % PLC(bgaq=1, bpip2=WILD, bca=WILD) ** CELL_MEMB
    )
    # Alias the complex Gaq:GDP
    Gaq_gdp = Gaq(bpar=None, bgdp=3, bgbg=None) ** CELL_MEMB % GDP(b=3) ** CELL_MEMB
    # Enhanced hydrolosis of GTP when Gaq is bound to PLC
    #   Gaq:GTP:PLC ---> Gaq:GDP + PLC
    Rule(
        "gtp_hydrolosis_plc",
        Gaq_gtp_PLC >> Gaq_gdp + PLC(bgaq=None, bpip2=WILD, bca=WILD) ** CYTOSOL,
        k_gtp_to_gdp_plc,
    )
    return


def observables():
    """Defines observables for the G-protein states.

    Defines 2 observables and 1 expression.

    Observables:
        totGaq = total amount of Gaq
        aGaq = amount of active Gaq (bound to GTP)

    Expressions:
        active_Gaq_ratio = ratio of aGaq to totGaq
    """
    alias_model_components()
    Observable("totGaq", Gaq())
    Observable(
        "aGaq", Gaq(bpar=WILD, bgdp=3, bgbg=WILD) ** CELL_MEMB % GTP(b=3) ** CELL_MEMB
    )
    alias_model_components()
    Expression("active_Gaq_ratio", aGaq / totGaq)
    return

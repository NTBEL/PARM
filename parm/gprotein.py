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
    alias_model_components()

    # Gaq binding activated-PAR2
    Parameter("kf_TAT_PAR2_A_bind_Gaq", defaults.KF_BIND / Vcyto.value)
    Parameter("kr_TAT_PAR2_A_bind_Gaq", defaults.KR_BIND)
    # Gaq release GDP
    Parameter("k_gdp_release", defaults.KR_BIND * 100)
    Parameter("k_gdp_bind", defaults.KF_BIND)
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
        TAT(b=1) ** EXTRACELLULAR % PAR2(state="A", bortho=1, bgaq=None) ** CELL_MEMB
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
        % PAR2(state="A", bortho=1, bgaq=2) ** CELL_MEMB
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
        % PAR2(state="A", bortho=1, bgaq=2) ** CELL_MEMB
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
        % PAR2(state="A", bortho=1, bgaq=2) ** CELL_MEMB
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
        % PAR2(state="A", bortho=1, bgaq=2) ** CELL_MEMB
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
    PAR2_i = PAR2(state="I", bortho=None, bgaq=None) ** CELL_MEMB
    # Alias the free Gprotein heterotrimer
    Gaq_gdp_Gbg = (
        Gaq(bpar=None, bgbg=3, bgdp=4) ** CELL_MEMB
        % GDP(b=3) ** CELL_MEMB
        % Gbg(b=4) ** CELL_MEMB
    )
    # Alias the complex 2AT:PAR2_I:Gaq:GDP:Gbg
    PAR2_i_Gaq_gdp_Gbg = (
        PAR2(state="I", bortho=None, bgaq=2) ** CELL_MEMB
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
    PAR2_a = PAR2(state="A", bortho=None, bgaq=None) ** CELL_MEMB
    # Alias the free Gprotein heterotrimer
    Gaq_gdp_Gbg = (
        Gaq(bpar=None, bgbg=3, bgdp=4) ** CELL_MEMB
        % GDP(b=3) ** CELL_MEMB
        % Gbg(b=4) ** CELL_MEMB
    )
    # Alias the complex 2AT:PAR2_A:Gaq:GDP:Gbg
    PAR2_a_Gaq_gdp_Gbg = (
        PAR2(state="A", bortho=None, bgaq=2) ** CELL_MEMB
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
        PAR2(state="A", bortho=None, bgaq=2) ** CELL_MEMB
        % Gaq(bpar=2, bgdp=3, bgbg=4) ** CELL_MEMB
        % GDP(b=3) ** CELL_MEMB
        % Gbg(b=4) ** CELL_MEMB
    )

    # Alias the complex  PAR2_A:Gaq:Gbg
    PAR2_a_Gaq_Gbg = (
        PAR2(state="A", bortho=None, bgaq=2) ** CELL_MEMB
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
        PAR2(state="A", bortho=None, bgaq=2) ** CELL_MEMB
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
        PAR2(state="A", bortho=None, bgaq=2) ** CELL_MEMB
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
    PAR2_a = PAR2(state="A", bortho=None, bgaq=None) ** CELL_MEMB
    # Gaq unbinds from PAR2
    Rule(
        "release_gaq_free_par2_a",
        PAR2_a_Gaq_gtp | Gaq_gtp + PAR2_a,
        kf_gaq_dissociation_free_PAR2_A,
        kr_gaq_dissociation_free_PAR2_A,
    )
    return


def gaq_hydrolyzes_gtp_to_gdp():
    # Hydrolosis of GTP bound to Gaq
    # 1. Autocatalysis rate for Gaq is ~0.8 1/min = 0.0133 1/s
    # Bernstein et al. https://doi.org/10.1016/0092-8674(92)90165-9
    # Also see Sprang https://dx.doi.org/10.1002%2Fbip.22836
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


def classic_activation_mechanism():
    gprotein_monomers()
    gprotein_initials()
    gprotein_activation_by_2at_bound_active_par2()
    return


def classic_activation_mechanism_with_constitutive_activity():
    classic_activation_mechanism()
    gprotein_activation_by_free_active_par2()
    return


def precoupled_activation_mechanism():
    gprotein_monomers()
    gprotein_initials()
    heterotrimer_precouples_free_inactive_par2()
    gprotein_activation_by_2at_bound_active_par2()
    return


def precoupled_activation_mechanism_with_constitutive_activity():
    gprotein_monomers()
    gprotein_initials()
    heterotrimer_precouples_free_inactive_par2()
    gprotein_activation_by_2at_bound_active_par2()
    gprotein_activation_by_free_active_par2()
    return


def addon_plc_enhances_gaq_hydrolosis_of_gtp_to_gdp():
    # 3. PLC binding enhanced conversion of GTP to GDP
    Parameter("k_gtp_to_gdp_plc", 1.33e-2 * 10)
    alias_model_components()
    Gaq_gtp_PLC = (
        Gaq(bpar=None, bgdp=3, bgbg=1) ** CELL_MEMB
        % GTP(b=3) ** CELL_MEMB
        % PLC(bgaq=1) ** CELL_MEMB
    )
    # Alias the complex Gaq:GDP
    Gaq_gdp = Gaq(bpar=None, bgdp=3, bgbg=None) ** CELL_MEMB % GDP(b=3) ** CELL_MEMB
    # Enhanced hydrolosis of GTP when Gaq is bound to PLC
    #   Gaq:GTP:PLC ---> Gaq:GDP + PLC
    Rule(
        "gtp_hydrolosis_plc",
        Gaq_gtp_PLC >> Gaq_gdp + PLC(bgaq=None, bpip2=None) ** CYTOSOL,
        k_gtp_to_gdp_plc,
    )
    return


def observables():
    alias_model_components()
    Observable("totGaq", Gaq())
    Observable(
        "aGaq", Gaq(bpar=WILD, bgdp=3, bgbg=WILD) ** CELL_MEMB % GTP(b=3) ** CELL_MEMB
    )
    alias_model_components()
    Expression("active_Gaq_ratio", aGaq / totGaq)
    return

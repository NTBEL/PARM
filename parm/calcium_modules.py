"""Defines functions for mechanisms in the calcium signaling pathway.

The full pathway including monomer and initial definitions is encoded by the
function:
    * gaq_activated_calcium_signaling

Measurement of the FRET ratio for TN-XXL sensing of free cytosolic Ca2+ is
encoded by the function:
    * fret_calcium_indicator_tnxxl

Other functions corresponding individual mechansitic elements or combinations
therof from the calcium signaling pathway are:
    * plc_binds_gaq_and_catalyzes_pip2_to_ip3
    * ip3_binds_ip3r
    * ip3r_transports_er_calcium_to_cytosol
    * ip3_degradation
    * cytosolic_calcium_inhibits_ip3r
    * calcium_binds_plc_and_enhances_pip2_hydrolysis
    * cytosolic_calcium_feedback (combines previous two functions)
    * cytosolic_calcium_buffering
    * calcium_extrusion_and_influx
    * regulation_of_cytosolic_calcium_concentration (combines previous two functions)

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
from sympy import Piecewise
import numpy as np

from . import defaults, units


def calcium_signal_monomers():
    """Declares the monomers that take part in intracellular PLC and IP3 calcium signaling pathway.

    Adds 6 monomers and 6 annotations.

    Monomers:
        * PLC - phospholipase C
        * PIP2 - Phosphatidylinositol bisphosphate
        * IP3 - Phosphatidylinositol trisphosphate
        * DAG - diacylglycerol
        * IP3R - IP3 receptor
        * Ca - calcium 2+
    """
    # Phospholipase C
    #   bgaq = binding site for Gaq
    #   bpip2 = binding site for PIP2
    #   bca = binding site for cytosolic calcium.
    Monomer("PLC", ["bgaq", "bpip2", "bca"])
    # PIP2
    Monomer("PIP2", ["b"])
    # IP3
    Monomer("IP3", ["b"])
    # DAG
    Monomer("DAG", ["b"])
    # IP3 receptor
    Monomer("IP3R", ["b1", "b2", "b3", "b4", "bcaer", "bcacyt"])
    # Monomer("IP3R_subunit", ["sub1", "sub2", "bip3", "bcaer", "bcacyt"])
    # Calcium 2+, loc: E = ER space, C = cytosol
    Monomer("Ca", ["b"])
    alias_model_components()

    # Annotations
    # ===========
    Annotation(PLC, "https://identifiers.org/uniprot:Q9NQ66", predicate="hasVersion")
    Annotation(PIP2, "https://identifiers.org/CHEBI:145879", predicate="hasVersion")
    Annotation(IP3, "https://identifiers.org/CHEBI:26034")
    Annotation(DAG, "https://identifiers.org/CHEBI:85722", predicate="hasVersion")
    Annotation(IP3R, "https://identifiers.org/uniprot:Q14643")
    Annotation(Ca, "https://identifiers.org/CHEBI:29108")
    return


def calcium_signal_initials():
    """Declares initial conditions for calcium signaling monomers.

    Adds 8 parameters.

    Parameters:
        * PLC_0 - initial amount of PLC in the cell membrane.
        * PIP2_0 - initial amount of PIP2 in the cell membrane.
        * IP3_0 - initial amount of IP3 in the cytosol.
        * DAG_0 - initial amount of DAG in the cell membrane.
        * IP3R_0 - initial amount of IP3R in the ER membrane.
        * Ca_E_0 - initial amount of Ca2+ in the ER lumen.
        * Ca_C_0 - initial amount of Ca2+ in the cytosol.
        * Ca_extra_0 - initial amount of Ca2+ in the extracellular space.
    """
    alias_model_components()
    # inactive PLC
    # From Falkenburger et al. 2010 https://dx.doi.org/10.1085%2Fjgp.200910344
    # tsA201 cells
    # Endogenous PLCB1 concentration: 3/micrometer^2
    # Overexpressed PLCB1 concentration: 3,000/micrometer^2
    # PLC total endogenous: 10/micrometer^2
    Parameter("PLC_0", 10 * SAcell.value)
    # PIP2
    # Basal no. of PIP2 molecules is 49997 as per Lemon et al. 2003 https://doi.org/10.1016/S0022-5193(03)00079-1
    # also free PIP2 of 5000 per micrometer^2 used by Falkenburger et al. 2013 https://doi.org/10.1085/jgp.201210887
    Parameter("PIP2_0", 5000 * SAcell.value)
    # Assume there is no IP3 or DAG initially.
    Parameter("IP3_0", 0)
    Parameter("DAG_0", 0)
    # IP3R
    # Set nominal as 1/micrometer^2
    Parameter("IP3R_0", 1 * SAer.value)
    # ER Ca2+ store
    # ER lumen of HEK-293 cells has between roughly 400-600 microM with an average
    # around around 530 +- 70 uM.
    # Foyouzi-Youssefi et al. https://doi.org/10.1073/pnas.97.11.5723 (Fig. 3C, control)
    Parameter("Ca_E_0", 530 * units.microM_to_molec_per_pL * Ver.value)
    # Initial concentration of Ca2+ in the cytosol expected to be around 100 nM.
    # 97 +- 5 nM from Tong et al. JBC 1999
    Parameter("Ca_C_0", 97 * units.nM_to_molec_per_pL * Vcyto.value)
    # In MH experiments the extracellular space is filled with ACSF with 3.1 mM
    # CaCl2, so extracellular Ca2+ should be around 3.1 mM.
    # Kang et al. https://doi.org/10.1021/acsnano.9b01993
    Parameter("Ca_extra_0", 3.1 * 1e3 * units.microM_to_molec_per_pL * Vextra.value)
    alias_model_components()

    Initial(PLC(bgaq=None, bpip2=None, bca=None) ** CELL_MEMB, PLC_0)
    Initial(PIP2(b=None) ** CELL_MEMB, PIP2_0)
    Initial(
        IP3R(b1=None, b2=None, b3=None, b4=None, bcaer=None, bcacyt=None) ** ER_MEMB,
        IP3R_0,
    )
    Initial(IP3(b=None) ** CYTOSOL, IP3_0)
    Initial(DAG(b=None) ** CELL_MEMB, DAG_0)
    Initial(Ca(b=None) ** ER_LUMEN, Ca_E_0)
    Initial(Ca(b=None) ** CYTOSOL, Ca_C_0)
    Initial(Ca(b=None) ** EXTRACELLULAR, Ca_extra_0)
    return


def plc_binds_gaq_and_catalyzes_pip2_to_ip3():
    """Defines Gaq dependent hydrolysis of PIP2 by PLC into IP3 and DAG.

    This function encodes Gaq binding to PLC which can then catalyzes the
    conversion of PIP2 into IP3 and DAG. There are 3 associated reactions:
        i) Gaq binds PLC:
            Gaq:GTP + PLC <---> Gaq:GTP:PLC
        ii) PIP2 binds to PLC:
            Gaq:GTP:PLC + PIP2 <---> Gaq:GTP:PLC:PIP2
        iii) PLC catalyzes conversion of PIP2 into IP3 and DAG:
            Gaq:GTP:PLC:PIP2 ---> Gaq:GTP:PLC + IP3 + DAG

    Adds 5 parameters.

    Parameters:
        * kf_PLC_bind_Gaq - forward rate constant for Gaq and PLC binding.
                Gaq:GTP + PLC ---> Gaq:GTP:PLC, kf
        * kr_TAT_PAR2_A_bind_Gaq - reverse rate constant for PAR2 binding to
            Gaq in the heterotrimer.
                Gaq:GTP:PLC ---> Gaq:GTP + PLC, kr
        * kf_PLC_bind_PIP2 - forward rate constant for PIP2 binding to PLC.
                Gaq:GTP:PLC + PIP2 ---> Gaq:GTP:PLC:PIP2, kf
        * kr_PLC_bind_PIP2 - rate constant for binding of GDP to Gaq after PAR2
            binding.
                Gaq:GTP:PLC:PIP2 ---> Gaq:GTP:PLC + PIP2, kr
        * kcat_PIP2_to_IP3 - catalytic rate constant for hydrolysis of PIP2 into
            IP3 and DAG by PLC.
                Gaq:GTP:PLC:PIP2 ---> Gaq:GTP:PLC + IP3 + DAG, kcat

    """

    alias_model_components()
    # PLC binding Gaq
    Parameter("kf_PLC_bind_Gaq", defaults.KF_BIND / Vcyto.value)
    Parameter("kr_PLC_bind_Gaq", defaults.KR_BIND)
    # Conversion of PIP2 to IP3
    Parameter("kf_PLC_bind_PIP2", defaults.KF_BIND / Vcyto.value)
    Parameter("kr_PLC_bind_PIP2", defaults.KR_BIND)
    # Nominal values for Ca bound PLC and Gaq_GTP from Flaherty et al.
    # Plos Comp Biol 2008 https://doi.org/10.1371/journal.pcbi.1000185
    #   PLCb4 22.85 1/s
    #   PLCb3 27.89 1/s
    # Note that Flaherty et al. assumes calcium must be bound to PLC to
    # hydrolyze PIP2.
    # We will assume that PLC can do the hydrolysis without Ca2+ but perhaps
    # less efficiently. We can separately encode a function that adds Ca2+
    # binding to PLC and that enhances the catalysis rate which is similar to
    # how Keizer and De Young (Biophys J. 61:549-660 1992) implement
    # calcium feedback on IP3 production.
    Parameter("kcat_PIP2_to_IP3", (22.85 + 27.89) / 4)
    alias_model_components()
    # PLC activation by binding Gaq:
    #    Gaq_A + PLC <---> Gaq_A:PLC
    #   Reusing the Gbg binding slot for PLC
    # Alias the complex Gaq:GTP
    Gaq_gtp = Gaq(bpar=None, bgdp=3, bgbg=None) ** CELL_MEMB % GTP(b=3) ** CELL_MEMB
    bind_complex(
        Gaq_gtp,
        "bgbg",
        PLC(bgaq=None, bpip2=None, bca=None) ** CELL_MEMB,
        "bgaq",
        [kf_PLC_bind_Gaq, kr_PLC_bind_Gaq],
    )
    # Conversion of PIP2 to IP3
    #    Gaq_A:PLC + PIP2 <---> Gaq_A:PLC:PIP2 ---> Gaq_A:PLC + IP3
    Gaq_gtp_PLC = (
        Gaq(bpar=None, bgdp=3, bgbg=1) ** CELL_MEMB
        % GTP(b=3) ** CELL_MEMB
        % PLC(bgaq=1, bpip2=None, bca=None) ** CELL_MEMB
    )
    Gaq_gtp_PLC_pip = (
        Gaq(bpar=None, bgdp=3, bgbg=1) ** CELL_MEMB
        % GTP(b=3) ** CELL_MEMB
        % PLC(bgaq=1, bpip2=2, bca=None) ** CELL_MEMB
        % PIP2(b=2) ** CELL_MEMB
    )
    Rule(
        "plc_binds_pip2",
        Gaq_gtp_PLC + PIP2(b=None) ** CELL_MEMB | Gaq_gtp_PLC_pip,
        kf_PLC_bind_PIP2,
        kr_PLC_bind_PIP2,
    )
    Rule(
        "plc_converts_pip2_to_ip3_and_dag",
        Gaq_gtp_PLC_pip
        >> Gaq_gtp_PLC + IP3(b=None) ** CYTOSOL + DAG(b=None) ** CELL_MEMB,
        kcat_PIP2_to_IP3,
    )

    return


def ip3_binds_ip3r():
    alias_model_components()
    # Binding of IP3 to IP3R
    # Note nominal values from Flaherty et al.
    # Plos Comp Biol 2008 https://doi.org/10.1371/journal.pcbi.1000185
    # for IP3 binidng to IP3R are:
    #   kf: 177.47 1 / (uM s)
    #   kr: 2.2 1/s
    R_o = 2e-7  # cm
    # IP3 diffuses in mammalian tissue at <= 10 micrometer^2/s
    # as per https://dx.doi.org/10.1126%2Fscisignal.aag1625
    D_ip3 = 10e-8  # cm^2/s
    # Assume the IP3 binding rate is diffusion-controlled by IP3 diffusion
    K_IP3_BIND = 4 * np.pi * D_ip3 * R_o * (1e-3) / (Vcyto.value * 1e-12)
    Parameter("kf_IP3_bind_IP3R", K_IP3_BIND)  # Diffusion controlled.
    Parameter("kr_IP3_bind_IP3R", 2.2)  # Nomimal value from Flaherty et al.
    alias_model_components()
    # Binding of IP3 to IP3R
    # Assume subunits are bound sequentially and that there is no cooperativity
    # between subunit binding.
    #   IP3R + IP3 <---> IP3R:IP3, subunit 1
    bind(
        IP3R(b1=None, b2=None, b3=None, b4=None, bcaer=None, bcacyt=WILD) ** ER_MEMB,
        "b1",
        IP3(b=None) ** CYTOSOL,
        "b",
        [kf_IP3_bind_IP3R, kr_IP3_bind_IP3R],
    )
    bind(
        IP3R(b1=ANY, b2=None, b3=None, b4=None, bcaer=None, bcacyt=WILD) ** ER_MEMB,
        "b2",
        IP3(b=None) ** CYTOSOL,
        "b",
        [kf_IP3_bind_IP3R, kr_IP3_bind_IP3R],
    )
    bind(
        IP3R(b1=ANY, b2=ANY, b3=None, b4=None, bcaer=None, bcacyt=WILD) ** ER_MEMB,
        "b3",
        IP3(b=None) ** CYTOSOL,
        "b",
        [kf_IP3_bind_IP3R, kr_IP3_bind_IP3R],
    )
    bind(
        IP3R(b1=ANY, b2=ANY, b3=ANY, b4=None, bcaer=None, bcacyt=WILD) ** ER_MEMB,
        "b4",
        IP3(b=None) ** CYTOSOL,
        "b",
        [kf_IP3_bind_IP3R, kr_IP3_bind_IP3R],
    )
    # #   IP3R + IP3 <---> IP3R:IP3, subunit 2
    # Rule(
    #     "bind_IP3_IPR3_sub2",
    #     IP3R(b1=1, b2=None, b3=None, b4=None, bcaer=None, bcacyt=None) ** ER_MEMB
    #     % IP3(b=1) ** CYTOSOL
    #     + IP3(b=None) ** CYTOSOL
    #     | IP3R(b1=1, b2=2, b3=None, b4=None, bcaer=None, bcacyt=None) ** ER_MEMB
    #     % IP3(b=1) ** CYTOSOL
    #     % IP3(b=2) ** CYTOSOL,
    #     kf_IP3_bind_IP3R,
    #     kr_IP3_bind_IP3R,
    # )
    #   IP3R + IP3 <---> IP3R:IP3, subunit 3
    # Rule(
    #     "bind_IP3_IPR3_sub3",
    #     IP3R(b1=1, b2=2, b3=None, b4=None, bcaer=None, bcacyt=None) ** ER_MEMB
    #     % IP3(b=1) ** CYTOSOL
    #     % IP3(b=2) ** CYTOSOL
    #     + IP3(b=None) ** CYTOSOL
    #     | IP3R(b1=1, b2=2, b3=3, b4=None, bcaer=None, bcacyt=None) ** ER_MEMB
    #     % IP3(b=1) ** CYTOSOL
    #     % IP3(b=2) ** CYTOSOL
    #     % IP3(b=3) ** CYTOSOL,
    #     kf_IP3_bind_IP3R,
    #     kr_IP3_bind_IP3R,
    # )
    # #   IP3R + IP3 <---> IP3R:IP3, subunit 4
    # Rule(
    #     "bind_IP3_IPR3_sub4",
    #     IP3R(b1=1, b2=2, b3=3, b4=None, bcaer=None, bcacyt=None) ** ER_MEMB
    #     % IP3(b=1) ** CYTOSOL
    #     % IP3(b=2) ** CYTOSOL
    #     % IP3(b=3) ** CYTOSOL
    #     + IP3(b=None) ** CYTOSOL
    #     | IP3R(b1=1, b2=2, b3=3, b4=4, bcaer=None, bcacyt=None) ** ER_MEMB
    #     % IP3(b=1) ** CYTOSOL
    #     % IP3(b=2) ** CYTOSOL
    #     % IP3(b=3) ** CYTOSOL
    #     % IP3(b=4) ** CYTOSOL,
    #     kf_IP3_bind_IP3R,
    #     kr_IP3_bind_IP3R,
    # )
    return


def ip3r_transports_er_calcium_to_cytosol():
    alias_model_components()
    # Assume IP3R is only open when all 4 subunits are bound by IP3 and that
    # it only flows ER-->Cytosol
    # Transport of Ca2+
    #  ER -> cytosol:
    # As a first estimate for the forward binding reactions of Ca2+ we will assume
    # that it is diffusion-controlled following Smoluchowski eqn.:
    #     kf = 4*pi*D*R_o
    # We will assume R_o is 2 nm and that D = D_Ca2+.
    # The diffusion coefficient of Ca2+ is 5.3 x 10^-6 cm^2/s; https://doi.org/10.1016/0143-4160(87)90027-3
    D_Ca = 5.3e-6  # cm^2/s
    R_o = 2e-7  # cm
    # (1e-3) term is for unit conversion from cm^3/(s*molecule) to 1/s*molecule*L
    # mL/s*molecule -> 10^-3 L/(s*molecule) and dividing by V converts to 1/(s*molecule/cell)
    K_CA_BIND = 4 * np.pi * D_Ca * R_o * (1e-3) / (Vcyto.value * 1e-12)
    Parameter("kf_erCa_bind_IP3R", K_CA_BIND)
    Parameter("kr_erCa_bind_IP3R", defaults.KR_BIND)
    # Effective IP3R channel permeability as per Lemon et al. 2003 https://doi.org/10.1016/S0022-5193(03)00079-1
    # is 575 1/s
    Parameter("k_tranport_erCa", 525)
    alias_model_components()
    # Transport of Ca2+ by activated IP3R
    #  ER -> cytosol:
    #    IP3R:IP3_4 + Ca_E <---> Ca_E:IP3R:IP3_4 ---> Ca_C + IP3R:IP3_4
    Rule(
        "bind_Ca_IPR3_er",
        IP3R(b1=ANY, b2=ANY, b3=ANY, b4=ANY, bcaer=None, bcacyt=WILD) ** ER_MEMB
        + Ca(b=None) ** ER_LUMEN
        | IP3R(b1=ANY, b2=ANY, b3=ANY, b4=ANY, bcaer=5, bcacyt=WILD) ** ER_MEMB
        % Ca(b=5) ** ER_LUMEN,
        kf_erCa_bind_IP3R,
        kr_erCa_bind_IP3R,
    )
    Rule(
        "transport_Ca_ER_CYTO",
        IP3R(b1=ANY, b2=ANY, b3=ANY, b4=ANY, bcaer=5, bcacyt=WILD) ** ER_MEMB
        % Ca(b=5) ** ER_LUMEN
        >> IP3R(b1=ANY, b2=ANY, b3=ANY, b4=ANY, bcaer=None, bcacyt=WILD) ** ER_MEMB
        + Ca(b=None) ** CYTOSOL,
        k_tranport_erCa,
    )
    return


def ip3_degradation():
    """1st order degradation of IP3.
    IP3 degradation is used to represent depletion or metabolim of IP3 after
    it's catalytic formation from PIP2 by PLC.

    Reactions:
        1. IP3 ---> None
    """
    # Depeletion/metabolism of IP3
    # 1.25 1/s as in Lemon et al. 2003 https://doi.org/10.1016/S0022-5193(03)00079-1
    Parameter("kdeg_ip3", 1.25)
    alias_model_components()
    # Metabolic consumption of IP3
    degrade(IP3(b=None) ** CYTOSOL, kdeg_ip3)
    return


def cytosolic_calcium_buffering():
    alias_model_components()
    Monomer("CalciumBuffer", ["b"])
    # Initial concentratio of calcium buffer from Flaherty et al. 2008
    # Plos Comp Biol 2008 https://doi.org/10.1371/journal.pcbi.1000185
    # is 5.05x10^-2 uM ( or 50.5 nM).
    Parameter("CalciumBuffer_0", 50.5 * units.nM_to_molec_per_pL * Vcyto.value)
    # Rate constants for Ca2+ binding to cytosolic buffer from Flaherty et al.:
    #   kf: 10 1/(uM s)
    #   kr : 7 1/s
    Parameter("kf_Ca_bind_buffer", 10 / (units.microM_to_molec_per_pL * Vcyto.value))
    Parameter("kr_Ca_bind_buffer", 7)
    alias_model_components()

    Initial(CalciumBuffer(b=None) ** CYTOSOL, CalciumBuffer_0)

    bind(
        Ca(b=None) ** CYTOSOL,
        "b",
        CalciumBuffer(b=None) ** CYTOSOL,
        "b",
        [kf_Ca_bind_buffer, kr_Ca_bind_buffer],
    )
    return


def calcium_extrusion_and_influx():
    """Defines reactions control extrusion and influx of Ca2+ from the cytosol and extracellular space.

    This function defines two first order reactions to control extrusion of Ca2+
    from the cytosol to the extracellular space and its influx from the
    extracellular space into the cytosol. This is based on the decay in the Ca2+
    FRET signal seen in Figure S3C of Kang et al. 2019
    (https://doi.org/10.1021/acsnano.9b01993) which we assume is primarily due
    to first-order excretion/extrusion of excess Ca2+ from the released ER store
    into the extracellular space by cell membrane ion channels. We assume that
    at steady-state in the abscence of agonist the extrusion and influx are
    balanced to maintain the resting cyctosolic concentration of free Ca2+.

    Reactions:
        1. Ca2+_CYTO ---> Ca2+_EXTRA
        2. Ca2+_EXTRA ---> Ca2+_CYTO

    Adds 2 parameters.

    Parameters:
        k_Ca_cyt_to_extra - 1st-order rate constant for calcium extrusion from
            the cytosol to the extracellular space.
        k_Ca_extra_to_cyt - 1st-order rate constant for calcium influx from the
            extracellular space to the cytosol.
    """
    # Cytosolic Ca2+ regulation

    # cytosol to extracellular space
    # From previous model fittings we get around 4 1/s.
    Parameter("k_Ca_cyt_to_extra", 4)  # 1/s
    alias_model_components()
    # extracellular space to cytosol
    # As a first estimate assume that the initial concentration of cytosolic
    # Ca2+ is its equilibrium value before any agonist is added and the rate out
    # of the cytosol equals the rate in. However, note that this doesn't take
    # into account any other reactions that can change the cytosolic calcium
    # concentration such as equlibration of Ca2+ buffer binding.
    Parameter(
        "k_Ca_extra_to_cyt", k_Ca_cyt_to_extra.value * Ca_C_0.value / Ca_extra_0.value
    )
    alias_model_components()

    # cytosol to extracellular space
    Rule(
        "Ca_cyt_to_extra",
        Ca(b=None) ** CYTOSOL >> Ca(b=None) ** EXTRACELLULAR,
        k_Ca_cyt_to_extra,
    )
    # extracellular space to cytosol.
    Rule(
        "Ca_extra_to_cyt",
        Ca(b=None) ** EXTRACELLULAR >> Ca(b=None) ** CYTOSOL,
        k_Ca_extra_to_cyt,
    )
    return


def calcium_cytosol_er_flux():
    """Defines reactions control flux of Ca2+ from the cytosol and ER lumen.

    This function defines two first order reactions to control transfer of Ca2+
    from the cytosol (back) to the ER lumen and its transfer from the ER lumen
    space into the cytosol. This function allows reloading of ER Ca2+ store
    after release. We assume that at steady-state in the abscence of agonist the
    extrusion and influx are balanced to maintain the resting cyctosolic
    concentration of free Ca2+.

    Reactions:
        1. Ca2+_CYTO ---> Ca2+V_ERL
        2. Ca2+_ERL ---> Ca2+_CYTO

    Adds 2 parameters.

    Parameters:
        k_Ca_cyt_to_er - 1st-order rate constant for calcium transfer from
            the cytosol to the ER lumen.
        k_Ca_er_to_cyt - 1st-order rate constant for calcium transfer from the
            ER lumen to the cytosol.
    """
    # Cytosolic Ca2+ regulation

    # cytosol to ER lumen
    Parameter("k_Ca_cyt_to_er", 4)  # 1/s
    alias_model_components()
    # ER lumen to cytosol
    # As a first estimate assume that the initial concentration of cytosolic
    # Ca2+ is its equilibrium value before any agonist is added and the rate out
    # of the cytosol equals the rate in. However, note that this doesn't take
    # into account any other reactions that can change the cytosolic calcium
    # concentration such as equlibration of Ca2+ buffer binding.
    Parameter("k_Ca_er_to_cyt", k_Ca_cyt_to_er.value * Ca_C_0.value / Ca_E_0.value)
    alias_model_components()

    # cytosol to ER lumen
    Rule(
        "Ca_cyt_to_er",
        Ca(b=None) ** CYTOSOL >> Ca(b=None) ** ER_LUMEN,
        k_Ca_cyt_to_er,
    )
    # ER lumen to cytosol.
    Rule(
        "Ca_er_to_cyt",
        Ca(b=None) ** ER_LUMEN >> Ca(b=None) ** CYTOSOL,
        k_Ca_er_to_cyt,
    )
    return


def calcium_extrusion_and_influx_mk():
    """Defines reactions control extrusion and influx of Ca2+ from the cytosol and extracellular space.

    This function defines two first order reactions with rate constant
    expressions that yield  Michaelis-Menten kinetics to control extrusion of
    Ca2+ from the cytosol to the extracellular space and its influx from the
    extracellular space into the cytosol. This is based on the decay in the Ca2+
    FRET signal seen in Figure S3C of Kang et al. 2019
    (https://doi.org/10.1021/acsnano.9b01993) which we assume is largely due to
    first-order excretion/extrusion of excess Ca2+ from the released ER store
    into the extracellular space by cell membrane ion channels. We assume that
    at steady-state in the abscence of agonist the extrusion and influx should
    be balanced to maintain the resting cyctosolic concentration of free Ca2+, but
    this must enforced when training model.

    Reactions:
        1. Ca2+_CYTO ---> Ca2+_EXTRA
        2. Ca2+_EXTRA ---> Ca2+_CYTO

    Adds 4 parameters and two expressions.

    Parameters:
        Vmax_Ca_cyt_to_extra
        Km_Ca_cyt_to_extra
        Vmax_Ca_extra_to_cyt
        Km_Ca_extra_to_cyt

    Expressions:
        k_Ca_cyt_to_extra - 1st-order rate constant for calcium extrusion from
            the cytosol to the extracellular space.
        k_Ca_extra_to_cyt - 1st-order rate constant for calcium influx from the
            extracellular space to the cytosol.
    """

    alias_model_components()
    # cytosol to extracellular space
    # Set a nominal max rate of 20 uM/s.
    Parameter("Vmax_Ca_cyt_to_extra", 10000 * units.nM_to_molec_per_pL * Vcyto.value)
    # Nominal value of 1 uM.
    Parameter("Km_Ca_cyt_to_extra", 0.15 * units.microM_to_molec_per_pL * Vcyto.value)
    # Parameter("n_Ca_cyt_to_extra", 2)
    # extracellular space to cytosol
    # Nominal value of 10 mM.
    Parameter("Km_Ca_extra_to_cyt", 10000 * units.microM_to_molec_per_pL * Ver.value)
    alias_model_components()
    # Nominal value of 10 uM/s.
    Parameter(
        "Vmax_Ca_extra_to_cyt",
        Vmax_Ca_cyt_to_extra.value
        * Ca_C_0.value
        * (Km_Ca_extra_to_cyt.value + Ca_extra_0.value)
        / ((Km_Ca_cyt_to_extra.value + Ca_C_0.value) * Ca_extra_0.value),
    )  # 10 * units.microM_to_molec_per_pL * Vcyto.value)
    # Some 'private' observables to monitor the ER and cytoslic calcium for
    # Michaelis-Menten rate.
    Observable("_Ca_extra", Ca(b=None) ** EXTRACELLULAR)
    Observable("_Ca_cyt", Ca(b=None) ** CYTOSOL)
    alias_model_components()
    # First order rate constants based on MK rate.
    Expression(
        "k_Ca_cyt_to_extra",
        (Vmax_Ca_cyt_to_extra) / (_Ca_cyt + Km_Ca_cyt_to_extra),
    )  # 1/s
    Expression(
        "k_Ca_extra_to_cyt", Vmax_Ca_extra_to_cyt / (_Ca_extra + Km_Ca_extra_to_cyt)
    )
    alias_model_components()
    # cytosol to extracellular space
    Rule(
        "Ca_cyt_to_extra",
        Ca(b=None) ** CYTOSOL >> Ca(b=None) ** EXTRACELLULAR,
        k_Ca_cyt_to_extra,
    )
    # extracellular space to cytosol.
    Rule(
        "Ca_extra_to_cyt",
        Ca(b=None) ** EXTRACELLULAR >> Ca(b=None) ** CYTOSOL,
        k_Ca_extra_to_cyt,
    )
    return


def calcium_cytosol_er_flux_mk():
    """Defines reactions control flux of Ca2+ from the cytosol and ER lumen.

    This function defines two first order reactions to control transfer of Ca2+
    from the cytosol (back) to the ER lumen and its transfer from the ER lumen
    space into the cytosol. This function allows reloading of ER Ca2+ store
    after release. We assume that at steady-state in the abscence of agonist the
    extrusion and influx are balanced to maintain the resting cyctosolic
    concentration of free Ca2+.

    Reactions:
        1. Ca2+_CYTO ---> Ca2+V_ERL
        2. Ca2+_ERL ---> Ca2+_CYTO

    Adds 2 parameters.

    Parameters:
        k_Ca_cyt_to_er - 1st-order rate constant for calcium transfer from
            the cytosol to the ER lumen.
        k_Ca_er_to_cyt - 1st-order rate constant for calcium transfer from the
            ER lumen to the cytosol.
    """
    # Cytosolic Ca2+ regulation

    # cytosol to ER lumen
    alias_model_components()
    Parameter("Vmax_Ca_cyt_to_er", 20 * units.microM_to_molec_per_pL * Vcyto.value)
    Parameter("Km_Ca_cyt_to_er", 0.65 * units.microM_to_molec_per_pL * Vcyto.value)
    # Parameter("n_Ca_cyt_to_er", 2)
    Parameter("Km_Ca_er_to_cyt", 1000 * units.microM_to_molec_per_pL * Ver.value)
    alias_model_components()
    Parameter(
        "Vmax_Ca_er_to_cyt",
        Vmax_Ca_cyt_to_er.value
        * Ca_C_0.value
        * (Km_Ca_er_to_cyt.value + Ca_E_0.value)
        / ((Km_Ca_cyt_to_er.value + Ca_C_0.value) * Ca_E_0.value),
    )
    Observable("__Ca_er", Ca(b=None) ** ER_LUMEN)
    Observable("__Ca_cyt", Ca(b=None) ** CYTOSOL)
    alias_model_components()
    Expression(
        "k_Ca_cyt_to_er", Vmax_Ca_cyt_to_er / (__Ca_cyt + Km_Ca_cyt_to_er)
    )  # 1/s
    Expression("k_Ca_er_to_cyt", Vmax_Ca_er_to_cyt / (__Ca_er + Km_Ca_er_to_cyt))
    alias_model_components()
    # ER lumen to cytosol
    # As a first estimate assume that the initial concentration of cytosolic
    # Ca2+ is its equilibrium value before any agonist is added and the rate out
    # of the cytosol equals the rate in. However, note that this doesn't take
    # into account any other reactions that can change the cytosolic calcium
    # concentration such as equlibration of Ca2+ buffer binding.

    # cytosol to ER lumen
    Rule(
        "Ca_cyt_to_er",
        Ca(b=None) ** CYTOSOL >> Ca(b=None) ** ER_LUMEN,
        k_Ca_cyt_to_er,
    )
    # ER lumen to cytosol.
    Rule(
        "Ca_er_to_cyt",
        Ca(b=None) ** ER_LUMEN >> Ca(b=None) ** CYTOSOL,
        k_Ca_er_to_cyt,
    )
    return


def calcium_extrusion_and_influx_single():
    """Defines reactions control extrusion and influx of Ca2+ from the cytosol and extracellular space.

    This function defines two first order reactions to control extrusion of Ca2+
    from the cytosol to the extracellular space and its influx from the
    extracellular space into the cytosol. This is based on the decay in the Ca2+
    FRET signal seen in Figure S3C of Kang et al. 2019
    (https://doi.org/10.1021/acsnano.9b01993) which we assume is primarily due
    to first-order excretion/extrusion of excess Ca2+ from the released ER store
    into the extracellular space by cell membrane ion channels. We assume that
    at steady-state in the abscence of agonist the extrusion and influx are
    balanced to maintain the resting cyctosolic concentration of free Ca2+.

    Reactions:
        1. Ca2+_CYTO ---> Ca2+_EXTRA
        2. Ca2+_EXTRA ---> Ca2+_CYTO

    Adds 2 parameters.

    Parameters:
        k_Ca_cyt_to_extra - 1st-order rate constant for calcium extrusion from
            the cytosol to the extracellular space.
        k_Ca_extra_to_cyt - 1st-order rate constant for calcium influx from the
            extracellular space to the cytosol.
    """
    # Cytosolic Ca2+ regulation

    # cytosol to extracellular space
    # From previous model fittings we get around 4 1/s.
    Parameter("k_Ca_cyt_to_extra", 5e-2)  # 1/s
    alias_model_components()
    # # free Ca2+ in the extracellular space
    Observable("_cytCa", Ca(b=None) ** CYTOSOL)
    alias_model_components()
    Expression("_k_Ca_cyt_to_extra", k_Ca_cyt_to_extra * (_cytCa - Ca_C_0) / _cytCa)
    alias_model_components()
    # cytosol to extracellular space
    Rule(
        "Ca_cyt_to_extra",
        Ca(b=None) ** CYTOSOL >> Ca(b=None) ** EXTRACELLULAR,
        _k_Ca_cyt_to_extra,
    )

    return

def calcium_extrusion_and_influx_single_o2():
    """Defines reactions control extrusion and influx of Ca2+ from the cytosol and extracellular space.

    This function defines two first order reactions to control extrusion of Ca2+
    from the cytosol to the extracellular space and its influx from the
    extracellular space into the cytosol. This is based on the decay in the Ca2+
    FRET signal seen in Figure S3C of Kang et al. 2019
    (https://doi.org/10.1021/acsnano.9b01993) which we assume is primarily due
    to first-order excretion/extrusion of excess Ca2+ from the released ER store
    into the extracellular space by cell membrane ion channels. We assume that
    at steady-state in the abscence of agonist the extrusion and influx are
    balanced to maintain the resting cyctosolic concentration of free Ca2+.

    Reactions:
        1. Ca2+_CYTO ---> Ca2+_EXTRA
        2. Ca2+_EXTRA ---> Ca2+_CYTO

    Adds 2 parameters.

    Parameters:
        k_Ca_cyt_to_extra - 1st-order rate constant for calcium extrusion from
            the cytosol to the extracellular space.
        k_Ca_extra_to_cyt - 1st-order rate constant for calcium influx from the
            extracellular space to the cytosol.
    """
    # Cytosolic Ca2+ regulation

    # cytosol to extracellular space
    # From previous model fittings we get around 4 1/s.
    Parameter("k_Ca_cyt_to_extra", 5e-2)  # 1/s
    alias_model_components()
    # # free Ca2+ in the extracellular space
    Observable("_cytCa", Ca(b=None) ** CYTOSOL)
    alias_model_components()
    Expression("_k_Ca_cyt_to_extra", k_Ca_cyt_to_extra * (_cytCa - Ca_C_0)**2 / _cytCa)
    alias_model_components()
    # cytosol to extracellular space
    Rule(
        "Ca_cyt_to_extra",
        Ca(b=None) ** CYTOSOL >> Ca(b=None) ** EXTRACELLULAR,
        _k_Ca_cyt_to_extra,
    )

    return

def calcium_extrusion_and_influx_single_mk():
    """Defines reactions control extrusion and influx of Ca2+ from the cytosol and extracellular space.

    This function defines two first order reactions with rate constant
    expressions that yield  Michaelis-Menten kinetics to control extrusion of
    Ca2+ from the cytosol to the extracellular space and its influx from the
    extracellular space into the cytosol. This is based on the decay in the Ca2+
    FRET signal seen in Figure S3C of Kang et al. 2019
    (https://doi.org/10.1021/acsnano.9b01993) which we assume is largely due to
    first-order excretion/extrusion of excess Ca2+ from the released ER store
    into the extracellular space by cell membrane ion channels. We assume that
    at steady-state in the abscence of agonist the extrusion and influx should
    be balanced to maintain the resting cyctosolic concentration of free Ca2+, but
    this must enforced when training model.

    Reactions:
        1. Ca2+_CYTO ---> Ca2+_EXTRA
        2. Ca2+_EXTRA ---> Ca2+_CYTO

    Adds 4 parameters and two expressions.

    Parameters:
        Vmax_Ca_cyt_to_extra
        Km_Ca_cyt_to_extra
        Vmax_Ca_extra_to_cyt
        Km_Ca_extra_to_cyt

    Expressions:
        k_Ca_cyt_to_extra - 1st-order rate constant for calcium extrusion from
            the cytosol to the extracellular space.

    """

    alias_model_components()
    # cytosol to extracellular space
    # Set a nominal max rate of 20 uM/s.
    Parameter("Vmax_Ca_cyt_to_extra", 20 * units.nM_to_molec_per_pL * Vcyto.value)
    # Nominal value of 1 uM.
    Parameter("Km_Ca_cyt_to_extra", 1. * units.microM_to_molec_per_pL * Vcyto.value)
    Parameter("n_Ca_cyt_to_extra", 2)
    alias_model_components()
    # Some 'private' observables to monitor the ER and cytoslic calcium for
    # Michaelis-Menten rate.
    Observable("_Ca_cyt", Ca(b=None) ** CYTOSOL)
    alias_model_components()
    # First order rate constants based on MK rate.
    Expression(
        "k_Ca_cyt_to_extra",
        (Vmax_Ca_cyt_to_extra * (_Ca_cyt - Ca_C_0)**n_Ca_cyt_to_extra) / (_Ca_cyt *(_Ca_cyt - Ca_C_0)**n_Ca_cyt_to_extra + Km_Ca_cyt_to_extra**n_Ca_cyt_to_extra),
    )  # 1/s
    alias_model_components()
    # cytosol to extracellular space
    Rule(
        "Ca_cyt_to_extra",
        Ca(b=None) ** CYTOSOL >> Ca(b=None) ** EXTRACELLULAR,
        k_Ca_cyt_to_extra,
    )
    return

def calcium_cytosol_er_flux_single():
    """Defines reactions control flux of Ca2+ from the cytosol and ER lumen.

    This function defines two first order reactions to control transfer of Ca2+
    from the cytosol (back) to the ER lumen and its transfer from the ER lumen
    space into the cytosol. This function allows reloading of ER Ca2+ store
    after release. We assume that at steady-state in the abscence of agonist the
    extrusion and influx are balanced to maintain the resting cyctosolic
    concentration of free Ca2+.

    Reactions:
        1. Ca2+_CYTO ---> Ca2+V_ERL
        2. Ca2+_ERL ---> Ca2+_CYTO

    Adds 2 parameters.

    Parameters:
        k_Ca_cyt_to_er - 1st-order rate constant for calcium transfer from
            the cytosol to the ER lumen.
        k_Ca_er_to_cyt - 1st-order rate constant for calcium transfer from the
            ER lumen to the cytosol.
    """
    # Cytosolic Ca2+ regulation

    # cytosol to ER lumen
    # Parameter("k_Ca_cyt_to_er", 2)  # 1/s
    # alias_model_components()
    # ER lumen to cytosol
    # As a first estimate assume that the initial concentration of cytosolic
    # Ca2+ is its equilibrium value before any agonist is added and the rate out
    # of the cytosol equals the rate in. However, note that this doesn't take
    # into account any other reactions that can change the cytosolic calcium
    # concentration such as equlibration of Ca2+ buffer binding.
    Parameter("k_Ca_er_to_cyt", 1e-3)
    # alias_model_components()
    # The Ca2+ in the ER Lumen
    Observable("_erCa", Ca(b=None) ** ER_LUMEN)
    alias_model_components()
    Expression("_k_Ca_er_to_cyt", k_Ca_er_to_cyt * (_erCa - Ca_E_0) / _erCa)
    alias_model_components()

    # cytosol to ER lumen
    # Rule(
    #     "Ca_cyt_to_er",
    #     Ca(b=None) ** CYTOSOL >> Ca(b=None) ** ER_LUMEN,
    #     k_Ca_cyt_to_er,
    # )
    # ER lumen to cytosol.
    Rule(
        "Ca_er_to_cyt",
        Ca(b=None) ** ER_LUMEN >> Ca(b=None) ** CYTOSOL,
        _k_Ca_er_to_cyt,
    )
    return



def regulation_of_cytosolic_calcium_concentration():
    """Reactions to regulate the concentration of free cytosolic calcium.

    This function combines cytosolic calcium buffering with extrusion and
    and influx, thus defining some key reactions in the regulation of the free
    cytosolic calcium concentration.

    Calls:
        * cytosolic_calcium_buffering
        * calcium_extrusion_and_influx
        * calcium_cytosol_er_flux
    """

    cytosolic_calcium_buffering()
    calcium_extrusion_and_influx_single()
    # calcium_cytosol_er_flux()
    return


def calcium_binds_plc_and_enhances_pip2_hydrolysis():
    """Defines positive feedback for cytosolic Ca2+ to enhance IP3 production.

    This function encodes cytosolic calcium binding to PLC to provide positive
    feedback by enhancing the production of IP3. It is assumed that the kcat
    for hydrolysis of PIP2 is higher when calcium is bound to PLC.

    See Keizer and De Young (Biophys J. 61:549-660 1992
    https://doi.org/10.1016/S0006-3495(92)81870-2) and Flaherty et al. (PloS
    Comput Biol 4(9):e1000185 2008
    https://doi.org/10.1371/journal.pcbi.1000185).

    Encodes the following reactions:
        1. Ca_C + PLC <---> PLC:Ca_C

    Adds 3 parameters and 1 expression:

    Parameters:
        * kf_PLC_bind_Ca - forward rate constant for Ca2+ binding to PLC.
        * kr_PLC_bind_Ca - reverse rate constant for
        * kcat_PIP2_to_IP3_Ca

    Expressions:
        * ef_ca_pip2_to_ip3
    """

    alias_model_components()
    # PLC binding cytosolic calcium and enhancement of PIP2 conversion

    # Nominal value from Flaherty et al. 2008
    # Plos Comp Biol 2008 https://doi.org/10.1371/journal.pcbi.1000185
    # Forward binding rate constant:
    #   PLCb4 20 1/(uM s)
    #   PLCb3 20 1/(uM s)
    Parameter("kf_PLC_bind_Ca", (20 / units.microM_to_molec_per_pL) / Vcyto.value)
    # Reverse binding rate constant
    #   PLCb4 8 1/s
    #   PLCb3 8 1/s
    Parameter("kr_PLC_bind_Ca", 8)
    # catalytic rate:
    #   PLCb4 22.85 1/s
    #   PLCb3 27.89 1/s
    Parameter("kcat_PIP2_to_IP3_Ca", (22.85 + 27.89) / 2)
    alias_model_components()
    Expression(
        "ef_ca_pip2_to_ip3", kcat_PIP2_to_IP3_Ca / kcat_PIP2_to_IP3
    )  # ef = enhancement factor
    PLC_Ca = PLC(bgaq=None, bpip2=None, bca=ANY) ** CELL_MEMB
    Gaq_gtp_PLC = (
        Gaq(bpar=None, bgdp=3, bgbg=1) ** CELL_MEMB
        % GTP(b=3) ** CELL_MEMB
        % PLC(bgaq=1, bpip2=None, bca=None) ** CELL_MEMB
    )
    Gaq_gtp_PLC_Ca = (
        Gaq(bpar=None, bgdp=3, bgbg=1) ** CELL_MEMB
        % GTP(b=3) ** CELL_MEMB
        % PLC(bgaq=1, bpip2=None, bca=4) ** CELL_MEMB
        % Ca(b=4) ** CYTOSOL
    )
    # Alias the complex Gaq:GTP
    Gaq_gtp = Gaq(bpar=None, bgdp=3, bgbg=None) ** CELL_MEMB % GTP(b=3) ** CELL_MEMB
    # PLC binds cytosolic Ca2+
    bind(
        PLC(bgaq=WILD, bpip2=None, bca=None) ** CELL_MEMB,
        "bca",
        Ca(b=None) ** CYTOSOL,
        "b",
        [kf_PLC_bind_Ca, kr_PLC_bind_Ca],
    )
    # Gaq:GTP binds PLC:Ca
    bind_complex(
        Gaq_gtp,
        "bgbg",
        PLC(bgaq=None, bpip2=None, bca=ANY) ** CELL_MEMB,
        "bgaq",
        [kf_PLC_bind_Gaq, kr_PLC_bind_Gaq],
    )
    # catalyze(PLC(bgaq=ANY, bca=ANY)**CELL_MEMB, 'bpip2',
    #          PIP2()**CELL_MEMB, 'b',
    #          (IP3(b=None)**CYTOSOL + DAG()**CYTOSOL),
    #          [kf_PLC_bind_PIP2, kr_PLC_bind_PIP2, kcat_PIP2_to_IP3_Ca])
    Rule(
        "plc_binds_pip2_ca",
        PLC(bgaq=ANY, bpip2=None, bca=ANY) ** CELL_MEMB + PIP2(b=None) ** CELL_MEMB
        | PLC(bgaq=ANY, bpip2=2, bca=ANY) ** CELL_MEMB % PIP2(b=2) ** CELL_MEMB,
        kf_PLC_bind_PIP2,
        kr_PLC_bind_PIP2,
    )
    Rule(
        "plc_converts_pip2_to_ip3_and_dag_ca",
        PLC(bgaq=ANY, bpip2=2, bca=ANY) ** CELL_MEMB % PIP2(b=2) ** CELL_MEMB
        >> PLC(bgaq=ANY, bpip2=None, bca=ANY) ** CELL_MEMB
        + IP3(b=None) ** CYTOSOL
        + DAG(b=None) ** CELL_MEMB,
        kcat_PIP2_to_IP3_Ca,
    )

    return


def cytosolic_calcium_inhibits_ip3r():
    """Defines negative feedback from cytosolic calcium that inhibits IP3R.

    The negative calcium feedback on IP3R is modeled by having calcium
    competitively inhibit the IP3 binding sites at IP3R. This mechanistic choice
    is based in part on the negative feedback mechanism implemented by Keizer
    and De Young (Biophys J. 61:549-660 1992
    https://doi.org/10.1016/S0006-3495(92)81870-2) as well as that proposed by
    Taylor and Tovey (Cold Spring Harb Perspect Biol 2010;2:a004010
    https://dx.doi.org/10.1101%2Fcshperspect.a004010, see Figure 1A). Note that
    Keizer and De Young implement the feedback as an allosteric modulation of
    the IP3-IP3R dissociation constant when Ca2+ is bound to a sub-unit.
    Competitive inhibition should have a similar effect of shifting up the EC50
    for IP3 binding while approximating the inhibited receptor that has calcium
    bound in the inhibitory binding site as described by Taylor and Tovey.

    Calcium can bind to each of the binding spots corresonding to the IP3R
    subunits at the IP3 binding sites (b1-b4).

    Adds 2 parameters and 1 expression.

    Parameters:
        * kf_cytCa_bind_IP3R - forward rate constant for calcium binding to
            IP3R.
        * Kd_cytCa_bind_IP3R - dissociation constant for caclium binding to
            IP3R.

    Expressions:
        * kr_cytCa_bind_IP3R - reverse rate constant for calcium binding to
            IP3R.

    """
    alias_model_components()
    # As a first estimate for the forward binding reactions of Ca2+ we will assume
    # that it is diffusion-controlled following Smoluchowski eqn.:
    #     kf = 4*pi*D*R_o
    # We will assume R_o is 2 nm and that D = D_Ca2+.
    # The diffusion coefficient of Ca2+ is 5.3 x 10^-6 cm^2/s; https://doi.org/10.1016/0143-4160(87)90027-3
    D_Ca = 5.3e-6  # cm^2/s
    R_o = 2e-7  # cm
    # (1e-3) term is for unit conversion from cm^3/(s*molecule) to 1/s*molecule*L
    # mL/s*molecule -> 10^-3 L/(s*molecule) and dividing by V converts to 1/(s*molecule/cell)
    K_CA_BIND = 4 * np.pi * D_Ca * R_o * (1e-3) / (Vcyto.value * 1e-12)
    Parameter("kf_cytCa_bind_IP3R", K_CA_BIND)
    Parameter("Kd_cytCa_bind_IP3R", 10.0 * units.microM_to_molec_per_pL * Vcyto.value)
    alias_model_components()
    Expression("kr_cytCa_bind_IP3R", kf_cytCa_bind_IP3R * Kd_cytCa_bind_IP3R / V_C)
    alias_model_components()
    # Binding of Ca to IP3R
    # Assume subunits are bound sequentially and that there is no cooperativity
    # between subunit binding.
    #   IP3R + Ca <---> IP3R:IP3, subunit 1
    bind(
        IP3R(b1=None, b2=None, b3=None, b4=None, bcaer=None, bcacyt=None) ** ER_MEMB,
        "b1",
        Ca(b=None) ** CYTOSOL,
        "b",
        [kf_cytCa_bind_IP3R, kr_cytCa_bind_IP3R],
    )
    bind(
        IP3R(b1=ANY, b2=None, b3=None, b4=None, bcaer=None, bcacyt=None) ** ER_MEMB,
        "b2",
        Ca(b=None) ** CYTOSOL,
        "b",
        [kf_cytCa_bind_IP3R, kr_cytCa_bind_IP3R],
    )
    bind(
        IP3R(b1=ANY, b2=ANY, b3=None, b4=None, bcaer=None, bcacyt=None) ** ER_MEMB,
        "b3",
        Ca(b=None) ** CYTOSOL,
        "b",
        [kf_cytCa_bind_IP3R, kr_cytCa_bind_IP3R],
    )
    bind(
        IP3R(b1=ANY, b2=ANY, b3=ANY, b4=None, bcaer=None, bcacyt=None) ** ER_MEMB,
        "b4",
        Ca(b=None) ** CYTOSOL,
        "b",
        [kf_cytCa_bind_IP3R, kr_cytCa_bind_IP3R],
    )
    return

def cytosolic_calcium_enhances_ip3r_calcium_transports():
    alias_model_components()
    # As a first estimate for the forward binding reactions of Ca2+ we will assume
    # that it is diffusion-controlled following Smoluchowski eqn.:
    #     kf = 4*pi*D*R_o
    # We will assume R_o is 2 nm and that D = D_Ca2+.
    # The diffusion coefficient of Ca2+ is 5.3 x 10^-6 cm^2/s; https://doi.org/10.1016/0143-4160(87)90027-3
    D_Ca = 5.3e-6  # cm^2/s
    R_o = 2e-7  # cm
    # (1e-3) term is for unit conversion from cm^3/(s*molecule) to 1/s*molecule*L
    # mL/s*molecule -> 10^-3 L/(s*molecule) and dividing by V converts to 1/(s*molecule/cell)
    K_CA_BIND = 4 * np.pi * D_Ca * R_o * (1e-3) / (Vcyto.value * 1e-12)
    Parameter("kf_cytCa_bind_IP3R", K_CA_BIND)
    Parameter("Kd_cytCa_bind_IP3R", 10.0 * units.microM_to_molec_per_pL * Vcyto.value)
    alias_model_components()
    Expression("kr_cytCa_bind_IP3R", kf_cytCa_bind_IP3R * Kd_cytCa_bind_IP3R / V_C)
    alias_model_components()
    bind(
        IP3R(b1=WILD, b2=WILD, b3=WILD, b4=WILD, bcaer=None, bcacyt=None) ** ER_MEMB,
        "bcacyt",
        Ca(b=None) ** CYTOSOL,
        "b",
        [kf_cytCa_bind_IP3R, kr_cytCa_bind_IP3R],
    )
    # Effective IP3R channel permeability as per Lemon et al. 2003 https://doi.org/10.1016/S0022-5193(03)00079-1
    # is 575 1/s
    Parameter("k_tranport_erCa_Ca", 575*10)
    alias_model_components()
    Rule(
        "transport_Ca_ER_CYTO_Ca",
        IP3R(b1=ANY, b2=ANY, b3=ANY, b4=ANY, bcaer=5, bcacyt=ANY) ** ER_MEMB
        % Ca(b=5) ** ER_LUMEN
        >> IP3R(b1=ANY, b2=ANY, b3=ANY, b4=ANY, bcaer=None, bcacyt=ANY) ** ER_MEMB
        + Ca(b=None) ** CYTOSOL,
        k_tranport_erCa_Ca,
    )
    return

def cytosolic_calcium_feedback():
    """Defines feedback from cytosolic Ca2+ on IP3 production and IP3R opening.

    This function combines the mechanism for positive feedback on IP3 production
    and negative feedback on IP3R gating.

    Calls:
        * calcium_binds_plc_and_enhances_pip2_hydrolysis
        * cytosolic_calcium_inhibits_ip3r

    """
    # Positive feedback on hydrolysis of PIP2.
    calcium_binds_plc_and_enhances_pip2_hydrolysis()
    # Negative feedback (inhibition) of IP3R gating.
    cytosolic_calcium_inhibits_ip3r()
    return


def cytosolic_calcium_positive_feedback():
    """Defines feedback from cytosolic Ca2+ on IP3 production.

    This function wraps the mechanism for positive feedback on IP3 production
    by Ca2+ binding to PLC.

    Calls:
        * calcium_binds_plc_and_enhances_pip2_hydrolysis
        * cytosolic_calcium_inhibits_ip3r

    """
    # Positive feedback on hydrolysis of PIP2.
    calcium_binds_plc_and_enhances_pip2_hydrolysis()
    #cytosolic_calcium_enhances_ip3r_calcium_transports()
    return


def gaq_activated_calcium_signaling():
    """Combines mechanistic elements to define the full calcium signaling pathway.

    Calls:
        * calcium_signal_monomers
        * calcium_signal_initials
        * plc_binds_gaq_and_catalyzes_pip2_to_ip3
        * ip3_binds_ip3r
        * ip3r_transports_er_calcium_to_cytosol
        * regulation_of_cytosolic_calcium_concentration
        * cytosolic_calcium_feedback
        * ip3_degradation

    """
    calcium_signal_monomers()
    calcium_signal_initials()
    plc_binds_gaq_and_catalyzes_pip2_to_ip3()
    ip3_binds_ip3r()
    ip3r_transports_er_calcium_to_cytosol()
    regulation_of_cytosolic_calcium_concentration()
    cytosolic_calcium_feedback()
    ip3_degradation()
    return


def gaq_activated_calcium_signaling_simplified():
    """Combines mechanistic elements to define the calcium signaling pathway.

    This simplified version of the calcium pathway does not include calcium
    buffering in the cytosol or cytosolic calcium feedback mechanisms.

    Calls:
        * calcium_signal_monomers
        * calcium_signal_initials
        * plc_binds_gaq_and_catalyzes_pip2_to_ip3
        * ip3_binds_ip3r
        * ip3r_transports_er_calcium_to_cytosol
        * calcium_extrusion_and_influx
        * ip3_degradation

    """
    calcium_signal_monomers()
    calcium_signal_initials()
    plc_binds_gaq_and_catalyzes_pip2_to_ip3()
    ip3_binds_ip3r()
    ip3r_transports_er_calcium_to_cytosol()
    cytosolic_calcium_positive_feedback()
    calcium_extrusion_and_influx_single()
    ip3_degradation()
    return


def gaq_activated_calcium_signaling_without_influx_extrusion():
    """Combines mechanistic elements to define the calcium signaling pathway.

    This simplified version of the calcium pathway does not include calcium
    buffering in the cytosol or cytosolic calcium feedback mechanisms.

    Calls:
        * calcium_signal_monomers
        * calcium_signal_initials
        * plc_binds_gaq_and_catalyzes_pip2_to_ip3
        * ip3_binds_ip3r
        * ip3r_transports_er_calcium_to_cytosol
        * cytosolic_calcium_buffering
        * cytosolic_calcium_feedback
        * ip3_degradation

    """
    calcium_signal_monomers()
    calcium_signal_initials()
    plc_binds_gaq_and_catalyzes_pip2_to_ip3()
    ip3_binds_ip3r()
    ip3r_transports_er_calcium_to_cytosol()
    cytosolic_calcium_buffering()
    cytosolic_calcium_feedback()
    ip3_degradation()
    return


def gaq_activated_calcium_signaling_minimal():
    """Combines mechanistic elements to define a minimalistic version of the calcium signaling pathway.

    This simplified version of the calcium pathway does not include calcium
    buffering, cytosolic calcium feedback mechanisms, or pump reactions to
    maintain calcium homeostasis.

    Calls:
        * calcium_signal_monomers
        * calcium_signal_initials
        * plc_binds_gaq_and_catalyzes_pip2_to_ip3
        * ip3_binds_ip3r
        * ip3r_transports_er_calcium_to_cytosol
        * ip3_degradation

    """
    calcium_signal_monomers()
    calcium_signal_initials()
    plc_binds_gaq_and_catalyzes_pip2_to_ip3()
    ip3_binds_ip3r()
    ip3r_transports_er_calcium_to_cytosol()
    ip3_degradation()
    return


def observables():
    """Defines observables for the calcium signaling pathway.

    Defines 7 observables and 3 expressions.

    Observables:
        totPIP2 = total amount of PIP2
        totIP3 = total amount of IP3
        totIP3R = total amount of IP3R
        aIP3R = active IP3R with IP3 bound at each subunit binding site
        erCa = amount of free Ca2+ in the ER lumen.
        cytoCa = amount of free Ca2+ in the cytosol.
        bound_cytoCa = amount of cytosolic Ca2+ bound to something
        buffered_cytoCa = amount of cytosolic Ca2+ bound to the calcium buffer protein.
        PLC_Ca = amount of PLC with Ca2+ bound to the calcium binding site.

    Expressions:
        active_IP3R_ratio = ratio of aIP3R to totIP3R
        Ca_num_to_nM = a conversion factor to convert the amount of cytosolic
            Ca2+ to nM concentrations.
        cytoCa_nM = concentration of free cytosolic Ca2+ in nM.

    """
    alias_model_components()
    Observable("totPIP2", PIP2())
    Observable("totIP3", IP3())

    Observable("totIP3R", IP3R())
    # Active IP3R (i.e., all 4 subunits bound by IP3)
    Observable(
        "aIP3R",
        IP3R(b1=1, b2=2, b3=3, b4=4, bcaer=WILD, bcacyt=WILD) ** ER_MEMB
        % IP3(b=1) ** CYTOSOL
        % IP3(b=2) ** CYTOSOL
        % IP3(b=3) ** CYTOSOL
        % IP3(b=4) ** CYTOSOL,
    )
    # The Ca2+ in the ER Lumen
    Observable("erCa", Ca(b=None) ** ER_LUMEN)
    # free Ca2+ in the Cytosol
    Observable("cytoCa", Ca(b=None) ** CYTOSOL)
    alias_model_components()

    Expression("active_IP3R_ratio", aIP3R / totIP3R)
    Expression("Ca_num_to_nM", 1 / (Vcyto * units.nM_to_molec_per_pL))
    alias_model_components()
    Expression("cytoCa_nM", cytoCa * Ca_num_to_nM)
    Expression("erCa_uM", erCa / (Ver * units.microM_to_molec_per_pL))
    # Bound cytosolic calcium.
    Observable("bound_cytoCa", Ca(b=ANY) ** CYTOSOL)
    # Buffered caclcium.
    # try:
    #     cba = kf_Ca_bind_buffer.value
    #     Observable("buffered_cytoCa", CalciumBuffer(b=ANY) ** CYTOSOL)
    # except:
    #     pass
    # Calcium bound PLC
    Observable("PLC_Ca", PLC(bca=ANY))

    return


def fret_calcium_indicator_tnxxl():
    """Defines parameters and expressions to get the FRET ratio for TN-XXL calcium sensing.

    TN-XXL is a genetically encoded ratiometric FRET sensor for the cytosolic
    level of Ca2+. It is what is used to measure the caclium response  of HEK293
    cells and its changes after Molecular Hyperthermia in Kang et al. 2019
    (https://doi.org/10.1021/acsnano.9b01993).

    The FRET ratio as a function of Ca2+ concentration is modeled using an
    empirically fitted Hill equation with parameters determined by Mank et al.
    (Nature Methods, 5:801-811, 2008, https://doi.org/10.1038/nmeth.1243).

    Adds 4 parameters and 5 expressions.

    Parameters:
        * Ca_C_resting - resting amount of free cytosolic Ca2+ used when
            computing the background FRET ratio (in case it is not the same
            as the initial amount of free cytosolic Ca2+).
        * Kd_cytCa_bind_TNXXL - the effective dissociation constant for
            calcium binding to TN-XXL.
        * Rmax - the maximum value of the FRET ratio.
        * HillCoeff_TNXXL - the Hill Coefficient for the FRET ratio Hill curve.

    Expressions:
        * _Ca_C_resting - defines the resting amount of Ca2+ to use for the
            FRET computation. It's a piecewise function that is equal
            to Ca_C_resting when Ca_C_resting is different than Ca_C_0, or
            equal to Ca_C_0 otherwise.
        * Ca_num_to_microM - conversion factor to convert the amount of
            cytosolic Ca2+ into a uM concentration.
        * Frc_base - the FRET ratio at the resting amount of free cytosolic
            Ca2+ (relative to zero Ca2+).
        * Frc_curr - the FRET ratio at the current amount of free cytosolic
            Ca2+ (relative to zero Ca2+).
        * FRET - the FRET ratio of the current Ca2+ amount relative to the
            resting amount of Ca2+. This corresponds to the experimentally
            observable FRET ratio.
    """
    alias_model_components()
    # Initial concentration of Ca2+ in the cytosol expected to be around 100 nM.
    Parameter("Ca_C_resting", Ca_C_0.value)
    alias_model_components()
    Expression(
        "_Ca_C_resting",
        Piecewise(
            (Ca_C_resting, Ca_C_0 > Ca_C_resting),
            (Ca_C_resting, Ca_C_0 < Ca_C_resting),
            (Ca_C_0, True),
        ),
    )
    # The maximum FRET ratio, deltaR/R, for TN-XXL is 2.3 at 39 microM Ca2+,
    # the effective Kd for Ca2+ binding to TN-XXL FRET reporter is
    #  Kd = 800 nM,and the Hill-Coefficient is 1.5, https://doi.org/10.1038/nmeth.1243
    Parameter("Kd_cytCa_bind_TNXXL", 800e-3)  # microM
    Parameter("Rmax", 2.3)
    Parameter("HillCoeff_TNXXL", 1.5)
    alias_model_components()
    # Get the FRET signal

    # Compute the FRET ratio change relative to zero (i.e., Rmin) using the Hill equation,
    #    (R-Rmin)/Rmin = Rmax*[Ca2+]**h / (Kd + [Ca2+]**h) ,
    #       where Rmax is maximum FRET ratio change at saturation, h is the
    #       Hill Coefficient, and Kd is effective dissociation constant.
    Expression("Ca_num_to_microM", 1 / (Vcyto * units.microM_to_molec_per_pL))
    alias_model_components()
    # FRET ratio change for baseline concentration relative to zero - dR/R = (Rb-Rmin)/Rmin
    Expression(
        "Frc_base",
        -Rmax
        * (_Ca_C_resting * Ca_num_to_microM) ** HillCoeff_TNXXL
        / (Kd_cytCa_bind_TNXXL + (_Ca_C_resting * Ca_num_to_microM) ** HillCoeff_TNXXL),
    )
    # FRET ratio change for current concentration relative to zero - dR/R = (Rc-Rmin)/Rmin
    Expression(
        "Frc_curr",
        Rmax
        * (cytoCa * Ca_num_to_microM) ** HillCoeff_TNXXL
        / (Kd_cytCa_bind_TNXXL + (cytoCa * Ca_num_to_microM) ** HillCoeff_TNXXL),
    )
    alias_model_components()
    # Exp. FRET ratio change which is relative to the baseline - dR/R = (Rc-Rb)/Rb
    Expression("FRET", (Frc_curr + Frc_base) / (-Frc_base + 1))
    return


def degradation_of_excess_cytosolic_calcium():
    alias_model_components()
    Parameter("k_cytoCa_loss", 4)
    alias_model_components()

    Expression(
        "k_cytoCa_loss_piecewise",
        Piecewise(
            (k_cytoCa_loss, cytoCa > Ca_C_resting),
            (0, True),
        ),
    )
    alias_model_components()
    degrade(Ca(b=None) ** CYTOSOL, k_cytoCa_loss_piecewise)
    return

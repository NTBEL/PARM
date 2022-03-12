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
import numpy as np

from . import defaults, units


def calcium_signal_monomers():
    # Phospholipase C
    Monomer("PLC", ["bgaq", "bpip2"])
    # PIP2
    Monomer("PIP2", ["b"])
    # IP3
    Monomer("IP3", ["b"])
    # DAG
    Monomer("DAG")
    # IP3 receptor
    Monomer("IP3R", ["b1", "b2", "b3", "b4", "bcaer", "bcacyt"])
    # Calcium 2+, loc: E = ER space, C = cytosol
    Monomer("Ca", ["b"])
    alias_model_components()

    # Annotations
    # ===========
    Annotation(PLC, "https://identifiers.org/uniprot:Q9NQ66")
    Annotation(PIP2, "https://identifiers.org/CHEBI:145879", predicate="hasVersion")
    Annotation(DAG, "https://identifiers.org/CHEBI:85722", predicate="hasVersion")
    Annotation(IP3R, "https://identifiers.org/uniprot:Q14643")
    Annotation(Ca, "https://identifiers.org/CHEBI:29108")
    return


def calcium_signal_initials():
    alias_model_components()
    # inactive PLC
    # From Falkenburger et al. 2010 https://dx.doi.org/10.1085%2Fjgp.200910344
    # tsA201 cells
    # Endogenous PLCB1 concentration: 3/micrometer^2
    # Overexpressed PLCB1 concentration: 3,000/micrometer^2
    # PLC total endogenous: 10/micrometer^2
    Parameter("PLC_0", 3 * SAcell.value)
    # PIP2
    # Basal no. of PIP2 molecules is 49997 as per Lemon et al. 2003 https://doi.org/10.1016/S0022-5193(03)00079-1
    # also free PIP2 of 5000 per micrometer^2 used by Falkenburger et al. 2013 https://doi.org/10.1085/jgp.201210887
    # For nominal value will start with Lemon et al. value.
    Parameter("PIP2_0", 49997)
    # IP3R
    Parameter("IP3R_0", defaults.SPC * Ver.value)
    # ER Ca2+ store
    # ER lumen of HEK-293 cells has between roughly 400-600 microM with an average
    # around 525 microM as reported in
    # Foyouzi-Youssefi et al. https://doi.org/10.1073/pnas.97.11.5723 (Fig. 3C, control)
    Parameter("Ca_E_0", 525 * units.microM_to_molec_per_pL * Ver.value)
    # Initial concentration of Ca2+ in the cytosol expected to be around 100 nM.
    Parameter("Ca_C_0", 100 * units.nM_to_molec_per_pL * Vcyto.value)
    # In MH experiments the extracellular space is filled with ACSF with 3.1 mM
    # CaCl2, so extracellular Ca2+ should be around 3.1 mM.
    # Kang et al. https://doi.org/10.1021/acsnano.9b01993
    Parameter("Ca_extra_0", 3.1 * 1e3 * units.microM_to_molec_per_pL * Vextra.value)
    alias_model_components()

    Initial(PLC(bgaq=None, bpip2=None) ** CELL_MEMB, PLC_0)
    Initial(PIP2(b=None) ** CELL_MEMB, PIP2_0)
    Initial(
        IP3R(b1=None, b2=None, b3=None, b4=None, bcaer=None, bcacyt=None) ** ER_MEMB,
        IP3R_0,
    )
    Initial(Ca(b=None) ** ER_LUMEN, Ca_E_0)
    Initial(Ca(b=None) ** CYTOSOL, Ca_C_0)
    Initial(Ca(b=None) ** EXTRACELLULAR, Ca_extra_0)
    return


def plc_binds_gaq_and_catalyzes_pip2_to_ip3():
    alias_model_components()
    # PLC binding Gaq
    Parameter("kf_PLC_bind_Gaq", defaults.KF_BIND / Vcyto.value)
    Parameter("kr_PLC_bind_Gaq", defaults.KR_BIND)
    # Conversion of PIP2 to IP3
    Parameter("kf_PLC_bind_PIP2", defaults.KF_BIND / Vcyto.value)
    Parameter("kr_PLC_bind_PIP2", defaults.KR_BIND)
    Parameter("kcat_PIP2_to_IP3", defaults.KCAT)
    alias_model_components()
    # PLC activation by binding Gaq:
    #    Gaq_A + PLC <---> Gaq_A:PLC
    #   Reusing the Gbg binding slot for PLC
    # Alias the complex Gaq:GTP
    Gaq_gtp = Gaq(bpar=None, bgdp=3, bgbg=None) ** CELL_MEMB % GTP(b=3) ** CELL_MEMB
    bind_complex(
        Gaq_gtp,
        "bgbg",
        PLC(bgaq=None, bpip2=None) ** CELL_MEMB,
        "bgaq",
        [kf_PLC_bind_Gaq, kr_PLC_bind_Gaq],
    )
    # Conversion of PIP2 to IP3
    #    Gaq_A:PLC + PIP2 <---> Gaq_A:PLC:PIP2 ---> Gaq_A:PLC + IP3
    Gaq_gtp_PLC = (
        Gaq(bpar=None, bgdp=3, bgbg=1) ** CELL_MEMB
        % GTP(b=3) ** CELL_MEMB
        % PLC(bgaq=1, bpip2=None) ** CELL_MEMB
    )
    Gaq_gtp_PLC_pip = (
        Gaq(bpar=None, bgdp=3, bgbg=1) ** CELL_MEMB
        % GTP(b=3) ** CELL_MEMB
        % PLC(bgaq=1, bpip2=2) ** CELL_MEMB
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
        Gaq_gtp_PLC_pip >> Gaq_gtp_PLC + IP3(b=None) ** CYTOSOL + DAG() ** CELL_MEMB,
        kcat_PIP2_to_IP3,
    )

    return


def ip3_binds_ip3r():
    alias_model_components()
    # Binding of IP3 to IP3R
    R_o = 2e-7  # cm
    # IP3 diffuses in mammalian at <= 10 micrometer^2/s https://dx.doi.org/10.1126%2Fscisignal.aag1625
    D_ip3 = 10e-8  # cm^2/s
    # Assume the IP3 binding rate is diffusion-controlled by IP3 diffusion
    K_IP3_BIND = 4 * np.pi * D_ip3 * R_o * (1e-3) / (Vcyto.value * 1e-12)
    Parameter("kf_IP3_bind_IP3R", K_IP3_BIND)
    Parameter("kr_IP3_bind_IP3R", defaults.KR_BIND)
    alias_model_components()
    # Binding of IP3 to IP3R
    # Assume subunits are bound sequentially and that there is no cooperativity
    # between subunit binding.
    #   IP3R + IP3 <---> IP3R:IP3, subunit 1
    bind(
        IP3R(b2=None, b3=None, b4=None, bcaer=None, bcacyt=None) ** ER_MEMB,
        "b1",
        IP3(b=None) ** CYTOSOL,
        "b",
        [kf_IP3_bind_IP3R, kr_IP3_bind_IP3R],
    )
    #   IP3R + IP3 <---> IP3R:IP3, subunit 2
    Rule(
        "bind_IP3_IPR3_sub2",
        IP3R(b1=1, b2=None, b3=None, b4=None, bcaer=None, bcacyt=None) ** ER_MEMB
        % IP3(b=1) ** CYTOSOL
        + IP3(b=None) ** CYTOSOL
        | IP3R(b1=1, b2=2, b3=None, b4=None, bcaer=None, bcacyt=None) ** ER_MEMB
        % IP3(b=1) ** CYTOSOL
        % IP3(b=2) ** CYTOSOL,
        kf_IP3_bind_IP3R,
        kr_IP3_bind_IP3R,
    )
    #   IP3R + IP3 <---> IP3R:IP3, subunit 3
    Rule(
        "bind_IP3_IPR3_sub3",
        IP3R(b1=1, b2=2, b3=None, b4=None, bcaer=None, bcacyt=None) ** ER_MEMB
        % IP3(b=1) ** CYTOSOL
        % IP3(b=2) ** CYTOSOL
        + IP3(b=None) ** CYTOSOL
        | IP3R(b1=1, b2=2, b3=3, b4=None, bcaer=None, bcacyt=None) ** ER_MEMB
        % IP3(b=1) ** CYTOSOL
        % IP3(b=2) ** CYTOSOL
        % IP3(b=3) ** CYTOSOL,
        kf_IP3_bind_IP3R,
        kr_IP3_bind_IP3R,
    )
    #   IP3R + IP3 <---> IP3R:IP3, subunit 4
    Rule(
        "bind_IP3_IPR3_sub4",
        IP3R(b1=1, b2=2, b3=3, b4=None, bcaer=None, bcacyt=None) ** ER_MEMB
        % IP3(b=1) ** CYTOSOL
        % IP3(b=2) ** CYTOSOL
        % IP3(b=3) ** CYTOSOL
        + IP3(b=None) ** CYTOSOL
        | IP3R(b1=1, b2=2, b3=3, b4=4, bcaer=None, bcacyt=None) ** ER_MEMB
        % IP3(b=1) ** CYTOSOL
        % IP3(b=2) ** CYTOSOL
        % IP3(b=3) ** CYTOSOL
        % IP3(b=4) ** CYTOSOL,
        kf_IP3_bind_IP3R,
        kr_IP3_bind_IP3R,
    )
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
    # is 525 1/s
    Parameter("k_tranport_erCa", 525)
    alias_model_components()
    # Transport of Ca2+ by activated IP3R
    #  ER -> cytosol:
    #    IP3R:IP3_4 + Ca_E <---> Ca_E:IP3R:IP3_4 ---> Ca_C + IP3R:IP3_4
    Rule(
        "bind_Ca_IPR3_er",
        IP3R(b1=1, b2=2, b3=3, b4=4, bcaer=None, bcacyt=None) ** ER_MEMB
        % IP3(b=1) ** CYTOSOL
        % IP3(b=2) ** CYTOSOL
        % IP3(b=3) ** CYTOSOL
        % IP3(b=4) ** CYTOSOL
        + Ca(b=None) ** ER_LUMEN
        | IP3R(b1=1, b2=2, b3=3, b4=4, bcaer=5, bcacyt=None) ** ER_MEMB
        % IP3(b=1) ** CYTOSOL
        % IP3(b=2) ** CYTOSOL
        % IP3(b=3) ** CYTOSOL
        % IP3(b=4) ** CYTOSOL
        % Ca(b=5) ** ER_LUMEN,
        kf_erCa_bind_IP3R,
        kr_erCa_bind_IP3R,
    )
    Rule(
        "transport_Ca_ER_CYTO",
        IP3R(b1=1, b2=2, b3=3, b4=4, bcaer=5, bcacyt=None) ** ER_MEMB
        % IP3(b=1) ** CYTOSOL
        % IP3(b=2) ** CYTOSOL
        % IP3(b=3) ** CYTOSOL
        % IP3(b=4) ** CYTOSOL
        % Ca(b=5) ** ER_LUMEN
        >> IP3R(b1=1, b2=2, b3=3, b4=4, bcaer=None, bcacyt=None) ** ER_MEMB
        % IP3(b=1) ** CYTOSOL
        % IP3(b=2) ** CYTOSOL
        % IP3(b=3) ** CYTOSOL
        % IP3(b=4) ** CYTOSOL
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


def cytosolic_calcium_regulation():
    """Reactions to regulate the concentration of free cytosolic calcium.
    The regulation of free cytosolic calcium is maintained by two first order
    reactions: one for transfer of free cytosolic calcium to
    the extracellular space and one for transfer of extracellular calcium into to the
    cytosol. The cytosol-to-extracellular reaction represent extrusion of
    calcium from the cell while the extracellular-to-ctyosol reaction represents
    influx of calcium from the extracellular space. The two rates cancel one
    another out so that the resting free cytosolic calcium concentration is
    maintained in the abscence of any perturbations that
    allow release from the ER store.

    Reactions:
        1. Ca2+_CYTO ---> Ca2+_EXTRA
        2. Ca2+_EXTRA ---> Ca2+_CYTO
    """
    # Cytosolic Ca2+ regulation

    # cytosol to extracellular space
    # From previous model fittings we get around 4 1/s.
    Parameter("k_Ca_cyt_to_extra", 4)  # 1/s
    alias_model_components()
    # extracellular space to cytosol
    # Assume that at equilibrium before any agonist is added the rate out of the
    # cytosol equals the rate in. We'll enforce this by setting the extra to cyt
    # rate constant with an Expression.
    # k_Ca_cyt_to_extra*Ca_C_0/(Ca_extra_0 * Vextra**2)
    Expression("k_Ca_extra_to_cyt", k_Ca_cyt_to_extra * Ca_C_0 / Ca_extra_0)
    alias_model_components()
    # Regulation of Cytosolic Ca2+ --
    # Here, we assume that the decay in the FRET signal is primarily due to
    # first-order excretion of excess Ca2+ from the released ER store into the
    # extracellular space by cell membrane ion channels. We assume both the
    # excretion from the cytosol to the extracellular space and the reverse
    # influx from the extracellular space to the cytosol are first order reactions
    # that are at equilibrium before any agonist is added.

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


def gaq_activated_calcium_signaling():
    calcium_signal_monomers()
    calcium_signal_initials()
    plc_binds_gaq_and_catalyzes_pip2_to_ip3()
    ip3_binds_ip3r()
    ip3r_transports_er_calcium_to_cytosol()
    cytosolic_calcium_regulation()
    ip3_degradation()
    return


def observables():
    # The maximum FRET ratio, deltaR/R, for TN-XXL is 2.3 at 39 microM Ca2+,
    # the effective Kd for Ca2+ binding to TN-XXL FRET reporter is
    #  Kd = 800 nM,and the Hill-Coefficient is 1.5, https://doi.org/10.1038/nmeth.1243
    Parameter("Kd_cytCa_bind_TNXXL", 800e-3)  # microM
    Parameter("Rmax", 2.3)
    Parameter("HillCoeff_TNXXL", 1.5)
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
    # Ca2+ in the Cytosol
    Observable("cytoCa", Ca(b=None) ** CYTOSOL)
    alias_model_components()

    Expression("active_IP3R_ratio", aIP3R / totIP3R)
    Expression("Ca_num_to_nM", 1 / (Vcyto * units.nM_to_molec_per_pL))
    alias_model_components()
    Expression("cytoCa_nM", cytoCa * Ca_num_to_nM)
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
        * (Ca_C_0 * Ca_num_to_microM) ** HillCoeff_TNXXL
        / (Kd_cytCa_bind_TNXXL + (Ca_C_0 * Ca_num_to_microM) ** HillCoeff_TNXXL),
    )
    # FRET ratio change for current concentration relative to zero - dR/R = (Rc-Rmin)/Rmin
    Expression(
        "Frc_curr",
        Rmax
        * (cytoCa * Ca_num_to_microM) ** HillCoeff_TNXXL
        / (Kd_cytCa_bind_TNXXL + (cytoCa * Ca_num_to_microM) ** HillCoeff_TNXXL),
    )
    alias_model_components()
    Expression("Fr_diff", Frc_curr + Frc_base)
    alias_model_components()
    # Exp. FRET ratio change which is relative to the baseline - dR/R = (Rc-Rb)/Rb
    Expression("FRET", (Frc_curr + Frc_base) / (-Frc_base + 1))

    return

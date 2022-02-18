"""PARM: PAR2 Activation-driven calcium Release Model
PARM is a mechanistic mass action model of PAR2 activation and downstream
calcium signaling via the phospholipase C and IP3 pathway. The model was
designed to mathematically model the underlying signaling dynamics
relevant to the inactivation of PAR2 by Molecular Hyperthermia as
described in:
Kang et al.,  Transient Photoinactivation of Cell Membrane Protein Activity
without Genetic Modification by Molecular Hyperthermia, ACS Nano 2019, 13, 11,
12487–12499 https://doi.org/10.1021/acsnano.9b01993
PAR2 activation and subsequent G-protein activation are modeled after the
Classical GPCR/G-protein activation model (e.g., see Fig 2A of Oliveira et al.
https://doi.org/10.3389/fnagi.2019.00089) in which G-protein heterotrimers
only interact with the receptor after receptor-activation. The model also
assumes that all four subunits of the IP3 receptor, IP3R, must be bound by IP3
before calcium can bind the receptor and be translocated between the ER lumen
and cytosol which is consistent with work by Alzayady et al
https://doi.org/10.1126/scisignal.aad6281. However,  feedback to either reduce
or enhance the IP3R calcium release is not included. Additionally, the
cytosolic calcium level maintenance is modeled by a unidirectional degradation
step that is active when the cytosolic calcium concentration goes  above the
starting value, approximating a lumped process for the regulation of Ca2+
concentration in the cytosol after the ER store is released (e.g., activation
of SERCA to pump Ca2+ back into the lumen, or activation of cell membrane ion
channels to release excess Ca2+ into the extracellular space). However, the
initial cytosolic concentration of Ca2+ is left at zero for puposes of the
differential equations in order to avoid having to define piecewise functions
for cytosolic calcium regulation, but the ctyosol is implicity assumed to have
a  baseline constant concentration of Ca2+ (nominally 100 nM) which is included
when computing the FRET ratio. Spefically, in this version of the model the
FRET ratio for cytosolic Ca2+ is estimated by using the Hill equation with
parameters from the dose-repsonse curve of TN-XXL to Ca2+.
The full set of interactions and sequence of rules included in the model are as
follows:
  1. Zero-order synthesis of PAR2:
      None ---> PAR2_I
  2. First order degradations of PAR2:
      i) free inactive PAR2:
         PAR2_I ---> None
     ii) bound inactive PAR2:
         TAT:PAR2_I ---> None
    iii) active PAR2:
         TAT:PAR2_A ---> None
     iv) denatured PAR2:
         PAR2_D ---> None
  Note: In the abscence of agonist and with no denatured PAR2 we expect 1 and
        2.i to cancel each other so that the net rate of change of PAR2_I is
        zero.
  3. Two-state receptor activation of PAR2 by 2AT:
      2AT + PAR2_I <---> TAT:PAR2_I <---> TAT:PAR2_A
  4. Gaq activation by activated-PAR2:  | Note: G-proteins are not pre-coupled to PAR2.
      i) G protein heterotrimer binds activated PAR2:
         PAR2_A + Gaq:GDP:Gbg <---> PAR2_A:Gaq:GDP:Gbg
     ii) GDP preferentially unbinds from Gaq:
         PAR2_A:Gaq:GDP:Gbc <---> PAR2_A:Gaq:Gbc + GDP
    iii) GTP preferentially binds Gaq:
         PAR2_A:Gaq:Gbc + GTP <---> PAR2_A:Gaq_A:GTP:Gbg
     iv) Gbg dissociates from Gaq (i.e., heterotrimer dissociation):
         PAR2_A:Gaq:GTP:Gbc ---> PAR2_A:Gaq:GTP + Gbc
      v) Gaq:GTP dissociates from PAR2, Gaq is now active (G protein dissociation from the receptor):
         PAR2_A:Gaq:GTP ---> PAR2_A + Gaq:GTP
  5. Hydrolosis of GTP by Gaq (inactivation of Gaq)
       a) Slow hydrolosis by Gaq alone
           Gaq:GTP ---> Gaq:GDP
       b) RGS enhanced hydrolosis
           Gaq:GTP + RGS <---> Gaq:GTP:RGS ---> Gaq:GDP + RGS
       c) PLC enhanced hydrolosis
           Gaq:GTP:PLC ---> Gaq:GDP + PLC
  6. Recombination of G protein heterotrimer
     Gaq:GDP + Gbg ---> Gaq:GDP:Gbq
  7. PLC activation by binding Gaq:
      Gaq_A:GTP + PLC <---> Gaq_A:GTP:PLC
  8. Conversion of PIP2 to IP3
      Gaq_A:GTP:PLC + PIP2 <---> Gaq_A:GTP:PLC:PIP2 ---> Gaq_A:PLC + IP3
  9. Binding of IP3 to IP3R:  | IP3R is only activated when all 4 subunits are bound.
     i) IP3R + IP3 <---> IP3R:IP3, subunit 1
    ii) IP3R + IP3 <---> IP3R:IP3, subunit 2
   iii) IP3R + IP3 <---> IP3R:IP3, subunit 3
    iv) IP3R + IP3 <---> IP3R:IP3, subunit 4
 10. Transport of Ca2+ by activated IP3R
    i) ER to cytosol:
       IP3R:IP3_4 + Ca_E <---> Ca_E:IP3R:IP3_4 ---> Ca_C + IP3R:IP3_4
   ii) Reverse, cytosol to ER:  | Assuming the transport is not just one way.
       IP3R:IP3_4 + Ca_C <---> Ca_C:IP3R:IP3_4 ---> Ca_E + IP3R:IP3_4
 11. Cytosolic calcium regualation
    i) First-order excretion of calcium to the extracellular space:
       Ca_C ---> Ca_EXTRA
   ii) First-order influx of extracellular calcium into the cytosol:
       Ca_EXTRA --> Ca_C
 Note: In the abscence of agonist we expect 11.i and 11.ii to cancel each other
       out so that the net rate of change of the cytosolic calcium is zero.
 12. Degradation of IP3
       IP3 ---> None
Unless otherwise noted the units used are:
    Volume : pL
    Area: micrometer^2
    Distance: micrometer
    Volume Concentration: number ( convert concentration by factor of V*N_A)
    Area Concentration: number/micrometer^2
    Forward bimolecular association rate constants (kf) : 1/(s*number)
    Reverse of bimolecular association (dissociation) rate constants (kr) : 1/s
    Catalytic rate constants (kcat) :  1/s
    Other unimolecular rate constants: 1/s
    Dissociation (bimolecular) constants: Kd : number
    Binding (bimolecular) constants : Kb : 1/number
"""

# PySB components
from pysb import (
    Model,
    Monomer,
    Parameter,
    Initial,
    Rule,
    Observable,
    Expression,
    Annotation,
    Compartment,
    ANY,
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

# NumPy
import numpy as np

# Avogadro's Number from scipy
from scipy.constants import N_A
from sympy.functions.elementary.miscellaneous import Max

# Conversion factors for concentration units.
# microMolar to number/pL
microM_to_num_per_pL = 1e-6 * N_A * 1e-12
# nanoMolar to number/pL
nM_to_num_per_pL = 1e-9 * N_A * 1e-12
# Cubic micron to picoliter:
#               cubic-micron to mL * mL to L * L to pL
cubicmicron_to_pL = (1e-4) ** 3 * 1e-3 * 1e12

Model()

# Cellular volume, 10^-12 L as assumed in
# Albeck et al. https://doi.org/10.1371/journal.pbio.0060299
Parameter("Vcell", 1)
# Cell-membrane surface area
# Note: Hek cells ~706.5 micrometer^2 https://doi.org/10.1038/s41598-017-07813-5
# But we'll define the effective surface area to be a function of the radius
# which corresponds to the volume.
Rcell = (3 / 4 * np.pi * Vcell.value / cubicmicron_to_pL) ** (1 / 3)
# Get the cell surface area using the radius.
Parameter("SAcell", 4 * np.pi * Rcell ** 2)
# Effective cell-membrane thickness
# Assume 10 nm (0.01 micron) as in https://github.com/RuleWorld/BNGTutorial/blob/master/CBNGL/LR_comp.bngl
Parameter("CMthickness", 0.01)
# Effective volume of the cell-membrane
Parameter("Vcm", SAcell.value * CMthickness.value * cubicmicron_to_pL)

# Volume of the extracellular space
# The following BNGL examples use 1000x the cell volume:
#   https://github.com/RuleWorld/BNGTutorial/blob/master/CBNGL/LRR_comp.bngl
#   https://github.com/RuleWorld/BNGTutorial/blob/master/CBNGL/LR_comp.bngl
# However, we'll just set the effective extracellular reaction volume to be on par
# (pun not intended, but acknowledged) with the cellular volume to avoid overly large
# numbers of agonist molecules.
Parameter("Vextra", Vcell.value)

# Volume of the ER lumen/cisternal space.
# It is often >10% of cell volume according Alberts et al. https://www.ncbi.nlm.nih.gov/books/NBK26841/ .
# Lemon et al. 2003 use ratio of 0.185 ER lumen/cytosol
# Politi et al. 2006 also use ratio of 0.185 ER lumen/cytosol
# We can also use 0.185
Parameter("Ver", Vcell.value * 0.185)
# Assume 10x the cell membrane surface area
# Parameter("SAer", SAcell.value*10)
Rer = (3 / 4 * np.pi * Ver.value / cubicmicron_to_pL) ** (1 / 3)
Parameter("SAer", 4 * np.pi * Rer ** 2)
# Effective thickness of ER membrane
# Assume 10 nm is still reasonable
Parameter("ERMthickness", 0.01)
# Effective Volume of the ER membrane
Parameter("Verm", ERMthickness.value * SAer.value * cubicmicron_to_pL)

# Default forward, reverse, and catalytic rates:
# KF_BIND is equivalent to kf/(Vcell*1e-12) / N_A for kf in 1/(M*s)
KF_BIND = 1e-6  # 1/(number*s) Default forward binding rate (for cell volume of 1 pL) from Aldridge et al. https://doi.org/10.1038/ncb1497
KR_BIND = 1e-3  # Default dissociation rate from Albeck et al. https://doi.org/10.1371/journal.pbio.0060299
KCAT = 10  # "average enzyme" from Bar-Even et al. https://doi.org/10.1021/bi2002289

# Default signaling protein concentration range
# of 1 nM to 1 microM range assumed by Albeck et al. https://doi.org/10.1371/journal.pbio.0060299
# for which they reference Wu and Pollard https://doi.org/10.1126/science.1113230
# Also see Aldridge et al. https://doi.org/10.1038/ncb1497
# Set to 100 nM or 0.1 microM
SPC = 0.1 * microM_to_num_per_pL

# As a first estimate for the forward binding reactions of Ca2+ we will assume
# that it is diffusion-controlled following Smoluchowski eqn.:
#     kf = 4*pi*D*R_o
# We will assume R_o is 2 nm and that D = D_Ca2+.
# The diffusion coefficient of Ca2+ is 5.3 x 10^-6 cm^2/s; https://doi.org/10.1016/0143-4160(87)90027-3
D_Ca = 5.3e-6  # cm^2/s
R_o = 2e-7  # cm
# (1e-3) term is for unit conversion from cm^3/s*number to 1/s*(number/L)
# mL/s*number -> 10^-3 L/(s*number) and dividing by V converts to 1/(s*number)
K_CA_BIND = 4 * np.pi * D_Ca * R_o * (1e-3) / (Vcell.value * 1e-12)

# Ion channel transport rate: up to 1e8 ions/s https://www.ncbi.nlm.nih.gov/books/NBK26910/
K_ION_CHANNEL = 1e8

# IP3 diffuses in mammalian at <= 10 micrometer^2/s https://dx.doi.org/10.1126%2Fscisignal.aag1625
D_ip3 = 10e-8  # cm^2/s
# Assume the IP3 binding rate is diffusion-controlled by IP3 diffusion
K_IP3_BIND = 4 * np.pi * D_ip3 * R_o * (1e-3) / (Vcell.value * 1e-12)

# Default molecule degradation rate.
K_DEGRADE = 1  # 1/s

# Default unidirection conversion rate.
K_CONVERT = 1  # 1/s

# Compartments
# ============
# Note: "For  elementary  bimolecular  reactions  with  the  rate  constant  given
#  in  units  of volume/time,  the  rate  constant  is  divided  by the  volume
#  of  the  reactant  compartment  with  the  highest  dimension"
# -- https://www.informs-sim.org/wsc09papers/087.pdf
# Although the above says "volume/time" for units of the forward binding rate,
# as far as I can tell it is actually volume/time*number so that dividing by
# the compartment volume would give you 1/time*number (e.g., 1/s*number),
# corresponding to concentrations given in number per cell.
# Also see BNGL example: https://github.com/RuleWorld/BNGTutorial/blob/master/CBNGL/LRR_comp.bngl

# Since we are already converting concentrations to numbers per cell and
# the bimolecular binding rates are already relative to the cell volume
# of 1 pL we can scale all compartment sizes relative to the cell/cytosol
# volume.
# When comparment volume scaling is applied it should yield:
#       KF_BIND/(Vcompartment/Vcell) = KF_BIND*Vcell/Vcomparment
# so that KF_BIND*Vcell is the binding rate in pL/s*number and then
# division by Vcompartment returns the compartment-specific scaled binding
# rate in 1/s*number.
Parameter("V_EXTRA", Vextra.value / Vcell.value)
Compartment("EXTRACELLULAR", dimension=3, size=V_EXTRA)
# Cell Membrane
Parameter("V_CM", Vcm.value / Vcell.value)
Compartment("CELL_MEMB", dimension=2, parent=EXTRACELLULAR, size=V_CM)
# Cytosol
Parameter("V_C", Vcell.value / Vcell.value)
Compartment("CYTOSOL", dimension=3, parent=CELL_MEMB, size=V_C)
#  ER membrane
Parameter("V_ERM", Verm.value / Vcell.value)
Compartment("ER_MEMB", dimension=2, parent=CYTOSOL, size=V_ERM)
# ER lumen volume
Parameter("V_ERL", Ver.value / Vcell.value)
Compartment("ER_LUMEN", dimension=3, parent=ER_MEMB, size=V_ERL)

# Monomers
# ========
# PAR2 agonist 2AT
# Note: this was modeled as Trypsin ('Tryp') in an earlier version of the model
# which is a canonical protease activator of PAR2 (via N-terminal cleavage),
# but in the experiments of Kang et al. PAR2 is actually activated by the
# agonist 2AT, Kang et al. https://doi.org/10.1021/acsnano.9b01993)
Monomer("TAT", ["b"])
# PAR2, states: I = inactive, A = active, D = denatured
Monomer("PAR2", ["bortho", "bgaq", "state"], {"state": ["I", "A", "D"]})
# G-alpha_q G-protein unit, states: I = inactive, A = active
Monomer("Gaq", ["bpar", "bgbg", "bgdp"])
# G-beta-gamma G-protein units
Monomer("Gbg", ["b"])
# GDP
Monomer("GDP", ["b"])
# GTP
Monomer("GTP", ["b"])
# Regulator of G protein Signaling (RGS)
Monomer("RGS", ["b"])

# Phospholipase C
Monomer("PLC", ["bgaq", "bpip2"])
# PIP2
Monomer("PIP2", ["b"])
# IP3
Monomer("IP3", ["b"])
# IP3 receptor
Monomer("IP3R", ["b1", "b2", "b3", "b4", "bcaer", "bcacyt"])
# Calcium 2+, loc: E = ER space, C = cytosol
Monomer("Ca", ["b", "loc"], {"loc": ["E", "C"]})

# Annotations
# ===========
Annotation(PAR2, "https://identifiers.org/uniprot:P55085")
Annotation(Gaq, "https://identifiers.org/uniprot:P50148")
Annotation(PLC, "https://identifiers.org/uniprot:Q9NQ66")
Annotation(IP3R, "https://identifiers.org/uniprot:Q14643")

# Initial conditions
# ==================
# PAR2 agonist 2AT, 330 nM for Fig 2D data from Kang et al. https://doi.org/10.1021/acsnano.9b01993
C_2AT = 330  # nM
V_2AT = 50e6  # Volume of agonist added to wells is 50 microL
Vwell = 150e6  # Looks like the total well volume was 150 microL (100 microL ACSF + 50 microL agonist in ACSF)
nM_2AT_to_num = nM_to_num_per_pL * (V_2AT / Vwell) * Vextra.value
# nM_2AT_to_molec = 1e-9 * V_2AT * N_A
Parameter("TAT_0", C_2AT * nM_2AT_to_num)
Initial(TAT(b=None) ** EXTRACELLULAR, TAT_0)
# total PAR2
# From Falkenburger et al. 2010 https://dx.doi.org/10.1085%2Fjgp.200910344
# # tsA201 cells
# Endogenous receptor density: 1/micrometer^2
# Overexpressed receptor density: 3,000/micrometer^2
# From Brinkerhoff et al. 2008: Receptor concentration 2e3 to 2e4 /cell
Parameter("PAR2_0", 1 * SAcell.value)
Parameter("f_denature", 0.0)  # By default no PAR2 has been denatured.
Expression("PAR2_0_D", f_denature * PAR2_0)
Expression("PAR2_0_I", Max(PAR2_0 - PAR2_0_D, 0))
Initial(PAR2(state="I", bortho=None, bgaq=None) ** CELL_MEMB, PAR2_0_I)
Initial(PAR2(state="D", bortho=None, bgaq=None) ** CELL_MEMB, PAR2_0_D)
# inactive G-protein heterotrimer Gaq-GDP:Gbg (the beta and gamma units are modeled as a single unit)
# From Falkenburger et al. 2010 https://dx.doi.org/10.1085%2Fjgp.200910344
# Endogenous G-protein density: 40/micrometer^2
# Overexpressed G-protein density: 3,000/micrometer^2
# Brinkerhoff et al. 2008: G-protein concentration 1e4 /cell
Parameter("Gaq_0", 40 * SAcell.value)
# Alias the free Gprotein heterotrimer
Gaq_gdp_Gbg = (
    Gaq(bpar=None, bgbg=3, bgdp=4) ** CELL_MEMB
    % GDP(b=3) ** CELL_MEMB
    % Gbg(b=4) ** CELL_MEMB
)
Initial(Gaq_gdp_Gbg, Gaq_0)
# GTP
# Physiolocal concentration of GTP in mammalian cells is generally 468 +/- 224 microM
# and for human cells it is 305 microM
# as per Traut https://doi.org/10.1007/bf00928361
Parameter("GTP_0", 305 * microM_to_num_per_pL * Vcell.value)
Initial(GTP(b=None) ** CYTOSOL, GTP_0)
# GDP
# Physiolocal concentration of GDP in human cells is 36 microM
# as per Traut https://doi.org/10.1007/BF00928361
Parameter("GDP_0", 36 * microM_to_num_per_pL * Vcell.value)
Initial(GDP(b=None) ** CYTOSOL, GDP_0)
# RGS
# For RAW 264.7 Cell model between 0.008 and 0.012 microM as
# per Maurya and Subramaniam 2007 https://doi.org/10.1529/biophysj.106.097469
# Assume 0.010 microM is a reasonable starting point.
Parameter("RGS_0", 0.010 * microM_to_num_per_pL * Vcell.value)
Initial(RGS(b=None) ** CYTOSOL, RGS_0)

# inactive PLC
# From Falkenburger et al. 2010 https://dx.doi.org/10.1085%2Fjgp.200910344
# tsA201 cells
# Endogenous PLCB1 concentration: 3/micrometer^2
# Overexpressed PLCB1 concentration: 3,000/micrometer^2
# PLC total endogenous: 10/micrometer^2
Parameter("PLC_0", 3 * SAcell.value)
Initial(PLC(bgaq=None, bpip2=None) ** CELL_MEMB, PLC_0)
# PIP2
# Basal no. of PIP2 molecules is 49997 as per Lemon et al. 2003 https://doi.org/10.1016/S0022-5193(03)00079-1
# also free PIP2 of 5000 per micrometer^2 used by Falkenburger et al. 2013 https://doi.org/10.1085/jgp.201210887
# For nominal value will start with Lemon et al. value.
Parameter("PIP2_0", 49997)  # Had to convert to area concentration.
Initial(PIP2(b=None) ** CELL_MEMB, PIP2_0)
# IP3R
Parameter("IP3R_0", SPC * Ver.value)
Initial(
    IP3R(b1=None, b2=None, b3=None, b4=None, bcaer=None, bcacyt=None) ** ER_MEMB, IP3R_0
)
# ER Ca2+ store
# ER lumen of HEK-293 cells has between roughly 400-600 microM with an average
# around 525 microM as reported in
# Foyouzi-Youssefi et al. https://doi.org/10.1073/pnas.97.11.5723 (Fig. 3C, control)
Parameter("Ca_0", 525 * microM_to_num_per_pL * Ver.value)
Initial(Ca(loc="E", b=None) ** ER_LUMEN, Ca_0)
# Initial concentration of Ca2+ in the cytosol expected to be around 100 nM.
Parameter("Ca_C_0", 100 * nM_to_num_per_pL * Vcell.value)
Initial(Ca(loc="E", b=None) ** CYTOSOL, Ca_C_0)
# In MH experiments the extracellular space is filled with ACSF with 3.1 mM
# CaCl2, so extracellular Ca2+ should be around 3.1 mM.
# Kang et al. https://doi.org/10.1021/acsnano.9b01993
Parameter("Ca_extra_0", 3.1 * 1e3 * microM_to_num_per_pL * Vextra.value)

# Kinetic Parameters
# ==================
# PAR2 synthesis
# rate of receptor synthesis from
# Yi et al. 2003 PNAS https://doi.org/10.1073/pnas.1834247100
# used in yeast G-protein cycle model is: 4 number/s
Parameter("k_PAR2_synthesis", 4)
# Rate constant for degradation of PAR_I - Use an expression to enforce the
# assumption that at equilibrium the net rate of PAR2_I change due to synthesis
# and degradation is zero in the abscence of agonist or any PAR2 denaturation by
# MH.
Expression("k_PAR2_I_degradation", k_par2_synthesis / PAR2_0)
# rate constant for ligand-bound receptor degradation from
# Yi et al. 2003 PNAS https://doi.org/10.1073/pnas.1834247100
# used in yeast G-protein cycle model is: 4x10^-3 1/s
Parameter("k_PAR2_bound_degradation", 4e-3)
# As a default we can set the denatured PAR2 degradation to have the same
# rate constant as bound PAR2 degradation, but we'll assume it could be
# different.
Parameter("k_PAR2_denatured_degradation", 4e-3)
# PAR2 activation by 2AT
# Note: Ca2+ signal Max. FRET Dose-Response for 2AT activation of PAR2
# has EC50 = 101.7 +- 28.7 nM, Kang et al. https://doi.org/10.1021/acsnano.9b01993
# Binding forward rate constant for ligand-receptor binding in the yeast
# G-protein cycle model of  Yi et al. 2003 PNAS https://doi.org/10.1073/pnas.1834247100
# is 2x10^6 1/(M*s). Converting for our model that would be:
#    2x10^6 1/(M*s) / (10^-12 L) / N_A = 3x10^-6 1/(number*s) = 3*KF_BIND
# We'll use that value as our initial estimate.
Parameter("kf_PAR2_bind_TAT", KF_BIND * 3)
# PAR2 agonists in HEK 293T cells - LeSarge et al. https://doi.org/10.1021/acsmedchemlett.9b00094
#   2f-LIGRLO(Sulfo-Cy5)-NH2 has Kd = 430 nM with EC50 = 296 nM
#   Isox-Cha-Chg-ARK(Sulfo-Cy5)-NH2 has Kd = 38 nM with EC50 = 16 nM
# Since 2AT has EC50 = 101.8 nM in Hek 293 cells its probably safe to
# assume that the Kd for 2AT is somewhere between those two compounds.
# 142 = (430-38)/(296-16) *101.5
Parameter("Kd_PAR2_bind_TAT", 142 * nM_to_num_per_pL * Vcell.value)
Expression("kr_PAR2_bind_TAT", Kd_PAR2_bind_TAT * kf_PAR2_bind_TAT)
Parameter("k_activate_PAR2", K_CONVERT * 10)
Parameter("k_inactivate_PAR2", K_CONVERT / 10)
# Gaq binding activated-PAR2
Parameter("kf_PAR2_bind_Gaq", KF_BIND)
Parameter("kr_PAR2_bind_Gaq", KR_BIND)
# Gaq release GDP
Parameter("k_gdp_release", KR_BIND * 100)
Parameter("k_gdp_bind", KF_BIND)
# Gaq bind GTP
Parameter("k_gtp_bind", KF_BIND)
Parameter("k_gtp_release", KR_BIND / 10)
# Gbg dissociates from Gaq
Parameter("k_gbg_release", K_CONVERT)
# Gaq:GTP dissociates from PAR2
Parameter("k_gaq_release", K_CONVERT)
# Hydrolosis of GTP bound to Gaq
# 1. Autocatalysis rate for Gaq is ~0.8 1/min = 0.0133 1/s
# Bernstein et al. https://doi.org/10.1016/0092-8674(92)90165-9
# Also see Sprang https://dx.doi.org/10.1002%2Fbip.22836
Parameter("k_gtp_to_gdp_auto", 1.33e-2)
# 2. RGS binding and enhanced conversion of GTP to GDP
Parameter("kf_rgs_bind_gaq", KF_BIND)
Parameter("kr_rgs_bind_gaq", KR_BIND)
# Enhanced catalysis rate is for Gaq catalysis by RGS4
# is 100x higher as per Chidiac and Ross https://doi.org/10.1074/jbc.274.28.19639
# Also see Sprang https://dx.doi.org/10.1002%2Fbip.22836
Parameter("k_gtp_to_gdp_rgs", k_gtp_to_gdp_auto.value * 100)
# 3. PLC binding enhanced conversion of GTP to GDP
Parameter("k_gtp_to_gdp_plc", k_gtp_to_gdp_rgs.value / 2)
# Free Gaq:GDP recombines with Gbg
Parameter("k_gaq_gdp_binds_gbg", K_CONVERT)
# Parameter('k_gaq_gdp_unbinds_gbg', KR_BIND)

# PLC binding Gaq
Parameter("kf_PLC_bind_Gaq", KF_BIND)
Parameter("kr_PLC_bind_Gaq", KR_BIND)
# Conversion of PIP2 to IP3
Parameter("kf_PLC_bind_PIP2", KF_BIND)
Parameter("kr_PLC_bind_PIP2", KR_BIND)
Parameter("kcat_PIP2_to_IP3", KCAT)
# Binding of IP3 to IP3R
Parameter("kf_IP3_bind_IP3R", K_IP3_BIND)
Parameter("kr_IP3_bind_IP3R", KR_BIND)
# Transport of Ca2+
#  ER -> cytosol:
Parameter("kf_erCa_bind_IP3R", K_CA_BIND)
Parameter("kr_erCa_bind_IP3R", KR_BIND)
# Effective IP3R channel permeability as per Lemon et al. 2003 https://doi.org/10.1016/S0022-5193(03)00079-1
# is 525 1/s
Parameter("kcat_tranport_erCa", 525)
#  cytosol -> ER:
# Parameter('kf_cytCa_bind_IP3R', K_CA_BIND)
# Parameter('kr_cytCa_bind_IP3R', KR_BIND)
# Parameter('kf_cytCa_bind_IP3R', KF_BIND/10)
# Parameter('kr_cytCa_bind_IP3R', KR_BIND*10)
# Parameter('kcat_tranport_cytCa', K_ION_CHANNEL)

# Cytosolic Ca2+ regulation
# cytosol to extracellular space
# From previous model fittings we get around 4 1/s.
Parameter("k_Ca_cyt_to_extra", 4)  # 1/s
# extracellular space to cytosol
# Assume that at equilibrium before any agonist is added the rate out of the
# cytosol equals the rate in. We'll enforce this by setting the extra to cyt
# rate constant with an Expression.
Expression("k_Ca_extra_to_cyt", k_Ca_cyt_to_extra * Ca_C_0 / Ca_extra_0)

# Depeletion/metabolism of IP3
# 1.25 1/s as in Lemon et al. 2003 https://doi.org/10.1016/S0022-5193(03)00079-1
Parameter("kdeg_ip3", 1.25)

# Rules
# =====

# Alias the TAT:PAR2 complexes
tat_PAR2_i = (
    TAT(b=1) ** EXTRACELLULAR % PAR2(state="I", bortho=1, bgaq=None) ** CELL_MEMB
)
tat_PAR2_a = (
    TAT(b=1) ** EXTRACELLULAR % PAR2(state="A", bortho=1, bgaq=None) ** CELL_MEMB
)

# PAR2 synthesis
synthesize(PAR2(state="I", bortho=None, bgaq=None) ** CELL_MEMB, k_par2_synthesis)

# PAR2 degradation
degrade(PAR2(state="I", bortho=None, bgaq=None) ** CELL_MEMB, k_PAR2_I_degradation)
degrade(tat_PAR2_i, k_PAR2_bound_degradation)
degrade(tat_PAR2_a, k_PAR2_bound_degradation)
degrade(
    PAR2(state="D", bortho=None, bgaq=None) ** CELL_MEMB, k_par2_denatured_degradation
)

# 2-step activation of PAR2 by 2AT agonist:
#    2AT + PAR2_I <---> TAT:PAR2_I
Rule(
    "tat_bind_PAR2",
    TAT(b=None) ** EXTRACELLULAR + PAR2(state="I", bortho=None, bgaq=None) ** CELL_MEMB
    | tat_PAR2_i,
    kf_PAR2_bind_TAT,
    kr_PAR2_bind_TAT,
)
#    TAT:PAR2_I <---> TAT:PAR2_A
Rule("tat_activate_PAR2", tat_PAR2_i | tat_PAR2_a, k_activate_PAR2, k_inactivate_PAR2)

# Gaq activation by activated-PAR2:
#    PAR2_A + Gaq_I <---> PAR2_A:Gaq_I ---> PAR2_A + Gaq_A
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
    kf_PAR2_bind_Gaq,
    kr_PAR2_bind_Gaq,
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
    tat_PAR2_a_Gaq_gtp_Gbg >> tat_PAR2_a_Gaq_gtp + Gbg(b=None) ** CELL_MEMB,
    k_gbg_release,
)
# Alias the complex Gaq:GTP
Gaq_gtp = Gaq(bpar=None, bgdp=3, bgbg=None) ** CELL_MEMB % GTP(b=3) ** CELL_MEMB
# Gaq unbinds from PAR2
Rule("release_gaq", tat_PAR2_a_Gaq_gtp >> Gaq_gtp + tat_PAR2_a, k_gaq_release)
# Alias the complex Gaq:GDP
Gaq_gdp = Gaq(bpar=None, bgdp=3, bgbg=None) ** CELL_MEMB % GDP(b=3) ** CELL_MEMB
# Gaq can (slowly) hydolyze GTP to GDP
Rule("gtp_hydrolosis_auto", Gaq_gtp >> Gaq_gdp, k_gtp_to_gdp_auto)
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
# The Inactivated Gaq (Gaq:GDP) can reassociate the Beta-Gamma subunits to
# reform the heterotrimer.
Rule(
    "heterotrimer_reassociation",
    Gaq_gdp + Gbg(b=None) ** CELL_MEMB >> Gaq_gdp_Gbg,
    k_gaq_gdp_binds_gbg,
)

# PLC activation by binding Gaq:
#    Gaq_A + PLC <---> Gaq_A:PLC
#   Reusing the Gbg binding slot for PLC
bind_complex(
    Gaq_gtp, "bgbg", PLC() ** CELL_MEMB, "bgaq", [kf_PLC_bind_Gaq, kr_PLC_bind_Gaq]
)
# Conversion of PIP2 to IP3
#    Gaq_A:PLC + PIP2 <---> Gaq_A:PLC:PIP2 ---> Gaq_A:PLC + IP3
Gaq_gtp_PLC = (
    Gaq(bpar=None, bgdp=3, bgbg=1) ** CELL_MEMB
    % GTP(b=3) ** CELL_MEMB
    % PLC(bgaq=1) ** CELL_MEMB
)
catalyze_complex(
    Gaq_gtp_PLC,
    "bpip2",
    PIP2() ** CELL_MEMB,
    "b",
    IP3(b=None) ** CYTOSOL,
    [kf_PLC_bind_PIP2, kr_PLC_bind_PIP2, kcat_PIP2_to_IP3],
)
# Enhanced hydrolosis of GTP when Gaq is bound to PLC
#   Gaq:GTP:PLC ---> Gaq:GDP + PLC
Rule(
    "gtp_hydrolosis_plc",
    Gaq_gtp_PLC >> Gaq_gdp + PLC(bgaq=None, bpip2=None) ** CYTOSOL,
    k_gtp_to_gdp_plc,
)

# Binding of IP3 to IP3R - IP3R is activated when all 4 subunits are bound
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
    + Ca(loc="E", b=None) ** ER_LUMEN
    | IP3R(b1=1, b2=2, b3=3, b4=4, bcaer=5, bcacyt=None) ** ER_MEMB
    % IP3(b=1) ** CYTOSOL
    % IP3(b=2) ** CYTOSOL
    % IP3(b=3) ** CYTOSOL
    % IP3(b=4) ** CYTOSOL
    % Ca(loc="E", b=5) ** ER_LUMEN,
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
    % Ca(loc="E", b=5) ** ER_LUMEN
    >> IP3R(b1=1, b2=2, b3=3, b4=4, bcaer=None, bcacyt=None) ** ER_MEMB
    % IP3(b=1) ** CYTOSOL
    % IP3(b=2) ** CYTOSOL
    % IP3(b=3) ** CYTOSOL
    % IP3(b=4) ** CYTOSOL
    + Ca(loc="E", b=None) ** CYTOSOL,
    kcat_tranport_erCa,
)

# Regulation of Cytosolic Ca2+ --
# Here, we assume that the decay in the FRET signal is primarily due to
# first-order excretion of excess Ca2+ from the released ER store into the
# extracellular space by cell membrane ion channels. We assume both the
# excretion from the cytosol to the extracellular space and the reverse
# influx from the extracellular space to the cytosol are first order reactions.
# cytosol to extracellular space
Rule(
    "Ca_cyt_to_extra",
    Ca(loc="E", b=None) ** CYTOSOL >> Ca(loc="E", b=None) ** EXTRACELLULAR,
    k_Ca_cyt_to_extra,
)
# extracellular space to cytosol.
Rule(
    "Ca_extra_to_cyt",
    Ca(loc="E", b=None) ** EXTRACELLULAR >> Ca(loc="E", b=None) ** CYTOSOL,
    k_Ca_extra_to_cyt,
)

# Metabolic consumption of IP3
degrade(IP3(b=None) ** CYTOSOL, kdeg_ip3)

# Observables
# ===========
# Total amounts of each monomer
Observable("totTAT", TAT())
Observable("totPAR2", PAR2())
Observable("totGaq", Gaq())
Observable("totGbg", Gbg())
Observable("totGDP", GDP())
Observable("totGTP", GTP())
Observable("totRGS", RGS())
Observable("totPLC", PLC())
Observable("totPIP2", PIP2())
Observable("totIP3", IP3())
Observable("totIP3R", IP3R())
Observable("totCa", Ca())
# Inactive PAR2
Observable("iPAR2", PAR2(state="I"))
# Active PAR2
Observable("aPAR2", PAR2(state="A"))
# Denatured PAR2
Observable("dPAR2", PAR2(state="D"))
# Ro
Expression("Ro", iPAR2 + aPAR2)
# RL
Observable("LR", tat_PAR2_i)
Expression("occupancy_ratio", LR / Ro)
Expression("active_ratio", aPAR2 / Ro)
Observable(
    "aGaq_i", Gaq(bpar=None, bgdp=3, bgbg=None) ** CELL_MEMB % GTP(b=3) ** CELL_MEMB
)
Observable("aGaq_ii", tat_PAR2_a_Gaq_gtp)
Observable("aGaq_iii", Gaq_gtp_RGS)
Observable("aGaq_iv", Gaq_gtp_PLC)
Expression("aGaq", aGaq_i + aGaq_ii + aGaq_iii + aGaq_iv)
Expression("active_G_ratio", aGaq / totGaq)
# Active IP3R (i.e., all 4 subunits bound by IP3)
Observable(
    "aIP3R",
    IP3R(b1=1, b2=2, b3=3, b4=4) ** ER_MEMB
    % IP3(b=1) ** CYTOSOL
    % IP3(b=2) ** CYTOSOL
    % IP3(b=3) ** CYTOSOL
    % IP3(b=4) ** CYTOSOL,
)
Expression("active_IP3R_ratio", aIP3R / totIP3R)
# Fully inactive IP3R (i.e., no IP3 bound)
Observable("iIP3R", IP3R(b1=None, b2=None, b3=None, b4=None))
# The Ca2+ in the ER Lumen
Observable("erCa", Ca(loc="E", b=None) ** ER_LUMEN)
# Ca2+ in the Cytosol
Observable("cytoCa", Ca(loc="E", b=None) ** CYTOSOL)
Expression("Ca_num_to_nM", 1 / (Vcell * nM_to_num_per_pL))
Expression("cytoCa_nM", (cytoCa + Ca_C_0) * Ca_num_to_nM)
# Get the FRET signal
# The maximum FRET ratio, deltaR/R, for TN-XXL is 2.3 at 39 microM Ca2+,
# the effective Kd for Ca2+ binding to TN-XXL FRET reporter is
#  Kd = 800 nM,and the Hill-Coefficient is 1.5, https://doi.org/10.1038/nmeth.1243
Parameter("Kd_cytCa_bind_TNXXL", 800e-3)  # microM
Parameter("Rmax", 2.3)
Parameter("HillCoeff_TNXXL", 1.5)
# Compute the FRET ratio change relative to zero (i.e., Rmin) using the Hill equation,
#    (R-Rmin)/Rmin = Rmax*[Ca2+]**h / (Kd + [Ca2+]**h) ,
#       where Rmax is maximum FRET ratio change at saturation, h is the
#       Hill Coefficient, and Kd is effective dissociation constant.
Expression("Ca_num_to_microM", 1 / (Vcell * microM_to_num_per_pL))
# FRET ratio change for baseline concentration relative to zero - dR/R = (Rb-Rmin)/Rmin
Expression(
    "Frc_base",
    Rmax
    * (Ca_C_0 * Ca_num_to_microM) ** HillCoeff_TNXXL
    / (Kd_cytCa_bind_TNXXL + (Ca_C_0 * Ca_num_to_microM) ** HillCoeff_TNXXL),
)
# FRET ratio change for current concentration relative to zero - dR/R = (Rc-Rmin)/Rmin
Expression(
    "Frc_curr",
    Rmax
    * ((cytoCa + Ca_C_0) * Ca_num_to_microM) ** HillCoeff_TNXXL
    / (Kd_cytCa_bind_TNXXL + ((cytoCa + Ca_C_0) * Ca_num_to_microM) ** HillCoeff_TNXXL),
)
# Exp. FRET ratio change which is relative to the baseline - dR/R = (Rc-Rb)/Rb
Expression("FRET", (Frc_curr - Frc_base) / (Frc_base + 1))
# print(Frc_base.get_value())

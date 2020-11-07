'''PARM: PAR2 Activation-driven calcium Release Model

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
https://doi.org/10.3389/fnagi.2019.00089). The model also assumes that all four
subunits of the IP3 receptor, IP3R, must be bound by IP3 before calcium can
bind the receptor and be translocated between the ER lumen and cytosol which is
consistent with work by Alzayady et al
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

  1. PAR2 activation by 2AT:
      2AT + PAR2_I <---> TAT:PAR2_I ---> TAT + PAR2_A
  2. Gaq activation by activated-PAR2:  | Note: G-proteins are not pre-assembled on PAR2.
      i) G protein heterotrimer binds activated PAR2:
         PAR2_A + Gaq:GDP:Gbg <---> PAR2_A:Gaq:GDP:Gbg
     ii) GDP unbinds from Gaq:
         PAR2_A:Gaq:GDP:Gbc ---> PAR2_A:Gaq:Gbc + GDP
    iii) GTP binds Gaq:
         PAR2_A:Gaq:Gbc + GTP ---> PAR2_A:Gaq_A:GTP:Gbg
     iv) Gbg dissociates from Gaq (i.e., heterotrimer dissociation):
         PAR2_A:Gaq:GTP:Gbc ---> PAR2_A:Gaq:GTP + Gbc
      v) Gaq:GTP dissociates from PAR2 (G protein dissociation from the receptor):
         PAR2_A:Gaq:GTP ---> PAR2_A + Gaq:GTP
  3. Hydrolosis of GTP by Gaq
       a) Slow hydrolosis by Gaq alone
           Gaq:GTP ---> Gaq:GDP
       b) RGS enhanced hydrolosis
           Gaq:GTP + RGS <---> Gaq:GTP:RGS ---> Gaq:GDP + RGS
  4. Recombination of G protein heterotrimer
    Gaq:GDP + Gbg ---> Gaq:GDP:Gbq
  5. PLC activation by binding Gaq:
      Gaq_A:GTP + PLC <---> Gaq_A:GTP:PLC
  6. Conversion of PIP2 to IP3
      Gaq_A:GTP:PLC + PIP2 <---> Gaq_A:GTP:PLC:PIP2 ---> Gaq_A:PLC + IP3
  7. Binding of IP3 to IP3R:  | IP3R is only activated when all 4 subunits are bound.
     i) IP3R + IP3 <---> IP3R:IP3, subunit 1
    ii) IP3R + IP3 <---> IP3R:IP3, subunit 2
   iii) IP3R + IP3 <---> IP3R:IP3, subunit 3
    iv) IP3R + IP3 <---> IP3R:IP3, subunit 4
  8. Transport of Ca2+ by activated IP3R
    i) ER to cytosol:
       IP3R:IP3_4 + Ca_E <---> Ca_E:IP3R:IP3_4 ---> Ca_C + IP3R:IP3_4
   ii) Reverse, cytosol to ER:  | Assuming the transport is not just one way.
       IP3R:IP3_4 + Ca_C <---> Ca_C:IP3R:IP3_4 ---> Ca_E + IP3R:IP3_4
  9. Degradation of Cytosolic Calcium
       Ca_C ---> None,  if Ca_C > Ca_C_0
  10. Degradation of IP3
       IP3 ---> None


Unless otherwise noted the units used are:
    Volume : pL
    Area: micrometer^2
    Distance: micrometer
    Volume Concentration: number ( convert concentration by factor of V*N_A)
    Area Concentration: number
    Forward bimolecular association rate constants (kf) : 1/(s*number)
    Reverse of bimolecular association (dissociation) rate constants (kr) : 1/s
    Catalytic rate constants (kcat) :  1/s
    Other unimolecular rate constants: 1/s
    Dissociation (bimolecular) constants: Kd : number
    Binding (bimolecular) constants : Kb : 1/number
'''

# PySB components
from pysb import Model, Monomer, Parameter, Initial, Rule, Observable, Expression, Annotation, Compartment, ANY
# PySB macros
from pysb.macros import bind, bind_complex, catalyze, catalyze_complex, catalyze_state, degrade
# NumPy
import numpy as np
# Avogadro's Number from scipy
from scipy.constants import N_A

Model()

# Cellular volume, 10^-12 L as assumed in
# Albeck et al. https://doi.org/10.1371/journal.pbio.0060299
Parameter('Vcell', 1)
# Cell-membrane surface area
# Hek cells ~706.5 micrometer^2 https://doi.org/10.1038/s41598-017-07813-5
Parameter("SAcell", 706.5)

# Volume of the extracellular space
# We'll just assume twice the cell volume
Parameter("Vextra", 2*Vcell.value)

# Volume of the ER lumen/cisternal space.
# It is often >10% of cell volume according Alberts et al. https://www.ncbi.nlm.nih.gov/books/NBK26841/ .
# but for simplicity we will assume it is 15% of cell volume.
Parameter("Ver", Vcell.value * 0.15) # L
# Assume 10x the cell membrane surface area
Parameter("SAer", SAcell.value*10)

# Conversion factors for concentration units.
# microMolar to number/pL
microM_to_num_per_pL = 1e-6*N_A*1e-12
# nanoMolar to number/pL
nM_to_num_per_pL = 1e-9*N_A*1e-12


# Default forward, reverse, and catalytic rates:
KF_BIND = 1e-6 # 1/(molec*s) Default forward binding rate from Aldridge et al. https://doi.org/10.1038/ncb1497
KR_BIND = 1e-3 # Default dissociation rate from Albeck et al. https://doi.org/10.1371/journal.pbio.0060299
KCAT = 10 # "average enzyme" from Bar-Even et al. https://doi.org/10.1021/bi2002289

# Default signaling protein concentration range
# of 1 nM to 1 microM range assumed by Albeck et al. https://doi.org/10.1371/journal.pbio.0060299
# for which they reference Wu and Pollard https://doi.org/10.1126/science.1113230
# Also see Aldridge et al. https://doi.org/10.1038/ncb1497
# Set to 100 nM or 0.1 microM
SPC = 0.1*microM_to_num_per_pL

# As a first estimate for the forward binding reactions of Ca2+ we will assume
# that it is diffusion-controlled following Smoluchowski eqn.:
#     kf = 4*pi*D*R_o
# We will assume R_o is 2 nm and that D = D_Ca2+.
# The diffusion coefficient of Ca2+ is 5.3 x 10^-6 cm^2/s; https://doi.org/10.1016/0143-4160(87)90027-3
D_Ca = 5.3e-6 # cm^2/s
R_o = 2e-7 # cm
# (1e-3) term is for unit conversion from cm^3/s*number to 1/s*(number/L)
# mL/s*molec -> 10^-3 L/(s*number) and dividing by Ver converts to 1/(s*number)
K_CA_BIND = 4*np.pi*D_Ca*R_o*(1e-3)/(Ver.value*1e-12)

# Ion channel transport rate: up to 1e8 ions/s https://www.ncbi.nlm.nih.gov/books/NBK26910/
K_ION_CHANNEL = 1e8

# IP3 diffuses in mammalian at <= 10 micrometer^2/s https://dx.doi.org/10.1126%2Fscisignal.aag1625
D_ip3 = 10e-8 # cm^2/s
K_IP3_BIND = 4*np.pi*D_ip3*R_o*(1e-3)/(Vcell.value*1e-12)

# Default molecule degradation rate.
K_DEGRADE = 1 # 1/s

# Default unidirection conversion rate.
K_CONVERT = 1 # 1/s

# Compartments
# ============
# Since we are converting concentrations to numbers we can just leave
# the compartment sizes as 1, which is the default.
Compartment('EXTRACELLULAR', dimension=3)
# Cell Membrane
Compartment('CELL_MEMB', dimension=2, parent=EXTRACELLULAR)
# Cytosol
Compartment('CYTOSOL', dimension=3, parent=CELL_MEMB)
#  ER membrane
Compartment('ER_MEMB', dimension=2, parent=CYTOSOL)
# ER lumen volume
Compartment('ER_LUMEN', dimension=3, parent=ER_MEMB)

# Monomers
# ========
# PAR2 agonist 2AT
# Note: this was modeled as Trypsin ('Tryp') in an earlier version of the model
# which is a canonical protease activator of PAR2 (via N-terminal cleavage),
# but in the experiments of Kang et al. PAR2 is actually activated by the
# agonist 2AT, Kang et al. https://doi.org/10.1021/acsnano.9b01993)
Monomer('TAT', ['b'])
# PAR2, states: I = inactive (unbound), A = active (bound)
Monomer('PAR2', ['btat', 'bgaq','state'], {'state': ['I','A']})
# G-alpha_q G-protein unit, states: I = inactive, A = active
Monomer('Gaq', ['bpar','bgbg','bgdp'])
# G-beta-gamma G-protein units
Monomer('Gbg', ['b'])
# GDP
Monomer('GDP', ['b'])
# GTP
Monomer('GTP', ['b'])
# Regulator of G protein Signaling (RGS)
Monomer('RGS', ['b'])

# Phospholipase C
Monomer('PLC', ['bgaq','bpip2'])
# PIP2
Monomer('PIP2', ['b'])
# IP3
Monomer('IP3', ['b'])
# IP3 receptor
Monomer('IP3R', ['b1', 'b2', 'b3', 'b4', 'bcaer', 'bcacyt'])
# Calcium 2+, loc: E = ER space, C = cytosol
Monomer('Ca',['b', 'loc'],{'loc': ['E', 'C']})
# FRET reporter TN-XXL
Monomer('TNXXL', ['bca'])

# Annotations
# ===========
Annotation(PAR2, 'https://identifiers.org/uniprot:P55085')
Annotation(Gaq, 'https://identifiers.org/uniprot:P50148')
Annotation(PLC, 'https://identifiers.org/uniprot:Q9NQ66')
Annotation(IP3R, 'https://identifiers.org/uniprot:Q14643')

# Initial conditions
# ==================
# PAR2 agonist 2AT, 330 nM for Fig 2D data from Kang et al. https://doi.org/10.1021/acsnano.9b01993
C_2AT = 330 # nM
V_2AT = 50e6 # Volume of agonist added to wells is 50 microL
Vwell = 150e6 # Looks like the total well volume was 150 microL (100 microL ACSF + 50 microL agonist in ACSF)
nM_2AT_to_num = 1e-9 * N_A * 1e-12 * (V_2AT / Vwell) * Vextra.value
#nM_2AT_to_molec = 1e-9 * V_2AT * N_A
Parameter('TAT_0', C_2AT*nM_2AT_to_num)
Initial(TAT(b=None)**EXTRACELLULAR, TAT_0)
# inactive PAR2
# Endogenous receptor density of 1/micrometer^2 as in
# Falkenburger et al. 2010 https://dx.doi.org/10.1085%2Fjgp.200910344
Parameter('PAR2_0', 1*SAcell.value)
Initial(PAR2(state='I', btat=None,bgaq=None)**CELL_MEMB, PAR2_0)
# inactive G-protein heterotrimer Gaq-GDP:Gbg (the beta and gamma units are modeled as a single unit)
# G-protein density of 40/micrometer^2 as in
# Falkenburger et al. 2010 https://dx.doi.org/10.1085%2Fjgp.200910344
Parameter('Gaq_0', 40*SAcell.value)
# Alias the free Gprotein heterotrimer
Gaq_gdp_Gbg = Gaq(bpar=None, bgbg=3, bgdp=4)**CELL_MEMB % GDP(b=3)**CELL_MEMB % Gbg(b=4)**CELL_MEMB
Initial(Gaq_gdp_Gbg, Gaq_0)
# GTP
# Physiolocal concentration of GTP in mammalian cells is 468 +/- 224 microM
# as per Traut https://doi.org/10.1007/bf00928361
Parameter('GTP_0', 468*microM_to_num_per_pL*Vcell.value)
Initial(GTP(b=None)**CYTOSOL, GTP_0)
# RGS
# For RAW 264.7 Cell model between 0.008 and 0.012 microM as
# per Maurya and Subramaniam 2007 https://doi.org/10.1529/biophysj.106.097469
# Assume 0.010 microM is a reasonable starting point.
Parameter('RGS_0', 0.010*microM_to_num_per_pL*Vcell.value)
Initial(RGS(b=None)**CYTOSOL, RGS_0)

# inactive PLC
# Endogenous PLCB1 concentration of 3/micrometer^2 as in
# Falkenburger et al. 2010 https://dx.doi.org/10.1085%2Fjgp.200910344
Parameter('PLC_0', 3*SAcell.value)
Initial(PLC(bgaq=None, bpip2=None)**CELL_MEMB, PLC_0)
# PIP2
# Basal no. of PIP2 molecules is 49997 as per Lemon et al. 2003 https://doi.org/10.1016/S0022-5193(03)00079-1
# also free PIP2 of 5000 per micrometer^2 used by Falkenburger et al. 2013 https://doi.org/10.1085/jgp.201210887
# For nominal value will start with Lemon et al. value.
Parameter('PIP2_0', 49997) # Had to convert to area concentration.
Initial(PIP2(b=None)**CELL_MEMB, PIP2_0)
# IP3R
Parameter('IP3R_0', SPC*Ver.value)
Initial(IP3R(b1=None, b2=None, b3=None, b4=None, bcaer=None, bcacyt=None)**ER_MEMB, IP3R_0)
# ER Ca2+ store
# ER lumen of HEK-293 cells has between roughly 400-600 microM with an average
# around 525 microM as reported in
# Foyouzi-Youssefi et al. https://doi.org/10.1073/pnas.97.11.5723 (Fig. 3C, control)
Parameter('Ca_0', 525*microM_to_num_per_pL*Ver.value)
Initial(Ca(loc='E', b=None)**ER_LUMEN, Ca_0)
Parameter('Ca_C_0', 100*nM_to_num_per_pL*Vcell.value)
#Initial(Ca(loc='E', b=None)**CYTOSOL, Ca_C_0)

# Kinetic Parameters
# ==================
# PAR2 activation by 2AT
# Ca2+ signal Max. FRET Dose-Response for 2AT activation of PAR2
# has EC50 = 101.7 +- 28.7 nM, Kang et al. https://doi.org/10.1021/acsnano.9b01993
Parameter('kf_PAR2_bind_TAT', KF_BIND)
Parameter('kr_PAR2_bind_TAT', KR_BIND)
Parameter('kcat_activate_PAR2', KCAT)
# Gaq binding activated-PAR2
Parameter('kf_PAR2_bind_Gaq', KF_BIND)
Parameter('kr_PAR2_bind_Gaq', KR_BIND)
# Gaq release GDP
Parameter('k_gdp_release', K_CONVERT)
# Gaq bind GTP
Parameter('k_gtp_bind', K_CONVERT)
# Gbg dissociates from Gaq
Parameter('k_gbg_release', K_CONVERT)
# Gaq:GTP dissociates from PAR2
Parameter('k_gaq_release', K_CONVERT)
# Hydrolosis of GTP bound to Gaq
# 1. Autocatalysis rate for Gaq is ~0.8 1/min = 0.0133 1/s
# Bernstein et al. https://doi.org/10.1016/0092-8674(92)90165-9
# Also see Sprang https://dx.doi.org/10.1002%2Fbip.22836
Parameter('k_gtp_to_gdp_auto', 1.33e-2)
# 2. RGS binding and enhanced conversion of GTP to GDP
Parameter('kf_rgs_bind_gaq', KF_BIND)
Parameter('kr_rgs_bind_gaq', KR_BIND)
# Enhanced catalysis rate is for Gaq catalysis by RGS4
# is 100x higher as per Chidiac and Ross https://doi.org/10.1074/jbc.274.28.19639
# Also see Sprang https://dx.doi.org/10.1002%2Fbip.22836
Parameter('k_gtp_to_gdp_rgs', k_gtp_to_gdp_auto.value*100)
# Free Gaq:GDP recombines with Gbg
Parameter('k_gaq_gdp_binds_gbg', K_CONVERT)

# PLC binding Gaq
Parameter('kf_PLC_bind_Gaq', KF_BIND)
Parameter('kr_PLC_bind_Gaq', KR_BIND)
# Conversion of PIP2 to IP3
Parameter('kf_PLC_bind_PIP2', KF_BIND)
Parameter('kr_PLC_bind_PIP2', KR_BIND)
Parameter('kcat_PIP2_to_IP3', KCAT)
# Binding of IP3 to IP3R
Parameter('kf_IP3_bind_IP3R', K_IP3_BIND)
Parameter('kr_IP3_bind_IP3R', KR_BIND)
# Transport of Ca2+
#  ER -> cytosol:
Parameter('kf_erCa_bind_IP3R', K_CA_BIND)
Parameter('kr_erCa_bind_IP3R', KR_BIND)
# Effective IP3R channel permeability as per Lemon et al. 2003 https://doi.org/10.1016/S0022-5193(03)00079-1
# is 525 1/s
Parameter('kcat_tranport_erCa', 525)
#  cytosol -> ER:
#Parameter('kf_cytCa_bind_IP3R', K_CA_BIND)
#Parameter('kr_cytCa_bind_IP3R', KR_BIND)
#Parameter('kf_cytCa_bind_IP3R', KF_BIND/10)
#Parameter('kr_cytCa_bind_IP3R', KR_BIND*10)
#Parameter('kcat_tranport_cytCa', K_ION_CHANNEL)

# Depletion of Cytosolic Ca2+
# Base rate
Parameter('kdeg_cytCa', K_DEGRADE) # 1/s

# Depeletion/metabolism of IP3
# 1.25 1/s as in Lemon et al. 2003 https://doi.org/10.1016/S0022-5193(03)00079-1
Parameter('kdeg_ip3', 1.25)

# Rules
# =====
# PAR2 activation by 2AT:
#    2AT + PAR2_I <---> TAT:PAR2_I
#    TAT:PAR2_I ---> TAT:PAR2_A
# Alias the TAT:PAR2 complexes
tat_PAR2_i = TAT(b=1)**EXTRACELLULAR % PAR2(state='I', btat=1, bgaq=None)**CELL_MEMB
tat_PAR2_a = TAT(b=1)**EXTRACELLULAR % PAR2(state='A', btat=1, bgaq=None)**CELL_MEMB
Rule('tat_bind_PAR2', TAT(b=None)**EXTRACELLULAR + PAR2(state='I', btat=None, bgaq=None)**CELL_MEMB
     | tat_PAR2_i, kf_PAR2_bind_TAT,kr_PAR2_bind_TAT)
Rule('tat_activate_PAR2', tat_PAR2_i >> tat_PAR2_a, kcat_activate_PAR2)
# Gaq activation by activated-PAR2:
#    PAR2_A + Gaq_I <---> PAR2_A:Gaq_I ---> PAR2_A + Gaq_A
tat_PAR2_a_Gaq_gdp_Gbg = (TAT(b=1)**EXTRACELLULAR %
                          PAR2(state='A', btat=1, bgaq=2)**CELL_MEMB %
                           Gaq(bpar=2, bgdp=3, bgbg=4)**CELL_MEMB %
                           GDP(b=3)**CELL_MEMB % Gbg(b=4)**CELL_MEMB)

Rule('par2_bind_gaq', tat_PAR2_a + Gaq_gdp_Gbg | tat_PAR2_a_Gaq_gdp_Gbg,
     kf_PAR2_bind_Gaq,kr_PAR2_bind_Gaq)
tat_PAR2_a_Gaq_Gbg = (TAT(b=1)**EXTRACELLULAR %
                          PAR2(state='A', btat=1, bgaq=2)**CELL_MEMB %
                           Gaq(bpar=2, bgdp=None, bgbg=4)**CELL_MEMB %
                           Gbg(b=4)**CELL_MEMB)
Rule('gaq_releases_gdp',tat_PAR2_a_Gaq_gdp_Gbg >> tat_PAR2_a_Gaq_Gbg +
     GDP(b=None)**CYTOSOL, k_gdp_release)
tat_PAR2_a_Gaq_gtp_Gbg = (TAT(b=1)**EXTRACELLULAR %
                          PAR2(state='A', btat=1, bgaq=2)**CELL_MEMB %
                           Gaq(bpar=2, bgdp=3, bgbg=4)**CELL_MEMB %
                           GTP(b=3)**CELL_MEMB % Gbg(b=4)**CELL_MEMB)
Rule('gaq_binds_gtp', tat_PAR2_a_Gaq_Gbg + GTP(b=None)**CYTOSOL >>
    tat_PAR2_a_Gaq_gtp_Gbg, k_gtp_bind)
tat_PAR2_a_Gaq_gtp = (TAT(b=1)**EXTRACELLULAR %
                          PAR2(state='A', btat=1, bgaq=2)**CELL_MEMB %
                           Gaq(bpar=2, bgdp=3, bgbg=None)**CELL_MEMB %
                           GTP(b=3)**CELL_MEMB)
Rule('release_gbg', tat_PAR2_a_Gaq_gtp_Gbg >> tat_PAR2_a_Gaq_gtp + Gbg(b=None)**CELL_MEMB, k_gbg_release)
Gaq_gtp = (Gaq(bpar=None, bgdp=3, bgbg=None)**CELL_MEMB % GTP(b=3)**CELL_MEMB)
Rule('release_gaq', tat_PAR2_a_Gaq_gtp >> Gaq_gtp, k_gaq_release)
Gaq_gdp = (Gaq(bpar=None, bgdp=3, bgbg=None)**CELL_MEMB % GDP(b=3)**CELL_MEMB)
Rule('gtp_hydrolosis_auto', Gaq_gtp >> Gaq_gdp, k_gtp_to_gdp_auto)
Gaq_gtp_RGS = (Gaq(bpar=None, bgdp=3, bgbg=1)**CELL_MEMB % GTP(b=3)**CELL_MEMB
               % RGS(b=1)**CYTOSOL)
Rule('gaq_gtp_binds_rgs', Gaq_gtp + RGS(b=None)**CYTOSOL | Gaq_gtp_RGS, kf_rgs_bind_gaq, kr_rgs_bind_gaq)
Rule('gtp_hydrolosis_rgs', Gaq_gtp_RGS >> Gaq_gdp + RGS(b=None)**CYTOSOL, k_gtp_to_gdp_rgs)

Rule('heterotrimer_reassociation', Gaq_gdp + Gbg(b=None)**CELL_MEMB >> Gaq_gdp_Gbg, k_gaq_gdp_binds_gbg)


# PLC activation by binding Gaq:
#    Gaq_A + PLC <---> Gaq_A:PLC
#   Reusing the Gbg binding slot for PLC
bind_complex(Gaq_gtp, 'bgbg', PLC()**CELL_MEMB, 'bgaq', [kf_PLC_bind_Gaq,kr_PLC_bind_Gaq])
# Conversion of PIP2 to IP3
#    Gaq_A:PLC + PIP2 <---> Gaq_A:PLC:PIP2 ---> Gaq_A:PLC + IP3
Gaq_gtp_PLC = (Gaq(bpar=None, bgdp=3, bgbg=1)**CELL_MEMB % GTP(b=3)**CELL_MEMB
               % PLC(bgaq=1)**CELL_MEMB)
catalyze_complex(Gaq_gtp_PLC, 'bpip2', PIP2()**CELL_MEMB, 'b', IP3(b=None)**CYTOSOL,
                 [kf_PLC_bind_PIP2,kr_PLC_bind_PIP2,kcat_PIP2_to_IP3])
# Binding of IP3 to IP3R - IP3R is activated when all 4 subunits are bound
#   IP3R + IP3 <---> IP3R:IP3, subunit 1
bind(IP3R(b2=None,b3=None,b4=None,bcaer=None,bcacyt=None)**ER_MEMB, 'b1', IP3(b=None)**CYTOSOL, 'b', [kf_IP3_bind_IP3R,kr_IP3_bind_IP3R])
#   IP3R + IP3 <---> IP3R:IP3, subunit 2
#bind_complex(IP3R(b1=1,b2=None,b3=None,b4=None,bcaer=None,bcacyt=None)**ER_MEMB % IP3(b=1)**CYTOSOL, 'b2', IP3(b=None)**CYTOSOL, 'b',
#             [kf_IP3_bind_IP3R,kr_IP3_bind_IP3R])
Rule('bind_IP3_IPR3_sub2', IP3R(b1=1,b2=None,b3=None,b4=None,bcaer=None,bcacyt=None)**ER_MEMB
     % IP3(b=1)**CYTOSOL + IP3(b=None)**CYTOSOL |
     IP3R(b1=1,b2=2,b3=None,b4=None,bcaer=None,bcacyt=None)**ER_MEMB %
     IP3(b=1)**CYTOSOL % IP3(b=2)**CYTOSOL, kf_IP3_bind_IP3R, kr_IP3_bind_IP3R)
#   IP3R + IP3 <---> IP3R:IP3, subunit 3
#bind_complex(IP3R(b1=1, b2=50, b3=None,b4=None,bcaer=None,bcacyt=None)**ER_MEMB % IP3(b=1)**CYTOSOL % IP3(b=50)**CYTOSOL, 'b3', IP3(b=None)**CYTOSOL, 'b',
#            [kf_IP3_bind_IP3R,kr_IP3_bind_IP3R])
Rule('bind_IP3_IPR3_sub3',
     IP3R(b1=1,b2=2,b3=None,b4=None,bcaer=None,bcacyt=None)**ER_MEMB %
     IP3(b=1)**CYTOSOL % IP3(b=2)**CYTOSOL + IP3(b=None)**CYTOSOL |
     IP3R(b1=1,b2=2,b3=3,b4=None,bcaer=None,bcacyt=None)**ER_MEMB %
     IP3(b=1)**CYTOSOL % IP3(b=2)**CYTOSOL % IP3(b=3)**CYTOSOL,
     kf_IP3_bind_IP3R, kr_IP3_bind_IP3R)
#   IP3R + IP3 <---> IP3R:IP3, subunit 4
#bind_complex(IP3R(b1=1, b2=50,b3=50,b4=None)**ER_MEMB % IP3(b=1)**CYTOSOL % IP3(b=50)**CYTOSOL % IP3(b=50)**CYTOSOL, 'b4',
#             IP3(b=None)**CYTOSOL, 'b', [kf_IP3_bind_IP3R,kr_IP3_bind_IP3R])
Rule('bind_IP3_IPR3_sub4',
     IP3R(b1=1,b2=2,b3=3,b4=None,bcaer=None,bcacyt=None)**ER_MEMB %
     IP3(b=1)**CYTOSOL % IP3(b=2)**CYTOSOL % IP3(b=3)**CYTOSOL +
     IP3(b=None)**CYTOSOL |
     IP3R(b1=1,b2=2,b3=3,b4=4,bcaer=None,bcacyt=None)**ER_MEMB %
     IP3(b=1)**CYTOSOL % IP3(b=2)**CYTOSOL % IP3(b=3)**CYTOSOL
     % IP3(b=4)**CYTOSOL,
     kf_IP3_bind_IP3R, kr_IP3_bind_IP3R)
# Transport of Ca2+ by activated IP3R
#  ER -> cytosol:
#    IP3R:IP3_4 + Ca_E <---> Ca_E:IP3R:IP3_4 ---> Ca_C + IP3R:IP3_4
#catalyze_complex(IP3R(b1=1, b2=2, b3=3, b4=4,bcaer=None,bcacyt=None)**ER_MEMB % IP3(b=1)**CYTOSOL % IP3(b=2)**CYTOSOL
#                  % IP3(b=3)**CYTOSOL % IP3(b=4)**CYTOSOL, 'bcaer', Ca(loc='E', b=None)**ER_LUMEN, 'b', Ca(loc='C', b=None)**CYTOSOL,
#                  [kf_erCa_bind_IP3R, kr_erCa_bind_IP3R, kcat_tranport_erCa])
Rule('bind_Ca_IPR3_er',
     IP3R(b1=1,b2=2,b3=3,b4=4,bcaer=None,bcacyt=None)**ER_MEMB %
     IP3(b=1)**CYTOSOL % IP3(b=2)**CYTOSOL % IP3(b=3)**CYTOSOL %
     IP3(b=4)**CYTOSOL + Ca(loc='E', b=None)**ER_LUMEN |
     IP3R(b1=1,b2=2,b3=3,b4=4,bcaer=5,bcacyt=None)**ER_MEMB %
     IP3(b=1)**CYTOSOL % IP3(b=2)**CYTOSOL % IP3(b=3)**CYTOSOL
     % IP3(b=4)**CYTOSOL % Ca(loc='E',b=5)**ER_LUMEN,
     kf_erCa_bind_IP3R, kr_erCa_bind_IP3R)
Rule('transport_Ca_ER_CYTO',
     IP3R(b1=1,b2=2,b3=3,b4=4,bcaer=5,bcacyt=None)**ER_MEMB %
     IP3(b=1)**CYTOSOL % IP3(b=2)**CYTOSOL % IP3(b=3)**CYTOSOL
     % IP3(b=4)**CYTOSOL % Ca(loc='E',b=5)**ER_LUMEN >>
     IP3R(b1=1,b2=2,b3=3,b4=4,bcaer=None,bcacyt=None)**ER_MEMB %
     IP3(b=1)**CYTOSOL % IP3(b=2)**CYTOSOL % IP3(b=3)**CYTOSOL %
     IP3(b=4)**CYTOSOL + Ca(loc='E', b=None)**CYTOSOL, kcat_tranport_erCa)
#Rule('transport_Ca_ER_CYTO', Ca(loc='E', b=None)**ER_LUMEN >> Ca(loc='E', b=None)**CYTOSOL, kcat_tranport_erCa)

#  Reverse, cytosol -> ER:
#    IP3R:IP3_4 + Ca_C <---> Ca_C:IP3R:IP3_4 ---> Ca_E + IP3R:IP3_4
#catalyze_complex(IP3R(b1=1, b2=2, b3=3, b4=4,bcaer=None,bcacyt=None)**ER_MEMB % IP3(b=1)**CYTOSOL % IP3(b=2)**CYTOSOL
#                  % IP3(b=3)**CYTOSOL % IP3(b=4)**CYTOSOL, 'bcacyt', Ca(loc='C', b=None)**CYTOSOL, 'b', Ca(loc='E', b=None)**ER_LUMEN,
#                  [kf_cytCa_bind_IP3R, kr_cytCa_bind_IP3R, kcat_tranport_cytCa])
# Rule('bind_Ca_IPR3_cyto',
#      IP3R(b1=1,b2=2,b3=3,b4=4,bcaer=None,bcacyt=None)**ER_MEMB %
#      IP3(b=1)**CYTOSOL % IP3(b=2)**CYTOSOL % IP3(b=3)**CYTOSOL %
#      IP3(b=4)**CYTOSOL + Ca(loc='E', b=None)**CYTOSOL |
#      IP3R(b1=1,b2=2,b3=3,b4=4,bcaer=None,bcacyt=5)**ER_MEMB %
#      IP3(b=1)**CYTOSOL % IP3(b=2)**CYTOSOL % IP3(b=3)**CYTOSOL
#      % IP3(b=4)**CYTOSOL % Ca(loc='E',b=5)**CYTOSOL,
#      kf_erCa_bind_IP3R, kr_erCa_bind_IP3R)
# Rule('transport_Ca_CYTO_ER',
#      IP3R(b1=1,b2=2,b3=3,b4=4,bcaer=None,bcacyt=5)**ER_MEMB %
#      IP3(b=1)**CYTOSOL % IP3(b=2)**CYTOSOL % IP3(b=3)**CYTOSOL
#      % IP3(b=4)**CYTOSOL % Ca(loc='E',b=5)**CYTOSOL >>
#      IP3R(b1=1,b2=2,b3=3,b4=4,bcaer=None,bcacyt=None)**ER_MEMB %
#      IP3(b=1)**CYTOSOL % IP3(b=2)**CYTOSOL % IP3(b=3)**CYTOSOL %
#      IP3(b=4)**CYTOSOL + Ca(loc='E', b=None)**ER_LUMEN, kcat_tranport_cytCa)



# Degradation of Cytosolic Ca2+ --
# This term was added to help fit the decay of FRET signal, presumably
# representing a lumped process for the regulation of Ca2+ concentration in the
# cytosol after the ER store is released (e.g., activation of
# SERCA to pump Ca2+ back into the lumen, or activation of cell membrane ion
# channels to release excess Ca2+ into the extracellular space).
degrade(Ca(loc='E', b=None)**CYTOSOL, kdeg_cytCa)

degrade(IP3(b=None)**CYTOSOL, kdeg_ip3)

# Observables
# ===========
Observable('iPAR2', PAR2(state='I'))
Observable('aPAR2', PAR2(state='A'))
Observable('aIP3R', IP3R(b1=1, b2=2, b3=3, b4=4)**ER_MEMB % IP3(b=1)**CYTOSOL % IP3(b=2)**CYTOSOL % IP3(b=3)**CYTOSOL % IP3(b=4)**CYTOSOL)
Observable('iIP3R', IP3R(b1=None,b2=None,b3=None,b4=None))
Observable('erCa', Ca(loc='E', b=None)**ER_LUMEN)
Observable('cytoCa', Ca(loc='E', b=None)**CYTOSOL)
# Get the FRET signal
# The maximum FRET ratio, deltaR/R, for TN-XXL is 2.3 at 39 microM Ca2+,
# the effective Kd for Ca2+ binding to TN-XXL FRET reporter is
#  Kd = 800 nM,and the Hill-Coefficient is 1.5, https://doi.org/10.1038/nmeth.1243
Parameter('Kd_cytCa_bind_TNXXL', 800e-3) # microM
Parameter('Rmax', 2.3)
Parameter('HillCoeff_TNXXL', 1.5)
# Compute the FRET ratio change relative to zero (i.e., Rmin) using the Hill equation,
#    (R-Rmin)/Rmin = Rmax*[Ca2+]**h / (Kd + [Ca2+]**h) ,
#       where Rmax is maximum FRET ratio change at saturation, h is the
#       Hill Coefficient, and Kd is effective dissociation constant.
Expression('Ca_num_to_microM', 1/(Vcell*microM_to_num_per_pL))
# FRET ratio change for baseline concentration relative to zero - dR/R = (Rb-Rmin)/Rmin
Expression('Frc_base', Rmax*(Ca_C_0*Ca_num_to_microM)**HillCoeff_TNXXL / (Kd_cytCa_bind_TNXXL + (Ca_C_0*Ca_num_to_microM)**HillCoeff_TNXXL))
# FRET ratio change for current concentration relative to zero - dR/R = (Rc-Rmin)/Rmin
Expression('Frc_curr', Rmax*((cytoCa+Ca_C_0)*Ca_num_to_microM)**HillCoeff_TNXXL / (Kd_cytCa_bind_TNXXL + ((cytoCa+Ca_C_0)*Ca_num_to_microM)**HillCoeff_TNXXL))
# Exp. FRET ratio change which is relative to the baseline - dR/R = (Rc-Rb)/Rb
Expression('FRET', (Frc_curr - Frc_base)/(Frc_base + 1))
print(Frc_base.get_value())

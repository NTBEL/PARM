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
 https://doi.org/10.3389/fnagi.2019.00089). In this model the initial cytosolic
 calcium is set to zero, so only calcium released from the ER is registered by
 the FRET calcium sensor TN-XXL. The interaction between Ca2+ and TN-XXL is
 explicitly included as a binding reaction with the effective Kd from the
 dose-repsonse curve of TN-XXL to Ca2+. The FRET ratio for cytosolic Ca2+ is
 therefore estimated as the fraction of occupied TN-XXL molecules multiplied by
 the maximum FRET ratio. The model also assumes that all four subunits of the
 IP3 receptor, IP3R, must be bound by IP3 before calcium can bind the receptor
 and be translocated between the ER lumen and cytosol which is consistent with
 work by Alzayady et al https://doi.org/10.1126/scisignal.aad6281. However,
 feedback to either reduce or enhance the IP3R calcium release is not included.
 Additionally, the cytosolic calcium level maintenance is modeled by a
 unidirectional degradation step, approximating a lumped process for the
 regulation of Ca2+ concentration in the cytosol after the ER store is released
 (e.g., activation of SERCA to pump Ca2+ back into the lumen, or activation of
 cell membrane ion channels to release excess Ca2+ into the extracellular
 space).

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
  9. Binding of Calcium to the TN-XXL FRET reporter:
       TNXXL + Ca <---> TNXXL:Ca
  10. Degradation of Cytosolic Calcium
       Ca_C ---> None,  if Ca_C > Ca_C_0
  11. Degradation of IP3
       IP3 ---> None


Unless otherwise noted the units used are:
    Volume : L
    Area: L^(2/3)
    Distance: L^(1/3)
    Volume Concentration: molec (M*N_A*V) - molec is short for molecules
    Area Concentration: molec
    Forward bimolecular association rate constants (kf) : 1/(s*(molec))
    Reverse of bimolecular association (dissociation) rate constants (kr) : 1/s
    Catalytic rate constants (kcat) :  1/s
    Other unimolecular rate constants: 1/s
    Dissociation (bimolecular) constants: Kd : molec
    Binding (bimolecular) constants : Kb : 1/(molec)
'''

from pysb import Model, Monomer, Parameter, Initial, Rule, Observable, Expression, Annotation, Compartment, ANY
from pysb.macros import bind, bind_complex, catalyze, catalyze_complex, catalyze_state, degrade
import numpy as np

Model()

# Cellular volume, 10^-12 L as assumed in
# Albeck et al. https://doi.org/10.1371/journal.pbio.0060299
Vcell = 1.e-12
# Volume of the ER lumen/cisternal space.
# It is often >10% of cell volume according Aleberts et al. https://www.ncbi.nlm.nih.gov/books/NBK26841/ .
# but for simplicity we will assume it is 10% of cell volume.
Ver = Vcell * 0.1

# Avogadro's number
N_A = 6.02214e23 # molec/mol

# Conversion factors for concentration units.
# microMolar to molec
microM_to_molec = 1e-6*N_A*Vcell
# nanoMolar to molec
nM_to_molec = 1e-9*N_A*Vcell

# Default forward, reverse, and catalytic rates:
KF_BIND = 1e-6 # 1/(molec*s) Default forward binding rate from Aldridge et al. https://doi.org/10.1038/ncb1497
KR_BIND = 1e-3 # Default dissociation rate from Albeck et al. https://doi.org/10.1371/journal.pbio.0060299
KCAT = 10 # "average enzyme" from Bar-Even et al. https://doi.org/10.1021/bi2002289

# Default signaling protein concentration, which is midway between
# the 1 nM to 1 microM range assumed by Albeck et al. https://doi.org/10.1371/journal.pbio.0060299
# for which they reference Wu and Pollard https://doi.org/10.1126/science.1113230
# Also see Aldridge et al. https://doi.org/10.1038/ncb1497
# Midpoint is 0.5 microM.
SPC = 0.5*microM_to_molec

# As a first estimate for the forward binding reactions of Ca2+ we will assume
# that it is diffusion-controlled following Smoluchowski eqn.:
#     kf = 4*pi*D*R_o
# We will assume R_o is 2 nm and that D = D_Ca2+.
# The diffusion coefficient of Ca2+ is 5.3 x 10^-6 cm^2/s; https://doi.org/10.1016/0143-4160(87)90027-3
D_Ca = 5.3e-6 # cm^2/s
R_o = 2e-7 # cm
# (1e-3) term is for unit conversion from cm^3/s*molec to 1/s*(molec/L)
# mL/s*molec -> 10^-3 L/s*molec and dividing by Vcell converts to 1/(s*molec)
K_CA_BIND = 4*np.pi*D_Ca*R_o*(1e-3)/Vcell

# Ion channel transport rate: up to 1e8 ions/s https://www.ncbi.nlm.nih.gov/books/NBK26910/
K_ION_CHANNEL = 1e8

# IP3 diffuses in mammalian at <= 10 micrometer^2/s https://dx.doi.org/10.1126%2Fscisignal.aag1625
D_ip3 = 10e-8 # cm^2/s
K_IP3_BIND = 4*np.pi*D_ip3*R_o*(1e-3)/Vcell

# Compartments
# ============
# Since we have already converted initial concentrations into molec units
# we need to scale the volumes used by compartments by Vcell.
# Assume extracellular volume is spherical with volume twice that of the
# cell volume.
Parameter('V_EXTRA', 2)
Compartment('EXTRACELLULAR', dimension=3, size=V_EXTRA)
# Cell Membrane
Parameter('SA_CM', 4*np.pi*(3/(4*np.pi))**1.5)
Compartment('CELL_MEMB', dimension=2, size=SA_CM, parent=EXTRACELLULAR)
# Cytosol
Parameter('V_C', 1)
Compartment('CYTOSOL', dimension=3, size=V_C, parent=CELL_MEMB)
#  ER membrane
Parameter('SA_ER', 4*np.pi*(3*Ver/(Vcell*4*np.pi))**1.5)
Compartment('ER_MEMB', dimension=2, size=SA_ER, parent=CYTOSOL)
#  internal ER lumen space
Parameter('V_ER', Ver/Vcell)
Compartment('ER_LUMEN', dimension=3, size=V_ER, parent=ER_MEMB)

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
V_2AT = 50e-6 # Volume of agonist added to wells is 50 microL
Vwell = 150e-6 # Looks like the total well volume was 150 microL (100 microL ACSF + 50 microL agonist in ACSF)
nM_2AT_to_molec = 1e-9 * V_2AT * N_A * (Vcell/Vwell)
#nM_2AT_to_molec = 1e-9 * V_2AT * N_A
Parameter('TAT_0', C_2AT*nM_2AT_to_molec)
Initial(TAT(b=None)**EXTRACELLULAR, TAT_0)
# inactive PAR2
Parameter('PAR2_0', SPC*V_C.value/SA_CM.value) # Had to convert to area concentration.
Initial(PAR2(state='I', btat=None,bgaq=None)**CELL_MEMB, PAR2_0)
# inactive G-protein heterotrimer Gaq-GDP:Gbg (the beta and gamma units are modeled as a single unit)
Parameter('Gaq_0', SPC*V_C.value/SA_CM.value)
# Alias the free Gprotein heterotrimer
Gaq_gdp_Gbg = Gaq(bpar=None, bgbg=3, bgdp=4)**CELL_MEMB % GDP(b=3)**CELL_MEMB % Gbg(b=4)**CELL_MEMB
Initial(Gaq_gdp_Gbg, Gaq_0)
# GTP
# Physiolocal concentration of GTP in mammalian cells is 468 +/- 224 microM
# as per Traut https://doi.org/10.1007/bf00928361
Parameter('GTP_0', 468*microM_to_molec)
Initial(GTP(b=None)**CYTOSOL, GTP_0)
# RGS
Parameter('RGS_0', SPC)
Initial(RGS(b=None)**CYTOSOL, RGS_0)

# inactive PLC
Parameter('PLC_0', SPC*V_C.value/SA_CM.value) # Had to convert to area concentration.
Initial(PLC(bgaq=None, bpip2=None)**CELL_MEMB, PLC_0)
# PIP2
Parameter('PIP2_0', SPC*V_C.value/SA_CM.value) # Had to convert to area concentration.
Initial(PIP2(b=None)**CELL_MEMB, PIP2_0)
# IP3R
Parameter('IP3R_0', SPC*V_ER.value/SA_ER.value) # Had to convert to area concentration.
Initial(IP3R(b1=None, b2=None, b3=None, b4=None, bcaer=None, bcacyt=None)**ER_MEMB, IP3R_0)
# ER Ca2+ store
# ER lumen of HEK-293 cells has between roughly 400-600 microM with an average
# around 525 microM as reported in
# Foyouzi-Youssefi et al. https://doi.org/10.1073/pnas.97.11.5723 (Fig. 3C, control)
Parameter('Ca_0', 525*microM_to_molec)
Initial(Ca(loc='E', b=None)**ER_LUMEN, Ca_0)
#Parameter('Ca_C_0', 100*nM_to_molec)
#Initial(Ca(loc='E', b=None)**CYTOSOL, Ca_C_0)
# TN-XXL
Parameter('TNXXL_0', 100e-3*microM_to_molec)
Initial(TNXXL(bca=None)**CYTOSOL, TNXXL_0)

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
Parameter('k_gdp_release', 1)
# Gaq bind GTP
Parameter('k_gtp_bind', 1)
# Gbg dissociates from Gaq
Parameter('k_gbg_release', 1)
# Gaq:GTP dissociates from PAR2
Parameter('k_gaq_release', 1)
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
Parameter('k_gaq_gdp_binds_gbg', 1)

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
#Parameter('kf_erCa_bind_IP3R', KF_BIND)
#Parameter('kr_erCa_bind_IP3R', KR_BIND)
Parameter('kcat_tranport_erCa', K_ION_CHANNEL)
#  cytosol -> ER:
Parameter('kf_cytCa_bind_IP3R', K_CA_BIND)
Parameter('kr_cytCa_bind_IP3R', KR_BIND)
#Parameter('kf_cytCa_bind_IP3R', KF_BIND/10)
#Parameter('kr_cytCa_bind_IP3R', KR_BIND*10)
Parameter('kcat_tranport_cytCa', K_ION_CHANNEL)
# Ca2+ binding to TN-XXL FRET reporter
#  Kd = 800 nM, https://doi.org/10.1038/nmeth.1243
Parameter('Kd_cytCa_bind_TNXXL', 800*nM_to_molec)
Parameter('kf_cytCa_bind_TNXXL', K_CA_BIND)
Expression('kr_cytCa_bind_TNXXL', kf_cytCa_bind_TNXXL*Kd_cytCa_bind_TNXXL)
#Parameter('kr_cytCa_bind_TNXXL', )
# Depletion of Cytosolic Ca2+
Parameter('kdeg_cytCa', 1) # 1/s
# Depeletion/metabolism of IP3
Parameter('kdeg_ip3', 1)
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


# Binding of Calcium to the TN-XXL FRET reporter
#  TNXXL + Ca <---> TNXXL:Ca
bind(TNXXL(bca=None)**CYTOSOL, 'bca', Ca(loc='E',b=None)**CYTOSOL, 'b',
     (kf_cytCa_bind_TNXXL,kr_cytCa_bind_TNXXL))


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
Observable('cytoIP3', IP3(b=ANY)**CYTOSOL)

# Get the FRET signal
Observable('fret_complex', TNXXL(bca=1)**CYTOSOL % Ca(b=1,loc='E')**CYTOSOL)
# The maximum FRET ratio, deltaR/R, for TN-XXL is 2.3 at 39 microM Ca2+,
# see Fig S1 C https://doi.org/10.1038/nmeth.1243
Parameter('max_FRET_ratio', 2.3)
Expression('FRET', max_FRET_ratio*fret_complex/TNXXL_0)
#Observable('FRET', fret_signal)

'''PARM: PAR2 Activation-driven calcium Release Model
Unless otherwise noted the units used are:
    Concentrations: microM
    Forward bimolecular association rate constants (kf) : 1/(s*microM)
    Reverse of bimolecular association (dissociation) rate constants (kr) : 1/s
    Catalytic rate constants (kcat) :  1/s
    Dissociation (bimolecular) constants: Kd : microM
    Binding (bimolecular) constants : Kb : 1/microM
    Volume : micron^3
    Distance: micron

Revelant equations:
    Relationship of equilibrium dissociation constant and bimolecular association
    rate constants:
        Kd = kr/kf
    Relationship of equilibrium binding constant and bimolecular association
    rate constants:
        Kb = kf/kr
    Relationship of equilibrium binding and dissociation constants for bimolecular
    association:
        Kb = 1/Kd
    Diffusion-controlled rate for bimolecular reaction
    (Smoluchowski's diffusion controlled encounter):
        kf = 4*pi*D*R_o/V,
            D = D_1 + D_2 is the apparent diffusion coefficient,
            R_o = R_1 + R_2 is the interaction radius,
            V is the total volume.
    Catalytic efficiency for Michaelis-Menten kinetics:
        keff = kcat/Km,
            Km is the Michaelis constant which is
                Km = (kcat + kr)/kf,
            and can be approximately equal to the ES dissociation constant (Kd)
            when kcat << kr.



'''

from pysb import Model, Monomer, Parameter, Initial, Rule, Observable, Expression, Annotation, Compartment
from pysb.macros import bind, bind_complex, catalyze, catalyze_complex
import numpy as np

Model()

# Cellular volume, value (10^-12 L) as assumed in
# Albeck et al. https://doi.org/10.1371/journal.pbio.0060299
Vcell = 10.
# Volume of the ER lumen/cisternal space.
# It is often >10% of cell volume according Aleberts et al. https://www.ncbi.nlm.nih.gov/books/NBK26841/ .
# but for simplicity we will assume it is 10% of cell volume.
Ver = Vcell * 0.1
# Default forward, reverse, and catalytic rates:
KF = 1e-6 #
KR = 1e-3 # Default dissociation rate from Albeck et al. https://doi.org/10.1371/journal.pbio.0060299
KCAT = 10 # "average enzyme" from Bar-Even et al. https://doi.org/10.1021/bi2002289

# Default signaling protein concentration, which is midway between
# the 1 nM to 1 microM range assumed by Albeck et al. https://doi.org/10.1371/journal.pbio.0060299
# for which they reference Wu and Pollard https://doi.org/10.1126/science.1113230
SPC = 0.5

# Compartments
# ============
# Cytosol
# Parameter('SA_CON', 4*np.pi*(6*Vcell/(4*np.pi))**1.5)
# Compartment('CONTAINER', dimension=2, size=SA_CM)
# Parameter('V_EXTRA', 2*Vcell)
# Compartment('EXTRACELLULAR', dimension=3, size=V_C)
# Parameter('SA_CM', 4*np.pi*(3*Vcell/(4*np.pi))**1.5)
# Compartment('CELL_MEMB', dimension=2, size=SA_CM)
# Parameter('V_C', Vcell)
# Compartment('CYTOSOL', dimension=3, size=V_C)
# #  ER membrane
# Parameter('SA_ER', 4*np.pi*(3*Ver/(4*np.pi))**1.5)
# Compartment('ER_MEMB', dimension=2, size=SA_ER, parent=CELL_MEMB)
# #  internal ER lumen space
# Parameter('V_ER', Ver)
# Compartment('ER_LUMEN', dimension=3, size=V_ER, parent=ER_MEMB)

# Monomers
# ========
# PAR2 agonist 2AT  Note: this was modeled as Trypsin or Tryp in an earlier
# version of model, but in experiments PAR2 is actually activated by the
# agonist 2AT, Kang et al. https://doi.org/10.1021/acsnano.9b01993)
Monomer('TAT', ['b'])
# PAR2, states: I = inactive (uncleaved), A = active (cleaved)
Monomer('PAR2', ['b','state'], {'state': ['I','A']})
# G-alpha_q protein, states: I = inactive, A = active
Monomer('Gaq', ['b','state'], {'state': ['I', 'A']})
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
# PAR2 agonist 2AT, Kang et al. https://doi.org/10.1021/acsnano.9b01993
Parameter('TAT_0', 330e-3)
Initial(TAT(b=None), TAT_0)
# inactive PAR2
Parameter('PAR2_0', SPC)
Initial(PAR2(state='I', b=None), PAR2_0)
# inactive Gaq
Parameter('Gaq_0', SPC)
Initial(Gaq(state='I', b=None), Gaq_0)
# inactive PLC
Parameter('PLC_0', SPC)
Initial(PLC(bgaq=None, bpip2=None), PLC_0)
# PIP2
Parameter('PIP2_0', SPC)
Initial(PIP2(b=None), PIP2_0)
# IP3R
Parameter('IP3R_0', SPC)
Initial(IP3R(b1=None, b2=None, b3=None, b4=None, bcaer=None, bcacyt=None), IP3R_0)
# ER Ca2+ store
# ER lumen of HEK-293 cells has between roughly 400-600 microM with an average
# around 525 microM as reported in
# Foyouzi-Youssefi et al. https://doi.org/10.1073/pnas.97.11.5723 (Fig. 3C, control)
Parameter('Ca_0', 525)
Initial(Ca(loc='E', b=None), Ca_0)
# TN-XXL
Parameter('TNXXL_0', 0.08)
Initial(TNXXL(bca=None), TNXXL_0)

# Kinetic Parameters
# ==================
# PAR2 activation by 2AT
# Ca2+ signal Max. FRET Dose-Response for 2AT activation of PAR2
# has EC50 = 101.7 +- 28.7 nM, Kang et al. https://doi.org/10.1021/acsnano.9b01993
Parameter('kf_PAR2_bind_TAT', KF)
Parameter('kr_PAR2_bind_TAT', KR)
Parameter('kcat_activate_PAR2', KCAT)
# Gaq activation by activated-PAR2
Parameter('kf_PAR2_bind_Gaq', KF)
Parameter('kr_PAR2_bind_Gaq', KR)
Parameter('kcat_activate_Gaq', KCAT)
# PLC binding Gaq
Parameter('kf_PLC_bind_Gaq', KF)
Parameter('kr_PLC_bind_Gaq', KR)
# Conversion of PIP2 to IP3
Parameter('kf_PLC_bind_PIP2', KF)
Parameter('kr_PLC_bind_PIP2', KR)
Parameter('kcat_PIP2_to_IP3', KCAT)
# Binding of IP3 to IP3R
Parameter('kf_IP3_bind_IP3R', KF)
Parameter('kr_IP3_bind_IP3R', KR)
# Transport of Ca2+
#  ER -> cytosol:
Parameter('kf_erCa_bind_IP3R', KF)
Parameter('kr_erCa_bind_IP3R', KR)
Parameter('kcat_tranport_erCa', KCAT)
#  cytosol -> ER:
Parameter('kf_cytCa_bind_IP3R', KF)
Parameter('kr_cytCa_bind_IP3R', KR)
Parameter('kcat_tranport_cytCa', KCAT)
# Ca2+ binding to TN-XXL FRET reporter
#  Kd = 800 nM, https://doi.org/10.1038/nmeth.1243
Parameter('Kd_cytCa_bind_TNXXL', 800e-3)
# Diffusion coefficient of Ca2+ is 5.3 x 10^-6 cm^2/s ; https://doi.org/10.1016/0143-4160(87)90027-3
# Diffusion-controlled rate: kf = 4*pi*D*R_o/V
Parameter('kf_cytCa_bind_TNXXL', 2)
Expression('kr_cytCa_bind_TNXXL', kf_cytCa_bind_TNXXL*Kd_cytCa_bind_TNXXL)
#Parameter('kr_cytCa_bind_TNXXL', )

# Rules
# =====
# PAR2 activation by 2AT:
#    2AT + PAR2_I <---> TAT:PAR2_I ---> TAT + PAR2_A
catalyze(TAT(), 'b', PAR2(state='I'), 'b', PAR2(state='A'),
         (kf_PAR2_bind_TAT,kr_PAR2_bind_TAT, kcat_activate_PAR2))
# Gaq activation by activated-PAR2:
#    PAR2_A + Gaq_I <---> PAR2_A:Gaq_I ---> PAR2_A + Gaq_A
catalyze(PAR2(state='A'), 'b', Gaq(state='I'), 'b', Gaq(state='A'),
         (kf_PAR2_bind_Gaq,kr_PAR2_bind_Gaq,kcat_activate_Gaq))
# PLC activation by binding Gaq:
#    Gaq_A + PLC <---> Gaq_A:PLC
bind(Gaq(state='A'), 'b', PLC(), 'bgaq', (kf_PLC_bind_Gaq,kr_PLC_bind_Gaq))
# Conversion of PIP2 to IP3
#    Gaq_A:PLC + PIP2 <---> Gaq_A:PLC:PIP2 ---> Gaq_A:PLC + IP3
catalyze_complex(PLC(bgaq=1) % Gaq(b=1, state='A'), 'bpip2', PIP2(), 'b', IP3(),
                (kf_PLC_bind_PIP2,kr_PLC_bind_PIP2,kcat_PIP2_to_IP3))
# Binding of IP3 to IP3R - IP3R is activated when all 4 subunits are bound
#   IP3R + IP3 <---> IP3R:IP3, subunit 1
bind(IP3R(), 'b1', IP3(), 'b', (kf_IP3_bind_IP3R,kr_IP3_bind_IP3R))
#   IP3R + IP3 <---> IP3R:IP3, subunit 2
bind_complex(IP3R(b1=1) % IP3(b=1), 'b2', IP3(), 'b',
             (kf_IP3_bind_IP3R,kr_IP3_bind_IP3R))
#   IP3R + IP3 <---> IP3R:IP3, subunit 3
bind_complex(IP3R(b1=1, b2=2) % IP3(b=1) % IP3(b=2), 'b3', IP3(), 'b',
             (kf_IP3_bind_IP3R,kr_IP3_bind_IP3R))
#   IP3R + IP3 <---> IP3R:IP3, subunit 4
bind_complex(IP3R(b1=1, b2=2,b3=3) % IP3(b=1) % IP3(b=2) % IP3(b=3), 'b4',
             IP3(), 'b', (kf_IP3_bind_IP3R,kr_IP3_bind_IP3R))
# Transport of Ca2+ by activated IP3R
#  ER -> cytosol:
#    IP3R:IP3_4 + Ca_E <---> Ca_E:IP3R:IP3_4 ---> Ca_C + IP3R:IP3_4
catalyze_complex(IP3R(b1=1, b2=2, b3=3, b4=4) % IP3(b=1) % IP3(b=2)
                 % IP3(b=3) % IP3(b=4), 'bcaer', Ca(loc='E'), 'b', Ca(loc='C'),
                 (kf_erCa_bind_IP3R, kr_erCa_bind_IP3R, kcat_tranport_erCa))
#  Reverse, cytosol -> ER:
#    IP3R:IP3_4 + Ca_C <---> Ca_C:IP3R:IP3_4 ---> Ca_E + IP3R:IP3_4
catalyze_complex(IP3R(b1=1, b2=2, b3=3, b4=4) % IP3(b=1) % IP3(b=2)
                 % IP3(b=3) % IP3(b=4), 'bcacyt', Ca(loc='C'), 'b', Ca(loc='E'),
                 (kf_cytCa_bind_IP3R, kr_cytCa_bind_IP3R, kcat_tranport_cytCa))

# Binding of Calcium to the TN-XXL FRET reporter
bind(TNXXL(), 'bca', Ca(loc='C'), 'b',
     (kf_cytCa_bind_TNXXL,kr_cytCa_bind_TNXXL))

# Observables
# ===========
Observable('iPAR2', PAR2(state='I'))
Observable('aPAR2', PAR2(state='A'))
Observable('cytoCa', Ca(loc='C', b=None))
Observable('FRETGen', TNXXL()%Ca())

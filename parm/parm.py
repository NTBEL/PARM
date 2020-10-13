'''PARM: PAR2 Activation-driven calcium Release Model
'''

from pysb import Model, Monomer, Parameter, Initial, Rule, Observable, Expression, Annotation
from pysb.macros import bind, bind_complex, catalyze, catalyze_complex

Model()

# Default forward, reverse, and catalytic rates:
KF = 1e-6
KR = 1e-3
KC = 1

# Compartments
# ============
# Cytosol
#Compartment('cytosol', dimension=3, size=1.0)
#  ER membrane
#Compartment('ER_MEMB', dimension=2, size=1.0, parent=cytosol)
#  internal ER space
#Compartment('ER', dimension=3, size=0.25, parent=ER_MEMB)

# Mononers
# ========
# Trypsin
Monomer('Tryp', ['b'])
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
# Tyrpsin
Parameter('Tryp_0', 0.5)
Initial(Tryp(b=None), Tryp_0)
# inactive PAR2
Parameter('PAR2_0', 0.05)
Initial(PAR2(state='I', b=None), PAR2_0)
# inactive Gaq
Parameter('Gaq_0', 0.05)
Initial(Gaq(state='I', b=None), Gaq_0)
# inactive PLC
Parameter('PLC_0', 0.06)
Initial(PLC(bgaq=None, bpip2=None), PLC_0)
# PIP2
Parameter('PIP2_0', 0.06)
Initial(PIP2(b=None), PIP2_0)
# IP3R
Parameter('IP3R_0', 0.01)
Initial(IP3R(b1=None, b2=None, b3=None, b4=None, bcaer=None, bcacyt=None), IP3R_0)
# ER Ca2+ store
Parameter('Ca_0', 0.07)
Initial(Ca(loc='E', b=None), Ca_0)
# TN-XXL
Parameter('TNXXL_0', 0.08)
Initial(TNXXL(bca=None), TNXXL_0)

# Kinetic Parameters
# ==================
# PAR2 activation/cleavage by Trypsin
Parameter('kf_PAR2_bind_Tryp', KF)
Parameter('kr_PAR2_bind_Tryp', KR)
Parameter('kc_activate_PAR2', KC)
# Gaq activation by activated-PAR2
Parameter('kf_PAR2_bind_Gaq', KF)
Parameter('kr_PAR2_bind_Gaq', KR)
Parameter('kc_activate_Gaq', KC)
# PLC binding Gaq
Parameter('kf_PLC_bind_Gaq', KF)
Parameter('kr_PLC_bind_Gaq', KR)
# Conversion of PIP2 to IP3
Parameter('kf_PLC_bind_PIP2', KF)
Parameter('kr_PLC_bind_PIP2', KR)
Parameter('kc_PIP2_to_IP3', KC)
# Binding of IP3 to IP3R
Parameter('kf_IP3_bind_IP3R', KF)
Parameter('kr_IP3_bind_IP3R', KR)
# Transport of Ca2+
#  ER -> cytosol:
Parameter('kf_erCa_bind_IP3R', KF)
Parameter('kr_erCa_bind_IP3R', KR)
Parameter('kc_tranport_erCa', KC)
#  cytosol -> ER:
Parameter('kf_cytCa_bind_IP3R', KF)
Parameter('kr_cytCa_bind_IP3R', KR)
Parameter('kc_tranport_cytCa', KC)
# Ca2+ binding to TN-XXL FRET reporter
#  Kd = 800 nM, https://doi.org/10.1038/nmeth.1243
Parameter('Kd_cytCa_bind_TNXXL', 800e3)
Parameter('kf_cytCa_bind_TNXXL', 2)
Expression('kr_cytCa_bind_TNXXL', kf_cytCa_bind_TNXXL*Kd_cytCa_bind_TNXXL)
#Parameter('kr_cytCa_bind_TNXXL', )

# Rules
# =====
# PAR2 activation/cleavage by Trypsin:
#    Tryp + PAR2_I <---> Tryp:PAR2_I ---> Tryp + PAR2_A
catalyze(Tryp(), 'b', PAR2(state='I'), 'b', PAR2(state='A'),
         (kf_PAR2_bind_Tryp,kr_PAR2_bind_Tryp, kc_activate_PAR2))
# Gaq activation by activated-PAR2:
#    PAR2_A + Gaq_I <---> PAR2_A:Gaq_I ---> PAR2_A + Gaq_A
catalyze(PAR2(state='A'), 'b', Gaq(state='I'), 'b', Gaq(state='A'),
         (kf_PAR2_bind_Gaq,kr_PAR2_bind_Gaq,kc_activate_Gaq))
# PLC activation by binding Gaq:
#    Gaq_A + PLC <---> Gaq_A:PLC
bind(Gaq(state='A'), 'b', PLC(), 'bgaq', (kf_PLC_bind_Gaq,kr_PLC_bind_Gaq))
# Conversion of PIP2 to IP3
#    Gaq_A:PLC + PIP2 <---> Gaq_A:PLC:PIP2 ---> Gaq_A:PLC + IP3
catalyze_complex(PLC(bgaq=1) % Gaq(b=1, state='A'), 'bpip2', PIP2(), 'b', IP3(),
                (kf_PLC_bind_PIP2,kr_PLC_bind_PIP2,kc_PIP2_to_IP3))
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
                 (kf_erCa_bind_IP3R, kr_erCa_bind_IP3R, kc_tranport_erCa))
#  Reverse, cytosol -> ER:
#    IP3R:IP3_4 + Ca_C <---> Ca_C:IP3R:IP3_4 ---> Ca_E + IP3R:IP3_4
catalyze_complex(IP3R(b1=1, b2=2, b3=3, b4=4) % IP3(b=1) % IP3(b=2)
                 % IP3(b=3) % IP3(b=4), 'bcacyt', Ca(loc='C'), 'b', Ca(loc='E'),
                 (kf_cytCa_bind_IP3R, kr_cytCa_bind_IP3R, kc_tranport_cytCa))

# Binding of Calcium to the TN-XXL FRET reporter
bind(TNXXL(), 'bca', Ca(loc='C'), 'b',
     (kf_cytCa_bind_TNXXL,kr_cytCa_bind_TNXXL))

# Observables
# ===========
Observable('iPAR2', PAR2(state='I'))
Observable('aPAR2', PAR2(state='A'))
Observable('cytoCa', Ca(loc='C', b=None))
Observable('FRETGen', TNXXL()%Ca())

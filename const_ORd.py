from math import exp, sqrt

const_d = {

  # Initial values of variables
  'v_'           : -87.5          , # Membrane voltage
  'dt_'          : 0.001          , # Time step (ms)
  'm_'           : 0.0            , # Na activation
  'hf_'          : 1.0            , # Na inactivation
  'hs_'          : 1.0            , #
  'hsp_'         : 1.0            , #
  'j_'           : 1.0            , # Na inactivation
  'jp_'          : 1.0            , #
  'mL_'          : 0.0            , # NaL activation
  'hL_'          : 1.0            , # NaL inactivation
  'hLp_'         : 1.0            , #
  'a_'           : 0.0            , #
  'iF_'          : 1.0            , #
  'iS_'          : 1.0            , #
  'ap_'          : 0.0            , #
  'iFp_'         : 1.0            , #
  'iSp_'         : 1.0            , #
  'd_'           : 0.0            , #
  'ff_'          : 1.0            , #
  'fs_'          : 1.0            , #
  'fcaf_'        : 1.0            , #
  'fcas_'        : 1.0            , #
  'jca_'         : 1.0            , #
  'ffp_'         : 1.0            , #
  'fcafp_'       : 1.0            , #
  'nca_'         : 0.0            , #
  'xrf_'         : 0.0            , #
  'xrs_'         : 0.0            , #
  'xs1_'         : 0.0            , # Slowly Activating K time-dependant activation
  'xs2_'         : 0.0            , # Slowly Activating K time-dependant activation
  'xk1_'         : 1.0            , #
  'nai_'         : 7.0            , # Intracellular Na (mM)
  'nass_'        : 7.0            , # Subspace Na (mM)
  'ki_'          : 145.0          , # Intracellular K (mM)
  'kss_'         : 145.0          , # Subspace Na (mM)
  'cai_'         : 0.0001         , # Intracellular Ca (mM)
  'cass_'        : 0.0001         , # Subspace Ca (mM)
  'Jrelnp_'      : 0.0            , #
  'Jrelp_'       : 0.0            , #
  'CaMKt_'       : 0.0            , # calmodulin-dependent protein kinase II total?
  'st_'          : 0.0            , # Stimulus current (uA/cm^2)

  # buffer paramaters
  'BSRmax_'      : 0.047          ,
  'KmBSR_'       : 0.00087        ,
  'BSLmax_'      : 1.124          ,
  'KmBSL_'       : 0.0087         ,
  'cmdnmax_'     : 0.05           ,
  'kmcmdn_'      : 0.00238        ,
  'trpnmax_'     : 0.07           ,
  'kmtrpn_'      : 0.0005         ,
  'csqnmax_'     : 10.0           ,
  'kmcsqn_'      : 0.8            ,

  # CaMK parameters
  'CaMKo_'       : 0.05           ,
  'aCaMK_'       : 0.05           ,
  'bCaMK_'       : 0.00068        ,
  'KmCaM_'       : 0.0015         ,
  'KmCaMK_'      : 0.15           ,

  # Terms for Solution of Conductance and Reversal Potential
  'R_'    : 8314 ,     # Universal Gas Constant (J/kmol*K)
  'frdy_' : 96485,     # Faraday's Constant (C/mol)
  'temp_' : 310  ,     # Temperature (K)

  # Extracellular ion concentration
  'nao_'   : 140.0          , # Extracellular Na (mM)
  'cao_'   : 1.8            , # Extracellular Ca (mM)
  'ko_'    : 5.4            , # Extracellular K (mM)

  # Cell Geometry
  'l_'       : 0.01         , # Length of the cell (cm)
  'a_'       : 0.0011       , # Radius of the cell (cm)
  'pi_'      : 3.141592     , # Pi
  'vcell_'   : 3.80132711e-5, # Cell volume (uL)                1000*pi_*a_*a_*l_
  'ageo_'    : 7.671769e-4  , # Geometric membrane area (cm^2)  2*pi_*a_*a_+2*pi_*a_*l_
  'acap_'    : 1.5343538e-5 , # Capacitive membrane area (cm^2) ageo_*2
  'vmyo_'    : 2.58490244e-5, # Myoplasm volume (uL)            vcell_*0.68
  'vmito_'   : 9.88344843e-6, # Mitochondria volume (uL)        vcell_*0.26
  'vsr_'     : 2.28079627e-6, # SR volume (uL)                  vcell_*0.06
  'vnsr_'    : 2.09833257e-6, # NSR volume (uL)                 vcell_*0.0552
  'vjsr_'    : 1.82463701e-7, # JSR volume (uL)                 vcell_*0.0048
  'volrss_'  : 7.60265422e-7, # Restricted Subspace (uL)        vcell_*0.02

  # Ion Valences
  'zna_' : 1,       # Na valence
  'zk_'  : 1,       # K valence
  'zca_' : 2,       # Ca valence

  # Slowly Activating Potassium Current
  'prnak_' : 0.01833, # Na/K Permiability Ratio

  # Q10 for conductance of ion currents
  'Q10NA_'   : 1.5,
  'Q10CAL_'  : 2.96,
  'Q10CAT_'  : 2.5,
  'Q10KR_'   : 1.94,
  'Q10KS_'   : 2.2,
  'Q10K1_'   : 1.5,
  'Q10TO_'   : 6.14,
  'Q10NACA_' : 2.2,
  'Q10NAK_'  : 1.87,
  'Q10REL_'  : 1.68,

  # Shift in steady-state (in)activation curves [mV/K]
  'SINFM_' : 0.8,
  'SINFHJ_' : 0.7,
  'SINFXR_' : -1.16,

  # Q10 for time constant of (in)activation
  'Q10TAUMHJ_' : 2.79,
  'Q10TAUD_' : 2.52,
  'Q10TAUF_' : 2.82,
  'Q10TAUXR_' : 6.25,
  'Q10TAUXS_' : 2.58,
}

from math import exp, sqrt

const_d = {

  # Initial values of variables
  'v_'           : -86.24         , # Membrane voltage
  'dt_'          : 0.001          , # Time step (ms)
  'm_'           : 0.000642       , # Na activation
  'h_'           : 0.996271       , # Na inactivation
  'j_'           : 0.997286       , # Na inactivation
  'xr_'          : 0.0001096      , # Rapidly Activating K time-dependant activation
  'b_'           : 0.000880231    , # T-Type Ca Voltage dependant activation gate
  'g_'           : 0.994305       , # T-Type Ca Voltage dependant inactivation gate
  'xs1_'         : 0.004761       , # Slowly Activating K time-dependant activation
  'xs2_'         : 0.03291        , # Slowly Activating K time-dependant activation
  'nai_'         : 10.71          , # Intracellular Na (mM)
  'naiss_'       : 10.71          , # Subspace Na (mM)
  'cai_'         : 0.0000821      , # Intracellular Ca (mM)
  'caiss_'       : 0.0000672      , # Subspace Ca (mM)
  'ki_'          : 137.6          , # Intracellular K (mM)
  'kiss_'        : 137.6          , # Subspace K (mM)
  'nsr_'         : 2.553          , # NSR Ca Concentration (mM)
  'jsr_'         : 2.479          , # JSR Ca Concentration (mM)
  'csqn_'        : 7.56           , # Calsequestrin Buffered Ca Concentration (mM)
  'ilca_'        : 0.0            , # Ca current through L-type Ca channel (uA/uF)
  'ilcana_'      : 0.0            , # Na current through L-type Ca channel (uA/uF)
  'ilcak_'       : 0.0            , # K current through L-type Ca channel (uA/uF)
  'irel_'        : 0.0            , # SR Ca Release
  'itr_'         : 0.0            , # Translocation current of Ca ions from NSR to JSR (mM/ms)
  'inacass_'     : 0.0            , # NaCa exchanger current in the Subspace (uA/uF)
  'fcasc_'       : 1.0            , #
  'fmode0_'      : 0.22           , #
  'ltypeCzero_'  : 0.948          , #
  'ltypeCone_'   : 0.052          , #
  'ltypeCtwo_'   : 0.0            , #
  'ltypeCthree_' : 0.0            , #
  'ltypeIVf_'    : 0.0            , #
  'ltypeIVs_'    : 0.0            , #
  'ltypeO_'      : 0.0            , #
  'ryrCone_'     : 0.89           , #
  'ryrCtwo_'     : 0.021          , #
  'ryrCthree_'   : 0.003          , #
  'ryrCfour_'    : 0.001          , #
  'ryrOone_'     : 0.0            , #
  'ryrIone_'     : 0.082          , #
  'ryrItwo_'     : 0.002          , #
  'ryrIthree_'   : 0.001          , #
  'ryrIfour_'    : 0.0            , #
  'ryrIfive_'    : 0.0            , #
  'sponrel_'     : 0.0            , #
  'it_'          : 0.0            , # Total current (uA/cm^2)
  'st_'          : 0.0            , # Stimulus current (uA/cm^2)

  # Terms for Solution of Conductance and Reversal Potential
  'R_'    : 8314 ,     # Universal Gas Constant (J/kmol*K)
  'frdy_' : 96485,     # Faraday's Constant (C/mol)
  'temp_' : 310.0  ,     # Temperature (K)

  # Extracellular ion concentration
  'nao_'   : 130.0          , # Extracellular Na (mM)
  'cao_'   : 1.8            , # Extracellular Ca (mM)
  'ko_'    : 4.5            , # Extracellular K (mM)

  # Cell Geometry
  'l_'       : 0.01         , # Length of the cell (cm)
  'a_'       : 0.0011       , # Radius of the cell (cm)
  'pi_'      : 3.141592     , # Pi
  'vcell_'   : 3.80132711e-5, # Cell volume (uL)                1000*pi_*a_*a_*l_
  'ageo_'    : 7.671769e-4  , # Geometric membrane area (cm^2)  2*pi_*a_*a_+2*pi_*a_*l_
  'acap_'    : 1.5343538e-5 , # Capacitive membrane area (cm^2) ageo_*2
  'vmyo_'    : 2.58490244e-5, # Myoplasm volume (uL)            vcell_*0.68
  'vmito_'   : 9.12318507e-6, # Mitochondria volume (uL)        vcell_*0.24
  'vsr_'     : 2.28079627e-6, # SR volume (uL)                  vcell_*0.06
  'vnsr_'    : 2.09833257e-6, # NSR volume (uL)                 vcell_*0.0552
  'vjsr_'    : 1.82463701e-7, # JSR volume (uL)                 vcell_*0.0048
  'volrss_'  : 7.60265422e-7, # Restricted Subspace (uL)        vcell_*0.02

  # Ion Valences
  'zna_' : 1,       # Na valence
  'zk_'  : 1,       # K valence
  'zca_' : 2,       # Ca valence

  # JSR Ca Ion Concentration Changes
  'kmcsqn_'  : 0.8,           # Equalibrium constant of buffering for CSQN (mM)
  'betaad_'  : 0  ,           # Set this Value to 1 to Introduce Beta Adrenergic Effects
  'csqnmut_' : 0  ,           # Set this Value to 1 to Introdue CASQ2 Mutation Effects
  'csqnbar_' : 10 ,           # (csqnmut_ == 0) ? 10:4

  # Translocation of Ca Ions from NSR to JSR
  'tautr_' : 120,             # Time constant of Ca transfer from NSR to JSR (ms)
  'taudiffss_' : 0.1,

  # NSR Ca Ion Concentration Changes
  'nsrbar_' : 15,
  'kmup_'   : 0.00092,
  'iupbar_' : 0.017325    , # Max. current through iup channel (mM/ms)
  'kleak_'  : 0.00875 / 15, # 0.00875 / nsrbar_ # Rate constant of Ca leakage from NSR (ms^-1)


  # Slowly Activating Potassium Current
  'prnak_' : 0.01833, # Na/K Permiability Ratio

  # Plateau Potassium Current
  'gkp_' : 0.000552,    # Channel Conductance of Plateau K Current (mS/uF)

  # Sodium-Calcium Exchanger V-S
  'c1_' : 0.000225   ,   # NaCa Scaling factor
  'c2_' : 0.0001     ,   # Half-saturation concentration of NaCa exhanger (mM)
  'c3_' : 1.5        ,   # NaCa exhanger Extracellular Na Dependence
  'gammanaca_' : 0.15,   # Position of energy barrier controlling voltage dependance of inaca

  # Sarcolemmal Ca Pump
  'ibarpca_' : 1,        # Max. Ca current through sarcolemmal Ca pump (uA/uF)

  'spontau_' : 125     , # Time constant of Spontaneous Release

  'pna_' : 0.000008265 , # Permiability of membrane to Na (cm/s)
  'pca_' : 0.006615    , # Permiability of membrane to Ca (cm/s)
  'pk_' : 0.000002363  , # Permiability of membrane to K (cm/s)
  'kmca_' : 0.0006     ,
  'gacai_' : 0.01      , # Activity coefficient of Ca
  'gacao_' : 0.341     , # Activity coefficient of Ca
  'ganai_' : 0.75      , # Activity coefficient of Na
  'ganao_' : 0.75      , # Activity coefficient of Na
  'gaki_' : 0.75       , # Activity coefficient of K
  'gako_' : 0.75       , # Activity coefficient of K
  'gcab_' : 0.003016   , # Max. conductance of Ca background (mS/uF)
  'gnab_' : 0.0002     , # Max. conductance of Na background (mS/uF)
  'gcat_' : 0.005      , # Max. Conductance of the T-type Ca channel (mS/uF)
  'gna_'  : 16         , # Max. Conductance of the Na Channel (mS/uF)
  'gkr_'   : 0.031368*sqrt(4.5/5.4), #0.031368*sqrt(ko_/5.4),
  'gki_'   : 0.5625*(sqrt(4.5/5.4)), #0.5625*(sqrt(ko_/5.4)),

  'kmbsr_' : 0.00087,
  'kmbsl_':0.127,
  'kmko_' : 1.5,
  'kmnai_' : 10.,
  'bslbar_':2.124,
  'bsrbar_' : 0.047,
  'ibarnak_' : 2.25,
  'cmdnbar_' : 0.050,    # Max. [Ca] buffered in CMDN (mM)
  'trpnbar_' : 0.070,    # Max. [Ca] buffered in TRPN (mM)
  'kmcmdn_'  : 0.00238,  # Equalibrium constant of buffering for CMDN (mM)
  'kmtrpn_'  : 0.0005,   # Equalibrium constant of buffering for TRPN (mM)

  'sigma_' : (exp(130.0/67.3)-1)/7,  #(exp(nao_/67.3)-1)/7

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
  'SINFXR_' : 1.16,

  # Q10 for time constant of (in)activation
  'Q10TAUMHJ_' : 2.79,
  'Q10TAUD_' : 2.52,
  'Q10TAUF_' : 2.82,
  'Q10TAUBG_' : 2.5,
  'Q10TAUXR_' : 6.25,
  'Q10TAUXS_' : 2.58,
  'Q10OTHERS_' : 3.0,
}

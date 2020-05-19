


const_d = {

# Initial values of variables */
  'v_'  : -81.2,       # Initial Voltage (mv) */
  'it_' : 0,
 # 't_' : 0,           # Time (ms) */
 #'udt_' : 0.01,     # Time step (ms) */
  'st_' : 0,        # Stimulus */
  'nai_' : 11.2,       # Initial Intracellular Na (mM) */
  'ki_' : 139,       # Initial Intracellular K (mM) */
  'cai_' : 0.000102,  # Initial Intracellular Ca (mM) */
  'ina_': 67.6,
  'inab_': 0.1,
  'iki_': 7.18,
  'inak_': 0.257,
  'inaca_':-3.21,
  'ikr_': 0.00000670,
  'iks_': 0.00034969,
  'ito_':0.0000370,
  'ikur_':0.00000000489,
  'ikach_':0.0261,
  'ilcatot_':-0.00192,
  # Initial Gate Conditions */
  'm_' : 0.00291,
  'h_' : 0.965,
  'j_' : 0.978,
  'd_' : 0.000137,
  'f_' : 0.999837,
  'xs_' : 0.0187,
  'xr_' : 0.0000329,
  'ato_' : 0.0304,
  'iito_' : 0.999,
  'uakur_' : 0.00496,
  'uikur_' : 0.999,
  
  'fca_' : 0.775,
  'ireljsrol_':0,
  # Initial values of variable coefficients
  'fibro_' : 1.0, # Brugada modification parameter
 
  # Initial Conditions */
  'jsr_' : 1.49,
  'nsr_' : 1.49,
  'trpn_' : 0.0118,
  'cmdn_' : 0.00205,
 #'csqn_' : 6.51,
  'dt_' : 0.01,
 # 'utsc_' : 50,
  'urel_' : 0.00,
  'vrel_' : 1.00,
  'wrel_' : 0.999,
  'yach_' : 2.54e-2,

# Constants*/

  'l_' : 0.01,       # Length of the cell (cm) */
  'a_' : 0.0008,     # Radius of the cell (cm) */
  'pi_' : 3.141592,  # Pi */
  'R_' : 8314,      # Universal Gas Constant (J/kmol*K) */
  'frdy_' : 96485, # Faraday's Constant (C/mol) */
  'temp_':310, #temperature */
  'nao_' : 140,      # Initial Extracellular Na (mM) */
  'ko_' : 4.5,       # Initial Extracellular K (mM) */
  'cao_' : 1.8,      # Initial Extracellular Ca (mM) */
  'APD_change_h_':1.0,
  'APD_change_j_':1.0,
  'APD_change_m_':1.0,
  'APD_change_gna_':1.0,  
  'APD_change_inak_':1.0,
  'APD_change_inab_':1.0,
  'APD_change_ina_':1.0,
  'APD_change_ikur_':3.0,
  'APD_change_ilca_':0.5,
  'thickness_':1.0,
  'sw_it_':1.0,

  # Ion Valences */
  'zna_' : 1,  # Na valence */
  'zk_' : 1,   # K valence */
  'zca_' : 2,  # Ca valence */
  'kmup_' : 0.00092,    # Half-saturation concentration of iup (mM) */
  'iupbar_' : 0.005,  # Max. current through iup channel (mM/ms) */
  'nsrbar_' : 15,
  'grelbarjsrol_' : 30, # Rate constant of Ca release from JSR due to overload (ms^-1)*/
  'csqnbar_' : 10,      # Max. [Ca] buffered in CSQN (mM)*/
  'kmcsqn_' : 0.8,      # Equalibrium constant of buffering for CSQN (mM)*/
  'tautr_' : 180,  # Time constant of Ca transfer from NSR to JSR (ms)*/
  'cmdnbar_' : 0.050,   # Max. [Ca] buffered in CMDN (mM) */
  'trpnbar_' : 0.070,   # Max. [Ca] buffered in TRPN (mM) */
  'kmcmdn_' : 0.00238,  # Equalibrium constant of buffering for CMDN (mM) */
  'kmtrpn_' : 0.0005,   # Equalibrium constant of buffering for TRPN (mM) */
  'gcalbar_' : 0.1238,
  'ach_' : 0.0, # Acetylcholine concentration */
  'prnak_' : 0.01833,  # Na/K Permiability Ratio */
  'kmnancx_' : 87.5,  # Na saturation constant for NaCa exchanger */
  'ksatncx_' : 0.1,   # Saturation factor for NaCa exchanger */
  'kmcancx_' : 1.38,  # Ca saturation factor for NaCa exchanger */
  'gammas_' : 0.35,  # Position of energy barrier controlling voltage dependance of inaca */
  'ibarnak_' : 1.0933,   # Max. current through Na-K pump (uA/uF) */
  'kmnai_' : 10,    # Half-saturation concentration of NaK pump (mM) */
  'kmko_' : 1.5,    # Half-saturation concentration of NaK pump (mM) */
  'ibarpca_' : 0.275, # Max. Ca current through sarcolemmal Ca pump (uA/uF) */
  'kmpca_' : 0.0005, # Half-saturation concentration of sarcolemmal Ca pump (mM) */
  'ekr_' : -89.19,
  'eks_' : -89.19,
  'eki_' : -89.19,
  'ekach_' : -89.19,
  'erevto_' : -89.19,
  'gki_' :0.08367,
  'kmnai_' : 10,
    

  }

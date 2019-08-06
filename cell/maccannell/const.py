from math import exp, sqrt, log

const_d = {
    
    # Initial values of variables    
    'v_'           : 0.0               , # membrane potential (mV)
    'dt_'          : 0.001             , # time step (ms)
    'temp_'        : 310.15            , # temprature (K)
    'sw_it_'       : 1.0               , # ion current switch 
    'st_'          : 0.0               , # stimulus current (uA/cm^2)
    'it_'          : 0.0               , # total ion current (uA/cm^2)
    'k_i_'         : 129.4349          , # fibroblast intracellular K+ concentration (mmol/L)
    'na_i_'        : 8.5547            , # fibroblast intracellular Na+ concentration (mmol/L)
    'r_'           : 0.0               , # activation parameter
    's_'           : 0.0               , # inactivation parameter
    'ikv_'         : 0.0               , # time and voltage-dependent K+ current (uA/cm^2)
    'ik1_'         : 0.0               , # inward rectifying K+ current (uA/cm^2)
    'inak_'        : 0.0               , # Na+-K+ pump current (uA/cm^2)
    'inab_'        : 0.0               , # Na+ background current (uA/cm^2)
    'ik_net_'      : 0.0               , #
    'ina_net_'     : 0.0               , #
    
    # Constants
    'R'            : 8.314472          , # gas constant (J.mol^-1.K^-1)
    'F'            : 96.4853415        , # Faraday constant (C/mmol)
    'G_KV'         : 0.25              , # maximum conductance of I_KV (mS/uF)
    'G_K1'         : 0.4822            , # maximum conductance of I_K1 (mS/uF)
    'G_NAB'        : 9.5e-3            , # leak conductance of I_Nab (mS/uF)
    'P_NAK'        : 2.002             , # maximum current of I_NaK (uA/uF)
    'V_REV'        : -150              , # fibroblast reversal potential for I_NaK (mV)
    'K_O'          : 5.4               , # extracellular K+ concentration (mmol/L)
    'NA_O'         : 140.0             , # extracellular Na+ concentration (mmol/L)
    'K_MK'         : 1.0               , # binding constant (mmol/L)
    'K_MNA'        : 11.0              , # binding constant (mmol/L)
    'ACAP'         : 1.5343538e-5      , # capacitive membrane area (cm^2)
    'V_MYO'        : 2.58490244e-5     , # myoplasm volume (uL)
    
}

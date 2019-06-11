from math import exp, sqrt

const_d = {

  # Initial values of variables
  'v_'           : -87.169816169406  , # Membrane voltage
  'temp_'        : 308               , # Temperature (K)
  'dt_'          : 0.001             , # Time step (ms)
  'm_'           : 0.001075453357    , # Na activation
  'h_'           : 0.990691306716    , # Na inactivation
  'j_'           : 0.993888937283    , # Na inactivation
  'dyad_'        : 1.716573130685    , #
  'c1_'          : 0.000018211252    , #
  'c2_'          : 0.979322592773    , #
  'xi1ca_'       : 0.001208153482    , #
  'xi1ba_'       : 0.000033616596    , #
  'xi2ca_'       : 0.004173008466    , #
  'xi2ba_'       : 0.015242594688    , #
  'xr_'          : 0.007074239331    , # Rapidly Activating K time-dependant activation
  'cai_'         : 0.256752008084    , # Intracellular Ca (mM)
  'xs1_'         : 0.048267587131    , # Slowly Activating K time-dependant activation
  'xs2_'         : 0.105468807033    , # Slowly Activating K time-dependant activation
  'xtos_'        : 0.00364776906     , #
  'xtof_'        : 0.003643592594    , #
  'ytos_'        : 0.174403618112    , #
  'ytof_'        : 0.993331326442    , #
  'nai_'         : 11.441712311614   , # Intracellular Na (mM)
  'submem_'      : 0.226941113355    , #
  'nsr_'         : 104.450004990523  , #
  'jsr_'         : 97.505463697266   , #
  'kj_'          : 50.0              , #
  'xir_'         : 0.006679257264    , #
  'tropi_'       : 22.171689894953   , #
  'trops_'       : 19.864701949854   , #
  'it_'          : 0.0               , # total ion current
  'st_'          : 0.0               , # Stimulus current (uA/cm^2)
  'xina_'        : 0.0               , #
  'xik1_'        : 0.0               , #
  'xikr_'        : 0.0               , #
  'xiks_'        : 0.0               , #
  'xito_'        : 0.0               , #
  'xitof_'       : 0.0               , #
  'xitos_'       : 0.0               , #
  'xiNaCa_'      : 0.0               , #
  'xica_'        : 0.0               , #
  'xiNaK_'       : 0.0               , #

  # Initial values of variable coefficients
  'c_brugada_'   : 1.4               ,

  # Initial values of constants
  'vth_'         : 0.0               , #
  's1t_'         : 0.00195           , #
  's2t_'         : 0.00105763        , #
  's6_'          : 8.0               , #
  'taupo_'       : 1.0               , #
  'taua_'        : 100.0             , #
  'cao_'         : 1.8               , #
  'cat_'         : 3.0               , #
  'tca_'         : 78.0329           , #
  'cpt_'         : 6.09365           , #
  'vy_'          : 40.0              , #
  'sy_'          : 4.0               , #
  'vx_'          : 40.0              , #
  'sx_'          : 3.0               , #
  'vyr_'         : 40.0              , #
  'syr_'         : 11.32             , #
  'tau3_'        : 3.0               , #
  'k2_'          : 1.03615e-4        , #
  'r1_'          : 0.3               , #
  'r2_'          : 3.0               , #
  'k1t_'         : 0.00413           , #
  'k2t_'         : 0.00224           , #
  'vup_'         : 0.4               , #
  'cup_'         : 0.5               , #
  'gleak_'       : 0.00002069        , #
  'ki_'          : 140.0             , #
  'bcal_'        : 24.0              , #
  'xkcal_'       : 7.0               , #
  'srmax_'       : 47.0              , #
  'srkd_'        : 0.6               , #
  'bmem_'        : 15.0              , #
  'kmem_'        : 0.3               , #
  'bsar_'        : 42.0              , #
  'ksar_'        : 13.0              , #
  'xkon_'        : 0.0327            , #
  'btrop_'       : 70.0              , #
  'xkoff_'       : 0.0196            , #
  'taud_'        : 4.0               , #
  'cstar_'       : 90.0              , #
  'av_'          : 11.3              , #
  'bv_'          : 977.0             , #
  'pca_'         : 0.00054           , #
  'F_'           : 96.4853415        , #
  'R_'           : 8.314472          , #
  'ax_'          : 0.3576            , #
  'ay_'          : 0.05              , #
  'gryr_'        : 2.58079           , #
  'taur_'        : 30.0              , #
  'gbarsr_'      : 26841.8           , #
  'gdyad_'       : 9000.0            , #
  'taups_'       : 0.5               , #
  'gca_'         : 182.0             , #
  'xkdna_'       : 0.3               , #
  'nao_'         : 136.0             , #
  'xmcao_'       : 1.3               , #
  'xmnao_'       : 87.5              , #
  'xmnai_'       : 12.3              , #
  'xmcai_'       : 0.0036            , #
  'gNaCa_'       : 0.84              , #
  'sigma_'       : 0.934910849       , #
  'gNaK_'        : 1.5               , #
  'gna_'         : 12.0              , #
  'xkmnai_'      : 12.0              , #
  'ko_'          : 5.4               , #
  'xkmko_'       : 1.5               , #
  'wca_'         : 8.0               , #
  'gkix_'        : 0.3               , #
  'gtos_'        : 0.04              , #
  'gtof_'        : 0.11              , #
  'gkr_'         : 0.0125            , #
  'prNaK_'       : 0.01833           , #
  'ki_'          : 140.0             , #
  'gks_'         : 0.1386            , #

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
  'Q10TAUBG_' : 2.5,
  'Q10TAUXR_' : 6.25,
  'Q10TAUXS_' : 2.58,
  'Q10OTHER_' : 2.50,
}

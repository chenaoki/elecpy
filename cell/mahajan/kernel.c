int _maskPosi;
int _maskNega;
int _NegaPosi;
int _NegaNega;

_v = v;

// voltage trapping to avoid error
_maskPosi =  ( _v < 0.0001 ) * 1;
_maskPosi *= ( _v >= 0 ) * 1;
_maskNega = ( _maskPosi == 0 ) * 1;
_v = _maskPosi * 0.0001 + _maskNega * _v;
_maskPosi =  ( _v > -0.0001 ) * 1;
_maskPosi *= ( _v <= 0 ) * 1;
_maskNega = ( _maskPosi == 0 ) * 1;
_v = _maskPosi * -0.0001 + _maskNega * _v;

_maskPosi =  ( _v < (-47.1300+0.0001) ) * 1;
_maskPosi *= ( _v >= -47.1300 ) * 1;
_maskNega = ( _maskPosi == 0 ) * 1;
_v = _maskPosi * (-47.1300+0.0001)+ _maskNega * _v;
_maskPosi =  ( _v > (-47.1300 -0.0001) ) * 1;
_maskPosi *= ( _v <= -47.1300 ) * 1;
_maskNega = ( _maskPosi == 0 ) * 1;
_v = _maskPosi * (-47.13-0.0001) + _maskNega * _v;

_maskPosi =  ( _v < (-30.0+0.0001) ) * 1;
_maskPosi *= ( _v >= -30.0 ) * 1;
_maskNega = ( _maskPosi == 0 ) * 1;
_v = _maskPosi * (-30.0+0.0001)+ _maskNega * _v;
_maskPosi =  ( _v > (-30.0 -0.0001) ) * 1;
_maskPosi *= ( _v <= -30.0 ) * 1;
_maskNega = ( _maskPosi == 0 ) * 1;
_v = _maskPosi * (-30.0-0.0001) + _maskNega * _v;

_maskPosi =  ( _v < (-7.0+0.0001) ) * 1;
_maskPosi *= ( _v >= -7.0 ) * 1;
_maskNega = ( _maskPosi == 0 ) * 1;
_v = _maskPosi * (-7.0+0.0001)+ _maskNega * _v;
_maskPosi =  ( _v > (-7.0 -0.0001) ) * 1;
_maskPosi *= ( _v <= -7.0 ) * 1;
_maskNega = ( _maskPosi == 0 ) * 1;
_v = _maskPosi * (-7.0-0.0001) + _maskNega * _v;

_maskPosi =  ( _v < (-10.0+0.0001) ) * 1;
_maskPosi *= ( _v >= -10.0 ) * 1;
_maskNega = ( _maskPosi == 0 ) * 1;
_v = _maskPosi * (-10.0+0.0001)+ _maskNega * _v;
_maskPosi =  ( _v > (-10.0 -0.0001) ) * 1;
_maskPosi *= ( _v <= -10.0 ) * 1;
_maskNega = ( _maskPosi == 0 ) * 1;
_v = _maskPosi * (-10.0-0.0001) + _maskNega * _v;

// reversible potential
double _ena = (1.00000/({F_}/({R_}*temp)))*log({nao_}/nai);
double _ek =  (1.00000/({F_}/({R_}*temp)))*log({ko_}/{ki_});
double _eks = (1.00000/({F_}/({R_}*temp)))*log(({ko_}+{prNaK_}*{nao_})/({ki_}+{prNaK_}*nai));

// calculate conductance
double _gna = {gna_} * pow({Q10NA_}, (temp-{temp_})/10.0);
double _gkix = {gkix_} * pow({Q10K1_}, (temp-{temp_})/10.0);
double _gkr = {gkr_} * pow({Q10KR_}, (temp-{temp_})/10.0);
double _gks = {gks_} * pow({Q10KS_}, (temp-{temp_})/10.0);
double _gtos = {gtos_} * pow({Q10TO_}, (temp-{temp_})/10.0)*1.4;
double _gtof = {gtof_} * pow({Q10TO_}, (temp-{temp_})/10.0)*1.4;
double _gNaCa = {gNaCa_} * pow({Q10NACA_}, (temp-{temp_})/10.0);
double _gca = {gca_} * pow({Q10CAL_}, (temp-{temp_})/10.0)/1.4;
double _gNaK = {gNaK_} * pow({Q10NAK_}, (temp-{temp_})/10.0);

// INa current
_maskPosi = (_v <= -40.0) * 1;
_maskNega = (_maskPosi == 0) * 1;
double _ah = ((0.135000*exp((80.0000+_v)/(-6.80000))) * _maskPosi + 0 * _maskNega) * pow({Q10TAUMHJ_}, (temp-{temp_})/10.0);
double _bh = ((3.56000*exp(0.0790000*_v)+310000.*exp(0.350000*_v)) * _maskPosi +
         (1.00000/(0.130000*(1.00000+exp((_v+10.6600)/(-11.1000))))) * _maskNega) * pow({Q10TAUMHJ_}, (temp-{temp_})/10.0);
double _aj = ((((-127140.*exp(0.244400*_v)-3.47400e-05*exp(-0.0439100*_v))*1.00000*(_v+37.7800))/(1.00000+exp(0.311000*(_v+79.2300)))) * _maskPosi +
         0 * _maskNega) * pow({Q10TAUMHJ_}, (temp-{temp_})/10.0);
double _bj = (((0.121200*exp(-0.0105200*_v))/(1.00000+exp(-0.137800*(_v+40.1400)))) * _maskPosi +
         ((0.300000*exp(-2.53500e-07*_v))/(1.00000+exp(-0.100000*(_v+32.0000)))) * _maskNega) * pow({Q10TAUMHJ_}, (temp-{temp_})/10.0);
double _am = ((0.320000*1.00000*(_v+47.1300))/(1.00000-exp(-0.100000*(_v+47.1300)))) * pow({Q10TAUMHJ_}, (temp-{temp_})/10.0);
double _bm = (0.0800000*exp(-_v/11.0000)) * pow({Q10TAUMHJ_}, (temp-{temp_})/10.0);

// IKs current
double _xs1ss = 1.00000/(1.00000+exp(-(_v-1.50000)/16.7000));
double _tauxs1 = (1.00000/((7.19000e-05*(_v+30.0000))/(1.00000-exp(-0.148000*(_v+30.0000)))+(0.000131000*(_v+30.0000))/(exp(0.0687000*(_v+30.0000))-1.00000))) / pow({Q10TAUXS_}, (temp-{temp_})/10.0);
double _xs2ss = _xs1ss;
double _tauxs2 =  4.00000*_tauxs1;
double _gksx = 1.00000+0.800000/(1.00000+pow(0.500000/cai, 3.00000));

// IK1 current
double _aki = 1.02000/(1.00000+exp(0.238500*(_v-_ek-59.2150)));
double _bki = (0.491240*exp(0.0803200*(_v-_ek+5.47600))+1.00000*exp(0.0617500*(_v-_ek-594.310)))/(1.00000+exp(-0.514300*(_v-_ek+4.75300)));
double _xkin = _aki/(_aki+_bki);

// IKr current
double _xkrinf = 1.00000/(1.00000+exp(-(_v+50.0000)/7.50000));
double _xkrv1 = ((0.00138000*1.00000*(_v+7.00000))/(1.00000-exp(-0.123000*(_v+7.00000))));
double _xkrv2 = ((0.000610000*1.00000*(_v+10.0000))/(exp( 0.145000*(_v+10.0000))-1.00000));
double _taukr = (1.00000/(_xkrv1+_xkrv2)) / pow({Q10TAUXR_}, (temp-{temp_})/10.0);
double _rg = 1.00000/(1.00000+exp((_v+33.0000)/22.4000));

// Na-K pump current
double _fNaK = 1.00000/(1.00000+0.124500*exp(-0.100000*_v*({F_}/({R_}*temp)))+0.0365000*{sigma_}*exp(-_v*({F_}/({R_}*temp))));

// Ito current
double _rt1 = - (_v+3.00000)/15.0000;
double _rt2 = (_v+33.5000)/10.0000;
double _rt3 = (_v+60.0000)/10.0000;
double _rt4 = ((-_v/30.0000)*_v)/30.0000;
double _rt5 = (_v+33.5000)/10.0000;
// Itos
double _rsinf = 1.00000/(1.00000+exp(_rt2));
double _xtosinf = 1.00000/(1.00000+exp(_rt1));
double _ytosinf = 1.00000/(1.00000+exp(_rt2));
double _txs = (9.00000/(1.00000+exp(-_rt1))+0.500000) / pow({Q10OTHER_}, (temp-{temp_})/10.0);
double _tys = (3000.00/(1.00000+exp(_rt3))+30.0000) / pow({Q10OTHER_}, (temp-{temp_})/10.0);
// Itof
double _xtofinf = _xtosinf;
double _ytofinf = _ytosinf;
double _txf = (3.50000*exp(_rt4)+1.50000) / pow({Q10OTHER_}, (temp-{temp_})/10.0);
double _tyf = (20.0000/(1.00000+exp(_rt5))+20.0000) / pow({Q10OTHER_}, (temp-{temp_})/10.0);

// L-type Ca current
double _Pr = 1.00000 - 1.00000/(1.00000+exp(-(_v+{vy_})/{sy_}));
double _Ps = 1.00000/(1.00000+exp(-(_v+{vyr_})/{syr_}));
double _recov = 10.0000+4954.00*exp(_v/15.6000);
double _tauca1 = {tca_}/(1.00000+pow(dyad/{cpt_}, 4.00000))+0.100000;
//T _tauca1 = ({tca_}+0.1*pow((1+dyad/{cpt_}), 4.00000))/(1.00000+pow(dyad/{cpt_}, 4.00000));
double _tauca2 = ((_recov-_tauca1)*_Pr+_tauca1);
double _tauba = ((_recov-450.000)*_Pr+450.000);
double _poinf = 1.00000/(1.00000+exp(-(_v-{vth_})/{s6_}));
double _alpha = _poinf/{taupo_};
double _beta = (1.00000-_poinf)/{taupo_};
double _fca = 1.00000/(1.00000+pow({cat_}/dyad, 3.00000));
double _s1 =  0.0182688*_fca;
double _k1 =  0.0241680*_fca;
double _s2 = _s1*({r1_}/{r2_})*({k2_}/_k1);
double _poi = 1.00000/(1.00000+exp(-(_v+{vx_})/{sx_}));
double _k3 = (1.00000 - _poi)/{tau3_};
double _k3t = _k3;
double _k5 = (1.00000-_Ps)/_tauca2;
double _k5t = (1.00000-_Ps)/_tauba;
double _k6 = (_fca*_Ps)/_tauca2;
double _k6t = _Ps/_tauba;
double _k4 = _k3*(_alpha/_beta)*(_k1/{k2_})*(_k5/_k6);
double _k4t = _k3t*(_alpha/_beta)*({k1t_}/{k2t_})*(_k5t/_k6t);
double _po = 1.00000 - (c1+c2+xi1ca+xi2ca+xi1ba+xi2ba);

// Diffusive flux
double _jd = (submem - cai)/{taud_};

// SR leak fluc
double _jleak = {gleak_}*(nsr*nsr/(nsr*nsr+{kj_}*{kj_}))*(nsr*16.6670-cai);

// Nonlinear buffering
double _spxs = ({srmax_}*{srkd_})/(({srkd_}+submem)*({srkd_}+submem));
double _bpxs = ({bcal_}*{xkcal_})/(({xkcal_}+submem)*({xkcal_}+submem));
double _mempxs = ({bmem_}*{kmem_})/(({kmem_}+submem)*({kmem_}+submem));
double _sarpxs = ({bsar_}*{ksar_})/(({ksar_}+submem)*({ksar_}+submem));
double _dcsib = 1.00000/(1.00000+_bpxs+_spxs+_mempxs+_sarpxs);
double _spxi = ({srmax_}*{srkd_})/(({srkd_}+cai)*({srkd_}+cai));
double _bpxi = ({bcal_}*{xkcal_})/(({xkcal_}+cai)*({xkcal_}+cai));
double _mempxi = ({bmem_}*{kmem_})/(({kmem_}+cai)*({kmem_}+cai));
double _sarpxi = ({bsar_}*{ksar_})/(({ksar_}+cai)*({ksar_}+cai));
double _dciib = 1.00000/(1.00000+_bpxi+_spxi+_mempxi+_sarpxi);
double _xbi = {xkon_}*cai*({btrop_}-tropi) - {xkoff_}*tropi;
double _xbs = {xkon_}*submem*({btrop_}-trops) - {xkoff_}*trops;

// NaCa exchange flux
double _aloss = 1.00000/(1.00000+pow({xkdna_}/submem, 3.00000));
double _csm = submem/1000.00;
double _zw3 = pow(nai, 3.00000)*{cao_}*exp(_v*0.350000*({F_}/({R_}*temp))) - pow({nao_}, 3.00000)*_csm*exp(_v*(0.350000-1.00000)*({F_}/({R_}*temp)));
double _zw4 = 1.00000+0.200000*exp(_v*(0.350000-1.00000)*{F_}/({R_}*temp));
double _yz1 = {xmcao_}*pow(nai, 3.00000)+pow({xmnao_}, 3.00000)*_csm;
double _yz2 = pow({xmnai_}, 3.00000)*{cao_}*(1.00000+_csm/{xmcai_});
double _yz3 = {xmcai_}*pow({nao_}, 3.00000)*(1.00000+pow(nai/{xmnai_}, 3.00000));
double _yz4 = pow(nai, 3.00000)*{cao_}+pow({nao_}, 3.00000)*_csm;
double _zw8 = _yz1+_yz2+_yz3+_yz4;
double _jNaCa = (_gNaCa*_aloss*_zw3)/(_zw4*_zw8);

// SERCA(uptake) pump
double _jup = ({vup_}*cai*cai)/(cai*cai+{cup_}*{cup_});

// L-type Ca current flux
double _za =  _v*2.00000*({F_}/({R_}*temp));
_maskPosi = (fabs(_za)<0.00100000) * 1;
_maskNega = (_maskPosi == 0) * 1;
double _rxa = ((4.00000*{pca_}*{F_}*({F_}/({R_}*temp))*(_csm*exp(_za)-0.341000*{cao_}))/(2.00000*({F_}/({R_}*temp)))) * _maskPosi +
         ((4.00000*{pca_}*_v*{F_}*({F_}/({R_}*temp))*(_csm*exp(_za)-0.341000*{cao_}))/(exp(_za) - 1.00000)) * _maskNega;
double _rxa2 = fabs(_rxa);
double _jca = _gca*_po*_rxa;

// Eautions for Ca cycling
_maskPosi = (jsr>50.0000) * 1;
_maskPosi *= (jsr<{cstar_}) * 1;
_maskNega = (_maskPosi == 0) * 1;
_NegaPosi = _maskNega * (jsr >={cstar_}) * 1;
_NegaNega = _maskNega * (_NegaPosi == 0) * 1;
double _Qr0 = ((jsr - 50.0000)/1.00000) * _maskPosi + ({av_}*jsr-{bv_}) * _NegaPosi + 0 * _NegaNega;
double _Qr = (nsr*_Qr0)/{cstar_};
double _sparkV = exp(-{ay_}*(_v+30.0000))/(1.00000+exp(-{ay_}*(_v+30.0000)));
double _sparkrate = ({gryr_}/1.00000)*_po*_rxa2*_sparkV;

// calculate currents
_xina = _gna*h*j*m*m*m*(_v-_ena);
_xiks = _gks*_gksx*xs1*xs2*(_v-_eks);
_xik1 = _gkix*pow(({ko_}/5.40000), 1.0/2)*_xkin*(_v-_ek);
_xikr = _gkr*pow(({ko_}/5.40000), 1.0/2)*xr*_rg*(_v-_ek);
_xiNaK = (((_gNaK*_fNaK*nai)/(nai+{xkmnai_}))*{ko_})/({ko_}+{xkmko_});
_xitos = _gtos*xtos*(ytos+0.500000*_rsinf)*(_v-_ek);
_xitof = _gtof*xtof*ytof*(_v-_ek);
_xito = _xitos+_xitof;
_xica =  2.00000*{wca_}*_jca;
_xiNaCa = {wca_}*_jNaCa;
double _xirp = (((_po*_Qr*_rxa2*{gbarsr_})/1.00000)*exp(-{ax_}*(_v+30.0000)))/(1.00000+exp(-{ax_}*(_v+30.0000)));
double _xicap = _po*{gdyad_}*_rxa2;
double _xiryr = _xirp+_xicap;

// calculate difference
double _djsr = (nsr - jsr)/{taua_};
double _dh = (_ah*(1.00000-h) - _bh*h)*1.4;
double _dj = (_aj*(1.00000-j) - _bj*j)*1.4;
double _dm = _am*(1.00000-m) - _bm*m;
double _dxs1 = (_xs1ss - xs1)/_tauxs1;
double _dxs2 = (_xs2ss - xs2)/_tauxs2;
double _dxr = (_xkrinf - xr)/_taukr;
double _dc2 = ((_beta*c1+_k5*xi2ca+_k5t*xi2ba) - (_k6+_k6t+_alpha)*c2) * pow({Q10OTHER_}, (temp-{temp_})/10.0);
double _dxi2ca = ((_k3*xi1ca+_k6*c2) - (_k5+_k4)*xi2ca) * pow({Q10OTHER_}, (temp-{temp_})/10.0);
double _dxi2ba = ((_k3t*xi1ba+_k6t*c2) - (_k5t+_k4t)*xi2ba) * pow({Q10OTHER_}, (temp-{temp_})/10.0);
double _dc1 = ((_alpha*c2+{k2_}*xi1ca+{k2t_}*xi1ba+{r2_}*_po) - (_beta+{r1_}+{k1t_}+_k1)*c1) * pow({Q10OTHER_}, (temp-{temp_})/10.0);
double _dxi1ca = ((_k1*c1+_k4*xi2ca+_s1*_po) - (_k3+{k2_}+_s2)*xi1ca) * pow({Q10OTHER_}, (temp-{temp_})/10.0);
double _dxi1ba = (({k1t_}*c1+_k4t*xi2ba+ {s1t_}*_po) - (_k3t+{k2t_}+{s2t_})*xi1ba) * pow({Q10OTHER_}, (temp-{temp_})/10.0);
double _dxtos = (_xtosinf - xtos)/_txs;
double _dxtof = (_xtofinf - xtof)/_txf;
double _dytos = (_ytosinf - ytos)/_tys;
double _dytof = (_ytofinf - ytof)/_tyf;
double _dcai = _dciib*(((_jd-_jup)+_jleak)-_xbi);
double _dtropi = _xbi;
double _dtrops = _xbs;
double _dnsr = -xir+_jup-_jleak;
double _dxir = _sparkrate*_Qr - (xir*(1.00000-({taur_}*_dnsr)/nsr))/{taur_};
double _ddyad = _xiryr - (dyad-submem)/{taups_};
double _dsubmem = _dcsib*(50.0000*(((xir-_jd)-_jca)+_jNaCa)-_xbs);
double _dnai = -(_xina+3.00000*_xiNaK+3.00000*_xiNaCa)/({wca_}*1000.00);
_it = _xina+_xik1+_xikr+_xiks+_xito+_xiNaCa+_xica+_xiNaK+st;

// rewrite status
_jsr = jsr + _djsr * dt;
_h = h + _dh * dt;
_j = j + _dj * dt;
_m = m + _dm * dt;
_xs1 = xs1 + _dxs1 * dt;
_xs2 = xs2 + _dxs2 * dt;
_xr = xr + _dxr * dt;
_c2 = c2 + _dc2 * dt;
_xi2ca = xi2ca + _dxi2ca * dt;
_xi2ba = xi2ba + _dxi2ba * dt;
_c1 = c1 + _dc1 * dt;
_xi1ca = xi1ca + _dxi1ca * dt;
_xi1ba = xi1ba + _dxi1ba * dt;
_xtos = xtos + _dxtos * dt;
_xtof = xtof + _dxtof * dt;
_ytos = ytos + _dytos * dt;
_ytof = ytof + _dytof * dt;
_cai = cai + _dcai * dt;
_tropi = tropi + _dtropi * dt;
_trops = trops + _dtrops * dt;
_nsr = nsr + _dnsr * dt;
_xir = xir + _dxir * dt;
_dyad = dyad + _ddyad * dt;
_submem = submem + _dsubmem * dt;
_nai = nai + _dnai * dt;
_v -= _it * dt;

_dt = dt;
_st = st;
_temp = temp;

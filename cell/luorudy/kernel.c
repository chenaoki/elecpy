int _maskPosi;
int _maskNega;
int _maskPosiPosi;
int _maskPosiNega;
int _maskA;
int _maskA_bar;
int _maskB;
int _maskB_bar;

_v = v;

// voltage trapping to avoid error
_maskPosi =  ( _v < 0.00001 ) * 1;
_maskPosi *= ( _v >= 0 ) * 1;
_maskNega = ( _maskPosi == 0 ) * 1;
_v = _maskPosi * 0.00001 + _maskNega * _v;
_maskPosi =  ( _v > -0.00001 ) * 1;
_maskPosi *= ( _v <= 0 ) * 1;
_maskNega = ( _maskPosi == 0 ) * 1;
_v = _maskPosi * -0.00001 + _maskNega * _v;

_maskPosi =  ( _v < -30.0+0.001 ) * 1;
_maskPosi *= ( _v >= -30.0 ) * 1;
_maskNega = ( _maskPosi == 0 ) * 1;
_v = _maskPosi * (-30.0+0.001)+ _maskNega * _v;
_maskPosi =  ( _v > -30.0 -0.001 ) * 1;
_maskPosi *= ( _v <= -30.0 ) * 1;
_maskNega = ( _maskPosi == 0 ) * 1;
_v = _maskPosi * (-30.0-0.001) + _maskNega * _v;

// comp_ina ()
double _shift_inf_m = {SINFM_}*(temp-{temp_});
double _shift_inf_hj = {SINFHJ_}*(temp-{temp_});

double _ena = (({R_}*temp)/{frdy_})*log({nao_}/nai);
_maskPosi = (_v >= -95.0 + _shift_inf_m) * 1;
_maskNega = (_maskPosi == 0) * 1;
double _lv = (_v - _shift_inf_m)*_maskPosi + (-95.0 - _shift_inf_m)*_maskNega;
double _am = (0.32*(_lv+47.13)/(1-exp(-0.1*(_lv+47.13))))*pow({Q10TAUMHJ_}, (temp-{temp_})/10.0);
double _bm = (0.08*exp(-_lv/11))*pow({Q10TAUMHJ_}, (temp-{temp_})/10.0);

_maskPosi = (_v < -40+_shift_inf_hj)*1;
_maskNega = (_maskPosi == 0)*1;
_maskPosiPosi = (_v >= -135.0+_shift_inf_hj)*1;
_maskPosiNega = (_maskPosiPosi == 0)*1;
double _ah = (0.135*exp(-(80+((_v-_shift_inf_hj)*_maskPosiPosi-135.0*_maskPosiNega))/6.8)*_maskPosi + 0)*pow({Q10TAUMHJ_}, (temp-{temp_})/10.0);
double _bh = ((3.56*exp(0.079*((_v-_shift_inf_hj)*_maskPosiPosi-135.0*_maskPosiNega))+310000*exp(0.35*((_v-_shift_inf_hj)*_maskPosiPosi-135.0*_maskPosiNega)))*_maskPosi + (1/(0.13*(1+exp(-(_v-_shift_inf_hj+10.66)/11.1)))) * _maskNega)*pow({Q10TAUMHJ_}, (temp-{temp_})/10.0);

_maskPosiPosi = (_v >= -240.0+_shift_inf_hj)*1;
_maskPosiNega = (_maskPosiPosi == 0)*1;
double _aj = (((-127140*exp(0.2444*((_v-_shift_inf_hj)*_maskPosiPosi-240.0*_maskPosiNega))-0.00003474*exp(-0.04391*((_v-_shift_inf_hj)*_maskPosiPosi-240.0*_maskPosiNega)))*((((_v-_shift_inf_hj)*_maskPosiPosi-240.0*_maskPosiNega)+37.78)/(1+exp(0.311*(((_v-_shift_inf_hj)*_maskPosiPosi-240.0*_maskPosiNega)+79.23))))) * _maskPosi + 0)*pow({Q10TAUMHJ_}, (temp-{temp_})/10.0);
double _bj = (((0.1212*exp(-0.01052*((_v-_shift_inf_hj)*_maskPosiPosi-240.0*_maskPosiNega)))/(1+exp(-0.1378*(((_v-_shift_inf_hj)*_maskPosiPosi-240.0*_maskPosiNega)+40.14))))*_maskPosi + ((0.3*exp(-0.0000002535*(_v-_shift_inf_hj)))/(1+exp(-0.1*((_v-_shift_inf_hj)+32))))*_maskNega)*pow({Q10TAUMHJ_}, (temp-{temp_})/10.0);

double _mtau = 1 / (_am+_bm);
double _htau = 1 / (_ah+_bh);
double _jtau = 1 / (_aj+_bj);
double _mss  = _am * _mtau;
double _hss  = _ah * _htau;
double _jss  = _aj * _jtau;
_m = _mss - (_mss-m)*exp(-dt/_mtau);
_h = _hss - (_hss-h)*exp(-dt/_htau);
_j = _jss - (_jss-j)*exp(-dt/_jtau);
double _gna = {gna_}*pow({Q10NA_}, (temp-{temp_})/10.0);
double _ina = _gna*_m*_m*_m*_h*_j*(_v-_ena);

// comp_ltypesc
//   +- comp_rates
double _alphaequation = 0.925*exp(_v/30);
double _betaequation = 0.39*exp(-_v/40);
double _alpha0 = 4*_alphaequation;
double _alpha1 = 3*_alphaequation;
double _alpha2 = 2*_alphaequation;
double _alpha3 = 1*_alphaequation;
double _beta0 = 1*_betaequation;
double _beta1 = 2*_betaequation;
double _beta2 = 3*_betaequation;
double _beta3 = 4*_betaequation;
double _gammmaf = 0.245*exp(_v/10);
double _gammmas = 0.005*exp(-_v/40);
double _phif = 0.02*exp(_v/500);
double _phis = 0.03*exp(-_v/280);
double _lambdaf = 0.035*exp(-_v/300);
double _lambdas = 0.0011*exp(_v/500);
double _omgegaf = (_lambdaf*_beta3*_gammmaf)/(_phif*_alpha3);
double _omgegas = (_lambdas*_beta3*_gammmas)/(_phis*_alpha3);
double _omgegasf = (_lambdas*_phif)/(_lambdaf);
double _omgegafs = _phis;
double _caon = 4/(1+1/caiss);
double _caoff = 0.01;
double _mode0on = 0.0008;
double _mode0off = 0.00018;

//   +- comp_ltypeca
double _Czero1 = ltypeCzero;
double _Cone1 = ltypeCone;
double _Ctwo1 = ltypeCtwo;
double _Cthree1 = ltypeCthree;
double _IVf1 = ltypeIVf;
double _IVs1 = ltypeIVs;
double _O1 = 1-(_Czero1+_Cone1+_Ctwo1+_Cthree1+_IVf1+_IVs1);
double _dltypeCzero = dt*(_Cone1*_beta0-_Czero1*(_alpha0));
double _dltypeCone = dt*(_Ctwo1*_beta1+_Czero1*_alpha0-_Cone1*(_alpha1+_beta0));
double _dltypeCtwo = dt*(_Cone1*_alpha1+_Cthree1*_beta2-_Ctwo1*(_beta1+_alpha2));
double _dltypeCthree = dt*(_O1*_beta3+_IVf1*_omgegaf+_IVs1*_omgegas+_Ctwo1*_alpha2-_Cthree1*(_gammmaf+_gammmas+_alpha3+_beta2));
double _dltypeO = dt*(_IVf1*_lambdaf+_IVs1*_lambdas+_Cthree1*_alpha3-_O1*(_phif+_phis+_beta3));
double _dltypeIVf = dt*(_O1*_phif+_Cthree1*_gammmaf+_IVs1*_omgegasf-_IVf1*(_omgegaf+_lambdaf+_omgegafs));
double _dltypeIVs = dt*(_O1*_phis+_Cthree1*_gammmas+_IVf1*_omgegafs-_IVs1*(_omgegas+_lambdas+_omgegasf));
double _Czero2 = ltypeCzero+_dltypeCzero;
double _Cone2 = ltypeCone+_dltypeCone;
double _Ctwo2 = ltypeCtwo+_dltypeCtwo;
double _Cthree2 = ltypeCthree+_dltypeCthree;
double _IVf2 = ltypeIVf+_dltypeIVf;
double _IVs2 = ltypeIVs+_dltypeIVs;
_ltypeCzero = _Czero2;
_ltypeCone = _Cone2;
_ltypeCtwo = _Ctwo2;
_ltypeCthree = _Cthree2;
_ltypeIVf = _IVf2;
_ltypeIVs = _IVs2;
_ltypeO = 1-(_Czero2+_Cone2+_Ctwo2+_Cthree2+_IVf2+_IVs2);

//   +- comp_Cacon
double _buffersss = 1/(1+({bsrbar_}*{kmbsr_})/(({kmbsr_}+caiss)*({kmbsr_}+caiss))+({bslbar_}*{kmbsl_})/(({kmbsl_}+caiss)*({kmbsl_}+caiss)));
double _idiffca = (caiss-cai)/({taudiffss_}/pow({Q10OTHERS_}, (temp-{temp_})/10.0));
double _dcaresspace = dt*_buffersss*(-(ilca-2*inacass)*{acap_}/({volrss_}*{frdy_}*2)+(irel*({vjsr_}/{volrss_}))-_idiffca);
_caiss = caiss+_dcaresspace;

double _idiffna = (naiss-nai)/({taudiffss_}/pow({Q10OTHERS_}, (temp-{temp_})/10.0));
double _dnaresspace = dt*(-(3*inacass+ilcana)*{acap_}/({volrss_}*{frdy_})-_idiffna);
_naiss = naiss+_dnaresspace;

double _idiffk = (kiss-ki)/({taudiffss_}/pow({Q10OTHERS_}, (temp-{temp_})/10.0));
double _dkresspace = dt*(-ilcak*{acap_}/({volrss_}*{frdy_})-_idiffk);
_kiss = kiss+_dkresspace;

// comp_ltypesc
double _ibarca = {pca_}*{zca_}*{zca_}*((_v*{frdy_}*{frdy_})/({R_}*temp))*(({gacai_}*_caiss*exp(({zca_}*(_v-0)*{frdy_})/({R_}*temp))-{gacao_}*{cao_})/(exp(({zca_}*(_v-0)*{frdy_})/({R_}*temp))-1));
double _ibarna = {pna_}*{zna_}*{zna_}*((_v*{frdy_}*{frdy_})/({R_}*temp))*(({ganai_}*_naiss*exp(({zna_}*_v*{frdy_})/({R_}*temp))-{ganao_}*{nao_})/(exp(({zna_}*_v*{frdy_})/({R_}*temp))-1));
double _ibark = {pk_}*{zk_}*{zk_}*((_v*{frdy_}*{frdy_})/({R_}*temp))*(({gaki_}*_kiss*exp(({zk_}*_v*{frdy_})/({R_}*temp))-{gako_}*{ko_})/(exp(({zk_}*_v*{frdy_})/({R_}*temp))-1));
_fcasc = _caoff/(_caon+_caoff)-((_caoff/(_caon+_caoff))-fcasc)*exp(-dt/((1/(_caon+_caoff))/pow({Q10OTHERS_}, (temp-{temp_})/10.0)));
_fmode0 = _mode0off/(_mode0on+_mode0off)-((_mode0off/(_mode0on+_mode0off))-fmode0)*exp(-dt/(1/(_mode0on+_mode0off))/pow({Q10OTHERS_}, (temp-{temp_})/10.0));
_ilca = _ibarca*_fcasc*_fmode0*ltypeO*pow({Q10CAL_}, (temp-{temp_})/10.0);
_ilcana = _ibarna*_fcasc*_fmode0*ltypeO*pow({Q10CAL_}, (temp-{temp_})/10.0);
_ilcak = _ibark*_fcasc*_fmode0*ltypeO*pow({Q10CAL_}, (temp-{temp_})/10.0);

// comp_icat
double _bss = 1/(1+exp(-(_v+14)/10.8));
double _taub = (3.7+6.1/(1+exp((_v+25)/4.5)))/pow({Q10TAUBG_}, (temp-{temp_})/10.0);
double _gss = 1/(1+exp((_v+60)/5.6));
_maskPosi = (_v <= 0)*1;
_maskNega = (_v > 0)*1;
double _taug = ((-0.875*_v+12)*_maskPosi + 12*_maskNega)/pow({Q10TAUBG_}, (temp-{temp_})/10.0);
_b = _bss - ( _bss - b )*exp(-dt/_taub);
_g = _gss - ( _gss - g )*exp(-dt/_taug);

double _eca = ({R_}*temp/(2*{frdy_}))*log({cao_}/cai);
double _gcat = {gcat_}*pow({Q10CAT_}, (temp-{temp_})/10.0);
double _icat = _gcat*_b*_b*_g*(_v-_eca);

// comp_ikr
double _shift_inf_xr = -{SINFXR_}*(temp-{temp_});
double _ekr = (({R_}*temp)/{frdy_})*log({ko_}/ki);
_maskPosi = (_v >= -5000.0+_shift_inf_xr)*1;
_maskNega = (_maskPosi == 0)*1;
double _xrss = (1/(1+exp(-(_v+21.5)/7.5)))*_maskPosi + 0;

_maskA = (_v >= -5000.0+_shift_inf_xr)*1;
_maskA_bar = (_maskA == 0) * 1;
_maskB = (_v <= 4500.0+_shift_inf_xr)*1;
_maskB_bar = (_maskB == 0) * 1;
double _tauxr = ((1*(1/(0.00138*(_v-_shift_inf_xr+14.2)/(1-exp(-0.123*(_v-_shift_inf_xr+14.2)))+0.00061*(_v-_shift_inf_xr+38.9)/(exp(0.145*(_v-_shift_inf_xr+38.9))-1))))*_maskA*_maskB+0.16*_maskA*_maskB_bar+0.3*_maskA_bar*_maskB)/pow({Q10TAUXR_}, (temp-{temp_})/10.0);
_xr = _xrss-(_xrss-xr)*exp(-dt/_tauxr);
double _r = 1/(1+exp((_v+9)/18.4));
double _gkr = {gkr_}*pow({Q10KR_}, (temp-{temp_})/10.0);
double _ikr = _gkr*_xr*_r*(_v-_ekr);

// comp_iks
double _gks = (0.3031*(1+0.6/(1+pow((0.000038/cai),1.4))))*pow({Q10KS_}, (temp-{temp_})/10.0);
double _eks = (({R_}*temp)/{frdy_})*log(({ko_}+{prnak_}*{nao_})/(ki+{prnak_}*nai));
double _xs1ss = 1/(1+exp(-(_v-1.5)/16.7));
double _xs2ss = _xs1ss;
double _tauxs1 = (1/(0.0000719*(_v+30)/(1-exp(-0.148*(_v+30)))+0.000131*(_v+30)/(exp(0.0687*(_v+30))-1)))/pow({Q10TAUXS_}, (temp-{temp_})/10.0);
double _tauxs2 = 4*_tauxs1;

_xs1 = _xs1ss-(_xs1ss-xs1)*exp(-dt/_tauxs1);
_xs2 = _xs2ss-(_xs2ss-xs2)*exp(-dt/_tauxs2);
double _iks = _gks*_xs1*_xs2*(_v-_eks);

// comp_iki
double _eki = (({R_}*temp)/{frdy_})*log({ko_}/ki);
double _aki = 1.02/(1+exp(0.2385*(_v-_eki-59.215)));
double _bki = (0.49124*exp(0.08032*(_v-_eki+5.476))+exp(0.06175*(_v-_eki-594.31)))/(1+exp(-0.5143*(_v-_eki+4.753)));

double _kin = _aki/(_aki+_bki);
double _gki = {gki_}*pow({Q10K1_}, (temp-{temp_})/10.0);
double _iki = _gki*_kin*(_v-_eki);

// comp_ikp
double _ekp = _eki;
double _kp = 1/(1+exp((7.488-_v)/5.98));
double _gkp = {gkp_}*pow({Q10K1_}, (temp-{temp_})/10.0);
double _ikp = _gkp*_kp*(_v-_ekp);

// comp_inaca
// Calculates Na-Ca Exchanger Current
double _inaca = (0.8*{c1_}*exp(({gammanaca_}-1)*_v*{frdy_}/({R_}*temp))*((exp(_v*{frdy_}/({R_}*temp))*nai*nai*nai*{cao_} -1.5*{nao_}*{nao_}*{nao_}*cai)/(1+{c2_}*exp(({gammanaca_}-1)*_v*{frdy_}/({R_}*temp))*(exp(_v*{frdy_}/({R_}*temp))*nai*nai*nai*{cao_}+{c3_}*{nao_}*{nao_}*{nao_}*cai))))*pow({Q10NACA_}, (temp-{temp_})/10.0);
_inacass = (0.2*{c1_}*exp(({gammanaca_}-1)*_v*{frdy_}/({R_}*temp))*((exp(_v*{frdy_}/({R_}*temp))*_naiss*_naiss*_naiss*{cao_}-1.5*{nao_}*{nao_}*{nao_}*_caiss)/(1+{c2_}*exp(({gammanaca_}-1)*_v*{frdy_}/({R_}*temp))*(exp(_v*{frdy_}/({R_}*temp))*_naiss*_naiss*_naiss*{cao_}+{c3_}*{nao_}*{nao_}*{nao_}*_caiss))))*pow({Q10NACA_}, (temp-{temp_})/10.0);

// comp_inak
double _fnak = 1/(1+0.1245*exp((-0.1*_v*{frdy_})/({R_}*temp))+0.0365*{sigma_}*exp((-_v*{frdy_})/({R_}*temp)));
double _inak = ({ibarnak_}*_fnak*(1/(1+pow({kmnai_}/nai,2.0)))*({ko_}/({ko_}+{kmko_})))*pow({Q10NAK_}, (temp-{temp_})/10.0);

// comp_ipca
double _ipca = ({ibarpca_}*(cai/(0.0005+cai))/(1+exp((-cai+0.00012)/0.00001)))*pow({Q10K1_}, (temp-{temp_})/10.0);

// comp_icab
double _ecan = (({R_}*temp)/(2*{frdy_}))*log({cao_}/cai);
double _icab = ({gcab_}*(_v-_ecan))*pow({Q10K1_}, (temp-{temp_})/10.0);

//def comp_inab
double _inab = ({gnab_}*(_v-_ena))*pow({Q10K1_}, (temp-{temp_})/10.0);

//def comp_it
double _naoint = _ina+_inab+_ilcana+3*_inak+3*_inaca+3*_inacass;
double _koint = _ikr+_iks+_iki+_ikp+_ilcak-2*_inak;
double _caiont = _ilca+_icab+_ipca+_icat-2*_inaca-2*_inacass;
_it = _naoint+_koint+_caiont+st;

//def conc_nai
double _dnai = -dt*(((_naoint-3*_inacass-_ilcana)*{acap_})/({vmyo_}*{zna_}*{frdy_})-_idiffna*({volrss_}/{vmyo_}));
_nai = _dnai + nai;

//def conc_ki
double _dki = -dt*(((_koint-_ilcak+st)*{acap_})/({vmyo_}*{zk_}*{frdy_})-_idiffk*({volrss_}/{vmyo_}));
_ki = _dki + ki;

//def conc_nsr
double _ileak = {kleak_}*nsr*pow({Q10OTHERS_}, (temp-{temp_})/10.0);
double _iup = {iupbar_}*cai/(cai+{kmup_})*pow({Q10OTHERS_}, (temp-{temp_})/10.0);
double _dnsr = dt*(_iup-_ileak-itr*{vjsr_}/{vnsr_});
_nsr = nsr+_dnsr;

//def comp_ryrrates
double C1_C2 = 1750*_caiss;
double C2_C3 = 5600*_caiss;
double C3_C4 = 5600*_caiss;
double C4_O1 = 5600*_caiss;
double I1_I2 = 1750*_caiss;
double I2_I3 = 5600*_caiss;
double I3_I4 = 5600*_caiss;
double I4_I5 = 5600*_caiss;
double C2_C1 = 5;
double C3_C2 = 2.625;
double C4_C3 = 1;
double O1_C4 = 6.25;
double I2_I1 = 5;
double I3_I2 = 2.625;
double I4_I3 = 1;
double I5_I4 = 6.25;
double C1_I1 = 0.4*_caiss;
double C2_I2 = 1.2*_caiss;
double C3_I3 = 2.8*_caiss;
double C4_I4 = 5.2*_caiss;
double O1_I5 = 8.4*_caiss;
double I1_C1 = 0.01/(1+pow(({csqnbar_}*0.75/csqn),9.0));
double I2_C2 = 0.001/(1+pow(({csqnbar_}*0.75/csqn),9.0));
double I3_C3 = 0.0001/(1+pow(({csqnbar_}*0.75/csqn),9.0));
double I4_C4 = 0.00001/(1+pow(({csqnbar_}*0.75/csqn),9.0));
double I5_O1 = 0.000001/(1+pow(({csqnbar_}*0.75/csqn),9.0));

//def comp_ryr
double _ryrCone1 = ryrCone;
double _ryrCtwo1 = ryrCtwo;
double _ryrCthree1 = ryrCthree;
double _ryrCfour1 = ryrCfour;
double _ryrOone1 = ryrOone;
double _ryrIone1 = ryrIone;
double _ryrItwo1 = ryrItwo;
double _ryrIthree1 = ryrIthree;
double _ryrIfour1 = ryrIfour;
double _ryrIfive1 = 1-(_ryrCone1+_ryrCtwo1+_ryrCthree1+_ryrCfour1+_ryrOone1+_ryrIone1+_ryrItwo1+_ryrIthree1+_ryrIfour1);
double _dryrCone = dt*(C2_C1*_ryrCtwo1+I1_C1*_ryrIone1-_ryrCone1*(C1_C2+C1_I1));
double _dryrCtwo = dt*(C1_C2*_ryrCone1+C3_C2*_ryrCthree1+I2_C2*_ryrItwo1-_ryrCtwo1*(C2_C1+C2_C3+C2_I2));
double _dryrCthree = dt*(C2_C3*_ryrCtwo1+C4_C3*_ryrCfour1+I3_C3*_ryrIthree1-_ryrCthree1*(C3_C2+C3_C4+C3_I3));
double _dryrCfour = dt*(O1_C4*_ryrOone1+C3_C4*_ryrCthree1+I4_C4*_ryrIfour1-_ryrCfour1*(C4_C3+C4_O1+C4_I4));
double _dryrOone = dt*(C4_O1*_ryrCfour1+I5_O1*_ryrIfive1-_ryrOone1*(O1_I5+O1_C4));
double _dryrIone = dt*(C1_I1*_ryrCone1+I2_I1*_ryrItwo1-_ryrIone1*(I1_C1+I1_I2));
double _dryrItwo = dt*(C2_I2*_ryrCtwo1+I1_I2*_ryrIone1+I3_I2*_ryrIthree1-_ryrItwo1*(I2_C2+I2_I1+I2_I3));
double _dryrIthree = dt*(C3_I3*_ryrCthree1+I2_I3*_ryrItwo1+I4_I3*_ryrIfour1-_ryrIthree1*(I3_C3+I3_I2+I3_I4));
double _dryrIfour = dt*(C4_I4*_ryrCfour1+I3_I4*_ryrIthree1+I5_I4*_ryrIfive1-_ryrIfour1*(I4_C4+I4_I3+I4_I5));
double _ryrCone2 = ryrCone     +_dryrCone;
double _ryrCtwo2 = ryrCtwo     +_dryrCtwo;
double _ryrCthree2 = ryrCthree +_dryrCthree;
double _ryrCfour2 = ryrCfour   +_dryrCfour;
double _ryrOone2 = ryrOone     +_dryrOone;
double _ryrIone2 = ryrIone     +_dryrIone;
double _ryrItwo2 = ryrItwo     +_dryrItwo;
double _ryrIthree2 = ryrIthree +_dryrIthree;
double _ryrIfour2 = ryrIfour   +_dryrIfour;
_ryrCone   = _ryrCone2;
_ryrCtwo   = _ryrCtwo2;
_ryrCthree = _ryrCthree2;
_ryrCfour  = _ryrCfour2;
_ryrOone   = _ryrOone2;
_ryrIone   = _ryrIone2;
_ryrItwo   = _ryrItwo2;
_ryrIthree = _ryrIthree2;
_ryrIfour  = _ryrIfour2;
_ryrIfive = 1-(_ryrCone2+_ryrCtwo2+_ryrCthree2+_ryrCfour2+_ryrOone2+_ryrIone2+_ryrItwo2+_ryrIthree2+_ryrIfour2);

//def conc_jsr
double _expArg = (5.3-jsr)/0.001;
_maskPosi = (_expArg <= 709)*1;
_maskNega = (_maskPosi == 0)*1;
double _xpArg = _expArg*_maskPosi + 709*_maskNega;
double _sponrelss = 25/(1+exp(_expArg));
_sponrel = _sponrelss-(_sponrelss-sponrel)*exp(-dt/({spontau_}/pow({Q10OTHERS_}, (temp-{temp_})/10.0)));
double _gradedrel = (1/(1+exp((20+_ilca)/6))-0.034445);
double _vgainofrel = (1+exp(((0.015*_ibarca)+1.25)/0.75));
double _vg = _gradedrel/_vgainofrel;
double _grel = 250*(_vg+_sponrel)*pow({Q10REL_}, (temp-{temp_})/10.0);
_irel = _grel*(ryrOone)*(jsr-_caiss);
_csqn = {csqnbar_}*(jsr/(jsr+{kmcsqn_}));

double _buffersjsr = 1/(1+({csqnbar_}*{kmcsqn_}/(pow({kmcsqn_}+jsr,2.0))));
double _djsr = dt*_buffersjsr*(itr-_irel);
_jsr = jsr+_djsr;

//def calc_itr
_itr = (_nsr-_jsr)/({tautr_}/pow({Q10OTHERS_}, (temp-{temp_})/10.0));

//def conc_cai
double _buffersmyo = 1/(1+({trpnbar_}*{kmtrpn_})/(({kmtrpn_}+cai)*({kmtrpn_}+cai))+({cmdnbar_}*{kmcmdn_})/(({kmcmdn_}+cai)*({kmcmdn_}+cai)));
double _dcai = -dt*_buffersmyo*((((_caiont-_ilca+2*_inacass)*{acap_})/({vmyo_}*{zca_}*{frdy_}))+((_iup-_ileak)*({vnsr_}/{vmyo_}))-_idiffca*({volrss_}/{vmyo_}));
_cai = cai+_dcai;

// Membrane voltage update
_v -= _it*dt;

// unmodified return variables
_dt = dt;
_st = st;
_temp = temp;


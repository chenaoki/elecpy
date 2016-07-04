T _maskPosi; 
T _maskNega; 

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
T _ena = (({R_}*{temp_})/{frdy_})*log({nao_}/nai);
T _am = 0.32*(_v+47.13)/(1-exp(-0.1*(_v+47.13)));
T _bm = 0.08*exp(-_v/11);
_maskPosi = (_v < -40)*1;
_maskNega = (_v >=-40)*1;
T _ah = 0.135*exp((80+_v)/-6.8)*_maskPosi + 0;
T _bh = (3.56*exp(0.079*_v)+310000*exp(0.35*_v))*_maskPosi + (1/(0.13*(1+exp((_v+10.66)/-11.1)))) * _maskNega;
T _aj = ((-127140*exp(0.2444*_v)-0.00003474*exp(-0.04391*_v))*((_v+37.78)/(1+exp(0.311*(_v+79.23))))) * _maskPosi + 0;
T _bj = ((0.1212*exp(-0.01052*_v))/(1+exp(-0.1378*(_v+40.14))))*_maskPosi + ((0.3*exp(-0.0000002535*_v))/(1+exp(-0.1*(_v+32))))*_maskNega;
_h = _ah/(_ah+_bh)-( (_ah/(_ah+_bh)) - h )*exp(-dt/(1/(_ah+_bh)));
_j = _aj/(_aj+_bj)-( (_aj/(_aj+_bj)) - j )*exp(-dt/(1/(_aj+_bj)));
_m = _am/(_am+_bm)-( (_am/(_am+_bm)) - m )*exp(-dt/(1/(_am+_bm)));
T _ina = {gna_}*_m*_m*_m*_h*_j*(_v-_ena);

// comp_ltypesc
//   +- comp_rates
T _alphaequation;
T _betaequation;
T _alpha0;
T _alpha1;
T _alpha2;
T _alpha3;
T _beta0;
T _beta1;
T _beta2;
T _beta3;
T _gammmaf;
T _gammmas;
T _phif;
T _phis;
T _lambdaf;
T _lambdas;
T _omgegaf;
T _omgegas;
T _omgegasf;
T _omgegafs;
T _caon;
T _caoff;
T _mode0on;
T _mode0off;

//if({betaad_}==0)
_alphaequation = 0.925*exp(_v/30);
_betaequation = 0.39*exp(_v/-40);
_alpha0 = 4*_alphaequation;
_alpha1 = 3*_alphaequation;
_alpha2 = 2*_alphaequation;
_alpha3 = 1*_alphaequation;
_beta0 = 1*_betaequation;
_beta1 = 2*_betaequation;
_beta2 = 3*_betaequation;
_beta3 = 4*_betaequation;
_gammmaf = 0.245*exp(_v/10);
_gammmas = 0.005*exp(_v/-40);
_phif = 0.02*exp(_v/500);
_phis = 0.03*exp(_v/-280);
_lambdaf = 0.035*exp(_v/-300);
_lambdas = 0.0011*exp(_v/500);
_omgegaf = (_lambdaf*_beta3*_gammmaf)/(_phif*_alpha3);
_omgegas = (_lambdas*_beta3*_gammmas)/(_phis*_alpha3);
_omgegasf = (_lambdas*_phif)/(_lambdaf);
_omgegafs = _phis;
_caon = 4/(1+1/caiss);
_caoff = 0.01;
_mode0on = 0.0008;
_mode0off = 0.00018;

//if({betaad_}==1)
//  _alphaequation = 0.925*exp(_v/30);
//  _betaequation = 0.39*exp(_v/-40);
//  _alpha0 = 4*_alphaequation;
//  _alpha1 = 3*_alphaequation;
//  _alpha2 = 2*_alphaequation;
//  _alpha3 = 1*_alphaequation;
//  _beta0 = 1*_betaequation;
//  _beta1 = 2*_betaequation;
//  _beta2 = 3*_betaequation;
//  _beta3 = 4*_betaequation;
//  _gammmaf = 0.245*exp(_v/10);
//  _gammmas = 0.005*exp(_v/-40);
//  _phif = 0.02*exp(_v/500);
//  _phis = 0.014*exp(_v/-20);
//  _lambdaf = 0.015*exp(_v/-300);
//  _lambdas = 0.0022*exp(_v/500);
//  _omgegaf = (_lambdaf*_beta3*_gammmaf)/(_phif*_alpha3);
//  _omgegas = (_lambdas*_beta3*_gammmas)/(_phis*_alpha3);
//  _omgegasf = (_lambdas*_phif)/(_lambdaf);
//  _omgegafs = _phis;
//  _caon = 5/(1+1/caiss);
//  _caoff = 0.001;
//  _mode0on = 0.0004;
//  _mode0off = 0.0004545;

//   +- comp_ltypeca
T _Czero1 = ltypeCzero;
T _Cone1 = ltypeCone;
T _Ctwo1 = ltypeCtwo;
T _Cthree1 = ltypeCthree;
T _IVf1 = ltypeIVf;
T _IVs1 = ltypeIVs;
T _O1 = 1-(_Czero1+_Cone1+_Ctwo1+_Cthree1+_IVf1+_IVs1);
T _dltypeCzero = dt*(_Cone1*_beta0-_Czero1*(_alpha0));
T _dltypeCone = dt*(_Ctwo1*_beta1+_Czero1*_alpha0-_Cone1*(_alpha1+_beta0));
T _dltypeCtwo = dt*(_Cone1*_alpha1+_Cthree1*_beta2-_Ctwo1*(_beta1+_alpha2));
T _dltypeCthree = dt*(_O1*_beta3+_IVf1*_omgegaf+_IVs1*_omgegas+_Ctwo1*_alpha2-_Cthree1*(_gammmaf+_gammmas+_alpha3+_beta2));
T _dltypeO = dt*(_IVf1*_lambdaf+_IVs1*_lambdas+_Cthree1*_alpha3-_O1*(_phif+_phis+_beta3));
T _dltypeIVf = dt*(_O1*_phif+_Cthree1*_gammmaf+_IVs1*_omgegasf-_IVf1*(_omgegaf+_lambdaf+_omgegafs));
T _dltypeIVs = dt*(_O1*_phis+_Cthree1*_gammmas+_IVf1*_omgegafs-_IVs1*(_omgegas+_lambdas+_omgegasf));
T _Czero2 = ltypeCzero+_dltypeCzero;
T _Cone2 = ltypeCone+_dltypeCone;
T _Ctwo2 = ltypeCtwo+_dltypeCtwo;
T _Cthree2 = ltypeCthree+_dltypeCthree;
//O2 = ltypeO+_dltypeO;
T _IVf2 = ltypeIVf+_dltypeIVf;
T _IVs2 = ltypeIVs+_dltypeIVs;
_ltypeCzero = _Czero2;
_ltypeCone = _Cone2;
_ltypeCtwo = _Ctwo2;
_ltypeCthree = _Cthree2;
_ltypeIVf = _IVf2;
_ltypeIVs = _IVs2;
_ltypeO = 1-(_Czero2+_Cone2+_Ctwo2+_Cthree2+_IVf2+_IVs2);

//   +- comp_Cacon
//bsr = {bsrbar_}*(caiss/(caiss+{kmbsr_}));
//bsl = {bslbar_}*(caiss/(caiss+{kmbsl_}));
T _buffersss = 1/(1+({bsrbar_}*{kmbsr_})/(({kmbsr_}+caiss)*({kmbsr_}+caiss))+({bslbar_}*{kmbsl_})/(({kmbsl_}+caiss)*({kmbsl_}+caiss)));
T _idiffca = (caiss-cai)/{taudiffss_};
T _dcaresspace = dt*_buffersss*(-(ical-2*inacass)*{acap_}/({volrss_}*{frdy_}*2)+(irel*({vjsr_}/{volrss_}))-_idiffca);
_caiss = caiss+_dcaresspace;

T _idiffna = (naiss-nai)/{taudiffss_};
T _dnaresspace = dt*(-(3*inacass+ilcana)*{acap_}/({volrss_}*{frdy_})-_idiffna);
_naiss = naiss+_dnaresspace;

T _idiffk = (kiss-ki)/{taudiffss_};
T _dkresspace = dt*(-ilcak*{acap_}/({volrss_}*{frdy_})-_idiffk);
_kiss = kiss+_dkresspace;

// comp_ltypesc
T _ibarca = {pca_}*{zca_}*{zca_}*((_v*{frdy_}*{frdy_})/({R_}*{temp_}))*(({gacai_}*_caiss*exp(({zca_}*(_v-0)*{frdy_})/({R_}*{temp_}))-{gacao_}*{cao_})/(exp(({zca_}*(_v-0)*{frdy_})/({R_}*{temp_}))-1));
T _ibarna = {pna_}*{zna_}*{zna_}*((_v*{frdy_}*{frdy_})/({R_}*{temp_}))*(({ganai_}*_naiss*exp(({zna_}*_v*{frdy_})/({R_}*{temp_}))-{ganao_}*{nao_})/(exp(({zna_}*_v*{frdy_})/({R_}*{temp_}))-1));
T _ibark = {pk_}*{zk_}*{zk_}*((_v*{frdy_}*{frdy_})/({R_}*{temp_}))*(({gaki_}*_kiss*exp(({zk_}*_v*{frdy_})/({R_}*{temp_}))-{gako_}*{ko_})/(exp(({zk_}*_v*{frdy_})/({R_}*{temp_}))-1));
_fcasc = _caoff/(_caon+_caoff)-((_caoff/(_caon+_caoff))-fcasc)*exp(-dt/(1/(_caon+_caoff)));
_fmode0 = _mode0off/(_mode0on+_mode0off)-((_mode0off/(_mode0on+_mode0off))-fmode0)*exp(-dt/(1/(_mode0on+_mode0off)));
_ical = _ibarca*_fcasc*_fmode0*ltypeO;
_ilcana = _ibarna*_fcasc*_fmode0*ltypeO;
_ilcak = _ibark*_fcasc*_fmode0*ltypeO;

T _ilca = _ical;

// comp_icat
T _bss = 1/(1+exp(-(_v+14)/10.8));
T _taub = 3.7+6.1/(1+exp((_v+25)/4.5));
T _gss = 1/(1+exp((_v+60)/5.6));
_maskPosi = (_v <= 0)*1;
_maskNega = (_v > 0)*1;
T _taug = (-0.875*_v+12)*_maskPosi + 12*_maskNega;
_b = _bss - ( _bss - b )*exp(-dt/_taub);
_g = _gss - ( _gss - g )*exp(-dt/_taug);

T _eca = ({R_}*{temp_}/(2*{frdy_}))*log({cao_}/cai);
T _icat = {gcat_}*_b*_b*_g*(_v-_eca);

// comp_ikr
T _ekr = (({R_}*{temp_})/{frdy_})*log({ko_}/ki);
T _xrss = 1/(1+exp(-(_v+21.5)/7.5));
T _tauxr = 1*(1/(0.00138*(_v+14.2)/(1-exp(-0.123*(_v+14.2)))+0.00061*(_v+38.9)/(exp(0.145*(_v+38.9))-1)));
_xr = _xrss-(_xrss-xr)*exp(-dt/_tauxr);

T _r = 1/(1+exp((_v+9)/18.4));
T _ikr = {gkr_}*_xr*_r*(_v-_ekr);

// comp_iks
T _gks;
T _eks;
T _xs1ss;
T _xs2ss;
T _tauxs1;
T _tauxs2;
//if({betaad_}==0)
_gks = 0.3031*(1+0.6/(1+pow((0.000038/cai),1.4)));
_eks = (({R_}*{temp_})/{frdy_})*log(({ko_}+{prnak_}*{nao_})/(ki+{prnak_}*nai));
_xs1ss = 1/(1+exp(-(_v-1.5)/16.7));
_xs2ss = _xs1ss;
_tauxs1 = 1/(0.0000719*(_v+30)/(1-exp(-0.148*(_v+30)))+0.000131*(_v+30)/(exp(0.0687*(_v+30))-1));
_tauxs2 = 4*_tauxs1;

//if({betaad_}==1)
//_gks = 1.5*0.3031*(1+0.6/(1+pow((0.000038/cai),1.4)));
//_eks = (({R_}*{temp_})/{frdy_})*log(({ko_}+{prnak_}*{nao_})/(ki+{prnak_}*nai));
//_xs1ss = 1/(1+exp(-(_v-1.5+4)/16.7));
//_xs2ss = _xs1ss;
//_tauxs1 = 1/(0.0000719*((_v+4)+30)/(1-exp(-0.148*((_v+4)+30)))+0.000131*((_v+4)+30)/(exp(0.0687*((_v+4)+30))-1));
//_tauxs2 = 4*_tauxs1;

_xs1 = _xs1ss-(_xs1ss-xs1)*exp(-dt/_tauxs1);
_xs2 = _xs2ss-(_xs2ss-xs2)*exp(-dt/_tauxs2);
T _iks = _gks*_xs1*_xs2*(_v-_eks);

// comp_iki
T _eki = (({R_}*{temp_})/{frdy_})*log({ko_}/ki);
T _aki = 1.02/(1+exp(0.2385*(_v-_eki-59.215)));
//if({betaad_}==0)
T _bki = (0.49124*exp(0.08032*(_v-_eki+5.476))+exp(0.06175*(_v-_eki-594.31)))/(1+exp(-0.5143*(_v-_eki+4.753)));
//if({betaad_}==1)
//T _bki = 2+(0.49124*exp(0.08032*(_v-_eki+5.476))+exp(0.06175*(_v-_eki-594.31)))/(1+exp(-0.5143*(_v-_eki+4.753)));

T _kin = _aki/(_aki+_bki);
T _iki = {gki_}*_kin*(_v-_eki);

// comp_ikp
T _ekp = _eki;
T _kp = 1/(1+exp((7.488-_v)/5.98));
T _ikp = {gkp_}*_kp*(_v-_ekp);

// comp_inaca
// Calculates Na-Ca Exchanger Current
T _inaca = 0.8*{c1_}*exp(({gammanaca_}-1)*_v*{frdy_}/({R_}*{temp_}))*((exp(_v*{frdy_}/({R_}*{temp_}))*nai*nai*nai*{cao_} -1.5*{nao_}*{nao_}*{nao_}*cai)/(1+{c2_}*exp(({gammanaca_}-1)*_v*{frdy_}/({R_}*{temp_}))*(exp(_v*{frdy_}/({R_}*{temp_}))*nai*nai*nai*{cao_}+{c3_}*{nao_}*{nao_}*{nao_}*cai)));
_inacass = 0.2*{c1_}*exp(({gammanaca_}-1)*_v*{frdy_}/({R_}*{temp_}))*((exp(_v*{frdy_}/({R_}*{temp_}))*_naiss*_naiss*_naiss*{cao_}-1.5*{nao_}*{nao_}*{nao_}*_caiss)/(1+{c2_}*exp(({gammanaca_}-1)*_v*{frdy_}/({R_}*{temp_}))*(exp(_v*{frdy_}/({R_}*{temp_}))*_naiss*_naiss*_naiss*{cao_}+{c3_}*{nao_}*{nao_}*{nao_}*_caiss)));

// comp_inak
T _fnak = 1/(1+0.1245*exp((-0.1*_v*{frdy_})/({R_}*{temp_}))+0.0365*{sigma_}*exp((-_v*{frdy_})/({R_}*{temp_})));
T _inak;
//if({betaad_}==0)
  _inak = {ibarnak_}*_fnak*(1/(1+pow({kmnai_}/nai,2)))*({ko_}/({ko_}+{kmko_}));
//if( {betaad_}==1)
//  _inak = 1.3*{ibarnak_}*_fnak*(1/(1+pow({kmnai_}/nai,2)))*({ko_}/({ko_}+{kmko_}));

// comp_ipca
T _ipca = {ibarpca_}*(cai/(0.0005+cai))/(1+exp((-cai+0.00012)/0.00001));

// comp_icab
T _ecan = (({R_}*{temp_})/(2*{frdy_}))*log({cao_}/cai);
T _icab = {gcab_}*(_v-_ecan);

//def comp_inab
T _inab = {gnab_}*(_v-_ena);

//def comp_it
T _naoint = _ina+_inab+_ilcana+3*_inak+3*_inaca+3*_inacass;
T _koint = _ikr+_iks+_iki+_ikp+_ilcak-2*_inak;
T _caiont = _ilca+_icab+_ipca+_icat-2*_inaca-2*_inacass;
_it = _naoint+_koint+_caiont+st;

//def conc_nai
T _dnai = -dt*(((_naoint-3*_inacass-_ilcana)*{acap_})/({vmyo_}*{zna_}*{frdy_})-_idiffna*({volrss_}/{vmyo_}));
_nai = _dnai + nai;

//def conc_ki
T _dki;
//if(stimtime>=0 && stimtime<0.5)
_dki = -dt*(((_koint-_ilcak+st)*{acap_})/({vmyo_}*{zk_}*{frdy_})-_idiffk*({volrss_}/{vmyo_}));
//else
//_dki = -dt*(((_koint-_ilcak)*{acap_})/({vmyo_}*{zk_}*{frdy_})-_idiffk*({volrss_}/{vmyo_}));
_ki = _dki + ki;

//def conc_nsr
T _ileak = {kleak_}*nsr;
//if({betaad_}==0)
T _iup = {iupbar_}*cai/(cai+{kmup_});
//if({betaad_}==1)
//T _iup = 1.2*{iupbar_}*cai/(cai+{kmup_});
T _dnsr = dt*(_iup-_ileak-itr*{vjsr_}/{vnsr_});
_nsr = nsr+_dnsr;

//def comp_ryrrates
T C1_C2 = 1750*_caiss;
T C2_C3 = 5600*_caiss;
T C3_C4 = 5600*_caiss;
T C4_O1 = 5600*_caiss;
T I1_I2 = 1750*_caiss;
T I2_I3 = 5600*_caiss;
T I3_I4 = 5600*_caiss;
T I4_I5 = 5600*_caiss;
T C2_C1 = 5;
T C3_C2 = 2.625;
T C4_C3 = 1;
T O1_C4 = 6.25;
T I2_I1 = 5;
T I3_I2 = 2.625;
T I4_I3 = 1;
T I5_I4 = 6.25;
T C1_I1 = 0.4*_caiss;
T C2_I2 = 1.2*_caiss;
T C3_I3 = 2.8*_caiss;
T C4_I4 = 5.2*_caiss;
T O1_I5 = 8.4*_caiss;
T I1_C1 = 0.01/(1+pow(({csqnbar_}*0.75/csqn),9));
T I2_C2 = 0.001/(1+pow(({csqnbar_}*0.75/csqn),9));
T I3_C3 = 0.0001/(1+pow(({csqnbar_}*0.75/csqn),9));
T I4_C4 = 0.00001/(1+pow(({csqnbar_}*0.75/csqn),9));
T I5_O1 = 0.000001/(1+pow(({csqnbar_}*0.75/csqn),9));

//def comp_ryr
T _ryrCone1 = ryrCone;
T _ryrCtwo1 = ryrCtwo;
T _ryrCthree1 = ryrCthree;
T _ryrCfour1 = ryrCfour;
T _ryrOone1 = ryrOone;
T _ryrIone1 = ryrIone;
T _ryrItwo1 = ryrItwo;
T _ryrIthree1 = ryrIthree;
T _ryrIfour1 = ryrIfour;
T _ryrIfive1 = 1-(_ryrCone1+_ryrCtwo1+_ryrCthree1+_ryrCfour1+_ryrOone1+_ryrIone1+_ryrItwo1+_ryrIthree1+_ryrIfour1);
T _dryrCone = dt*(C2_C1*_ryrCtwo1+I1_C1*_ryrIone1-_ryrCone1*(C1_C2+C1_I1));
T _dryrCtwo = dt*(C1_C2*_ryrCone1+C3_C2*_ryrCthree1+I2_C2*_ryrItwo1-_ryrCtwo1*(C2_C1+C2_C3+C2_I2));
T _dryrCthree = dt*(C2_C3*_ryrCtwo1+C4_C3*_ryrCfour1+I3_C3*_ryrIthree1-_ryrCthree1*(C3_C2+C3_C4+C3_I3));
T _dryrCfour = dt*(O1_C4*_ryrOone1+C3_C4*_ryrCthree1+I4_C4*_ryrIfour1-_ryrCfour1*(C4_C3+C4_O1+C4_I4));
T _dryrOone = dt*(C4_O1*_ryrCfour1+I5_O1*_ryrIfive1-_ryrOone1*(O1_I5+O1_C4));
T _dryrIone = dt*(C1_I1*_ryrCone1+I2_I1*_ryrItwo1-_ryrIone1*(I1_C1+I1_I2));
T _dryrItwo = dt*(C2_I2*_ryrCtwo1+I1_I2*_ryrIone1+I3_I2*_ryrIthree1-_ryrItwo1*(I2_C2+I2_I1+I2_I3));
T _dryrIthree = dt*(C3_I3*_ryrCthree1+I2_I3*_ryrItwo1+I4_I3*_ryrIfour1-_ryrIthree1*(I3_C3+I3_I2+I3_I4));
T _dryrIfour = dt*(C4_I4*_ryrCfour1+I3_I4*_ryrIthree1+I5_I4*_ryrIfive1-_ryrIfour1*(I4_C4+I4_I3+I4_I5));
//T _dryrIfive = dt*(O1_I5*_ryrOone1+I4_I5*_ryrIfour1-_ryrIfive1*(I5_O1+I5_I4));
T _ryrCone2 = ryrCone     +_dryrCone;
T _ryrCtwo2 = ryrCtwo     +_dryrCtwo;
T _ryrCthree2 = ryrCthree +_dryrCthree;
T _ryrCfour2 = ryrCfour   +_dryrCfour;
T _ryrOone2 = ryrOone     +_dryrOone;
T _ryrIone2 = ryrIone     +_dryrIone;
T _ryrItwo2 = ryrItwo     +_dryrItwo;
T _ryrIthree2 = ryrIthree +_dryrIthree;
T _ryrIfour2 = ryrIfour   +_dryrIfour;
//T _ryrIfive2 = ryrIfive   +_dryrIfive;
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
T _expArg = (5.3-jsr)/.001;
_maskPosi = (_expArg <= 709)*1;
_maskNega = (_expArg >  709)*1;
T _xpArg = _expArg*_maskPosi + 709*_maskNega;
T _sponrelss = 25/(1+exp(_expArg));
T _sponrel = _sponrelss-(_sponrelss-_sponrel)*exp(-dt/{spontau_});
T _gradedrel = (1/(1+exp((20+_ical)/6))-0.034445);
T _vgainofrel = (1+exp(((0.015*_ibarca)+1.25)/0.75));
T _vg = _gradedrel/_vgainofrel;
T _grel = 250*(_vg+_sponrel);
_irel = _grel*(ryrOone)*(jsr-_caiss);
_csqn = {csqnbar_}*(jsr/(jsr+{kmcsqn_}));

T _buffersjsr = 1/(1+({csqnbar_}*{kmcsqn_}/(pow({kmcsqn_}+jsr,2))));
T _djsr = dt*_buffersjsr*(itr-_irel);
_jsr = jsr+_djsr;

//def calc_itr
_itr = (_nsr-_jsr)/{tautr_};

//def conc_cai
//trpn = {trpnbar_}*(cai/(cai+{kmtrpn_}));
//cmdn = {cmdnbar_}*(cai/(cai+{kmcmdn_}));
T _buffersmyo = 1/(1+({trpnbar_}*{kmtrpn_})/(({kmtrpn_}+cai)*({kmtrpn_}+cai))+({cmdnbar_}*{kmcmdn_})/(({kmcmdn_}+cai)*({kmcmdn_}+cai)));
T _dcai = -dt*_buffersmyo*((((_caiont-_ilca+2*_inacass)*{acap_})/({vmyo_}*{zca_}*{frdy_}))+((_iup-_ileak)*({vnsr_}/{vmyo_}))-_idiffca*({volrss_}/{vmyo_}));
_cai = cai+_dcai;

// Membrane voltage update
_v -= _it*dt;

// unmodified return variables
_dt = dt;
_st = st;

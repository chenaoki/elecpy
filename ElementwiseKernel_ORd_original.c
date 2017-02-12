#define CELLTYPE_EPI

// compute reversal potentials
T ENa = {R_}*temp/{frdy_}*log({nao_}/nai);
T EK  = {R_}*temp/{frdy_}*log({ko_}/ki);
T EKs = {R_}*temp/{frdy_}*log(({ko_}+{prnak_}*{nao_})/(ki+{prnak_}*nai));

// CaMK valuables?
T CaMKb = {CaMKo_}*(1.0-CaMKt)/(1.0+{KmCaM_}/cass);
T CaMKa = CaMKb+CaMKt;
T vffrt = v*{frdy_}*{frdy_}/({R_}*temp);
T vfrt  = v*{frdy_}/({R_}*temp);

//Na ion current
T mss   = 1.0/(1.0+exp((-(v+39.57))/9.871));
T tm    = 1.0/(6.765*exp((v+11.64)/34.77)+8.552*exp(-(v+77.42)/5.955));
_m      = mss-(mss-m)*exp(-dt/tm);
T hss   = 1.0/(1+exp((v+82.90)/6.086));
T thf   = 1.0/(1.432e-5*exp(-(v+1.196)/6.285)+6.149*exp((v+0.5096)/20.27));
T ths   = 1.0/(0.009794*exp(-(v+17.95)/28.05)+0.3343*exp((v+5.730)/56.66));
T Ahf   = 0.99;
T Ahs   = 1.0-Ahf;
_hf     = hss-(hss-hf)*exp(-dt/thf);
_hs     = hss-(hss-hs)*exp(-dt/ths);
T h     = Ahf*_hf+Ahs*_hs;
T jss   = hss;
T tj    = 2.038+1.0/(0.02136*exp(-(v+100.6)/8.281)+0.3052*exp((v+0.9941)/38.45));
_j      = jss-(jss-j)*exp(-dt/tj);
T hssp  = 1.0/(1+exp((v+89.1)/6.086));
T thsp  = 3.0*ths;
_hsp    = hssp-(hssp-hsp)*exp(-dt/thsp);
T hp    = Ahf*_hf+Ahs*_hsp;
T tjp   = 1.46*tj;
_jp     = jss-(jss-jp)*exp(-dt/tjp);
T GNa   = 75; //conductance
T fINap = (1.0/(1.0+{KmCaMK_}/CaMKa));
T INa    = GNa*(v-ENa)*_m*_m*_m*((1.0-fINap)*h*_j+fINap*hp*_jp);

// NaL ion current
T mLss   = 1.0/(1.0+exp((-(v+42.85))/5.264));
T tmL    = tm;
_mL      = mLss-(mLss-mL)*exp(-dt/tmL);
T hLss   = 1.0/(1.0+exp((v+87.61)/7.488));
T thL    = 200.0;
_hL      = hLss-(hLss-hL)*exp(-dt/thL);
T hLssp  = 1.0/(1.0+exp((v+93.81)/7.488));
T thLp   = 3.0*thL;
_hLp     = hLssp-(hLssp-hLp)*exp(-dt/thLp);
T GNaL   = 0.0075;
// celltype (endo = 0, epi = 1, M = 2)
#ifdef CELLTYPE_EPI
GNaL    *= 0.6;
#endif
T fINaLp = (1.0/(1.0+{KmCaMK_}/CaMKa));
T INaL    = GNaL*(v-ENa)*_mL*((1.0-fINaLp)*_hL+fINaLp*_hLp);

//Ito ion current
T ass     = 1.0/(1.0+exp((-(v-14.34))/14.82));
T ta      = 1.0515/(1.0/(1.2089*(1.0+exp(-(v-18.4099)/29.3814)))+3.5/(1.0+exp((v+100.0)/29.3814)));
_a        = ass-(ass-a)*exp(-dt/ta);
T iss     = 1.0/(1.0+exp((v+43.94)/5.711));
T delta_epi;
#ifdef CELLTYPE_EPI
delta_epi = 1.0-(0.95/(1.0+exp((v+70.0)/5.0)));
#else
delta_epi=1.0;
#endif
T tiF         = 4.562+1/(0.3933*exp((-(v+100.0))/100.0)+0.08004*exp((v+50.0)/16.59));
T tiS         = 23.62+1/(0.001416*exp((-(v+96.52))/59.05)+1.780e-8*exp((v+114.1)/8.079));
tiF          *= delta_epi;
tiS          *= delta_epi;
T AiF         = 1.0/(1.0+exp((v-213.6)/151.2));
T AiS         = 1.0-AiF;
_iF           = iss-(iss-iF)*exp(-dt/tiF);
_iS           = iss-(iss-iS)*exp(-dt/tiS);
T c           = AiF*_iF+AiS*_iS;
T assp        = 1.0/(1.0+exp((-(v-24.34))/14.82));
_ap           = assp-(assp-ap)*exp(-dt/ta);
T dti_develop = 1.354+1.0e-4/(exp((v-167.4)/15.89)+exp(-(v-12.23)/0.2154));
T dti_recover = 1.0-0.5/(1.0+exp((v+70.0)/20.0));
T tiFp        = dti_develop*dti_recover*tiF;
T tiSp        = dti_develop*dti_recover*tiS;
_iFp          = iss-(iss-iFp)*exp(-dt/tiFp);
_iSp          = iss-(iss-iSp)*exp(-dt/tiSp);
T ip          = AiF*_iFp+AiS*_iSp;
T Gto         = 0.02;
#ifdef CELLTYPE_EPI
Gto*=4.0;
#endif
#ifdef CSELLTYPE_M
Gto*=4.0;
#endif
T fItop       = (1.0/(1.0+{KmCaMK_}/CaMKa));
T Ito          = Gto*(v-EK)*((1.0-fItop)*_a*c+fItop*_ap*ip);

// Ca ion current
T dss     = 1.0/(1.0+exp((-(v+3.940))/4.230));
T td      = 0.6+1.0/(exp(-0.05*(v+6.0))+exp(0.09*(v+14.0)));
_d        = dss-(dss-d)*exp(-dt/td);
T fss     = 1.0/(1.0+exp((v+19.58)/3.696));
T tff     = 7.0+1.0/(0.0045*exp(-(v+20.0)/10.0)+0.0045*exp((v+20.0)/10.0));
T tfs     = 1000.0+1.0/(0.000035*exp(-(v+5.0)/4.0)+0.000035*exp((v+5.0)/6.0));
T Aff     = 0.6;
T Afs     = 1.0-Aff;
_ff       = fss-(fss-ff)*exp(-dt/tff);
_fs       = fss-(fss-fs)*exp(-dt/tfs);
T f       = Aff*_ff+Afs*_fs;
T fcass   = fss;
T tfcaf   = 7.0+1.0/(0.04*exp(-(v-4.0)/7.0)+0.04*exp((v-4.0)/7.0));
T tfcas   = 100.0+1.0/(0.00012*exp(-v/3.0)+0.00012*exp(v/7.0));
T Afcaf   = 0.3+0.6/(1.0+exp((v-10.0)/10.0));
T Afcas   = 1.0-Afcaf;
_fcaf     = fcass-(fcass-fcaf)*exp(-dt/tfcaf);
_fcas     = fcass-(fcass-fcas)*exp(-dt/tfcas);
T fca     = Afcaf*_fcaf+Afcas*_fcas;
T tjca    = 75.0;
_jca      = fcass-(fcass-jca)*exp(-dt/tjca);
T tffp    = 2.5*tff;
_ffp      = fss-(fss-ffp)*exp(-dt/tffp);
T fp      = Aff*_ffp+Afs*_fs;
T tfcafp  = 2.5*tfcaf;
_fcafp    = fcass-(fcass-fcafp)*exp(-dt/tfcafp);
T fcap    = Afcaf*_fcafp+Afcas*_fcas;
T Kmn     = 0.002;
T k2n     = 1000.0;
T km2n    = _jca*1.0;
T anca    = 1.0/(k2n/km2n+pow(1.0+Kmn/cass,4.0));
_nca      = anca*k2n/km2n-(anca*k2n/km2n-nca)*exp(-km2n*dt);
T PhiCaL  = 4.0*vffrt*(cass*exp(2.0*vfrt)-0.341*{cao_})/(exp(2.0*vfrt)-1.0);
T PhiCaNa = 1.0*vffrt*(0.75*nass*exp(1.0*vfrt)-0.75*{nao_})/(exp(1.0*vfrt)-1.0);
T PhiCaK  = 1.0*vffrt*(0.75*kss*exp(1.0*vfrt)-0.75*{ko_})/(exp(1.0*vfrt)-1.0);
T PCa     = 0.0001;
#ifdef CELLTYPE_EPI
PCa*=1.2;
#endif
#ifdef CELLTYPE_M
PCa*=2.5;
#endif
T PCap    = 1.1*PCa;
T PCaNa   = 0.00125*PCa;
T PCaK    = 3.574e-4*PCa;
T PCaNap  = 0.00125*PCap;
T PCaKp   = 3.574e-4*PCap;
T fICaLp  = (1.0/(1.0+{KmCaMK_}/CaMKa));
T ICaL     = (1.0-fICaLp)*PCa*PhiCaL*_d*(f*(1.0-_nca)+_jca*fca*_nca)+fICaLp*PCap*PhiCaL*_d*(fp*(1.0-_nca)+_jca*fcap*_nca);
T ICaNa    = (1.0-fICaLp)*PCaNa*PhiCaNa*_d*(f*(1.0-_nca)+_jca*fca*_nca)+fICaLp*PCaNap*PhiCaNa*_d*(fp*(1.0-_nca)+_jca*fcap*_nca);
T ICaK     = (1.0-fICaLp)*PCaK*PhiCaK*_d*(f*(1.0-_nca)+_jca*fca*_nca)+fICaLp*PCaKp*PhiCaK*_d*(fp*(1.0-_nca)+_jca*fcap*_nca);

// IKr ion current
T xrss = 1.0/(1.0+exp((-(v+8.337))/6.789));
T txrf = 12.98+1.0/(0.3652*exp((v-31.66)/3.869)+4.123e-5*exp((-(v-47.78))/20.38));
T txrs = 1.865+1.0/(0.06629*exp((v-34.70)/7.355)+1.128e-5*exp((-(v-29.74))/25.94));
T Axrf = 1.0/(1.0+exp((v+54.81)/38.21));
T Axrs = 1.0-Axrf;
_xrf   = xrss-(xrss-xrf)*exp(-dt/txrf);
_xrs   = xrss-(xrss-xrs)*exp(-dt/txrs);
T xr   = Axrf*_xrf+Axrs*_xrs;
T rkr  = 1.0/(1.0+exp((v+55.0)/75.0))*1.0/(1.0+exp((v-10.0)/30.0));
T GKr  = 0.046;
#ifdef CELLTYPE_EPI
GKr*=1.3;
#endif
#ifdef CELLTYPE_M
GKr*=0.8;
#endif
T IKr  = GKr*sqrt({ko_}/5.4)*xr*rkr*(v-EK);

// IKs ion current
T xs1ss = 1.0/(1.0+exp((-(v+11.60))/8.932));
T txs1  = 817.3+1.0/(2.326e-4*exp((v+48.28)/17.80)+0.001292*exp((-(v+210.0))/230.0));
_xs1    = xs1ss-(xs1ss-xs1)*exp(-dt/txs1);
T xs2ss = xs1ss;
T txs2  = 1.0/(0.01*exp((v-50.0)/20.0)+0.0193*exp((-(v+66.54))/31.0));
_xs2    = xs2ss-(xs2ss-xs2)*exp(-dt/txs2);
T KsCa  = 1.0+0.6/(1.0+pow(3.8e-5/cai,1.4));
T GKs   = 0.0034;
#ifdef CELLTYPE_EPI
GKs*=1.4;
#endif
T IKs   = GKs*KsCa*_xs1*_xs2*(v-EKs);

// IK1 ion current
T xk1ss = 1.0/(1.0+exp(-(v+2.5538*{ko_}+144.59)/(1.5692*{ko_}+3.8115)));
T txk1  = 122.2/(exp((-(v+127.2))/20.36)+exp((v+236.8)/69.33));
_xk1    = xk1ss-(xk1ss-xk1)*exp(-dt/txk1);
T rk1   = 1.0/(1.0+exp((v+105.8-2.6*{ko_})/9.493));
T GK1   = 0.1908;
#ifdef CELLTYPE_EPI
GK1*=1.2;
#endif
#ifdef CELLTYPE_M
GK1*=1.3;
#endif
T IK1   = GK1*sqrt({ko_})*rk1*_xk1*(v-EK);

//NaCa ion current
T kna1    = 15.0;
T kna2    = 5.0;
T kna3    = 88.12;
T kasymm  = 12.5;
T wna     = 6.0e4;
T wca     = 6.0e4;
T wnaca   = 5.0e3;
T kcaon   = 1.5e6;
T kcaoff  = 5.0e3;
T qna     = 0.5224;
T qca     = 0.1670;
T hca     = exp((qca*v*{frdy_})/({R_}*temp));
T hna     = exp((qna*v*{frdy_})/({R_}*temp));
T h1      = 1+nai/kna3*(1+hna);
T h2      = (nai*hna)/(kna3*h1);
T h3      = 1.0/h1;
T h4      = 1.0+nai/kna1*(1+nai/kna2);
T h5      = nai*nai/(h4*kna1*kna2);
T h6      = 1.0/h4;
T h7      = 1.0+{nao_}/kna3*(1.0+1.0/hna);
T h8      = {nao_}/(kna3*hna*h7);
T h9      = 1.0/h7;
T h10     = kasymm+1.0+{nao_}/kna1*(1.0+{nao_}/kna2);
T h11     = {nao_}*{nao_}/(h10*kna1*kna2);
T h12     = 1.0/h10;
T k1      = h12*{cao_}*kcaon;
T k2      = kcaoff;
T k3p     = h9*wca;
T k3pp    = h8*wnaca;
T k3      = k3p+k3pp;
T k4p     = h3*wca/hca;
T k4pp    = h2*wnaca;
T k4      = k4p+k4pp;
T k5      = kcaoff;
T k6      = h6*cai*kcaon;
T k7      = h5*h2*wna;
T k8      = h8*h11*wna;
T x1      = k2*k4*(k7+k6)+k5*k7*(k2+k3);
T x2      = k1*k7*(k4+k5)+k4*k6*(k1+k8);
T x3      = k1*k3*(k7+k6)+k8*k6*(k2+k3);
T x4      = k2*k8*(k4+k5)+k3*k5*(k1+k8);
T E1      = x1/(x1+x2+x3+x4);
T E2      = x2/(x1+x2+x3+x4);
T E3      = x3/(x1+x2+x3+x4);
T E4      = x4/(x1+x2+x3+x4);
T KmCaAct = 150.0e-6;
T allo    = 1.0/(1.0+pow(KmCaAct/cai,2.0));
T JncxNa  = 3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
T JncxCa  = E2*k2-E1*k1;
T Gncx    = 0.0008;
#ifdef CELLTYPE_EPI
Gncx*=1.1;
#endif
#ifdef CELLTYPE_M
Gncx*=1.4;
#endif
T INaCa_i = 0.8*Gncx*allo*({zna_}*JncxNa+{zca_}*JncxCa);

h1         = 1+nass/kna3*(1+hna);
h2         = (nass*hna)/(kna3*h1);
h3         = 1.0/h1;
h4         = 1.0+nass/kna1*(1+nass/kna2);
h5         = nass*nass/(h4*kna1*kna2);
h6         = 1.0/h4;
h7         = 1.0+{nao_}/kna3*(1.0+1.0/hna);
h8         = {nao_}/(kna3*hna*h7);
h9         = 1.0/h7;
h10        = kasymm+1.0+{nao_}/kna1*(1+{nao_}/kna2);
h11        = {nao_}*{nao_}/(h10*kna1*kna2);
h12        = 1.0/h10;
k1         = h12*{cao_}*kcaon;
k2         = kcaoff;
k3p        = h9*wca;
k3pp       = h8*wnaca;
k3         = k3p+k3pp;
k4p        = h3*wca/hca;
k4pp       = h2*wnaca;
k4         = k4p+k4pp;
k5         = kcaoff;
k6         = h6*cass*kcaon;
k7         = h5*h2*wna;
k8         = h8*h11*wna;
x1         = k2*k4*(k7+k6)+k5*k7*(k2+k3);
x2         = k1*k7*(k4+k5)+k4*k6*(k1+k8);
x3         = k1*k3*(k7+k6)+k8*k6*(k2+k3);
x4         = k2*k8*(k4+k5)+k3*k5*(k1+k8);
E1         = x1/(x1+x2+x3+x4);
E2         = x2/(x1+x2+x3+x4);
E3         = x3/(x1+x2+x3+x4);
E4         = x4/(x1+x2+x3+x4);
KmCaAct    = 150.0e-6;
allo       = 1.0/(1.0+pow(KmCaAct/cass,2.0));
JncxNa     = 3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
JncxCa     = E2*k2-E1*k1;
T INaCa_ss = 0.2*Gncx*allo*({zna_}*JncxNa+{zca_}*JncxCa);
T INaCa    = INaCa_i+INaCa_ss;

// NaK ion current
T k1p    = 949.5;
T k1m    = 182.4;
T k2p    = 687.2;
T k2m    = 39.4;
k3p      = 1899.0;
T k3m    = 79300.0;
k4p      = 639.0;
T k4m    = 40.0;
T Knai0  = 9.073;
T Knao0  = 27.78;
T delta  = -0.1550;
T Knai   = Knai0*exp((delta*v*{frdy_})/(3.0*{R_}*temp));
T Knao   = Knao0*exp(((1.0-delta)*v*{frdy_})/(3.0*{R_}*temp));
T Kki    = 0.5;
T Kko    = 0.3582;
T MgADP  = 0.05;
T MgATP  = 9.8;
T Kmgatp = 1.698e-7;
T H      = 1.0e-7;
T eP     = 4.2;
T Khp    = 1.698e-7;
T Knap   = 224.0;
T Kxkur  = 292.0;
T P      = eP/(1.0+H/Khp+nai/Knap+ki/Kxkur);
T a1     = (k1p*pow(nai/Knai,3.0))/(pow(1.0+nai/Knai,3.0)+pow(1.0+ki/Kki,2.0)-1.0);
T b1     = k1m*MgADP;
T a2     = k2p;
T b2     = (k2m*pow({nao_}/Knao,3.0))/(pow(1.0+{nao_}/Knao,3.0)+pow(1.0+{ko_}/Kko,2.0)-1.0);
T a3     = (k3p*pow({ko_}/Kko,2.0))/(pow(1.0+{nao_}/Knao,3.0)+pow(1.0+{ko_}/Kko,2.0)-1.0);
T b3     = (k3m*P*H)/(1.0+MgATP/Kmgatp);
T a4     = (k4p*MgATP/Kmgatp)/(1.0+MgATP/Kmgatp);
T b4     = (k4m*pow(ki/Kki,2.0))/(pow(1.0+nai/Knai,3.0)+pow(1.0+ki/Kki,2.0)-1.0);
x1       = a4*a1*a2+b2*b4*b3+a2*b4*b3+b3*a1*a2;
x2       = b2*b1*b4+a1*a2*a3+a3*b1*b4+a2*a3*b4;
x3       = a2*a3*a4+b3*b2*b1+b2*b1*a4+a3*a4*b1;
x4       = b4*b3*b2+a3*a4*a1+b2*a4*a1+b3*b2*a1;
E1       = x1/(x1+x2+x3+x4);
E2       = x2/(x1+x2+x3+x4);
E3       = x3/(x1+x2+x3+x4);
E4       = x4/(x1+x2+x3+x4);
T JnakNa = 3.0*(E1*a3-E2*b3);
T JnakK  = 2.0*(E4*b1-E3*a1);
T Pnak   = 30;
#ifdef CELLTYPE_EPI
Pnak*=0.9;
#endif
#ifdef CELLTYPE_M
Pnak*=0.7;
#endif
T INaK   = Pnak*({zna_}*JnakNa+{zk_}*JnakK);

// Kb ion current
T xkb = 1.0/(1.0+exp(-(v-14.48)/18.34));
T GKb = 0.003;
#ifdef CELLTYPE_EPI
GKb*=0.6;
#endif
T IKb = GKb*xkb*(v-EK);

// Nab ion current
T PNab = 3.75e-10;
T INab = PNab*vffrt*(nai*exp(vfrt)-{nao_})/(exp(vfrt)-1.0);

// Cab ion current
T PCab = 2.5e-8;
T ICab = PCab*4.0*vffrt*(cai*exp(2.0*vfrt)-0.341*{cao_})/(exp(2.0*vfrt)-1.0);

// pCa ion current
T GpCa = 0.0005;
T IpCa = GpCa*cai/(0.0005+cai);


// Membrane voltage update
_it = INa+INaL+Ito+ICaL+ICaNa+ICaK+IKr+IKs+IK1+INaCa+INaK+INab+IKb+IpCa+ICab+st;
_v -= dt*_it;


//calculate fluxes, buffers and concentrations
_CaMKt      = CaMKt + dt*({aCaMK_}*CaMKb*(CaMKb+CaMKt)-{bCaMK_}*CaMKt);
T JdiffNa   = (nass-nai)/2.0;
T JdiffK    = (kss-ki)/2.0;
T Jdiff     = (cass-cai)/0.2;

T bt        = 4.75;
T a_rel     = 0.5*bt;
T Jrel_inf  = a_rel*(-ICaL)/(1.0+pow(1.5/cajsr,8.0));
#ifdef CELLTYPE_M
Jrel_inf*=1.7;
#endif
T tau_rel   = bt/(1.0+0.0123/cajsr);
T maskPosi;
T maskNega;
maskPosi    = (tau_rel > 0.005) * 1;
maskNega    = (tau_rel <= 0.005) * 1;
tau_rel     = maskNega * 0.005 + maskPosi * tau_rel;
_Jrelnp     = Jrel_inf-(Jrel_inf-Jrelnp)*exp(-dt/tau_rel);
T btp       = 1.25*bt;
T a_relp    = 0.5*btp;
T Jrel_infp = a_relp*(-ICaL)/(1.0+pow(1.5/cajsr,8.0));
#ifdef CELLTYPE_M
Jrel_infp*=1.7;
#endif
T tau_relp  = btp/(1.0+0.0123/cajsr);
maskPosi    = (tau_relp > 0.005) * 1;
maskNega    = (tau_relp <= 0.005) * 1;
tau_relp    = maskNega * 0.005 + maskPosi * tau_relp;
_Jrelp      = Jrel_infp-(Jrel_infp-Jrelp)*exp(-dt/tau_relp);
T fJrelp    = (1.0/(1.0+{KmCaMK_}/CaMKa));
T Jrel      = (1.0-fJrelp)*_Jrelnp+fJrelp*_Jrelp;

T Jupnp     = 0.004375*cai/(cai+0.00092);
T Jupp      = 2.75*0.004375*cai/(cai+0.00092-0.00017);
#ifdef CELLTYPE_EPI
Jupnp*=1.3;
Jupp*=1.3;
#endif
T fJupp     = (1.0/(1.0+{KmCaMK_}/CaMKa));
T Jleak     = 0.0039375*cansr/15.0;
T Jup       = (1.0-fJupp)*Jupnp+fJupp*Jupp-Jleak;

T Jtr       = (cansr-cajsr)/100.0;

_nai        = nai + dt*(-(INa+INaL+3.0*INaCa_i+3.0*INaK+INab)*{acap_}/({frdy_}*{vmyo_})+JdiffNa*{volrss_}/{vmyo_});
_nass       = nass + dt*(-(ICaNa+3.0*INaCa_ss)*{acap_}/({frdy_}*{volrss_})-JdiffNa);

_ki         = ki + dt*(-(Ito+IKr+IKs+IK1+IKb+st-2.0*INaK)*{acap_}/({frdy_}*{vmyo_})+JdiffK*{volrss_}/{vmyo_});
_kss        = kss + dt*(-(ICaK)*{acap_}/({frdy_}*{volrss_})-JdiffK);

T Bcai;
#ifdef CELLTYPE_EPI
Bcai=1.0/(1.0+1.3*{cmdnmax_}*{kmcmdn_}/pow({kmcmdn_}+cai,2.0)+{trpnmax_}*{kmtrpn_}/pow({kmtrpn_}+cai,2.0));
#else
Bcai=1.0/(1.0+{cmdnmax_}*{kmcmdn_}/pow({kmcmdn_}+cai,2.0)+{trpnmax_}*{kmtrpn_}/pow({kmtrpn_}+cai,2.0));
#endif
_cai        = cai + dt*(Bcai*(-(IpCa+ICab-2.0*INaCa_i)*{acap_}/(2.0*{frdy_}*{vmyo_})-Jup*{vnsr_}/{vmyo_}+Jdiff*{volrss_}/{vmyo_}));

T Bcass     = 1.0/(1.0+{BSRmax_}*{KmBSR_}/pow({KmBSR_}+cass,2.0)+{BSLmax_}*{KmBSL_}/pow({KmBSL_}+cass,2.0));
_cass       = cass + dt*(Bcass*(-(ICaL-2.0*INaCa_ss)*{acap_}/(2.0*{frdy_}*{volrss_})+Jrel*{vjsr_}/{volrss_}-Jdiff));

_cansr      = cansr + dt*(Jup-Jtr*{vjsr_}/{vnsr_});

T Bcajsr    = 1.0/(1.0+{csqnmax_}*{kmcsqn_}/pow({kmcsqn_}+cajsr,2.0));
_cajsr      = cajsr + dt*(Bcajsr*(Jtr-Jrel));


// unmodified return variables
_dt = dt;
_st = st;
_temp = temp;

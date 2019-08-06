double rtf_   = ({R})*temp/({F});
double e_naf_  = log(({NA_O})/na_i)*(rtf_); // fibroblast Nernst potential for Na+
double e_kf_ = log(({K_O})/k_i)*(rtf_);   // fibroblast Nernst potential for K+

// I_Kv
double tau_r_ = 20.3+138.*exp(-pow(((v+20.)/25.9), 2.0));
double tau_s_ = 1574.+5268.*exp(-pow(((v+23.)/22.7), 2.0));
double r_inf_ = pow(1+exp(-(v+20.)/11),-1.0);
double s_inf_ = pow(1+exp(-(v+23.)/7),-1.0);
_r = ( (r-r_inf_) / (tau_r_) ) * dt;
_s = ( (s-s_inf_) / (tau_s_) ) * dt;
_ikv = ({G_KV})*r*s*(v-e_kf_);
 
// I_K1
double alpha_ = 0.1 * pow( 1.0+exp(0.06*(v-e_kf_-200.0)),-1.0);
double beta_ = ( 3.0*exp(2.0e-4*(v-e_kf_+100.0))+exp(0.1*(v-e_kf_-10.0)))/(1.0+exp(-0.5*(v-e_kf_)));
_ik1 = ({G_K1})*(alpha_/(alpha_+beta_))*(v-e_kf_);

// I_NaK
_inak = double(({P_NAK}));
_inak *= double(({K_O})/(({K_O})+({K_MK})));
_inak *= pow(na_i,1.5)/( pow(na_i,1.5)+pow(({K_MNA}),1.5) );
_inak *= (v-({V_REV})) / (v+200.0);

// // I_Nab
_inab = ({G_NAB}) * (v-e_naf_);

// total current 
_ina_net = _inab + 3*_inak;
_ik_net  = _ikv + _ik1 - 2*_inak;
_it = sw_it * (_ina_net + _ik_net + st);

// ion concentration update
double dki_ = -dt*(_ik_net)*(({ACAP})/(({V_MYO})*({F})));
_k_i = k_i + dki_;
double dnai_ = -dt*(_ina_net)*(({ACAP})/(({V_MYO})*({F})));
_na_i = na_i + dnai_;

// constant variables
_v = v;
_dt = dt;
_st = st;
_temp = temp;
_sw_it = sw_it;

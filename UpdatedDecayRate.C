#include "UpdatedDecayRate.h"
#include <cstddef>

//______________________________________________________________________________
double decayRates::dGdwD(FFModel pFF,double VcbEff, double w){

  double Gw = 1.;
  double f0 = 1.;
  
  if(pFF.model=="CLN")
      FFfunctionsCLN(w,pFF.FFpars,Gw,f0);
  else if(pFF.model=="BGL")
      FFfunctionsBGL(w,pFF.FFpars,Gw,f0);
  else{
      cout << " \t Unknown FF model in dGdw! \n";
      return 0;
  }
    
  double GF  = Constants::GF;
  double pi3 = pow(TMath::Pi(),3);
  double mB  = Mass::Bs/1000;
  double mD  = Mass::Ds/1000;
  
  return GF*GF *mD*mD*mD /48./pi3 *(mB+mD)*(mB+mD)
               *sqrt((w*w-1.)*(w*w-1.)*(w*w-1.))
               *VcbEff*VcbEff * Gw*Gw;

}

//_________________________________________________________________________________________
void decayRates::FFfunctionsCLN(double w, vector<parameter> FFpar, double& Gw, double& f0){
                                                                                   
  double rho2 = 0;
  double G0   = 0;
  for(auto p : FFpar){
    if     (p.name=="G0")    G0   = p.value;
    else if(p.name=="rho2P") rho2 = p.value;
	else { cout << "!!! ERROR FFfunctionsCLN for B->D, no valid parameters!" << endl;
	       Gw=0; f0=0.; return; }
  }

  double z = (sqrt(w+1.) - sqrt(2.))/(sqrt(w+1.) + sqrt(2.));
  Gw = G0*(1. - 8.*rho2*z + (51.*rho2 - 10.)*z*z - (252.*rho2 - 84.)*z*z*z);
  
  double r     =  Mass::Ds/Mass::Bs;
  double fplus = (1+r)/2./sqrt(r)*Gw;
  double w1    = w - 1.;
  f0 = fplus* (4.*r/(1.+r)/(1.+r))*(1+w)/2. * 1.0036*
 			  (1. - 0.0068*w1 + 0.0017*w1*w1 - 0.0013*w1*w1*w1);
 
  return;
}

//_________________________________________________________________________________________
void decayRates::FFfunctionsBGL(double w, vector<parameter> FFpar, double& Gw, double& f0){
    
    //coeff of the series
    const int maxorder = 4;
    double dp[maxorder] = {0.};
    double ep[maxorder] = {0.};
    double G0=1;

    //link now the parameters to the coeff of the series
    int count_d=0; int count_e=0;
    for(auto p : FFpar){
        if(p.name.Contains("d")) {
		  dp[count_d] = p.value; ++count_d;
        }else if(p.name.Contains("e")) {
		  ep[count_e] = p.value; ++count_e;
		}else if(p.name == "G0")  G0 = p.value;
        else { cout << "!!! ERROR FFfunctionsBGL for B->D, no valid parameters!" << endl;
		  Gw=0.; f0=0.; return; };
    }

	double z = (sqrt(w+1.) - sqrt(2.))/(sqrt(w+1.) + sqrt(2.));
    
    double mB  = Mass::Bs/1000;
    double mD  = Mass::Ds/1000;
    double r = mD/mB;
    double tplus = (mB+mD)*(mB+mD);
    double tminus = (mB-mD)*(mB-mD);
    double mBv[5]={6.32847,6.91947,7.030,7.280,0}; //new poles
   // double mBv[5]={6.329,6.920,7.020,7.280,0}; //old poles
    double P_p  = P1pm(z,tplus,tminus,mBv);
       
     //outer function
     double nI    = 2.6;
     double chiT0 = 5.131e-4;//Gambino paper 1606.08030
     double kp    = 8.*r*r/mB*sqrt(8.*nI/3./TMath::Pi()/chiT0);
     double phi_p = kp*(1.+z)*(1.+z)*sqrt(1.-z)/pow( (1.+r)*(1.-z) + 2*sqrt(r)*(1.+z), 5);

     //make the first coefficient linked to G0
     double phi_p0 = kp/pow(1+r +2*sqrt(r),5);
     double P_p0   = P1pm(0,tplus,tminus,mBv);
     double d0 = (1.+r)/2./sqrt(r)*phi_p0*P_p0;
     //build the series
     double Sum_z = d0; int order = 1;
	 //cout << " ============== " << endl;
     //cout << " d0 = " << d0 << endl;
     for(auto cdp : dp){
	   //cout << " d" << order << " = " << cdp << endl;
	   Sum_z += cdp*pow(z, order);  ++order;
	 }
     double f_p = Sum_z/phi_p/P_p;
     //here it is G(w)
     Gw = G0*2*sqrt(r)/(1+r)*f_p;
  
     //now let's do f0:
     double mBa[5]={6.70347,7.122,0,0,0};
     double P_z = P1pm(z,tplus,tminus,mBa);
   
     double chiL0 = 6.204e-3;
     double k0    = r*(1.-r*r)*sqrt(8.*nI/TMath::Pi()/chiL0);
     double phi_z = k0*(1.-z*z)*sqrt(1.-z)/pow( (1.+r)*(1.-z) + 2*sqrt(r)*(1.+z), 4);
  
	 //make the strong contraints on f0: f0(q2=0) = f+(q2=0)
     double wmax = (mB*mB + mD*mD) / (2.*mB*mD);
     double zq20 = (sqrt(wmax+1.) - sqrt(2.))/(sqrt(wmax+1.) + sqrt(2.));
     //build f+ at q2=0
	 double P_pq20 = P1pm(zq20,tplus,tminus,mBv);
	 double phi_pq20 = kp*(1.+zq20)*(1.+zq20)*sqrt(1.-zq20)/pow( (1.+r)*(1.-zq20) + 2*sqrt(r)*(1.+zq20), 5);
	 double Sum_zq20 = d0; order = 1;
	 for(auto cdp : dp){ Sum_zq20 += cdp*pow(zq20, order);  ++order;}
     double f_pq20 = Sum_zq20/phi_pq20/P_pq20;
     //build the "truncated f0 series": that without one coeff
	 double P_zq20 = P1pm(zq20,tplus,tminus,mBa);
     double phi_zq20 = k0*(1.-zq20*zq20)*sqrt(1.-zq20)/pow( (1.+r)*(1.-zq20) + 2*sqrt(r)*(1.+zq20), 4);
	 double Sum_f0q20 = 0;
     for(int ij=1; ij<maxorder; ij++){ Sum_f0q20 += ep[ij]*pow(zq20, ij);}
	
     double e0 = f_pq20*P_zq20*phi_zq20 - Sum_f0q20;
	 //double e2 = (f_pq20*P_zq20*phi_zq20 - f_0q20)/zq20/zq20;
    
	 double Sum_f0 = 0;
	 //to apply strong constraints on f0, uncomment below
	 //(note: e0 must be a fixed parameter in the config file)
     //my realtion:
	 //ep[0] = e0;
	 //Sneha relation:
     //ep[0] = 6.500331203665146*dp[0] + 0.3922131457342012*dp[1] +
	 //       0.023665125186233833*dp[2] + 0.0014278923492779892*dp[3] -
	 //		 0.060337409502004426*ep[1] - 0.003640602985412573*ep[2] -
	 //		 0.00021966455316505822*ep[3];

     for(int ij=0; ij<maxorder; ij++){
	   //cout << " e" << ij << " = " << ep[ij] << endl;
	   Sum_f0 += ep[ij]*pow(z, ij);
	 }
  
	 f0 = Sum_f0/phi_z/P_z;
 
	 return;

}

//______________________________________________________________________________
double decayRates::P1pm(double z, double tplus, double tminus, double* mBc){

  double zp[5];
  
  int nRes = 5;
  if(mBc[4]==0) nRes=4;
  if(mBc[3]==0) nRes=3;
  if(mBc[2]==0) nRes=2;
  
  for(int i=0;i<nRes;++i)
    zp[i]=(sqrt(tplus - mBc[i]*mBc[i]) - sqrt(tplus - tminus))/(sqrt(tplus - mBc[i]*mBc[i]) + sqrt(tplus - tminus));

  double f=1;
  for(int i=0; i<nRes; ++i) f *= (z-zp[i])/(1.-z*zp[i]);

  return f;

}


//____________________________________________________________________________________________________________________
void decayRates::FFfunctionsCLN(double w, vector<parameter> FFpar, double& hw, double& R1w, double& R2w, double& R0w){

    double F1   = 0;
    double rho2 = 0;
    double R1   = 0;
    double R2   = 0;
    double R0   = 0;
    
    for(auto p : FFpar){
      if     (p.name == "F1")    F1   = p.value;
      else if(p.name == "rho2V") rho2 = p.value;
      else if(p.name == "R1")    R1   = p.value;
      else if(p.name == "R2")    R2   = p.value;
	  else if(p.name == "R0")    R0   = p.value;
      else { cout << "!!! ERROR FFfunctionsCLN for B->D*, no valid parameters!" << endl;
		hw=0; R1w=1; R2w=1; R0w=1; return;}
    }

    double z = (sqrt(w+1.) - sqrt(2.))/(sqrt(w+1.) + sqrt(2.));

    hw  = F1*(1. - 8.*rho2*z + (53.*rho2 - 15.)*z*z - (231.*rho2 - 91.)*z*z*z);
    R0w = R0 - 0.11*(w-1.) + 0.01*(w-1.)*(w-1.);
    R1w = R1 - 0.12*(w-1.) + 0.05*(w-1.)*(w-1.);
    R2w = R2 + 0.11*(w-1.) - 0.06*(w-1.)*(w-1.);
    
    return;

}

//________________________________________________________________________________________________________________
void decayRates::FFfunctionsBGL(double w, vector<parameter> FFpar, double& A1, double& A2, double& V, double& A0){

  //these are the coeff of the 3 series
  double bp[5] = {0.};
  double cp[5] = {0.};
  double ap[5] = {0.};
  double hp[5] = {0.};

  //link now the parameters to the coeff of the series
  int count_b=0;
  int count_c=0;
  int count_a=0;
  int count_h=0;
  double F1=0;
  for(auto p : FFpar){
      if     (p.name.Contains("b")) {
        bp[count_b] = p.value; ++count_b;
      }
      else if(p.name.Contains("c")) {
        cp[count_c] = p.value; ++count_c;
      }
      else if(p.name.Contains("a")) {
        ap[count_a] = p.value; ++count_a;
      }
	  else if(p.name.Contains("h")) {
		hp[count_h] = p.value; ++count_h;
	  }
      else if(p.name == "F1") {
		F1 = p.value;
      }
      else { cout << "!!! ERROR FFfunctionsBGL for B->D*, no valid parameters!" << endl;
		A1=1; A2=1; V=1; A0=1; return;
      };
  }
  
  double z  = (sqrt(w+1.) - sqrt(2.))/(sqrt(w+1.) + sqrt(2.));
  double mB = Mass::Bs/1000;
  double mD = Mass::DsS/1000;
  double r  = mD/mB;
  double tplus  = (mB+mD)*(mB+mD);
  double tminus = (mB-mD)*(mB-mD);

  double mBv[5]={6.32847,6.91947,7.030,7.280,7.365};//new poles
  double mBa[5]={6.73847,6.750  ,7.145,7.150,0};//new poles
  //double mBv[5]={6.329,6.920,7.020,7.280,0};//old poles
  //double mBa[5]={6.739,6.750,7.145,7.150,0};//old poles
  
  double mBp[5]={6.27447,6.8712,7.250,0,0};
  
  double nI = 2.6;
  double chiT1m = 5.131e-4;
  double chiT1p = 3.894e-4;
  double chi0m = 19.421e-3;

  double b0 = 2.*sqrt(mB*mD)*P1pm(0,tplus,tminus,mBa)*
              4.*r/mB/mB*sqrt(nI/3./TMath::Pi()/chiT1p)/
              pow( (1.+r) + 2.*sqrt(r), 4)*F1;

  double r_phi_F1f_0 = 1./mB/sqrt(2.)/((1.+r) + 2.*sqrt(r));
  double c0 = (mB-mD)* r_phi_F1f_0* b0;

  double phi_f = 4.*r/mB/mB * sqrt(nI/3./TMath::Pi()/chiT1p)*
                 (1.+z)*sqrt((1.-z)*(1.-z)*(1.-z))/
                 pow( (1.+r)*(1.-z) + 2.*sqrt(r)*(1.+z), 4);

  double phi_g = sqrt(nI/3./TMath::Pi()/chiT1m)*
                  (16.*r*r*(1.+z)*(1.+z)/sqrt(1.-z))/
                  pow( (1.+r)*(1.-z) + 2.*sqrt(r)*(1.+z), 4);

  double phi_F1 = 4.*r/mB/mB/mB * sqrt(nI/6./TMath::Pi()/chiT1p)*
                  (1.+z)*sqrt((1.-z)*(1.-z)*(1.-z)*(1.-z)*(1.-z))/
                  pow( (1.+r)*(1.-z) + 2.*sqrt(r)*(1.+z), 5);
  
  double phi_F2 = 8.*r*r * sqrt(2.*nI/TMath::Pi()/chi0m)*
				  (1.+z)*(1.+z)/sqrt(1.-z)/
				  pow( (1.+r)*(1.-z) + 2.*sqrt(r)*(1.+z), 4);

  double p1p  = P1pm(z,tplus,tminus,mBa);

  double sumf = b0; int order = 1;
  //cout << "b0 = " << b0 << endl;
  for(auto cbp : bp){ sumf += cbp*pow(z, order); ++order;}
  double fBGL = sumf/p1p/phi_f;

  double p1m  = P1pm(z,tplus,tminus,mBv);
  double sumg = 0; order =0;
  for(auto cap : ap){ sumg += cap*pow(z, order); ++order;}
  double gBGL = sumg/p1m/phi_g;

  double sumF1 = c0; order =1;
  for(auto ccp : cp){ sumF1 += ccp*pow(z, order); ++order;}
  double F1BGL = sumF1/p1p/phi_F1;
  
  double p0m  = P1pm(z,tplus,tminus,mBp);
  double sumF2 = 0; order =0;
  //to apply strong constraints on F2 uncomment below
  //(note: h0 must be a fixed parameter)
  // Sneha relation:
  //hp[0] = 4.310226659869955*b0 + 1.395860907385888*cp[0] + 0.07319749527752503*cp[1] + 0.0038384005788494402*cp[2] - 0.05243896070891928*hp[1] - 0.002749844600231581*hp[2] - 0.00014419899294717772*hp[3];
  
  for(auto hap : hp){
	//cout << "h" << order << " = " << hap << endl;
	sumF2 += hap*pow(z, order); ++order;
  }
  double F2BGL = sumF2/p0m/phi_F2;

  double q2  = mB*mB + mD*mD - 2.*mB*mD*w;
  
  double lambda = mB*mB*mB*mB + mD*mD*mD*mD + q2*q2 - 2.0*(mB*mB*mD*mD + mB*mB*q2 + mD*mD*q2);
  
  V  = (mB+mD)/2.0*gBGL;
  A1 = fBGL/(mB+mD);
  A2 = (mB+mD)/lambda*( (mB*mB-mD*mD-q2)*fBGL - 2.0*mD*F1BGL );
  A0 = F2BGL/2.;
 
  return;
}

//Decay width differential in w, Bs->Ds* mu nu
//Vcb,etaEw,F1,rho2,b2,R1,R2,model
//______________________________________________________________________________
double decayRates::dGdwDst(FFModel pFF, double VcbEff, double w){

  double GF    = Constants::GF;
  double pi3   = pow(TMath::Pi(),3);
    
  double mB    = Mass::Bs/1000;
  double mD    = Mass::DsS/1000;
  double r = mD/mB;

  double hw  = 1;
  double R1w = 1;
  double R2w = 1;
  double R0w = 1;

  if(pFF.model=="CLN")
      FFfunctionsCLN(w,pFF.FFpars,hw,R1w,R2w,R0w);
  else if(pFF.model=="BGL"){
	double A1{1}, A2{1}, V{1}, A0{1};
    FFfunctionsBGL(w,pFF.FFpars,A1,A2,V,A0);
    double R = 2.*sqrt(r)/(1.+r);
	hw = 2./(w+1.)*A1/R;
	R1w = V *R/hw;
	R2w = A2*R/hw;
  }
  else{
      cout << " \t Unknown FF model in dGdw! \n";
      return 0;
  }
 
  double XiFw2 = hw*hw * sqrt(w*w-1.) *(w+1.)*(w+1.)
    * ( 2. * (1. - 2.*w*r + r*r)/(1.-r)/(1.-r) * (1. + R1w*R1w *(w -1.)/(w+1.) )  +
    ( 1. + (1. - R2w)*(w-1.)/(1.-r)) * ( 1. + (1. - R2w)*(w-1.)/(1.-r)) );

  return GF*GF *mD*mD*mD /48./pi3 *(mB-mD)*(mB-mD) *XiFw2 *VcbEff*VcbEff;

}


//______________________________________________________________________________
double decayRates::dGdwdAnglesDst(FFModel pFF, double VcbEff, double w, double ctl, double ctd, double chi){

  double GF    = Constants::GF;
  double pi4   = pow(TMath::Pi(),4);

  double mB    = Mass::Bs/1000;
  double mD    = Mass::DsS/1000;
  double r = mD/mB;

  double hw  = 1.;
  double R1w = 1.;
  double R2w = 1.;
  double R0w = 1.;

  if(pFF.model=="CLN")
	FFfunctionsCLN(w,pFF.FFpars,hw,R1w,R2w,R0w);
  else if(pFF.model=="BGL"){
	double A1{1}, A2{1}, V{1}, A0{1};
	FFfunctionsBGL(w,pFF.FFpars,A1,A2,V,A0);
	double R = 2.*sqrt(r)/(1.+r);
	hw = 2./(w+1.)*A1/R;
	R1w = V *R/hw;
	R2w = A2*R/hw;
  }
  else{
	cout << " \t Unknown FF model in dGdw! \n";
	return 0;
  }

  double arg  = 1. - 2.*w*r+ r*r;
  double qfac = (arg<0.)?  0. : sqrt(arg)/(1.-r);
  double hP = qfac*(1. - sqrt( (w-1.)/(w+1.) )*R1w);
  double hM = qfac*(1. + sqrt( (w-1.)/(w+1.) )*R1w);
  double h0 = 1. + (w-1.)/(1.-r)*(1. - R2w);
  double rp = 2.*sqrt(mB*mD)/(mB+mD);
  double hfactor = -mB*rp*(1.-r*r)*(w+1)/2./sqrt(1.-2.*w*r + r*r) * hw;
  hP*=hfactor;
  hM*=hfactor;
  h0*=hfactor;

  //angular functions for the Pi0
  double fPP = (1.-ctd*ctd) * (1.-ctl)*(1.-ctl);
  double fMM = (1.-ctd*ctd) * (1.+ctl)*(1.+ctl);
  double f00 =  4. * ctd*ctd  * (1.-ctl*ctl);
  double fPM = -2. * (1.-ctd*ctd) * (1.-ctl*ctl) * cos(2.*chi);
  double fP0 = -4. * sqrt(1.-ctd*ctd)*ctd * sqrt(1.-ctl*ctl) *(1.-ctl) *cos(chi);
  double fM0 =  4. * sqrt(1.-ctd*ctd)*ctd * sqrt(1.-ctl*ctl) *(1.+ctl) *cos(chi);

  double angularTermPi =
    hP*hP * fPP +
    hM*hM * fMM +
    h0*h0 * f00 +
    hP*hM * fPM +
    hP*h0 * fP0 +
    hM*h0 * fM0 ;

  //angular functions for the gamma
  ctl = -ctl;
    
  fPP = 2*(1. + ctd*ctd)*(1./4. + 1./4.*ctl*ctl + 1./2.*ctl);
  fMM = 2*(1. + ctd*ctd)*(1./4. + 1./4.*ctl*ctl - 1./2.*ctl);
  f00 = (1. - ctd*ctd)*2*(1. - ctl*ctl);
  fPM = (1. - ctd*ctd)*(1. - ctl*ctl)*(2*cos(chi)*cos(chi) - 1.);
  fP0 = (2*ctd*sqrt(1. - ctd*ctd))*cos(chi)*(ctl + 1.)*sqrt(1. - ctl*ctl);
  fM0 = (2*ctd*sqrt(1. - ctd*ctd))*cos(chi)*(ctl - 1.)*sqrt(1. - ctl*ctl);

  double angularTermGamma =
    hP*hP * fPP +
    hM*hM * fMM +
    h0*h0 * f00 +
    hP*hM * fPM +
    hP*h0 * fP0 +
    hM*h0 * fM0 ;

  double angularTerm = 0.06*angularTermPi + 0.94* angularTermGamma;
 
  //Put all together, here the full rate with all factors
  double Norm = 3.*GF*GF/4./(256.*pi4)*mB *mD*mD* VcbEff*VcbEff;

  return Norm * sqrt(w*w-1.)*(1.-2.*w*r + r*r) * angularTerm;

}

//______________________________________________________________________________
double decayRates::dGdwdctdDst(FFModel pFF, double VcbEff, double w, double ctd){

  double GF    = Constants::GF;
  double pi4   = pow(TMath::Pi(),4);

  double mB    = Mass::Bs/1000;
  double mD    = Mass::DsS/1000;
  double r = mD/mB;

  double hw  = 1.;
  double R1w = 1.;
  double R2w = 1.;
  double R0w = 1.;

  if(pFF.model=="CLN")
	FFfunctionsCLN(w,pFF.FFpars,hw,R1w,R2w,R0w);
  else if(pFF.model=="BGL"){
	double A1{1}, A2{1}, V{1}, A0{1};
	FFfunctionsBGL(w,pFF.FFpars,A1,A2,V,A0);
	double R = 2.*sqrt(r)/(1.+r);
	hw = 2./(w+1.)*A1/R;
	R1w = V *R/hw;
	R2w = A2*R/hw;
  }
  else{
	cout << " \t Unknown FF model in dGdw! \n";
	return 0;
  }

  double arg  = 1. - 2.*w*r+ r*r;
  double qfac = (arg<0.)?  0. : sqrt(arg)/(1.-r);
  double hP = qfac*(1. - sqrt( (w-1.)/(w+1.) )*R1w);
  double hM = qfac*(1. + sqrt( (w-1.)/(w+1.) )*R1w);
  double h0 = 1. + (w-1.)/(1.-r)*(1. - R2w);
  double rp = 2.*sqrt(mB*mD)/(mB+mD);
  double hfactor = -mB*rp*(1.-r*r)*(w+1)/2./sqrt(1.-2.*w*r + r*r) * hw;
  hP*=hfactor;
  hM*=hfactor;
  h0*=hfactor;

  //angular functions for the Pi0
  double fPP = (1.-ctd*ctd) * 8./3;
  double fMM = (1.-ctd*ctd) * 8./3;
  double f00 =  4. * ctd*ctd  * 4./3;

  double angularTermPi =
    hP*hP * fPP +
    hM*hM * fMM +
    h0*h0 * f00 ;

  //angular functions for the gamma
  fPP = 2*(1. + ctd*ctd)*2./3;
  fMM = 2*(1. + ctd*ctd)*2./3;
  f00 = (1. - ctd*ctd)*2*4./3;

  double angularTermGamma =
    hP*hP * fPP +
    hM*hM * fMM +
    h0*h0 * f00 ;

  double angularTerm = 0.06*angularTermPi + 0.94* angularTermGamma;
 
  //Put all together, here the full rate with all factors
  double Norm = 3.*GF*GF/4./(256.*pi4)*mB *mD*mD* VcbEff*VcbEff;

  return Norm * sqrt(w*w-1.)*(1.-2.*w*r + r*r) * angularTerm;

}


//______________________________________________________________________________
double decayRates::calculateGamma(FFModel pFF, double VcbEff, bool isDst){
    
    if(pFF.FFpars.size()==0){
       cout << "No FF PARAMETERS to calculateBR!" << endl;
       return 0;
    }
    
    double wmin = 1.;
    double mB    = Mass::Bs/1000;
    double mD    = isDst? Mass::DsS/1000 : Mass::Ds/1000;
    double wmax  = (mB*mB + mD*mD) / (2.*mB*mD);
    
    decayRates * thisModel = new decayRates(pFF,VcbEff,isDst);
    TF1 * f = new TF1("f",thisModel,&decayRates::TF_dGdw, wmin, wmax, 0);

    double gamma = f->Integral(wmin,wmax);
    
    return gamma;

}

//______________________________________________________________________________
double decayRates::calculateBR(FFModel pFF, double VcbEff, double tauB, bool isDst){
    
    double Gamma = calculateGamma(pFF,VcbEff, isDst);
    
    double norm = Constants::hbar/tauB;
    double BR = Gamma/norm;
    
    return BR;
    
}

//______________________________________________________________________________
double decayRates::pperp_w_cos(double w, double cos, double massD){
    
  return massD*sqrt(w*w - 1.)*sqrt(1. - cos*cos);
    
}

//______________________________________________________________________________
double decayRates::Pm1(double M, double m1, double m2){
    
  double sum = M*M - (m1+m2)*(m1+m2);
  double dif = M*M - (m1-m2)*(m1-m2);
    
  return sqrt(sum*dif)/(2.*M);
}

//______________________________________________________________________________
double decayRates::Em1(double M, double m1, double m2){
    
  return (M*M - m2*m2 + m1*m1)/(2.*M);
    
}

//______________________________________________________________________________
double decayRates::pperp_w_ctd_cos(double beta, double ctd, double cos, double M, double m1, double m2){
    
  double gamma = 1./sqrt(1.-beta*beta);
    
  double p = Pm1(M,m1,m2);
  double E = Em1(M,m1,m2);
    
  double pTransf_perp =  p* sqrt(1. -ctd*ctd);

  double p_parallel = p* ctd;
  double pTransf_parallel = beta*gamma*E + gamma*p_parallel;
   
  double pTransf = sqrt(pTransf_parallel*pTransf_parallel + pTransf_perp*pTransf_perp);
    
  double sinalpha = sqrt(1. - cos*cos);
        
  return pTransf * sinalpha;
    
}

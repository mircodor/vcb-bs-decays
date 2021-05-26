#include "decayRates.h"

//______________________________________________________________________________
double decayRates::dGdwD(FFModel pFF,double VcbEff, double w){

  double Gw = 1.;
    
  if(pFF.model=="CLN")
      FFfunctionsCLN(w,pFF.FFpars,Gw);
  else if(pFF.model=="BGL")
      FFfunctionsBGL(w,pFF.FFpars,Gw);
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

//______________________________________________________________________________
void decayRates::FFfunctionsCLN(double w, vector<parameter> FFpar, double& Gw){
                                                                                   
  double rho2 = 0;
  double G0   = 0;
  for(auto p : FFpar){
    if     (p.name=="G0")    G0   = p.value;
    else if(p.name=="rho2P") rho2 = p.value;
    else { cout << "!!! ERROR FFfunctionsCLN for B->D, no valid parameters!" << endl;  Gw=0; return;}
  }

  double z = (sqrt(w+1.) - sqrt(2.))/(sqrt(w+1.) + sqrt(2.));
  Gw  = G0*(1. - 8.*rho2*z + (51.*rho2 - 10.)*z*z - (252.*rho2 - 84.)*z*z*z);

  return;
}

//______________________________________________________________________________

void decayRates::FFfunctionsBGL(double w, vector<parameter> FFpar, double& Gw){
    
    //da ricontrollare!!!
    
    //coeff of the series
    double dp[10] = {0.};
    double G0=1;

    //link now the parameters to the coeff of the series
    int count_d=0; double F1;
    for(auto p : FFpar){
        if     (p.name.Contains("d")) {
                dp[count_d] = p.value; ++count_d;
        }
        else if(p.name == "G0") G0 = p.value;
        else { cout << "!!! ERROR FFfunctionsBGL for B->D, no valid parameters!" << endl;
                Gw = 0.; return; };
    }

     double z = (sqrt(w+1.) - sqrt(2.))/(sqrt(w+1.) + sqrt(2.));
    
     double mB  = Mass::Bs/1000;
     double mD  = Mass::Ds/1000;
     double r = mD/mB;
     double tplus = (mB+mD)*(mB+mD);
     double tminus = (mB-mD)*(mB-mD);
     double mBv[4]={6.329,6.920,7.020,7.280};
     double Pvscale = 1.;//2.52733;
     double P_p  = Pvscale*P1pm(z,tplus,tminus,mBv);
     
     //outer function
     //double Phi_p = 1.1213; //Belle paper and 1503.07237
     double nI    = 2.6;
     double chiT0 = 5.131e-4;//Gambino paper 1606.08030
     double Phi_p = 8.*r*r/mB*sqrt(8.*nI/3./TMath::Pi()/chiT0);
     double phi_p = Phi_p*(1.+z)*(1.+z)*sqrt(1.-z)/pow( (1.+r)*(1.-z) + 2*sqrt(r)*(1.+z), 5);

     //series and p_p
     double phi_p0 = Phi_p/pow(1+r +2*sqrt(r),5);
     double P_p0 = Pvscale*P1pm(0,tplus,tminus,mBv);

     double d0 = (1.+r)/2./sqrt(r)*phi_p0*P_p0;

     double Sum_z = d0; int order = 1;
     for(auto cdp : dp){ Sum_z += cdp*pow(z, order);  ++order;}
  
     double f_p = Sum_z/phi_p/P_p;

     Gw = G0*2*sqrt(r)/(1+r)*f_p;
    
    return;

}

//______________________________________________________________________________
double decayRates::P1pm(double z, double tplus, double tminus, double* mBc){

  double zp[4];

  int nRes =4;
  if(mBc[3]==0) nRes=3;

  for(int i=0;i<nRes;++i)
    zp[i]=(sqrt(tplus - mBc[i]*mBc[i]) - sqrt(tplus - tminus))/(sqrt(tplus - mBc[i]*mBc[i]) + sqrt(tplus - tminus));

  double f=1;
  for(int i=0; i<nRes; ++i) f *= (z-zp[i])/(1.-z*zp[i]);

  return f;

}


//______________________________________________________________________________
void decayRates::FFfunctionsCLN(double w, vector<parameter> FFpar, double& hw, double& R1w, double& R2w){

    double F1   = 0;
    double rho2 = 0;
    double R1   = 0;
    double R2   = 0;
    
    for(auto p : FFpar){
      if     (p.name == "F1")    F1 = p.value;
      else if(p.name == "rho2V") rho2 = p.value;
      else if(p.name == "R1")    R1   = p.value;
      else if(p.name == "R2")    R2   = p.value;
      else { cout << "!!! ERROR FFfunctionsCLN for B->D*, no valid parameters!" << endl;
          hw=0; R1w=1; R2w=1; return;}
    }

    double z = (sqrt(w+1.) - sqrt(2.))/(sqrt(w+1.) + sqrt(2.));

    hw  = F1*(1. - 8.*rho2*z + (53.*rho2 - 15.)*z*z - (231.*rho2 - 91.)*z*z*z) ;
    R1w = R1 - 0.12*(w-1.) + 0.05*(w-1.)*(w-1.);
    R2w = R2 + 0.11*(w-1.) - 0.06*(w-1.)*(w-1.);

    return;

}

//______________________________________________________________________________
void decayRates::FFfunctionsBGL(double w, vector<parameter> FFpar, double& hw, double& R1w, double& R2w){

  //these are the coeff of the 3 series
  double bp[10] = {0.};
  double cp[10] = {0.};
  double ap[10] = {0.};

  //link now the parameters to the coeff of the series
  int count_b=0;
  int count_c=0;
  int count_a=0;
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
      else if(p.name == "F1") {
          F1 = p.value;
      }
      else { cout << "!!! ERROR FFfunctionsBGL for B->D*, no valid parameters!" << endl;
          hw=0; R1w=1; R2w=1; return;
      };
  }
  
  double z = (sqrt(w+1.) - sqrt(2.))/(sqrt(w+1.) + sqrt(2.));
  double mB  = Mass::Bs/1000;
  double mD  = Mass::DsS/1000;
  double r = mD/mB;
  double tplus = (mB+mD)*(mB+mD);
  double tminus = (mB-mD)*(mB-mD);

  double mBv[4]={6.329,6.920,7.020,7.280};
  double mBa[4]={6.739,6.750,7.145,7.150};

  double nI=2.6;
  double chiT1m = 5.131e-4;
  double chiT1p = 3.894e-4;
  double Pvscale = 1;
  double Pascale = 1;
  //double Pvscale =  2.52733;
  //double Pascale =  2.02159;

  double b0 = 2.*sqrt(mB*mD)*Pascale*P1pm(0,tplus,tminus,mBa)*
              4.*r/mB/mB*sqrt(nI/3./TMath::Pi()/chiT1p)/
              pow( (1.+r) + 2.*sqrt(r), 4)*F1;

  double r_phi_F1f_0 = 1./mB/sqrt(2.)/((1.+r) + 2.*sqrt(r));
  double c0 = (mB-mD)* r_phi_F1f_0* b0;

  //cout << "c0 " << c0 << endl;
  double phi_f = 4.*r/mB/mB * sqrt(nI/3./TMath::Pi()/chiT1p)*
                 (1.+z)*sqrt((1.-z)*(1.-z)*(1.-z))/
                 pow( (1.+r)*(1.-z) + 2.*sqrt(r)*(1.+z), 4);

  double phi_g = sqrt(nI/3./TMath::Pi()/chiT1m)*
                  (16.*r*r*(1.+z)*(1.+z)/sqrt(1.-z))/
                  pow( (1.+r)*(1.-z) + 2.*sqrt(r)*(1.+z), 4);

  double phi_F1 = 4.*r/mB/mB/mB * sqrt(nI/6./TMath::Pi()/chiT1p)*
                  (1.+z)*sqrt((1.-z)*(1.-z)*(1.-z)*(1.-z)*(1.-z))/
                  pow( (1.+r)*(1.-z) + 2.*sqrt(r)*(1.+z), 5);

  double p1p  = Pascale*P1pm(z,tplus,tminus,mBa);

  double sumf = b0; int order = 1;
  for(auto cbp : bp){ sumf += cbp*pow(z, order); ++order;}
  double fBGL = sumf/p1p/phi_f;

  double p1m  = Pvscale*P1pm(z,tplus,tminus,mBv);
  double sumg = 0; order =0;
  for(auto cap : ap){ sumg += cap*pow(z, order); ++order;}
  double gBGL = sumg/p1m/phi_g;

  double sumF1 = c0; order =1;
  for(auto ccp : cp){ sumF1 += ccp*pow(z, order); ++order;}
  double F1BGL = sumF1/p1p/phi_F1;

  hw  = fBGL/sqrt(mB*mD)/(1.+w);
  R1w = (w+1.) *mB*mD* gBGL/fBGL;
  R2w = (w-r)/(w-1.) - F1BGL/(mB*(w-1.)*fBGL);

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

  if(pFF.model=="CLN")
      FFfunctionsCLN(w,pFF.FFpars,hw,R1w,R2w);
  else if(pFF.model=="BGL")
      FFfunctionsBGL(w,pFF.FFpars,hw,R1w,R2w);
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

  if(pFF.model=="CLN")
      FFfunctionsCLN(w,pFF.FFpars,hw,R1w,R2w);
  else if(pFF.model=="BGL")
      FFfunctionsBGL(w,pFF.FFpars,hw,R1w,R2w);
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

  //double angularTerm = 0.95*angularTermPi + 0.05*angularTermGamma;
  //double angularTerm = angularTermPi;
  //if(isBs==1)
  double angularTerm = 0.06*angularTermPi + 0.94* angularTermGamma;
  //angularTerm = angularTermGamma;
  

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

  if(pFF.model=="CLN")
      FFfunctionsCLN(w,pFF.FFpars,hw,R1w,R2w);
  else if(pFF.model=="BGL")
      FFfunctionsBGL(w,pFF.FFpars,hw,R1w,R2w);
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

  //double angularTerm = 0.95*angularTermPi + 0.05*angularTermGamma;
  //double angularTerm = angularTermPi;
  //if(isBs==1)
  double angularTerm = 0.06*angularTermPi + 0.94* angularTermGamma;
  //angularTerm = angularTermGamma;
  

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

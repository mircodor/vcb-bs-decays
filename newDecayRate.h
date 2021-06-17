#ifndef newdecayrates_h
#define newdecayrates_h

#include "Tools.h"
#include "TF1.h"

class FFModel
{
 public :
  TString model;
  unsigned int nPars;
  vector<parameter> FFpars;
  
  void initNFFpars(unsigned int i){
    nPars = i;
    FFpars.clear();
    FFpars.reserve(nPars);
    for(unsigned int jp=0; jp<nPars; ++jp){
      parameter p;
      p.inum = jp;
      p.name = "no name";
      p.value=-99999.999;
      p.error=0.;
      p.min=0;
      p.max=0;
      p.gconst=false;
      FFpars.push_back(p);
    }
  }
    
  void setNFFpars(){
    nPars = FFpars.size();
  }
  
  void print(){
    cout << " Form factor model: " << model << endl;
    for(auto p : FFpars) p.print();
  }
    
    FFModel() {}
    ~FFModel() {};
    
};



class decayRates {
private :
    FFModel _pFF;
    double  _VcbEff;
    double  _tauB;
    bool    _isDst;
    static double P1pm(double z, double tplus, double tminus, double* mBc);
    
    double dGdwD  (FFModel pFF,double VcbEff, double w);
    double dGdwDst(FFModel pFF,double VcbEff, double w);
    double dGdwdctdDst   (FFModel pFF,double VcbEff, double w, double ctd);
    double dGdwdAnglesDst(FFModel pFF,double VcbEff, double w, double ctl, double ctd, double chi);
    
    double calculateGamma(FFModel pFF, double VcbEff, bool isDst);
    double calculateBR   (FFModel pFF, double VcbEff, double tauB, bool isDst);
    double pperp_w_cos    (double w, double cos, double massD);
    double pperp_w_ctd_cos(double beta, double ctd, double cos, double M, double m1, double m2);
    static double Em1(double M, double m1, double m2);
    static double Pm1(double M, double m1, double m2);
    
public :
	static void FFfunctionsCLN(double w, vector<parameter> FFpar, double& Gw);
	static void FFfunctionsCLN(double w, vector<parameter> FFpar, double& hw, double& R1w, double& R2w);
	static void FFfunctionsBGL(double w, vector<parameter> FFpar, double& Gw);
	static void FFfunctionsBGL(double w, vector<parameter> FFpar, double& A1, double& A2, double& V);
  
    double TF_dGdw      (double* w, double *) { return _isDst? dGdwDst(_pFF,_VcbEff,w[0]) : dGdwD(_pFF,_VcbEff,w[0]) ; };
    double TF_dGdwdctd  (double* w, double *) { return _isDst? dGdwdctdDst(_pFF,_VcbEff,w[0],w[1]) : 0; };
    double TF_dGdAngle  (double* w, double *p) { return _isDst? dGdwdAnglesDst(_pFF,_VcbEff,p[0],w[0],w[1],w[2]) : 0; };
   
    double Eval_dGdw      (double w) { return _isDst? dGdwDst(_pFF,_VcbEff,w) : dGdwD(_pFF,_VcbEff,w); };
    double Eval_dGdwdctd  (double w, double ctd){ return _isDst? dGdwdctdDst(_pFF,_VcbEff,w,ctd) : 0; };
    double Eval_dGdwdAngle(double w, double ctl, double ctd, double chi){ return _isDst? dGdwdAnglesDst(_pFF,_VcbEff,w,ctl,ctd,chi) : 0; };
    
    double Eval_Gamma(){ return calculateGamma(_pFF,_VcbEff,_isDst);    };
    double Eval_BR()   { return calculateBR(_pFF,_VcbEff,_tauB,_isDst); };
    double TF_Pperp(double* w, double *p){ return _isDst? pperp_w_ctd_cos(w[0], w[1], w[2], p[0], p[1], p[2]) :
            pperp_w_cos(w[0], w[1], p[0]); };
    
    double Eval_Pperp(double w, double cos, double M){ return pperp_w_cos(w,cos,M); };
    double Eval_Pperp(double beta, double ctd, double cos, double M, double m1, double m2){ return _isDst?
        pperp_w_ctd_cos(beta, ctd, cos, M, m1, m2) : 0; };
    
    
    void SetFFModel(FFModel pFF)  { _pFF = pFF;      };
    void SetVcbEff (double VcbEff){ _VcbEff = VcbEff;};
    void SetTauB   (double tauB)  { _tauB = tauB;    };
    void SetIsDst  (bool isDst)   { _isDst = isDst;  };
    
    FFModel GetFFModel(){ return _pFF ;  };
    double GetVcbEff()  { return _VcbEff;};
    double GetTauB()    { return _tauB;  };
    bool   GetIsDst()   { return _isDst; };
    
    decayRates(FFModel pFF, double VcbEff, double tauB, bool isDst) :
        _pFF{pFF}, _VcbEff{VcbEff}, _tauB{tauB}, _isDst{isDst} {}
    
    decayRates(FFModel pFF, double VcbEff, bool isDst) :
    _pFF{pFF}, _VcbEff{VcbEff}, _tauB{TauPs::Bs}, _isDst{isDst} {}
    
    decayRates(FFModel pFF, bool isDst) :
        _pFF{pFF}, _VcbEff{Constants::VcbRef*Constants::etaEW}, _tauB{TauPs::Bs}, _isDst{isDst} {}
    
    
    ~decayRates();

};

#endif

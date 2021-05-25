#ifndef TOOLS_H
#define TOOLS_H
#include "TMath.h"
#include "Riostream.h"
#include "TString.h"
#include <vector>

class parameter{
    
 public:
  unsigned int inum;
  TString  name;
  double   value;
  double   error;
  double   min;
  double   max;
  bool     gconst;
    
  void     print() {
    std::cout << "\t"
          <<  inum  << "\t" << name  << "\t"
          <<  value << "\t" << error << "\t"
          <<  min   << "\t" << max   << "\t"
          << gconst << std::endl;
  }
};


// -----------------------------------------
// Useful constants
// -----------------------------------------
namespace Constants
{
  const double C         = TMath::C()*1e-9;
  const double GF        = 1.1663787e2; //GeV-2 (actually mutiplied by e7)
  const double hbar      = 65.8211928; //GeV * ps (actually in unit of e-14)
  const double VcbRef    = 39.25e-3; //HFLAV combined exclusive, error 0.56e-3
  const double etaEW     = 1.0066; //error 0.0055
}


// -----------------------------------------
// Particle masses (MeV/c2)
// -----------------------------------------
namespace Mass
{
  const double E   = 0.510998910;
  const double Mu  = 105.658367;
  const double Tau = 1776.84;
  const double K   = 493.677;
  const double Pi  = 139.57018;
  const double Pi0 = 134.9766;
  const double D0  = 1864.84;
  const double DS  = 2010.27;
  const double D   = 1869.62;
  const double Ds  = 1968.49;
  const double DsS = 2112.1;
  const double P   = 938.272013;
  const double Phi = 1019.455;
  const double B   = 5279.17;
  const double Bd  = 5279.5;
  const double Bs  = 5366.3;
  const double K0  = 4976.14;
  const double Lm  = 1115.683;
}


// -----------------------------------------
// Particle lifetimes (ps) as of PDG 2015
// -----------------------------------------
namespace TauPs
{
  const double D0  = 0.4101;
  const double D   = 1.040;
  const double Ds  = 0.500;
  const double B   = 1.638;
  const double Bd  = 1.519;
  const double Bs  = 1.510;
  const double K0s = 89.54;
  const double K0l = 5116.;
}

// -----------------------------------------
// Errors on particle lifetimes (ps) as of PDG 2015
// -----------------------------------------
namespace ErrTauPs
{
  const double D0  = 0.0015;
  const double D   = 0.007;
  const double Ds  = 0.007;
  const double B   = 0.004;
  const double Bd  = 0.004;
  const double Bs  = 0.005;
  const double K0s = 0.040;
  const double K0l = 21.00;
}


// -----------------------------------------
// Particle MC Codes
// -----------------------------------------
namespace PdgCode
{
  int Gamma= 22;
  int E    = 11;
  int NuE  = 12;
  int Mu   = 13;
  int NuMu = 14;
  int Tau  = 15;
  int NuTau= 16;
  int K    = 321;
  int Pi   = 211;
  int Pi0  = 111;
  int D0   = 421;
  int D    = 411;
  int DS   = 413;
  int Ds   = 431;
  int DsS  = 433;
  int P    = 2212;
  int Phi  = 333;
  int B    = 521;
  int Bd   = 511;
  int Bs   = 531;
  int Bc   = 541;
  int K0   = 311;
  int K0s  = 310;
  int K0l  = 130;
  int Lm   = 3122;
  int Lc   = 4122;
  int Lb   = 5122;
}




#endif // TOOLS_H



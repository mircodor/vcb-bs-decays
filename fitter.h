#ifndef fitter_h                                                                        
#define fitter_h

#include "Tools.h"
#include "decayRates.h"

#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLegend.h"
#include "TMatrixD.h"

class MCcand{
 public:
  double wvar;
  double ctl;
  double ctd;
  double chi;
  double pperp;
  double den;
  int signalCategory;
    
    void print() {
        std::cout
            << " MC cand, signal type " << signalCategory << "\n"
            << " w     = " << wvar << "\n"
            << " ctl   = " << ctl  << "\n"
            << " ctd   = " << ctd  << "\n"
            << " chi   = " << chi  << "\n"
            << " pperp = " << pperp  << "\n"
            << " den   = " << den  << "\n"
            << "---------------------\n"
            << std::endl;
    }
};

class pdfComp         
{                 
 public:             
  int      id;            
  double   yield;       
  TH1D*    pdf;
  TString  name;      
  double   color;         
  TString  legend;      
  void     print() {    
    cout << " Component " << id << "\t" << name << "\t" << legend << endl;
  }                                                                        
};                          

vector<MCcand>    mccand;
vector<parameter> otherPars;
vector<parameter> fitPars;
vector<parameter> parConst;
vector<vector<double> > covParConst(100,vector<double>(100));

decayRates* decFitDs;
decayRates* decFitDsS;
decayRates* decRefDs;
decayRates* decRefDsS;

pdfComp fitComp[4];


std::vector<TString> theoryInputs;


double _Bmass, _Dsmass, _DsSmass;
double _chi2, _ndf;
TH1D * hData;
TH1D* hacc[2];

bool _rate1D = false;

void   SetAllPars(double* );
void   SetAllPars();
void   FillHistogram();
void   fcn_tot(int &, double *, double &, double *, int );
void   calculateYields();

class fitter
{
 private:
  bool Configure(TString filename);
  bool DoFit(double strategy=2, bool useHesse=true, bool useMinos=false);
  void PrintConfigInfo();
  bool FillMCcandidates(TString filename);
  bool Plot(TCanvas *);
  bool ReadConfigFile(TString filename);
  void SetDataAndBkg();

 public: 
  bool RunFitter(TString filename);
    
    fitter() {}
  
  ~fitter() {
    fitPars.clear();
  };

};

#endif 

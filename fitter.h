#ifndef fitter_h                                                                        
#define fitter_h

#include "Tools.h"
#include "UpdatedDecayRate.h"

#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLegend.h"
#include "TMatrixD.h"
#include "TGraphErrors.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TDecompChol.h"

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


class extInputs {                                                                                                                                 
 public:
  int dim;
  std::vector<TString> name;
  std::vector<double>  x;
  std::vector<double>  x_err;
  std::vector<double>  val;
  std::vector<double>  val_err;
  map<TString,double>  corr;
  map<TString,double>  cov;
  TMatrixD corr_matrix;                        
  TMatrixD cov_matrix;                            
  double cov_inv_matrix[30][30];                   
  void print() {
	cout << " dimension: " << dim << endl;
	for(int i =0; i<dim; ++i)
	  cout << "\t" << name[i] << "\t" << x[i] << "\t" << x_err[i] << "\t" << val[i] << "\t" << val_err[i] << endl;
	cout << " covariance matrix " << endl;
	cov_matrix.Print();
  }

  extInputs() {};                                                           
  ~extInputs() {};                      
};                                



vector<MCcand>    mccand;
vector<parameter> otherPars;
vector<parameter> fitPars;
vector<parameter> parConst;
vector<vector<double> > covParConst(100,vector<double>(100));
double fitFullMatrix[100][100];

decayRates* decFitDs;
decayRates* decFitDsS;
decayRates* decRefDs;
decayRates* decRefDsS;

pdfComp fitComp[4];
pdfComp fitBand[4];

std::vector<TString> theoryInputs;

double _Bmass{Mass::Bs}, _Dsmass{Mass::Ds}, _DsSmass{Mass::DsS};
double _chi2, _ndf;
double _normHistLHCb;
TH1D * hData;
TH1D* hacc[2];
std::map<TString,extInputs> _extInputs;

const bool _rate1D = false;
const TString dirInputs = "./new_inputs_modified/";

void   SetAllPars(double* );
void   SetAllPars();
void   FillHistogram();
void   fcn_tot(int &, double *, double &, double *, int );
void   calculateYields();
double FFfunctions(TString xname, double x, double xerr);

class fitter
{
 private:
  bool Configure(TString filename);
  bool DoFit(double strategy=2, bool useHesse=true, bool useMinos=false);
  void PrintConfigInfo();
  bool FillMCcandidates(TString filename);
  bool DoProjection();
  bool ReadConfigFile(TString filename);
  bool SetDataAndBkg();
  bool SetExtInputs();
  void DrawResiduals(TH1D* hD, TH1D* hF, TH1D* hp);
  void DrawFFErrorBand(TString xname, std::vector<TGraphErrors*>& gr);
 public: 
  bool RunFitter(TString filename);
    
    fitter() {}
  
  ~fitter() {
    fitPars.clear();
  };

};

#endif 

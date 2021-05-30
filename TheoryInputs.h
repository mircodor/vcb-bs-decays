#include "TMatrixD.h"
#include "decayRates.h"
#include "Tools.h"

std::vector<double> w;
std::vector<double> wlow;
std::vector<double> werr;
std::vector<double> f;
std::vector<double> ferr;
double covParDs[30][30];

std::vector<double> wHPQCD;
std::vector<double> fHPQCD;
std::vector<double> fHPQCDerr;

std::vector<double> wMILC;
std::vector<double> fMILC;
std::vector<double> fMILCerr;

std::vector<double> wLCSRDs;
std::vector<double> fLCSRDs;
std::vector<double> fLCSRDserr;

std::vector<double> wLCSRDsS;
std::vector<double> fLCSRDsS;
std::vector<double> fLCSRDsSerr;
std::vector<double> V;
std::vector<double> A1;
std::vector<double> A2;
std::vector<double> Verr;
std::vector<double> A1err;
std::vector<double> A2err;
double covParDsS[30][30];


std::vector<double> wLHCb;
std::vector<double> wLHCblow;
std::vector<double> wLHCberr;
std::vector<double> fLHCb;
std::vector<double> fLHCberr;
double covParLHCb[30][30];

void readLHCbData(TString whichData = "LHCb") {

  w.clear(); wlow.clear(); werr.clear();
  f.clear(); ferr.clear();
  
  TFile * fin = TFile::Open("inputs/"+whichData+"_data.root");
  TH1D * h = (TH1D*)fin->Get(whichData+"_data");
  for(int i = 1; i <= h->GetNbinsX(); i++) {
    w.push_back(h->GetBinCenter(i));
    wlow.push_back(h->GetBinLowEdge(i));
    werr.push_back(h->GetBinWidth(i));
    f.push_back(h->GetBinContent(i));
    ferr.push_back(h->GetBinError(i));
  }
  //fin->Close();
  TFile * finCov = TFile::Open("inputs/"+whichData+"_CovInv.root");
  TH2D * hCovInv = (TH2D*)finCov->Get("hLHCbCovInv");
  //finCov->Close();
  
  for(int i = 1; i <= hCovInv->GetNbinsX(); i++)
    for(int j = 1; j <= hCovInv->GetNbinsY(); j++) {
      covParLHCb[i-1][j-1] = hCovInv->GetBinContent(i,j);
      //cout << i-1 << " " << j-1 << " "<< covParLHCb[i-1][j-1] << endl; 
    }

  if(whichData == "LHCb") {
    wLHCb = w;
    wLHCblow = wlow;
    wLHCberr = werr;
    fLHCb = f;
    fLHCberr = ferr;
  }
    
  
}


//____________________________________________________________
void readTheoryDataDsS(TString whichData = "LCSRDsS") {


  std::vector<TString> parName;
  if(whichData == "LCSRDsS") {
    parName.push_back("vs");
    parName.push_back("a1s");
    parName.push_back("a2s");
  }

  w.clear(); f.clear(); ferr.clear();
  V.clear(); A1.clear(); A2.clear();

  TString name;
  double corr, wi, fi, ferri;
  std::map<TString, double> corrs;

  //read w points, FF and uncertainty
  std::ifstream infile ( "inputs/"+whichData+"_values.txt" );
  while(infile >> name >> wi >> fi >> ferri) {
    w.push_back(wi);
    f.push_back(fi);
    ferr.push_back(ferri);          
    if(name.Contains("vs")) {
      V.push_back(fi);
      Verr.push_back(ferri);
    }
    if(name.Contains("a1s")) {
      A1.push_back(fi);
      A1err.push_back(ferri);
    }
    if(name.Contains("a2s")) {
      A2.push_back(fi);
      A2err.push_back(ferri);
    }
  }

  infile.close();
  
  ///read correlation matrix
  TMatrixD corrMatrix(w.size(), w.size());
  infile.open( "inputs/"+whichData+"_corr.txt" );

  while(infile >> name >> corr) 
    corrs.insert(std::make_pair(name,corr));
  
  infile.close();
  
  for(long unsigned int i = 1; i <= w.size(); i++){
    for(long unsigned int j = 1; j <= w.size(); j++) {
      unsigned long int itmp = i, jtmp = j, par1 = 0, par2 = 0;
      if(i>4) {itmp = i-4; par1 = 1;}
      if(j>4) {jtmp = j-4; par2 = 1;}
      if(i>8) {itmp = i-8; par1 = 2;}
      if(j>8) {jtmp = j-8; par2 = 2;}

      if(i<j)  {
	double c = corrs[Form("%s%lu_%s%lu",parName[par1].Data(),itmp,parName[par2].Data(),jtmp)];
	if(c==0) c = corrs[Form("%s%lu_%s%lu",parName[par2].Data(),jtmp,parName[par1].Data(),itmp)];
	corrMatrix(i-1,j-1) = c;
      }
      else if(i>j)
	corrMatrix(i-1,j-1) = corrMatrix(j-1,i-1);
      else if(i==j)
	corrMatrix(i-1,j-1) = 1;
    }
  }
  //cout << "CORRELATION MATRIX "+whichData << endl;
  //corrMatrix.Print();

  ///build covariance matrix
  TMatrixD covMatrix(w.size(), w.size());
  for(long unsigned int i = 0; i < w.size(); i++){
    for(long unsigned int j = 0; j < w.size(); j++) {
      covMatrix(i,j) = corrMatrix(i,j)*ferr[i]*ferr[j];
    }
  }

  //cout << "COVARIANCE MATRIX LCSRDsS" << endl;
  //covMatrix.Print();


  TMatrixD covMatrixInv(w.size(), w.size());
  covMatrixInv = covMatrix.InvertFast();

  //cout << "INVERSE COVARIANCE MATRIX LCSRDsS" << endl;
  //covMatrixInv.Print();

  for(long unsigned int i = 0; i < w.size(); i++){
    for(long unsigned int j = 0; j < w.size(); j++) {
      covParDsS[i][j] = covMatrixInv(i,j);
    }
  }

  if(whichData == "LCSRDsS") {
    wLCSRDsS = w;
    fLCSRDsS = f;
    fLCSRDsSerr = ferr;
  }

}

//____________________________________________________________
void readTheoryDataDs(TString whichData= "HPQCD") {

  ///some setup and declarations
  TString parName;
  if(whichData == "MILC")  parName = "fMILCs";
  if(whichData == "HPQCD") parName = "fHPs";
  if(whichData == "LCSRDs")  parName = "ffps";


  w.clear(); f.clear(); ferr.clear();

  TString name;
  double corr, wi, fi, ferri;
  std::map<TString,double> corrs;


  //read w points, FF and uncertainty
  std::ifstream infile ( "inputs/"+whichData+"_values.txt" );
  while(infile >> wi >> fi >> ferri) {
    w.push_back(wi);
    f.push_back(fi);
    ferr.push_back(ferri);      
    //cout << whichData << " " << wi << " " << fi << " " << endl; 
  }


  infile.close();
  
  ///read correlation matrix
  TMatrixD corrMatrix(w.size(), w.size());
  infile.open( "inputs/"+whichData+"_corr.txt" );

  while(infile >> name >> corr) 
    corrs.insert(std::make_pair(name,corr));
  
  infile.close();
  
  for(long unsigned int i = 1; i <= w.size(); i++){
    for(long unsigned int j = 1; j <= w.size(); j++) {
      if(i<j) 
	corrMatrix(i-1,j-1) = corrs[Form("%s%lu_%s%lu",parName.Data(),i,parName.Data(),j)];
      else if(i>j)
	corrMatrix(i-1,j-1) = corrMatrix(j-1,i-1);
      else if(i==j)
	corrMatrix(i-1,j-1) = 1;

    }
  }
  //cout << "CORRELATION MATRIX "+whichData << endl;
  //corrMatrix.Print();

  ///build covariance matrix
  TMatrixD covMatrix(w.size(), w.size());
  for(long unsigned int i = 0; i < w.size(); i++){
    for(long unsigned int j = 0; j < w.size(); j++) {
      covMatrix(i,j) = corrMatrix(i,j)*ferr[i]*ferr[j];
    }
  }

  //cout << "COVARIANCE MATRIX MILC" << endl;
  //covMatrix.Print();


  TMatrixD covMatrixInv(w.size(), w.size());
  covMatrixInv = covMatrix.InvertFast();

  //  cout << "INVERSE COVARIANCE MATRIX MILC" << endl;
  //covMatrixInv.Print();

  for(long unsigned int i = 0; i < w.size(); i++){
    for(long unsigned int j = 0; j < w.size(); j++) {
      covParDs[i][j] = covMatrixInv(i,j);
    }
  }

  if(whichData == "HPQCD") {
    wHPQCD = w;
    fHPQCD = f;
    fHPQCDerr = ferr;
  }

  if(whichData == "MILC") {
    wMILC = w;
    fMILC = f;
    fMILCerr = ferr;
  }

  if(whichData == "LCSRDs") {
    wLCSRDs = w;
    fLCSRDs = f;
    fLCSRDserr = ferr;
  }

}

//____________________________________________________________
double FFfunctionsCLN(double w, vector<parameter> FFpar){       
  double rho2 = 0;                
  double G0 = 0;                         
  double r     = Mass::Ds/Mass::Bs;          

  for(auto p : FFpar){               
    if     (p.name.Contains("rho2P")) rho2 = p.value;          
    if     (p.name.Contains("G0"))    G0   = p.value;
  }                  
  double z = (sqrt(w+1.) - sqrt(2.))/(sqrt(w+1.) + sqrt(2.));   
  double Gw = G0*(1. - 8.*rho2*z + (51.*rho2 - 10.)*z*z - (252.*rho2 - 84.)*z*z*z);      
  
  return (1+r)/sqrt(4*r)*Gw;     
} 

//____________________________________________________________
void FFfunctionsCLN_DsS(double w, vector<parameter> FFpar, double& V, double& A1, double& A2) {
  double rho2 = 0;
  double R1 = 0;
  double R2 = 0;
  double F1 = 0;
  double mB = Mass::Bs/1000;
  double mDsS = Mass::DsS/1000;
  double r = mDsS/mB;
  double q2 = mB*mB + mDsS*mDsS - 2*w*mB*mDsS;

  for(auto p : FFpar){               
    if     (p.name.Contains("rho2V")) rho2 = p.value;          
    if     (p.name.Contains("R1"))    R1   = p.value;          
    if     (p.name.Contains("R2"))    R2   = p.value;          
    if     (p.name.Contains("F1"))    F1   = p.value;          
  }                  
  double z = (sqrt(w+1.) - sqrt(2.))/(sqrt(w+1.) + sqrt(2.));   

  double hw = F1*(1. - 8.*rho2*z + (53.*rho2 - 15.)*z*z - (231.*rho2 - 91.)*z*z*z);
  double R1w =  R1 - 0.12*(w-1.) + 0.05*(w-1.)*(w-1.);  
  double R2w =  R2 + 0.11*(w-1.) - 0.06*(w-1.)*(w-1.);  
  
  double fw = hw * sqrt(mB*mDsS)*(w+1.);
  double gw = R1w*fw/(w+1)/mB/mDsS;
  double lambda = pow(mB,4) + pow(mDsS,4) + pow(q2,2) - 2.0*(mB*mB*mDsS*mDsS + mB*mB*q2 + mDsS*mDsS*q2);
  double F1w = ( (w-r)/(w-1)-R2w ) * mB*(w-1)*fw;

  V = (mB+mDsS)/2.0*gw;
  A1 = 1/(mB+mDsS)*fw;
  A2 = (mB+mDsS)/lambda*( (mB*mB-mDsS*mDsS-q2)*fw - 2.0*mDsS*F1w );
}

//____________________________________________________________

double dGdwDst(double * x, double * par) {
  double w    = x[0];              
  double F1   = par[0];                
  double rho2 = par[1];                
  double R1   = par[2];             
  double R2   = par[3];               
  double normLHCb = par[4];
  double mB = Mass::Bs/1000;
  double mDsS = Mass::DsS/1000;
  double r = mDsS/mB;

  double z = (sqrt(w+1)-sqrt(2))/(sqrt(w+1)+sqrt(2));
  double hw  = F1*(1. - 8.*rho2*z + (53.*rho2 - 15.)*z*z - (231.*rho2 - 91.)*z*z*z) ;
  double R1w = R1 - 0.12*(w-1.) + 0.05*(w-1.)*(w-1.);       
  double R2w = R2 + 0.11*(w-1.) - 0.06*(w-1.)*(w-1.);           
  double XiFw2 =  hw*hw * sqrt(w*w-1.) *(w+1.)*(w+1.) * ( 2. * (1. - 2.*w*r + r*r)/(1.-r)/(1.-r) * (1. + R1w*R1w *(w -1.)/(w+1.) )  +  ( 1. + (1. - R2w)*(w-1.)/(1.-r)) * ( 1. + (1. - R2w)*(w-1.)/(1.-r)) );      
  
  return XiFw2*normLHCb;

}

void addTheoryInputDs(TString whichData, vector<parameter> FFpar, double& _chi2, double& _ndf){

  readTheoryDataDs(whichData);
  
  for(long unsigned int i = 0; i < w.size(); i++) {
    for(long unsigned int j = 0; j < w.size(); j++) {
      if(whichData == "MILC") {
	if(i!=j) continue;
	_chi2 += (FFfunctionsCLN(w[i],FFpar)-f[i])*(FFfunctionsCLN(w[j],FFpar)-f[j])/ferr[i]/ferr[j];
	//_chi2 += (FFfunctionsCLN(w[i],FFpar)-f[i])*(FFfunctionsCLN(w[j],FFpar)-f[j])*covParDs[i][j];
	//cout << i << " " << j << " " << FFfunctionsCLN(w[i],FFpar)-f[i] << " " << FFfunctionsCLN(w[j],FFpar)-f[j] << " " << covPar[i][j] << " " << _chi2 << endl;
      }
      if(whichData == "HPQCD") {
	//if(i!=j) continue;
	//_chi2 += (FFfunctionsCLN(w[i],FFpar)-f[i])*(FFfunctionsCLN(w[j],FFpar)-f[j])/ferr[i]/ferr[j];
	_chi2 += (FFfunctionsCLN(w[i],FFpar)-f[i])*(FFfunctionsCLN(w[j],FFpar)-f[j])*covParDs[i][j];
      }
      
      //cout << i << " " << j << " " << FFfunctionsCLN(w[i],FFpar) << " " << f[i] << " " << ferr[i] << " " << _chi2 << endl;

      }
  }
  _ndf+=w.size();


}

//____________________________________________________________
void addTheoryInputDsS(TString whichData, vector<parameter> FFpar, double& _chi2, double& _ndf){

  readTheoryDataDsS(whichData);
  
  double Vwi, A1wi, A2wi, Vwj, A1wj, A2wj;
  for(long unsigned int i = 0; i < w.size(); i++) {
    for(long unsigned int j = 0; j < w.size(); j++) {
      
      //if(i > 7 || j > 7) continue;
      FFfunctionsCLN_DsS(w[i],FFpar,Vwi,A1wi,A2wi);
      FFfunctionsCLN_DsS(w[j],FFpar,Vwj,A1wj,A2wj);
      double di, dj, ffi, ffj;
      if(i <= 3) { ffi = Vwi; }
      else if( i > 3 && i <= 7) {ffi = A1wi;}
      else if( i > 7) {ffi = A2wi;}
      if(j <= 3) {ffj = Vwj;}
      else if( j > 3 && j <= 7) {ffj = A1wj;}
      else if( j > 7) {ffj = A2wj;}
      di = f[i];
      dj = f[j];



      
      //if(i!=j) continue;
      //_chi2 += (ffi-di)*(ffj-dj)/ferr[i]/ferr[j];

      _chi2 += (ffi-di)*(ffj-dj)*covParDsS[i][j];
      //cout << i << " " << j  << " " << w[i] << " " << w[j] << " " << ffi << " " << di << " " << ffj << " " <<  dj << " " << ferr[i] << " " << ferr[j] << " " << covParDsS[i][j] << " " << _chi2 << endl;      
    }
  }
  _ndf+=w.size();
}


//____________________________________________________________
//Aggiungere altro vector di paramtri come argomenti (quelli Other) e aggiungere un parametro di norm nel config
//poi passare questo parametro alla funzione come fattore di normalizzazione al posto di 0.65
void addLHCbInputDsS(TString whichData, vector<parameter> FFpar, vector<parameter> otherPars, double& _chi2, double& _ndf){ 

  readLHCbData();

  TF1 * fff = new TF1("fff", dGdwDst, 1, 1.4667, 5);
  double rho2(0), R1(0), R2(0), F1(0), normLHCb(0);
  for(auto p : FFpar) {               
    if     (p.name.Contains("F1"))    F1   = p.value;          
    if     (p.name.Contains("rho2V")) rho2 = p.value;          
    if     (p.name.Contains("R1"))    R1   = p.value;          
    if     (p.name.Contains("R2"))    R2   = p.value;          
  }
  for(auto p : otherPars) {
    if     (p.name.Contains("normLHCb"))  normLHCb   = p.value;          
  }

  fff->SetParameter(0,F1);
  fff->SetParameter(1,rho2);
  fff->SetParameter(2,R1);
  fff->SetParameter(3,R2);
  fff->SetParameter(4,normLHCb);
  //cout << F1 << " " << rho2 << " " << R1 << " " << R2 << " " << normLHCb << endl;
  
  for(long unsigned int i = 0; i < w.size(); i++) {           
    for(long unsigned int j = 0; j < w.size(); j++) {

      double ffi = fff->Integral(wlow[i],wlow[i]+werr[i])/werr[i];
      double ffj = fff->Integral(wlow[j],wlow[j]+werr[j])/werr[j];
      //if(i != j)continue;
      //_chi2 += (ffi-f[i])*(ffj-f[j])/ferr[i]/ferr[j];

      _chi2 += (ffi-f[i])*(ffj-f[j])/covParLHCb[i][j];
      //cout << i << " " << j << " " << ffi << " " << f[i] << " " << ffj << " " << f[j] << " " << ferr[i] << " " << ferr[j] << " " << 1./(ferr[i]*ferr[j]) << " " << covParLHCb[i][j] << " " << _chi2 << endl;
    }
  }
  
  _ndf+=w.size();
  
}

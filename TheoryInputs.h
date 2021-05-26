#include "TMatrixD.h"
#include "decayRates.h"
#include "Tools.h"

std::vector<double> w;
std::vector<double> f;
std::vector<double> ferr;
double covPar[30][30];

std::vector<double> wHPQCD;
std::vector<double> fHPQCD;
std::vector<double> fHPQCDerr;
double covParHPQCD[30][30];

std::vector<double> wMILC;
std::vector<double> fMILC;
std::vector<double> fMILCerr;
double covParMILC[30][30];


void readTheoryDataDs(TString whichData= "HPQCD") {

  ///some setup and declarations
  TString parName;
  if(whichData == "MILC")  parName  = "fMILCs";
  if(whichData == "HPQCD") parName  = "fHPs";
  if(whichData == "LCSR")  parName  = "ffps";

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

  //cout << "INVERSE COVARIANCE MATRIX MILC" << endl;
  //covMatrixInv.Print();

  for(long unsigned int i = 0; i < w.size(); i++){
    for(long unsigned int j = 0; j < w.size(); j++) {
      covPar[i][j] = covMatrixInv(i,j);
    }
  }

  if(whichData == "HPQCD") {
    wHPQCD = w;
    fHPQCD = f;
    fHPQCDerr = ferr;
    for(long unsigned int i = 0; i < w.size(); i++){               
      for(long unsigned int j = 0; j < w.size(); j++) {                  
	covParHPQCD[i][j] = covPar[i][j];             
      }                    
    }                
  }

  if(whichData == "MILC") {
    wMILC = w;
    fMILC = f;
    fMILCerr = ferr;
    for(long unsigned int i = 0; i < w.size(); i++){               
      for(long unsigned int j = 0; j < w.size(); j++) {                  
	covParMILC[i][j] = covPar[i][j];             
      }                    
    }                
  }
}

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

void addTheoryInputDs(TString whichData, vector<parameter> FFpar, double& _chi2, double& _ndf){

  readTheoryDataDs(whichData);
  
  for(long unsigned int i = 0; i < w.size(); i++) {
    for(long unsigned int j = 0; j < w.size(); j++) {
      if(whichData == "MILC") {
	//if(i!=j) continue;
	//_chi2 += (FFfunctionsCLN(w[i],FFpar)-f[i])*(FFfunctionsCLN(w[j],FFpar)-f[j])/ferr[i]/ferr[j];
	_chi2 += (FFfunctionsCLN(w[i],FFpar)-f[i])*(FFfunctionsCLN(w[j],FFpar)-f[j])*covPar[i][j];
      }
      else if(whichData == "HPQCD") {
	//if(i!=j) continue;
	//_chi2 += (FFfunctionsCLN(w[i],FFpar)-f[i])*(FFfunctionsCLN(w[j],FFpar)-f[j])/ferr[i]/ferr[j];
	_chi2 += (FFfunctionsCLN(w[i],FFpar)-f[i])*(FFfunctionsCLN(w[j],FFpar)-f[j])*covPar[i][j];
      }
      
      //cout << i << " " << j << " " << FFfunctionsCLN(w[i],FFpar) << " " << f[i] << " " << ferr[i] << " " << _chi2 << endl;
      //cout << i << " " << j << " " << FFfunctionsCLN(w[i],FFpar)-f[i] << " " << FFfunctionsCLN(w[j],FFpar)-f[j] << " " << covPar[i][j] << " " << _chi2 << endl;
      }
  }
  _ndf+=w.size();


}

#include "fitter.h"
#include "plotStyle.C"
#include <TVirtualFitter.h>                        
#include <TFitter.h>  
#include <TMinuit.h>
#include <TRandom3.h>
#include <TGraphErrors.h>
#include <TF2.h>
#include <TF3.h>
#include <TMultiGraph.h>
#include <map>

using namespace std;

bool fitter::RunFitter(TString filename)  {
  
    if (!Configure(filename)) return false;
	if (!DoFit(1)) return false;
    if (!DoProjection()) return false;
    return true;
  
}


bool fitter::Configure(TString filename) {
    
  cout << "==========================================================\n";
  cout << " Configuration file read: \n";
  cout << "----------------------------------------------------------\n";

  //Take data and all histograms for the fit (hardcoded)
  if(!SetDataAndBkg()) return false;

  //Read the configuration file to set the fit
  if(!ReadConfigFile(filename)) return false;
    
  //Set external inputs
  if(!SetExtInputs()) return false;
  
  //Print all information read from the config
  PrintConfigInfo();
    
  //Store the candidates for the templates
  if(!FillMCcandidates(filename)) return false;
    
  return true;
    
}

bool fitter::SetDataAndBkg(){
       
    TFile * fileData = TFile::Open("LHCb_input_data_mc.root");
	if(!fileData){
	  cout << "File to take p_perp data to fit not found!" << endl;
	  return false;
	}
    hData   = (TH1D*) fileData->Get("h_lhcb_data");
    fitComp[2].pdf = (TH1D*) fileData->Get("h_phys_bkg");
    fitComp[3].pdf = (TH1D*) fileData->Get("h_comb_bkg");
    
    hData->SetDirectory(0);
    fitComp[2].pdf->SetDirectory(0);
    fitComp[3].pdf->SetDirectory(0);
 
    fileData->Close();
   
    fitComp[2].yield = fitComp[2].pdf->Integral();
    fitComp[3].yield = fitComp[3].pdf->Integral();
    
    return true;
}


bool fitter::ReadConfigFile(TString filename){
    
  string line;
  unsigned int np=1;
  map< pair<TString,TString>, double > parCorr;
  ifstream fin(filename);
  while(fin>>line) {
    if(line=="FF-Model-Ref-P")  {
      FFModel FF;
      fin >> FF.model  >> np;
      for (unsigned int i=0; i<np; ++i) {
        parameter p;
        fin >> p.name >> p.value;
        FF.FFpars.push_back(p);
      }
        decRefDs = new decayRates(FF,false);
    }
    else if(line=="FF-Model-Ref-V"){
      FFModel FF;
      fin >> FF.model  >> np;
      for (unsigned int i=0; i<np; ++i) {
        parameter p;
        fin >> p.name >> p.value;
        FF.FFpars.push_back(p);
      }
        decRefDsS = new decayRates(FF,true);
    }
    else if(line=="FF-Model-Fit-P") {
      FFModel FF;
      fin >> FF.model >> np;
      for (unsigned int i=0; i<np; ++i) {
          parameter p;
          p.inum=i;
          fin >> p.name >> p.value >> p.error >> p.min >> p.max >> p.gconst;
          FF.FFpars.push_back(p);
      }
      decFitDs = new decayRates(FF,false);
    }
    else if(line=="FF-Model-Fit-V") {
      FFModel FF;
      fin >> FF.model >> np;
      for (unsigned int i=0; i<np; ++i) {
          parameter p;
          p.inum=i;
          fin >> p.name >> p.value >> p.error >> p.min >> p.max >> p.gconst;
          FF.FFpars.push_back(p);
      }
        decFitDsS = new decayRates(FF,true);
    }
    else if(line=="Other-Parameters") {
      fin >> np;
      for (unsigned int i=0; i<np; ++i) {
        parameter p;
        p.inum=i;
        fin >> p.name >> p.value >> p.error >> p.min >> p.max >> p.gconst;
        otherPars.push_back(p);
      }
    }
    else if(line=="Correlations") {
      fin >> np;
      for (unsigned int i=0; i<np; ++i) {
          TString par1, par2;
          double value;
          fin >> par1 >> par2 >> value;
          parCorr[make_pair(par1, par2)] = value;
          parCorr[make_pair(par2, par1)] = value;
      }
    }
    else if(line=="Theory-Inputs") {
      fin >> np;
      for (unsigned int i=0; i<np; ++i) {
		TString model;
		fin >> model;
		//cout << model << endl;
		theoryInputs.push_back(model);
      }
    }
    else continue;
  }
  fin.close();
    
  //Store the different parameters in a single vector, to be used for the fit.
  SetAllPars();
    
  unsigned int nparConst=0;
  for(auto p : fitPars) {
        if(!p.gconst) continue;
        parameter pg;
        pg = p;
        pg.inum = nparConst;
        parConst.push_back(pg);
        ++nparConst;
  }
  if(parConst.size() == 0) return true;
     
  TMatrixD covMatrixStat(nparConst,nparConst);
  for(unsigned int i=0; i<nparConst; ++i){
     for(unsigned int j=0; j<nparConst; ++j){
          auto it = parCorr.find(make_pair(parConst[i].name,parConst[j].name));
          if(it != parCorr.end()) {
              double rho = it->second;
              covMatrixStat(i,j) = rho*parConst[i].error*parConst[j].error;
          }
          else covMatrixStat(i,j) = 0;
          if(i==j) covMatrixStat(i,j) = parConst[i].error*parConst[j].error;
      }
   }
      
    //cout << "--------------------------------------- \n";
    //cout << " Covariance matrix: ";
    //covMatrixStat.Print();

    TMatrixD covMatrixStatInv(nparConst,nparConst);
    if(nparConst>1) {
        covMatrixStatInv = covMatrixStat.InvertFast();
        for(unsigned int i=0; i<nparConst; ++i){
           for(unsigned int j=0; j<nparConst; ++j){
               covParConst[i][j] = covMatrixStatInv(i,j);
            }
          }
        }
        else {
          for(unsigned int i=0; i<nparConst; ++i){
             for(unsigned int j=0; j<nparConst; ++j){
                 covParConst[i][j] = 1/covMatrixStat(i,j);
            }
          }
     }
    

     return true;
    
}


bool fitter::SetExtInputs(){
  
  //if LCSR is present for both Bs -> Ds and Bs -> Ds*, read it once with full correlation matrix
  if(std::find(theoryInputs.begin(),theoryInputs.end(), "LCSRDs") != theoryInputs.end() &&
	 std::find(theoryInputs.begin(),theoryInputs.end(), "LCSRDsS") != theoryInputs.end()) {
	theoryInputs.erase(std::remove(theoryInputs.begin(), theoryInputs.end(), "LCSRDs"), theoryInputs.end());
	theoryInputs.erase(std::remove(theoryInputs.begin(), theoryInputs.end(), "LCSRDsS"), theoryInputs.end());
	theoryInputs.push_back("LCSR");
  }
  //read external inputs and store them
  for(auto model : theoryInputs) {
	if(model != "LHCb-PAPER-2019-046") {
	  TString n;
	  double x, y, err;
	  int i = 0;
	  ifstream ftmp("inputs/"+model+"_values.txt");
	  if(!ftmp.is_open()){
		cout << model << " input file not found!" << endl; return false;
	  }
	  while(ftmp >> n >> x >> y >> err) {
		_extInputs[model].name.push_back(n);
		_extInputs[model].x.push_back(x);
		_extInputs[model].x_err.push_back(0);
		_extInputs[model].val.push_back(y);
		_extInputs[model].val_err.push_back(err);
		i++;
	  }
	  _extInputs[model].dim = i;
	  _extInputs[model].corr_matrix.ResizeTo(i,i);
	  _extInputs[model].cov_matrix.ResizeTo(i,i);
	  
	  
	  //read correlation
	  double c;
	  ifstream ftmpc("inputs/"+model+"_corr.txt");
	  while(ftmpc >> n >> c) {
		_extInputs[model].corr[n] = c;
	  }
	  
	  //build covariance and its inverse
	  for(int i = 0; i < _extInputs[model].dim; i++) {
		for(int j = 0; j < _extInputs[model].dim; j++) {
		  double rho = _extInputs[model].corr[Form("%s_%s",_extInputs[model].name[i].Data(),_extInputs[model].name[j].Data())];
		  if(rho == 0) rho = _extInputs[model].corr[Form("%s_%s",_extInputs[model].name[j].Data(),_extInputs[model].name[i].Data())];
		  if(i == j) rho = 1;
		  _extInputs[model].cov_matrix(i,j) = rho*_extInputs[model].val_err[i]*_extInputs[model].val_err[j];
		  if(i!=j)
			_extInputs[model].cov_matrix(j,i) = _extInputs[model].cov_matrix(i,j);
		}
	  }
	  ftmp.close();
	  ftmpc.close();
	}
	else { //for the LHCb input, the covariance matrix is directly provided
	  //and we need to take into account bin widths (xerr)
	  TString n;
	  double x, y, xerr, err;
	  int i = 0;
	  ifstream ftmp("inputs/"+model+"_values.txt");
	  if(!ftmp.is_open()){
		cout << model << " input file not found!" << endl; return false;
	  }
	  while(ftmp >> n >> x >> xerr >> y >> err) {
		_extInputs[model].name.push_back(n);
		_extInputs[model].x.push_back(x);
		_extInputs[model].x_err.push_back(xerr);
		_extInputs[model].val.push_back(y);
		_extInputs[model].val_err.push_back(err);
		i++;
	  }
	  _extInputs[model].dim = i;
	  _extInputs[model].corr_matrix.ResizeTo(i,i);
	  _extInputs[model].cov_matrix.ResizeTo(i,i);
	  
	  //read covariance
	  double c;
	  ifstream ftmpc("inputs/"+model+"_cov.txt");
	  while(ftmpc >> n >> c) {
		_extInputs[model].cov[n] = c;
	  }
	  
	  //build covariance and its inverse
	  for(int i = 0; i < _extInputs[model].dim; i++) {
		for(int j = 0; j < _extInputs[model].dim; j++) {
		  double cov = _extInputs[model].cov[Form("%s_%s",_extInputs[model].name[i].Data(),_extInputs[model].name[j].Data())];
		  if(cov == 0) cov = _extInputs[model].cov[Form("%s_%s",_extInputs[model].name[j].Data(),_extInputs[model].name[i].Data())];
		  _extInputs[model].cov_matrix(i,j) = cov;
		  if(i!=j)
			_extInputs[model].cov_matrix(j,i) = _extInputs[model].cov_matrix(i,j);
		}
	  }
	  ftmp.close();
	  ftmpc.close();
	}
	
	///inverse covariance matrix
	TMatrixD cov_inv_matrix(_extInputs[model].dim,_extInputs[model].dim);
	cov_inv_matrix = _extInputs[model].cov_matrix.InvertFast();
	//cov_inv_matrix.Print();
	for(int i = 0; i < _extInputs[model].dim; i++)
	  for(int j = 0; j < _extInputs[model].dim; j++)
		_extInputs[model].cov_inv_matrix[i][j] = cov_inv_matrix(i,j);
	
  }
  
  return true;
}


void fitter::PrintConfigInfo(){
 
    FFModel FFModelFitDs  =  decFitDs->GetFFModel();
    FFModel FFModelFitDsS =  decFitDsS->GetFFModel();
    FFModel FFModelRefDs  =  decRefDs->GetFFModel();
    FFModel FFModelRefDsS =  decRefDsS->GetFFModel();
    
    cout << " REFERENCE FF models of the templates: \n";
    cout << " - Bs->Ds :  " << FFModelRefDs.model << ", parameters: " << FFModelRefDs.FFpars.size() << endl;
    for(auto p : FFModelRefDs.FFpars) cout << "\t " << p.name << "\t " << p.value << "\n";
    cout << " - Bs->Ds*:  " << FFModelRefDsS.model << ", parameters: " << FFModelRefDsS.FFpars.size() << endl;
    for(auto p : FFModelRefDsS.FFpars) cout << "\t " << p.name << "\t " << p.value << "\n";
    cout << " Corresponding to BR(Bs->Ds): " <<  decRefDs->Eval_BR() << ", BR(Bs->Ds*):" << decRefDsS->Eval_BR() << endl;
    cout << "--------------------------------------------" << endl;
    cout << " FIT FF models for the templates: \n";
    cout << " - Bs->Ds :  " << FFModelFitDs.model << ", parameters: " << FFModelFitDs.FFpars.size() << endl;
    cout << " - Bs->Ds*:  " << FFModelFitDsS.model << ", parameters: " << FFModelFitDsS.FFpars.size() << endl;
    cout << " Initial pars give BR(Bs->Ds): " <<  decFitDs->Eval_BR() << ", BR(Bs->Ds*):" << decFitDsS->Eval_BR() << endl;
    cout << "--------------------------------------------" << endl;
    cout << " Fit parameters: " << fitPars.size() << "\n";
    cout << "\tnum\tname\tvalue\terror\tmin\tmax\tgconst \n";
    for(auto p : fitPars) p.print();
    cout << "--------------------------------------------\n";
    cout << " Gaussian constrained parameters: ";
    cout << parConst.size() << endl;
    if(parConst.size() == 0) return;
    for(auto p : parConst) p.print();
    cout << "--------------------------------------------\n";
    cout << " External inputs: " << theoryInputs.size() << endl;
    for(auto model : theoryInputs) cout << "\t " << model;
    cout << "\n They are: " << endl;
    for(auto model : theoryInputs) { cout << "\t " << model << endl; _extInputs[model].print(); }
    cout << "--------------------------------------------\n" << endl;
   
	return;
}

bool fitter::FillMCcandidates(TString filename){
 
  cout << "---------------------------------------------------------" << endl;
  cout << " Fill signal templates " << endl;
    
  string line;
  ifstream fin(filename);
  TString infile, treename;
  while(fin>>line) {
    if(line=="Templates")  {
      fin >> infile >> treename;
    }
    else continue;
  }
  fin.close();
    
  TFile* file = TFile::Open(infile);
  cout << " File to get acceptances and candidates for templates: ";
  if(!file) { cout << "NOT FOUND!!!" << endl; return false; }
  else cout << infile << endl;
  
  hacc[0] = (TH1D*) file->Get("h_sig_Ds");
  hacc[1] = (TH1D*) file->Get("h_sig_DsS");
  hacc[0]->SetDirectory(0);
  hacc[1]->SetDirectory(0);
  
  double min{-999}, max{-999};
  for(int i=0; i<2; ++i){
	cout << " Get acceptance from histo: ";
	if(!hacc[i]){ cout << "NOT FOUND!!!" << endl; return false; }
	else cout << hacc[0]->GetName() << ", which has min and max values:" << endl;
	
	hacc[i]->GetMinimumAndMaximum(min,max);
	cout << " component: " << i << "\t " << min << "\t " << max << endl;
	
	fitComp[i].pdf = (TH1D*) hacc[i]->Clone(Form("h_fitComp_%i",i));
	fitComp[i].pdf->SetDirectory(0);
	fitComp[i].pdf->Reset();
  }
    
  TTree* tree = (TTree*) file->Get(treename);
  cout << " Take tree: ";
  if(!tree) { cout << "NOT FOUND!!!" << endl; return false; }
  else cout << treename << ", with entries " << tree->GetEntries() << endl;
 
  int signalCategory{-999};
  double w{-999}, ctl{-999}, ctd{-999}, chi{-999}, pperp_t{-999}, pperp{-999}, den{-999}, den1D=-999, den4D=-999;
  tree->SetBranchAddress("signalCategory",&signalCategory);
  tree->SetBranchAddress("w",&w);
  tree->SetBranchAddress("ctl",&ctl);
  tree->SetBranchAddress("ctd",&ctd);
  tree->SetBranchAddress("chi",&chi);
  tree->SetBranchAddress("pperp",&pperp);
  tree->SetBranchAddress("den",&den);
  tree->SetBranchAddress("den1D",&den1D);
  tree->SetBranchAddress("den4D",&den4D);
    
  for(int i=0; i<tree->GetEntries(); ++i){
    tree->GetEntry(i);
	
	MCcand mc;
	
	mc.signalCategory = signalCategory;
	mc.pperp = pperp;
	mc.wvar  = w;
	mc.ctl   = ctl;
	mc.ctd   = ctd;
	mc.chi   = chi;
	
	mc.den   = den;
	if(signalCategory==1){
	  if(_rate1D) mc.den= den1D;
	  else mc.den=den4D;
	}
	
	mccand.push_back(mc);
	
  }
      
  calculateYields();
  cout << " Yields of P and V components: " << fitComp[0].yield << " " << fitComp[1].yield << endl;
	
  file->Close();
    
  return true;
  
    
}

bool fitter::DoFit(double strategy, bool useHesse, bool useMinos) {
  cout << "==========================================================" << endl; 
  cout << " Fitting " << endl;
  cout << "----------------------------------------------------------" << endl;

  TMinuit *fitter = new TMinuit(fitPars.size());
  fitter->SetFCN(fcn_tot); 
  int ierflg(0);
  
  for(auto p : fitPars) fitter->mnparm(p.inum,p.name,p.value,p.error,p.min,p.max,ierflg);
  
  double arglist[2];  
  arglist[0] = strategy;   
  fitter->mnexcm("SET STRATEGY", arglist, 1, ierflg);

  double fmin(0.), fedm(0.), errdef(1.);  
  int npari(0), nparx(0), istat(0), n(0); 
  arglist[0] = 5000.;  
  arglist[1] = 0.001;
  
  while (istat!=3 && n++<3) {
    fitter->mnexcm("MIGRAD", arglist ,2, ierflg);
    fitter->mnstat(fmin, fedm, errdef, npari, nparx, istat);  
  }
  cout << "MIGRAD status = " << istat << endl;  

  if (istat!=3){
    cout << "Fit failed after 3 attempts" << endl;
    return false;
  }

  if(useHesse){          
    fitter->mnexcm("HESSE", arglist ,2, ierflg);     
    fitter->mnstat(fmin, fedm, errdef, npari, nparx, istat);    
    cout << "HESSE status = " << istat << endl;
  }                      
  if (istat!=3){           
    cout << "Hesse failed" << endl;     
    return false;
  }                                                           
  
  if(useMinos){                                     
    fitter->mnexcm("MINOS", arglist ,2, ierflg);          
    fitter->mnstat(fmin, fedm, errdef, npari, nparx, istat);           
    cout << "MINOS status = " << istat << endl;        
  }                                 

  if (istat!=3){  
    cout << "Minos failed" << endl;       
    return false;
  }                                        
  _ndf -= fitter->GetNumFreePars();           
  cout << "chi2/ndf = " << _chi2 << "/" << _ndf;           
  cout << ", prob = " << TMath::Prob(_chi2,_ndf) << endl; 

  cout << "----------------------------------------------------------" << endl;
  cout << " Fit results give BR(Bs->Ds): " <<  decFitDs->Eval_BR() << ", BR(Bs->Ds*):" << decFitDsS->Eval_BR() << endl;

  return true;

}

void fcn_tot(int &, double *, double &f, double *p, int ) {

  SetAllPars(p);
    
  FillHistogram();            

  _chi2 = 0.; _ndf = 0;
  for(int x=1; x<=hData->GetNbinsX(); ++x) {
        
    double data  = hData->GetBinContent(x);
    if(data==0) continue;
    double sigma = hData->GetBinError(x);       
    double sd2 = sigma*sigma;
    
    int ncomp = sizeof(fitComp)/sizeof(fitComp[0]);

    double mc(0.), smc2(0.);
    for(int comp = 0; comp < ncomp; comp++) {
      mc += fitComp[comp].pdf->GetBinContent(x);
      double sj = fitComp[comp].pdf->GetBinError(x);
        smc2 += sj*sj;
    }
    if(mc==0) continue;// da ricontrollare
      
    _chi2 += (data-mc)*(data-mc)/(sd2 + smc2);
    _ndf++;
  }
    
   //consider now gauss-constrained parameters
   if(parConst.size()!=0){

     vector<double> pg;
     for(auto p : fitPars) {
       if(!p.gconst) continue;
       pg.push_back(p.value);
     }
     if(pg.size()!=parConst.size()){ cout << "wrong setting of gauss constraints in the chi2 func!" << endl; }
     
     for(unsigned int i=0; i<parConst.size(); ++i) {
       for(unsigned int j=0; j<parConst.size(); ++j) {
           _chi2 += (pg[i] - parConst[i].value)*covParConst[i][j]*(pg[j] - parConst[j].value);
        }
      }
     _ndf+=parConst.size();
   }

   //add external inputs to chi2
  for(auto model : theoryInputs) {
	extInputs inputs = _extInputs[model];
	for(int i=0; i<inputs.dim; i++){
	  for(int j=0; j<inputs.dim; j++){
		_chi2 += (inputs.val[i] - FFfunctions(inputs.name[i],inputs.x[i],inputs.x_err[i]))*
			  (inputs.val[j] - FFfunctions(inputs.name[j],inputs.x[j],inputs.x_err[j]))*inputs.cov_inv_matrix[i][j];
		
		//cout << inputs.val[i] << "\t" << inputs.val[j] << "\t" << FFfunctions(inputs.name[i],inputs.x[i],inputs.x_err[i]) << "\t" <<
		//FFfunctions(inputs.name[j],inputs.x[j],inputs.x_err[j]) << "\t" << inputs.cov_inv_matrix[i][j] << endl;
	  }
	}
	_ndf+=inputs.dim;
  }
 
  f = _chi2;
  
}

double FFfunctions(TString xname, double w, double werr) {
  
  double FFvalue = 0;
  
  if(xname.Contains("f")){
	double Gw = 0;
	double r =  _Dsmass/_Bmass;
	FFModel ffmodel =  decFitDs->GetFFModel();
	if     (ffmodel.model == "CLN") decayRates::FFfunctionsCLN(w, ffmodel.FFpars, Gw);
	else if(ffmodel.model == "BGL") decayRates::FFfunctionsBGL(w, ffmodel.FFpars, Gw);
	else{ cout << "Wrong FF model for decFitDs in FFfunctions!"; return 1e20; }
	FFvalue = (1+r)/2./sqrt(r)*Gw;
  }
  else {
	double mB = _Bmass; double mDsS = _DsSmass;
	if(xname.Contains("bin")){
	  double norm =1.;
	  for(auto p : otherPars) if (p.name == "normLHCb") norm = p.value;
	  double wmin = 1.;
	  double wmax = (mB*mB + mDsS*mDsS) / (2.*mB*mDsS);
	  TF1 * f = new TF1("f",decFitDsS,&decayRates::TF_dGdw, wmin, wmax, 0);
	  FFvalue = norm*f->Integral(w-werr,w+werr);
	  f->Delete();
	}
	else{
	  double hw{0}, R1w{0}, R2w{0};
	  FFModel ffmodel = decFitDsS->GetFFModel();
	  if     (ffmodel.model == "CLN") decayRates::FFfunctionsCLN(w, ffmodel.FFpars, hw, R1w, R2w);
	  else if(ffmodel.model == "BGL"){
		double wp = 0;
		//if(xname.Contains("a2") && w==1) wp = 0.001;
		decayRates::FFfunctionsBGL(w+wp, ffmodel.FFpars, hw, R1w, R2w);
	  }
	  else { cout << "Wrong FF model for decFitDs in FFfunctions!"; return 1e20; }
	
	  /*
	  double fw = hw * sqrt(mB*mDsS)*(w+1.);
	  double gw = R1w*fw/(w+1)/mB/mDsS;
	  double q2 = mB*mB + mDsS*mDsS - 2*w*mB*mDsS;
	  double lambda = pow(mB,4) + pow(mDsS,4) + pow(q2,2) - 2.0*(mB*mB*mDsS*mDsS + mB*mB*q2 + mDsS*mDsS*q2);
	  double r = mDsS/mB;
	  double F1w = ( (w-r)/(w-1)-R2w ) * mB*(w-1)*fw;
	
	  double V = (mB+mDsS)/2.0*gw;
	  double A1 = 1/(mB+mDsS)*fw;
	  double A2 = (mB+mDsS)/lambda*( (mB*mB-mDsS*mDsS-q2)*fw - 2.0*mDsS*F1w );
	   */
	  
	  double r = mDsS/mB;
	  double R = 2.*sqrt(r)/(1.+r);
	  double A1 = hw/2.*(w+1.)*R;
	  double V  = R1w*hw/R;
	  double A2 = R2w*hw/R;
	
	  if     (xname.Contains("v")) { FFvalue = V;  if(FFvalue!=FFvalue) { cout << " w = " << w << ", V  = " <<  V  << endl;} }
	  else if(xname.Contains("a1")){ FFvalue = A1; if(FFvalue!=FFvalue) { cout << " w = " << w << ", A1 = " <<  A1 << endl;} }
	  else if(xname.Contains("a2")){ FFvalue = A2; if(FFvalue!=FFvalue) { cout << " w = " << w << ", A2 = " <<  A2 << endl;} }
	  else{
		cout << "Wrong name for x point for the external input in FFfunction!" << endl;
		return 1e20;
	  }
	}
  }
  
  
  return FFvalue;
}

void FillHistogram() { 

  fitComp[0].pdf->Reset();
  fitComp[1].pdf->Reset();

  double rateDs  = decFitDs->Eval_Gamma();
  double rateDsS = decFitDsS->Eval_Gamma();

  double num(0);

  for(auto cand : mccand) {
      if(cand.signalCategory == 0)   num = decFitDs->Eval_dGdw(cand.wvar) / rateDs;
	  else if(cand.signalCategory == 1) {
		if(_rate1D) num = decFitDsS->Eval_dGdw(cand.wvar) / rateDsS;
		else num = decFitDsS->Eval_dGdwdAngle(cand.wvar, cand.ctl, cand.ctd, cand.chi) / rateDsS;
	  }
   
    double wFF = num/cand.den;

    if (wFF!=wFF) {
        cout << "wFF is NaN, this should never happen!" << endl;
        printf("num: %f den: %f rateDs: %f rateDsS: %f \n", num, cand.den, rateDs, rateDsS);
		FFModel FFModelFitDs  =  decFitDs->GetFFModel();
		FFModel FFModelFitDsS =  decFitDsS->GetFFModel();
		cout << " - Bs->Ds :  " << FFModelFitDs.model << ", parameters: " << FFModelFitDs.FFpars.size() << endl;
	    for(auto p : FFModelFitDs.FFpars) p.print();
	    cout << " - Bs->Ds*:  " << FFModelFitDsS.model << ", parameters: " << FFModelFitDsS.FFpars.size() << endl;
	    for(auto p : FFModelFitDsS.FFpars) p.print();
	
        continue;
    }
    
    fitComp[cand.signalCategory].pdf->Fill(cand.pperp,wFF);
  }
  
  for(int i =0; i<2; ++i) fitComp[i].pdf->Multiply(hacc[i]);
    
  calculateYields();
  
  return;
}


void SetAllPars(double* p) {
    
    FFModel FFModelFitDs  =  decFitDs->GetFFModel();
    FFModel FFModelFitDsS =  decFitDsS->GetFFModel();
    
  for(unsigned int i = 0; i < fitPars.size(); i++){
     
     fitPars[i].value = p[i];
    
     if(i<FFModelFitDs.FFpars.size())  { FFModelFitDs.FFpars[i].value = p[i]; }
    
     else if(i>=FFModelFitDs.FFpars.size() && i<FFModelFitDs.FFpars.size()+FFModelFitDsS.FFpars.size()) { FFModelFitDsS.FFpars[i-FFModelFitDs.FFpars.size()].value = p[i]; }
    
     else if(i>=FFModelFitDs.FFpars.size()+FFModelFitDsS.FFpars.size()){ otherPars[i-FFModelFitDs.FFpars.size()-FFModelFitDsS.FFpars.size()].value = p[i]; }
  }
   
    decFitDs->SetFFModel(FFModelFitDs);
    decFitDsS->SetFFModel(FFModelFitDsS);
    
  return;

}

void SetAllPars() {
    
    FFModel FFModelFitDs  =  decFitDs->GetFFModel();
    FFModel FFModelFitDsS =  decFitDsS->GetFFModel();
    
    for(unsigned int i=0; i<FFModelFitDs.FFpars.size()+FFModelFitDsS.FFpars.size()+otherPars.size(); ++i){
        if(i<FFModelFitDs.FFpars.size())  {
            fitPars.push_back(FFModelFitDs.FFpars[i]);
            fitPars[i].inum = i;
        }
        else if(i>=FFModelFitDs.FFpars.size() && i<FFModelFitDs.FFpars.size()+FFModelFitDsS.FFpars.size()) {
            fitPars.push_back(FFModelFitDsS.FFpars[i-FFModelFitDs.FFpars.size()]);
            fitPars[i].inum = i; }
        else if(i>=FFModelFitDs.FFpars.size()+FFModelFitDsS.FFpars.size()){
            fitPars.push_back(otherPars[i-FFModelFitDs.FFpars.size()-FFModelFitDsS.FFpars.size()]);
            fitPars[i].inum = i; }
    }
}




void calculateYields() {

    double Vcb{0}, etaEW{0}, tauBs{0}, fsfdBrDs{0}, BrD{0}, BrDst{0},
        BrBdD{0}, BrBdDst{0}, NrefD{0}, NrefDst{0}, effRD{0}, effRDst{0},
        physBkg{0}, combBkg{0}, normLHCb{0};
    
    for(auto p : otherPars){
        if      (p.name == "Vcb")     Vcb = p.value;
        else if (p.name == "etaEW")   etaEW = p.value;
        else if (p.name == "tauBs")   tauBs = p.value;
        else if (p.name == "fsfdBrDs")fsfdBrDs = p.value;
        else if (p.name == "BrD")     BrD = p.value;
        else if (p.name == "BrDst")   BrDst = p.value;
        else if (p.name == "BrBdD")   BrBdD = p.value;
        else if (p.name == "BrBdDst") BrBdDst = p.value;
        else if (p.name == "NrefD")   NrefD = p.value;
        else if (p.name == "NrefDst") NrefDst = p.value;
        else if (p.name == "effRD")   effRD = p.value;
        else if (p.name == "effRDst") effRDst = p.value;
        else if (p.name == "physBkg") physBkg = p.value;
        else if (p.name == "combBkg") combBkg = p.value;
        else if (p.name == "normLHCb")normLHCb = p.value;
        else { cout << "ERROR: calculateYields: unkown parameter name!!! " << endl; return; }
    }
    
    double VcbEff = Vcb*etaEW;
    decFitDs ->SetTauB(tauBs);
    decFitDs ->SetVcbEff(VcbEff);
    decFitDsS->SetTauB(tauBs);
    decFitDsS->SetVcbEff(VcbEff);
    
    double BrBsDs  = decFitDs->Eval_BR();
    double BrBsDsS = decFitDsS->Eval_BR();
    
    
    fitComp[0].yield = effRD  * BrBsDs/ BrBdD  *  fsfdBrDs / BrD * NrefD ;
    fitComp[1].yield = effRDst* BrBsDsS/BrBdDst * fsfdBrDs / BrD / BrDst * NrefDst ;
    fitComp[2].yield = physBkg;
    fitComp[3].yield = combBkg;

    for(auto c : fitComp) c.pdf->Scale(c.yield/c.pdf->Integral());
    
    return;
}


void fitter::DrawResiduals(TH1D* hD, TH1D* hF, TH1D* hp){
  
  for(int i = 1; i <= hD->GetNbinsX() ; i++) {
	double pull = (hD->GetBinContent(i)-hF->GetBinContent(i))/
	sqrt(hD->GetBinError(i)*hD->GetBinError(i) +
		 hF->GetBinError(i)*hF->GetBinError(i));
	hp->SetBinContent(i,pull);
  }
  
  hp->SetLineColor(kBlack);
  hp->SetMarkerColor(kBlack);
  hp->GetYaxis()->SetNdivisions(505);
  hp->GetYaxis()->SetLabelSize(0.13);
  hp->GetYaxis()->SetRangeUser(-5,5);
  hp->GetXaxis()->SetLabelSize(0);
  hp->GetYaxis()->SetTitle("Pulls");
  hp->GetYaxis()->SetTitleSize(0.15);
  hp->GetYaxis()->SetTitleOffset(0.4);
  hp->GetXaxis()->SetTitleSize(0);
  
  TLine * l1 = new TLine(hD->GetXaxis()->GetXmin(),-3,hD->GetXaxis()->GetXmax(),-3);
  TLine * l2 = new TLine(hD->GetXaxis()->GetXmin(),+3,hD->GetXaxis()->GetXmax(),+3);
  TLine * l3 = new TLine(hD->GetXaxis()->GetXmin(),0,hD->GetXaxis()->GetXmax(),0);
  l1->SetLineStyle(kDashed);
  l2->SetLineStyle(kDashed);
  l3->SetLineStyle(kDashed);
  
  hp->GetYaxis()->SetRangeUser(-5,5);
  hp->SetFillColorAlpha(40, 0.35);
  hp->Draw("F");
  l1->Draw("SAME");
  l2->Draw("SAME");
  l3->Draw("SAME");
  
}


bool fitter::DoProjection() {
  
	plotStyle();
	FFModel FFModelFitDs  =  decFitDs->GetFFModel();
	FFModel FFModelFitDsS =  decFitDsS->GetFFModel();
  
	TCanvas* cperp = new TCanvas("cperp", "cperp", 800,1200);
	TPad* upperPad; TPad* lowerPad;
	upperPad = new TPad("plot_perp", "plot_perp", .005, .2525, .995, .995);
	lowerPad = new TPad("residuals_perp", "residuals_perp", .005, .005, .995, .2475);
	upperPad->SetBottomMargin(0.18);
	upperPad->SetLeftMargin(0.2);
	lowerPad->SetLeftMargin(0.2);
	upperPad->SetRightMargin(0.1);
	lowerPad->SetRightMargin(0.1);
	upperPad->SetTopMargin(0.1);
	lowerPad->Draw();
	upperPad->Draw();
	upperPad->cd();
	hData->SetMarkerStyle(20);
	hData->SetMarkerColor(kBlack);
	hData->SetLineColor(kBlack);
	hData->SetMarkerSize(0.8);
    hData->GetYaxis()->SetRangeUser(0,26.5e3);
	hData->GetYaxis()->SetNdivisions(505);
    hData->GetYaxis()->SetMaxDigits(3);
	hData->Draw("PE");
    TH1D * hFit = (TH1D*) hData->Clone("hFit");
    hFit->Reset();
    int colors[] = {46,38,14,12};
    double shades[]= {0.8,0.4,1,1};
    int ic=0;
	for(auto c : fitComp) {
	  c.pdf->SetLineColor(colors[ic]);
	  c.pdf->SetFillColorAlpha(colors[ic],shades[ic]);
	  c.pdf->SetFillStyle(3001);
	  c.pdf->SetLineWidth(2);
	  ++ic;
	  hFit->Add(c.pdf);
	}
    hFit->SetLineWidth(2);
	hFit->SetLineColor(9);
	fitComp[1].pdf->Draw("HISTSAME][");
	fitComp[0].pdf->Draw("HISTSAME][");
	fitComp[2].pdf->Draw("HISTSAME][");
	fitComp[3].pdf->Draw("HISTSAME][");
	hFit->Draw("HISTSAME][");
    hData->Draw("PESAME");
  
	TLegend * leg;
	leg = new TLegend(0.250,0.75,0.35,0.9);
	leg->SetFillStyle(0);
	leg->SetTextSize(0.040);
	leg->SetBorderSize(0);
	leg->SetTextFont(132);
	leg->AddEntry(hData,"LHCb, PRD101 (2020) 072004","lpe");
	leg->AddEntry(hFit,"Fit","l");
	leg->Draw("SAME");
    
	lowerPad->cd();
	TH1D * hpull = (TH1D*)hData->Clone("hpull");
	hpull->Reset();
    DrawResiduals(hData,hFit,hpull);
    cperp->SaveAs("fit_projection_pperp_Ds_"+FFModelFitDs.model+"_DsS_"+FFModelFitDsS.model+".pdf");
	cperp->SaveAs("fit_projection_pperp_Ds_"+FFModelFitDs.model+"_DsS_"+FFModelFitDsS.model+".C");
  
  
  TString fmodels[]  ={"HPQCD_Ds","MILC","LCSRDs"};
  const int Nfmodels = sizeof(fmodels)/sizeof(fmodels[0]);
  int fcolor[]={1,2,36};
  int marker[]={21,22,23};
  TGraphErrors* grfplus[Nfmodels];
  for(int im=0; im<Nfmodels; ++im){
	grfplus[im] = new TGraphErrors();
	grfplus[im]->SetMarkerColor(fcolor[im]);
	grfplus[im]->SetLineColor(fcolor[im]);
	grfplus[im]->SetMarkerStyle(marker[im]);
	grfplus[im]->SetMarkerSize(0.8);
	TString n;
	double x, y, err;
	int ip = 0;
	ifstream ftmp("inputs/"+fmodels[im]+"_values.txt");
	if(!ftmp.is_open()){
	  cout << fmodels[im] << " input file not found!" << endl; return false;
	}
	while(ftmp >> n >> x >> y >> err) {
	  grfplus[im]->SetPoint(ip,x,y);
	  grfplus[im]->SetPointError(ip,0,err);
	  ip++;
	}
  }
  TGraph* grfplusFit = new TGraph();
  grfplusFit->SetLineWidth(2);
  grfplusFit->SetLineColor(9);
  for(int p=0;p<1000;++p){
	double w = 1.0000001 + 1.4/1000*p;
	grfplusFit->SetPoint(p,w,FFfunctions("f",w,0));
  }
  TCanvas* cf = new TCanvas("cf","cf",800,800);
  cf->SetLeftMargin(0.2);
  cf->SetRightMargin(0.1);
  TMultiGraph* gf = new TMultiGraph();
  gf->Add(grfplusFit,"l");
  for(int im=0; im<Nfmodels; ++im) gf->Add(grfplus[im],"p");
  gf->GetHistogram()->GetXaxis()->SetRangeUser(0.9,2.5);
  gf->GetHistogram()->GetYaxis()->SetRangeUser(0,1.8);
  gf->GetYaxis()->SetTitle("#it{f}_{+}(#it{w})");
  gf->GetXaxis()->SetTitle("#it{w}");
  gf->Draw("a");
  TLegend * legf;
  legf = new TLegend(0.3,0.75,0.5,0.95);
  legf->SetFillStyle(0);
  legf->SetTextSize(0.040);
  legf->SetBorderSize(0);
  legf->SetTextFont(132);
  TString fmodelsName[]  ={"HPQCD, PRD101 (2020) 7, 074513","MILC, PRD85 (2012) 114502","LCSR, EPJC80 (2020) 4, 347"};
  for(int im=0; im<Nfmodels; ++im) legf->AddEntry(grfplus[im],fmodelsName[im],"pe");
  legf->AddEntry(grfplusFit,"Fit","l");
  legf->Draw("SAME");
  cf->SaveAs("fit_projection_fplus_Ds_"+FFModelFitDs.model+"_DsS_"+FFModelFitDsS.model+".pdf");
  cf->SaveAs("fit_projection_fplus_Ds_"+FFModelFitDs.model+"_DsS_"+FFModelFitDsS.model+".C");
  
  
  
  TString VAmodels[] ={"HPQCD_DsS","LCSRDsS"};
  const int NVAmodels = sizeof(VAmodels)/sizeof(VAmodels[0]);
  int colorsVA[]={1,36};
  int markersVA[]={21,23};
  TGraphErrors* grDst[3][NVAmodels];
  for(int im=0; im<NVAmodels; ++im){
	for(int i=0;i<3;++i) {
	  grDst[i][im] = new TGraphErrors();
	  grDst[i][im]->SetLineWidth(2);
	  grDst[i][im]->SetMarkerColor(colorsVA[im]);
	  grDst[i][im]->SetLineColor(colorsVA[im]);
	  grDst[i][im]->SetMarkerStyle(markersVA[im]);
	  grDst[i][im]->SetMarkerSize(0.8);
	}
	TString n;
	double x, y, err;
	int ipv{0},ipa1{0},ipa2{0};
	ifstream ftmp("inputs/"+VAmodels[im]+"_values.txt");
	if(!ftmp.is_open()){
	  cout << VAmodels[im] << " input file not found!" << endl; return false;
	}
	while(ftmp >> n >> x >> y >> err) {
	  //cout << n << x << y << err<< endl;
	  if     (n.Contains("v") ){ grDst[0][im]->SetPoint(ipv,x,y);  grDst[0][im]->SetPointError(ipv,0,err);  ipv++; }
	  else if(n.Contains("a1")){ grDst[1][im]->SetPoint(ipa1,x,y); grDst[1][im]->SetPointError(ipa1,0,err); ipa1++;}
	  else if(n.Contains("a2")){ grDst[2][im]->SetPoint(ipa2,x,y); grDst[2][im]->SetPointError(ipa2,0,err); ipa2++;}
	}
  }
  TGraph* grDstFit[3];
  TString titley[]={"#it{V}(#it{w})","#it{A}_{1}(#it{w})","#it{A}_{2}(#it{w})"};
  for(int i=0;i<3;++i){
	grDstFit[i] = new TGraph();
	grDstFit[i]->SetLineWidth(2);
	grDstFit[i]->SetLineColor(9);
  }
  for(int p=0;p<1000;++p){
	double w = 1.0000001 + 1.2/1000*p;
	grDstFit[0]->SetPoint(p,w,FFfunctions("v" ,w,0));
	grDstFit[1]->SetPoint(p,w,FFfunctions("a1",w,0));
	grDstFit[2]->SetPoint(p,w,FFfunctions("a2",w,0));
  }
  
  TGraphErrors* grDstHPQCD = new TGraphErrors();
  grDstHPQCD->SetLineWidth(2);
  grDstHPQCD->SetMarkerColor(4);
  grDstHPQCD->SetLineColor(4);
  grDstHPQCD->SetMarkerStyle(20);
  grDstHPQCD->SetMarkerSize(0.8);
  double mB = _Bmass; double mDsS = _DsSmass;
  double hpqcdDsS  = 1./(mB+mDsS)*0.902* sqrt(mB*mDsS)*(2.);
  double ehpqcdDsS = 1./(mB+mDsS)*0.013 * sqrt(mB*mDsS)*(2.);
  grDstHPQCD->SetPoint(0,1,hpqcdDsS);
  grDstHPQCD->SetPointError(0,0,ehpqcdDsS);
  
  TMultiGraph* gVA[3];
  TCanvas* cVA[3];
  TLegend * legVA[3];
  TString funcname[]={"V","A1","A2"};
  for(int cp=0;cp<3;++cp){
	cVA[cp] = new TCanvas(Form("c_%i",cp),Form("c_%i",cp),800,800);
	cVA[cp]->SetLeftMargin(0.2);
	cVA[cp]->SetRightMargin(0.1);
	gVA[cp] = new TMultiGraph();
	gVA[cp]->Add(grDstFit[cp],"l");
	for(int m=0; m<NVAmodels; m++) gVA[cp]->Add(grDst[cp][m],"p");
	if(cp==1) gVA[cp]->Add(grDstHPQCD,"p");
    gVA[cp]->GetHistogram()->GetXaxis()->SetRangeUser(0.9,2.5);
	gVA[cp]->GetHistogram()->GetYaxis()->SetRangeUser(0,1.5);
	if(cp==0) gVA[cp]->GetHistogram()->GetYaxis()->SetRangeUser(0,2.1);
	gVA[cp]->GetYaxis()->SetTitle(titley[cp]);
	gVA[cp]->GetXaxis()->SetTitle("#it{w}");
	gVA[cp]->Draw("a");
	legVA[cp] = new TLegend(0.35,0.75,0.55,0.95);
	legVA[cp]->SetFillStyle(0);
	legVA[cp]->SetTextSize(0.040);
	legVA[cp]->SetBorderSize(0);
	legVA[cp]->SetTextFont(132);
	if(cp==0) gVA[cp]->GetHistogram()->GetYaxis()->SetRangeUser(0,2.1);
	if(cp==1) legVA[cp]->AddEntry(grDstHPQCD,"HPQCD, PRD99 (2019) 114512","pe");
	legVA[cp]->AddEntry(grDst[cp][0],"HPQCD, arXiv:2105.11433","pe");
	legVA[cp]->AddEntry(grDst[cp][1],"LCSR, EPJC80 (2020) 4, 347","pe");
	legVA[cp]->AddEntry(grDstFit[cp],"Fit","l");
	legVA[cp]->Draw("SAME");
	
	cVA[cp]->SaveAs("fit_projection_"+funcname[cp]+"_Ds_"+FFModelFitDs.model+"_DsS_"+FFModelFitDsS.model+".pdf");
	cVA[cp]->SaveAs("fit_projection_"+funcname[cp]+"_Ds_"+FFModelFitDs.model+"_DsS_"+FFModelFitDsS.model+".C");
  }
  
  
  double wBins[] = {1.000, 1.1087, 1.1688, 1.2212, 1.2717, 1.3226, 1.3814, 1.4667};
  TH1D* hw = new TH1D("h_w","h_w",7,wBins);//wLHCb
  TH1D* hwFit = (TH1D*) hw->Clone("h_wFit");
  TString n;
  double x, y, xerr, yerr;
  int ip = 1;
  ifstream ftmp("inputs/LHCb-PAPER-2019-046_values.txt");
  if(!ftmp.is_open()){
	cout << "LHCb w data input file not found!" << endl; return false;
  }
  while(ftmp >> n >> x >> xerr>> y >> yerr) {
	hw->SetBinContent(ip,y/hw->GetBinWidth(ip));
	hw->SetBinError(ip,yerr/hw->GetBinWidth(ip));
	hwFit->SetBinContent(ip,FFfunctions("bin",x,xerr)/hwFit->GetBinWidth(ip));
	hwFit->SetBinError(ip,0);
	ip++;
  }
  
  TCanvas* cw = new TCanvas("cw", "cw", 800,1200);
  TPad* upperPadw; TPad* lowerPadw;
  upperPadw = new TPad("plot_w", "plot_w", .005, .2525, .995, .995);
  lowerPadw = new TPad("residuals_w", "residuals_w", .005, .005, .995, .2475);
  upperPadw->SetBottomMargin(0.18);
  upperPadw->SetLeftMargin(0.2);
  lowerPadw->SetLeftMargin(0.2);
  upperPadw->SetRightMargin(0.1);
  lowerPadw->SetRightMargin(0.1);
  upperPadw->SetTopMargin(0.1);
  lowerPadw->Draw();
  upperPadw->Draw();
  upperPadw->cd();
  hw->SetMarkerStyle(20);
  hw->SetMarkerColor(kBlack);
  hw->SetLineColor(kBlack);
  hw->SetMarkerSize(0.8);
  hw->GetYaxis()->SetRangeUser(0.5,3.5);
  hw->GetYaxis()->SetTitle("(1/#it{N})(d#it{N}/d#it{w})");
  hw->GetXaxis()->SetTitle("#it{w}");
  hw->Draw("PE");
  hwFit->SetLineWidth(2);
  hwFit->SetLineColor(9);
  hwFit->Draw("HISTSAME][");
  hw->Draw("PESAME");
  
  TLegend * legw;
  legw = new TLegend(0.25,0.75,0.50,0.90);
  legw->SetFillStyle(0);
  legw->SetTextSize(0.040);
  legw->SetBorderSize(0);
  legw->SetTextFont(132);
  legw->AddEntry(hw,"LHCb, JHEP12 (2020) 144","lpe");
  legw->AddEntry(hwFit,"Fit","l");
  legw->Draw("SAME");
  
  lowerPadw->cd();
  TH1D * hpullw = (TH1D*) hw->Clone("hpullw");
  hpullw->Reset();
  DrawResiduals(hw,hwFit,hpullw);
  
  cw->SaveAs("fit_projection_wLHCb_Ds_"+FFModelFitDs.model+"_DsS_"+FFModelFitDsS.model+".pdf");
  cw->SaveAs("fit_projection_wLHCb_Ds_"+FFModelFitDs.model+"_DsS_"+FFModelFitDsS.model+".C");
  
  return true;
	
}


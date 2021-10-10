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
	if (!DoFit(1,true,false)) return false;//arguments: strategy, useHesse, useMinos
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
      fin >> np >> dirInputs;
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
      
    //cout << "--------------------------------------------------- \n";
    //cout << " Covariance matrix gaussian constrained parameters: ";
    //covMatrixStat.Print();

    TMatrixD covMatrixStatInv(nparConst,nparConst);
    if(nparConst>1) {
	  if(covMatrixStat.Determinant()==0) { cout << "cannot invert Gauss const covariance, det = 0" << endl; return false;}
        //covMatrixStatInv = covMatrixStat.Invert();
		covMatrixStatInv = covMatrixStat;
	  TDecompLU lu(covMatrixStatInv);
	  lu.SetTol(1.0e-24);
	  int nr = 0;
	  while(!lu.Invert(covMatrixStatInv) && nr < 10) {
		lu.SetMatrix(covMatrixStatInv);
		lu.SetTol(0.1*lu.GetTol());
		nr++;
	  }
	  cout << "Gauss constr cov matrix inversion sucessfully done\n" << endl;
	  //covMatrixStatInv.Print();
	  //cout  << "check the inverse: A*A-1 = 1" << endl;
	  //TMatrixD identity = covMatrixStat*covMatrixStatInv;
	  //identity.Print();
	  
        for(unsigned int i=0; i<nparConst; ++i){
           for(unsigned int j=0; j<nparConst; ++j){
               covParConst[i][j] = covMatrixStatInv(i,j);
            }
          }
        }
        else {
          for(unsigned int i=0; i<nparConst; ++i){
             for(unsigned int j=0; j<nparConst; ++j){
                 covParConst[i][j] = 1./covMatrixStat(i,j);
            }
          }
     }
    //cout << "Matrix inversion done" << endl;

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
	cout << "Setting ext inputs for " << model << endl;
	if(model != "LHCb-PAPER-2019-046") {
	  TString n;
	  double x, y, err;
	  int i = 0;
	  ifstream ftmp(dirInputs+"/"+model+"_values.txt");
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
	  ifstream ftmpc(dirInputs+"/"+model+"_corr.txt");
	  while(ftmpc >> n >> c) {
		_extInputs[model].corr[n] = c;
	  }
	  
	  //build covariance and its inverse
	  for(int i = 0; i < _extInputs[model].dim; i++) {
		for(int j = 0; j <=i; j++) {
		  double rho = 0;
		  if(i == j) rho = 1;
		  else {
			if(_extInputs[model].corr.find(Form("%s_%s",_extInputs[model].name[i].Data(),_extInputs[model].name[j].Data()))!=_extInputs[model].corr.end())
			  rho = _extInputs[model].corr[Form("%s_%s",_extInputs[model].name[i].Data(),_extInputs[model].name[j].Data())];
			else if (_extInputs[model].corr.find(Form("%s_%s",_extInputs[model].name[j].Data(),_extInputs[model].name[i].Data()))!=_extInputs[model].corr.end())
			  rho = _extInputs[model].corr[Form("%s_%s",_extInputs[model].name[j].Data(),_extInputs[model].name[i].Data())];
			//if(_extInputs[model].name[j].Contains("fzero") || _extInputs[model].name[j].Contains("a0") ||
			//   _extInputs[model].name[i].Contains("fzero") || _extInputs[model].name[i].Contains("a0")) rho=0;
		  }
		  _extInputs[model].cov_matrix(i,j) = rho*_extInputs[model].val_err[i]*_extInputs[model].val_err[j];
		  _extInputs[model].cov_matrix(j,i) = _extInputs[model].cov_matrix(i,j);
		}
	  }
	  ftmp.close();
	  ftmpc.close();
	}
	else { //for the LHCb input, the covariance matrix is directly provided
	  //and we need to take into account bin widths (xerr)
	  _normHistLHCb=0;
	  TString n;
	  double x, y, xerr, err;
	  int i = 0;
	  ifstream ftmp(dirInputs+"/"+model+"_values.txt");
	  if(!ftmp.is_open()){
		cout << model << " input file not found!" << endl; return false;
	  }
	  while(ftmp >> n >> x >> xerr >> y >> err) {
		_extInputs[model].name.push_back(n);
		_extInputs[model].x.push_back(x);
		_extInputs[model].x_err.push_back(xerr);
		_extInputs[model].val.push_back(y);
		_extInputs[model].val_err.push_back(err);
		_normHistLHCb+=y;
		i++;
	  }
	  _extInputs[model].dim = i;
	  _extInputs[model].corr_matrix.ResizeTo(i,i);
	  _extInputs[model].cov_matrix.ResizeTo(i,i);
	  
	  //read covariance
	  double c;
	  ifstream ftmpc(dirInputs+"/"+model+"_cov.txt");
	  while(ftmpc >> n >> c) {
		_extInputs[model].cov[n] = c;
	  }
	  
	  //build covariance and its inverse
	  for(int i = 0; i < _extInputs[model].dim; i++) {
		for(int j = 0; j <=i ; j++) {
		  double cov = 0;
		  if(_extInputs[model].cov.find(Form("%s_%s",_extInputs[model].name[i].Data(),_extInputs[model].name[j].Data()))!=_extInputs[model].cov.end())
			cov = _extInputs[model].cov[Form("%s_%s",_extInputs[model].name[i].Data(),_extInputs[model].name[j].Data())];
		  else if(_extInputs[model].cov.find(Form("%s_%s",_extInputs[model].name[j].Data(),_extInputs[model].name[i].Data()))!=_extInputs[model].cov.end())
			cov = _extInputs[model].cov[Form("%s_%s",_extInputs[model].name[j].Data(),_extInputs[model].name[i].Data())];
		  
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
	TMatrixD tobeInverted = _extInputs[model].cov_matrix;
	cout << "Covariance input matrix inversion..." << endl;
	double det = tobeInverted.Determinant();
	cout << " matrix determinant = " << det << endl;
	if(det==0 || det!=det) { cout << "cannot invert covariance for model " << model << endl; return false;}
	TDecompLU lu(tobeInverted);
	lu.SetTol(1.0e-24);
	int nr = 0;
	while(!lu.Invert(tobeInverted) && nr < 10) {
	  lu.SetMatrix(tobeInverted);
	  lu.SetTol(0.1*lu.GetTol());
	  nr++;
	}
	cout << " matrix inversion done" << endl;
	//tobeInverted.Print();
	//cout  << "check the inverse: A*A-1 = 1" << endl;
	//TMatrixD identity = _extInputs[model].cov_matrix*tobeInverted;
	//identity.Print();
	
	for(int i = 0; i < _extInputs[model].dim; i++)
	  for(int j = 0; j < _extInputs[model].dim; j++)
		_extInputs[model].cov_inv_matrix[i][j] = tobeInverted(i,j);
	
	cout << "-------------------------------------\n" << endl;
	
  }
  
  return true;
}


void fitter::PrintConfigInfo(){
 
    FFModel FFModelFitDs  =  decFitDs->GetFFModel();
    FFModel FFModelFitDsS =  decFitDsS->GetFFModel();
    FFModel FFModelRefDs  =  decRefDs->GetFFModel();
    FFModel FFModelRefDsS =  decRefDsS->GetFFModel();
    
    cout << " REFERENCE FF models of the templates: \n";
    cout << " - Bs->Ds :  model " << FFModelRefDs.model << ", check is Ds*? " << decRefDs->GetIsDst() <<  ", VcbEff = " << decRefDs->GetVcbEff() << ", tauB = " << decRefDs->GetTauB()<< ", parameters: " << FFModelRefDs.FFpars.size() <<endl;
    for(auto p : FFModelRefDs.FFpars) cout << "\t " << p.name << "\t " << p.value << "\n";
    cout << " - Bs->Ds*:  model " << FFModelRefDsS.model << ", check is Ds*? " << decRefDsS->GetIsDst()<< ", VcbEff = " << decRefDsS->GetVcbEff() << ", tauB = " << decRefDsS->GetTauB()<< ", parameters: " << FFModelRefDsS.FFpars.size() << endl;
    for(auto p : FFModelRefDsS.FFpars) cout << "\t " << p.name << "\t " << p.value << "\n";
    cout << " Corresponding to BR(Bs->Ds): " <<  decRefDs->Eval_BR() << ", BR(Bs->Ds*):" << decRefDsS->Eval_BR() << endl;
    cout << "--------------------------------------------" << endl;
    cout << " FIT FF models for the templates: \n";
    cout << " - Bs->Ds :  model " <<  FFModelFitDs.model << ", check is Ds*? " << decFitDs->GetIsDst() <<  ", VcbEff = " << decFitDs->GetVcbEff() << ", tauB = " << decFitDs->GetTauB()<< ", parameters: " << FFModelFitDs.FFpars.size() << endl;
	for(auto p : FFModelFitDs.FFpars) cout << "\t " << p.name << "\t " << p.value << "\n";
    cout << " - Bs->Ds*:  model " <<  FFModelFitDsS.model << ", check is Ds*? " << decFitDsS->GetIsDst()<<  ", VcbEff = " << decFitDsS->GetVcbEff() << ", tauB = " << decFitDsS->GetTauB()<< ", parameters: " << FFModelFitDsS.FFpars.size() << endl;
    for(auto p : FFModelFitDsS.FFpars) cout << "\t " << p.name << "\t " << p.value << "\n";
    cout << " Initial pars give BR(Bs->Ds): " <<  decFitDs->Eval_BR() << ", BR(Bs->Ds*):" << decFitDsS->Eval_BR() << endl;
    cout << "--------------------------------------------" << endl;
    cout << " Fit parameters: " << fitPars.size() << "\n";
    cout << "\tnum\tname\tvalue\terror\tmin\tmax\tgconst \n";
    for(auto p : fitPars) p.print();
    cout << "--------------------------------------------\n";
    cout << " Gaussian constrained parameters: ";
    cout << parConst.size() << endl;
    if(parConst.size() != 0) for(auto p : parConst) p.print();
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
  int migradstat = istat;

  if (istat!=3){
    cout << "Fit failed after 3 attempts" << endl;
    return false;
  }

  int hessestat = -999;
  if(useHesse){          
    fitter->mnexcm("HESSE", arglist ,2, ierflg);     
    fitter->mnstat(fmin, fedm, errdef, npari, nparx, istat);    
    cout << "HESSE status = " << istat << endl;
	hessestat = istat;
	if (istat!=3)
	  cout << "!!! -- HESSE FAILED -- !!!" << endl;
  }                                                           
  
  int minosstat = -999;
  if(useMinos){                                     
    fitter->mnexcm("MINOS", arglist ,2, ierflg);          
    fitter->mnstat(fmin, fedm, errdef, npari, nparx, istat);
    cout << "MINOS status = " << istat << endl;
	minosstat = istat;
	if (istat!=3)
	  cout << "!!! -- MINOS FAILED -- !!!" << endl;
  }
  
  int npf = fitter->GetNumFreePars();
  _ndf -= npf;
  
  cout << "chi2/ndf = " << _chi2 << "/" << _ndf;
  cout << ", prob = " << TMath::Prob(_chi2,_ndf) << endl;
  cout << "----------------------------------------------------------" << endl;
  cout << " Fit results give BR(Bs->Ds): " <<  decFitDs->Eval_BR() << ", BR(Bs->Ds*):" << decFitDsS->Eval_BR() << endl;
  
  // expand covariance matrix to account for fixed parameters
  double fitMatrix[fitPars.size()][fitPars.size()];
  fitter->mnemat(&fitMatrix[0][0],fitPars.size());
    
  for (unsigned int i=0; i<fitPars.size(); ++i) {
    int ioff(0);                                                                          
    for(unsigned int k=0; k<i; ++k)                                                       
      if (fitPars[k].error==0.) ioff++;                                                       
    for (unsigned int j=0; j<fitPars.size(); ++j) {                                                
      if (fitPars[i].error==0. || fitPars[j].error==0.) {                                         
        fitFullMatrix[i][j] = 0.;
        continue;                                                                         
      }                                                                                   
      int joff(0);                                                                        
      for(unsigned int k=0; k<j; ++k)                                                     
        if (fitPars[k].error==0.) joff++;
		fitFullMatrix[i][j] = fitMatrix[i-ioff][j-joff];
    }
  }
  
  //print out results:
  FFModel FFModelFitDs  =  decFitDs->GetFFModel();
  FFModel FFModelFitDsS =  decFitDsS->GetFFModel();
  ofstream fresult("results_fit_Ds_"+FFModelFitDs.model+"_DsS_"+FFModelFitDsS.model+".txt");
  fresult << " Fit results " << endl;
  fresult << "----------------------------------------------------------" << endl;
  fresult << " MIGRAD status = " << migradstat << endl;
  fresult << " HESSE  status = "; if(useHesse) fresult << hessestat << endl; else fresult << "not used" << endl;
  fresult << " MINOS  status = "; if(useMinos) fresult << minosstat << endl; else fresult << "not used" << endl;
  fresult << "----------------------------------------------------------" << endl;
  for(auto p : fitPars)
	fresult << p.inum << "\t" << p.name << "\t\t" << p.value << "\t+-\t" << sqrt(fitFullMatrix[p.inum][p.inum]) << endl;
  fresult << "----------------------------------------------------------" << endl;
  fresult << " chi2/ndf = " << _chi2 << "/" << _ndf;
  fresult << ", prob = " << TMath::Prob(_chi2,_ndf) << endl;
  fresult << "----------------------------------------------------------" << endl;
  fresult << " Covariance matrix elements: " << endl;
  for(unsigned int i = 0; i<fitPars.size(); ++i){
	for(unsigned int j = 0; j<=i; ++j)
	  fresult << i << "\t " << j << "\t " << fitFullMatrix[i][j] << endl;
  }
  fresult.close();
  
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
	//cout << " ================ " << endl;
	//cout << model << endl;
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
 
  //cout << " chi2 = " << _chi2 << endl;
  f = _chi2;
  
}



//______________________________________________
void fitter::DrawFFErrorBand(TString xname, std::vector<TGraphErrors*>& gr) {

  TColor * col = new TColor(); 
  //Color_t c1 = col->GetColor(227,74,51);                                                                                  
  //Color_t c2 = col->GetColor(253,187,132);      
  //Color_t c3 = col->GetColor(254,232,200);      
  Color_t c1 = col->GetColor(49,163,84);     
  Color_t c2 = col->GetColor(161,217,155);      
  Color_t c3 = col->GetColor(229,245,224);      


  double FFvalue = 0;
  double w = 0;
  FFModel ffmodelDs =  decFitDs->GetFFModel();
  FFModel ffmodelDsS =  decFitDsS->GetFFModel();
  //quantiles of gaussian
  const int dim = 6;
  double quantiles[dim];
  double sum[dim] = {0.15865, 1-0.15865, 0.0227, 1-0.0227, 0.0015, 1-0.0015};
  //TGraphErrors
  gr.resize(3);
  for(int i = 0; i < 3; i++) gr[i] = new TGraphErrors();

  vector<parameter> allPars;
  for(auto p : fitPars) {
    if(p.error == 0) continue;
    allPars.push_back(p);
  }

  vector<parameter> newPars;
  vector<parameter> tmpFFPars;
  newPars = fitPars;
  //for (auto p : allPars) p.print();
  //get the fit covariance matrix only for free pars
  TMatrixD cov(allPars.size(),allPars.size());
  for (unsigned int i=0; i<fitPars.size(); ++i) {
	int ioff(0);
	for(unsigned int k=0; k<i; ++k)
	  if (fitPars[k].error==0.) ioff++;
	for (unsigned int j=0; j<fitPars.size(); ++j) {
	  if (fitPars[i].error==0. || fitPars[j].error==0.)continue;
	  int joff(0);
	  for(unsigned int k=0; k<j; ++k)
		if (fitPars[k].error==0.) joff++;
	   cov(i-ioff,j-joff) = fitFullMatrix[i][j];
	}
  }
  //cov.Print();
  

  if(xname.Contains("proj")) {

    std::vector<TH1D*> hg;
    for(int i = 0; i < 20; i++) {
      TH1D * hgtmp =  new TH1D(Form("hg%d",i),Form("hg%d",i),10000,0,50e3);
      hg.push_back(hgtmp);
    }
    for(int i = 0; i < 4; i++) {
      fitBand[i] = fitComp[i];
    }

    double Vcb(1), etaEW(1);
    double effRD(1), BrBsDs(1), BrBdD(1), fsfdBrDs(1), BrD(1), NrefD(1);
    double effRDst(1), BrBsDsS(1), BrBdDst(1), BrDst(1), NrefDst(1);
    double physBkg(1), combBkg(1);

    double rateDs  = decFitDs->Eval_Gamma();
    double rateDsS = decFitDsS->Eval_Gamma(); 

    //need new models to draw the bars, to avoid modifying the origianl one
    decayRates * newModelDs = new decayRates(decFitDs->GetFFModel(),false);
    FFModel newFFModelDs;
    newFFModelDs.model = ffmodelDs.model;
    newFFModelDs.initNFFpars(ffmodelDs.FFpars.size());	
    newModelDs->SetTauB(decFitDs->GetTauB());                                 
    newModelDs->SetVcbEff(decFitDs->GetVcbEff());    

    decayRates * newModelDsS = new decayRates(decFitDsS->GetFFModel(),true);
    FFModel newFFModelDsS;
    newFFModelDsS.model = ffmodelDsS.model;
    newFFModelDsS.initNFFpars(ffmodelDsS.FFpars.size());	
    newModelDsS->SetTauB(decFitDsS->GetTauB());                                 
    newModelDsS->SetVcbEff(decFitDsS->GetVcbEff());    

      for(int i = 0; i < 100; i++) { //nExtractions per point
	fitBand[0].pdf->Reset();
	fitBand[1].pdf->Reset();
	tmpFFPars.clear();
	TVectorD shift(allPars.size());
	TDecompChol tdc(cov);
	if ( !tdc.Decompose() ) {cout << "Something is wrong with decomposition!" << endl; return;}
	TMatrixD u = tdc.GetU(); u.T();
	TRandom3 ran(i+1);
	
	for (long unsigned int m=0; m<allPars.size(); ++m) shift(m) = ran.Gaus();
	shift = u*shift;
	
	for (long unsigned int m=0; m<allPars.size(); ++m) {
	  for (long unsigned int n=0; n<newPars.size(); ++n)
	    //if(newPars[n].inum == allPars[m].inum) newPars[n].value = allPars[m].value + shift(m);
	    if(newPars[n].inum == allPars[m].inum) newPars[n].value = allPars[m].value;
	}
	
	for (long unsigned int m=0; m<newPars.size(); ++m) {
	  if(newPars[m].name == "Vcb")           Vcb  = newPars[m].value;
	  else if(newPars[m].name == "etaEW")    etaEW  = newPars[m].value;
	  else if(newPars[m].name == "effRD")    effRD  = newPars[m].value;
	  else if(newPars[m].name == "BrBsDs")   BrBsDs  = newPars[m].value;
	  else if(newPars[m].name == "BrBdD")    BrBdD  = newPars[m].value;
	  else if(newPars[m].name == "fsfdBrDs") fsfdBrDs  = newPars[m].value;
	  else if(newPars[m].name == "BrD")      BrD  = newPars[m].value;
	  else if(newPars[m].name == "NrefD")    NrefD  = newPars[m].value;
	  else if(newPars[m].name == "effRDst")  effRDst  = newPars[m].value;
	  else if(newPars[m].name == "BrBsDsS")  BrBsDsS  = newPars[m].value;
	  else if(newPars[m].name == "BrBdDst")  BrBdDst  = newPars[m].value;
	  else if(newPars[m].name == "BrDst")    BrDst  = newPars[m].value;
	  else if(newPars[m].name == "NrefDst")  NrefDst  = newPars[m].value;
	  else if(newPars[m].name == "physBkg")  physBkg  = newPars[m].value;
	  else if(newPars[m].name == "combBkg")  combBkg  = newPars[m].value;
	}
	for (long unsigned int m=0; m<newPars.size(); ++m) {
	  if(m < ffmodelDs.FFpars.size())  tmpFFPars.push_back(newPars[m]);
	  if(m >= ffmodelDs.FFpars.size() && m < ffmodelDsS.FFpars.size()+ffmodelDs.FFpars.size()) tmpFFPars.push_back(newPars[m]);
	}

	for(long unsigned int m = 0; m < tmpFFPars.size(); ++m) {
	  if(m < ffmodelDs.FFpars.size()) newFFModelDs.FFpars[m] = tmpFFPars[m]; 
	  if(m >= ffmodelDs.FFpars.size() && m < ffmodelDsS.FFpars.size()+ffmodelDs.FFpars.size()) newFFModelDsS.FFpars[m-ffmodelDs.FFpars.size()] = tmpFFPars[m]; 
	}

	
	newModelDs->SetFFModel(newFFModelDs);
	newModelDs->SetVcbEff(Vcb*etaEW);
	newModelDsS->SetFFModel(newFFModelDsS);
	newModelDsS->SetVcbEff(Vcb*etaEW);

	double num(0);
	for(auto cand : mccand) {                                                                                                     
	  if(cand.signalCategory == 0)   num = newModelDs->Eval_dGdw(cand.wvar) / rateDs;                                             
	  else if(cand.signalCategory == 1) {                                                                                   
	    if(_rate1D) num = newModelDsS->Eval_dGdw(cand.wvar) / rateDsS;                                                    
	    else num = newModelDsS->Eval_dGdwdAngle(cand.wvar, cand.ctl, cand.ctd, cand.chi) / rateDsS;                       
	  }
	  
	  double wFF = num/cand.den;                                                                                                  
	  
	  if (wFF!=wFF) {                                                                                                             
	    cout << "wFF is NaN, this should never happen!" << endl;                                                                
	    printf("num: %f den: %f rateDs: %f rateDsS: %f \n", num, cand.den, rateDs, rateDsS);                                    
	    FFModel FFModelFitDs  =  newModelDs->GetFFModel();                                                                
	    FFModel FFModelFitDsS =  newModelDsS->GetFFModel();                                                               
	    cout << " - Bs->Ds :  " << FFModelFitDs.model << ", parameters: " << FFModelFitDs.FFpars.size() << endl;        
	    for(auto p : FFModelFitDs.FFpars) p.print();                                                                        
	    cout << " - Bs->Ds*:  " << FFModelFitDsS.model << ", parameters: " << FFModelFitDsS.FFpars.size() << endl;          
	    for(auto p : FFModelFitDsS.FFpars) p.print();                                                                       
	    
	    continue;                                                                                                               
	  }                                                                                                                           


	  fitBand[cand.signalCategory].pdf->Fill(cand.pperp,wFF);                                                                     
	}

	for(int i =0; i<2; ++i) fitBand[i].pdf->Multiply(hacc[i]);            
	

	BrBsDs  = newModelDs->Eval_BR(); 
	BrBsDsS = newModelDsS->Eval_BR();
	
	fitBand[0].yield = effRD  * BrBsDs/ BrBdD  *  fsfdBrDs / BrD * NrefD ;                                                      
	fitBand[1].yield = effRDst* BrBsDsS/BrBdDst * fsfdBrDs / BrD / BrDst * NrefDst ;                                            
	fitBand[2].yield = physBkg;                                                                                                 
	fitBand[3].yield = combBkg;                       
	
	for(auto c : fitBand)  c.pdf->Scale(c.yield/c.pdf->Integral()); 
	
	TH1D * hFit = (TH1D*) hData->Clone("hFit"); 
	hFit->Reset(); 
	for(auto c : fitBand) hFit->Add(c.pdf);
	
	for(int k = 0; k < hData->GetNbinsX(); k++) {
	  hg[k]->Fill(hFit->GetBinContent(k+1));
	}
	hFit->Delete();

	
      }
      
      for(int k = 0; k < hData->GetNbinsX(); k++) {  
	double w(0), werr(0);
	w = hData->GetBinCenter(k+1);
	werr = hData->GetBinWidth(k+1)/2.0;
	hg[k]->GetQuantiles(dim,quantiles,sum);
	gr[0]->SetPoint(k,w,hg[k]->GetMean());
	gr[0]->SetPointError(k,werr,hg[k]->GetMean()-quantiles[0]);
	gr[1]->SetPoint(k,w,hg[k]->GetMean());
	gr[1]->SetPointError(k,werr,hg[k]->GetMean()-quantiles[2]);
	gr[2]->SetPoint(k,w,hg[k]->GetMean());
	gr[2]->SetPointError(k,werr,hg[k]->GetMean()-quantiles[4]);
    
	//hg->Delete();
      }
    

  }
  else if(xname.Contains("fplus") || xname.Contains("fzero")){
	
	double Gw{0}, f0{0};
    double r =  _Dsmass/_Bmass;

    //start the toys
    int wPoints = 400;
    for(int k = 0; k < wPoints; k++) { //w points
	  if(k!=0)
		w = 1.0000001 + 1.33/wPoints*k;
	  else
		w = 1.5465;
		
      TH1D * hg = new TH1D("hg","hg",10000,0,5);
      
      for(int i = 0; i < 1000; i++) { //nExtractions per point
		tmpFFPars.clear();
		TVectorD shift(allPars.size());
		TDecompChol tdc(cov);
		if ( !tdc.Decompose() ) {cout << "Something is wrong with decomposition!" << endl; return;}
		TMatrixD u = tdc.GetU(); u.T();
		TRandom3 ran(i+1);

		for (long unsigned int m=0; m<allPars.size(); ++m) shift(m) = ran.Gaus();
		shift = u*shift;
		for (long unsigned int m=0; m<allPars.size(); ++m) {
		  for (long unsigned int n=0; n<newPars.size(); ++n)
			if(newPars[n].inum == allPars[m].inum) newPars[n].value = allPars[m].value + shift(m);
		}
		
		for (long unsigned int m=0; m<newPars.size(); ++m)
		  if(m < ffmodelDs.FFpars.size())  tmpFFPars.push_back(newPars[m]);
	
		if     (ffmodelDs.model == "CLN") decayRates::FFfunctionsCLN(w, tmpFFPars, Gw, f0);
		else if(ffmodelDs.model == "BGL") decayRates::FFfunctionsBGL(w, tmpFFPars, Gw, f0);
		else{ cout << "Wrong FF model for decFitDs in FFfunctions!"; return; }

		if(xname.Contains("fplus")) FFvalue = (1+r)/2./sqrt(r)*Gw;
		else FFvalue = f0;
		
		hg->Fill(FFvalue);
      }

      hg->GetQuantiles(dim,quantiles,sum);
      gr[0]->SetPoint(k,w,hg->GetMean());
      gr[0]->SetPointError(k,0,hg->GetMean()-quantiles[0]);
      gr[1]->SetPoint(k,w,hg->GetMean());
      gr[1]->SetPointError(k,0,hg->GetMean()-quantiles[2]);
      gr[2]->SetPoint(k,w,hg->GetMean());
      gr[2]->SetPointError(k,0,hg->GetMean()-quantiles[4]);
      
      hg->Delete();

    }
  }
  else {
    double mB = _Bmass; double mDsS = _DsSmass;
    double wBins[] = {1.000, 1.1087, 1.1688, 1.2212, 1.2717, 1.3226, 1.3814, 1.4667}; 

    //need new model to draw the bars, to avoid modifying the origianl one
    decayRates * newModel = new decayRates(decFitDsS->GetFFModel(),true);
    FFModel newFFModel;
    newFFModel.model = ffmodelDsS.model;
    newFFModel.initNFFpars(ffmodelDsS.FFpars.size());	
    newModel->SetTauB(decFitDsS->GetTauB());                                 
    newModel->SetVcbEff(decFitDsS->GetVcbEff());    
    
    if(xname.Contains("bin")){ 
      double norm  = 1.;
      double Vcb   = 1.;
      double etaEW = 1.;
      double wmin = 1.;
      double wmax = (mB*mB + mDsS*mDsS) / (2.*mB*mDsS);
      double w(0), werr(0);
      for(long unsigned int k = 0; k < 7; k++) { //w points
		w    = (wBins[k+1]+wBins[k])/2.;
		werr = (wBins[k+1]-wBins[k])/2.;

		TH1D * hg = new TH1D("hg","hg",10000,0,5);

		for(int i = 0; i < 1000; i++) { //nExtractions per point
	  
		  tmpFFPars.clear();
		  TVectorD shift(allPars.size());
		  TDecompChol tdc(cov);
		  if ( !tdc.Decompose() ) {cout << "Something is wrong with decomposition!" << endl; return;}
		  TMatrixD u = tdc.GetU(); u.T();
		  TRandom3 ran(i+1);
	  
		  for (long unsigned int m=0; m<allPars.size(); ++m) shift(m) = ran.Gaus();
		  shift = u*shift;
		  
		  for (long unsigned int m=0; m<allPars.size(); ++m) {
			for (long unsigned int n=0; n<newPars.size(); ++n)
			  if(newPars[n].inum == allPars[m].inum) newPars[n].value = allPars[m].value + shift(m);
		  }
		  
		  for (long unsigned int m=0; m<newPars.size(); ++m) {
			if(newPars[m].name == "normLHCb") norm  = newPars[m].value;
			else if(newPars[m].name == "Vcb") Vcb  = newPars[m].value;
			else if(newPars[m].name == "etaEW") etaEW  = newPars[m].value;
			//only the DsS FFs
			if(m >= ffmodelDs.FFpars.size() && m < ffmodelDsS.FFpars.size()+ffmodelDs.FFpars.size())  {
			  tmpFFPars.push_back(newPars[m]);
			}
		  }

	  for(long unsigned int m = 0; m < tmpFFPars.size(); ++m) newFFModel.FFpars[m] = tmpFFPars[m];
	  
	  ///set extracted parameters in the new FF model
	  newModel->SetFFModel(newFFModel);
	  newModel->SetVcbEff(Vcb*etaEW);

	  TF1 * f = new TF1("f",newModel,&decayRates::TF_dGdw, wmin, wmax, 0);
	  double yFF = norm*f->Integral(w-werr,w+werr);
	  double intTot = norm*f->Integral(wmin,wmax);
	  FFvalue = yFF/intTot*_normHistLHCb;
	  hg->Fill(FFvalue/(werr*2.0));
	  f->Delete();
		}

	hg->GetQuantiles(dim,quantiles,sum);
	gr[0]->SetPoint(k,w,hg->GetMean());
	gr[0]->SetPointError(k,werr,hg->GetMean()-quantiles[0]);
	gr[1]->SetPoint(k,w,hg->GetMean());
	gr[1]->SetPointError(k,werr,hg->GetMean()-quantiles[2]);
	gr[2]->SetPoint(k,w,hg->GetMean());
	gr[2]->SetPointError(k,werr,hg->GetMean()-quantiles[4]);
	
	hg->Delete();
      }
    }
    else {
	  double A1{0}, A2{0}, V{0}, A0{0};
      
      //start the toys
      int wPoints = 400;
      for(int k = 0; k < wPoints; k++) { //w points
		if(k!=0)
		  w = 1.0000001 + 1.2/wPoints*k;
		else
		  w = 1.4672;
	
	TH1D * hg = new TH1D("hg","hg",10000,0,5);
	
	for(int i = 0; i < 1000; i++) { //nExtractions per point
	  tmpFFPars.clear();
	  TVectorD shift(allPars.size());
	  TDecompChol tdc(cov);                                                                                                                        
	  if ( !tdc.Decompose() ) {cout << "Something is wrong with decomposition!" << endl; return;}
	  TMatrixD u = tdc.GetU(); u.T();                                                                                                              
	  TRandom3 ran(i+1);                                                                                                                            
	  
	  for (long unsigned int m=0; m<allPars.size(); ++m) shift(m) = ran.Gaus();    
	  shift = u*shift;
	  
	  for (long unsigned int m=0; m<allPars.size(); ++m) {
		for (long unsigned int n=0; n<newPars.size(); ++n)
		  if(newPars[n].inum == allPars[m].inum) newPars[n].value = allPars[m].value + shift(m);
	  }

	  for (long unsigned int m=0; m<newPars.size(); ++m) {
	    if(m >= ffmodelDs.FFpars.size() && m < ffmodelDsS.FFpars.size()+ffmodelDs.FFpars.size())   tmpFFPars.push_back(newPars[m]);
	  }

	  if     (newFFModel.model == "CLN")  {
		double hw{0}, R1w{0}, R2w{0}, R0w{0};
		decayRates::FFfunctionsCLN(w, tmpFFPars, hw, R1w, R2w, R0w);
		double r = mDsS/mB;
		double R = 2.*sqrt(r)/(1.+r);
		A1 = hw/2.*(w+1.)*R;
		A2 = R2w*hw/R;
		V  = R1w*hw/R;
		A0 = R0w*hw/R;
	  }
	  else if(newFFModel.model == "BGL"){
		decayRates::FFfunctionsBGL(w, tmpFFPars, A1, A2, V, A0);
	  }
	
	  if     (xname.Contains("v")) { FFvalue = V;  if(FFvalue!=FFvalue) { cout << " w = " << w << ", V  = " <<  V  << endl;} }
	  else if(xname.Contains("a1")){ FFvalue = A1; if(FFvalue!=FFvalue) { cout << " w = " << w << ", A1 = " <<  A1 << endl;} }
	  else if(xname.Contains("a2")){ FFvalue = A2; if(FFvalue!=FFvalue) { cout << " w = " << w << ", A2 = " <<  A2 << endl;} }
	  else if(xname.Contains("a0")){ FFvalue = A0; if(FFvalue!=FFvalue) { cout << " w = " << w << ", A0 = " <<  A0 << endl;} }
	  else{
		cout << "Wrong name for x point for the external input in FFfunction!" << endl;
		return;
	  }
	  hg->Fill(FFvalue);
      }

      hg->GetQuantiles(dim,quantiles,sum);
      gr[0]->SetPoint(k,w,hg->GetMean());
      gr[0]->SetPointError(k,0,hg->GetMean()-quantiles[0]);
      gr[1]->SetPoint(k,w,hg->GetMean());
      gr[1]->SetPointError(k,0,hg->GetMean()-quantiles[2]);
      gr[2]->SetPoint(k,w,hg->GetMean());
      gr[2]->SetPointError(k,0,hg->GetMean()-quantiles[4]);
      
      hg->Delete();
      
      }
    }
  }
  
  gr[0]->SetLineColor(c1);
  gr[0]->SetMarkerColor(c1);
  gr[0]->SetFillColor(c1);
  gr[0]->SetMarkerSize(0);
  gr[1]->SetLineColor(c2);
  gr[1]->SetLineColor(c2);
  gr[1]->SetMarkerColor(c2);
  gr[1]->SetFillColor(c2);
  gr[1]->SetMarkerSize(0);
  gr[2]->SetLineColor(c3);
  gr[2]->SetMarkerColor(c3);
  gr[2]->SetFillColor(c3);
  gr[2]->SetMarkerSize(0);
  
  /*
  gr[0]->SetLineColorAlpha(kOrange,0.2);
  gr[0]->SetMarkerColorAlpha(kOrange,0.2);
  gr[0]->SetFillColorAlpha(kOrange,0.2);
  gr[0]->SetMarkerSize(0);
  gr[1]->SetLineColorAlpha(kOrange,0.5);
  gr[1]->SetMarkerColorAlpha(kOrange,0.5);
  gr[1]->SetFillColorAlpha(kOrange,0.5);
  gr[1]->SetMarkerSize(0);
  gr[2]->SetLineColorAlpha(kOrange,0.8);
  gr[2]->SetMarkerColorAlpha(kOrange,0.8);
  gr[2]->SetFillColorAlpha(kOrange,0.8);
  gr[2]->SetMarkerSize(0);
  */
  return;
}


double FFfunctions(TString xname, double w, double werr) {
  
  double FFvalue = 0;
  
  if(xname.Contains("fplus") || xname.Contains("fzero")){
	double Gw{0}, f0{0};
	double r =  _Dsmass/_Bmass;
	FFModel ffmodel =  decFitDs->GetFFModel();
	if     (ffmodel.model == "CLN") decayRates::FFfunctionsCLN(w, ffmodel.FFpars, Gw, f0);
	else if(ffmodel.model == "BGL") decayRates::FFfunctionsBGL(w, ffmodel.FFpars, Gw, f0);
	else{ cout << "Wrong FF model for decFitDs in FFfunctions!"; return 1e20; }
	FFvalue = (1+r)/2./sqrt(r)*Gw;
	if(xname.Contains("fzero")) FFvalue = f0;
  }
  else {
	double mB = _Bmass; double mDsS = _DsSmass;
	if(xname.Contains("bin")){
	  double norm =1.;
	  for(auto p : otherPars) if (p.name == "normLHCb") norm = p.value;
	  double wmin = 1.;
	  double wmax = (mB*mB + mDsS*mDsS) / (2.*mB*mDsS);
	  TF1 * f = new TF1("f",decFitDsS,&decayRates::TF_dGdw, wmin, wmax, 0);
	  double yFF = norm*f->Integral(w-werr,w+werr);
	  double intTot = norm*f->Integral(wmin,wmax);
	  FFvalue = yFF/intTot*_normHistLHCb;
	  f->Delete();
	}
	else{
	  double V{0}, A1{0}, A2{0}, A0{0};
	  FFModel ffmodel = decFitDsS->GetFFModel();
	  if     (ffmodel.model == "CLN")  {
		double hw{0}, R1w{0}, R2w{0}, R0w{0};
		decayRates::FFfunctionsCLN(w, ffmodel.FFpars, hw, R1w, R2w, R0w);
		double r = mDsS/mB;
		double R = 2.*sqrt(r)/(1.+r);
		A1 = hw/2.*(w+1.)*R;
		A2 = R2w*hw/R;
		V  = R1w*hw/R;
		A0 = R0w*hw/R;
	  }
	  else if(ffmodel.model == "BGL"){
		decayRates::FFfunctionsBGL(w, ffmodel.FFpars, A1, A2, V, A0);
	  }
	  else { cout << "Wrong FF model for decFitDs in FFfunctions!"; return 1e20; }
	
	  if     (xname.Contains("v")) { FFvalue = V;  if(FFvalue!=FFvalue) { cout << " w = " << w << ", V  = " <<  V  << endl;} }
	  else if(xname.Contains("a1")){ FFvalue = A1; if(FFvalue!=FFvalue) { cout << " w = " << w << ", A1 = " <<  A1 << endl;} }
	  else if(xname.Contains("a2")){ FFvalue = A2; if(FFvalue!=FFvalue) { cout << " w = " << w << ", A2 = " <<  A2 << endl;} }
	  else if(xname.Contains("a0")){ FFvalue = A0; if(FFvalue!=FFvalue) { cout << " w = " << w << ", A0 = " <<  A0 << endl;} }
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
   

	TCanvas* cperp = new TCanvas("cperp", "cperp", 800,800);
	
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

	hData->Draw("PE");
	for(int i = 1; i <= hData->GetNbinsX(); i++) {
	  double error = sqrt(hData->GetBinError(i)*hData->GetBinError(i) + hFit->GetBinError(i)*hFit->GetBinError(i));
	  hData->SetBinError(i,error);
	}

	fitComp[1].pdf->Draw("HISTSAME][");
	fitComp[0].pdf->Draw("HISTSAME][");
	fitComp[2].pdf->Draw("HISTSAME][");
	fitComp[3].pdf->Draw("HISTSAME][");

	/*
	std::vector<TGraphErrors*> gfit;
	DrawFFErrorBand("proj",gfit);
	TMultiGraph* gfitproj = new TMultiGraph();
	gfitproj->Add(gfit[2],"p2][");
	gfitproj->Add(gfit[1],"p2][");
	gfitproj->Add(gfit[0],"p2][");
	gfitproj->Draw("a][");
	*/
	hFit->Draw("HISTSAME][");
	hData->Draw("PESAME");
  
	TLegend * leg;
	leg = new TLegend(0.22,0.75,0.45,0.9);
	leg->SetFillStyle(0);
	leg->SetTextSize(0.045);
	leg->SetBorderSize(0);
	leg->SetTextFont(132);
	leg->AddEntry(hData,"LHCb, PRD101 (2020) 072004","lpe");
	leg->AddEntry(hFit,"Fit","l");
	leg->Draw("SAME");
    
	lowerPad->cd();
	TH1D * hpull = (TH1D*)hData->Clone("hpull");
	hpull->Reset();
    DrawResiduals(hData,hFit,hpull);
    for(int i = 1; i < hpull->GetNbinsX(); i++) {
      hpull->SetBinContent(i, (hData->GetBinContent(i)-hFit->GetBinContent(i))/hData->GetBinError(i));
    }
    cperp->SaveAs("fit_projection_pperp_Ds_"+FFModelFitDs.model+"_DsS_"+FFModelFitDsS.model+".pdf");
    cperp->SaveAs("fit_projection_pperp_Ds_"+FFModelFitDs.model+"_DsS_"+FFModelFitDsS.model+".C");

  TString fmodels[]  ={"HPQCD_Ds","MILC","LCSRDs"};
  const int Nfmodels = sizeof(fmodels)/sizeof(fmodels[0]);
  int fcolor[]={1,2,13};
  int marker[]={21,22,23};
  const int nFFDs = FFModelFitDs.model=="CLN"? 1 : 2;

  TGraphErrors* grfplus[nFFDs][Nfmodels];
  TGraphErrors* grfpluspull[nFFDs][Nfmodels];
  for(int im=0; im<Nfmodels; ++im){
	for(int i=0;i<nFFDs;++i) {
	  grfplus[i][im] = new TGraphErrors();
	  grfplus[i][im]->SetMarkerColor(fcolor[im]);
	  grfplus[i][im]->SetLineColor(fcolor[im]);
	  grfplus[i][im]->SetMarkerStyle(marker[im]);
	  grfplus[i][im]->SetMarkerSize(1.1);
	  grfpluspull[i][im] = new TGraphErrors();
	  grfpluspull[i][im]->SetMarkerColor(fcolor[im]);
	  grfpluspull[i][im]->SetLineColor(fcolor[im]);
	  grfpluspull[i][im]->SetMarkerStyle(marker[im]);
	  grfpluspull[i][im]->SetMarkerSize(1);
	}
	
	TString n;
	double x, y, err;
	int ip = 0; int ip0 = 0;
	ifstream ftmp(dirInputs+"/"+fmodels[im]+"_values.txt");
	if(!ftmp.is_open()){
	  cout << fmodels[im] << " input file not found!" << endl; return false;
	}
	while(ftmp >> n >> x >> y >> err) {
	  if(n.Contains("fplus")){
		grfplus[0][im]->SetPoint(ip,x,y);
		grfplus[0][im]->SetPointError(ip,0,err);
		grfpluspull[0][im]->SetPoint(ip,x,(y-FFfunctions("fplus",x,0))/err);
		grfpluspull[0][im]->SetPointError(ip,0,0);
		ip++;
	  }else if(n.Contains("fzero") && nFFDs>1){
		grfplus[1][im]->SetPoint(ip0,x,y);
		grfplus[1][im]->SetPointError(ip0,0,err);
		grfpluspull[1][im]->SetPoint(ip0,x,(y-FFfunctions("fzero",x,0))/err);
		grfpluspull[1][im]->SetPointError(ip0,0,0);
		ip0++;
	  }
	}
  }
  
  TGraph* grfplusFit[nFFDs];
  TCanvas* cf[nFFDs];
  TPad*    upperPadcf[nFFDs];
  TPad*    lowerPadcf[nFFDs];
  TMultiGraph* gf[nFFDs];
  TLegend * legf[nFFDs];
  TLine * l1gf[nFFDs]; TLine * l2gf[nFFDs]; TLine * l3gf[nFFDs];
  TString titleyf[]={"#it{f}_{+}(#it{w})","#it{f}_{0}(#it{w})"};
  TString namef[]={"fplus","fzero"};
  TString fmodelsName[]  ={"HPQCD, PRD101 (2020) 7, 074513","MILC, PRD85 (2012) 114502","LCSR, EPJC80 (2020) 4, 347"};
  TMultiGraph* gfpull[nFFDs];
  for(int i=0;i<nFFDs;++i){
	grfplusFit[i] = new TGraph();
	grfplusFit[i]->SetLineWidth(2);
	grfplusFit[i]->SetLineColor(9);
  }
  
  for(int p=0;p<1000;++p){
	double w = 1.0 + 1.33/1000*p;
	grfplusFit[0]->SetPoint(p,w,FFfunctions("fplus",w,0));
	if(nFFDs>1)
	  grfplusFit[1]->SetPoint(p,w,FFfunctions("fzero",w,0));
  }
  
  for(int i=0;i<nFFDs;++i){
	cf[i] = new TCanvas(Form("cf_%i",i),Form("cf_%i",i),800,800);
	cf[i]->cd();
	//cf->SetLeftMargin(0.2);
	//cf->SetRightMargin(0.1);
	upperPadcf[i] = new TPad(Form("upperPadcf_%i",i), Form("upperPadcf_%i",i),   .005, .2525, .995, .995);
	lowerPadcf[i] = new TPad(Form("lowerPadcf_%i",i), Form("lowerPadcf_%i",i),   .005, .005, .995, .2475);
	upperPadcf[i]->SetLeftMargin(0.2);
	upperPadcf[i]->SetRightMargin(0.1);
	lowerPadcf[i]->SetLeftMargin(0.2);
	lowerPadcf[i]->SetRightMargin(0.1);
	upperPadcf[i]->Draw();
	lowerPadcf[i]->Draw();
	upperPadcf[i]->cd();

	gf[i] = new TMultiGraph();
	  
	std::vector<TGraphErrors*> gfCL;
	if(i==0) {DrawFFErrorBand("fplus",gfCL);
	  cout << "f+(q2=0) = " << gfCL[0]->GetPointY(0) << " +- " << gfCL[0]->GetErrorY(0) << ", w = " << gfCL[0]->GetPointX(0) << endl;
	}
	else { DrawFFErrorBand("fzero",gfCL);
	  cout << "f0(q2=0) = " << gfCL[0]->GetPointY(0) << " +- " << gfCL[0]->GetErrorY(0) << ", w = " << gfCL[0]->GetPointX(0) << endl;
	}
	  
	gf[i]->Add(gfCL[2],"p");
	gf[i]->Add(gfCL[1],"p");
	gf[i]->Add(gfCL[0],"p");
	for(int im=0; im<Nfmodels; ++im) gf[i]->Add(grfplus[i][im],"p");
	gf[i]->GetXaxis()->SetLimits(0.9,2.4);
	gf[i]->SetMinimum(0);
	gf[i]->SetMaximum(1.8);
	gf[i]->GetYaxis()->SetTitle(titleyf[i]);
	gf[i]->GetXaxis()->SetTitle("#it{w}");
	gf[i]->Draw("a");
	legf[i] = new TLegend(0.39,0.75,0.59,0.95);
	legf[i]->SetFillStyle(0);
	legf[i]->SetTextSize(0.045);
	legf[i]->SetBorderSize(0);
	legf[i]->SetTextFont(132);
	
	for(int im=0; im<Nfmodels; ++im) legf[i]->AddEntry(grfplus[i][im],fmodelsName[im],"pe");
	//legf->AddEntry(grfplusFit,"Fit","l");
	legf[i]->Draw("SAME");

	lowerPadcf[i]->cd();
	gfpull[i] = new TMultiGraph();
	for(int im=0; im<Nfmodels; ++im) gfpull[i]->Add(grfpluspull[i][im],"p");
	gfpull[i]->GetXaxis()->SetLimits(0.9,2.4);
	gfpull[i]->GetHistogram()->GetYaxis()->SetNdivisions(505);
	gfpull[i]->GetHistogram()->GetYaxis()->SetLabelSize(0.13);
	gfpull[i]->GetHistogram()->GetYaxis()->SetRangeUser(-5,5);
	gfpull[i]->GetHistogram()->GetXaxis()->SetLabelSize(0);
	gfpull[i]->GetHistogram()->GetYaxis()->SetTitle("Pulls");
	gfpull[i]->GetHistogram()->GetYaxis()->SetTitleSize(0.15);
	gfpull[i]->GetHistogram()->GetYaxis()->SetTitleOffset(0.4);
	gfpull[i]->GetHistogram()->GetXaxis()->SetTitleSize(0);

	l1gf[i] = new TLine(gfpull[i]->GetHistogram()->GetXaxis()->GetXmin(),-3,gfpull[i]->GetHistogram()->GetXaxis()->GetXmax(),-3);
	l2gf[i] = new TLine(gfpull[i]->GetHistogram()->GetXaxis()->GetXmin(),+3,gfpull[i]->GetHistogram()->GetXaxis()->GetXmax(),+3);
	l3gf[i] = new TLine(gfpull[i]->GetHistogram()->GetXaxis()->GetXmin(),0,gfpull[i]->GetHistogram()->GetXaxis()->GetXmax(),0);
	l1gf[i]->SetLineStyle(kDashed);
	l2gf[i]->SetLineStyle(kDashed);
	l3gf[i]->SetLineStyle(kDashed);
	gfpull[i]->Draw("a");
	l1gf[i]->Draw("SAME");
	l2gf[i]->Draw("SAME");
	l3gf[i]->Draw("SAME");

	cf[i]->SaveAs("fit_projection_"+namef[i]+"_Ds_"+FFModelFitDs.model+"_DsS_"+FFModelFitDsS.model+".pdf");
	cf[i]->SaveAs("fit_projection_"+namef[i]+"_Ds_"+FFModelFitDs.model+"_DsS_"+FFModelFitDsS.model+".C");
  }
  
  TString VAmodels[] ={"HPQCD_DsS","LCSRDsS"};
  const int NVAmodels = sizeof(VAmodels)/sizeof(VAmodels[0]);
  int colorsVA[]={1,13};
  int markersVA[]={21,23};
  const int nFFDsS = FFModelFitDsS.model=="CLN"? 3 : 4;
  TGraphErrors* grDst[nFFDsS][NVAmodels];
  TGraphErrors* grDstpull[nFFDsS][NVAmodels];
  for(int im=0; im<NVAmodels; ++im){
	for(int i=0;i<nFFDsS;++i) {
	  grDst[i][im] = new TGraphErrors();
	  grDst[i][im]->SetLineWidth(2);
	  grDst[i][im]->SetMarkerColor(colorsVA[im]);
	  grDst[i][im]->SetLineColor(colorsVA[im]);
	  grDst[i][im]->SetMarkerStyle(markersVA[im]);
	  grDst[i][im]->SetMarkerSize(1.1);
	  grDstpull[i][im] = new TGraphErrors();
	  grDstpull[i][im]->SetLineWidth(2);
	  grDstpull[i][im]->SetMarkerColor(colorsVA[im]);
	  grDstpull[i][im]->SetLineColor(colorsVA[im]);
	  grDstpull[i][im]->SetMarkerStyle(markersVA[im]);
	  grDstpull[i][im]->SetMarkerSize(1);
	}
	
	TString n;
	double x, y, err;
	int ipv{0},ipa1{0},ipa2{0},ipa0{0};
	ifstream ftmp(dirInputs+"/"+VAmodels[im]+"_values.txt");
	if(!ftmp.is_open()){
	  cout << VAmodels[im] << " input file not found!" << endl; return false;
	}
	while(ftmp >> n >> x >> y >> err) {
	  //cout << n << x << y << err<< endl;
	  if     (n.Contains("v") ){ 
	    grDst[0][im]->SetPoint(ipv,x,y);  grDst[0][im]->SetPointError(ipv,0,err);  
	    grDstpull[0][im]->SetPoint(ipv,x,(y-FFfunctions("v",x,0))/err);  grDstpull[0][im]->SetPointError(ipv,0,0);
	    ipv++; 
	  }
	  else if(n.Contains("a1")){ 
	    grDst[1][im]->SetPoint(ipa1,x,y); grDst[1][im]->SetPointError(ipa1,0,err); 
	    grDstpull[1][im]->SetPoint(ipa1,x,(y-FFfunctions("a1",x,0))/err); grDstpull[1][im]->SetPointError(ipa1,0,0); 
	    ipa1++;
	  }
	  else if(n.Contains("a2")){ 
	    grDst[2][im]->SetPoint(ipa2,x,y); grDst[2][im]->SetPointError(ipa2,0,err); 
	    grDstpull[2][im]->SetPoint(ipa2,x,(y-FFfunctions("a2",x,0))/err); grDstpull[2][im]->SetPointError(ipa2,0,0); 
	    ipa2++;
	  }
	  else if(n.Contains("a0") && nFFDsS>3 ){
		grDst[3][im]->SetPoint(ipa0,x,y); grDst[3][im]->SetPointError(ipa0,0,err);
		grDstpull[3][im]->SetPoint(ipa0,x,(y-FFfunctions("a0",x,0))/err); grDstpull[3][im]->SetPointError(ipa0,0,0);
		ipa0++;
	  }
	}
  }
  
  TGraph* grDstFit[nFFDsS];
  TString titley[]={"#it{V}(#it{w})","#it{A}_{1}(#it{w})","#it{A}_{2}(#it{w})","#it{A}_{0}(#it{w})"};
  for(int i=0;i<nFFDsS;++i){
	grDstFit[i] = new TGraph();
	grDstFit[i]->SetLineWidth(2);
	grDstFit[i]->SetLineColor(9);
  }
  for(int p=0;p<1000;++p){
	double w = 1. + 1.2/1000*p;
	grDstFit[0]->SetPoint(p,w,FFfunctions("v" ,w,0));
	grDstFit[1]->SetPoint(p,w,FFfunctions("a1",w,0));
	grDstFit[2]->SetPoint(p,w,FFfunctions("a2",w,0));
	if (nFFDsS>3)
	  grDstFit[3]->SetPoint(p,w,FFfunctions("a0",w,0));
  }
  
  TGraphErrors* grDstHPQCD = new TGraphErrors();
  grDstHPQCD->SetLineWidth(2);
  grDstHPQCD->SetMarkerColor(4);
  grDstHPQCD->SetLineColor(4);
  grDstHPQCD->SetMarkerStyle(20);
  grDstHPQCD->SetMarkerSize(1.1);
  TGraphErrors* grDstHPQCDpull = new TGraphErrors();
  grDstHPQCDpull->SetLineWidth(2);
  grDstHPQCDpull->SetMarkerColor(4);
  grDstHPQCDpull->SetLineColor(4);
  grDstHPQCDpull->SetMarkerStyle(20);
  grDstHPQCDpull->SetMarkerSize(1);
  double mB = _Bmass; double mDsS = _DsSmass;
  double hpqcdDsS  = 1./(mB+mDsS)*0.902* sqrt(mB*mDsS)*(2.);
  double ehpqcdDsS = 1./(mB+mDsS)*0.013 * sqrt(mB*mDsS)*(2.);
  grDstHPQCD->SetPoint(0,1,hpqcdDsS);
  grDstHPQCD->SetPointError(0,0,ehpqcdDsS);
  grDstHPQCDpull->SetPoint(0,1,(hpqcdDsS-FFfunctions("a1",1,0))/ehpqcdDsS);
  grDstHPQCDpull->SetPointError(0,0,0);
 
  TMultiGraph* gVA[nFFDsS];
  TMultiGraph* gVApull[nFFDsS];
  TCanvas* cVA[nFFDsS];
  TPad * upperPadVA[nFFDsS];
  TPad * lowerPadVA[nFFDsS];
  TLegend * legVA[nFFDsS];
  TString funcname[]={"V","A1","A2","A0"};
  for(int cp=0;cp<nFFDsS;++cp){
        std::vector<TGraphErrors*> gVACL;
        cVA[cp] = new TCanvas(Form("c_%i",cp),Form("c_%i",cp),800,800);
	//cVA[cp]->SetLeftMargin(0.2);
	//cVA[cp]->SetRightMargin(0.1);
	upperPadVA[cp] = new TPad(Form("upperPad_%i",cp), Form("upperPad_%i",cp),   .005, .2525, .995, .995);   
	lowerPadVA[cp] = new TPad(Form("lowerPad_%i",cp), Form("lowerPad_%i",cp),   .005, .005, .995, .2475);   
	upperPadVA[cp]->SetLeftMargin(0.2);
	upperPadVA[cp]->SetRightMargin(0.1);
	lowerPadVA[cp]->SetLeftMargin(0.2);
	lowerPadVA[cp]->SetRightMargin(0.1);
	upperPadVA[cp]->Draw();
	lowerPadVA[cp]->Draw();
	upperPadVA[cp]->cd();
	gVA[cp] = new TMultiGraph();
	gVApull[cp] = new TMultiGraph();
	
	if(funcname[cp] == "V"){
	  DrawFFErrorBand("v",gVACL);
	  cout << "V(q2=0) = " << gVACL[0]->GetPointY(0) << " +- " << gVACL[0]->GetErrorY(0) << ", w = " << gVACL[0]->GetPointX(0) << endl;
	}
	if(funcname[cp] == "A1"){
	  DrawFFErrorBand("a1",gVACL);
	  cout << "A1(q2=0) = " << gVACL[0]->GetPointY(0) << " +- " << gVACL[0]->GetErrorY(0) << ", w = " << gVACL[0]->GetPointX(0) << endl;
	}
	if(funcname[cp] == "A2"){
	  DrawFFErrorBand("a2",gVACL);
	  cout << "A2(q2=0) = " << gVACL[0]->GetPointY(0) << " +- " << gVACL[0]->GetErrorY(0) << ", w = " << gVACL[0]->GetPointX(0) << endl;
	}
	if(funcname[cp] == "A0"){
	  DrawFFErrorBand("a0",gVACL);
	  cout << "A0(q2=0) = " << gVACL[0]->GetPointY(0) << " +- " << gVACL[0]->GetErrorY(0) << ", w = " << gVACL[0]->GetPointX(0) << endl;
	}

	gVA[cp]->Add(gVACL[2],"p");
	gVA[cp]->Add(gVACL[1],"p");
	gVA[cp]->Add(gVACL[0],"p");
  
	//gVA[cp]->Add(grDstFit[cp],"l");
	for(int m=0; m<NVAmodels; m++) {
	  gVA[cp]->Add(grDst[cp][m],"p");
	  gVApull[cp]->Add(grDstpull[cp][m],"p");	
	}
	if(cp==1) {
	  gVA[cp]->Add(grDstHPQCD,"p");
	  gVApull[cp]->Add(grDstHPQCDpull,"p");
	}
	gVA[cp]->GetXaxis()->SetLimits(0.9,2.3);
	gVA[cp]->GetHistogram()->GetYaxis()->SetRangeUser(0,1.5);
	gVApull[cp]->GetXaxis()->SetLimits(0.9,2.3);
	gVApull[cp]->GetHistogram()->GetYaxis()->SetRangeUser(-5,5);
	if(cp==0) gVA[cp]->GetHistogram()->GetYaxis()->SetRangeUser(0,2.1);
	gVA[cp]->GetYaxis()->SetTitle(titley[cp]);
	gVA[cp]->GetXaxis()->SetTitle("#it{w}");
	gVA[cp]->Draw("a");
	legVA[cp] = new TLegend(0.45,0.75,0.65,0.95);
	legVA[cp]->SetFillStyle(0);
	legVA[cp]->SetTextSize(0.045);
	legVA[cp]->SetBorderSize(0);
	legVA[cp]->SetTextFont(132);
	if(cp==0) gVA[cp]->GetHistogram()->GetYaxis()->SetRangeUser(0,2.1);
	if(cp==1) legVA[cp]->AddEntry(grDstHPQCD,"HPQCD, PRD99 (2019) 114512","pe");
	legVA[cp]->AddEntry(grDst[cp][0],"HPQCD, arXiv:2105.11433","pe");
	legVA[cp]->AddEntry(grDst[cp][1],"LCSR, EPJC80 (2020) 4, 347","pe");
	//legVA[cp]->AddEntry(grDstFit[cp],"Fit","l");
	legVA[cp]->Draw("SAME");
	
	lowerPadVA[cp]->cd();
	gVApull[cp]->GetXaxis()->SetLimits(0.9,2.3);
	gVApull[cp]->GetHistogram()->GetYaxis()->SetNdivisions(505);
	gVApull[cp]->GetHistogram()->GetYaxis()->SetLabelSize(0.13);
	gVApull[cp]->GetHistogram()->GetYaxis()->SetRangeUser(-5,5);
	gVApull[cp]->GetHistogram()->GetXaxis()->SetLabelSize(0);
	gVApull[cp]->GetHistogram()->GetYaxis()->SetTitle("Pulls");
	gVApull[cp]->GetHistogram()->GetYaxis()->SetTitleSize(0.15);
	gVApull[cp]->GetHistogram()->GetYaxis()->SetTitleOffset(0.4);
	gVApull[cp]->GetHistogram()->GetXaxis()->SetTitleSize(0);
	
	TLine * l1gf = new TLine(gVApull[cp]->GetHistogram()->GetXaxis()->GetXmin(),-3,gVApull[cp]->GetHistogram()->GetXaxis()->GetXmax(),-3);
	TLine * l2gf = new TLine(gVApull[cp]->GetHistogram()->GetXaxis()->GetXmin(),+3,gVApull[cp]->GetHistogram()->GetXaxis()->GetXmax(),+3);
	TLine * l3gf = new TLine(gVApull[cp]->GetHistogram()->GetXaxis()->GetXmin(),0,gVApull[cp]->GetHistogram()->GetXaxis()->GetXmax(),0);
	l1gf->SetLineStyle(kDashed);
	l2gf->SetLineStyle(kDashed);
	l3gf->SetLineStyle(kDashed);
	gVApull[cp]->Draw("a");
	l1gf->Draw("SAME");
	l2gf->Draw("SAME");
	l3gf->Draw("SAME");

	cVA[cp]->SaveAs("fit_projection_"+funcname[cp]+"_Ds_"+FFModelFitDs.model+"_DsS_"+FFModelFitDsS.model+".pdf");
	cVA[cp]->SaveAs("fit_projection_"+funcname[cp]+"_Ds_"+FFModelFitDs.model+"_DsS_"+FFModelFitDsS.model+".C");
  }
  
  
  double wBins[] = {1.000, 1.1087, 1.1688, 1.2212, 1.2717, 1.3226, 1.3814, 1.4667};
  TH1D* hw = new TH1D("h_w","h_w",7,wBins);//wLHCb
  TH1D* hwFit = (TH1D*) hw->Clone("h_wFit");
  TString n;
  double x, y, xerr, yerr;
  int ip = 1;
  ifstream ftmp(dirInputs+"/"+"LHCb-PAPER-2019-046_values.txt");
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
  
  TCanvas* cw = new TCanvas("cw", "cw", 800,800);
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
  hw->SetMarkerSize(1.1);
  hw->GetYaxis()->SetRangeUser(0.5,3.5);
  hw->GetYaxis()->SetTitle("(1/#it{N})(d#it{N}/d#it{w})");
  hw->GetXaxis()->SetTitle("#it{w}");

  std::vector<TGraphErrors*> gwbin;
  DrawFFErrorBand("bin",gwbin);
  TMultiGraph* gw = new TMultiGraph();

  gw->Add(gwbin[2],"p2][");
  gw->Add(gwbin[1],"p2][");
  gw->Add(gwbin[0],"p2][");
  gw->GetHistogram()->GetYaxis()->SetRangeUser(0,4);
  gw->GetXaxis()->SetLimits(1.0,1.4667);
  gw->Draw("a][");
  //hw->Draw("PESAME");
  //gwbin[2]->Draw("SAMEP2");
  //gwbin[1]->Draw("SAMEP2");
  //gwbin[0]->Draw("SAMEP2");
  hw->Draw("PESAME");
  hwFit->SetLineWidth(2);
  hwFit->SetLineColor(9);
  //hwFit->Draw("HISTSAME][");
  hw->Draw("PESAME");
  
  TLegend * legw;
  legw = new TLegend(0.23,0.75,0.48,0.90);
  legw->SetFillStyle(0);
  legw->SetTextSize(0.045);
  legw->SetBorderSize(0);
  legw->SetTextFont(132);
  legw->AddEntry(hw,"LHCb, JHEP12 (2020) 144","lpe");
  //legw->AddEntry(hwFit,"Fit","l");
  legw->Draw("SAME");

  gPad->RedrawAxis(); 
  
  lowerPadw->cd();
  TH1D * hpullw = (TH1D*) hw->Clone("hpullw");
  hpullw->Reset();
  DrawResiduals(hw,hwFit,hpullw);
  
  cw->SaveAs("fit_projection_wLHCb_Ds_"+FFModelFitDs.model+"_DsS_"+FFModelFitDsS.model+".pdf");
  cw->SaveAs("fit_projection_wLHCb_Ds_"+FFModelFitDs.model+"_DsS_"+FFModelFitDsS.model+".C");
  
  return true;
	
}


#include "fitter.h"
#include <TVirtualFitter.h>                        
#include <TFitter.h>  
#include <TMinuit.h>
#include <TRandom3.h>
#include <TGraphErrors.h>
#include <TF2.h>
#include <TF3.h>
#include <map>
#include "TheoryInputs.h"

using namespace std;

bool fitter::RunFitter(TString filename)  {
  
    if (!Configure(filename)) return false;
    if (!DoFit(1)) return false;
    
    TCanvas* c1 = new TCanvas("c1","c1",600,1800);
    c1->Divide(theoryInputs.size()+1,1);
    Plot(c1);
     
    return true;
}


bool fitter::Configure(TString filename) {
    
  cout << "==========================================================\n";
  cout << " Configuration file read: \n";
  cout << "----------------------------------------------------------\n";

  //Take data and all histograms for the fit (hardocded)
  SetDataAndBkg();
  
  //Read the configuration file to set the fit
  if(!ReadConfigFile(filename)) return false;
    
  //Print all information read from the config
  PrintConfigInfo();
    
  //Store the candidates for the templates
  if(!FillMCcandidates(filename)) return false;
    
  return true;
    
}

void fitter::SetDataAndBkg(){
       
    TFile * fileData = TFile::Open("LHCb_input_data_mc.root");
    hData   = (TH1D*) fileData->Get("h_lhcb_data");
    fitComp[2].pdf = (TH1D*) fileData->Get("h_phys_bkg");
    fitComp[3].pdf = (TH1D*) fileData->Get("h_comb_bkg");
    
    hData->SetDirectory(0);
    fitComp[2].pdf->SetDirectory(0);
    fitComp[3].pdf->SetDirectory(0);
 
    fileData->Close();
   
    fitComp[2].yield = fitComp[2].pdf->Integral();
    fitComp[3].yield = fitComp[3].pdf->Integral();
    
    return;
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
	cout << " Initial parameters give BR(Bs->Ds): " <<  decFitDs->Eval_BR() << ", BR(Bs->Ds*):" << decFitDsS->Eval_BR() << endl;
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
    cout << " Theory inputs: ";
    for(long unsigned int i = 0; i < theoryInputs.size(); i++) cout << theoryInputs[i] << " ";
    cout << endl;
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
	else cout << hacc[i]->GetName() << ", which has min and max values:" << endl;
	
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
  
  FFModel FFModelFitDs  =  decFitDs->GetFFModel();
  FFModel FFModelFitDsS =  decFitDsS->GetFFModel();
  cout << " Fitted FF models are: \n";
  cout << " - Bs->Ds :  " << FFModelFitDs.model << ", parameters: " << FFModelFitDs.FFpars.size() << endl;
  for(auto p : FFModelFitDs.FFpars) cout << "\t " << p.name << "\t " << p.value << "\n";
  cout << " - Bs->Ds*:  " << FFModelFitDsS.model << ", parameters: " << FFModelFitDsS.FFpars.size() << endl;
  for(auto p : FFModelFitDsS.FFpars) cout << "\t " << p.name << "\t " << p.value << "\n";
  cout << " Which give BR(Bs->Ds): " <<  decFitDs->Eval_BR() << ", BR(Bs->Ds*):" << decFitDsS->Eval_BR() << endl;

  double negErr, posErr;              
  for(auto p : fitPars)
    fitter->mnpout(p.inum,p.name,p.value,p.error,negErr,posErr,ierflg);

  fitter->Delete();  

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



   //add theory inputs to chi2
   for(unsigned int i = 0; i < theoryInputs.size(); i++) 
     addTheoryInputDs(theoryInputs[i],decFitDs->GetFFModel().FFpars,_chi2,_ndf);

 
  f = _chi2;
  
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
        physBkg{0}, combBkg{0};
    
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



  bool fitter::Plot(TCanvas *c) {
	
	c->cd(1);
	TPad*    upperPad = new TPad("upperPad", "upperPad",   .005, .2525, .995, .995);
	TPad*    lowerPad = new TPad("lowerPad", "lowerPad",   .005, .005,  .995, .2475);
	upperPad->SetBottomMargin(0.18);
	upperPad->SetLeftMargin(0.15);
	lowerPad->SetLeftMargin(0.15);
	upperPad->SetRightMargin(0.2);
	lowerPad->SetRightMargin(0.2);
	upperPad->SetTopMargin(0.1);
	lowerPad->Draw();
	upperPad->Draw();
	upperPad->cd();
	hData->SetMarkerStyle(20);
	hData->SetMarkerColor(kBlack);
	hData->SetLineColor(kBlack);
	hData->SetMarkerSize(0.7);
	hData->Draw("PE");
	fitComp[0].pdf->SetLineColor(kRed);
	fitComp[0].pdf->Draw("HISTSAME][");
	fitComp[1].pdf->SetLineColor(kGreen);
	fitComp[1].pdf->Draw("HISTSAME][");
	fitComp[2].pdf->SetLineColor(kCyan);
	fitComp[2].pdf->Draw("HISTSAME][");
	fitComp[3].pdf->SetLineColor(kOrange);
	fitComp[3].pdf->Draw("HISTSAME][");
	
	TH1D * hFit = (TH1D*)fitComp[0].pdf->Clone("hFit");
	hFit->Add(fitComp[1].pdf);
	hFit->Add(fitComp[2].pdf);
	hFit->Add(fitComp[3].pdf);
	hFit->SetLineColor(kBlue);
	hFit->Draw("HISTSAME][");
	
	
	TH1D * hpull = (TH1D*)hData->Clone("hpull");
	hpull->Reset();
	
	for(int i = 1; i <= hData->GetNbinsX() ; i++) {
	  double pull = (hData->GetBinContent(i)-hFit->GetBinContent(i))/sqrt(hData->GetBinError(i)*hData->GetBinError(i) +
									      hFit->GetBinError(i)*hFit->GetBinError(i));
	  cout << i << " " << pull << endl;
	  hpull->SetBinContent(i,pull);
	}
	
	lowerPad->cd();
	hpull->SetLineColor(kBlack);
	hpull->SetMarkerColor(kBlack);
	hpull->GetYaxis()->SetNdivisions(505);
	hpull->GetYaxis()->SetLabelSize(0.13);
	hpull->GetYaxis()->SetRangeUser(-5,5);
	hpull->GetXaxis()->SetLabelSize(0);
	hpull->GetYaxis()->SetTitle("Pulls");
	hpull->GetYaxis()->SetTitleSize(0.15);
	hpull->GetYaxis()->SetTitleOffset(0.4);
	hpull->GetXaxis()->SetTitleSize(0);
	
	
	TLine * l1 = new TLine(0.2,-3,hData->GetXaxis()->GetXmax(),-3);
	TLine * l2 = new TLine(0.2,+3,hData->GetXaxis()->GetXmax(),+3);
	TLine * l3 = new TLine(0.2,0,hData->GetXaxis()->GetXmax(),0);
	l1->SetLineStyle(kDashed);
	l2->SetLineStyle(kDashed);
	l3->SetLineStyle(kDashed);
	
	hpull->GetYaxis()->SetRangeUser(-5,5);
	hpull->SetFillColorAlpha(kBlue, 0.35);
	hpull->Draw("F");
	l1->Draw("SAME");
	l2->Draw("SAME");
	l3->Draw("SAME");
	
	upperPad->cd();
	TLegend * leg;
	leg = new TLegend(0.60,0.65,0.80,0.85);
	
	leg->SetFillStyle(0);
	leg->SetTextSize(0.050);
	leg->SetBorderSize(0);
	leg->SetTextFont(132);
	leg->AddEntry(hData,"Back. sub. data","lpe");
	//leg->AddEntry(fitComp[0].pdf,"Fit","F");
	
	leg->Draw("SAME");
	

	//////////
	std::vector<TGraphErrors*> gr;
	for(long unsigned int k = 0; k < theoryInputs.size(); k++) {
	  std::vector<double> wtmp;
	  std::vector<double> ftmp;
	  std::vector<double> ftmperr;
	  if(theoryInputs[k] == "MILC") {
	    wtmp = wMILC;
	    ftmp = fMILC;
	    ftmperr = fMILCerr;
	  }
	  if(theoryInputs[k] == "HPQCD") {
	    wtmp = wHPQCD;
	    ftmp = fHPQCD;
	    ftmperr = fHPQCDerr;
	  }
	  
	  TGraphErrors * gr1 = new TGraphErrors(wtmp.size());
	  for(long unsigned int i = 0; i < wtmp.size(); i++) {
	    gr1->SetPoint(i,wtmp[i],ftmp[i]);
	    gr1->SetPointError(i,0,ftmperr[i]);
	  }
	  gr.push_back(gr1);
	}
	
	TGraphErrors * gr2 = new TGraphErrors(100);
	for(long unsigned int i = 0; i < 100; i++) {
	  double step = (1.2 - 1)/100;
	  //decFitDs->GetFFModel().FFpars
	  gr2->SetPoint(i,1+i*step,FFfunctionsCLN(1+i*step,decFitDs->GetFFModel().FFpars));
	  gr2->SetPointError(i,0,0);
	}
	

	gr2->SetMarkerStyle(20);
	gr2->SetMarkerSize(0.5);
	gr2->GetYaxis()->SetRangeUser(0.8,1.4);

	for(long unsigned int i = 0; i < gr.size(); i++) {
	  c->cd(2+i);
	  gr2->Draw("APE");
	  gr[i]->SetMarkerStyle(21+i);
	  gr[i]->SetMarkerColor(2+i);
	  gr[i]->Draw("PESAME");
	}
  
	c->SaveAs("figs/projection_fit_simul_4D.pdf");
	
	return true;
	
  }


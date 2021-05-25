#include "Tools.h"
#include "decayRates.h"

void generateSignalTemplates(TString filename = "config/fit_testCLN.cfg",
                             unsigned int nentries = 5e5,
                             bool rateDS_4D = true,
                             bool _res2D = true)
{
 
    cout << "---------------------------------------------------------" << endl;
    cout << " Will generate signal templates with entries " << nentries << endl;
    
    decayRates* decRefDs;
    decayRates* decRefDsS;
    
    string line;
    unsigned int np=1;
    ifstream fin(filename);
    TString outfile, outtree;
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
	  else if(line=="Templates")  {
		fin >> outfile >> outtree;
	  }
      else continue;
    }
    fin.close();
    
    FFModel FFModelRefDs = decRefDs->GetFFModel();
    FFModel FFModelRefDsS = decRefDsS->GetFFModel();
    cout << " REFERENCE FF models of the templates: \n";
    cout << " - Bs->Ds :  " << FFModelRefDs.model << ", parameters: " << FFModelRefDs.FFpars.size() << endl;
    for(auto p : FFModelRefDs.FFpars) cout << "\t " << p.name << "\t " << p.value << "\n";
    cout << " - Bs->Ds*:  " << FFModelRefDsS.model << ", parameters: " << FFModelRefDsS.FFpars.size() << endl;
    for(auto p : FFModelRefDsS.FFpars) cout << "\t " << p.name << "\t " << p.value << "\n";
    
    double _Bmass = Mass::Bs*1e-3;
    double _Dsmass = Mass::Ds*1e-3;
    double _DsSmass = Mass::DsS*1e-3;
    double _Mumass = Mass::Mu*1e-3;
    double wmin = 1.0;
    double wmax  = (_Bmass*_Bmass + _Dsmass*_Dsmass - _Mumass*_Mumass) / (2.*_Bmass*_Dsmass);
    double wmaxS = (_Bmass*_Bmass + _DsSmass*_DsSmass - _Mumass*_Mumass) / (2.*_Bmass*_DsSmass);
     
    double rateRefDs  = decRefDs->Eval_Gamma();
    double rateRefDsS = decRefDsS->Eval_Gamma();
    cout << " Calculated reference Gamma for P and V decays: " << rateRefDs << ", " << rateRefDsS << endl;
    cout << " which give BRs: " << decRefDs->Eval_BR() << ", " << decRefDsS->Eval_BR()<< endl;
     
    TF1 * fDs = new TF1("fDs",decRefDs,&decayRates::TF_dGdw, wmin, wmax, 0);
    fDs->SetNpx(10000);
    TF1 * fDsS = new TF1("fDsS",decRefDsS,&decayRates::TF_dGdw, wmin, wmaxS, 0);
    fDsS->SetNpx(10000);
    
    TF2 * fDsS2D = new TF2("fDsS2D",decRefDsS,&decayRates::TF_dGdwdctd, wmin, wmaxS,-1.,1.,0);
    fDsS2D->SetNpx(1000); fDsS2D->SetNpy(100);
    TF3 * fDsSAngle = new TF3("fDsSAngle",decRefDsS,&decayRates::TF_dGdAngle,
                                  -1,1,-1,1,0.,2*TMath::Pi(),1);
   // fDsSAngle->SetNpx(100); fDsSAngle->SetNpy(100); fDsSAngle->SetNpz(100);
          
    cout << "fDs integral " << fDs->Integral(wmin,wmax) << endl;
    if(rateDS_4D)
        cout << "fDsS integral " << fDsS->Integral(wmin,wmaxS) << endl;
    else
        cout << "fDsS integral " << fDsS2D->Integral(wmin,wmaxS,-1,1) << endl;
    
    TH2D* hres2D[2];
    TH1D* hres[2];
    TH1D* hacc[2];
    TFile * fileData = TFile::Open("LHCb_input_data_mc.root");
    hres2D[0] = (TH2D*) fileData->Get("h_res2D_Ds");
    hres2D[1] = (TH2D*) fileData->Get("h_res2D_DsS");
    hres[0] = (TH1D*) fileData->Get("h_res_Ds");
    hres[1] = (TH1D*) fileData->Get("h_res_DsS");
    hacc[0] = (TH1D*) fileData->Get("h_sig_Ds");
    hacc[1] = (TH1D*) fileData->Get("h_sig_DsS");

    TH1D* hden[2];
    for(int i=0;i<2;++i) {
        hres2D[i]->SetDirectory(0);
        hres[i]->SetDirectory(0);
        hacc[i]->SetDirectory(0);
        hden[i] = (TH1D*) hacc[i]->Clone(Form("hden_%i",1));
        hden[i]->Reset();
        hden[i]->SetDirectory(0);
    }
    fileData->Close();
    
    int signalCategory;
    double w{-999}, ctl{-999}, ctd{-999}, chi{-999}, pperp_t{-999}, pperp{-999}, den{-999}, den1D{-999}, den4D{-999};
    TFile* file = new TFile(outfile.Data(),"RECREATE");
    TTree* tree = new TTree(outtree.Data(),outtree.Data());
    tree->Branch("signalCategory",&signalCategory);
    tree->Branch("w",&w);
    tree->Branch("ctl",&ctl);
    tree->Branch("ctd",&ctd);
    tree->Branch("chi",&chi);
    tree->Branch("pperp_t",&pperp_t);
    tree->Branch("pperp",&pperp);
    tree->Branch("den",&den);
    tree->Branch("den1D",&den1D);
    tree->Branch("den4D",&den4D);
    
    TRandom3 rnd(345);
    for(unsigned int i = 0; i < nentries; i++) {
   
       if(i%50000==0 && i!=0) cout << i << endl;
       //generate for BsDs
       signalCategory = 0;
       w=-999; ctl=-999; ctd=-999; chi=-999; pperp_t=-999; pperp=-999; den=-999, den1D=-999, den4D=-999;
       double cos = rnd.Uniform(0,1); //cos alpha extraction
       w        = fDs->GetRandom();
       pperp_t  = decRefDs->Eval_Pperp(w,cos,_Dsmass);
       den      = decRefDs->Eval_dGdw(w) / rateRefDs;
        
       if(_res2D){
           TH1D * htmp0 = (TH1D*) hres2D[0]->ProjectionX("hTmp0",
                                                         hres2D[0]->GetYaxis()->FindBin(pperp_t),
                                                         hres2D[0]->GetYaxis()->FindBin(pperp_t));
           pperp = pperp_t + htmp0->GetRandom();
           htmp0->Delete();
       }else
         pperp = pperp_t + hres[0]->GetRandom();
        
       hden[0]->Fill(pperp);
       tree->Fill();
       
        
       //generate for BsDsS
       signalCategory = 1;
       w=-999; ctl=-999; ctd=-999; chi=-999; pperp_t=-999; pperp=-999; den=-999, den1D=-999, den4D=-999;
       cos   = rnd.Uniform(0,1); //cos alpha extraction
        if(rateDS_4D){
            w = fDsS->GetRandom();
            fDsSAngle->SetParameter(0,w);
            fDsSAngle->GetRandom3(ctl,ctd,chi);
            den4D   = decRefDsS->Eval_dGdwdAngle(w,ctl,ctd,chi) / rateRefDs;
		    den1D   = decRefDsS->Eval_dGdw(w) / rateRefDsS;
        }else{
            fDsS2D->GetRandom2(w,ctd);
            den1D   = decRefDsS->Eval_dGdw(w) / rateRefDsS;
        }
       
       double beta = sqrt(w*w-1)/w;
       pperp_t = decRefDsS->Eval_Pperp(beta,ctd,cos,_DsSmass,_Dsmass,0);
       if(_res2D){
           TH1D * htmp1 = (TH1D*) hres2D[1]->ProjectionX("hTmp1",
                                                          hres2D[1]->GetYaxis()->FindBin(pperp_t),
                                                          hres2D[1]->GetYaxis()->FindBin(pperp_t));
           pperp = pperp_t + htmp1->GetRandom();
           htmp1->Delete();
       }else{
          pperp = pperp_t + hres[1]->GetRandom();
       }
        
       hden[1]->Fill(pperp);
       tree->Fill();
        
     }
    
    cout << "Calculate acceptance histograms, which has min and max values: " << endl;
    double min{-999}, max{-999};
     for(int i=0; i<2; ++i){
       hacc[i]->Scale(1./hacc[i]->Integral());
       hden[i]->Scale(1./hden[i]->Integral());
       hacc[i]->Divide(hden[i]);
       hacc[i]->GetMinimumAndMaximum(min,max);
       cout << " component: " << i << "\t " << min << "\t " << max << endl;
     }
       
     for(int i=0; i<2; ++i){
         hacc[i]->Write();
         hres[i]->Write();
         hres2D[i]->Write();
         tree->Write();
         fDs->Write();
         fDsS->Write();
         fDsS2D->Write();
         fDsSAngle->Write();
     }
     file->Close();
       
     return;

    
 return;
}


void runGeneration(TString configFile = "config/fit_testCLN.cfg")
{
  gROOT->ProcessLine(".L ./decayRates.C+");
  generateSignalTemplates(configFile);
}


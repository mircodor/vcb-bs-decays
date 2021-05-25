void RunFit(TString configFile = "config/fit_testCLN.cfg")
{

  gROOT->ProcessLine(".L ./decayRates.C+");
  gROOT->ProcessLine(".L ./fitter.C+");
  gROOT->ProcessLine(".x ./lhcbStyle.C");
  gROOT->ProcessLine("fitter fit;");
  gROOT->ProcessLine(Form("fit.RunFitter(\"%s\");",
			  configFile.Data()));

}

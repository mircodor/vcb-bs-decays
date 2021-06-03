void RunFit(TString configFile = "config/fit_testBGL.cfg")
{

  //gROOT->ProcessLine(".L ./lhcbStyle.C");
  gROOT->ProcessLine(".L ./decayRates.C+");
  gROOT->ProcessLine(".L ./fitter.C+");
  gROOT->ProcessLine("fitter fit;");
  gROOT->ProcessLine(Form("fit.RunFitter(\"%s\");",
			  configFile.Data()));

}

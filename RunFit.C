void RunFit(TString configFile = "config/fit_testBGL.cfg")
{

  //gROOT->ProcessLine(".L ./plotStyle.C");
  gROOT->ProcessLine(".L ./UpdatedDecayRate.C+");
  gROOT->ProcessLine(".L ./fitter.C+");
  gROOT->ProcessLine("fitter fit;");
  gROOT->ProcessLine(Form("fit.RunFitter(\"%s\");",
			  configFile.Data()));

}

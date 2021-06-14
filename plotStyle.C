#include "TROOT.h"
#include "TStyle.h"
#include "TText.h"
#include "TLatex.h"


void plotStyle(){

  Int_t black  = 1;
  Int_t red    = 2;
  Int_t green  = 3;
  Int_t blue   = 4;
  Int_t yellow = 5; 
  Int_t magenta= 6;
  Int_t cyan   = 7;
  Int_t purple = 9;
  


  Int_t plotFont        = 132;  
  Double_t plotWidth    = 2.00; 
  Double_t plotTSize    = 0.06; 
  

  gROOT->SetStyle("Plain"); 
  TStyle *plotStyle= new TStyle("plotStyle","Plots style");
  

  plotStyle->SetFillColor(1);
  plotStyle->SetFillStyle(1001);  
  plotStyle->SetFrameFillColor(0);
  plotStyle->SetFrameBorderMode(0);
  plotStyle->SetPadBorderMode(0);
  plotStyle->SetPadColor(0);
  plotStyle->SetCanvasBorderMode(0);
  plotStyle->SetCanvasColor(0);
  plotStyle->SetStatColor(0);
  plotStyle->SetLegendBorderSize(0);

  plotStyle->SetPalette(1);

  int colors[8] = {0,5,7,3,6,2,4,1};
  plotStyle->SetPalette(8,colors);


  plotStyle->SetPaperSize(20,26);
  plotStyle->SetPadTopMargin(0.05);
  plotStyle->SetPadRightMargin(0.05); 
  plotStyle->SetPadBottomMargin(0.16);
  plotStyle->SetPadLeftMargin(0.14);
  

  plotStyle->SetTextFont(plotFont);
  plotStyle->SetTextSize(plotTSize);
  plotStyle->SetLabelFont(plotFont,"x");
  plotStyle->SetLabelFont(plotFont,"y");
  plotStyle->SetLabelFont(plotFont,"z");
  plotStyle->SetLabelSize(plotTSize,"x");
  plotStyle->SetLabelSize(plotTSize,"y");
  plotStyle->SetLabelSize(plotTSize,"z");
  plotStyle->SetTitleFont(plotFont);
  plotStyle->SetTitleFont(plotFont,"x");
  plotStyle->SetTitleFont(plotFont,"y");
  plotStyle->SetTitleFont(plotFont,"z");
  plotStyle->SetTitleSize(1.2*plotTSize,"x");
  plotStyle->SetTitleSize(1.2*plotTSize,"y");
  plotStyle->SetTitleSize(1.2*plotTSize,"z");

  plotStyle->SetLineWidth(plotWidth);
  plotStyle->SetFrameLineWidth(plotWidth);
  plotStyle->SetHistLineWidth(plotWidth);
  plotStyle->SetFuncWidth(plotWidth);
  plotStyle->SetGridWidth(plotWidth);
  plotStyle->SetLineStyleString(2,"[12 12]"); 
  plotStyle->SetMarkerStyle(20);
  plotStyle->SetMarkerSize(1.0);

  plotStyle->SetLabelOffset(0.010,"X");
  plotStyle->SetLabelOffset(0.010,"Y");

  plotStyle->SetOptStat(0);  
  plotStyle->SetStatFormat("6.3g"); 
  plotStyle->SetOptTitle(0);
  plotStyle->SetOptFit(0);

  plotStyle->SetTitleOffset(1.0,"X");
  plotStyle->SetTitleOffset(1.2,"Y");
  plotStyle->SetTitleOffset(1.2,"Z");
  plotStyle->SetTitleFillColor(0);
  plotStyle->SetTitleStyle(0);
  plotStyle->SetTitleBorderSize(0);
  plotStyle->SetTitleFont(plotFont,"title");
  plotStyle->SetTitleX(0.0);
  plotStyle->SetTitleY(1.0); 
  plotStyle->SetTitleW(1.0);
  plotStyle->SetTitleH(0.05);
  

  plotStyle->SetStatBorderSize(0);
  plotStyle->SetStatFont(plotFont);
  plotStyle->SetStatFontSize(0.05);
  plotStyle->SetStatX(0.9);
  plotStyle->SetStatY(0.9);
  plotStyle->SetStatW(0.25);
  plotStyle->SetStatH(0.15);


  plotStyle->SetPadTickX(0);
  plotStyle->SetPadTickY(0);


  plotStyle->SetNdivisions(505,"x");
  plotStyle->SetNdivisions(510,"y");
  
  gROOT->SetStyle("plotStyle");
  gROOT->ForceStyle();

  TText *plotLabel = new TText();
  plotLabel->SetTextFont(plotFont);
  plotLabel->SetTextColor(1);
  plotLabel->SetTextSize(plotTSize);
  plotLabel->SetTextAlign(12);

  TLatex *plotLatex = new TLatex();
  plotLatex->SetTextFont(plotFont);
  plotLatex->SetTextColor(1);
  plotLatex->SetTextSize(plotTSize);
  plotLatex->SetTextAlign(12);

  cout << "-------------------------" << endl;  
  cout << "Set plot Style" << endl;
  cout << "-------------------------" << endl;  
  
}

                 

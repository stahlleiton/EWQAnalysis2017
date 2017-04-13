#include "GetAccEff.h"
#include "HiMuonTree.h"
#include "HiMETTree.h"



void DrawAccEff_DYZ (
    const string muIDType="Tight",
    const string beamDir="Pbp"
    ) {

  string inputFile=Form("AccEff_%s_%s.root",muIDType.c_str(),beamDir.c_str());
  TFile finput(inputFile.c_str());

  // Define efficiency histograms
  // CombineAccEff == true:  acc*eff
  // num[0]: normal
  // num[1]: veto eff
  TH1D *hden_eta, *hnum_eta[2];
  TGraphAsymmErrors *heff_eta[2];

  // Load efficiency histograms
  string suffix_ = beamDir + "_" + muIDType;
  string suffix;
  hden_eta = static_cast<TH1D*>(finput.Get(Form("hden_eta_%s",suffix_.c_str())));
  for (int i=0; i<2; i++) {
    if (i==0) suffix = suffix_ + "_normal";
    if (i==1) suffix = suffix_ + "_veto";
    hnum_eta[i] = static_cast<TH1D*>(finput.Get(Form("hnum_eta_%s",suffix.c_str())));
    heff_eta[i] = static_cast<TGraphAsymmErrors*>(finput.Get(Form("heff_eta_%s",suffix.c_str())));
    heff_eta[i]->GetXaxis()->SetTitle("#eta_{lab}");
    heff_eta[i]->GetYaxis()->SetTitle("Efficiency");
  }


  const int color[] = {kBlack, kRed+1, kOrange+2, kSpring+1, kGreen+1, kAzure+7, kBlue+1, kViolet, kViolet+3, kPink};

  for (unsigned int i=0; i<2; i++) {
    SetGraphStyle(heff_eta[i],color[i],kOpenCircle,1,-2.5,2.5,0,1.2);
  }

  // Draw graphs
  TLegend leg(0.67,0.85,0.95,0.95);
  leg.SetFillColor(0);
  leg.SetFillStyle(4000);
  leg.SetBorderSize(0);
  leg.SetMargin(0.15);

  leg.AddEntry(heff_eta[0],"DY, Z vetoed","lp");
  leg.AddEntry(heff_eta[1],"DY, Z left","lp");
  
  TCanvas canv("canv","canv",600,600);
  canv.SetTopMargin(0.04);
  canv.SetBottomMargin(0.1);
  canv.SetRightMargin(0.04);
  canv.SetLeftMargin(0.13);

  heff_eta[0]->Draw("ap");
  heff_eta[1]->Draw("p");
  heff_eta[0]->GetYaxis()->SetTitleOffset(1.3);
  canv.Update();
  leg.Draw();
  canv.SaveAs(Form("heff_eta_%s_%s.png",muIDType.c_str(),beamDir.c_str()));
  canv.SaveAs(Form("heff_eta_%s_%s.pdf",muIDType.c_str(),beamDir.c_str()));
  canv.Clear();


}



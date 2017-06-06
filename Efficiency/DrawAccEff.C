#include "GetAccEff.h"
#include "HiMuonTree.h"
#include "HiMETTree.h"



void DrawAccEff (
    const string beamDir="Pbp",
    const int CombineAccEff=0
    ) {

  string inputFile=Form("AccEff_%s.root",beamDir.c_str());
  TFile finput(inputFile.c_str());

  // Define efficiency histograms
  //[0]: mu+
  //[1]: mu-
  // CombineAccEff == 0: acc*eff
  // CombineAccEff == 1: eff, 2: 5 TeV W analysis setting eff
  TH1D *hden_eta[2], *hnum_eta[2];
  TH1D *hden_pt[2][nbins_eta], *hnum_pt[2][nbins_eta];
  TGraphAsymmErrors *heff_eta[2];
  TGraphAsymmErrors *heff_pt[2][nbins_eta];

  // CombineAccEff == 1: acc, 2: 5 TeV W analysis setting acc
  TH1D *hden_acc_eta[2], *hnum_acc_eta[2];
  TH1D *hden_acc_pt[2][nbins_eta], *hnum_acc_pt[2][nbins_eta];
  TGraphAsymmErrors *heff_acc_eta[2];
  TGraphAsymmErrors *heff_acc_pt[2][nbins_eta];

  // Load efficiency histograms
  string suffix = beamDir;
  for (int i=0; i<2; i++) {
    if (i==0) suffix = beamDir + "_muPlus";
    if (i==1) suffix = beamDir + "_muMinus";
    hden_eta[i] = static_cast<TH1D*>(finput.Get(Form("hden_eta_%s",suffix.c_str())));
    hnum_eta[i] = static_cast<TH1D*>(finput.Get(Form("hnum_eta_%s",suffix.c_str())));
    heff_eta[i] = static_cast<TGraphAsymmErrors*>(finput.Get(Form("heff_eta_%s",suffix.c_str())));
    heff_eta[i]->GetXaxis()->SetTitle("#eta_{lab}");
    heff_eta[i]->GetYaxis()->SetTitle("Efficiency");
    if (CombineAccEff>0) {
      hden_acc_eta[i] = static_cast<TH1D*>(finput.Get(Form("hden_acc_eta_%s",suffix.c_str())));
      hnum_acc_eta[i] = static_cast<TH1D*>(finput.Get(Form("hnum_acc_eta_%s",suffix.c_str())));
      heff_acc_eta[i] = static_cast<TGraphAsymmErrors*>(finput.Get(Form("heff_acc_eta_%s",suffix.c_str())));
      heff_acc_eta[i]->GetXaxis()->SetTitle("#eta_{lab}");
      heff_acc_eta[i]->GetYaxis()->SetTitle("Acceptance");
    }
    for (int j=0; j<nbins_eta; j++) {
      hden_pt[i][j] = static_cast<TH1D*>(finput.Get(Form("hden_pt_%.0f-%.0f_%s",bins_eta[j]*10,bins_eta[j+1]*10,suffix.c_str())));
      hnum_pt[i][j] = static_cast<TH1D*>(finput.Get(Form("hnum_pt_%.0f-%.0f_%s",bins_eta[j]*10,bins_eta[j+1]*10,suffix.c_str())));
      heff_pt[i][j] = static_cast<TGraphAsymmErrors*>(finput.Get(Form("heff_pt_%.0f-%.0f_%s",bins_eta[j]*10,bins_eta[j+1]*10,suffix.c_str())));
      heff_pt[i][j]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
      heff_pt[i][j]->GetYaxis()->SetTitle("Efficiency");
      if (CombineAccEff>0) {
        hden_acc_pt[i][j] = static_cast<TH1D*>(finput.Get(Form("hden_acc_pt_%.0f-%.0f_%s",bins_eta[j]*10,bins_eta[j+1]*10,suffix.c_str())));
        hnum_acc_pt[i][j] = static_cast<TH1D*>(finput.Get(Form("hnum_acc_pt_%.0f-%.0f_%s",bins_eta[j]*10,bins_eta[j+1]*10,suffix.c_str())));
        heff_acc_pt[i][j] = static_cast<TGraphAsymmErrors*>(finput.Get(Form("heff_acc_pt_%.0f-%.0f_%s",bins_eta[j]*10,bins_eta[j+1]*10,suffix.c_str())));
        heff_acc_pt[i][j]->GetXaxis()->SetTitle("p_{T}");
        heff_acc_pt[i][j]->GetYaxis()->SetTitle("Acceptance");
      }
    }
  }

  // Calculate efficiency
/*  for (int i=0; i<2; i++) {
    CheckUnderFlow(hnum_eta[i],hden_eta[i]);
    heff_eta[i]->Divide(hnum_eta[i],hden_eta[i]);
    if (CombineAccEff>0) {
      CheckUnderFlow(hnum_acc_eta[i],hden_acc_eta[i]);
      heff_acc_eta[i]->Divide(hnum_acc_eta[i],hden_acc_eta[i]);
    }
    for (int j=0; j<nbins_eta; j++) {
      CheckUnderFlow(hnum_pt[i][j],hden_pt[i][j]);
      heff_pt[i][j]->Divide(hnum_pt[i][j],hden_pt[i][j]);
      if (CombineAccEff>0) {
        CheckUnderFlow(hnum_acc_pt[i][j],hden_acc_pt[i][j]);
        heff_acc_pt[i][j]->Divide(hnum_acc_pt[i][j],hden_acc_pt[i][j]);
      }
    }
  }
*/

  const int color[] = {kBlack, kRed+1, kOrange+2, kSpring+1, kGreen+1, kAzure+7, kBlue+1, kViolet, kViolet+3, kPink};

  for (unsigned int i=0; i<2; i++) {
    SetGraphStyle(heff_eta[i],color[i],kOpenCircle,1,-2.5,2.5,0,1.2);
    if (CombineAccEff>0) {
      SetGraphStyle(heff_acc_eta[i],color[i],kOpenCircle,1,-2.5,2.5,0,1.2);
    }
    for (unsigned int j=0; j<nbins_eta; j++) {
      SetGraphStyle(heff_pt[i][j],color[i],kOpenCircle,1,0,0,0,1.2);
      if (CombineAccEff>0) {
        SetGraphStyle(heff_acc_pt[i][j],color[i],kOpenCircle,1,0,0,0,1.2);
      }
    }
  }

  // Draw graphs
  TLegend leg(0.67,0.85,0.95,0.95);
  leg.SetFillColor(0);
  leg.SetFillStyle(4000);
  leg.SetBorderSize(0);
  leg.SetMargin(0.15);

  leg.AddEntry(heff_eta[0],"W+","lp");
  leg.AddEntry(heff_eta[1],"W-","lp");
  
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
  canv.SaveAs(Form("heff_eta_%s.png",beamDir.c_str()));
  canv.SaveAs(Form("heff_eta_%s.pdf",beamDir.c_str()));
  canv.Clear();

  if (CombineAccEff>0) {
    heff_acc_eta[0]->Draw("ap");
    heff_acc_eta[1]->Draw("p");
    heff_acc_eta[0]->GetYaxis()->SetTitleOffset(1.3);
    canv.Update();
    leg.Draw();
    canv.SaveAs(Form("heff_acc_eta_%s.png",beamDir.c_str()));
    canv.SaveAs(Form("heff_acc_eta_%s.pdf",beamDir.c_str()));
    canv.Clear();
  }

/*  TLegend leg_pt(0.67,0.35,0.95,0.45);
  leg_pt.SetFillColor(0);
  leg_pt.SetFillStyle(4000);
  leg_pt.SetBorderSize(0);
  leg_pt.SetMargin(0.15);

  heff_pt[0][0]->Draw("ap");
  for (int j=0; j<nbins_eta; j++) {
    for (int i=0; i<2; i++) {
      leg_pt.AddEntry(heff_pt[i][j],Form("%s, #eta: [%.1f, %.1f]",(i==0?"W+":"W-"),bins_eta[j],bins_eta[j+1]),"lp");
      if (i==0) heff_pt[i][j]->Draw("ap");
      else heff_pt[i][j]->Draw("p");
    }
    leg_pt.Draw();
    canv.SaveAs(Form("heff_pt_eta_%.0f-%.0f_%s_%s.png",bins_eta[j]*10,bins_eta[j+1]*10,muIDType.c_str(),beamDir.c_str()));
    canv.SaveAs(Form("heff_pt_eta_%.0f-%.0f_%s_%s.pdf",bins_eta[j]*10,bins_eta[j+1]*10,muIDType.c_str(),beamDir.c_str()));
    canv.Clear();
    leg_pt.Clear();
  }

  if (CombineAccEff>0) {
    heff_acc_pt[0][0]->Draw("ap");
    for (int j=0; j<nbins_eta; j++) {
      for (int i=0; i<2; i++) {
        leg_pt.AddEntry(heff_pt[i][j],Form("%s, #eta: [%.1f, %.1f]",(i==0?"W+":"W-"),bins_eta[j],bins_eta[j+1]),"lp");
        if (i==0) heff_acc_pt[i][j]->Draw("ap");
        else heff_acc_pt[i][j]->Draw("p");
      }
      leg_pt.Draw();
      canv.SaveAs(Form("heff_acc_pt_eta_%.0f-%.0f_%s_%s.png",bins_eta[j]*10,bins_eta[j+1]*10,muIDType.c_str(),beamDir.c_str()));
      canv.SaveAs(Form("heff_acc_pt_eta_%.0f-%.0f_%s_%s.pdf",bins_eta[j]*10,bins_eta[j+1]*10,muIDType.c_str(),beamDir.c_str()));
      canv.Clear();
      leg_pt.Clear();
    }

  }
*/

}



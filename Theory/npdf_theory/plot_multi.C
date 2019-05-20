// to compile:
//   g++ plot_multi.C `lhapdf-config --cflags --libs` `root-config --cflags --libs` -lAfterImage -o plot_multi
//
// to run:
//   ./plot_multi CT14nlo,EPPS16nlo_CT14nlo_Pb208

#include "rebin.C"
#include "lhapdf_utils.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h"
#include "TH1F.h"
#include "TString.h"
#include "TROOT.h"

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <cstdlib>

#include "tdrstyle.C"
#include "CMS_lumi.C"

#define NCT10 53
#define NEPS09 30
#define NCT14 57
#define NEPPS16 97
#define NnCTEQ15 33

void plot_graphs(vector<TGraphAsymmErrors*> graphs, vector<string> names, vector<Color_t> colors, vector<int> fillstyles, double ymin, double ymax, const char* xtitle, const char* ytitle, const char* cname);
void plot(TString pdfnames, double lumi);

void plot_multi(const std::string type="EPPS16nlo_CT14nlo_Pb208", const double lumi=-1)
{
  plot(type.c_str(), lumi);
}

void plot(TString pdfnames, double lumi)
{
  map<string,bool> Use;
  Use["CT10nlo"]=false;
  Use["CT14nlo"]=false;
  Use["EPS09nlo"]=false;
  Use["EPPS16nlo_CT14nlo_Pb208"]=false;
  Use["nCTEQ15_208_82"]=false;

  TString tok;
  Ssiz_t from = 0;
  while (pdfnames.Tokenize(tok, from, ",")) {
    // Analyse tok
    if (tok=="CT14nlo") {
      Use["CT14nlo"] = true;
      cout << "using CT14nlo" << endl;

    }
    if (tok=="CT10nlo") {
      Use["CT10nlo"] = true;
      cout << "using CT10nlo" << endl;

    }
    if (tok=="EPS09nlo") {
      Use["EPS09nlo"] = true;
      cout << "using EPS09nlo" << endl;

    }
    if (tok=="nCTEQ15_208_82") {
      Use["nCTEQ15_208_82"] = true;
      cout << "using nCTEQ15nlo" << endl;

    }
    if (tok=="EPPS16nlo_CT14nlo_Pb208") {
      Use["EPPS16nlo_CT14nlo_Pb208"] = true;
      cout << "using EPPS16nlo_CT14nlo_Pb208" << endl;

    }
  }

  map<string,vector<string> > wpname, wmname;
  for (int i=0; i<NCT14; i++) {
    wpname["CT14nlo"].push_back(Form("MCFM/CT14nlo/W_only_nlo_CT14nlo_80___80___W1_nlo_CT14nlo_%d.root",i));
    wmname["CT14nlo"].push_back(Form("MCFM/CT14nlo/W_only_nlo_CT14nlo_80___80___W6_nlo_CT14nlo_%d.root",i));
  }

  for (int i=0; i<NCT10; i++) {
    wpname["CT10nlo"].push_back(Form("MCFM/CT10nlo/W_only_nlo_CT10nlo_80___80___W1_nlo_CT10nlo_%d_1.root",i));
    wmname["CT10nlo"].push_back(Form("MCFM/CT10nlo/W_only_nlo_CT10nlo_80___80___W6_nlo_CT10nlo_%d_1.root",i));
  }

  for (int i=0; i<NnCTEQ15; i++) {
    wpname["nCTEQ15_208_82"].push_back(Form("MCFM/nCTEQ15nlo/W_only_nlo_CT14nlo_80___80___W1_nlo_nCTEQ15nlo_0_%d.root",i));
    wmname["nCTEQ15_208_82"].push_back(Form("MCFM/nCTEQ15nlo/W_only_nlo_CT14nlo_80___80___W6_nlo_nCTEQ15nlo_0_%d.root",i));
  }
  for (uint i = 1; i < NCT14; i++) {
    wpname["nCTEQ15_208_82"].push_back(Form("MCFM/nCTEQ15nlo/W_only_nlo_CT14nlo_80___80___W1_nlo_nCTEQ15nlo_%d_0.root",i));
    wmname["nCTEQ15_208_82"].push_back(Form("MCFM/nCTEQ15nlo/W_only_nlo_CT14nlo_80___80___W6_nlo_nCTEQ15nlo_%d_0.root",i));
  }

  for (int i=0; i<NEPPS16; i++) {
    int i1 = (i<=40) ? 0 : i-40;
    wpname["EPPS16nlo_CT14nlo_Pb208"].push_back(Form("MCFM/EPPS16nlo/W_only_nlo_CT14nlo_80___80___W1_nlo_EPPS16nlo_%d_%d.root",i1,i));
    wmname["EPPS16nlo_CT14nlo_Pb208"].push_back(Form("MCFM/EPPS16nlo/W_only_nlo_CT14nlo_80___80___W6_nlo_EPPS16nlo_%d_%d.root",i1,i));
  }

  for (int i=0; i<NCT10; i++) {
    wpname["EPS09nlo"].push_back(Form("MCFM/EPS09nlo/W_only_nlo_CT10nlo_80___80___W1_nlo_EPS09_%d_%d.root",i,1));
    wmname["EPS09nlo"].push_back(Form("MCFM/EPS09nlo/W_only_nlo_CT10nlo_80___80___W6_nlo_EPS09_%d_%d.root",i,1));
  }
  for (int i=2; i<NEPS09+2; i++) {
    wpname["EPS09nlo"].push_back(Form("MCFM/EPS09nlo/W_only_nlo_CT10nlo_80___80___W1_nlo_EPS09_%d_%d.root",0,i));
    wmname["EPS09nlo"].push_back(Form("MCFM/EPS09nlo/W_only_nlo_CT10nlo_80___80___W6_nlo_EPS09_%d_%d.root",0,i));
  }

  map<string,Color_t> allcolors;
  allcolors["CT10nlo"] = kYellow+2;
  allcolors["EPS09nlo"] = kGreen+2;
  allcolors["CT14nlo"] = kRed;
  allcolors["nCTEQ15_208_82"] = kMagenta;
  allcolors["EPPS16nlo_CT14nlo_Pb208"] = kBlue;

  map<string,int> allstyles;
  allstyles["CT10nlo"] = 3003;
  allstyles["CT14nlo"] = 3004;
  allstyles["EPS09nlo"] = 3144;
  allstyles["nCTEQ15_208_82"] = 3007;
  allstyles["EPPS16nlo_CT14nlo_Pb208"] = 3005;

  std::vector< TGraphAsymmErrors* > gps, gms, gchs, ga1ps, ga1ms, ga3s; std::vector<string> names; std::vector<Color_t> colors; std::vector<int> fillstyles;
  std::map<std::string,bool>::iterator it;
  for (it = Use.begin(); it != Use.end(); it++)
    {
      std::vector< TH1F* > hps, hms, hchs, ha1ps, ha1ms, ha3s;
      if (!it->second) continue;
      string tag = it->first;

      // read and make the histos for all PDF members
      for (uint j=0; j<wpname[tag].size(); j++) {
        const char* fwp = wpname[tag][j].c_str();
        const char* fwm = wmname[tag][j].c_str();
        hps.push_back(rebin(fwp));
        hms.push_back(rebin(fwm));
        hchs.push_back(chasym(fwp,fwm));
        ha1ps.push_back(a1plus(fwp));
        ha1ms.push_back(a1plus(fwm));
        ha3s.push_back(a3(fwp,fwm));
      }

      // now produce the final graphs with uncertainties
      gps.push_back(getPDFHessianError(hps));
      gms.push_back(getPDFHessianError(hms));
      gchs.push_back(getPDFHessianError(hchs));
      ga1ps.push_back(getPDFHessianError(ha1ps));
      ga1ms.push_back(getPDFHessianError(ha1ms));
      ga3s.push_back(getPDFHessianError(ha3s));
      
      names.push_back(tag);
      colors.push_back(allcolors[tag]);
      fillstyles.push_back(allstyles[tag]);
    }

  if (lumi>0) { // add pseudo-data projection
    TH1F *hplus = rebin(wpname["EPPS16nlo_CT14nlo_Pb208"][0].c_str());
    hplus->Scale(lumi);
    for (int i=1; i<=hplus->GetNbinsX(); i++) hplus->SetBinError(i,sqrt(hplus->GetBinContent(i)));
      
    TH1F *hminus = rebin(wmname["EPPS16nlo_CT14nlo_Pb208"][0].c_str());
    hminus->Scale(lumi);
    for (int i=1; i<=hminus->GetNbinsX(); i++) hminus->SetBinError(i,sqrt(hminus->GetBinContent(i)));

    // TODO FIXME add systs
    double syst=0.01; // 0.02
    gps.push_back(hist2graph(hplus));
    gms.push_back(hist2graph(hminus));
    gchs.push_back(hist2graph(chasym(hplus,hminus),syst));
    ga1ps.push_back(hist2graph(a1plus(hplus),syst));
    ga1ms.push_back(hist2graph(a1plus(hminus),syst));
    ga3s.push_back(hist2graph(a3(hplus,hminus),syst));

    // names.push_back(Form("Projection (%.0f nb^{-1})",lumi));
    names.push_back("Data");
    colors.push_back(kBlack);
    fillstyles.push_back(0);
  }

  plot_graphs(gps, names, colors, fillstyles, 0, 200, "#eta_{cm}", "d#sigma(W^{+} #rightarrow l^{+}#nu / d#eta_{cm} [nb]","Wp");
  plot_graphs(gms, names, colors, fillstyles, 0, 200, "#eta_{cm}", "d#sigma(W^{-} #rightarrow l^{-}#nu / d#eta_{cm} [nb]","Wm");
  plot_graphs(gchs, names, colors, fillstyles, -0.2, 0.3, "#eta_{cm}", "(N^{+} - N^{-}) / (N^{+} + N^{-})","chasym");
  plot_graphs(ga1ps, names, colors, fillstyles, 0.8, 1.6, "#eta_{cm}", "N^{+} (+#eta_{cm}) / N^{+} (-#eta_{cm})","A1p");
  plot_graphs(ga1ms, names, colors, fillstyles, 0.6, 1.1, "#eta_{cm}", "N^{-} (+#eta_{cm}) / N^{-} (-#eta_{cm})","A1m");
  plot_graphs(ga3s, names, colors, fillstyles, 0.7, 1.25, "#eta_{cm}", "N (+#eta_{cm}) / N (-#eta_{cm})","A3");

}

void plot_graphs(vector<TGraphAsymmErrors*> graphs, vector<string> names, vector<Color_t> colors, vector<int> fillstyles, double ymin, double ymax, const char* xtitle, const char* ytitle, const char* cname)
{
  setTDRStyle();
  gROOT->SetStyle( "tdrStyle" );

  TCanvas *cch = new TCanvas(cname,"",600,600);
  cch->cd();
   
  // checks
  unsigned int nh = graphs.size();
  if (names.size() != nh) return;
  if (colors.size() != nh) return;
  if (fillstyles.size() != nh) return;

  for (uint i=0; i<nh; i++)
    {
      graphs[i]->SetFillColor(colors[i]);
      graphs[i]->SetFillStyle(fillstyles[i]);
      graphs[i]->SetLineColor(colors[i]);
      graphs[i]->SetMarkerColor(colors[i]);

      if (i==0) 
        {
          graphs[i]->GetYaxis()->SetRangeUser(ymin,ymax);
          if (!TString(cname).Contains("A")) graphs[i]->GetXaxis()->SetRangeUser(bins[0],bins[nbins]);
          else graphs[i]->GetXaxis()->SetRangeUser(bins2[0],bins[nbins2]);
          graphs[i]->GetYaxis()->SetTitle(ytitle);
          graphs[i]->GetXaxis()->SetTitle(xtitle);
          graphs[i]->Draw("a5");
        }
      else if (!TString(names[i]).Contains("Data")) graphs[i]->Draw("5");
      else graphs[i]->Draw("P");
    }

  double x1=0.55, y1=0.17, x2=0.88, y2=0.4;
  if (TString(cname)=="A3") {
    x1=0.55;y1=0.68;x2=0.88;y2=0.91;
  }
  if (TString(cname)=="A1p") {
    x1=0.2;y1=0.63;x2=0.55;y2=0.86;
  }
  if (TString(cname)=="A1m") {
    x1=0.2;y1=0.17;x2=0.55;y2=0.4;
  }

  TLegend *tlegch = new TLegend(x1,y1,x2,y2);
  tlegch->SetHeader("MCFM nlo");
  tlegch->SetBorderSize(0);
  tlegch->SetTextSize(0.04);
  for (uint i=0; i<nh; i++) {
    TString label(names[i]);
    if (label=="EPPS16nlo_CT14nlo_Pb208") label = "EPPS16nlo";
    if (label.Contains("Data")) tlegch->AddEntry(graphs[i],label.Data(),"lp");
    else tlegch->AddEntry(graphs[i],label.Data(),"lpf");
  }
  tlegch->Draw();

  extraText = "Projection";
  writeExtraText = true;
  CMS_lumi( cch, 111);//, 0 );

  cch->Print(Form("Plot/%s.pdf",cname));
  cch->SaveAs(Form("Plot/%s.root",cname));
  cch->SaveAs(Form("Plot/%s.C",cname));
}

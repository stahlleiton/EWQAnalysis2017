#include "Macros/CMS/CMS_lumi.C"
#include "Macros/CMS/tdrstyle.C"
#include "Macros/Utilities/bin.h"
#include "Macros/Utilities/initClasses.h"
#include "Macros/Utilities/resultsUtils.h"
#include "Macros/Utilities/texUtils.h"

#include <vector>
#include <map>
#include <string>
#include <sstream>
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TLegend.h"
#include "TLine.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TArrow.h"

////////////////
// PARAMETERS //
////////////////

#ifndef poiname_check
#define poiname_check
const char* poiname = "N_WToMu"; // for W Asymmetry (will correct automatically for efficiency)
#endif
string ylabel  = "";
string nameTag_base = "";    // 
bool addTheory = true;

typedef std::map< std::string, std::map< anabin<0>, TGraphAsymmErrors* >, std::greater<std::string> >         TGraphBinDict;
typedef std::map< std::string, std::map< anabin<0> , std::vector<double> > >       DoubleVecBinDict;
typedef std::map< std::string, std::map< anabin<0> , std::vector< anabin<0> > > >  BinVecBinDict;

//////////////////
// DECLARATIONS //
//////////////////

void printOptions();
void plotGraph(const TGraphBinDict& theGraphs, const std::string& xaxis, const std::string& outputDir,  const std::string type, 
               const std::vector<std::string>& infoLabel=std::vector<std::string>(), const bool& merge=false);
void plot(std::vector< anabin<0> > thecats, string xaxis, string workDirName);
void flipXAxis(TGraphAsymmErrors& g);
int color(int i);
int markerstyle(int i);
string nameTag;

class asymm_input {
public:
  double nPl;
  double nMi;
  double dnPl_stat;
  double dnMi_stat;
};


/////////////////////////////////////////////
// MAIN FUNCTIONS TO BE CALLED BY THE USER //
/////////////////////////////////////////////
/*
     theCats.push_back(anabin<0>(-2.4,-2.0,0.0,0.15,0,200));
     theCats.push_back(anabin<0>(-2.0,-1.5,0.0,0.15,0,200));
     theCats.push_back(anabin<0>(-1.5,-1.0,0.0,0.15,0,200));
     theCats.push_back(anabin<0>(-1.0,-0.5,0.0,0.15,0,200));
     theCats.push_back(anabin<0>(-0.5,+0.0,0.0,0.15,0,200));
     theCats.push_back(anabin<0>(+0.0,+0.5,0.0,0.15,0,200));
     theCats.push_back(anabin<0>(+0.5,+1.0,0.0,0.15,0,200));
     theCats.push_back(anabin<0>(+1.0,+1.5,0.0,0.15,0,200));
     theCats.push_back(anabin<0>(+1.5,+2.0,0.0,0.15,0,200));
     theCats.push_back(anabin<0>(+2.0,+2.4,0.0,0.15,0,200));
*/

void plotEta(string workDirName)
{
  std::string xaxis = "eta";
  std::vector<anabin<0>> theCats;

  theCats.push_back(anabin<0>(-2.4,2.4,0.0,0.15,0,200));

  nameTag = nameTag_base;
  plot(theCats,xaxis,workDirName);
};

void plotAll(string workDirName)
{
  plotEta(workDirName);
};

/////////////////////
// OTHER FUNCTIONS //
/////////////////////

void plot(std::vector< anabin<0> > thecats, std::string xaxis, std::string outputDir)
{ 
  TFile *f = new TFile(treeFileName(outputDir.c_str(),"DATA","","MET"));
  if (!f || !f->IsOpen()) {
    resultsEWQ2tree(outputDir.c_str(),"DATA");
    f = new TFile(treeFileName(outputDir.c_str(),"DATA","","MET"));
    if (!f) return;
  }
  TTree *tr = (TTree*) f->Get("fitresults");
  if (!tr) { f->Close(); delete f; return; }

  std::map<std::string, std::map< anabin<0>, asymm_input> > theVars_inputs;

  std::vector<double> x, ex, y, ey;
  std::map< std::string , std::map < std::string , float > > evtVar = 
    {
      { "MET",        {{"Min" , -99.}, {"Max" , -99.}, {"Val" , -99.}, {"Err" , -99.}} },
      { "Muon_Pt",    {{"Min" , -99.}, {"Max" , -99.}, {"Val" , -99.}, {"Err" , -99.}} },
      { "Muon_Eta",   {{"Min" , -99.}, {"Max" , -99.}, {"Val" , -99.}, {"Err" , -99.}} },
      { "Muon_Iso",   {{"Min" , -99.}, {"Max" , -99.}, {"Val" , -99.}, {"Err" , -99.}} },
      { "Muon_MT",    {{"Min" , -99.}, {"Max" , -99.}, {"Val" , -99.}, {"Err" , -99.}} },
      { "Centrality", {{"Min" , -99.}, {"Max" , -99.}, {"Val" , -99.}, {"Err" , -99.}} }
    };
  float val, err=0;
  char collSystem[5], charge[5], cutLabel[200];

  // Event Variables
  for (auto v : evtVar) {
    for (auto p : v.second) {
      tr->SetBranchAddress(Form("%s_%s", v.first.c_str(), p.first.c_str()) ,&(evtVar.at(v.first).at(p.first)));
    }
  }
  tr->SetBranchAddress(Form("%s_Val",poiname),&val);
  tr->SetBranchAddress(Form("%s_Err",poiname),&err);
  tr->SetBranchAddress("collSystem",collSystem);
  tr->SetBranchAddress("charge",charge);
  tr->SetBranchAddress("CutAndCount_WToMu",cutLabel);
  const char* token = Form("%s_%s", charge, collSystem);

  int ntr = tr->GetEntries();
  for (int i=0; i<ntr; i++) {
    tr->GetEntry(i);
    anabin<0> thebin(evtVar.at("Muon_Eta").at("Min"), evtVar.at("Muon_Eta").at("Max"), evtVar.at("Muon_Iso").at("Min"), evtVar.at("Muon_Iso").at("Max"), 0, 200);
    if (thebin==anabin<0>(-99.,-99.,-99.,-99.,-99.,-99.)) { std::cout << "[ERROR] The bin was not set properly!" << std::endl; return; }
    if ((TString(collSystem)=="pPb") && (TString(charge)=="Pl")) {
      theVars_inputs["pPb"][thebin].nPl = val;
      theVars_inputs["pPb"][thebin].dnPl_stat = err;
    } 
    else if ((TString(collSystem)=="pPb") && (TString(charge)=="Mi")) {
      theVars_inputs["pPb"][thebin].nMi = val;
      theVars_inputs["pPb"][thebin].dnMi_stat = err;
    }
    else if ((TString(collSystem)=="Pbp") && (TString(charge)=="Pl")) {
      theVars_inputs["Pbp"][thebin].nPl = val;
      theVars_inputs["Pbp"][thebin].dnPl_stat = err;
    }
    else if ((TString(collSystem)=="Pbp") && (TString(charge)=="Mi")) {
      theVars_inputs["Pbp"][thebin].nMi = val;
      theVars_inputs["Pbp"][thebin].dnMi_stat = err;
    }
  }
  BinVecBinDict    theBins;
  DoubleVecBinDict theVarsBinned_Asymm;
  DoubleVecBinDict theVarsBinned_Asymm_stat;
  TGraphBinDict    theGraphs_Asymm;
  DoubleVecBinDict theVarsBinned_A1Pl;
  DoubleVecBinDict theVarsBinned_A1Pl_stat;
  TGraphBinDict    theGraphs_A1Pl;
  DoubleVecBinDict theVarsBinned_A1Mi;
  DoubleVecBinDict theVarsBinned_A1Mi_stat;
  TGraphBinDict    theGraphs_A1Mi;
  DoubleVecBinDict theVarsBinned_A3;
  DoubleVecBinDict theVarsBinned_A3_stat;
  TGraphBinDict    theGraphs_A3;
  // initialize the maps
  for (auto colM : theVars_inputs) {
    std::string col = colM.first;
    for (auto it : thecats) {
      theBins[col][it] = std::vector< anabin<0> >();
      theVarsBinned_Asymm[col][it] = std::vector<double>();
      theVarsBinned_Asymm_stat[col][it] = std::vector<double>();
      theVarsBinned_A1Pl[col][it] = std::vector<double>();
      theVarsBinned_A1Pl_stat[col][it] = std::vector<double>();
      theVarsBinned_A1Mi[col][it] = std::vector<double>();
      theVarsBinned_A1Mi_stat[col][it] = std::vector<double>();
      theVarsBinned_A3[col][it] = std::vector<double>();
      theVarsBinned_A3_stat[col][it] = std::vector<double>();
    }
  }
  std::map< std::string, std::map< std::string, std::map< double, double > > > N , dN;
  std::map< std::string, double > norm;
  for (auto colM : theVars_inputs) {
    std::string col = colM.first;
    norm[col] = 0;
    for (auto it : colM.second) {
      anabin<0> thebin = it.first;
      if (!binok(thecats,xaxis,thebin)) continue;
      theBins.at(col).at(thebin).push_back(it.first);
      // Compute Charge Asymmetry
      asymm_input v = theVars_inputs.at(col).at(it.first);
      double asym = (v.nPl - v.nMi)/(v.nPl + v.nMi);
      double dasym_stat = 2.0*sqrt( pow((v.nMi*v.dnPl_stat), 2.0) + pow((v.nPl*v.dnMi_stat), 2.0) )/(pow((v.nPl + v.nMi), 2.0));
      theVarsBinned_Asymm.at(col).at(thebin).push_back(asym);
      theVarsBinned_Asymm_stat.at(col).at(thebin).push_back(dasym_stat);
      // Compute the eta asymmetries
      double eta = (it.first.muetabin().high()+it.first.muetabin().low())/2.;
      N[col]["Pl"][eta] = v.nPl; N[col]["Mi"][eta] = v.nMi;
      dN[col]["Pl"][eta] = v.dnPl_stat; dN[col]["Mi"][eta] = v.dnMi_stat;
      norm.at(col) += (v.nPl + v.nMi);
    }
  }
  // make TGraphAsymmErrors
  for (auto colM : theVars_inputs) {
    int cnt=0;
    std::string col = colM.first;
    for (auto it : thecats) {
      int n = theBins.at(col).at(it).size();
      if(n==0) {
        cout << "Error, nothing found for category" << endl;
        theGraphs_Asymm[col][it] = NULL;
        continue;
      }
      if (!theGraphs_Asymm[col][it]) { theGraphs_Asymm.at(col).at(it) = new TGraphAsymmErrors(n); }
      theGraphs_Asymm.at(col).at(it)->SetName(Form("bin_%s_%i",col.c_str(),cnt));
      for (int i=0; i<n; i++) {
        double x=0, exl=0, exh=0, y=0, eyl=0, eyh=0;
        double low=0, high=0; 
        anabin<0> thebin = theBins.at(col).at(it).at(i);
        y = theVarsBinned_Asymm.at(col).at(it).at(i);
        if (xaxis=="eta") {
          low= thebin.muetabin().low();
          high = thebin.muetabin().high();
          x = (low+high)/2.;
          exh = (high-low)/2.;
          exl = (high-low)/2.;
        }
        eyl = fabs(theVarsBinned_Asymm_stat.at(col).at(it).at(i));
        eyh = eyl;
        theGraphs_Asymm.at(col).at(it)->SetPoint(i,x,y);
        theGraphs_Asymm.at(col).at(it)->SetPointError(i,exl,exh,eyl,eyh);
      }
      cnt++;
    }
  }
  for (auto colM : theVars_inputs) {
    std::string col = colM.first;
    for (auto it : colM.second) {
      anabin<0> thebin = it.first;
      if (!binok(thecats,xaxis,thebin)) continue;
      double eta = (it.first.muetabin().high()+it.first.muetabin().low())/2.;
      if (eta<0.) continue;
      if (col.find("pPb")!=std::string::npos) { eta = -1.0*eta; }
      double a1Pl = N.at(col).at("Pl").at(-1.*eta)/N.at(col).at("Pl").at(eta);
      double da1Pl_stat = sqrt( 
                               pow((N.at(col).at("Pl").at(eta)*dN.at(col).at("Pl").at(-1.*eta)), 2.0) + 
                               pow((N.at(col).at("Pl").at(-1.*eta)*dN.at(col).at("Pl").at(eta)), 2.0) 
                                ) / (pow(N.at(col).at("Pl").at(eta), 2.0));
      theVarsBinned_A1Pl.at(col).at(thebin).push_back(a1Pl);
      theVarsBinned_A1Pl_stat.at(col).at(thebin).push_back(da1Pl_stat);
      double a1Mi = N.at(col).at("Mi").at(-1.*eta)/N.at(col).at("Mi").at(eta);
      double da1Mi_stat = sqrt( 
                               pow((N.at(col).at("Mi").at(eta)*dN.at(col).at("Mi").at(-1.*eta)), 2.0) + 
                               pow((N.at(col).at("Mi").at(-1.*eta)*dN.at(col).at("Mi").at(eta)), 2.0) 
                                ) / (pow(N.at(col).at("Mi").at(eta), 2.0));
      theVarsBinned_A1Mi.at(col).at(thebin).push_back(a1Mi);
      theVarsBinned_A1Mi_stat.at(col).at(thebin).push_back(da1Mi_stat);
      double a3 = (N.at(col).at("Mi").at(-1.*eta)+N.at(col).at("Pl").at(-1.*eta))/(N.at(col).at("Mi").at(eta)+N.at(col).at("Pl").at(eta));
      double da3_stat = sqrt( 
                             pow( (N.at(col).at("Pl").at(eta)+N.at(col).at("Mi").at(eta)) , 2.0)*(pow( (dN.at(col).at("Pl").at(-1.*eta)) , 2.0) + pow( (dN.at(col).at("Mi").at(-1.*eta)) , 2.0) ) +
                             pow( (N.at(col).at("Pl").at(-1.*eta)+N.at(col).at("Mi").at(-1.*eta)) , 2.0)*(pow( (dN.at(col).at("Pl").at(eta)) , 2.0) + pow( (dN.at(col).at("Mi").at(eta)) , 2.0) )
                              ) / (pow( (N.at(col).at("Mi").at(eta)+N.at(col).at("Pl").at(eta)) , 2.0));
      theVarsBinned_A3.at(col).at(thebin).push_back(a3);
      theVarsBinned_A3_stat.at(col).at(thebin).push_back(da3_stat);
    }
  }
  for (auto colM : theVars_inputs) {
    int cnt=0;
    std::string col = colM.first;
    for (auto it : thecats) {
      int n = int(theBins.at(col).at(it).size());
      if(n==0) {
        cout << "Error, nothing found for category" << endl;
        theGraphs_A1Pl[col][it] = NULL;
        theGraphs_A1Mi[col][it] = NULL;
        theGraphs_A3[col][it] = NULL;
        continue;
      }
      if (!theGraphs_A1Pl[col][it]) { theGraphs_A1Pl.at(col).at(it) = new TGraphAsymmErrors(int(n/2.)); }
      if (!theGraphs_A1Mi[col][it]) { theGraphs_A1Mi.at(col).at(it) = new TGraphAsymmErrors(int(n/2.)); }
      if (!theGraphs_A3[col][it]) { theGraphs_A3.at(col).at(it) = new TGraphAsymmErrors(int(n/2.)); }
      theGraphs_A1Pl.at(col).at(it)->SetName(Form("bin_%s_%i",col.c_str(),cnt));
      theGraphs_A1Mi.at(col).at(it)->SetName(Form("bin_%s_%i",col.c_str(),cnt));
      theGraphs_A3.at(col).at(it)->SetName(Form("bin_%s_%i",col.c_str(),cnt));
      int i=0;
      for (int j=0; j<n; j++) {
        double x=0, exl=0, exh=0, y=0, eyl=0, eyh=0;
        double low=0, high=0; 
        anabin<0> thebin = theBins.at(col).at(it).at(j);
        double eta = (thebin.muetabin().high()+thebin.muetabin().low())/2.;
        if (eta<0.) continue;
        y = theVarsBinned_A1Pl.at(col).at(it).at(i);
        if (xaxis=="eta") {
          low= thebin.muetabin().low();
          high = thebin.muetabin().high();
          x = (low+high)/2.;
          exh = (high-low)/2.;
          exl = (high-low)/2.;
        }
        eyl = fabs(theVarsBinned_A1Pl_stat.at(col).at(it).at(i));
        eyh = eyl;
        theGraphs_A1Pl.at(col).at(it)->SetPoint(i,x,y);
        theGraphs_A1Pl.at(col).at(it)->SetPointError(i,exl,exh,eyl,eyh);
        y = theVarsBinned_A1Mi.at(col).at(it).at(i);
        eyl = fabs(theVarsBinned_A1Mi_stat.at(col).at(it).at(i));
        eyh = eyl;
        theGraphs_A1Mi.at(col).at(it)->SetPoint(i,x,y);
        theGraphs_A1Mi.at(col).at(it)->SetPointError(i,exl,exh,eyl,eyh);
        y = theVarsBinned_A3.at(col).at(it).at(i);
        eyl = fabs(theVarsBinned_A3_stat.at(col).at(it).at(i));
        eyh = eyl;
        theGraphs_A3.at(col).at(it)->SetPoint(i,x,y);
        theGraphs_A3.at(col).at(it)->SetPointError(i,exl,exh,eyl,eyh);
        i += 1;
      }
      cnt++;
    }
  }
  // Adding MCFM values
  if (addTheory) {
    for (auto it : thecats) {
      TFile* fh = TFile::Open("Data/chasym_8tev.root");
      if (fh && fh->Get("chasym_8tev_th")) {
        TH1F h = *((TH1F*)fh->Get("chasym_8tev_th"));
        int nbin = h.GetNbinsX();
        for (int ibin=0; ibin<nbin; ibin++) {
          double xv = h.GetBinCenter(ibin+1), yv = h.GetBinContent(ibin+1), xerr = h.GetBinWidth(ibin+1)/2., yerr = h.GetBinError(ibin+1);
          if (!theGraphs_Asymm["MCFM"][it]) { theGraphs_Asymm.at("MCFM").at(it) = new TGraphAsymmErrors(nbin); }
          theGraphs_Asymm.at("MCFM").at(it)->SetPoint(ibin, xv, yv);
          theGraphs_Asymm.at("MCFM").at(it)->SetPointError(ibin, xerr, xerr, yerr, yerr);
        }
      }
      if (fh) { fh->Close(); delete fh; }
    }
  }
  // plot
  std::vector< std::string > text;
  text.push_back(Form("%.1f GeV/c #leq p_{T}", 25.));
  if (strcmp(cutLabel,"")) { text.push_back( formatCut(cutLabel, varEWQLabel) ); }
  ylabel = "(N^{+}-N^{-})/(N^{+}+N^{-})";
  nameTag = "_asymm";
  plotGraph(theGraphs_Asymm, xaxis, outputDir, "asymm", text, true);
  ylabel = "N^{+}(-#eta_{lab})/N^{+}(#eta_{lab})";
  nameTag = "_a1Pl";
  plotGraph(theGraphs_A1Pl, xaxis, outputDir, "a1Pl", text, true);
  ylabel = "N^{-}(-#eta_{lab})/N^{-}(#eta_{lab})";
  nameTag = "_a1Mi";
  plotGraph(theGraphs_A1Mi, xaxis, outputDir, "a1Mi", text, true);
  ylabel = "N(-#eta_{lab})/N(#eta_{lab})";
  nameTag = "_a3";
  plotGraph(theGraphs_A3, xaxis, outputDir, "a3", text, true);
  f->Close(); delete f;
}

void plotGraph(const TGraphBinDict& theGraphs, const std::string& xaxis, const std::string& outputDir, const std::string type, const std::vector<std::string>& infoLabel, const bool& merge)
{
  setTDRStyle();
  std::vector< anabin<0> > theCats;
  std::map< std::string , TCanvas* > c1;
  std::map< std::string , TH1F*    > haxes;
  std::map< std::string , TLegend* > tleg;
  bool first = true;
  std::string nom = "all";
  int cnt=0;
  for(auto token : theGraphs) {
    std::string tok = token.first;
    if (!merge) { nom = tok; }
    if (c1.count(nom)==0) { 
      c1[nom] = new TCanvas(Form("c1_%s",nom.c_str()),Form("c1_%s",nom.c_str()),600,600); 
      c1.at(nom)->cd();
    }
    // the axes
    std::string xname; 
    if (xaxis=="eta") { xname = "\\eta"; }
    if (haxes.count(nom)==0) {
      haxes[nom]=NULL; TLine line;
      if (xaxis=="eta") {
        if (type=="asymm") {
          haxes.at(nom) = new TH1F(Form("haxes_%s", nom.c_str()), Form("haxes_%s", nom.c_str()), 1, -2.5, 2.5);
          line = TLine(-2.5,0.,2.5,0.);
        } else {
          haxes.at(nom) = new TH1F(Form("haxes_%s", nom.c_str()), Form("haxes_%s", nom.c_str()), 1, 0.0, 2.5);
          line = TLine(0.0,0.,2.5,0.);
        }
      }
      if (type=="asymm") { haxes.at(nom)->GetYaxis()->SetRangeUser(-0.4,0.4); }
      else { haxes.at(nom)->GetYaxis()->SetRangeUser(0.6,1.8); }
      haxes.at(nom)->GetYaxis()->SetTitle(ylabel.c_str());
      std::string xlabel;
      if (xaxis=="eta") { xlabel = "#eta_{lab}"; }
      haxes.at(nom)->GetXaxis()->SetTitle(xlabel.c_str());
      haxes.at(nom)->GetXaxis()->SetNdivisions(505);
      haxes.at(nom)->Draw();
      line.Draw();
    }
    // the legend
    if (merge && tleg.count(nom)==0) {
      double xshift=0.;
      tleg[nom] = new TLegend(0.21+xshift,0.18,0.50,0.28);
      tleg.at(nom)->SetBorderSize(0);
      tleg.at(nom)->SetTextSize(0.035);
    }
    // prepare for the printing of the result tables
    gSystem->mkdir(Form("Output/%s/tex/", outputDir.c_str()), kTRUE); 
    char texname[2048]; sprintf(texname, "Output/%s/tex/result_%s%s_%s.tex",outputDir.c_str(),xaxis.c_str(),nameTag.c_str(),tok.c_str());
    string yname("$\\asym (\\W)$");
    inittex(texname, xname.c_str(), yname);
    if (!merge) { cnt = 0; }
    for (auto it : theGraphs.at(tok)) {
      anabin<0> thebin = it.first;
      TGraphAsymmErrors* tg = it.second;
      if (tok.find("pPb")==std::string::npos&&(tg->GetX()[0]<0.)) { flipXAxis(*tg); }
      theCats.push_back(thebin);
      tg->SetMarkerStyle(markerstyle(cnt));
      tg->SetMarkerColor(color(cnt));
      tg->SetLineColor(color(cnt));
      tg->SetMarkerSize(1.0);
      tg->SetLineWidth(tg->GetLineWidth()*2);
      gStyle->SetEndErrorSize(5);
      tg->Draw("P");
      if (merge) {
        if (tok.find("pPb")==std::string::npos) { tleg.at(nom)->AddEntry(tg, (tok+" (flipped)").c_str(), "p"); }
        else { tleg.at(nom)->AddEntry(tg, tok.c_str(), "p"); }
      }
      // print tex
      ostringstream oss;
      oss.precision(1); oss.setf(ios::fixed);
      oss << "$" << it.first.muetabin().low() << "<#eta<" << it.first.muetabin().high() << "$, ";
      addline(texname,oss.str());
      printGraph(tg, texname);
      double x=0, dx=0, y=0, dy=0;
      if (xaxis=="eta") { dx = 10; }
      x = 2*dx*cnt + dx;
      y = 1;
      cnt++;
    }
    if (!merge || (first&&merge)) {
      TLatex tl;
      double tlx = 0.2; double tly = 0.85, dy = 0.045;
      tl.SetTextSize(0.03);
      if (infoLabel.size()>0) { for(auto f : infoLabel) { tl.DrawLatexNDC(tlx, tly, f.c_str()); tly-=dy; } }
      int iPos = 33, lumiId = 0;
      if (merge) { lumiId = 111; }
      else if (tok.find("pPb")!=std::string::npos) { lumiId = 109; } 
      else if (tok.find("Pbp")!=std::string::npos) { lumiId = 110; } 
      else if (tok.find("PA")!=std::string::npos) { lumiId = 111; }
      CMS_lumi( (TPad*) gPad, lumiId, iPos, "" );
    }
    closetex(texname);
    first = false;
  }
  if (merge) { tleg.at(nom)->Draw(); }
  for (auto& c : c1) {
    c.second->cd();
    c.second->Update();
    c.second->RedrawAxis();
    gSystem->mkdir(Form("Output/%s/plot/RESULT/root/", outputDir.c_str()), kTRUE); 
    c.second->SaveAs(Form("Output/%s/plot/RESULT/root/result_%s%s_%s.root",outputDir.c_str(), xaxis.c_str(), nameTag.c_str(), c.first.c_str()));
    gSystem->mkdir(Form("Output/%s/plot/RESULT/png/", outputDir.c_str()), kTRUE);
    c.second->SaveAs(Form("Output/%s/plot/RESULT/png/result_%s%s_%s.png",outputDir.c_str(), xaxis.c_str(), nameTag.c_str(), c.first.c_str()));
    gSystem->mkdir(Form("Output/%s/plot/RESULT/pdf/", outputDir.c_str()), kTRUE);
    c.second->SaveAs(Form("Output/%s/plot/RESULT/pdf/result_%s%s_%s.pdf",outputDir.c_str(), xaxis.c_str(), nameTag.c_str(), c.first.c_str()));
  }
  for (auto& l : tleg)  { if (l.second) delete l.second; }
  for (auto& h : haxes) { if (h.second) delete h.second; }
  for (auto& c : c1)    { if (c.second) delete c.second; }
}

void flipXAxis(TGraphAsymmErrors& g)
{
  int nbins = g.GetN();
  std::vector< std::vector< double > > val;
  for (int ibin=0; ibin<nbins; ibin++) {
    val.push_back(std::vector<double>({ g.GetX()[ibin], g.GetY()[ibin], g.GetErrorXlow(ibin), g.GetErrorXhigh(ibin), g.GetErrorYlow(ibin), g.GetErrorYhigh(ibin) }));
  }
  for (int ibin=0; ibin<nbins; ibin++) {
    int jbin = (nbins-ibin-1);
    g.SetPoint(ibin, val.at(ibin).at(0), val.at(jbin).at(1));
    g.SetPointError(ibin, val.at(ibin).at(2), val.at(ibin).at(3), val.at(jbin).at(4), val.at(jbin).at(5));
  }
};
  
int color(int i)
{
   if (i==0) return kRed+2;
   else if (i==1) return kBlue+2;
   else if (i==2) return kGreen+2;
   else if (i==3) return kCyan+2;
   else if (i==4) return kMagenta+2;
   else if (i==5) return kOrange+2;
   else return kBlack;
};

int markerstyle(int i)
{
   if (i==0) return kFullSquare;
   else if (i==1) return kFullCircle;
   else if (i==2) return kOpenCircle;
   else if (i==3) return kFullCross;
   else if (i==4) return kOpenSquare;
   else if (i==5) return kFullStar;
   else if (i==6) return kOpenStar;
   else return kOpenCross;
};

void printOptions()
{
   cout << 
     "nameTag_base = \"" << nameTag_base << "\"" << 
     endl;
};

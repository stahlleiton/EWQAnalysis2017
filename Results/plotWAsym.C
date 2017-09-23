// Auxiliary Headers
#include "Utilities/resultsUtils.h"
#include "resultsEWQ2tree.C"
// ROOT headers
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
// RooFit headers

// c++ headers
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <string>


////////////////
// PARAMETERS //
////////////////

/*
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
void flipXAxis(TGraphAsymmErrors& g);
int color(int i);
int markerstyle(int i);
string nameTag;
*/

void plot ( const std::string xaxis = "eta" , const std::string workDirName = "NominalCM" );

void plotEta(const string& workDirName)
{
  std::string xaxis = "eta";
  std::vector< anabin<0> > bins;
  //
  plot(xaxis, workDirName);
};


void plotAll(const string& workDirName)
{
  plotEta(workDirName);
};


/////////////////////
// OTHER FUNCTIONS //
/////////////////////

void plot(const std::string xaxis, const std::string workDirName)
{
  //
  const std::string metTag      = "METPF_RAW";
  const std::string dsTag       = "DATA";
  const std::string colTag      = "PA";
  const std::string thePoiNames = "all";
  //
  // --------------------------------------------------------------------------------- //
  //
  // Define the input file info
  const std::string CWD = getcwd(NULL, 0);
  const std::string inputDirPath = Form("%s/Tree/%s/%s/%s/%s", CWD.c_str(), workDirName.c_str(), metTag.c_str(), dsTag.c_str(), colTag.c_str());
  const std::string inputFileName = "tree_allvars.root";
  const std::string inputFilePath = Form("%s/%s", inputDirPath.c_str(), inputFileName.c_str());
  // Open the input file
  TFile inputFile(inputFilePath.c_str(), "READ");
  if (inputFile.IsOpen()==false || inputFile.IsZombie()==true) {
    std::cout << "[WARNING] The input file " << inputFilePath << " was not found, will create it!" << std::endl;
    if (!resultsEWQ2tree(workDirName, metTag, dsTag, colTag)) { return; };
    inputFile.ReOpen("READ");
    if (inputFile.IsOpen()==false || inputFile.IsZombie()==true) {
      std::cout << "[ERROR] The input file " << inputFilePath << " could not be re-created!" << std::endl;
    }
  }
  //
  // --------------------------------------------------------------------------------- //
  //
  // Define the tree info container
  TreeInfo info;
  // Initialize the tree info container
  iniResultsTreeInfo(info, thePoiNames);
  //
  // Extract the input tree
  TTree * tree = (TTree*) inputFile.Get("fitResults");
  if (tree==NULL) { std::cout << "[ERROR] The input tree fitResults was not found in " << inputFilePath << "" << std::endl; inputFile.Close(); return; }
  // Set the address of the tree branches
  setBranchAddress(*tree, info);
  //
  // Extract the input variables
  VarBinMap inputVar;
  bool useEtaCM = false;
  //
  const uint nEntries = tree->GetEntries();
  for (uint i = 0; i < nEntries; i++) {
    tree->GetEntry(i);
    //
    const std::string collSystem = *info.StrP.at("collSystem");
    const std::string charge     = *info.StrP.at("charge");
    //
    // Determine the bin
    double etaMin = info.Var.at("VAR_Muon_Eta").at("Min");
    double etaMax = info.Var.at("VAR_Muon_Eta").at("Max");
    if ((etaMin == -99.) || (etaMax == -99.)) { std::cout << "[ERROR] The bin was not set properly!" << std::endl; return; }
    // If the bins where fitted in the Center of Mass, change from LAB to CM
    if (info.Flag.at("useEtaCM")) {
      bool usepPb = false; if ( (collSystem == "PA") || (collSystem == "pPb") ) { usepPb = true; }
      etaMin = PA::EtaLABtoCM(etaMin, usepPb);
      etaMax = PA::EtaLABtoCM(etaMax, usepPb);
      useEtaCM = true;
    }
    // Round the bin boundaries to two decimals
    roundValue(etaMin, 2);
    roundValue(etaMax, 2);
    anabin<0> bin(etaMin, etaMax);
    //
    // Get the fit variables
    for (auto& v : info.Var) {
      if (v.first.find("POI_")==std::string::npos) continue;
      std::string name = v.first; name.erase(name.find("POI_"),4);
      inputVar[collSystem][bin][charge][name]["Val"] = v.second.at("Val");
      inputVar[collSystem][bin][charge][name]["Err_Stat_High"] = v.second.at("Err");
      inputVar[collSystem][bin][charge][name]["Err_Stat_Low"]  = v.second.at("Err");
      inputVar[collSystem][bin][charge][name]["Err_Syst_High"] = 0.0;
      inputVar[collSystem][bin][charge][name]["Err_Syst_Low"]  = 0.0;
    }
    // Get the dataset variables
    for (auto& v : info.Var) {
      if (v.first.find("VAR_")==std::string::npos) continue;
      std::string name = v.first; name.erase(name.find("VAR_"),4);
      inputVar[collSystem][bin][charge][name]["Val"] = v.second.at("Val");
      inputVar[collSystem][bin][charge][name]["RMS"] = v.second.at("Err");
      inputVar[collSystem][bin][charge][name]["Min"] = v.second.at("Min");
      inputVar[collSystem][bin][charge][name]["Max"] = v.second.at("Max");
    }
  }
  //
  // --------------------------------------------------------------------------------- //
  //
  // Initialize the output variables
  BinPentaMap var;
  //
  const std::vector< std::string > varType = { "Var" , "Err_Stat_High" , "Err_Stat_Low" , "Err_Syst_High" , "Err_Syst_Low" };
  //
  // Compute the Charge Asymmetry
  //
  if (!computeChargeAsymmetry(var, inputVar)) { return; }
  //
  // Compute the Forward-Backward ratio
  //
  if (!computeForwardBackwardRatio(var, inputVar)) { return; }
  //
  // Compute the Cross-Section
  //
  if (!computeCrossSection(var, inputVar)) { return; }
  //
  // --------------------------------------------------------------------------------- //
  //
  // Initialize the Output Graphs
  GraphQuadMap graph;
  iniResultsGraph(graph, var);
  //
  // Fill the Output Graphs
  //
  if (!fillResultsGraph(graph, var)) { return; }
  //
  // --------------------------------------------------------------------------------- //
  //
  // Set Style
  setStyle();
  //
  // Draw the Output Graphs
  drawGraph(graph, CWD);
  //
  //
  inputFile.Close();
  /*
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
  */
}
/*
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
*/
/*
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
*/

#if !defined(__CINT__) || defined(__MAKECINT__)
// Auxiliary Headers
// ROOT headers
#include "TROOT.h"
#include "TSystem.h"
#include "TClass.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TEfficiency.h"
#include "TVectorD.h"
#include "TH1.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TFrame.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TPaletteAxis.h"
#include "TMath.h"
// c++ headers
#include <dirent.h>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <string>
// CMS headers
#include "../Utilities/CMS/tdrstyle.C"
#include "../Utilities/CMS/CMS_lumi.C"

#endif


// ------------------ TYPE -------------------------------
using Unc1DVec_t   =  std::map< std::string , TVectorD >;
using Unc1DMap2_t  =  std::map< std::string , std::map< std::string , std::map< std::string , Unc1DVec_t > > >;
using Unc1DMap_t   =  std::map< std::string , std::map< std::string , Unc1DMap2_t > >;
using EffVec_t     =  std::map< std::string , std::vector< TEfficiency > >;
using EffMap2_t    =  std::map< std::string , std::map< std::string , std::map< std::string , EffVec_t > > >;
using EffMap_t     =  std::map< std::string , std::map< std::string , EffMap2_t > >;
using BinMap_t     =  std::map< std::string , std::vector< double > >;
using BinMapMap_t  =  std::map< std::string , BinMap_t >;
using KeyPVec_t    =  std::vector< TKey* >;



// ------------------ FUNCTION -------------------------------
void         makeDir             ( const std::string& dir );
bool         existDir            ( const std::string& dir );
bool         getObjectsFromFile  ( EffMap_t& eff , Unc1DMap_t& unc , const std::string& filePath );
void         formatEff1D         ( TGraphAsymmErrors& graph , const std::string& col , const std::string& var , const std::string& charge , const std::string& type );
void         formatEff1D         ( TEfficiency& eff , const std::string& col , const std::string& var , const std::string& charge , const std::string& type );
void         drawEff1D           ( const std::string& outDir , EffMap_t& effMap , Unc1DMap_t& unc , const bool isCutAndCount = false );
void         createEff1DTable    ( std::vector< std::string >& texTable , const std::vector< std::string >& colTyp , const std::vector< std::string >& colCor ,
                                   const std::vector< std::string >& colTitle1 , const std::vector< std::string >& colChg ,
                                   const EffMap2_t& effMap , const Unc1DMap2_t& uncMap , const std::string& col );
void         makeCorrEff1DTable  ( std::ofstream& file    , const EffMap_t& iEffMap  , const Unc1DMap_t& iUncMap , const std::string& col );
void         makeMCEff1DTable    ( std::ofstream& file    , const EffMap_t& iEffMap  , const Unc1DMap_t& iUncMap , const std::string& col , const bool isHFReweight );
void         makeTnPUnc1DTable   ( std::ofstream& file    , const EffMap_t& iEffMap  , const Unc1DMap_t& iUncMap , const std::string& col , const std::string& chg );
bool         printEff1DTables    ( const EffMap_t& effMap , const Unc1DMap_t& uncMap , const std::string& outDir , const bool isHFReweight );

void         formatLegendEntry   ( TLegendEntry& e );
void         formatDecayLabel    ( std::string& label , const std::string& inLabel, const std::string c="#" );
void         setStyle            ( );
KeyPVec_t    getK                ( TList* list );
const char*  sgn                 ( const double n ) { if (n >= 0.) { return "+"; } else { return ""; } }


// ------------------ GLOBAL ------------------------------- 
//
// Muon Charge
const std::vector< std::string > CHG_  = { "Plus" , "Minus" };


void printEfficiency(const std::string workDirName = "NominalCM", const uint applyHFCorr = 1)
{
  // Change the working directory
  const std::string CWD = getcwd(NULL, 0);
  const std::string mainDir = Form("%s/Output/%s/", CWD.c_str(), (workDirName+((applyHFCorr==1) ? "_WithHF" : ((applyHFCorr==2) ? "_WithNTrack" : ""))).c_str());
  gSystem->mkdir(mainDir.c_str(), kTRUE);
  gSystem->ChangeDirectory(mainDir.c_str());
  //
  // Define the working flags
  bool isCutAndCount = false; if (workDirName.find("CutAndCount")!=std::string::npos) { isCutAndCount = true; }
  //
  // Declare the efficiencies
  EffMap_t eff1D;
  // Declare the TnP uncertainties
  Unc1DMap_t unc1D;
  //
  // Extract information from the input file
  const std::string inputFileName = "efficiencyTnP.root";
  if (!getObjectsFromFile(eff1D, unc1D, inputFileName)) { return; }
  //
  // ------------------------------------------------------------------------------------------------------------------------
  //
  // Draw the Efficiencies
  drawEff1D(mainDir, eff1D, unc1D, isCutAndCount);
  //
  // ------------------------------------------------------------------------------------------------------------------------
  //
  // Create Tables
  //if (!printTableEff1D(mainDir, eff1D, unc1D, isCutAndCount)) { return; };
  if (!printEff1DTables(eff1D, unc1D, mainDir, applyHFCorr)) { return; }
  //
};


void makeDir(const std::string& dir)
{
  if (existDir(dir.c_str())==false){ 
    std::cout << "[INFO] DataSet directory: " << dir << " doesn't exist, will create it!" << std::endl;  
    gSystem->mkdir(dir.c_str(), kTRUE);
  }
};


bool existDir(const std::string& dir)
{
  bool exist = false;
  void * dirp = gSystem->OpenDirectory(dir.c_str());
  if (dirp) { gSystem->FreeDirectory(dirp); exist = true; }
  return exist;
};


bool getObjectsFromFile(EffMap_t& eff, Unc1DMap_t& unc, const std::string& filePath)
{
  // Open the input file
  TFile file(filePath.c_str(), "READ");
  if (file.IsZombie()) { std::cout << "[ERROR] The input file " << filePath << " was not found!" << std::endl; return false; }
  file.cd();
  //
  //gain time, do not add the objects in the list in memory
  Bool_t status = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  //
  const std::string mainDirName = gDirectory->GetListOfKeys()->First()->GetName();
  TDirectory* mainDir = file.GetDirectory(mainDirName.c_str());
  if (mainDir==NULL) { std::cout << "[ERROR] Dir: " << mainDirName << " was not found!" << std::endl; return false; }
  mainDir->cd();
  // loop over all keys in this directory
  for (const auto& v : getK(mainDir->GetListOfKeys())) {
    if (v->ReadObj()->IsA()->InheritsFrom(TDirectory::Class()) == false) continue;
    TDirectory* varDir = mainDir->GetDirectory(v->GetName());
    varDir->cd();
    for (const auto& s : getK(varDir->GetListOfKeys())) {
      if (s->ReadObj()->IsA()->InheritsFrom(TDirectoryFile::Class()) == false) continue;
      TDirectory* sampleDir = varDir->GetDirectory(s->GetName());
      sampleDir->cd();
      for (const auto& c : getK(sampleDir->GetListOfKeys())) {
        if (c->ReadObj()->IsA()->InheritsFrom(TDirectoryFile::Class()) == false) continue;
        TDirectory* colDir = sampleDir->GetDirectory(c->GetName());
        colDir->cd();
        for (const auto& t : getK(colDir->GetListOfKeys())) {
          if (t->ReadObj()->IsA()->InheritsFrom(TDirectoryFile::Class()) == false) continue;
          TDirectory* typeDir = colDir->GetDirectory(t->GetName());
          typeDir->cd();
          for (const auto& co : getK(typeDir->GetListOfKeys())) {
            if (co->ReadObj()->IsA()->InheritsFrom(TDirectoryFile::Class())) {
              TDirectory* corDir = typeDir->GetDirectory(co->GetName());
              corDir->cd();
              for (auto& ch : CHG_) {
                for (const auto& key : getK(corDir->GetListOfKeys())) {
                  if (std::string(key->GetName()).find(ch)!=std::string::npos) {
                    if (key->ReadObj()->IsA()->InheritsFrom(TEfficiency::Class())) {
                      eff[v->GetName()][s->GetName()][c->GetName()][ch][t->GetName()][co->GetName()].push_back(*((TEfficiency*)key->ReadObj()));
                      const uint i = ( eff.at(v->GetName()).at(s->GetName()).at(c->GetName()).at(ch).at(t->GetName()).at(co->GetName()).size() - 1);
                      eff.at(v->GetName()).at(s->GetName()).at(c->GetName()).at(ch).at(t->GetName()).at(co->GetName())[i].SetName(key->GetName());
                    }
                    if (key->ReadObj()->IsA()->InheritsFrom(TVectorD::Class())) {
                      unc[v->GetName()][s->GetName()][c->GetName()][ch][t->GetName()][co->GetName()].ResizeTo( ((TVectorD*)key->ReadObj())->GetNrows() );
                      unc.at(v->GetName()).at(s->GetName()).at(c->GetName()).at(ch).at(t->GetName()).at(co->GetName()) = *((TVectorD*)key->ReadObj());
                    }
                  }
                }
                if (eff.at(v->GetName()).at(s->GetName()).at(c->GetName()).at(ch).at(t->GetName()).count("NoCorr")>0) {
                  eff[v->GetName()][s->GetName()][c->GetName()][ch][t->GetName()]["NoCorrOnly"] = eff.at(v->GetName()).at(s->GetName()).at(c->GetName()).at(ch).at(t->GetName()).at("NoCorr");
                  std::string name = eff.at(v->GetName()).at(s->GetName()).at(c->GetName()).at(ch).at(t->GetName()).at("NoCorr")[0].GetName();
                  name.replace(name.find("NoCorr"), name.find("NoCorr")+6, "NoCorrOnly");
                  eff.at(v->GetName()).at(s->GetName()).at(c->GetName()).at(ch).at(t->GetName()).at("NoCorrOnly")[0].SetName(name.c_str());
                }
                if (eff.at(v->GetName()).at(s->GetName()).at(c->GetName()).at(ch).at(t->GetName()).count("HFCorr")>0) {
                  eff[v->GetName()][s->GetName()][c->GetName()][ch][t->GetName()]["HFCorrOnly"] = eff.at(v->GetName()).at(s->GetName()).at(c->GetName()).at(ch).at(t->GetName()).at("HFCorr");
                  std::string name = eff.at(v->GetName()).at(s->GetName()).at(c->GetName()).at(ch).at(t->GetName()).at("HFCorr")[0].GetName();
                  name.replace(name.find("HFCorr"), name.find("HFCorr")+6, "HFCorrOnly");
                  eff.at(v->GetName()).at(s->GetName()).at(c->GetName()).at(ch).at(t->GetName()).at("HFCorrOnly")[0].SetName(name.c_str());
                }
              }
              typeDir->cd();
            }
            if (co->ReadObj()->IsA()->InheritsFrom(TVectorD::Class())) {
              for (auto& ch : CHG_) {
                if (std::string(co->GetName()).find(ch)!=std::string::npos) {
                  std::string corrName = "";
                  if (std::string(co->GetName()).find("TnP_Stat")!=std::string::npos) { corrName = "TnP_Stat"; }
                  if (std::string(co->GetName()).find("TnP_Syst")!=std::string::npos) { corrName = "TnP_Syst"; }
                  if (std::string(co->GetName()).find("TnP_Tot" )!=std::string::npos) { corrName = "TnP_Tot";  }
                  if (std::string(co->GetName()).find("MC_Syst" )!=std::string::npos) { corrName = "MC_Syst";  }
                  if (corrName!="") {
                    unc[v->GetName()][s->GetName()][c->GetName()][ch][t->GetName()][corrName].ResizeTo( ((TVectorD*)co->ReadObj())->GetNrows() );
                    unc.at(v->GetName()).at(s->GetName()).at(c->GetName()).at(ch).at(t->GetName()).at(corrName) = *((TVectorD*)co->ReadObj());
                  }
                }
              }
            }
          }
          colDir->cd();
        }
        sampleDir->cd();
      }
      varDir->cd();
    }
    mainDir->cd();
  }
  // Close File
  file.Close("R");
  // return
  return true;
};


void formatEff1D(TEfficiency& eff, const std::string& col, const std::string& var, const std::string& charge, const std::string& type)
{
  // Set the format of all graphs
  if (eff.GetDimension() != 1) { std::cout << "[ERROR] formatEff1D only works for dimension == 1" << std::endl; return; }
  gPad->Update();
  auto graph = eff.GetPaintedGraph();
  const std::string objLbl = ( (type=="Total" || type=="Acceptance") ? "Gen #mu" : "Reco #mu" );
  // X-axis
  std::string xLabel = objLbl; if (charge=="Plus") xLabel += "^{+}"; if (charge=="Minus") xLabel += "^{-}";
  if (var=="Eta"  ) { xLabel += " #eta_{LAB}";    }
  if (var=="EtaCM") { xLabel += " #eta_{CM}";     }
  if (var=="Pt"   ) { xLabel += " p_{T} (GeV/c)"; }
  // Y-axis
  std::string yLabel = type;
  if (type!="Acceptance") { type + " Efficiency"; }
  // Set Axis Titles
  eff.SetTitle(Form(";%s;%s", xLabel.c_str(), yLabel.c_str()));
  if (graph) { formatEff1D(*graph, col, var, charge, type); }
  gPad->Update(); 
};


void formatEff1D(TGraphAsymmErrors& graph, const std::string& col, const std::string& var, const std::string& charge, const std::string& type)
{
  // General
  graph.SetMarkerColor(kBlue);
  graph.SetMarkerStyle(20);
  graph.SetMarkerSize(1.0);
  graph.SetFillStyle(1001);
  const std::string objLbl = ( (type=="Total" || type=="Acceptance") ? "Gen #mu" : "Reco #mu" );
  // X-axis
  std::string xLabel = objLbl; if (charge=="Plus") xLabel += "^{+}"; if (charge=="Minus") xLabel += "^{-}";
  if (var=="Eta"  ) { xLabel += " #eta_{LAB}";    }
  if (var=="EtaCM") { xLabel += " #eta_{CM}";     }
  if (var=="Pt"   ) { xLabel += " p_{T} (GeV/c)"; }
  graph.GetXaxis()->CenterTitle(kFALSE);
  graph.GetXaxis()->SetTitleOffset(0.9);
  graph.GetXaxis()->SetTitleSize(0.050);
  graph.GetXaxis()->SetLabelSize(0.035);
  double xMin, xMax, yDummy;
  graph.GetPoint(0, xMin, yDummy); xMin -= graph.GetErrorXlow(0);
  int n = (graph.GetN()-1);
  graph.GetPoint(n, xMax, yDummy); xMax += graph.GetErrorXhigh(n);
  xMin = (std::floor((xMin-0.1)*10.0)/10.0);
  xMax = (std::ceil((xMax+0.1)*10.0)/10.0);
  graph.GetXaxis()->SetLimits(xMin , xMax);
  // Y-axis
  std::string yLabel = type;
  if (type!="Acceptance") { type + " Efficiency"; }
  graph.GetYaxis()->CenterTitle(kFALSE);
  graph.GetYaxis()->SetTitleOffset(1.4);
  graph.GetYaxis()->SetTitleSize(0.04);
  graph.GetYaxis()->SetLabelSize(0.035);
  if (type!="Acceptance") { graph.GetYaxis()->SetRangeUser(0.79, 1.13); }
  else { graph.GetYaxis()->SetRangeUser(0.60, 0.90); }
  // Set Axis Titles
  graph.SetTitle(Form(";%s;%s", xLabel.c_str(), yLabel.c_str()));
};


void drawEff1D(const std::string& outDir, EffMap_t& effMap, Unc1DMap_t& uncMap, const bool isCutAndCount)
{
  // Set Style
  setStyle();
  // Draw all graphs
  for (auto& v : effMap) {
    for (auto& s : v.second) {
      for (auto& cl : s.second) {
        for (auto& ch : cl.second) {
          for (auto& t : ch.second) {
            for (auto& co : t.second) {
              const std::string var    = v.first;
              const std::string sample = s.first;
              const std::string col    = cl.first;
              const std::string charge = ch.first;
              const std::string type   = t.first;
              const std::string corr   = co.first;
              std::vector<TEfficiency>& eff = co.second;
              // Create Canvas
              TCanvas c("c", "c", 1000, 1000); c.cd();
              // Create the Text Info
              TLatex tex; tex.SetNDC(); tex.SetTextSize(0.028); float dy = 0;
              std::vector< std::string > textToPrint;
              std::string sampleLabel; formatDecayLabel(sampleLabel, (sample+"_"+charge));
              textToPrint.push_back(sampleLabel);
              if (var=="Eta" || var=="EtaCM") { textToPrint.push_back("p^{#mu}_{T} > 25 GeV/c"); }
              if (var=="Pt") { textToPrint.push_back("|#eta^{#mu}| < 2.4"); }
              if (isCutAndCount) { textToPrint.push_back("20 GeV/c #geq |#slash{E}_{T}| & 40 GeV/c^{2} #geq M_{T}^{#mu}"); }
              //
              // Declare the graph vector (for drawing with markers)
              std::vector< TGraphAsymmErrors > grVec;
              TLegend leg(0.2, 0.66, 0.4, 0.79);
              double xMin , xMax;
              // Draw graph
              if ( (corr=="NoCorrOnly") || (corr=="HFCorrOnly") )  {
                // Extract the Uncorrected Efficiency graph
                eff[0].Draw(); gPad->Update();
                grVec.push_back(*(eff[0].GetPaintedGraph()));
                // Format the Graphs
                for (auto& g : grVec) { formatEff1D(g, col, var, charge, type); }
                // Create Legend
                if (type=="Acceptance") {
                  if (corr=="NoCorrOnly") { formatLegendEntry(*leg.AddEntry(&grVec[0], "MC Acceptance", "pe")); }
                  if (corr=="HFCorrOnly") { formatLegendEntry(*leg.AddEntry(&grVec[0], "MC Acceptance HF Reweighted", "pe")); }
                }
                else {
                  if (corr=="NoCorrOnly") { formatLegendEntry(*leg.AddEntry(&grVec[0], "MC Truth Efficiency", "pe")); }
                  if (corr=="HFCorrOnly") { formatLegendEntry(*leg.AddEntry(&grVec[0], "MC Truth Efficiency HF Reweighted", "pe")); }
                }
                // Draw the graph
                grVec[0].Draw("ap");
                xMin = grVec[0].GetXaxis()->GetXmin(); xMax = grVec[0].GetXaxis()->GetXmax();
              }
              else if (corr=="HFCorr" && type=="Acceptance") {
                // Extract the Uncorrected Efficiency graph
                t.second.at("NoCorr")[0].Draw(); gPad->Update();
                grVec.push_back(*(t.second.at("NoCorr")[0].GetPaintedGraph()));
                // Extract the HF Corrected Efficiency graph
                eff[0].Draw(); gPad->Update();
                grVec.push_back(*(eff[0].GetPaintedGraph()));
                // Format the Graphs
                for (auto& g : grVec) { formatEff1D(g, col, var, charge, type); }
                grVec[0].SetMarkerColor(kRed);
                // Create Legend
                formatLegendEntry(*leg.AddEntry(&grVec[0], "MC Acceptance", "pe"));
                formatLegendEntry(*leg.AddEntry(&grVec[1], "MC Acceptance HF Reweighted", "pe"));
                // Draw the graph
                grVec[0].Draw("ap");
                grVec[1].Draw("samep");
                xMin = grVec[0].GetXaxis()->GetXmin(); xMax = grVec[0].GetXaxis()->GetXmax();
              }
              else if (corr=="HFCorr") {
                // Extract the Uncorrected Efficiency graph
                t.second.at("NoCorr")[0].Draw(); gPad->Update();
                grVec.push_back(*(t.second.at("NoCorr")[0].GetPaintedGraph()));
                // Extract the HF Corrected Efficiency graph
                eff[0].Draw(); gPad->Update();
                grVec.push_back(*(eff[0].GetPaintedGraph()));
                // Extract the TnP Corrected Efficiency graph
                t.second.at("TnP_Nominal")[0].Draw(); gPad->Update();
                grVec.push_back(*(t.second.at("TnP_Nominal")[0].GetPaintedGraph()));
                // Fill Corrected Efficiency graph with total TnP Uncertainties
                Unc1DVec_t& unc = uncMap.at(v.first).at(s.first).at(cl.first).at(ch.first).at(t.first);
                for (int j = 0; j < grVec[2].GetN(); j++) {
                  const double errorYlow  = std::sqrt( std::pow(grVec[2].GetErrorYlow(j)  , 2.0) + std::pow(unc.at("TnP_Tot")[j] , 2.0) );
                  const double errorYhigh = std::sqrt( std::pow(grVec[2].GetErrorYhigh(j) , 2.0) + std::pow(unc.at("TnP_Tot")[j] , 2.0) );
                  grVec[2].SetPointError(j, grVec[1].GetErrorXlow(j), grVec[1].GetErrorXhigh(j), errorYlow, errorYhigh);
                }
                // Format the Graphs
                for (auto& g : grVec) { formatEff1D(g, col, var, charge, type); }
                grVec[0].SetMarkerColor(kRed);
                grVec[1].SetMarkerColor(kGreen+2);
                // Create Legend
                formatLegendEntry(*leg.AddEntry(&grVec[0], "MC Truth Efficiency", "pe"));
                formatLegendEntry(*leg.AddEntry(&grVec[1], "MC Truth Efficiency HF Reweighted", "pe"));
                formatLegendEntry(*leg.AddEntry(&grVec[2], "Corrected Efficiency", "pe"));
                // Draw the graph
                grVec[0].Draw("ap");
                grVec[1].Draw("samep");
                grVec[2].Draw("samep");
                xMin = grVec[0].GetXaxis()->GetXmin(); xMax = grVec[0].GetXaxis()->GetXmax();
              }
              else if (corr=="NoCorr" && type!="Acceptance") {
                // Extract the Uncorrected Efficiency graph
                eff[0].Draw(); gPad->Update();
                grVec.push_back(*(eff[0].GetPaintedGraph()));
                // Extract the Corrected Efficiency graph
                t.second.at("TnP_Nominal")[0].Draw(); gPad->Update();
                grVec.push_back(*(t.second.at("TnP_Nominal")[0].GetPaintedGraph()));
                // Fill Corrected Efficiency graph with total TnP Uncertainties
                Unc1DVec_t& unc = uncMap.at(v.first).at(s.first).at(cl.first).at(ch.first).at(t.first);
                for (int j = 0; j < grVec[1].GetN(); j++) {
                  const double errorYlow  = std::sqrt( std::pow(grVec[1].GetErrorYlow(j)  , 2.0) + std::pow(unc.at("TnP_Tot")[j] , 2.0) );
                  const double errorYhigh = std::sqrt( std::pow(grVec[1].GetErrorYhigh(j) , 2.0) + std::pow(unc.at("TnP_Tot")[j] , 2.0) );
                  grVec[1].SetPointError(j, grVec[1].GetErrorXlow(j), grVec[1].GetErrorXhigh(j), errorYlow, errorYhigh);
                }
                // Format the Graphs
                for (auto& g : grVec) { formatEff1D(g, col, var, charge, type); }
                grVec[0].SetMarkerColor(kRed);
                // Create Legend
                formatLegendEntry(*leg.AddEntry(&grVec[0], "MC Truth Efficiency", "pe"));
                formatLegendEntry(*leg.AddEntry(&grVec[1], "Corrected Efficiency", "pe"));
                // Draw the graph
                grVec[0].Draw("ap");
                grVec[1].Draw("samep");
                xMin = grVec[0].GetXaxis()->GetXmin(); xMax = grVec[0].GetXaxis()->GetXmax();
              }
              else if (corr=="NoCorr" && type=="Acceptance") continue;
              else if (corr=="TnP_Nominal") {
                // Extract the Corrected Efficiency graph
                eff[0].Draw(); gPad->Update();
                auto graph = eff[0].GetPaintedGraph();
                grVec.push_back(*graph);
                // Fill Corrected Efficiency graph with total TnP Uncertainties
                Unc1DVec_t& unc = uncMap.at(v.first).at(s.first).at(cl.first).at(ch.first).at(t.first);
                const bool useMCSyst = (unc.count("MC_Syst")>0);
                //
                for (int l = 0; l < grVec[0].GetN(); l++) {
                  const double errorYlow  = std::sqrt( (grVec[0].GetErrorYlow(l)*grVec[0].GetErrorYlow(l)) + (unc.at("TnP_Tot")[l]*unc.at("TnP_Tot")[l]) );
                  const double errorYhigh = std::sqrt( (grVec[0].GetErrorYhigh(l)*grVec[0].GetErrorYhigh(l)) + (unc.at("TnP_Tot")[l]*unc.at("TnP_Tot")[l]) );
                  grVec[0].SetPointError(l, grVec[0].GetErrorXlow(l), grVec[0].GetErrorXhigh(l), errorYlow, errorYhigh);
                }
                // Fill graph for TnP Statistical Uncertainties
                grVec.push_back(*graph);
                for (int l = 0; l < grVec[1].GetN(); l++) {
                  grVec[1].SetPointError(l, grVec[1].GetErrorXlow(l)*0.5, grVec[1].GetErrorXhigh(l)*0.5, unc.at("TnP_Stat")[l], unc.at("TnP_Stat")[l]);
                }
                // Fill graph for TnP Systematic Uncertainties
                grVec.push_back(*graph);
                for (int l = 0; l < grVec[2].GetN(); l++) {
                  grVec[2].SetPointError(l, grVec[2].GetErrorXlow(l)*0.5, grVec[2].GetErrorXhigh(l)*0.5, unc.at("TnP_Syst")[l], unc.at("TnP_Syst")[l]);
                }
                // Fill graph for MC Statistical Uncertainties
                grVec.push_back(*graph);
                for (int l = 0; l < grVec[3].GetN(); l++) {
                  grVec[3].SetPointError(l, grVec[3].GetErrorXlow(l)*0.5, grVec[3].GetErrorXhigh(l)*0.5, grVec[3].GetErrorYlow(l), grVec[3].GetErrorYhigh(l));
                }
                // Fill graph for MC Systematic Uncertainties
                if (useMCSyst) {
                  grVec.push_back(*graph);
                  for (int l = 0; l < grVec[4].GetN(); l++) {
                    grVec[4].SetPointError(l, grVec[4].GetErrorXlow(l)*0.5, grVec[4].GetErrorXhigh(l)*0.5, unc.at("MC_Syst")[l], unc.at("MC_Syst")[l]);
                  }
                }
                // Create Legend
                formatLegendEntry(*leg.AddEntry(&grVec[0], "Corrected Efficiency", "pe"));
                formatLegendEntry(*leg.AddEntry(&grVec[3], "MC Statistical Uncertainty", "f"));
                if (useMCSyst) { formatLegendEntry(*leg.AddEntry(&grVec[4], "MC Systematic Uncertainty", "f")); }
                formatLegendEntry(*leg.AddEntry(&grVec[1], "TnP Statistical Uncertainty", "f"));
                formatLegendEntry(*leg.AddEntry(&grVec[2], "TnP Systematic Uncertainty", "f"));
                // Format the graphs
                for (auto& g : grVec) { formatEff1D(g, col, var, charge, type); }
                grVec[0].SetMarkerColor(kBlack);
                grVec[1].SetFillColor(kOrange);
                grVec[2].SetFillColor(kGreen+3);
                grVec[3].SetFillColor(kRed);
                if (useMCSyst) { grVec[4].SetFillColor(kBlue+2); }
                // Draw the graphs
                grVec[0].Draw("ap");
                grVec[2].Draw("same2");
                grVec[1].Draw("same2");
                grVec[3].Draw("same2");
                if (useMCSyst) { grVec[4].Draw("same2");; }
                grVec[0].Draw("samep");
                xMin = grVec[0].GetXaxis()->GetXmin(); xMax = grVec[0].GetXaxis()->GetXmax();
              }
              else {
                // Extract the Corrected Efficiency graph
                t.second.at("TnP_Nominal")[0].Draw(); gPad->Update();
                auto graph = t.second.at("TnP_Nominal")[0].GetPaintedGraph();
                grVec.push_back(*graph);
                // Fill Corrected Efficiency graph with the TnP Uncertainties
                Unc1DVec_t& unc = uncMap.at(v.first).at(s.first).at(cl.first).at(ch.first).at(t.first);
                grVec.push_back(*graph);
                for (int j = 0; j < grVec[1].GetN(); j++) {
                  grVec[1].SetPointError(j, grVec[1].GetErrorXlow(j)*0.8, grVec[1].GetErrorXhigh(j)*0.8, unc.at(corr)[j], unc.at(corr)[j]);
                }
                // Extract the Varied Efficiency graphs
                for (uint i = 0; i < eff.size(); i++) {
                  eff[i].Draw(); gPad->Update();
                  grVec.push_back(*(eff[i].GetPaintedGraph()));
                  for (int j = 0; j < grVec[i+2].GetN(); j++) {
                    grVec[i+2].SetPointError(j, grVec[i+2].GetErrorXlow(j), grVec[i+2].GetErrorXhigh(j), 0.0, 0.0);
                  }
                }
                // Format the graphs
                for (auto& g : grVec) { formatEff1D(g, col, var, charge, type); }
                grVec[0].SetMarkerColor(kAzure+10);
                grVec[1].SetFillColor(kRed);
                for (uint n=2; n<grVec.size(); n++) { grVec[n].SetMarkerColor(kBlack);  grVec[n].SetMarkerSize(0.0); grVec[n].SetLineColor(kBlack); }
                // Create Legend
                formatLegendEntry(*leg.AddEntry(&grVec[0], "Corrected Efficiency", "pe"));
                formatLegendEntry(*leg.AddEntry(&grVec[1], Form("%s Uncertainty", corr.c_str()), "f"));
                formatLegendEntry(*leg.AddEntry(&grVec[2], Form("%s Variation", corr.c_str()), "l"));
                // Draw the Graph
                grVec[0].Draw("ap");
                for (uint n=2; n<grVec.size(); n++) { grVec[n].Draw("samep"); }
                grVec[1].Draw("same2");
                grVec[0].Draw("samep");
                xMin = grVec[0].GetXaxis()->GetXmin(); xMax = grVec[0].GetXaxis()->GetXmax();
              }
              // Draw Legend and line
              TLine line( xMin , 1.0 ,  xMax , 1.0 ); line.SetLineStyle(2);
              line.Draw("same");
              leg.Draw("same");
              // Update
              c.Modified(); c.Update();
              //
              // Draw the text
              for (const auto& s: textToPrint) { tex.DrawLatex(0.22, 0.86-dy, s.c_str()); dy+=0.04; }
              c.Modified(); c.Update();
              // set the CMS style
              int option = 114;
              if (col.find("pPb")!=std::string::npos) option = 112;
              if (col.find("Pbp")!=std::string::npos) option = 113;
              CMS_lumi(&c, option, 33, "");
              c.Modified(); c.Update();
              // Create Output Directory
              const std::string plotDir = outDir + "TnPEfficiencyLATEX/" + var+"/" + sample+"/" + col+"/" + type;
              makeDir(plotDir + "/png/");
              makeDir(plotDir + "/pdf/");
              makeDir(plotDir + "/root/");
              // Save Canvas
              const std::string ename = "eff1D_" + v.first +"_"+ s.first +"_"+ cl.first +"_"+ ch.first +"_"+ t.first +"_"+ co.first;
              c.SaveAs(( plotDir + "/png/" + ename + ".png" ).c_str());
              c.SaveAs(( plotDir + "/pdf/" + ename + ".pdf" ).c_str());
              c.SaveAs(( plotDir + "/root/" + ename + ".root" ).c_str());
              // Clean up memory
              c.Clear(); c.Close();
            }
          }
        }
      }
    }
  }
};
 
 
void setStyle()
{
  // Set the CMS style
  setTDRStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  //
};


void formatLegendEntry(TLegendEntry& e)
{
  e.SetTextSize(0.028);
};


 void formatDecayLabel(std::string& label, const std::string& inLabel, const std::string c)
{
  label = "";
  if (inLabel.find("MC_DYToMuMu")!=std::string::npos) {
    label = "Z/"+c+"gamma* "+c+"rightarrow "+c+"mu^{+} + "+c+"mu^{-}";
    if (inLabel.find("M_10_30")!=std::string::npos) { label = "Z/"+c+"gamma* "+c+"rightarrow "+c+"mu^{+} + "+c+"mu^{-} [10 < M < 30]"; }
    if (inLabel.find("M_30_")  !=std::string::npos) { label = "Z/"+c+"gamma* "+c+"rightarrow "+c+"mu^{+} + "+c+"mu^{-} [M > 30] ";      }
  }
  if (inLabel.find("MC_ZToMuMu")!=std::string::npos) {
    label = "Z "+c+"rightarrow "+c+"mu^{+} + "+c+"mu^{-}";
    if (inLabel.find("M_10_30")!=std::string::npos) { label = "Z "+c+"rightarrow "+c+"mu^{+} + "+c+"mu^{-} [10 < M < 30]"; }
    if (inLabel.find("M_30_")  !=std::string::npos) { label = "Z "+c+"rightarrow "+c+"mu^{+} + "+c+"mu^{-} [M > 30] ";      }
  }
  if (inLabel.find("MC_WToMu")!=std::string::npos) {
     label = "W "+c+"rightarrow "+c+"mu + "+c+"nu_{"+c+"mu}";
    if (inLabel.find("Plus") !=std::string::npos) { label = "W^{+} "+c+"rightarrow "+c+"mu^{+} + "+c+"nu_{"+c+"mu}"; }
    if (inLabel.find("Minus")!=std::string::npos) { label = "W^{-} "+c+"rightarrow "+c+"mu^{-} + "+c+"bar{"+c+"nu}_{"+c+"mu}"; }
  }
  if (inLabel.find("MC_WToTau")!=std::string::npos) {
     label = "W "+c+"rightarrow "+c+"tau "+c+"rightarrow "+c+"mu + "+c+"nu_{"+c+"mu} + "+c+"bar{"+c+"nu}_{"+c+"tau}";
    if (inLabel.find("Plus") !=std::string::npos) { label = "W^{+} "+c+"rightarrow "+c+"tau^{+} "+c+"rightarrow "+c+"mu^{+} + "+c+"nu_{"+c+"mu} + "+c+"bar{"+c+"nu}_{"+c+"tau}"; }
    if (inLabel.find("Minus")!=std::string::npos) { label = "W^{-} "+c+"rightarrow "+c+"tau^{-} "+c+"rightarrow "+c+"mu^{-} + "+c+"bar{"+c+"nu}_{"+c+"mu} + "+c+"nu_{"+c+"tau}"; }
  }
  if (inLabel == "MC_QCDToMu") { label = "QCD "+c+"rightarrow "+c+"mu"; }
  if (inLabel == "MC_TTall") { label = "t + "+c+"bar{t} "+c+"rightarrow All "+c+"rightarrow "+c+"mu"; }
  if (c!="#") { while(label.find(" \\")!=std::string::npos) { label.replace(label.find(" "), 1, "~"); } }
};


KeyPVec_t getK(TList* list)
{
  TIter iter(list);
  KeyPVec_t out;
  while (TKey* key = (TKey*)iter()) { out.push_back(key); }
  return out;
};


void createEff1DTable(std::vector< std::string >& texTable, const std::vector< std::string >& colTyp, const std::vector< std::string >& colCor, const std::vector< std::string >& colTitle1,
                      const std::vector< std::string >& colChg, const EffMap2_t& effMap, const Unc1DMap2_t& uncMap, const std::string& col)
{
  //
  const uint nCol = colTyp.size();
  //
  texTable.push_back("  \\renewcommand{\\arraystretch}{1.5}");
  texTable.push_back(Form("  \\begin{tabular}{|c|*%dc|}", (nCol-1)));
  texTable.push_back("    \\hline");
  std::string tmp;
  tmp = ("    ");
  for (uint i = 0; i < nCol; i++) {
    tmp += colTitle1[i];
    if (i<(nCol-1)) { tmp += " & "; }
    else { tmp += "\\\\"; }
  }
  texTable.push_back(tmp);
  texTable.push_back("    \\hline\\hline");
  //
  const auto& nomEff = effMap.at(col).at("Plus").at("Total").at("NoCorr")[0];
  //
  for (int iBin = 1; iBin <= nomEff.GetTotalHistogram()->GetNbinsX(); iBin++) {
    tmp = ("    ");
    for (uint i = 0; i < nCol; i++) {
      const auto&  typ = colTyp[i];
      const auto&  cor = colCor[i];
      const auto&  chg = colChg[i];
      std::string val;
      if (cor=="NoCorr" || cor=="HFCorr") {
        const auto& eff = effMap.at(col).at(chg).at(typ).at(cor)[0];
        const double errLow = eff.GetEfficiencyErrorLow(iBin);
        const double errHi  = eff.GetEfficiencyErrorUp(iBin);
        const double rVal   = eff.GetEfficiency(iBin);
        if (std::abs(errLow-errHi)<0.00001) { val = Form("$%.3f \\pm %.3f$", rVal, errLow ); }
        else { val = Form("$%.3f^{+%.3f}_{-%.3f}$", rVal, errHi, errLow ); }
      }
      else if (cor=="TnP_Nominal") {
        const auto& eff = effMap.at(col).at(chg).at(typ).at(cor)[0];
        const double errLow = eff.GetEfficiencyErrorLow(iBin);
        const double errHi  = eff.GetEfficiencyErrorUp(iBin);
        const double rVal   = eff.GetEfficiency(iBin);
        const auto& tnpStat = uncMap.at(col).at(chg).at(typ).at("TnP_Stat");
        const auto& tnpSyst = uncMap.at(col).at(chg).at(typ).at("TnP_Syst");
        const double tnpErr = TMath::Sqrt(TMath::Power(tnpStat[iBin-1], 2.0) + TMath::Power(tnpSyst[iBin-1], 2.0));
        if (std::abs(errLow-errHi)<0.00001) { val = Form("$%.3f \\pm %.3f", rVal, errLow ); }
        else { val = Form("$%.3f^{+%.3f}_{-%.3f}", rVal, errHi, errLow ); }
        val += Form(" \\pm %.3f \\textrm{ (tnp)}$", tnpErr);
      }
      else {
        const double min = nomEff.GetTotalHistogram()->GetXaxis()->GetBinLowEdge(iBin);
        const double max = nomEff.GetTotalHistogram()->GetXaxis()->GetBinUpEdge(iBin);
        val = Form("%s%.2f , %s%.2f", sgn(min), min, sgn(max) , max);
      }
      tmp += val;
      if (i<(nCol-1)) { tmp += " & "; } else { tmp += "\\\\"; }
    }
    texTable.push_back(tmp);
    texTable.push_back("    \\hline");
  }
  //
  texTable.push_back("  \\end{tabular}");
  //
};


void makeCorrEff1DTable(std::ofstream& file, const EffMap_t& iEffMap, const Unc1DMap_t& iUncMap, const std::string& col)
{
  //
  bool useEtaCM = false; if (iEffMap.count("EtaCM")>0) { useEtaCM = true; }
  std::string varLbl = "EtaCM"; if (!useEtaCM) { varLbl = "Eta"; }
  if (iEffMap.count(varLbl)==0) { std::cout << "[ERROR] Variable " << varLbl << " was not found" << std::endl; return; }
  // Draw all graphs
  const auto& effMap = iEffMap.at(varLbl).at("MC_WToMuNu");
  const auto& uncMap = iUncMap.at(varLbl).at("MC_WToMuNu");
  //
  std::string pTCUT = "$p_{T} > 25$~GeV/c";
  //
  // Initialize the latex table
  std::vector< std::string > texTable;
  //
  // Determine number of columns
  const std::vector< std::string > colChg = { "" , "Minus" , "Plus" };
  const std::vector< std::string > colCor = { "" , "TnP_Nominal" , "TnP_Nominal" };
  const std::vector< std::string > colTyp = { "VAR" , "Total" , "Total" };
  std::vector< std::string >    colTitle1 = { "$\\eta_{LAB}$ Range" , "$\\mu^{-}$ Corrected Eff." , "$\\mu^{+}$ Corrected Eff." };
  if (useEtaCM) { colTitle1[0] = "$\\eta_{CM}$ Range"; }
  //
  texTable.push_back("\\begin{table}[h!]");
  texTable.push_back("  \\centering");
  //texTable.push_back("  \\resizebox{\\textwidth}{!}{");
  createEff1DTable(texTable, colTyp, colCor, colTitle1, colChg, effMap, uncMap, col);
  //texTable.push_back("  }");
  texTable.push_back(Form("  \\caption{%s}",
                          Form("Muon corrected efficiency as a function of the generated %s, derived from the $\\WToMuNu$ %s \\POWHEG samples separated in negative and positive charged muons. %sGenerated muons are required to be matched to reconstructed muons passing all analysis cuts. The muon efficiency has been corrected by applying the Tag and Probe scale factors event by event. The label \"tnp\" represent the TnP total uncertainty.",
                               (useEtaCM ? "muon $\\eta_{CM}$" : "$\\eta_{LAB}$"),
                               col.c_str(),
                               ( (col=="PA") ? " The \\pPb and \\Pbp MC samples are combined as described in \\sect{sec:CombiningBeamDirection}. " : "")
                               )
                          )
                     );
  texTable.push_back(Form("  \\label{tab:corrEfficiency_WToMu_%s}", col.c_str()));
  texTable.push_back("\\end{table}");
  //
  for (const auto& row : texTable) { file << row << std::endl; }
  file << std::endl; file << std::endl;
  //
};


void makeMCEff1DTable(std::ofstream& file, const EffMap_t& iEffMap, const Unc1DMap_t& iUncMap, const std::string& col, const bool isHFReweight)
{
  //
  bool useEtaCM = false; if (iEffMap.count("EtaCM")>0) { useEtaCM = true; }
  std::string varLbl = "EtaCM"; if (!useEtaCM) { varLbl = "Eta"; }
  if (iEffMap.count(varLbl)==0) { std::cout << "[ERROR] Variable " << varLbl << " was not found" << std::endl; return; }
  // Draw all graphs
  const auto& effMap = iEffMap.at(varLbl).at("MC_WToMuNu");
  const auto& uncMap = iUncMap.at(varLbl).at("MC_WToMuNu");
  //
  std::string pTCUT = "$p_{T} > 25$~GeV/c";
  //
  // Initialize the latex table
  std::vector< std::string > texTable;
  //
  // Determine number of columns
  const std::vector< std::string > colChg = (isHFReweight ? std::vector< std::string >({ "" , "Minus" , "Minus" , "Plus" , "Plus" }) : std::vector< std::string >({ "" , "Minus" , "Plus" }));
  const std::vector< std::string > colCor = (isHFReweight ? std::vector< std::string >({ "" , "NoCorr" , "HFCorr" , "NoCorr" , "HFCorr" }) : std::vector< std::string >({ "" , "NoCorr" , "NoCorr" }));
  const std::vector< std::string > colTyp = (isHFReweight ? std::vector< std::string >({ "VAR" , "Total" , "Total" , "Total" , "Total" }) : std::vector< std::string >({ "VAR" , "Total" , "Total" }));
  std::vector< std::string >    colTitle1 = (isHFReweight ? std::vector< std::string >({ "$\\eta_{LAB}$ Range" , "$\\mu^{-}$ Truth Eff." , "$\\mu^{-}$ Reweighted Eff." , "$\\mu^{+}$ Truth Eff." , "$\\mu^{+}$ Reweighted Eff." }) : std::vector< std::string >({ "$\\eta_{LAB}$ Range" , "$\\mu^{-}$ Truth Eff." , "$\\mu^{+}$ Truth Eff." }));
if (useEtaCM) { colTitle1[0] = "$\\eta_{CM}$ Range"; }
  //
  texTable.push_back("\\begin{table}[h!]");
  texTable.push_back("  \\centering");
  if (colChg.size()>3) { texTable.push_back("  \\resizebox{\\textwidth}{!}{"); }
  createEff1DTable(texTable, colTyp, colCor, colTitle1, colChg, effMap, uncMap, col);
  if (colChg.size()>3) { texTable.push_back("  }"); }
  texTable.push_back(Form("  \\caption{%s}",
                          Form("Muon truth efficiency as a function of the generated %s, derived from the $\\WToMuNu$ %s \\POWHEG samples separated in negative and positive charged muons. %sGenerated muons are required to be matched to reconstructed muons passing all analysis cuts.%s Results corresponds to \\eq{eq:MCTruthEfficiency}",
                               (useEtaCM ? "muon $\\eta_{CM}$" : "$\\eta_{LAB}$"),
                               col.c_str(),
                               ( (col=="PA") ? "The \\pPb and \\Pbp MC samples are combined as described in \\sect{sec:CombiningBeamDirection}. " : ""),
                               (isHFReweight ? " The event activity of the MC samples has been re-weighted." : "")
                               )
                          )
                     );
  texTable.push_back(Form("  \\label{tab:mcEfficiency_WToMu_%s}", col.c_str()));
  texTable.push_back("\\end{table}");
  //
  for (const auto& row : texTable) { file << row << std::endl; }
  file << std::endl; file << std::endl;
  //
};


void createUnc1DTable(std::vector< std::string >& texTable, const std::vector< std::string >& colTyp, const std::vector< std::string >& colCor, const std::vector< std::string >& colTitle1,
                      const std::vector< std::string >& colChg, const EffMap2_t& effMap, const Unc1DMap2_t& uncMap, const std::string& col)
{
  //
  const uint nCol = colTyp.size();
  //
  texTable.push_back("  \\renewcommand{\\arraystretch}{1.5}");
  texTable.push_back(Form("  \\begin{tabular}{|c|*%dc|}", (nCol-1)));
  texTable.push_back("    \\hline");
  std::string tmp;
  tmp = ("    ");
  for (uint i = 0; i < nCol; i++) {
    tmp += colTitle1[i];
    if (i<(nCol-1)) { tmp += " & "; }
    else { tmp += "\\\\"; }
  }
  texTable.push_back(tmp);
  texTable.push_back("    \\hline\\hline");
  //
  const auto& nomUnc = uncMap.at(col).at("Plus").at("Total").at("TnP_Stat");
  //
  for (int iBin = 1; iBin <= nomUnc.GetNrows(); iBin++) {
    tmp = ("    ");
    for (uint i = 0; i < nCol; i++) {
      const auto&  typ = colTyp[i];
      const auto&  cor = colCor[i];
      const auto&  chg = colChg[i];
      std::string val;
      if (cor=="TnP_Syst_STA") {
        const auto& nomEff = effMap.at(col).at("Plus").at("Total").at("TnP_Nominal")[0];
        const double tnpErr = 0.0060*nomEff.GetEfficiency(iBin);
        val = Form("%.4f", tnpErr );
      }
      else  if (cor=="TnP_Syst_PU") {
        const auto& nomEff = effMap.at(col).at("Plus").at("Total").at("TnP_Nominal")[0];
        const double tnpErr = 0.0034*nomEff.GetEfficiency(iBin);
        val = Form("%.4f", tnpErr );
      }
      else if (cor=="TnP_Stat_MuID" || cor=="TnP_Stat_Iso" || cor=="TnP_Stat_Trig") {
        double tnpErr = 0.0;
        if (uncMap.at(col).at(chg).at(typ).count(cor+"_p21_p24")>0) {
          for (const auto& cr : uncMap.at(col).at(chg).at(typ)) {
            if (cr.first.find(cor+"_")!=std::string::npos) {
              tnpErr += std::pow( cr.second[iBin-1] , 2.0 );
            }
          }
          tnpErr = std::sqrt( tnpErr );
        }
        else {
          tnpErr = uncMap.at(col).at(chg).at(typ).at(cor)[iBin-1];
        }
        val = Form("%.4f", tnpErr );
      }
      else if (typ=="Total") {
        const auto& tnpErr = uncMap.at(col).at(chg).at(typ).at(cor)[iBin-1];
        val = Form("%.4f", tnpErr );
      }
      else {
        const auto& nomEff = effMap.at(col).at("Plus").at("Total").at("NoCorr")[0];
        const double min = nomEff.GetTotalHistogram()->GetXaxis()->GetBinLowEdge(iBin);
        const double max = nomEff.GetTotalHistogram()->GetXaxis()->GetBinUpEdge(iBin);
        val = Form("%s%.2f , %s%.2f", sgn(min), min, sgn(max) , max);
      }
      tmp += val;
      if (i<(nCol-1)) { tmp += " & "; } else { tmp += "\\\\"; }
    }
    texTable.push_back(tmp);
    texTable.push_back("    \\hline");
  }
  //
  texTable.push_back("  \\end{tabular}");
  //
};


void makeTnPUnc1DTable(std::ofstream& file, const EffMap_t& iEffMap, const Unc1DMap_t& iUncMap, const std::string& col, const std::string& chg)
{
  //
  bool useEtaCM = false; if (iEffMap.count("EtaCM")>0) { useEtaCM = true; }
  std::string varLbl = "EtaCM"; if (!useEtaCM) { varLbl = "Eta"; }
  if (iUncMap.count(varLbl)==0) { std::cout << "[ERROR] Variable " << varLbl << " was not found" << std::endl; return; }
  // Draw all graphs
  const auto& uncMap = iUncMap.at(varLbl).at("MC_WToMuNu");
  const auto& effMap = iEffMap.at(varLbl).at("MC_WToMuNu");
  //
  std::string pTCUT = "$p_{T} > 25$~GeV/c";
  //
  // Initialize the latex table
  std::vector< std::string > texTable;
  //
  // Determine number of columns
  std::vector< std::string > colChg = { "" , chg , chg , chg , chg , chg , chg , chg , chg };
  std::vector< std::string > colCor = { "" , "TnP_Syst_MuID" , "TnP_Syst_Iso" , "TnP_Syst_Trig" , "TnP_Syst_BinIso" , "TnP_Syst_BinMuID" , "TnP_Syst_STA" , "TnP_Syst_PU" , "TnP_Syst" };
  std::vector< std::string > colTyp = { "VAR" , "Total" , "Total" , "Total" , "Total" , "Total" , "Total" , "Total" , "Total" };
  std::vector< std::string > colTitle1 = { "$\\eta_{LAB}$ Range" , "MuID" , "Iso" , "Trig" , "Iso Binned" , "MuID Binned" , "STA" , "PU" , "Total" };
  if (useEtaCM) { colTitle1[0] = "$\\eta_{CM}$ Range"; }
  //
  texTable.push_back("\\begin{table}[h!]");
  texTable.push_back("  \\centering");
  texTable.push_back("  \\resizebox{\\textwidth}{!}{");
  createUnc1DTable(texTable, colTyp, colCor, colTitle1, colChg, effMap, uncMap, col);
  texTable.push_back("  }");
  texTable.push_back(Form("  \\caption{%s}",
                          Form("Systematic uncertainties corresponding to the muon efficiency corrections derived applying the Tag and Probe scale factors event by event. Each source  for each source  The errors are shown as a function of the generated %s, derived from the %s %s \\POWHEG samples separated in negative and positive charged muons. %sGenerated muons are require to match reconstructed muons passing all analysis cuts.",
                               (useEtaCM ? "muon $\\eta_{CM}$" : "$\\eta_{LAB}$"),
                               (chg=="Plus" ? "$\\WToMuNuPl$" : "$\\WToMuNuMi$"),
                               col.c_str(),
                               ( (col=="PA") ? "The \\pPb and \\Pbp MC samples are combined as described in \\sect{sec:CombiningBeamDirection}. " : "")
                               )
                          )
                     );
  texTable.push_back(Form("  \\label{tab:tnpSystUncertainty_WToMu_%s_%s}", chg.c_str(), col.c_str()));
  texTable.push_back("\\end{table}");
  //
  for (const auto& row : texTable) { file << row << std::endl; }
  file << std::endl; file << std::endl;
  //
  //
  // Initialize the latex table
  texTable.clear(); colChg.clear(); colCor.clear(); colTyp.clear(); colTitle1.clear();
  //
  // Determine number of columns
  colChg = std::vector< std::string >({ "" , chg , chg , chg , chg });
  colCor = std::vector< std::string >({ "" , "TnP_Stat_MuID" , "TnP_Stat_Iso" , "TnP_Stat_Trig" , "TnP_Stat" });
  colTyp = std::vector< std::string >({ "VAR" , "Total" , "Total" , "Total" , "Total" });
  colTitle1 = std::vector< std::string >({ "$\\eta_{LAB}$ Range" , "MuID" , "Iso" , "Trig" , "Total" });
  if (useEtaCM) { colTitle1[0] = "$\\eta_{CM}$ Range"; }
  //
  texTable.push_back("\\begin{table}[h!]");
  texTable.push_back("  \\centering");
  //texTable.push_back("  \\resizebox{\\textwidth}{!}{");
  createUnc1DTable(texTable, colTyp, colCor, colTitle1, colChg, effMap, uncMap, col);
  //texTable.push_back("  }");
  texTable.push_back(Form("  \\caption{%s}",
                          Form("Statistical uncertainties corresponding to the muon efficiency corrections derived by applying the Tag and Probe scale factors event by event. Each source  for each source  The errors are shown as a function of the generated %s, derived from the %s %s \\POWHEG samples separated in negative and positive charged muons. %sGenerated muons are require to match reconstructed muons passing all analysis cuts.",
                               (useEtaCM ? "muon $\\eta_{CM}$" : "$\\eta_{LAB}$"),
                               (chg=="Plus" ? "$\\WToMuNuPl$" : "$\\WToMuNuMi$"),
                               col.c_str(),
                               ( (col=="PA") ? "The \\pPb and \\Pbp MC samples are combined as described in \\sect{sec:CombiningBeamDirection}. " : "")
                               )
                          )
                     );
  texTable.push_back(Form("  \\label{tab:tnpStatUncertainty_WToMu_%s_%s}", chg.c_str(), col.c_str()));
  texTable.push_back("\\end{table}");
  //
  for (const auto& row : texTable) { file << row << std::endl; }
  file << std::endl; file << std::endl;



};


bool printEff1DTables(const EffMap_t& effMap, const Unc1DMap_t& uncMap, const std::string& outDir, const bool isHFReweight)
{
  //
  std::cout << "[INFO] Filling the efficiency tables" << std::endl;
  //
  for (const auto& c : effMap.at("EtaCM").at("MC_WToMuNu")) {
    // Create Output Directory
    const std::string tableDir = outDir + "/Tables/Efficiency/" + c.first;
    makeDir(tableDir);
    // Create Output Files for MC Efficiency
    const std::string fileName_MC = "mcEfficiency_" + c.first;
    std::ofstream file_MC((tableDir + "/" + fileName_MC + ".tex").c_str());
    if (file_MC.is_open()==false) { std::cout << "[ERROR] File " << fileName_MC << " was not created!" << std::endl; return false; }
    //
    makeMCEff1DTable(file_MC, effMap, uncMap, c.first, isHFReweight);
    //
    // Create Output Files for TnP Efficiency
    const std::string fileName_TnPEff = "corrEfficiency_" + c.first;
    std::ofstream file_TnPEff((tableDir + "/" + fileName_TnPEff + ".tex").c_str());
    if (file_TnPEff.is_open()==false) { std::cout << "[ERROR] File " << fileName_TnPEff << " was not created!" << std::endl; return false; }
    //
    makeCorrEff1DTable(file_TnPEff, effMap, uncMap, c.first);
    //
    // Create Output Files for TnP Uncertainties
    const std::string fileName_TnPUnc = "tnpUncertainty_" + c.first;
    std::ofstream file_TnPUnc((tableDir + "/" + fileName_TnPUnc + ".tex").c_str());
    if (file_TnPUnc.is_open()==false) { std::cout << "[ERROR] File " << fileName_TnPUnc << " was not created!" << std::endl; return false; }
    //
    makeTnPUnc1DTable(file_TnPUnc, effMap, uncMap, c.first, "Plus");
    //
    makeTnPUnc1DTable(file_TnPUnc, effMap, uncMap, c.first, "Minus");
  }
  //
  return true;
};

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
#include "TVector.h"
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
using Unc1DVec_t   =  std::map< std::string , TVector >;
using Unc1DMap_t   =  std::map< std::string , std::map< std::string , std::map< std::string , std::map< std::string , std::map< std::string , Unc1DVec_t > > > > >;
using EffVec_t     =  std::map< std::string , std::vector< TEfficiency > >;
using EffMap_t     =  std::map< std::string , std::map< std::string , std::map< std::string , std::map< std::string , std::map< std::string , EffVec_t > > > > >;
using BinMap_t     =  std::map< std::string , std::vector< double > >;
using KeyPVec_t    =  std::vector< TKey* >;



// ------------------ FUNCTION -------------------------------
void         makeDir             ( const std::string& dir );
bool         existDir            ( const std::string& dir );
bool         getObjectsFromFile  ( EffMap_t& eff , Unc1DMap_t& unc , const std::string& filePath );
void         formatEff1D         ( TGraphAsymmErrors& graph , const std::string& var , const std::string& charge , const std::string& type );
void         formatEff1D         ( TEfficiency& eff , const std::string& var , const std::string& charge , const std::string& type );
void         drawEff1D           ( const std::string& outDir , EffMap_t& effMap , Unc1DMap_t unc );
bool         printTableEff1D     ( const std::string& outDir , const EffMap_t& effMap , const Unc1DMap_t& uncMap );
void         makeTable           ( std::ofstream& file , const std::vector< TEfficiency >& eff ,
                                   const std::vector< TVector >& tnpStat , const std::vector< TVector >& tnpSyst , 
                                   const std::string& type , const std::vector<std::string>& info , const bool& isAcceptance );
void         formatLegendEntry   ( TLegendEntry& e );
void         formatDecayLabel    ( std::string& label , const std::string& inLabel, const std::string c="#" );
void         setStyle            ( );
KeyPVec_t    getK                ( TList* list );
const char*  sgn                 ( const double n ) { if (n >= 0.) { return "+"; } else { return ""; } }


// ------------------ GLOBAL ------------------------------- 
// Kinematic Info
const BinMap_t  MU_BIN_RANGE_ = {
  { "Eta" , { -2.5 , 2.5 } },
  { "Pt"  , { 25., 200. } }
};
//
// Muon Charge
const std::vector< std::string > CHG_  = { "Plus" , "Minus" };

void printEfficiency(void)
{
  // Change the working directory
  const std::string CWD = getcwd(NULL, 0);
  const std::string mainDir = Form("%s/TnPEfficiency/", CWD.c_str());
  gSystem->mkdir(mainDir.c_str(), kTRUE);
  gSystem->ChangeDirectory(mainDir.c_str());
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
  // Set Style
  //setStyle();
  // Draw the Efficiencies
  //drawEff1D(mainDir, eff1D, unc1D);
  //
  // ------------------------------------------------------------------------------------------------------------------------
  //
  // Create Tables
  if (!printTableEff1D(mainDir, eff1D, unc1D)) { return; };
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
                    if (key->ReadObj()->IsA()->InheritsFrom(TVector::Class())) {
                      unc[v->GetName()][s->GetName()][c->GetName()][ch][t->GetName()][co->GetName()].ResizeTo( ((TVector*)key->ReadObj())->GetNrows() );
                      unc.at(v->GetName()).at(s->GetName()).at(c->GetName()).at(ch).at(t->GetName()).at(co->GetName()) = *((TVector*)key->ReadObj());
                    }
                  }
                }
                eff[v->GetName()][s->GetName()][c->GetName()][ch][t->GetName()]["NoCorrOnly"] = eff.at(v->GetName()).at(s->GetName()).at(c->GetName()).at(ch).at(t->GetName()).at("NoCorr");
                std::string name = eff.at(v->GetName()).at(s->GetName()).at(c->GetName()).at(ch).at(t->GetName()).at("NoCorr")[0].GetName();
                name.replace(name.find("NoCorr"), name.find("NoCorr")+6, "NoCorrOnly");
                eff.at(v->GetName()).at(s->GetName()).at(c->GetName()).at(ch).at(t->GetName()).at("NoCorrOnly")[0].SetName(name.c_str());
              }
              typeDir->cd();
            }
            if (co->ReadObj()->IsA()->InheritsFrom(TVector::Class())) {
              for (auto& ch : CHG_) {
                if (std::string(co->GetName()).find(ch)!=std::string::npos) {
                  std::string corrName = "";
                  if (std::string(co->GetName()).find("TnP_Stat")!=std::string::npos) { corrName = "TnP_Stat"; }
                  if (std::string(co->GetName()).find("TnP_Syst")!=std::string::npos) { corrName = "TnP_Syst"; }
                  if (std::string(co->GetName()).find("TnP_Tot" )!=std::string::npos) { corrName = "TnP_Tot";  }
                  if (corrName!="") {
                    unc[v->GetName()][s->GetName()][c->GetName()][ch][t->GetName()][corrName].ResizeTo( ((TVector*)co->ReadObj())->GetNrows() );
                    unc.at(v->GetName()).at(s->GetName()).at(c->GetName()).at(ch).at(t->GetName()).at(corrName) = *((TVector*)co->ReadObj());
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


void formatEff1D(TEfficiency& eff, const std::string& var, const std::string& charge, const std::string& type)
{
  // Set the format of all graphs
  if (eff.GetDimension() != 1) { std::cout << "[ERROR] formatEff1D only works for dimension == 1" << std::endl; return; }
  gPad->Update();
  auto graph = eff.GetPaintedGraph();
  const std::string objLbl = ( (type=="Total" || type=="Acceptance") ? "Gen #mu" : "Reco #mu" );
  // X-axis
  std::string xLabel = objLbl; if (charge=="Plus") xLabel += "^{+}"; if (charge=="Minus") xLabel += "^{-}";
  if (var=="Eta") { xLabel += " #eta"; }
  if (var=="Pt" ) { xLabel += " p_{T} (GeV/c)"; }
  // Y-axis
  std::string yLabel = type;
  if (type!="Acceptance") { type + " Efficiency"; }
  // Set Axis Titles
  eff.SetTitle(Form(";%s;%s", xLabel.c_str(), yLabel.c_str()));
  if (graph) { formatEff1D(*graph, var, charge, type); }
  gPad->Update(); 
};


void formatEff1D(TGraphAsymmErrors& graph, const std::string& var, const std::string& charge, const std::string& type)
{
  // General
  graph.SetMarkerColor(kBlue);
  graph.SetMarkerStyle(20);
  graph.SetMarkerSize(1.0);
  graph.SetFillStyle(1001);
  const std::string objLbl = ( (type=="Total" || type=="Acceptance") ? "Gen #mu" : "Reco #mu" );
  // X-axis
  std::string xLabel = objLbl; if (charge=="Plus") xLabel += "^{+}"; if (charge=="Minus") xLabel += "^{-}";
  if (var=="Eta") { xLabel += " #eta"; }
  if (var=="Pt" ) { xLabel += " p_{T} (GeV/c)"; }
  graph.GetXaxis()->CenterTitle(kFALSE);
  graph.GetXaxis()->SetTitleOffset(0.9);
  graph.GetXaxis()->SetTitleSize(0.050);
  graph.GetXaxis()->SetLabelSize(0.035);
  graph.GetXaxis()->SetLimits(MU_BIN_RANGE_.at(var)[0] , MU_BIN_RANGE_.at(var)[MU_BIN_RANGE_.at(var).size()-1]);
  // Y-axis
  std::string yLabel = type;
  if (type!="Acceptance") { type + " Efficiency"; }
  graph.GetYaxis()->CenterTitle(kFALSE);
  graph.GetYaxis()->SetTitleOffset(1.4);
  graph.GetYaxis()->SetTitleSize(0.04);
  graph.GetYaxis()->SetLabelSize(0.035);
  graph.GetYaxis()->SetRangeUser(0.79, 1.13);
  // Set Axis Titles
  graph.SetTitle(Form(";%s;%s", xLabel.c_str(), yLabel.c_str()));
};


void drawEff1D(const std::string& outDir, EffMap_t& effMap, Unc1DMap_t uncMap)
{
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
              Unc1DVec_t&       unc = uncMap.at(v.first).at(s.first).at(cl.first).at(ch.first).at(t.first);
              // Create Canvas
              TCanvas c("c", "c", 1000, 1000); c.cd();
              // Create the Text Info
              TLatex tex; tex.SetNDC(); tex.SetTextSize(0.028); float dy = 0;
              std::vector< std::string > textToPrint;
              std::string sampleLabel; formatDecayLabel(sampleLabel, sample);
              textToPrint.push_back(sampleLabel);
              if (var=="Eta") textToPrint.push_back("p^{#mu}_{T} > 25 GeV/c");
              if (var=="Pt" ) textToPrint.push_back("|#eta^{#mu}| < 2.4");
              // Declare the graph vector (for drawing with markers)
              std::vector< TGraphAsymmErrors > grVec;
              TLegend leg(0.2, 0.66, 0.4, 0.79);
              // Draw graph
              if (corr=="NoCorrOnly") {
                // Extract the Uncorrected Efficiency graph
                eff[0].Draw(); gPad->Update();
                grVec.push_back(*(eff[0].GetPaintedGraph()));
                // Format the Graphs
                for (auto& g : grVec) { formatEff1D(g, var, charge, type); }
                // Create Legend
                formatLegendEntry(*leg.AddEntry(&grVec[0], "MC Truth Efficiency", "pe"));
                // Draw the graph
                grVec[0].Draw("ap");
              }
              else if (corr=="NoCorr" && type=="Acceptance") {
                // Extract the Uncorrected Efficiency graph
                eff[0].Draw(); gPad->Update();
                grVec.push_back(*(eff[0].GetPaintedGraph()));
                // Format the Graphs
                for (auto& g : grVec) { formatEff1D(g, var, charge, type); }
                grVec[0].SetMarkerColor(kRed);
                // Create Legend
                formatLegendEntry(*leg.AddEntry(&grVec[0], "MC Acceptance", "pe"));
                // Draw the graph
                grVec[0].Draw("ap");
              }
              else if (corr=="NoCorr") {
                // Extract the Uncorrected Efficiency graph
                eff[0].Draw(); gPad->Update();
                grVec.push_back(*(eff[0].GetPaintedGraph()));
                // Extract the Corrected Efficiency graph
                t.second.at("TnP_Nominal")[0].Draw(); gPad->Update();
                grVec.push_back(*(t.second.at("TnP_Nominal")[0].GetPaintedGraph()));
                // Fill Corrected Efficiency graph with total TnP Uncertainties
                for (int l = 0; l < grVec[1].GetN(); l++) {
                  const double errorYlow  = std::sqrt( (grVec[1].GetErrorYlow(l)*grVec[1].GetErrorYlow(l)) + (unc.at("TnP_Tot")[l]*unc.at("TnP_Tot")[l]) );
                  const double errorYhigh = std::sqrt( (grVec[1].GetErrorYhigh(l)*grVec[1].GetErrorYhigh(l)) + (unc.at("TnP_Tot")[l]*unc.at("TnP_Tot")[l]) );
                  grVec[1].SetPointError(l, grVec[1].GetErrorXlow(l), grVec[1].GetErrorXhigh(l), errorYlow, errorYhigh);
                }
                // Format the Graphs
                for (auto& g : grVec) { formatEff1D(g, var, charge, type); }
                grVec[0].SetMarkerColor(kRed);
                // Create Legend
                formatLegendEntry(*leg.AddEntry(&grVec[0], "MC Truth Efficiency", "pe"));
                formatLegendEntry(*leg.AddEntry(&grVec[1], "Corrected Efficiency", "pe"));
                // Draw the graph
                grVec[0].Draw("ap");
                grVec[1].Draw("samep");
              }
              else if (corr=="TnP_Nominal") {
                // Extract the Corrected Efficiency graph
                eff[0].Draw(); gPad->Update();
                auto graph = eff[0].GetPaintedGraph();
                grVec.push_back(*graph);
                // Fill Corrected Efficiency graph with total TnP Uncertainties
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
                // Create Legend
                formatLegendEntry(*leg.AddEntry(&grVec[0], "Corrected Efficiency", "pe"));
                formatLegendEntry(*leg.AddEntry(&grVec[3], "MC Statistical Uncertainty", "f"));
                formatLegendEntry(*leg.AddEntry(&grVec[1], "TnP Statistical Uncertainty", "f"));
                formatLegendEntry(*leg.AddEntry(&grVec[2], "TnP Systematic Uncertainty", "f"));
                // Format the graphs
                for (auto& g : grVec) { formatEff1D(g, var, charge, type); }
                grVec[0].SetMarkerColor(kBlack);
                grVec[1].SetFillColor(kOrange);
                grVec[2].SetFillColor(kGreen+3);
                grVec[3].SetFillColor(kRed);
                // Draw the graphs
                grVec[2].Draw("a2");
                grVec[1].Draw("same2");
                grVec[3].Draw("same2");
                grVec[0].Draw("samep");
              }
              else {
                // Extract the Corrected Efficiency graph
                t.second.at("TnP_Nominal")[0].Draw(); gPad->Update();
                auto graph = t.second.at("TnP_Nominal")[0].GetPaintedGraph();
                grVec.push_back(*graph);
                // Fill Corrected Efficiency graph with the TnP Uncertainties
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
                for (auto& g : grVec) { formatEff1D(g, var, charge, type); }
                grVec[0].SetMarkerColor(kAzure+10);
                grVec[1].SetFillColor(kRed);
                for (uint n=2; n<grVec.size(); n++) { grVec[n].SetMarkerColor(kBlack);  grVec[n].SetMarkerSize(0.0); grVec[n].SetLineColor(kBlack); }
                // Create Legend
                formatLegendEntry(*leg.AddEntry(&grVec[0], "Corrected Efficiency", "pe"));
                formatLegendEntry(*leg.AddEntry(&grVec[1], Form("%s Uncertainty", corr.c_str()), "f"));
                formatLegendEntry(*leg.AddEntry(&grVec[2], Form("%s Variation", corr.c_str()), "l"));
                // Draw the Graph
                grVec[1].Draw("a2");
                for (uint n=2; n<grVec.size(); n++) { grVec[n].Draw("samep"); }
                grVec[1].Draw("same2");
                grVec[0].Draw("samep");
              }
              // Draw Legend and line
              TLine line(-2.4, 1.0, 2.4, 1.0); line.SetLineStyle(2);
              line.Draw("same");
              leg.Draw("same");
              // Update
              c.Modified(); c.Update();
              //
              if (corr=="NoCorr" || corr=="TnP_Nominal") {
                // Add Min, Max and Mean value of efficiency
                auto graph = eff[0].GetPaintedGraph();
                if (graph) {
                  const double min = TMath::MinElement(graph->GetN(), graph->GetY())*100.;
                  const double max = TMath::MaxElement(graph->GetN(), graph->GetY())*100.;
                  const double avg = graph->GetMean(2)*100.;
                  const double err = graph->GetRMS(2)*100.;
                  //textToPrint.push_back(Form("min: %.2f %% , max: %.2f %%", min, max));
                  //textToPrint.push_back(Form("mean: %.2f #pm %.2f %%", avg, err));
                }
              }
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


void makeTable(std::ofstream& file, const std::vector< TEfficiency >& eff, const std::vector< TVector >& tnpStat, const std::vector< TVector >& tnpSyst,
               const std::string& type, const std::vector<std::string>& info, const bool& isAcceptance)
{
  // Initialize the latex table
  std::vector< std::string > texTable;
  // 
  texTable.push_back("\\begin{table}[h!]");
  texTable.push_back("  \\centering");
  texTable.push_back("  \\begin{tabular}{*3c}");
  texTable.push_back("    \\hline");
  texTable.push_back("    \\multirow{2}{*}{Bin} & \\multicolumn{2}{c}{Efficiency}  \\\\");
  texTable.push_back("                          & $\\mu^{+}$      &     $\\mu^{-}$ \\\\[1ex]");
  texTable.push_back("    \\hline\\hline\\\\");
  //
  for (int i = 1; i <= eff[0].GetTotalHistogram()->GetNbinsX(); i++) {
    std::string bin;
    const double min = eff[0].GetTotalHistogram()->GetXaxis()->GetBinLowEdge(i);
    const double max = eff[0].GetTotalHistogram()->GetXaxis()->GetBinUpEdge(i);
    if (info[0]=="Pt" ) { bin = Form("$%s%.1f < \\pt < %s%.1f$", sgn(min), min, sgn(max) , max); }
    if (info[0]=="Eta") { bin = Form("$%s%.1f < \\eta < %s%.1f$", sgn(min), min, sgn(max) , max); }
    std::string eff_Pl, eff_Mi;
    if (type=="NoCorr") {
      eff_Pl = Form("$%.3f^{+%.3f}_{-%.3f}$", eff[0].GetEfficiency(i), eff[0].GetEfficiencyErrorUp(i), eff[0].GetEfficiencyErrorLow(i));
      eff_Mi = Form("$%.3f^{+%.3f}_{-%.3f}$", eff[1].GetEfficiency(i), eff[1].GetEfficiencyErrorUp(i), eff[1].GetEfficiencyErrorLow(i));
    }
    else {
      eff_Pl = Form("$%.3f^{+%.3f}_{-%.3f}", eff[0].GetEfficiency(i), eff[0].GetEfficiencyErrorUp(i), eff[0].GetEfficiencyErrorLow(i));
      eff_Pl += Form(" \\pm %.3f\\textrm{ (stat)} \\pm %.3f\\textrm{ (syst)}$", tnpStat[0][i-1], tnpSyst[0][i-1]);
      eff_Mi = Form("$%.3f^{+%.3f}_{-%.3f}", eff[1].GetEfficiency(i), eff[1].GetEfficiencyErrorUp(i), eff[1].GetEfficiencyErrorLow(i));
      eff_Mi += Form(" \\pm %.3f\\textrm{ (stat)} \\pm %.3f\\textrm{ (syst)}$", tnpStat[1][i-1], tnpSyst[1][i-1]);
    }
    texTable.push_back(("    " + bin + " & " + eff_Pl + " & " + eff_Mi + " \\\\[1ex]").c_str());
  }
  //
  texTable.push_back("    \\hline");
  texTable.push_back("  \\end{tabular}");
  std::string decay; formatDecayLabel(decay, info[1], "\\");
  std::string colInfo = " The \\pPb and \\Pbp MC samples are combined as described in \\sect{sec:CombiningBeamDirection}.";
  if (info[2]!="PA") { colInfo = ""; }
  std::string colType = (" "+info[2]);
  if (info[2]=="PA") { colType = ""; }
  std::string varInfo = ""; if (info[0]=="Pt") { varInfo = "\\pt"; }; if (info[0]=="Eta") { varInfo = "\\eta"; }
  std::string mcInfo = ""; if (type=="NoCorr") { mcInfo = "truth efficiency"; }; if (type=="TnP_Nominal") { mcInfo = "TnP corrected efficiency"; }
  std::string tnpInfo = " The labels \"stat\" and \"syst\" represent the TnP statisitical and systematic uncertainties, respectively.";
  if (type!="TnP_Nominal") { tnpInfo = ""; }
  if (isAcceptance) { texTable.push_back(Form("  \\label{tab:Acceptance_%s_%s_%s}", info[0].c_str(), info[1].c_str(), info[2].c_str())); }
  else { texTable.push_back(Form("  \\label{tab:Efficiency_%s_%s_%s_%s}", info[0].c_str(), info[1].c_str(), info[2].c_str(), type.c_str())); }
  if (isAcceptance) {
    texTable.push_back(Form("  \\caption{\\PW\\ acceptance as a function of the generated muon $%s$, derived from the $%s$%s \\POWHEG sample separated in negative and positive charged muons.%s Results corresponds to \\eq{eq:WAcceptance}.}", varInfo.c_str(), decay.c_str(), colType.c_str(), colInfo.c_str()));
  }
  else {
    texTable.push_back(Form("  \\caption{Muon %s as a function of the generated muon $%s$, derived from the $%s$%s \\POWHEG sample separated in negative and positive charged muons.%s%s Generated muons are require to match reconstructed muons passing all analysis cuts. Results corresponds to \\eq{eq:MCTruthEfficiency}.}", mcInfo.c_str(), varInfo.c_str(), decay.c_str(), colType.c_str(), tnpInfo.c_str(), colInfo.c_str()));
  }
  texTable.push_back("\\end{table}");
  //
  for (const auto& row : texTable) { file << row << std::endl; }
  file << endl; file << endl;
  //
};


bool printTableEff1D(const std::string& outDir, const EffMap_t& effMap, const Unc1DMap_t& uncMap)
{
  // Draw all graphs
  for (auto& v : effMap) {
    for (auto& s : v.second) {
      for (auto& cl : s.second) {
        // Create Output Directory
        const std::string tableDir = outDir + "TablesLATEX/" + v.first+"/" + s.first+"/" + cl.first;
        makeDir(tableDir);
        // Create Output Files
        const std::string fileName_NoCorr = "eff1D_" + v.first +"_"+ s.first +"_"+ cl.first +"_"+ "Total" +"_"+ "NoCorr";
        std::ofstream file_NoCorr((tableDir + "/" + fileName_NoCorr + ".tex").c_str());
        if (file_NoCorr.is_open()==false) { std::cout << "[ERROR] File " << fileName_NoCorr << " was not created!" << std::endl; return false; }
        const std::string fileName_Corr   = "eff1D_" + v.first +"_"+ s.first +"_"+ cl.first +"_"+ "Total" +"_"+ "TnP_Nominal";
        std::ofstream file_Corr((tableDir + "/" + fileName_Corr + ".tex").c_str());
        if (file_Corr.is_open()==false) { std::cout << "[ERROR] File " << fileName_Corr << " was not created!" << std::endl; return false; }
        // Efficiency
        const TEfficiency& eff_Pl_NoCorr = effMap.at(v.first).at(s.first).at(cl.first).at("Plus").at("Total").at("NoCorr")[0];
        const TEfficiency& eff_Mi_NoCorr = effMap.at(v.first).at(s.first).at(cl.first).at("Minus").at("Total").at("NoCorr")[0];
        const TEfficiency& eff_Pl_Corr   = effMap.at(v.first).at(s.first).at(cl.first).at("Plus").at("Total").at("TnP_Nominal")[0];
        const TEfficiency& eff_Mi_Corr   = effMap.at(v.first).at(s.first).at(cl.first).at("Minus").at("Total").at("TnP_Nominal")[0];
        // Uncertainty
        const TVector&     unc_Pl_Stat   = uncMap.at(v.first).at(s.first).at(cl.first).at("Plus").at("Total").at("TnP_Stat");
        const TVector&     unc_Mi_Stat   = uncMap.at(v.first).at(s.first).at(cl.first).at("Minus").at("Total").at("TnP_Stat");
        const TVector&     unc_Pl_Syst   = uncMap.at(v.first).at(s.first).at(cl.first).at("Plus").at("Total").at("TnP_Syst");
        const TVector&     unc_Mi_Syst   = uncMap.at(v.first).at(s.first).at(cl.first).at("Minus").at("Total").at("TnP_Syst");
        // Make table for uncorrected efficiency
        makeTable(file_NoCorr, {eff_Pl_NoCorr, eff_Mi_NoCorr}, {}, {}, "NoCorr", {v.first, s.first, cl.first}, false);
        // Make table for corrected efficiency
        makeTable(file_Corr, {eff_Pl_Corr, eff_Mi_Corr}, {unc_Pl_Stat, unc_Mi_Stat}, {unc_Pl_Syst, unc_Mi_Syst}, "TnP_Nominal", {v.first, s.first, cl.first}, false);
        // Make table with all the uncertainties splitted
        //makeTable(file_Corr, uncMap.at(v.first).at(s.first).at(cl.first).at("Plus").at("Total"), "Plus", {v.first, s.first, cl.first});
        //makeTable(file_Corr, uncMap.at(v.first).at(s.first).at(cl.first).at("Minus").at("Total"), "Minus", {v.first, s.first, cl.first});
        // Acceptance
        const TEfficiency& acc_Pl_NoCorr = effMap.at(v.first).at(s.first).at(cl.first).at("Plus").at("Acceptance").at("NoCorr")[0];
        const TEfficiency& acc_Mi_NoCorr = effMap.at(v.first).at(s.first).at(cl.first).at("Minus").at("Acceptance").at("NoCorr")[0];
        // Make table for uncorrected efficiency
        makeTable(file_NoCorr, {acc_Pl_NoCorr, acc_Mi_NoCorr}, {}, {}, "NoCorr", {v.first, s.first, cl.first}, true);
      }
    }
  }
  return true;
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

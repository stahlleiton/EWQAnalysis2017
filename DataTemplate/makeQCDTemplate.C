// Auxiliary Headers
#include "../Utilities/EVENTUTILS.h"
// ROOT headers
#include "TSystem.h"
#include "TSystemFile.h"
#include "TSystemDirectory.h"
#include "TROOT.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TLegendEntry.h"
// RooFit headers
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooAbsData.h"
// c++ headers
#include <dirent.h>
#include <iostream>
#include <fstream>
#include <map>
#include <tuple>
#include <vector>
// CMS headers
#include "../Utilities/CMS/tdrstyle.C"
#include "../Utilities/CMS/CMS_lumi.C"


// Type definitions
using RooRealVarMap_t = std::map<std::string, RooRealVar>;
using Content_t    = std::tuple< std::string , RooRealVarMap_t , RooRealVarMap_t >;
using Par_t        = std::pair< double , double >;
using Bin_t        = std::tuple< Par_t , Par_t >;
using Val_t        = std::tuple< Par_t , Par_t >;
using ContentMap_t = std::map< std::string , std::map< std::string , std::map< Bin_t , std::vector< Val_t > > > >;
using GraphMap_t   = std::map< std::string , std::map< std::string , std::map< Bin_t , TGraphErrors* > > >;
using FitMap_t     = std::map< std::string , std::map< std::string , std::map< Bin_t , TF1* > > >;
using FitMap2_t    = std::map< std::string , std::map< Bin_t , std::map< std::string , TF1* > > >;
using BinMap_t     =  std::map< std::string , std::vector< double > >;
using BinMapMap_t  =  std::map< std::string , BinMap_t >;

// Function declaration
bool     findFiles   ( const std::string&  dirPath , std::vector<std::string>& fileNames );
bool     readFiles   ( const std::vector<std::string>& fileNames, ContentMap_t& contentMap , bool& useEtaCM );
bool     getInfo     ( const std::string&  fileName , Content_t& content , bool& useEtaCM );
bool     getBinPars  ( const RooWorkspace& myws , RooRealVarMap_t& binPars );
bool     getFitPars  ( const RooWorkspace& myws , RooRealVarMap_t& fitPars , std::string& label );
void     findBinLimits ( BinMapMap_t& binMap , const ContentMap_t& contentMap , const bool& useEtaCM );
bool     fillGraph   ( const ContentMap_t& contentMap , GraphMap_t& graphMap );
void     formatGraph ( GraphMap_t& graphMap );
void     drawGraph   ( const std::string& outDir , const GraphMap_t& graphMap , const FitMap_t& fitMap , const bool& useEtaCM );
void     clearGraph  ( GraphMap_t& graphMap );
bool     fitGraph    ( const GraphMap_t& graphMap , FitMap_t& fitMap );
bool     computeMeanAndVariance ( FitMap_t& fitMap , const std::string& col , const bool& useEtaCM , const BinMapMap_t& binMap );
void     formatFit   ( FitMap_t& fitMap );
void     clearFit    ( FitMap_t& fitMap );
bool     fillExGraph ( const FitMap_t& fitMap , GraphMap_t& exGraphMap );
void     formatExGraph ( GraphMap_t& exGraphMap , const std::string& col , const bool& useEtaCM , const BinMapMap_t& binMap );
void     drawExGraph ( const std::string& outDir , const GraphMap_t& exGraphMap , const bool& useEtaCM );
void     getBinText  ( const Bin_t& bin , std::vector<std::string>& text , const int& chg=-1 , const bool& useEtaCM = true );
void     getFitText  ( const TF1& fit, std::vector<std::string>& text );
void     formatLegendEntry ( TLegendEntry& e );
void     updateYAxisRange ( TGraphErrors* graph , TF1* fit );
double   roundVal    ( const double& v, const uint& n=3 ) { return ( floor(v * pow(10, n)) / pow(10, n) ); }
Double_t pol1        ( Double_t *x , Double_t *par ) { return ( par[0] + par[1]*(x[0] - par[2]) ); }
bool     createInputFiles ( const std::string& dirPath , const FitMap_t& iniFitMap , const bool& useEtaCM , const BinMapMap_t& binMap );

// Global 
const std::vector< std::string > cutSelection = {
  "p^{#mu}_{T}>25 GeV/c, |#eta^{#mu}|<2.4, Iso^{#mu}<0.15", "#mu Tight ID , Drell-Yan Veto"
};
std::map< std::string, std::tuple< std::string > > yAxis;
const double isoPoint = 0.03;
const std::string useBin = "Inclusive"; // Can also  be Inclusive
std::string MODELNAME_ = "";

void makeQCDTemplate(
                     const std::string workDirName = "QCDTemplateCM",
                     const std::string METTag      = "METPF_RAW",
                     const std::string DSTag       = "DATA"
                     )
{
  //
  if (workDirName.find("QCDTemplate")==std::string::npos) { std::cout << "[ERROR] Invalid input workdirname " << workDirName << std::endl; return; }
  const std::string CWD = getcwd(NULL, 0);
  std::string preCWD = CWD; preCWD.erase(preCWD.find_last_of("/"), 100);
  //
  const std::vector< std::string > col = { "pPb" , "Pbp" , "PA" };
  for (const auto& c : col) {
    const std::string dirPath = Form("%s/Fitter/Output/%s/%s/%s/%s/result/", preCWD.c_str(), workDirName.c_str(), METTag.c_str(), DSTag.c_str(), c.c_str());
    const std::string outDir  = Form("%s/Output/%s/%s/%s/%s/", CWD.c_str(), workDirName.c_str(), METTag.c_str(), DSTag.c_str(), c.c_str());
    // Find all the workspaces
    std::vector<std::string> fileNames;
    if (!findFiles(dirPath, fileNames)) continue;
    // For checking if eta is at CM
    bool useEtaCM = false;
    // Extract the information from each file
    ContentMap_t contentMap;
    if (!readFiles(fileNames, contentMap, useEtaCM)) { return; }
    // Find Bin Boundaries
    BinMapMap_t  MU_BIN_RANGE;
    findBinLimits(MU_BIN_RANGE, contentMap, useEtaCM);
    // Fill the graphs
    GraphMap_t graphMap;
    if (!fillGraph(contentMap, graphMap)) { return; }
    formatGraph(graphMap);
    // Fit the graphs
    FitMap_t fitMap;
    if (!fitGraph(graphMap, fitMap)) { return; }
    computeMeanAndVariance(fitMap, c, useEtaCM, MU_BIN_RANGE);
    formatFit(fitMap);
    // Create the graphs with extrapolated variables
    GraphMap_t exGraphMap;
    if (!fillExGraph(fitMap, exGraphMap)) { return; }
    formatExGraph(exGraphMap, c, useEtaCM, MU_BIN_RANGE);
    // Set the CMS style
    setTDRStyle();
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    // Draw the graphs
    drawGraph(outDir, graphMap, fitMap, useEtaCM);
    // Delete the graphs
    clearGraph(graphMap);
    // Draw the extrapolated graphs
    drawExGraph(outDir, exGraphMap, useEtaCM);
    // Delete the extrapolated graphs
    clearGraph(exGraphMap);
    // Create new input files
    createInputFiles((outDir+"InputFiles/"), fitMap, useEtaCM, MU_BIN_RANGE);
    // Create the fits
    clearFit(fitMap);
  }
};


bool findFiles(const std::string& dirPath, std::vector<std::string>& fileNames)
{
  DIR *dpdf; struct dirent *epdf;
  // Open the working directory
  dpdf = opendir(dirPath.c_str());
  // Search for all the files inside this directory
  if (dpdf != NULL){
    while ( (epdf=readdir(dpdf)) ){
      if (strcmp(epdf->d_name,".")!=0 && strcmp(epdf->d_name,"..")!=0 ) {
        std::cout << "[INFO] Adding file: " << epdf->d_name << std::endl;
        fileNames.push_back(dirPath+std::string(epdf->d_name));
      }
    }
  } else {
    std::cout << "[ERROR] Working directory was not found!" << std::endl; return false;
  }
  return true;
};


bool readFiles(const std::vector<std::string>& fileNames, ContentMap_t& contentMap, bool& useEtaCM)
{
  if (fileNames.size()==0) { return false; }
  for (const auto& fileName : fileNames) {
    // Extract the content for the file
    Content_t content;
    if (!getInfo(fileName, content, useEtaCM)) { return false; }
    // Insert the content in the map
    std::string label = std::get<0>(content);
    RooRealVarMap_t binPars = std::get<1>(content);
    RooRealVarMap_t fitPars = std::get<2>(content);
    // Make bin
    bool ispPb = false; if (fileName.find("_PA_")!=std::string::npos || fileName.find("_pPb_")!=std::string::npos) { ispPb = true; }
    double etaMin = binPars.at("Muon_Eta").getMin(); if (useEtaCM) { etaMin = PA::EtaLABtoCM(etaMin, ispPb); }
    double etaMax = binPars.at("Muon_Eta").getMax(); if (useEtaCM) { etaMax = PA::EtaLABtoCM(etaMax, ispPb); }
    Par_t eta = std::make_pair(roundVal(etaMin), roundVal(etaMax));
    Par_t pt  = std::make_pair(roundVal(binPars.at("Muon_Pt").getMin()), roundVal(binPars.at("Muon_Pt").getMax()));
    Bin_t bin = std::make_tuple(eta, pt);
    // Make x var
    Par_t xPar = std::make_pair(binPars.at("Muon_Iso").getVal(), binPars.at("Muon_Iso").getError());
    // Make y var
    for (const auto& par : fitPars) {
      std::string parName = par.first;
      RooRealVar  parVar  = par.second;
      Par_t       yPar    = std::make_pair(parVar.getVal(), parVar.getError());
      Val_t       value   = std::make_tuple(xPar, yPar);
      contentMap[parName][label][bin].push_back(value);
    }
  }
  if (contentMap.size()==0) { std::cout << "[ERROR] The contentMap is empty!" << std::endl; return false; }
  return true;
};


bool getInfo(const std::string& fileName, Content_t& content, bool& useEtaCM)
{
  // Extract the workspace
  TFile *f = TFile::Open( fileName.c_str() );
  if (!f) { std::cout << "[Error] " << fileName << " not found" << std::endl; return false; }
  RooWorkspace *ws = (RooWorkspace*) f->Get("workspace");
  if (!ws) { std::cout << "[ERROR] Workspace not found in " << fileName << std::endl; f->Close(); delete f; return false; }
  // Extract the muon parameters used for binning
  RooRealVarMap_t binPars;
  if (!getBinPars(*ws, binPars)) { f->Close(); delete f; return false; }
  // Extract the MET fitted parameters
  std::string label;
  RooRealVarMap_t fitPars;
  if(!getFitPars(*ws, fitPars, label)) { f->Close(); delete f; return false; }
  // Organize the fit and bin information
  content = std::make_tuple(label, binPars, fitPars);
  // Check if we are using eta at Center of Mass
  if (ws->var("useEtaCM")!=NULL && ws->var("useEtaCM")->getVal()==1.0) { useEtaCM = true; }
  // Clean up
  f->Close(); delete f;
  return true;
};


bool getBinPars(const RooWorkspace& myws, RooRealVarMap_t& binPars)
{
  // Extract the parameters used for binning
  if (myws.var("Muon_Eta")) { binPars["Muon_Eta"] = *myws.var("Muon_Eta"); }
  else { std::cout << "[ERROR] " << "Muon_Eta" << " was not found!" << std::endl; return false; }
  if (myws.var("Muon_Pt")) { binPars["Muon_Pt"] = *myws.var("Muon_Pt"); }
  else { std::cout << "[ERROR] " << "Muon_Pt" << " was not found!" << std::endl; return false; }
  if (myws.var("Muon_Iso")) { binPars["Muon_Iso"] = *myws.var("Muon_Iso"); }
  else { std::cout << "[ERROR] " << "Muon_Iso" << " was not found!" << std::endl; return false; }
  // Compute the mean and RMS values for each bin parameter
  const std::string DSTAG = ( myws.obj("DSTAG") ? ((TObjString*)myws.obj("DSTAG"))->GetString().Data() : "" ); 
  if (DSTAG=="") { std::cout << "[ERROR] DSTAG was not found!" << std::endl; return false; }
  const std::string charge  = ( myws.obj("fitCharge") ? ((TObjString*)myws.obj("fitCharge"))->GetString().Data() : "" );
  if (charge=="") { std::cout << "[ERROR] fitCharge was not found!" << std::endl; return false; };
  const std::string dsName = "d" + charge + "_" + DSTAG;
  if (myws.data(dsName.c_str())) {
    for (auto& par : binPars) {
      par.second.setVal(myws.data(dsName.c_str())->meanVar(*myws.var(par.first.c_str()))->getVal());
      par.second.setError(myws.data(dsName.c_str())->rmsVar(*myws.var(par.first.c_str()))->getVal());
    }
  }
  else { std::cout << "[ERROR] " << dsName << " was not found!" << std::endl; return false; }
  return true;
};
 
 
bool getFitPars(const RooWorkspace& myws, RooRealVarMap_t& fitPars, std::string& label)
{
  // Extract the internal info
  const std::string channel = ( myws.obj("channel") ? ((TObjString*)myws.obj("channel"))->GetString().Data() : "" ); 
  if (channel=="") { std::cout << "[ERROR] channel was not found!" << std::endl; return false; }
  const std::string object  = ( myws.obj("fitObject") ? ((TObjString*)myws.obj("fitObject"))->GetString().Data() : "" );
  if (object=="") { std::cout << "[ERROR] fitObject was not found!" << std::endl; return false; }
  const std::string charge  = ( myws.obj("fitCharge") ? ((TObjString*)myws.obj("fitCharge"))->GetString().Data() : "" );
  if (charge=="") { std::cout << "[ERROR] fitCharge was not found!" << std::endl; return false; };
  const std::string beamDir = ( myws.obj("fitObject") ? ((TObjString*)myws.obj("fitSystem"))->GetString().Data() : "" );
  if (beamDir=="") { std::cout << "[ERROR] fitSystem was not found!" << std::endl; return false; }
  label = object + channel + charge + "_" + beamDir;
  // Extract the model name
  const std::string modelName = ( myws.obj(("Model_"+label).c_str()) ? ((TObjString*)myws.obj(("Model_"+label).c_str()))->GetString().Data() : "" );
  if (modelName=="") { std::cout << "[ERROR] " << ("Model_"+label) << " was not found!" << std::endl; return false; }
  // Extract the parameters of the model
  if (modelName=="MultiJetBkg") {
    if (myws.var(("Alpha_"+label).c_str())) { fitPars["Alpha"] = *myws.var(("Alpha_"+label).c_str()); }
    else { std::cout << "[ERROR] " << ("Alpha_"+label).c_str() << " was not found!" << std::endl; return false; }
    if (myws.var(("Beta_"+label).c_str())) { fitPars["Beta"] = *myws.var(("Beta_"+label).c_str()); }
    else { std::cout << "[ERROR] " << ("Beta_"+label).c_str() << " was not found!" << std::endl; return false; }
    if (myws.var(("x0_"+label).c_str())) { fitPars["x0"] = *myws.var(("x0_"+label).c_str()); }
    else { std::cout << "[ERROR] " << ("x0_"+label).c_str() << " was not found!" << std::endl; return false; }
  }
  if (modelName=="ModifiedRayleigh") {
    if (myws.var(("Sigma0_"+label).c_str())) { fitPars["Sigma0"] = *myws.var(("Sigma0_"+label).c_str()); }
    else { std::cout << "[ERROR] " << ("Sigma0_"+label).c_str() << " was not found!" << std::endl; return false; }
    if (myws.var(("Sigma1_"+label).c_str())) { fitPars["Sigma1"] = *myws.var(("Sigma1_"+label).c_str()); }
    else { std::cout << "[ERROR] " << ("Sigma1_"+label).c_str() << " was not found!" << std::endl; return false; }
    if (myws.var(("Sigma2_"+label).c_str())) { fitPars["Sigma2"] = *myws.var(("Sigma2_"+label).c_str()); }
    else { std::cout << "[ERROR] " << ("Sigma2_"+label).c_str() << " was not found!" << std::endl; return false; }
  }
  MODELNAME_ = modelName;
  return true;
};


void findBinLimits(BinMapMap_t& binMap, const ContentMap_t& contentMap, const bool& useEtaCM)
{
  // Clean the input bin Map
  binMap.clear();
  // Find the bin boundaries
  for (const auto& lbl : contentMap.begin()->second) {
    std::string col = "";
    if (lbl.first.find("pPb")!=std::string::npos) { col = "pPb"; }
    if (lbl.first.find("Pbp")!=std::string::npos) { col = "Pbp"; }
    if (lbl.first.find("PA" )!=std::string::npos) { col = "PA";  }
    double minEta = 999999. , maxEta = -999999. , minPt = 999999. , maxPt = -999999.;
    for (const auto& bin : lbl.second) {
      const double etaMin = std::get<0>(bin.first).first;
      const double etaMax = std::get<0>(bin.first).second;
      const double ptMin  = std::get<1>(bin.first).first;
      const double ptMax  = std::get<1>(bin.first).second;
      if (minEta > etaMin) { minEta = etaMin; }
      if (maxEta < etaMax) { maxEta = etaMax; }
      if (minPt  > ptMin ) { minPt  = ptMin;  }
      if (maxPt  < ptMax ) { maxPt  = ptMax;  }
    }
    binMap[col]["Pt"]  = { minPt , maxPt };
    if (useEtaCM) { binMap[col]["EtaCM_Inc"] = { minEta , maxEta }; }
    else          { binMap[col]["Eta_Inc"  ] = { minEta , maxEta }; }
    minEta = (std::floor((minEta-0.1)*10.0)/10.0);
    maxEta = (std::ceil ((maxEta+0.1)*10.0)/10.0);
    if (useEtaCM) { binMap[col]["EtaCM"] = { minEta , maxEta }; }
    else          { binMap[col]["Eta"  ] = { minEta , maxEta }; }
  };
};


bool fillGraph(const ContentMap_t& contentMap, GraphMap_t& graphMap)
{
  // Clean the input graph Map
  clearGraph(graphMap);
  // Fill the new graphs
  for (const auto& var : contentMap) {
    for (const auto& lbl : var.second) {
      for (const auto& bin : lbl.second) {
        const uint n = bin.second.size();
        Double_t x[n] , y[n] , ex[n] , ey[n];
        for (uint i = 0; i < n; i++) {
          x[i]  = std::get<0>(bin.second[i]).first;
          ex[i] = std::get<0>(bin.second[i]).second;
          y[i]  = std::get<1>(bin.second[i]).first;
          ey[i] = std::get<1>(bin.second[i]).second;
        }
        graphMap[var.first][lbl.first][bin.first] = new TGraphErrors(n, x, y, ex, ey);
        if (graphMap.at(var.first).at(lbl.first).at(bin.first)==NULL) {
          std::cout << "[ERROR] Element in graphMap is NULL!" << std::endl; return false;
        }
        std::string name = Form("graph_%s_%s_%.0f_eta_%.0f_%.0f_pt_%.0f",
                                var.first.c_str(),
                                lbl.first.c_str(), 
                                std::get<0>(bin.first).first*10., std::get<0>(bin.first).second*10.,
                                std::get<1>(bin.first).first*10., std::get<1>(bin.first).second*10.);
        graphMap.at(var.first).at(lbl.first).at(bin.first)->SetName(name.c_str());
      }
    }
  }
  return true;
};


void formatGraph(GraphMap_t& graphMap)
{
  yAxis["Alpha"]  = std::make_tuple( "#alpha"     );
  yAxis["Beta"]   = std::make_tuple( "#beta"      );
  yAxis["x0"]     = std::make_tuple( "x_{0}"      );
  yAxis["Sigma0"] = std::make_tuple( "#sigma_{0}" );
  yAxis["Sigma1"] = std::make_tuple( "#sigma_{1}" );
  yAxis["Sigma2"] = std::make_tuple( "#sigma_{2}" );
  // Set the format of all graphs
  for (auto& var : graphMap) {
    for (auto& lbl : var.second) {
      for (auto& graph : lbl.second) {
        if (graph.second) {
          // General
          graph.second->SetTitle("");
          graph.second->SetMarkerColor(kBlue);
          graph.second->SetMarkerStyle(20);
          // X-axis
          graph.second->GetXaxis()->SetTitle("Relative PF r=0.3 #beta=0 Isolation");
          graph.second->GetXaxis()->CenterTitle(kFALSE);
          graph.second->GetXaxis()->SetTitleOffset(1.2);
          graph.second->GetXaxis()->SetTitleSize(0.04);
          graph.second->GetXaxis()->SetLabelSize(0.035);
          //graph.second->GetXaxis()->SetNdivisions(204);
          graph.second->GetXaxis()->SetLimits(-0.1, roundVal(graph.second->GetXaxis()->GetXmax(),1));
          // Y-axis
          graph.second->GetYaxis()->SetTitle(std::get<0>(yAxis.at(var.first)).c_str());
          graph.second->GetYaxis()->CenterTitle(kFALSE);
          graph.second->GetYaxis()->SetTitleOffset(1.0);
          graph.second->GetYaxis()->SetTitleSize(0.04);
          graph.second->GetYaxis()->SetLabelSize(0.035);
          //graph.second->GetYaxis()->SetNdivisions(204);
          graph.second->GetYaxis()->SetRangeUser(-50., 50.);
        }
      }
    }
  }
};


void drawGraph(const std::string& outDir, const GraphMap_t& graphMap, const FitMap_t& fitMap, const bool& useEtaCM)
{
  // Make output directories
  makeDir(outDir+"Plot/Initial/png/");
  makeDir(outDir+"Plot/Initial/pdf/");
  makeDir(outDir+"Plot/Initial/root/");
  Double_t xl1=.52, yl1=0.74, xl2=xl1+.2, yl2=yl1+.150;
  // Set the format of all graphs
  for (auto& var : graphMap) {
    for (auto& lbl : var.second) {
      for (auto& graph : lbl.second) {
        if (graph.second) {
          TF1 *f = fitMap.at(var.first).at(lbl.first).at(graph.first);
          // Create Canvas
          TCanvas* c   = new TCanvas("c", "c", 1000, 1000); c->cd();
          // Create Legend
          TLegend *leg = new TLegend(xl1,yl1,xl2,yl2);
          formatLegendEntry(*leg->AddEntry(graph.second, "Data", "p"));
          if (f) formatLegendEntry(*leg->AddEntry(f, "Fit", "l"));
          // Create the Text Info
          TLatex *tex = new TLatex(); tex->SetNDC(); tex->SetTextSize(0.025); float dy = 0;
          std::vector< std::string > cutSelection;
          getBinText(graph.first, cutSelection, ((lbl.first.find("Pl")!=std::string::npos)*2-1), useEtaCM);
          std::vector< std::string > fitInfo;
          if (f) getFitText(*f, fitInfo);
          // Draw graph
          updateYAxisRange(graph.second, f);
          graph.second->Draw("ap");
          // Draw the fit
          TGraphErrors *tmp = NULL;
          if (f) {
            double x[]  = { f->GetParameter(2) } , y[]  = { f->GetParameter(0) };
            double ex[] = { f->GetParError(2)  } , ey[] = { f->GetParError(0)  };
            tmp = new TGraphErrors(1, x, y, ex, ey);
            tmp->Draw("samep"); tmp->SetMarkerColor(kRed); tmp->SetMarkerStyle(20); tmp->SetMarkerSize(1.5);
            f->Draw("same");
          }
          // Draw the text
          for (const auto& s: cutSelection) { tex->DrawLatex(0.20, 0.85-dy, s.c_str()); dy+=0.04; }
          if (f) for (const auto& s: fitInfo) { tex->DrawLatex(0.70, 0.84-dy, s.c_str()); dy+=0.04; }
          // Draw the legend
          leg->Draw("same");
          c->Modified(); c->Update();
          // set the CMS style
          int option = 111;
          if (lbl.first.find("pPb")!=std::string::npos) option = 109;
          if (lbl.first.find("Pbp")!=std::string::npos) option = 110;
          CMS_lumi(c, option, 33, "");
          c->Modified(); c->Update();
          // Save Canvas
          c->SaveAs(( outDir + "Plot/Initial/png/" + graph.second->GetName() + ".png" ).c_str());
          c->SaveAs(( outDir + "Plot/Initial/pdf/" + graph.second->GetName() + ".pdf" ).c_str());
          c->SaveAs(( outDir + "Plot/Initial/root/" + graph.second->GetName() + ".root" ).c_str());
          // Clean up memory
          c->Clear(); c->Close();
          delete c;
          delete leg;
          delete tex;
          if (tmp) delete tmp;
        }
      }
    }
  }
};


void updateYAxisRange(TGraphErrors* graph, TF1* fit)
{
  if (graph==NULL) return;
  // Find the max and min of graph
  double yMin = 9999999999. , yMax = -9999999999.;
  for (int i = 0; i < graph->GetN(); i++) {
    double x, y, yL, yH;
    graph->GetPoint(i, x, y);
    yL = y - graph->GetErrorY(i);
    yH = y + graph->GetErrorY(i);
    if (yH > yMax) { yMax = yH; }
    if (yL < yMin) { yMin = yL; }
  }
  // Find the max and min of TF1
  if (fit) {
    double fitMin = fit->GetMinimum(graph->GetXaxis()->GetXmin(), graph->GetXaxis()->GetXmax());
    double fitMax = fit->GetMaximum(graph->GetXaxis()->GetXmin(), graph->GetXaxis()->GetXmax());
    if (fitMax > yMax) yMax = fitMax;
    if (fitMin < yMin) yMin = fitMin;
  }
  // Set the up and down of y axis
  double fMax = 0.5 , fMin = 0.1;
  double yDown = (fMax*yMin - fMin*yMax)/(fMax - fMin);
  double yUp   = yDown + (yMax-yMin)/(fMax-fMin);
  // Update the y range
  graph->GetYaxis()->SetRangeUser(floor(yDown), ceil(yUp));
};


void getBinText(const Bin_t& bin, std::vector<std::string>& text, const int& chg, const bool& useEtaCM)
{
  // Clean the text
  text.clear();
  std::string sgn = (chg>0 ? "+" : "-");
  // Set Eta Bin Text
  if (std::get<0>(bin).first!=std::get<0>(bin).second) {
    if (useEtaCM) {
      text.push_back(Form("%.2f < #eta^{#mu%s}_{CM} < %.2f", std::get<0>(bin).first, sgn.c_str(), std::get<0>(bin).second));
    }
    else {
      text.push_back(Form("%.2f < #eta^{#mu%s}_{LAB} < %.2f", std::get<0>(bin).first, sgn.c_str(), std::get<0>(bin).second));
    }
  }
  // Set Pt Bin Text
  if (std::get<1>(bin).first!=std::get<1>(bin).second) {
    text.push_back(Form("p^{#mu%s}_{T} > %.1f GeV/c", sgn.c_str(), std::get<1>(bin).first));
  }
  // Set Quality Cuts Text
  text.push_back(Form("#mu%s Tight ID , Drell-Yan Veto", sgn.c_str()));
};


void formatLegendEntry(TLegendEntry& e)
{
  e.SetTextSize(0.035);
};


void clearGraph(GraphMap_t& graphMap)
{
  // Delete pointers to TGraph
  for (auto& var : graphMap) {
    for (auto& lbl : var.second) {
      for (auto& graph : lbl.second) {
        if (graph.second) delete graph.second; 
      }
    }
  }
  // Clear the graph map
  graphMap.clear();
};


bool fitGraph(const GraphMap_t& graphMap, FitMap_t& fitMap)
{
  // Clean the input fit Map
  clearFit(fitMap);
  // Set the extrapolation point
  double x0 = isoPoint;
  // Fit the graphs
  for (const auto& var : graphMap) {
    for (const auto& lbl : var.second) {
      for (const auto& graph : lbl.second) {
        std::string name = Form("fit_%s_%s_%.0f_eta_%.0f_%.0f_pt_%.0f",
                                var.first.c_str(),
                                lbl.first.c_str(), 
                                std::get<0>(graph.first).first*10., std::get<0>(graph.first).second*10.,
                                std::get<1>(graph.first).first*10., std::get<1>(graph.first).second*10.);
        fitMap[var.first][lbl.first][graph.first] = new TF1(name.c_str(), pol1, graph.second->GetXaxis()->GetXmin(), graph.second->GetXaxis()->GetXmax(), 3);
        if (fitMap.at(var.first).at(lbl.first).at(graph.first)==NULL) {
          std::cout << "[ERROR] Element in fitMap is NULL!" << std::endl; return false;
        }
        std::string v = std::get<0>(yAxis.at(var.first));
        fitMap.at(var.first).at(lbl.first).at(graph.first)->SetParNames(Form("%s", v.c_str()), "s", "Iso");
        // Estimate the initial parameters of the fit
        double x1, y1, x2, y2, ex1, ex2;
        graph.second->GetPoint(0, x1, y1);
        graph.second->GetPoint((graph.second->GetN()-1), x2, y2);
        ex1 = graph.second->GetErrorX(0);
        ex2 = graph.second->GetErrorX(1);
        double m  = ((y2 - y1)/(x2 - x1));
        double y0 = (y1 - m*(x1 - x0));
        fitMap.at(var.first).at(lbl.first).at(graph.first)->SetParameters(y0, m, x0);
        fitMap.at(var.first).at(lbl.first).at(graph.first)->SetParLimits(0, -50., 50.);
        fitMap.at(var.first).at(lbl.first).at(graph.first)->SetParLimits(1, -50., 50.);
        fitMap.at(var.first).at(lbl.first).at(graph.first)->FixParameter(2, x0);
        // Fit the graph
        graph.second->Fit(name.c_str(), "0QEFS", "", (floor((x1-ex1)*10.)/10.), (ceil((x2+ex2)*10.)/10.));
      }
    }
  }
  return true;
};


bool computeMeanAndVariance(FitMap_t& fitMap, const std::string& col, const bool& useEtaCM, const BinMapMap_t& binMap)
{
  const FitMap_t tmp = fitMap;
  // Fill the new graphs
  for (const auto& var : tmp) {
    for (const auto& lbl : var.second) {
      // Loop over all bins
      double sumW = 0.0 , sumY2 = 0.0 , sumY = 0.0;
      for (const auto& fit : lbl.second) {
        const double ptMin  = std::get<1>(fit.first).first;
        const double ptMax  = std::get<1>(fit.first).second;
        const double etaMin = std::get<0>(fit.first).first;
        const double etaMax = std::get<0>(fit.first).second;
        const double etaHW  = (etaMax - etaMin)/2.0;
        if (etaHW > 2.4) continue;
        const double y  = fit.second->GetParameter(0);
        const double ey = fit.second->GetParError(0);
        sumW  += (1.0 / ey);
        sumY  += (y / ey);
        sumY2 += ((y*y) / ey);
      }
      // Compute Weighted Mean
      const double wMean = ( sumY / sumW );
      // Compute Weighted Variance
      const double wVar = std::sqrt( ( sumY2 - (sumW*wMean*wMean) ) / ( sumW - 1.0 ) );
      // Add the results to the FitMap
      const double etaMin = binMap.at(col).at(useEtaCM ? "EtaCM" : "Eta")[0];
      const Par_t eta = std::make_pair(etaMin, etaMin+0.06);
      const Par_t pt  = std::make_pair(binMap.at(col).at("Pt")[0], binMap.at(col).at("Pt")[1]);
      const Bin_t bin = std::make_tuple(eta, pt);
      fitMap.at(var.first).at(lbl.first)[bin] = new TF1(*fitMap.at(var.first).at(lbl.first).begin()->second);
      fitMap.at(var.first).at(lbl.first).at(bin)->SetParameter(0, wMean);
      fitMap.at(var.first).at(lbl.first).at(bin)->SetParError(0, wVar);
    }
  }
  return true;
};


void formatFit(FitMap_t& fitMap)
{
  // Set the format of all fits
  for (auto& var : fitMap) {
    for (auto& lbl : var.second) {
      for (auto& fit : lbl.second) {
        if (fit.second) {
          // General
          fit.second->SetLineColor(kBlack);
          fit.second->SetLineWidth(3);
        }
      }
    }
  }
};


void getFitText(const TF1& fit, std::vector<std::string>& text)
{
  // Clean the text
  text.clear();
  // Set the Chi2 / NDF text
  text.push_back(Form("#chi^{2} / ndf     %.3f / %d", fit.GetChisquare(), fit.GetNDF()));
  // Set the parameters text
  for (int i = 0; i < fit.GetNpar(); i++) {
    if (fit.GetParError(i)>0.0) {
      text.push_back(Form("%s     %.3f #pm %.3f", fit.GetParName(i), fit.GetParameter(i), fit.GetParError(i)));
    }
  }
};


void clearFit(FitMap_t& fitMap)
{
  // Delete pointers to TGraph
  for (auto& var : fitMap) {
    for (auto& lbl : var.second) {
      for (auto& fit : lbl.second) {
        if (fit.second) delete fit.second; 
      }
    }
  }
  // Clear the graph map
  fitMap.clear();
};


bool fillExGraph(const FitMap_t& fitMap, GraphMap_t& exGraphMap)
{
  // Clean the input graph Map
  clearGraph(exGraphMap);
  // Fill the new graphs
  for (const auto& var : fitMap) {
    for (const auto& lbl : var.second) {
      int i = 0;
      const uint n = lbl.second.size();
      Double_t x[n] , y[n] , ex[n] , ey[n];
      double minPt , maxPt;
      for (const auto& fit : lbl.second) {
        const double ptMin  = std::get<1>(fit.first).first;  minPt = ptMin;
        const double ptMax  = std::get<1>(fit.first).second; maxPt = ptMax;
        const double etaMin = std::get<0>(fit.first).first;
        const double etaMax = std::get<0>(fit.first).second;
        const double etaMid = (etaMax + etaMin)/2.0;
        const double etaHW  = (etaMax - etaMin)/2.0;
        if (etaHW > 2.4) continue;
        x[i]  = etaMid;
        ex[i] = etaHW;
        y[i]  = fit.second->GetParameter(0);
        ey[i] = fit.second->GetParError(0);
        i++;
      }
      Par_t eta = std::make_pair(0., 0.);
      Par_t pt  = std::make_pair(minPt, maxPt);
      Bin_t bin = std::make_tuple(eta, pt);
      exGraphMap[var.first][lbl.first][bin] = new TGraphErrors(n, x, y, ex, ey);
      if (exGraphMap.at(var.first).at(lbl.first).at(bin)==NULL) {
        std::cout << "[ERROR] Element in exGraphMap is NULL!" << std::endl; return false;
      }
      std::string name = Form("exGraph_ETA_%s_%s", var.first.c_str(), lbl.first.c_str());
      exGraphMap.at(var.first).at(lbl.first).at(bin)->SetName(name.c_str());
    }
  }
  return true;
};


void formatExGraph(GraphMap_t& exGraphMap, const std::string& col, const bool& useEtaCM, const BinMapMap_t& binMap)
{
  // Set the format of all extrapolated graphs
  formatGraph(exGraphMap);
  for (auto& var : exGraphMap) {
    for (auto& lbl : var.second) {
      for (auto& graph : lbl.second) {
        if (graph.second) {
          // General
          // X-axis
          const double etaMin = binMap.at(col).at(useEtaCM ? "EtaCM" : "Eta")[0];
          const double etaMax = binMap.at(col).at(useEtaCM ? "EtaCM" : "Eta")[1];
          graph.second->GetXaxis()->SetLimits(etaMin, etaMax);
          if (useEtaCM) { graph.second->GetXaxis()->SetTitle("Leading Muon #eta_{CM}");  }
          else          { graph.second->GetXaxis()->SetTitle("Leading Muon #eta_{LAB}"); }
          // Y-axis
        }
      }
    }
  }
};


void drawExGraph(const std::string& outDir, const GraphMap_t& exGraphMap, const bool& useEtaCM)
{
  // Make output directories
  makeDir(outDir+"Plot/Final/png/");
  makeDir(outDir+"Plot/Final/pdf/");
  makeDir(outDir+"Plot/Final/root/");
  Double_t xl1=.52, yl1=0.74, xl2=xl1+.2, yl2=yl1+.150;
  // Set the format of all graphs
  for (auto& var : exGraphMap) {
    for (auto& lbl : var.second) {
      for (auto& graph : lbl.second) {
        if (graph.second) {
          // Create Canvas
          TCanvas* c   = new TCanvas("c", "c", 1000, 1000); c->cd();
          // Create Legend
          TLegend *leg = new TLegend(xl1,yl1,xl2,yl2);
          formatLegendEntry(*leg->AddEntry(graph.second, "Data", "p"));
          // Create the Text Info
          TLatex *tex = new TLatex(); tex->SetNDC(); tex->SetTextSize(0.025); float dy = 0;
          std::vector< std::string > cutSelection;
          getBinText(graph.first, cutSelection, ((lbl.first.find("Pl")!=std::string::npos)*2-1), useEtaCM);
          // Draw graph
          updateYAxisRange(graph.second, NULL);
          graph.second->Draw("ap");
          TGraphErrors* tmp = new TGraphErrors(*graph.second); tmp->Set(1);
          tmp->SetFillColor(kGreen+3); tmp->SetMarkerColor(kOrange);
          tmp->Draw("samep2");          
          // Draw the text
          for (const auto& s: cutSelection) { tex->DrawLatex(0.20, 0.85-dy, s.c_str()); dy+=0.04; }
          // Draw the legend
          leg->Draw("same");
          c->Modified(); c->Update();
          // set the CMS style
          int option = 111;
          if (lbl.first.find("pPb")!=std::string::npos) option = 109;
          if (lbl.first.find("Pbp")!=std::string::npos) option = 110;
          CMS_lumi(c, option, 33, "");
          c->Modified(); c->Update();
          // Save Canvas
          c->SaveAs(( outDir + "Plot/Final/png/" + graph.second->GetName() + ".png" ).c_str());
          c->SaveAs(( outDir + "Plot/Final/pdf/" + graph.second->GetName() + ".pdf" ).c_str());
          c->SaveAs(( outDir + "Plot/Final/root/" + graph.second->GetName() + ".root" ).c_str());
          // Clean up memory
          c->Clear(); c->Close();
          delete c;
          delete leg;
          delete tex;
          delete tmp;
        }
      }
    }
  }
};


bool createInputFiles(const std::string& dirPath , const FitMap_t& iniFitMap, const bool& useEtaCM, const BinMapMap_t& binMap)
{
  // Define inputFile as a matrix where the row index are the bins and column index are the headers
  std::map< std::string , std::map< Bin_t , std::vector< std::string > > > inputFile;
  // Initialize the fit map
  FitMap2_t fitMap;
  for (const auto& var : iniFitMap) {
    for (const auto& lbl : var.second) {
      for (const auto& fit : lbl.second) {
        fitMap[lbl.first][fit.first][var.first] = fit.second;
      }
    }
  }
  // Initialize the new input files
  std::string name = "InitialParam_MET_QCDToMu";
  std::map< std::string , std::map< std::string , uint > > COLMAP;
  for (const auto& lbl : fitMap) {
    std::string chg = "";
    if (lbl.first.find("Pl_")!=std::string::npos) { chg = "Pl"; }
    if (lbl.first.find("Mi_")!=std::string::npos) { chg = "Mi"; }
    std::string col = "";
    if (lbl.first.find("_pPb")!=std::string::npos) { col = "pPb"; }
    if (lbl.first.find("_Pbp")!=std::string::npos) { col = "Pbp"; }
    if (lbl.first.find("_PA") !=std::string::npos) { col = "PA";  }
    std::string inputFileName = name + "_" + col + ".csv";
    // Create the header
    std::vector< std::string > row;
    if (inputFile.count(inputFileName)==0) {
      std::string header = "";
      header = "Muon_Eta"; if (useEtaCM) { header = "Muon_EtaCM"; };  COLMAP[inputFileName][header] = row.size(); row.push_back(header);
      header = Form("Model_QCDToMu_%s", col.c_str()); COLMAP[inputFileName][header] = row.size(); row.push_back(header);
      if (chg!="") {
        for (const auto& var : lbl.second.begin()->second) {
          header = Form("%s_QCDToMuPl_%s", var.first.c_str(), col.c_str()); COLMAP[inputFileName][header] = row.size(); row.push_back(header);
          header = Form("%s_QCDToMuMi_%s", var.first.c_str(), col.c_str()); COLMAP[inputFileName][header] = row.size(); row.push_back(header);
        }
      }
      else {
        for (const auto& var : lbl.second.begin()->second) {
          header = Form("%s_QCDToMu_%s", var.first.c_str(), col.c_str()); COLMAP[inputFileName][header] = row.size(); row.push_back(header);
        }
      }
      Bin_t bin = std::make_tuple(std::make_pair(-99., -99.), std::make_pair(-99., -99.)); // header dummy bin
      inputFile[inputFileName][bin] = row;
    }
    std::map< std::string , uint > colIdx = ( (COLMAP.count(inputFileName)>0) ? COLMAP.at(inputFileName) : std::map< std::string , uint >() );
    // Fill the new input files
    for (const auto& bin : lbl.second) {
      // Check that pT range is consistent
      const double ptMin  = std::get<1>(bin.first).first;
      const double ptMax  = std::get<1>(bin.first).second;
      // Set the Eta Label
      const double etaMin = std::get<0>(bin.first).first;
      const double etaMax = std::get<0>(bin.first).second;
      if (std::abs(etaMax-etaMin) < 0.1) continue;
      // Initialize the values of the new bin row
      if (inputFile.at(inputFileName).count(bin.first)==0) {
        for (const auto& ele : row) { inputFile.at(inputFileName)[bin.first].push_back(""); }
      }
      std::string header = "";
      header = "Muon_Eta"; if (useEtaCM) { header = "Muon_EtaCM"; }
      inputFile.at(inputFileName).at(bin.first).at(colIdx[header]) = Form("%s%.2f_%s%.2f", (etaMin>0 ? "+" : "-"), std::abs(etaMin), (etaMax>0 ? "+" : "-"), std::abs(etaMax));
      // Set the Model Label
      header = Form("Model_QCDToMu_%s", col.c_str()); inputFile.at(inputFileName).at(bin.first).at(colIdx[header]) = MODELNAME_;
      const double etaMeanMin = binMap.at(col).at(useEtaCM ? "EtaCM" : "Eta")[0];
      const double etaIncMin  = binMap.at(col).at(useEtaCM ? "EtaCM_Inc" : "Eta_Inc")[0];
      const double etaIncMax  = binMap.at(col).at(useEtaCM ? "EtaCM_Inc" : "Eta_Inc")[1];
      // Set each variable Label
      if ((useBin=="Inclusive") && lbl.second.count(std::make_tuple(std::make_pair(etaIncMin,etaIncMax), std::make_pair(ptMin, ptMax)))>0) {
        for (const auto& var : lbl.second.at(std::make_tuple(std::make_pair(etaIncMin,etaIncMax), std::make_pair(ptMin, ptMax)))) {
          double val  = var.second->GetParameter(0);
          double err = var.second->GetParError(0);
          header = Form("%s_QCDToMu%s_%s", var.first.c_str(), chg.c_str(), col.c_str());  inputFile.at(inputFileName).at(bin.first).at(colIdx[header]) = Form("[%.4f;%.4f]", val, err);
        }
      }
      else if ((useBin=="Mean") && lbl.second.count(std::make_tuple(std::make_pair(etaMeanMin,etaMeanMin+0.06), std::make_pair(ptMin, ptMax)))>0) {
        for (const auto& var : lbl.second.at(std::make_tuple(std::make_pair(etaMin,etaMin+0.06), std::make_pair(ptMin, ptMax)))) {
          double val  = var.second->GetParameter(0);
          double err = var.second->GetParError(0);
          header = Form("%s_QCDToMu%s_%s", var.first.c_str(), chg.c_str(), col.c_str());  inputFile.at(inputFileName).at(bin.first).at(colIdx[header]) = Form("[%.4f;%.4f]", val, err);
        }
      }
      else {
        for (const auto& var : bin.second) {
          double val  = var.second->GetParameter(0);
          double err = var.second->GetParError(0);
          header = Form("%s_QCDToMu%s_%s", var.first.c_str(), chg.c_str(), col.c_str());  inputFile.at(inputFileName).at(bin.first).at(colIdx[header]) = Form("[%.4f;%.4f]", val, err);
        }
      }      
    }
  }
  // Print the new input files
  makeDir(dirPath);
  for (const auto& name : inputFile) {
    ofstream fout( dirPath + name.first );
    // Print all except inclusive bins
    for (const auto& row : name.second) {
      if (std::abs(std::get<0>(row.first).first-std::get<0>(row.first).second)<2.4 || std::abs(std::get<0>(row.first).second)==99.) {
        std::string rowLine = "";
        for (const auto& col : row.second) {
          rowLine += col + ",";
        }
        fout << rowLine << endl;
      }
    }
    // Print inclusive bins
    for (const auto& row : name.second) {
      if (std::abs(std::get<0>(row.first).first-std::get<0>(row.first).second)>2.4 && std::abs(std::get<0>(row.first).second)!=99.) {
        std::string rowLine = "";
        for (const auto& col : row.second) {
          rowLine += col + ",";
        }
        fout << rowLine << endl;
      }
    }
  }
  return true;
};

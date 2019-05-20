#if !defined(__CINT__) || defined(__MAKECINT__)
// Auxiliary Headers
#include "../../Utilities/HiEvtTree.h"
#include "../../Utilities/HiMuonTree.h"
#include "../../Utilities/EVENTUTILS.h"
// ROOT headers
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TFrame.h"
#include "TAxis.h"
#include "TPaletteAxis.h"
// c++ headers
#include <iostream>
#include <map>
#include <vector>
#include <string>
// CMS headers
#include "../../Utilities/CMS/tdrstyle.C"
#include "../../Utilities/CMS/CMS_lumi.C"

#endif


// ------------------ TYPE -------------------------------
using TH2DMap_t    =  std::map< std::string , std::map< std::string , std::map< std::string , TH2D > > >;
using BinMap_t     =  std::map< std::string , std::vector< double > >;
using BinMapMap_t  =  std::map< std::string , BinMap_t >;


// ------------------ FUNCTION -------------------------------
void     initHist2D ( TH2DMap_t& h2D    , const BinMapMap_t& binMap );
void     drawHist2D ( TH2DMap_t& h2DMap , const std::string& outDir );

// ------------------ GLOBAL ---------------------------------
//
// Collision System
const std::vector< std::string > COLL_ = { "pPb" , "Pbp" , "PA" };
//
// Muon Charge
const std::vector< std::string > CHG_  = { "Plus" , "Minus" };
//
// Efficiency Categories
const std::vector< std::string > VAR_  = {"EtaCM_xPb", "EtaCM_xp"};
//
// Input Files for analysis
const std::string path_MC = "/data_CMS/cms/stahl";
const std::map< std::string , std::vector< std::pair< std::string , double > > > inputFileMap_ = {
  {"MC_WToMuNu_Plus_pPb"      , { { Form("%s/%s", path_MC.c_str(), "POWHEG/HiEWQForest_Embedded_Official_POWHEG_CT14_EPPS16_WToMuNu_Plus_pPb_8160GeV_20171003.root")  , POWHEG::XSec.at("WToMuNu_Plus").at("pPb")  } } },
  {"MC_WToMuNu_Minus_pPb"     , { { Form("%s/%s", path_MC.c_str(), "POWHEG/HiEWQForest_Embedded_Official_POWHEG_CT14_EPPS16_WToMuNu_Minus_pPb_8160GeV_20171003.root") , POWHEG::XSec.at("WToMuNu_Minus").at("pPb") } } },
  {"MC_WToMuNu_Plus_Pbp"      , { { Form("%s/%s", path_MC.c_str(), "POWHEG/HiEWQForest_Embedded_Official_POWHEG_CT14_EPPS16_WToMuNu_Plus_Pbp_8160GeV_20171003.root")  , POWHEG::XSec.at("WToMuNu_Plus").at("Pbp")  } } },
  {"MC_WToMuNu_Minus_Pbp"     , { { Form("%s/%s", path_MC.c_str(), "POWHEG/HiEWQForest_Embedded_Official_POWHEG_CT14_EPPS16_WToMuNu_Minus_Pbp_8160GeV_20171003.root") , POWHEG::XSec.at("WToMuNu_Minus").at("Pbp") } } }
};
std::map< std::string , std::vector< std::string > > sampleType_;


void makeXPbvsEtaPlot(const std::string workDirName = "NominalCM", const bool remakeHist2D = false)
{
  //
  // Initialize the Kinematic Bin info
  BinMapMap_t  BIN;
  if (workDirName.find("NominalCM")!=std::string::npos) {
    BinMap_t  TMP_pPb = {
      { "EtaCM" , { -3.0 , -2.86 , -2.60 , -2.40 , -2.20 , -1.93 , -1.80 , -1.60 , -1.40 , -1.20 , -1.00 , -0.80 , -0.60 , -0.40 , -0.20 , 0.00 , 0.20 , 0.40 , 0.60 , 0.80 , 1.00 , 1.20 , 1.40 , 1.60 , 1.80, 1.93 , 2.1 } }
    };
    BinMap_t  TMP_Pbp = {
      { "EtaCM" , { -2.1 , -1.93 , -1.80 , -1.60 , -1.40 , -1.20 , -1.00 , -0.80 , -0.60 , -0.40 , -0.20 , 0.00 , 0.20 , 0.40 , 0.60 , 0.80 , 1.00 , 1.20 , 1.40 , 1.60 , 1.80 , 1.93 , 2.20 , 2.40 , 2.60 , 2.86 , 3.0 } }
    };
    const double xMin = 0.0001 , xMax = 100.0;
    const uint xBins = 200;
    const double xFactor = std::pow( (xMax/xMin) , (1.0/xBins) );
    for (uint i=0; i<xBins; i++) {
      TMP_pPb["xPb"].push_back(xMin*std::pow(xFactor, i));
      TMP_pPb["xp" ].push_back(xMin*std::pow(xFactor, i));
      TMP_Pbp["xPb"].push_back(xMin*std::pow(xFactor, i));
      TMP_Pbp["xp" ].push_back(xMin*std::pow(xFactor, i));
    }
    for (const auto& v : TMP_pPb) { BIN["PA"][v.first] = v.second; BIN["pPb"][v.first] = v.second; }
    for (const auto& v : TMP_Pbp) { BIN["Pbp"][v.first] = v.second; }
  }
  //
  // Change the working directory
  const std::string CWD = getcwd(NULL, 0);
  const std::string mainDir = Form("%s/Output/", CWD.c_str());
  gSystem->mkdir(mainDir.c_str(), kTRUE);
  gSystem->ChangeDirectory(mainDir.c_str());
  //
  // ------------------------------------------------------------------------------------------------------------------------
  //
  // Declare the histograms
  TH2DMap_t h2D;
  //
  // Initialize the 2D histograms
  initHist2D(h2D , BIN);
  //
  // ------------------------------------------------------------------------------------------------------------------------
  //
  // Extract the histograms
  bool histFound = true;
  if (!remakeHist2D) {
    const std::string fileName = "h2DInfo.root";
    TFile file((mainDir + fileName).c_str(), "READ"); file.cd();
    if (file.IsOpen()==true && file.IsZombie()==false) {
      for (const auto& col : COLL_) {
        for (const auto& chg : CHG_) {
          for (const auto& var : VAR_) {
            const std::string& histName = h2D.at(col).at(chg).at(var).GetName();
            if (file.Get(histName.c_str())) { h2D.at(col).at(chg).at(var) = *(TH2D*)file.Get(histName.c_str()); }
            else { std::cout << "[ERROR] Histogram " << histName << " was not found in " << fileName << std::endl; histFound = false; break; }
          }
        }
      }
    }
    else { std::cout << "[ERROR] File " << fileName << " was not found!" << std::endl; histFound = false; }
    file.Close();
  }
  //
  // ------------------------------------------------------------------------------------------------------------------------
  //
  // Fill the histograms
  if (remakeHist2D || !histFound) {
    //
    // ------------------------------------------------------------------------------------------------------------------------
    //
    // Extract all the samples
    std::map< std::string , Long64_t > nentries;
    std::map< std::string , std::unique_ptr< HiMuonTree > > muonTree;
    std::map< std::string , std::unique_ptr< HiMETTree  > > metTree;
    std::map< std::string , std::unique_ptr< HiEvtTree  > > evtTree;
    for (const auto & inputFile : inputFileMap_) {
      const std::string sample = inputFile.first;
      std::vector< std::string >  fileInfo;
      for (const auto& f : inputFile.second) { fileInfo.push_back(f.first); }
      //
      muonTree[sample] = std::unique_ptr<HiMuonTree>(new HiMuonTree());
      if (!muonTree.at(sample)->GetTree(fileInfo)) return;
      nentries[sample] = muonTree.at(sample)->GetEntries();
      //
      evtTree[sample] = std::unique_ptr<HiEvtTree>(new HiEvtTree());
      if (!evtTree.at(sample)->GetTree(fileInfo)) return;
      if (evtTree.at(sample)->GetEntries() != nentries.at(sample)) { std::cout << "[ERROR] Inconsistent number of entries between Event (" << 
          evtTree.at(sample)->GetEntries() << ") and Muon Tree (" << muonTree.at(sample)->GetEntries() << ") !" << std::endl; return; }
    }
    //
    // ------------------------------------------------------------------------------------------------------------------------
    //
    // Loop over the samples
    //
    for (auto & inputFile : inputFileMap_) {
      const std::string sample = inputFile.first;
      // Loop over the events
      int treeIdx = -1;
      double crossSection = -1.0;
      std::cout << "[INFO] Starting to process " << nentries.at(sample) << " nentries" << std::endl;
      //
      for (Long64_t jentry = 0; jentry < nentries.at(sample); jentry++) { //
        //
        // Get the entry in the trees
        if (muonTree.at(sample)->GetEntry(jentry)<0    ) { std::cout << "[ERROR] Muon Tree invalid entry!"     << std::endl; return; }
        if ( evtTree.at(sample)->GetEntry(jentry)<0    ) { std::cout << "[ERROR] Event Tree invalid entry!"    << std::endl; return; }
        //
        if (muonTree.at(sample)->Chain()->GetTreeNumber()!=treeIdx) {
          treeIdx = muonTree.at(sample)->Chain()->GetTreeNumber();
          std::cout << "[INFO] Processing Root File: " << inputFile.second[treeIdx].first << std::endl;
          // Get the MC Cross-Section
          crossSection = inputFile.second[treeIdx].second;
        }
        //
        loadBar(jentry, nentries.at(sample));
        // 
        // Check that the different tree agrees well
        if (muonTree.at(sample)->Event_Run()    != evtTree.at(sample)->run()) { std::cout << "[ERROR] Event Run does not agree!"     << std::endl; return; }
        if (muonTree.at(sample)->Event_Number() != evtTree.at(sample)->evt()) { std::cout << "[ERROR] Event Event does not agree!"   << std::endl; return; }
        //
        // Determine the collision system of the sample
        std::string col = "";
        if (sample.find("Pbp")!=std::string::npos) col = "Pbp"; // for Pbp
        if (sample.find("pPb")!=std::string::npos) col = "pPb"; // for pPb
        if (col=="") { std::cout << "[ERROR] Could not determine the collision system in the sample" << std::endl; return; }
        //
        // Determine the type of sample : i.e. MC_DYToMuMu
        std::string sampleType = sample;
        sampleType = sampleType.substr(0, (sampleType.find(col)-1));
        //
        // Get the Lumi re-weight for MC (global weight)
        const double lumi = ( (col=="pPb") ? PA::LUMI::Data_pPb : PA::LUMI::Data_Pbp );
        const double mcWeight = ( ( crossSection * lumi ) / muonTree.at(sample)->GetTreeEntries() );
        //
        // Muon Based Information
        //
        // Find the Muon associated to the W boson
        int muGenIdx = -1;
        for (uint iGenMu = 0; iGenMu < muonTree.at(sample)->Gen_Muon_Mom().size(); iGenMu++) {
          const auto& mom = muonTree.at(sample)->MuonMother(iGenMu);
          if (mom.pdg==24) { muGenIdx = iGenMu; break; }
        }
        if (muGenIdx<0) { std::cout << "[ERROR] Muon mother was not found!" << std::endl; return; }
        //
        // Apply PT and Eta Cut
        TLorentzVector muGenP4 = muonTree.at(sample)->Gen_Muon_Mom()[muGenIdx];
        if (muGenP4.Pt()<25. || std::abs(muGenP4.Eta())>2.4) continue;
        const int muGenChg = muonTree.at(sample)->Gen_Muon_Charge()[muGenIdx];
        //
        // PDF Based Information
        //
        const double xp  = (col=="pPb" ? evtTree.at(sample)->pdfX().first  : evtTree.at(sample)->pdfX().second);
        const double xPb = (col=="pPb" ? evtTree.at(sample)->pdfX().second : evtTree.at(sample)->pdfX().first );
        //
        // Fill the histograms
        //
        for (const auto& c : COLL_) {
          if (c!="PA" && c!=col) continue;
          const double etaCM = PA::EtaLABtoCM( ( (c=="PA" && col=="Pbp") ? -1.0 : 1.0 )*muGenP4.Eta()  , (c!="Pbp") );
          if ( ( (c=="PA" || c=="pPb") && (etaCM>-2.86 && etaCM<1.93) ) || ( (c=="Pbp") && (etaCM>-1.93 && etaCM<2.86) ) ) {
            if (muGenChg>0) {
              h2D.at(c).at("Plus").at("EtaCM_xPb").Fill( etaCM , xPb , mcWeight );
              h2D.at(c).at("Plus").at("EtaCM_xp" ).Fill( etaCM , xp  , mcWeight );
            }
            else if (muGenChg<0) {
              h2D.at(c).at("Minus").at("EtaCM_xPb").Fill( etaCM , xPb , mcWeight );
                h2D.at(c).at("Minus").at("EtaCM_xp" ).Fill( etaCM , xp  , mcWeight );
            }
          }
        }
        //
      }
    }
    //
    // Save the 2D histograms
    const std::string fileName = "h2DInfo.root";
    TFile file((mainDir + fileName).c_str(), "RECREATE"); file.cd();
    if (file.IsOpen()==true && file.IsZombie()==false) {
      for (const auto& col : COLL_) {
        for (const auto& chg : CHG_) {
          for (const auto& var : VAR_) {
            h2D.at(col).at(chg).at(var).Write();
          }
        }
      }
      file.Write();
    }
    else { std::cout << "[ERROR] File " << fileName << " failed to open!" << std::endl; }
    file.Close();
  }
  //
  // ------------------------------------------------------------------------------------------------------------------------
  //
  // Draw the 2D histograms
  drawHist2D(h2D, mainDir);
  //  
};


void initHist2D(TH2DMap_t& h2D, const BinMapMap_t& binMap)
{
  for (const auto& col : COLL_) {
    for (const auto& chg : CHG_) {
      for (const auto& var : VAR_) {
        const std::string name = "h2D_" + col +"_"+ chg +"_"+ var;
        if (var=="EtaCM_xPb") { h2D[col][chg][var] = TH2D(name.c_str(), name.c_str(),
                                                          (binMap.at(col).at("EtaCM").size()-1), binMap.at(col).at("EtaCM").data(),
                                                          (binMap.at(col).at("xPb").size()-1), binMap.at(col).at("xPb").data()); }
        if (var=="EtaCM_xp" ) { h2D[col][chg][var] = TH2D(name.c_str(), name.c_str(),
                                                          (binMap.at(col).at("EtaCM").size()-1), binMap.at(col).at("EtaCM").data(),
                                                          (binMap.at(col).at("xp").size()-1), binMap.at(col).at("xp").data()); }
        if (var=="EtaCM_xPb") { h2D.at(col).at(chg).at(var).SetTitle(";#eta^{#mu}_{CM};x_{Pb}"); }
        if (var=="EtaCM_xp" ) { h2D.at(col).at(chg).at(var).SetTitle(";#eta^{#mu}_{CM};x_{p}");  }

        std::string xLabel = "#eta^{#mu}_{CM}";
        std::string yLabel = "x_{Pb}";
        h2D.at(col).at(chg).at(var).Sumw2();
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
  TGaxis::SetMaxDigits(3); // to display powers of 10
  //
  // Set Palette
  gStyle->SetPalette(55);
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
  //
};


void drawHist2D(TH2DMap_t& h2DMap, const std::string& outDir)
{
  //
  // Set Style
  setStyle();
  TGaxis::SetMaxDigits(2); // to display powers of 10
  //
  for (const auto& c : h2DMap) {
    for (const auto& ch : c.second) {
      for (const auto& v : ch.second) {
        //
        const auto& col = c.first;
        const auto& chg = ch.first;
        const auto& var = v.first;
        auto& h2D = h2DMap.at(col).at(chg).at(var);
        //
        for (int i=0; i<=h2D.GetNbinsX(); i++) {
          for (int j=0; j<=h2D.GetNbinsY(); j++) {
            const auto& binWidthX = h2D.GetXaxis()->GetBinWidth(i);
            const auto& binValue  = h2D.GetBinContent(i, j);
            if (binValue>0.0 && binWidthX>0.0) { h2D.SetBinContent(i, j, (binValue/binWidthX)); }
          }
        }
	//
	// Create the Text Info
	TLatex tex; tex.SetNDC(); tex.SetTextSize(0.045); float dy = 0;
	std::vector< std::string > textToPrint;
	std::string sampleLabel = "W #rightarrow #mu + #nu_{#mu}";
	if (chg == "Plus") { sampleLabel = "W^{+} #rightarrow #mu^{+} + #nu_{#mu}"; }
	if (chg == "Minus") { sampleLabel = "W^{#font[122]{\55}} #rightarrow #mu^{#font[122]{\55}} + #bar{#nu}_{#mu}"; }
	textToPrint.push_back(sampleLabel);
	textToPrint.push_back("p^{#mu}_{T} > 25 GeV/c");
        //
        // Create Canvas
        TCanvas c("c", "c", 1000, 1000); c.cd();
        c.SetRightMargin(2.8);
        //
        // Format graph
        //
        // X-axis
        h2D.GetXaxis()->CenterTitle(kTRUE);
        h2D.GetXaxis()->SetTitleOffset(0.70);
        h2D.GetXaxis()->SetTitleSize(0.065);
        h2D.GetXaxis()->SetLabelSize(0.035);
        // Y-axis
        h2D.GetYaxis()->CenterTitle(kTRUE);
        h2D.GetYaxis()->SetTitleOffset(1.05);
        h2D.GetYaxis()->SetTitleSize(0.065);
        h2D.GetYaxis()->SetLabelSize(0.035);
        // Z axis
        h2D.GetZaxis()->SetLabelSize(0.030);
        // Set axis labels
        if (var=="EtaCM_xPb") { h2D.SetTitle(";#eta^{#mu}_{CM};x_{Pb}"); }
        if (var=="EtaCM_xp" ) { h2D.SetTitle(";#eta^{#mu}_{CM};x_{p}");  }
        //
        // Draw graph
        h2D.Draw("COLZ");
        //
	// Draw the text
	tex.SetTextSize(0.050); tex.DrawLatex(0.20, 0.77, textToPrint[0].c_str());
        tex.SetTextSize(0.050); tex.DrawLatex(0.60, 0.77, "POWHEG v2");
        tex.SetTextSize(0.058); tex.SetTextFont(61); tex.DrawLatex(0.20, 0.85, "CMS"); tex.SetTextFont(62);
        //tex.SetTextSize(0.044); tex.SetTextFont(52); tex.DrawLatex(0.33, 0.85, "Simulation Preliminary"); tex.SetTextFont(62);
        tex.SetTextSize(0.044); tex.SetTextFont(52); tex.DrawLatex(0.33, 0.85, "Simulation"); tex.SetTextFont(62);
	if (textToPrint.size()>1) { tex.SetTextSize(0.035); tex.DrawLatex(0.20, 0.70, textToPrint[1].c_str()); }
        c.Modified(); c.Update();
        //
        // set the CMS style
        int option = 1141;
        CMS_lumi(&c, option, 33, "", false, 0.6, false);
        // Set y axis log
        c.SetLogy(); c.Modified(); c.Update();
        // Set z axis log
        //c.SetLogz(); c.Modified(); c.Update();
        // Update
        c.Modified(); c.Update(); // Pure paranoia
        TPaletteAxis *palette = (TPaletteAxis*)h2D.GetListOfFunctions()->FindObject("palette");
        palette->SetX2NDC(0.93);
        c.Modified();
        //
        // Save canvas
        //
        // Create Output Directory
        const std::string plotDir = outDir + "Plot/" + col+"/" + var;
        makeDir(plotDir + "/png/");
        makeDir(plotDir + "/pdf/");
        makeDir(plotDir + "/root/");
        //
        // Save Canvas
        const std::string name = Form("h2D_WToMu%s_%s_%s", chg.c_str(), col.c_str(), var.c_str());
        c.SaveAs(( plotDir + "/png/"  + name + ".png"  ).c_str());
        c.SaveAs(( plotDir + "/pdf/"  + name + ".pdf"  ).c_str());
        c.SaveAs(( plotDir + "/root/" + name + ".root" ).c_str());
        //
        // Clean up memory
        c.Clear(); c.Close();
      }
    }
  }
};
        

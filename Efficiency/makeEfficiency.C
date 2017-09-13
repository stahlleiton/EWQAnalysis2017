#if !defined(__CINT__) || defined(__MAKECINT__)
// Auxiliary Headers
#include "../Utilities/HiMETTree.h"
#include "../Utilities/HiMuonTree.h"
#include "../Utilities/EVENTUTILS.h"
// ROOT headers
#include "TROOT.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TFrame.h"
#include "TLine.h"
#include "TBox.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TPaletteAxis.h"
#include "TMath.h"
// c++ headers
#include <dirent.h>
#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <chrono>
// CMS headers
#include "../Utilities/CMS/tdrstyle.C"
#include "../Utilities/CMS/CMS_lumi.C"

#endif


// ------------------ TYPE -------------------------------
using TH1DMap_t    =  std::map< std::string , std::map< std::string , std::map< std::string , std::map< std::string , std::map< std::string , std::tuple< TH1D , TH1D , double > > > > > >;
using TH2DMap_t    =  std::map< std::string , std::map< std::string , std::map< std::string , std::map< std::string , std::map< std::string , std::tuple< TH2D , TH2D , double > > > > > >;
using EffMap_t     =  std::map< std::string , std::map< std::string , std::map< std::string , std::map< std::string , std::map< std::string , TEfficiency > > > > >;
using VarMap_t     =  std::map< std::string , std::map< std::string , std::map< std::string , std::map< std::string , std::vector< double > > > > >;
using BinMap_t     =  std::map< std::string , std::vector< double > >;
using Bin2DMap_t   =  std::map< std::string , std::pair< std::vector< double > , std::vector< double > > >;
using FileInfo_t   =  std::vector< std::pair< std::string , double > >;


// ------------------ FUNCTION -------------------------------
void     initEff1D          ( TH1DMap_t& h , const BinMap_t& binMap );
bool     fillEff1D          ( TH1DMap_t& h , const bool& pass , const std::string& type , const VarMap_t& var , const double sfTnP=1.0 , const double evtWeight=1.0 );
bool     loadEff1D          ( EffMap_t& eff, const TH1DMap_t& h );
void     formatEff1D        ( TEfficiency& eff , const std::string& var , const std::string& charge , const std::string& type );
void     drawEff1D          ( const std::string& outDir, EffMap_t& effMap );
void     initEff2D          ( TH2DMap_t& h , const Bin2DMap_t& binMap );
bool     fillEff2D          ( TH2DMap_t& h , const bool& pass , const std::string& type , const VarMap_t& var , const double sfTnP=1.0 , const double evtWeight=1.0 );
bool     loadEff2D          ( EffMap_t& eff, const TH2DMap_t& h );
void     formatEff2D        ( TEfficiency& eff , const std::string& var , const std::string& charge , const std::string& type);
void     drawEff2D          ( const std::string& outDir, EffMap_t& effMap );
void     formatComparison2D ( TH2& histo , const std::string& var , const std::string& charge , const std::string& type );
void     compareEff2D       ( const std::string& outDir , EffMap_t& effMap , const std::string comp );
void     mergeEff           ( EffMap_t& eff );
void     writeEff           ( TFile& file , const EffMap_t& eff , const std::string& mainDirName );
void     saveEff            ( const std::string& outDir , const EffMap_t& eff1D , const EffMap_t& eff2D );
void     setGlobalWeight    ( TH1DMap_t& h , const double& weight , const std::string& sample , const std::string& col );
void     setGlobalWeight    ( TH2DMap_t& h , const double& weight , const std::string& sample , const std::string& col );
void     setStyle           ( );
void     formatLegendEntry  ( TLegendEntry& e );
void     formatDecayLabel   ( std::string& label , const std::string& inLabel );
void     clearEff           ( EffMap_t& effMap );
void     getStats2D         ( double& mean , double& rms , double& min , double& max , const TH2& hist );
void     fillErrorEff2D     ( TH2D& hist , const TEfficiency& eff , const std::string& var);
bool     checkWeights       ( const TH1& pass , const TH1& total );
void     redrawBorder       ( );


// ------------------ GLOBAL -------------------------------
// Kinematic Info
const BinMap_t  MU_BIN_ = {
  { "Eta" , { -2.4 , -2.2 , -2.0 , -1.8 , -1.6 , -1.4 , -1.2 , -1.0 , -0.8 , -0.6 , -0.4 , -0.2 , 0.0 , 0.2 , 0.4 , 0.6 , 0.8 , 1.0 , 1.2 , 1.4 , 1.6 , 1.8 , 2.0 , 2.2 , 2.4 } },
  { "Pt"  , { 25., 30., 35., 40., 50, 70., 200. } }
};
const Bin2DMap_t MU_BIN_2D_ = {
  { "Pt_Eta" , { MU_BIN_.at("Pt") , MU_BIN_.at("Eta")} }
};
//
// Trigger Info
const ushort triggerIndex_ = PA::HLT_PAL3Mu12;
//
// Collision System
const std::vector< std::string > COLL_ = { "pPb" , "Pbp" , "PA" };
//
// Muon Charge
const std::vector< std::string > CHG_  = { "Plus" , "Minus" };
//
// Efficiency Categories
const std::vector< std::string > effType = {"Acceptance" , "Identification" , "Trigger" , "Isolation" , "Offline" , "Event" , "DrellYanVeto" , "Total" , "Acceptance_Efficiency"};
//
// Input Files for analysis
const std::string path_MC = "root://cms-xrd-global.cern.ch//store/group/phys_heavyions/anstahll/EWQAnalysis2017/pPb2016/8160GeV/MC/Embedded/Official";
const std::map< std::string , std::vector< std::pair< std::string , double > > > inputFileMap_ = {
  {"MC_QCDToMu_pPb"          , { { Form("%s/%s", path_MC.c_str(), "PYTHIA/HiEWQForest_Embedded_Official_PYTHIA8_QCD_PtHat20_MuPt15_pPb_8160GeV_20170813.root")         , PYTHIA::XSec.at("QCDToMu").at("pPb")           } } },
  {"MC_QCDToMu_Pbp"          , { { Form("%s/%s", path_MC.c_str(), "PYTHIA/HiEWQForest_Embedded_Official_PYTHIA8_QCD_PtHat20_MuPt15_Pbp_8160GeV_20170813.root")         , PYTHIA::XSec.at("QCDToMu").at("Pbp")           } } },
  {"MC_TTall_pPb"            , { { Form("%s/%s", path_MC.c_str(), "POWHEG/HiEWQForest_Embedded_Official_POWHEG_CT14_EPPS16_TTall_pPb_8160GeV_20170813.root")           , POWHEG::XSec.at("TTall").at("pPb")             } } },
  {"MC_TTall_Pbp"            , { { Form("%s/%s", path_MC.c_str(), "POWHEG/HiEWQForest_Embedded_Official_POWHEG_CT14_EPPS16_TTall_Pbp_8160GeV_20170813.root")           , POWHEG::XSec.at("TTall").at("Pbp")             } } },
  {"MC_ZToMuMu_M_30_Inf_pPb" , { { Form("%s/%s", path_MC.c_str(), "POWHEG/HiEWQForest_Embedded_Official_POWHEG_CT14_EPPS16_DYtoMuMu_M_30_pPb_8160GeV_20170813.root")   , POWHEG::XSec.at("DYToMuMu_M_30_Inf").at("pPb") } } },
  {"MC_ZToMuMu_M_30_Inf_Pbp" , { { Form("%s/%s", path_MC.c_str(), "POWHEG/HiEWQForest_Embedded_Official_POWHEG_CT14_EPPS16_DYtoMuMu_M_30_Pbp_8160GeV_20170813.root")   , POWHEG::XSec.at("DYToMuMu_M_30_Inf").at("Pbp") } } },
  {"MC_DYToMuMu_M_10_30_pPb" , { { Form("%s/%s", path_MC.c_str(), "POWHEG/HiEWQForest_Embedded_Official_POWHEG_CT14_EPPS16_DYtoMuMu_M_10_30_pPb_8160GeV_20170813.root"), POWHEG::XSec.at("DYToMuMu_M_10_30").at("pPb")  } } },
  {"MC_DYToMuMu_M_10_30_Pbp" , { { Form("%s/%s", path_MC.c_str(), "POWHEG/HiEWQForest_Embedded_Official_POWHEG_CT14_EPPS16_DYtoMuMu_M_10_30_Pbp_8160GeV_20170813.root"), POWHEG::XSec.at("DYToMuMu_M_10_30").at("Pbp")  } } },
  {"MC_WToTauNu_Plus_pPb"    , { { Form("%s/%s", path_MC.c_str(), "POWHEG/HiEWQForest_Embedded_Official_POWHEG_CT14_EPPS16_WToTauNu_Plus_pPb_8160GeV_20170813.root")   , POWHEG::XSec.at("WToTauNu_Plus").at("pPb")     } } },
  {"MC_WToTauNu_Minus_pPb"   , { { Form("%s/%s", path_MC.c_str(), "POWHEG/HiEWQForest_Embedded_Official_POWHEG_CT14_EPPS16_WToTauNu_Minus_pPb_8160GeV_20170813.root")  , POWHEG::XSec.at("WToTauNu_Minus").at("pPb")    } } },
  {"MC_WToTauNu_Plus_Pbp"    , { { Form("%s/%s", path_MC.c_str(), "POWHEG/HiEWQForest_Embedded_Official_POWHEG_CT14_EPPS16_WToTauNu_Plus_Pbp_8160GeV_20170813.root")   , POWHEG::XSec.at("WToTauNu_Plus").at("Pbp")     } } },
  {"MC_WToTauNu_Minus_Pbp"   , { { Form("%s/%s", path_MC.c_str(), "POWHEG/HiEWQForest_Embedded_Official_POWHEG_CT14_EPPS16_WToTauNu_Minus_Pbp_8160GeV_20170813.root")  , POWHEG::XSec.at("WToTauNu_Minus").at("Pbp")    } } },
  {"MC_WToMuNu_Plus_pPb"     , { { Form("%s/%s", path_MC.c_str(), "POWHEG/HiEWQForest_Embedded_Official_POWHEG_CT14_EPPS16_WToMuNu_Plus_pPb_8160GeV_20170813.root")    , POWHEG::XSec.at("WToMuNu_Plus").at("pPb")      } } },
  {"MC_WToMuNu_Minus_pPb"    , { { Form("%s/%s", path_MC.c_str(), "POWHEG/HiEWQForest_Embedded_Official_POWHEG_CT14_EPPS16_WToMuNu_Minus_pPb_8160GeV_20170813.root")   , POWHEG::XSec.at("WToMuNu_Minus").at("pPb")     } } },
  {"MC_WToMuNu_Plus_Pbp"     , { { Form("%s/%s", path_MC.c_str(), "POWHEG/HiEWQForest_Embedded_Official_POWHEG_CT14_EPPS16_WToMuNu_Plus_Pbp_8160GeV_20170813.root")    , POWHEG::XSec.at("WToMuNu_Plus").at("Pbp")      } } },
  {"MC_WToMuNu_Minus_Pbp"    , { { Form("%s/%s", path_MC.c_str(), "POWHEG/HiEWQForest_Embedded_Official_POWHEG_CT14_EPPS16_WToMuNu_Minus_Pbp_8160GeV_20170813.root")   , POWHEG::XSec.at("WToMuNu_Minus").at("Pbp")     } } }
};
std::map< std::string , std::vector< std::string > > sampleType_;


void makeEfficiency(void)
{
  // Change the working directory
  const std::string CWD = getcwd(NULL, 0);
  const std::string mainDir = Form("%s/Efficiency/", CWD.c_str());
  gSystem->mkdir(mainDir.c_str(), kTRUE);
  gSystem->ChangeDirectory(mainDir.c_str());
  //
  // Create list of samples
  sampleType_["sample"] = {};
  for (const auto& f : inputFileMap_) { std::string tmp = f.first; tmp = tmp.substr(0, tmp.find_last_of("_"));
    if (std::find(sampleType_.at("sample").begin(), sampleType_.at("sample").end(), tmp)==sampleType_.at("sample").end()) { sampleType_.at("sample").push_back(tmp); }
  }
  //
  // For computing time optimization
  auto t1_ALL = std::chrono::high_resolution_clock::now();
  auto t1 = t1_ALL;
  //
  // ------------------------------------------------------------------------------------------------------------------------
  //
  // Declare the histograms for the efficiency
  TH1DMap_t h1D;   // Stores the total and passing histograms separately
  TH2DMap_t h2D;   // Stores the total and passing histograms separately
  //
  // Initialize the efficiencies
  initEff1D(h1D , MU_BIN_);
  initEff2D(h2D , MU_BIN_2D_);
  //
  // ------------------------------------------------------------------------------------------------------------------------
  //
  // Extract all the samples
  std::map< std::string , Long64_t > nentries;
  std::map< std::string , std::unique_ptr< HiMuonTree > > muonTree;
  std::map< std::string , std::unique_ptr< HiMETTree  > > metTree;
  for (const auto & inputFile : inputFileMap_) {
    const std::string sample   = inputFile.first;
    const FileInfo_t  fileInfo = inputFile.second;
    for (const auto& f : fileInfo) { std::cout << "[INFO] Loading File: " << f.first << std::endl; }
    //
    muonTree[sample] = std::unique_ptr<HiMuonTree>(new HiMuonTree());
    if (!muonTree.at(sample)->GetTree(fileInfo)) return;
    nentries[sample] = muonTree.at(sample)->GetEntries();
    //
    metTree[sample] = std::unique_ptr<HiMETTree>(new HiMETTree());
    if (!metTree.at(sample)->GetTree(fileInfo)) return;
    if (metTree.at(sample)->GetEntries() != nentries.at(sample)) { std::cout << "[ERROR] Inconsistent number of entries between MET (" << 
        metTree.at(sample)->GetEntries() << ") and Muon Tree (" << muonTree.at(sample)->GetEntries() << ") !" << std::endl; return; }
  }
  //
  // ------------------------------------------------------------------------------------------------------------------------
  //
  // Loop over the samples
  //
  for (auto & inputFile : inputFileMap_) {
    const std::string sample = inputFile.first;
    // Loop over the events
    for (Long64_t jentry = 0; jentry < nentries.at(sample); jentry++) {
      //
      // Get the entry in the trees
      if (metTree.at(sample)->GetEntry(jentry)<0) break;
      if (muonTree.at(sample)->GetEntry(jentry)<0) break;
      //
      if (jentry%200000==0) std::cout << sample << " : " << jentry << "/" << nentries[sample] << std::endl;
      if (jentry%200000==0) { 
        std::cout << "[INFO] Processing time: " << std::chrono::duration_cast<std::chrono::seconds>( std::chrono::high_resolution_clock::now() - t1 ).count() << " sec" << std::endl;
        t1 = std::chrono::high_resolution_clock::now();
      }
      // 
      // Check that the different tree agrees well
      if (muonTree.at(sample)->Event_Run()    != metTree.at(sample)->Event_Run()   ) { std::cout << "[ERROR] MET Run does not agree!"   << std::endl; return; }
      if (muonTree.at(sample)->Event_Number() != metTree.at(sample)->Event_Number()) { std::cout << "[ERROR] MET Event does not agree!" << std::endl; return; }
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
      const double mcWeight = ( ( muonTree.at(sample)->GetCrossSection() * lumi ) / muonTree.at(sample)->GetEntriesFast() );
      // Set the global weight only in the first event (i.e. once per sample)
      if (jentry==0) {
        setGlobalWeight(h1D, mcWeight, sampleType, col);
        setGlobalWeight(h2D, mcWeight, sampleType, col);
      }
      // Define the event weight (set to 1.0 by default)
      const double evtWeight = 1.0;
      //
      // Define the TnP scale factor (set to 1.0 by default)
      const double sfTnP = 1.0;
      //
      // Check Event Conditions
      //
      // Determine if the event pass the event filters
      const bool passEvent  = PA::passEventFilter(metTree.at(sample));
      // Determine if the event pass the Drell-Yan veto
      const bool passDYVeto = PA::passDrellYanVeto(muonTree.at(sample));
      //
      // Check Muon Conditions
      //
      // Loop over the generated muons
      for (ushort iGenMu = 0; iGenMu < muonTree.at(sample)->Gen_Muon_Mom().size(); iGenMu++) {
        // Check the content of the MC
        if (PA::checkGenMuon(iGenMu, sampleType, muonTree.at(sample)) == false) continue;
        // Check if the generated muon pass the kinematic cuts
        const bool isGoodGenMuon  = PA::selectGenMuon(iGenMu, muonTree.at(sample));
        // Extract the kinematic information of generated muon
        const double mu_Gen_Pt  = muonTree.at(sample)->Gen_Muon_Mom()[iGenMu].Pt();
        const double mu_Gen_Eta = muonTree.at(sample)->Gen_Muon_Mom()[iGenMu].Eta();
        // Determine the charge of the generated muon
        const short charge = muonTree.at(sample)->Gen_Muon_Charge()[iGenMu];
        const std::string chg = ( (charge < 0) ? "Minus" : "Plus" );
        // Fill the VarMap with the kinematic information
        VarMap_t genMuonVar;
        genMuonVar[sampleType][col][chg]["Eta"]    = { mu_Gen_Eta };
        genMuonVar[sampleType][col][chg]["Pt" ]    = { mu_Gen_Pt  };
        genMuonVar[sampleType][col][chg]["Pt_Eta"] = { mu_Gen_Pt , mu_Gen_Eta };
        // Initialize the boolean flags
        bool passIdentification = false;
        bool passTrigger        = false;
        bool passIsolation      = false;
        //
        // Check that the generated muon is within the analysis kinematic range
        if (isGoodGenMuon) {
          // Find the reconstructed muon matched to gen
          const short iPFMu = muonTree.at(sample)->Gen_Muon_PF_Idx()[iGenMu];
          if (iPFMu >= 0) {
            //
            // PF Muon was matched to generated muon
            //
            const short iRecoMu = muonTree.at(sample)->PF_Muon_Reco_Idx()[iPFMu];
            if (iRecoMu < 0) { std::cout << "[ERROR] Reco idx is negative" << std::endl; return; }
            //
            // Check if the reconstructed muon pass muon ID and kinematic cuts
            passIdentification = PA::isGoodMuon(iPFMu , muonTree.at(sample));
            //
            if (passIdentification) {
              // Extract the kinematic information of generated muon
              const double mu_PF_Pt  = muonTree.at(sample)->PF_Muon_Mom()[iPFMu].Pt();
              const double mu_PF_Eta = muonTree.at(sample)->PF_Muon_Mom()[iPFMu].Eta();
              // Fill the VarMap with the kinematic information
              VarMap_t pfMuonVar;
              pfMuonVar[sampleType][col][chg]["Eta"]    = { mu_PF_Eta };
              pfMuonVar[sampleType][col][chg]["Pt" ]    = { mu_PF_Pt  };
              pfMuonVar[sampleType][col][chg]["Pt_Eta"] = { mu_PF_Pt , mu_PF_Eta };
              //
              // Trigger Efficiency (Based on Identified muons)
              //
              // Check if the reconstructed muon is matched to the trigger
              passTrigger = PA::isTriggerMatched(triggerIndex_, iPFMu, muonTree.at(sample));
              //
              if (!fillEff1D(h1D, passTrigger, "Trigger", pfMuonVar, sfTnP, evtWeight)) { return; }
              if (!fillEff2D(h2D, passTrigger, "Trigger", pfMuonVar, sfTnP, evtWeight)) { return; }
              //
              if (passTrigger) {
                //
                // Isolation Efficiency (Based on Identified muons that passed the Trigger)
                //
                // Check if the reconstructed muon pass isolation cuts
                passIsolation = PA::isIsolatedMuon(iPFMu , muonTree.at(sample));
                //
                if (!fillEff1D(h1D, passIsolation, "Isolation", pfMuonVar, sfTnP, evtWeight)) { return; }
                if (!fillEff2D(h2D, passIsolation, "Isolation", pfMuonVar, sfTnP, evtWeight)) { return; }
                //
                if (passIsolation) {
                  //
                  // Event-Based Efficiency (Based on Identified and Isolated muons that passed the Trigger)
                  //
                  if (!fillEff1D(h1D, passEvent, "Event", pfMuonVar, sfTnP, evtWeight)) { return; }
                  if (!fillEff2D(h2D, passEvent, "Event", pfMuonVar, sfTnP, evtWeight)) { return; }
                  if (!fillEff1D(h1D, passDYVeto, "DrellYanVeto", pfMuonVar, sfTnP, evtWeight)) { return; }
                  if (!fillEff2D(h2D, passDYVeto, "DrellYanVeto", pfMuonVar, sfTnP, evtWeight)) { return; }
                }
              }
            }
          }
          //
          // Identification Efficiency (Based on Generated muons)
          //
          if (!fillEff1D(h1D, passIdentification, "Identification", genMuonVar, sfTnP, evtWeight)) { return; }
          if (!fillEff2D(h2D, passIdentification, "Identification", genMuonVar, sfTnP, evtWeight)) { return; }
          //
          if (!fillEff1D(h1D, (passIdentification && passTrigger && passIsolation), "Offline", genMuonVar, sfTnP, evtWeight)) { return; }
          if (!fillEff2D(h2D, (passIdentification && passTrigger && passIsolation), "Offline", genMuonVar, sfTnP, evtWeight)) { return; }
          //
          if (!fillEff1D(h1D, (passIdentification && passTrigger && passIsolation && passEvent && passDYVeto), "Total", genMuonVar, sfTnP, evtWeight)) { return; }
          if (!fillEff2D(h2D, (passIdentification && passTrigger && passIsolation && passEvent && passDYVeto), "Total", genMuonVar, sfTnP, evtWeight)) { return; }
        }
        //
        if (!fillEff1D(h1D, isGoodGenMuon, "Acceptance", genMuonVar, sfTnP, evtWeight)) { return; }
        if (!fillEff2D(h2D, isGoodGenMuon, "Acceptance", genMuonVar, sfTnP, evtWeight)) { return; }
        //
        if (!fillEff1D(h1D, (isGoodGenMuon && passIdentification && passTrigger && passIsolation && passEvent && passDYVeto), "Acceptance_Efficiency", genMuonVar, sfTnP, evtWeight)) { return; }
        if (!fillEff2D(h2D, (isGoodGenMuon && passIdentification && passTrigger && passIsolation && passEvent && passDYVeto), "Acceptance_Efficiency", genMuonVar, sfTnP, evtWeight)) { return; }

      }
    }
  }
  //
  // ------------------------------------------------------------------------------------------------------------------------
  //
  // Declare the efficiencies
  EffMap_t  eff1D; // Stores the efficiency
  EffMap_t  eff2D; // Stores the efficiency
  //
  // Load the efficiencies with the histograms
  if (!loadEff1D(eff1D, h1D)) { return; };
  if (!loadEff2D(eff2D, h2D)) { return; };
  //
  // ------------------------------------------------------------------------------------------------------------------------
  //
  // Merge Efficiencies
  mergeEff(eff1D);
  mergeEff(eff2D);
  //
  // ------------------------------------------------------------------------------------------------------------------------
  //
  // Set Style
  setStyle();
  // Draw the Efficiencies
  drawEff1D(mainDir, eff1D);
  drawEff2D(mainDir, eff2D);
  //
  // ------------------------------------------------------------------------------------------------------------------------
  //
  // Compare 2D histograms
  compareEff2D(mainDir, eff2D, "MC_ZToMuMu_M_30_Inf");
  compareEff2D(mainDir, eff2D, "PA");
  compareEff2D(mainDir, eff2D, "Minus");
  //
  // ------------------------------------------------------------------------------------------------------------------------
  //
  // Store the Efficiencies
  //
  saveEff(mainDir, eff1D, eff2D);
};


void initEff1D(TH1DMap_t& h, const BinMap_t& binMap)
{
  for (const auto& bins : binMap) {
    for (const auto& sample : sampleType_.at("sample")) {
      for (const auto& col : COLL_) {
        for (const auto& chg : CHG_) {
          for (const auto& type : effType) {
            Double_t bin[bins.second.size()];
            for (uint i = 0; i < bins.second.size(); i++) { bin[i] = bins.second[i]; }
            std::get<0>(h[bins.first][sample][col][chg][type]) = TH1D("Passed", "Passed", (bins.second.size()-1), bin);
            std::get<1>(h[bins.first][sample][col][chg][type]) = TH1D("Total" , "Total" , (bins.second.size()-1), bin);
            std::get<2>(h[bins.first][sample][col][chg][type]) = 1.0;
            //
            const std::string name = "eff1D_" + bins.first +"_"+ sample +"_"+ col +"_"+ chg +"_"+ type;
            std::get<0>(h.at(bins.first).at(sample).at(col).at(chg).at(type)).SetName((name+"_Passed").c_str());
            std::get<1>(h.at(bins.first).at(sample).at(col).at(chg).at(type)).SetName((name+"_Total").c_str());
            std::get<0>(h.at(bins.first).at(sample).at(col).at(chg).at(type)).Sumw2(kTRUE);
            std::get<1>(h.at(bins.first).at(sample).at(col).at(chg).at(type)).Sumw2(kTRUE);
          }
        }
      }
    }
  }
};


bool fillEff1D(TH1DMap_t& h, const bool& pass, const std::string& type, const VarMap_t& var, const double sfTnP, const double evtWeight)
{
  const std::string sample = var.begin()->first;
  const std::string col    = var.at(sample).begin()->first;
  const std::string charge = var.at(sample).at(col).begin()->first;
  if (sample.find(col)!=std::string::npos) { std::cout << "[ERROR] Sample name " << sample << " has wrong format!" << std::endl; return false; }
  for (auto& v : h) {
    if (var.at(sample).at(col).at(charge).count(v.first)==0) { std::cout << "[ERROR] Variable " << v.first << " was not loaded properly" << std::endl; return false; }
    for (auto& s : v.second) {
      if (s.first!=sample) continue;
      for (auto& c : s.second) {
        if ((c.first=="PA" && col!="Pbp") || (c.first!="PA" && c.first!=col)) continue;
        for (auto& ch : c.second) {
          if (ch.first!="Inc" && ch.first!=charge) continue;
          if ( (s.first.find("Minus")!=std::string::npos || s.first.find("Plus")!=std::string::npos) && (s.first.find(ch.first)==std::string::npos) ) continue;
          if (ch.second.count(type)==0) { std::cout << "[ERROR] Efficiency type " << type << " is not defined" << std::endl; return false; }
          double xVar = var.at(sample).at(col).at(charge).at(v.first)[0];
          if (v.first=="Eta" && c.first=="PA" && col=="Pbp") { xVar = -xVar; }
          if (xVar >= MU_BIN_.at(v.first)[0] && xVar <= MU_BIN_.at(v.first)[MU_BIN_.at(v.first).size()-1]) { // Don't include values outside of range
            // Fill the Pass histogram
            if (pass) { std::get<0>(ch.second.at(type)).Fill(xVar , (evtWeight*sfTnP)); }
            // Fill the total histogram
            std::get<1>(ch.second.at(type)).Fill(xVar , evtWeight);
          }
        }
      }
    }
  }
  return true;
};


bool loadEff1D(EffMap_t& eff, const TH1DMap_t& h)
{
  for (const auto& v : h) {
    for (const auto& s : v.second) {
      for (const auto& c : s.second) {
        for (const auto& ch : c.second) {
          for (const auto& t : ch.second) {
            const TH1D&   hPassed = std::get<0>(t.second);
            const TH1D&   hTotal  = std::get<1>(t.second);
            const double& weight  = std::get<2>(t.second);
            const std::string name = "eff1D_" + v.first +"_"+ s.first +"_"+ c.first +"_"+ ch.first +"_"+ t.first;
            if ( TEfficiency::CheckConsistency(hPassed, hTotal) ) {
              eff[v.first][s.first][c.first][ch.first][t.first] = TEfficiency(hPassed , hTotal);
              eff.at(v.first).at(s.first).at(c.first).at(ch.first).at(t.first).SetName(name.c_str());
              eff.at(v.first).at(s.first).at(c.first).at(ch.first).at(t.first).SetTitle("");
              // Set Global Weight
              eff.at(v.first).at(s.first).at(c.first).at(ch.first).at(t.first).SetWeight(weight);
              if ( checkWeights(hPassed, hTotal) ) {
                eff.at(v.first).at(s.first).at(c.first).at(ch.first).at(t.first).SetStatisticOption(TEfficiency::kFNormal);
              }
            }
            else {
              std::cout << "[ERROR] 1D Histograms for " << name << " are not consistent!" << std::endl; return false;
            }
          }
        }
      }
    }
  }
  return true;
};


void formatEff1D(TEfficiency& eff, const std::string& var, const std::string& charge, const std::string& type)
{
  // Set the format of all graphs
  if (eff.GetDimension() != 1) { std::cout << "[ERROR] formatEff1D only works for dimension == 1" << std::endl; return; }
  gPad->Update();
  auto graph = eff.GetPaintedGraph();
  if (graph) {
    // General
    graph->SetMarkerColor(kBlue);
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(1.0);
    const std::string objLbl = ( (type.find("Acceptance")!=std::string::npos || type=="Total" || type=="Identification" || type=="Offline") ? "Gen #mu" : "Reco #mu" );
    // X-axis
    std::string xLabel = objLbl; if (charge=="Plus") xLabel += "^{+}"; if (charge=="Minus") xLabel += "^{-}";
    if (var=="Eta") { xLabel += " #eta"; }
    if (var=="Pt" ) { xLabel += " p_{T} (GeV/c)"; }
    graph->GetXaxis()->CenterTitle(kFALSE);
    graph->GetXaxis()->SetTitleOffset(0.9);
    graph->GetXaxis()->SetTitleSize(0.050);
    graph->GetXaxis()->SetLabelSize(0.035);
    graph->GetXaxis()->SetLimits(MU_BIN_.at(var)[0] , MU_BIN_.at(var)[MU_BIN_.at(var).size()-1]);
    // Y-axis
    std::string yLabel = type + " Efficiency";
    if (type == "Acceptance") { yLabel = "Acceptance"; }
    if (type == "Acceptance_Efficiency") { yLabel = "Acceptance x Efficiency"; }
    graph->GetYaxis()->CenterTitle(kFALSE);
    graph->GetYaxis()->SetTitleOffset(1.4);
    graph->GetYaxis()->SetTitleSize(0.04);
    graph->GetYaxis()->SetLabelSize(0.035);
    graph->GetYaxis()->SetRangeUser(0.0, 1.4);
    // Set Axis Titles
    eff.SetTitle(Form(";%s;%s", xLabel.c_str(), yLabel.c_str()));
  }
  gPad->Update(); 
};


void drawEff1D(const std::string& outDir, EffMap_t& effMap)
{
  if (effMap.begin()->second.begin()->second.begin()->second.begin()->second.begin()->second.GetDimension() != 1) { std::cout << "[ERROR] drawEff1D only works for dimension == 1" << std::endl; return; }
  // Draw all graphs
  for (auto& v : effMap) {
    for (auto& s : v.second) {
      for (auto& c : s.second) {
        for (auto& ch : c.second) {
          for (auto& t : ch.second) {
            if (t.second.GetTotalHistogram() && t.second.GetTotalHistogram()->GetEntries()>0) {
              const std::string var    = v.first;
              const std::string sample = s.first;
              const std::string col    = c.first;
              const std::string charge = ch.first;
              const std::string type   = t.first;
              TEfficiency&      eff    = t.second;
              // Create Canvas
              TCanvas c("c", "c", 1000, 1000); c.cd();
              // Create the Text Info
              TLatex tex; tex.SetNDC(); tex.SetTextSize(0.028); float dy = 0;
              std::vector< std::string > textToPrint;
              std::string sampleLabel; formatDecayLabel(sampleLabel, sample);
              textToPrint.push_back(sampleLabel);
              if (var=="Eta") textToPrint.push_back("p^{#mu}_{T} > 25 GeV/c");
              if (var=="Pt" ) textToPrint.push_back("|#eta^{#mu}| < 2.4");
              // Draw graph
              eff.Draw("AP");
              c.Modified(); c.Update();
              // Format the Graph
              formatEff1D(eff, var, charge, type);
              c.Modified(); c.Update();
              // Add Min, Max and Mean value of efficiency
              auto graph = eff.GetPaintedGraph();
              if (graph) {
                const double min = TMath::MinElement(graph->GetN(), graph->GetY())*100.;
                const double max = TMath::MaxElement(graph->GetN(), graph->GetY())*100.;
                const double avg = graph->GetMean(2)*100.;
                const double err = graph->GetRMS(2)*100.;
                textToPrint.push_back(Form("min: %.2f %% , max: %.2f %%", min, max));
                textToPrint.push_back(Form("mean: %.2f #pm %.2f %%", avg, err));
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
              const std::string plotDir = outDir + "Efficiency1D/" + var+"/" + sample+"/" + col;
              makeDir(plotDir + "/png/");
              makeDir(plotDir + "/pdf/");
              makeDir(plotDir + "/root/");
              // Save Canvas
              c.SaveAs(( plotDir + "/png/" + eff.GetName() + ".png" ).c_str());
              c.SaveAs(( plotDir + "/pdf/" + eff.GetName() + ".pdf" ).c_str());
              c.SaveAs(( plotDir + "/root/" + eff.GetName() + ".root" ).c_str());
              // Clean up memory
              c.Clear(); c.Close();
            }
          }
        }
      }
    }
  }
};


void initEff2D(TH2DMap_t& h, const Bin2DMap_t& binMap)
{
  for (const auto& bins : binMap) {
    for (const auto& sample : sampleType_.at("sample")) {
      for (const auto& col : COLL_) {
        for (const auto& chg : CHG_) {
          for (const auto& type : effType) {
            Double_t binX[bins.second.first.size()];
            for (uint i = 0; i < bins.second.first.size(); i++) { binX[i] = bins.second.first[i]; }
            Double_t binY[bins.second.second.size()];
            for (uint i = 0; i < bins.second.second.size(); i++) { binY[i] = bins.second.second[i]; }

            std::get<0>(h[bins.first][sample][col][chg][type]) = TH2D("Passed", "Passed", (bins.second.first.size()-1), binX, (bins.second.second.size()-1), binY);
            std::get<1>(h[bins.first][sample][col][chg][type]) = TH2D("Total" , "Total" , (bins.second.first.size()-1), binX, (bins.second.second.size()-1), binY);
            std::get<2>(h[bins.first][sample][col][chg][type]) = 1.0;
            //
            const std::string name = "eff2D_" + bins.first +"_"+ sample +"_"+ col +"_"+ chg +"_"+ type;
            std::get<0>(h.at(bins.first).at(sample).at(col).at(chg).at(type)).SetName((name+"_Passed").c_str());
            std::get<1>(h.at(bins.first).at(sample).at(col).at(chg).at(type)).SetName((name+"_Total").c_str());
            std::get<0>(h.at(bins.first).at(sample).at(col).at(chg).at(type)).Sumw2(kTRUE);
            std::get<1>(h.at(bins.first).at(sample).at(col).at(chg).at(type)).Sumw2(kTRUE);
          }
        }
      }
    }
  }
};


bool fillEff2D(TH2DMap_t& h, const bool& pass, const std::string& type, const VarMap_t& var, const double sfTnP, const double evtWeight)
{
  const std::string sample = var.begin()->first;
  const std::string col    = var.at(sample).begin()->first;
  const std::string charge = var.at(sample).at(col).begin()->first;
  if (sample.find(col)!=std::string::npos) { std::cout << "[ERROR] Sample name " << sample << " has wrong format!" << std::endl; return false; }
  for (auto& v : h) {
    if (var.at(sample).at(col).at(charge).count(v.first)==0) { std::cout << "[ERROR] Variable " << v.first << " was not loaded properly" << std::endl; return false; }
    for (auto& s : v.second) {
      if (s.first!=sample) continue;
      for (auto& c : s.second) {
        if ((c.first=="PA" && col!="Pbp") || (c.first!="PA" && c.first!=col)) continue;
        for (auto& ch : c.second) {
          if (ch.first!="Inc" && ch.first!=charge) continue;
          if ( (s.first.find("Minus")!=std::string::npos || s.first.find("Plus")!=std::string::npos) && (s.first.find(ch.first)==std::string::npos) ) continue;
          if (ch.second.count(type)==0) { std::cout << "[ERROR] Efficiency type " << type << " is not defined" << std::endl; return false; }
          double xVar = var.at(sample).at(col).at(charge).at(v.first)[0];
          double yVar = var.at(sample).at(col).at(charge).at(v.first)[1];
          if (v.first=="Pt_Eta" && c.first=="PA" && col=="Pbp") { yVar = -yVar; }
          if (
              (xVar >= MU_BIN_2D_.at(v.first).first[0] && xVar <= MU_BIN_2D_.at(v.first).first[MU_BIN_2D_.at(v.first).first.size()-1]) &&
              (yVar >= MU_BIN_2D_.at(v.first).second[0] && yVar <= MU_BIN_2D_.at(v.first).second[MU_BIN_2D_.at(v.first).second.size()-1])
              ) { // Don't include values outside of range
            // Fill the Pass histogram
            if (pass) { std::get<0>(ch.second.at(type)).Fill(xVar , yVar , (evtWeight*sfTnP)); }
            // Fill the total histogram
            std::get<1>(ch.second.at(type)).Fill(xVar , yVar , evtWeight);
          }
        }
      }
    }
  }
  return true;
};


bool loadEff2D(EffMap_t& eff, const TH2DMap_t& h)
{
  for (const auto& v : h) {
    for (const auto& s : v.second) {
      for (const auto& c : s.second) {
        for (const auto& ch : c.second) {
          for (const auto& t : ch.second) {
            const TH2D&   hPassed = std::get<0>(t.second);
            const TH2D&   hTotal  = std::get<1>(t.second);
            const double& weight  = std::get<2>(t.second);
            const std::string name = "eff2D_" + v.first +"_"+ s.first +"_"+ c.first +"_"+ ch.first +"_"+ t.first;
            if ( TEfficiency::CheckConsistency(hPassed, hTotal) ) {
              eff[v.first][s.first][c.first][ch.first][t.first] = TEfficiency(hPassed , hTotal);
              eff.at(v.first).at(s.first).at(c.first).at(ch.first).at(t.first).SetName(name.c_str());
              eff.at(v.first).at(s.first).at(c.first).at(ch.first).at(t.first).SetTitle("");
              // Set Global Weight
              eff.at(v.first).at(s.first).at(c.first).at(ch.first).at(t.first).SetWeight(weight);
              if ( checkWeights(hPassed, hTotal) ) {
                eff.at(v.first).at(s.first).at(c.first).at(ch.first).at(t.first).SetStatisticOption(TEfficiency::kFNormal);
              }
            }
            else {
              std::cout << "[ERROR] 2D Histograms for " << name << " are not consistent!" << std::endl; return false;
            }
          }
        }
      }
    }
  }
  return true;
};


void formatEff2D(TEfficiency& eff, const std::string& var, const std::string& charge, const std::string& type)
{
  // Set the format of all histograms
  if (eff.GetDimension() != 2) { std::cout << "[ERROR] formatEff2D only works for dimension == 2" << std::endl; return; }
  gPad->Update();
  auto histo = eff.GetPaintedHistogram();
  if (histo) {
    // General
    const std::string objLbl = ( (type.find("Acceptance")!=std::string::npos || type=="Total" || type=="Identification" || type=="Offline") ? "Gen #mu" : "Reco #mu" );
    // X-axis
    std::string xVar = ""; if (var=="Pt_Eta") { xVar = "Pt"; }
    std::string xLabel = objLbl; if (charge=="Plus") xLabel += "^{+}"; if (charge=="Minus") xLabel += "^{-}";
    if (xVar=="Eta") { xLabel += " #eta"; }
    if (xVar=="Pt" ) { xLabel += " p_{T} (GeV/c)"; }
    histo->GetXaxis()->CenterTitle(kFALSE);
    histo->GetXaxis()->SetTitleOffset(1.2);
    histo->GetXaxis()->SetTitleSize(0.043);
    histo->GetXaxis()->SetLabelSize(0.035);
    const double xMin = MU_BIN_.at(xVar)[0];
    const double xMax = MU_BIN_.at(xVar)[MU_BIN_.at(xVar).size()-1];
    histo->GetXaxis()->SetRangeUser(xMin , xMax);
    // Y-axis
    std::string yVar = ""; if (var=="Pt_Eta") { yVar = "Eta"; }
    std::string yLabel = objLbl; if (charge=="Plus") yLabel += "^{+}"; if (charge=="Minus") yLabel += "^{-}";
    if (yVar=="Eta") { yLabel += " #eta"; }
    if (yVar=="Pt" ) { yLabel += " p_{T} (GeV/c)"; }
    histo->GetYaxis()->CenterTitle(kFALSE);
    histo->GetYaxis()->SetTitleOffset(1.1);
    histo->GetYaxis()->SetTitleSize(0.045);
    histo->GetYaxis()->SetLabelSize(0.035);
    const double yMin = MU_BIN_.at(yVar)[0];
    const double yMax = MU_BIN_.at(yVar)[MU_BIN_.at(yVar).size()-1] + (MU_BIN_.at(yVar)[MU_BIN_.at(yVar).size()-1] - MU_BIN_.at(yVar)[0])*0.3;
    histo->GetYaxis()->SetLimits(yMin , yMax);
    histo->SetAxisRange(yMin , yMax , "Y");
    // Z-axis
    histo->GetZaxis()->SetLabelSize(0.025);
    histo->GetZaxis()->SetRangeUser(0.0 , 1.0);
    auto palette = (TPaletteAxis*)histo->GetListOfFunctions()->FindObject("palette");
    palette->SetX1NDC(0.91);
    palette->SetX2NDC(0.94);
    palette->SetY1NDC(0.12);
    palette->SetY2NDC(0.92);
    // Set Axis Titles
    eff.SetTitle(Form(";%s;%s;", xLabel.c_str(), yLabel.c_str()));
  }
  gPad->Update();
};


void drawEff2D(const std::string& outDir, EffMap_t& effMap)
{
  if (effMap.begin()->second.begin()->second.begin()->second.begin()->second.begin()->second.GetDimension() != 2) { std::cout << "[ERROR] drawEff2D only works for dimension == 2" << std::endl; return; }
  // Draw all graphs
  for (auto& v : effMap) {
    for (auto& s : v.second) {
      for (auto& c : s.second) {
        for (auto& ch : c.second) {
          for (auto& t : ch.second) {
            if (t.second.GetTotalHistogram() && t.second.GetTotalHistogram()->GetEntries()>0) {
              const std::string var    = v.first;
              const std::string sample = s.first;
              const std::string col    = c.first;
              const std::string charge = ch.first;
              const std::string type   = t.first;
              TEfficiency&      eff    = t.second;
              // Create Canvas
              TCanvas c("c", "c", 1000, 1000); c.cd();
              c.SetRightMargin(1.5);
              // Create the Text Info
              TLatex tex; tex.SetNDC(); tex.SetTextSize(0.028); float dy = 0;
              std::vector< std::string > textToPrint;
              std::string zLabel = type + " Efficiency";
              if (type == "Acceptance") { zLabel = "Acceptance"; }
              if (type == "Acceptance_Efficiency") { zLabel = "Acceptance x Efficiency"; }
              textToPrint.push_back(zLabel);
              std::string sampleLabel; formatDecayLabel(sampleLabel, sample);
              textToPrint.push_back(sampleLabel);
              // Draw graph
              eff.Draw("colz");
              c.Modified(); c.Update();
              // Format the Graph
              formatEff2D(eff, var, charge, type);
              c.Modified(); c.Update();
              // Add Min, Max and Mean value of efficiency
              auto histo = eff.GetPaintedHistogram();
              if (histo) {
                double avg = 0. , err = 0. , min = 0. , max = 0.;
                getStats2D(avg, err, min, max, *histo);
                textToPrint.push_back(Form("min: %.2f %% , max: %.2f %%", (min*100.), (max*100.)));
                textToPrint.push_back(Form("mean: %.2f #pm %.2f %%", (avg*100.), (err*100.)));
                if (min==max) { min -= 0.01; max += 0.01; } if (max>1.0) { max = 1.0; }
                histo->GetZaxis()->SetRangeUser( (std::floor(min*10.)/10.) , (std::ceil(max*10.)/10.) );
              }
              // Draw the text
              for (const auto& s: textToPrint) { tex.DrawLatex(0.22, 0.87-dy, s.c_str()); dy+=0.04; }
              c.Modified(); c.Update();
              // Draw the efficiency relative errors with text
              TH2D hist;                        
              fillErrorEff2D(hist, eff, var);
              hist.SetMarkerSize(0.7);
              hist.Draw("SAME,TEXT");
              c.Modified(); c.Update();
              // set the CMS style
              int option = 114;
              if (col.find("pPb")!=std::string::npos) option = 112;
              if (col.find("Pbp")!=std::string::npos) option = 113;
              CMS_lumi(&c, option, 33, "");
              c.Modified(); c.Update();
              // Create Output Directory
              const std::string plotDir = outDir + "Efficiency2D/" + var+"/" + sample+"/" + col;
              makeDir(plotDir + "/png/");
              makeDir(plotDir + "/pdf/");
              makeDir(plotDir + "/root/");
              // Save Canvas
              c.SaveAs(( plotDir + "/png/" + eff.GetName() + ".png" ).c_str());
              c.SaveAs(( plotDir + "/pdf/" + eff.GetName() + ".pdf" ).c_str());
              c.SaveAs(( plotDir + "/root/" + eff.GetName() + ".root" ).c_str());
              // Clean up memory
              c.Clear(); c.Close();
            }
          }
        }
      }
    }
  }
};


void formatComparisonEff2D(TH2& histo, const std::string& var, const std::string& charge, const std::string& type)
{
  // Set the format of all histograms
  // General
  const std::string objLbl = ( (type.find("Acceptance")!=std::string::npos || type=="Total" || type=="Identification" || type=="Offline") ? "Gen #mu" : "Reco #mu" );
  // X-axis
  std::string xVar = ""; if (var=="Pt_Eta") { xVar = "Pt"; }
  std::string xLabel = objLbl; if (charge=="Plus") xLabel += "^{+}"; if (charge=="Minus") xLabel += "^{-}";
  if (xVar=="Eta") { xLabel += " #eta"; }
  if (xVar=="Pt" ) { xLabel += " p_{T} (GeV/c)"; }
  histo.GetXaxis()->CenterTitle(kFALSE);
  histo.GetXaxis()->SetTitleOffset(1.2);
  histo.GetXaxis()->SetTitleSize(0.043);
  histo.GetXaxis()->SetLabelSize(0.035);
  const double xMin = MU_BIN_.at(xVar)[0];
  const double xMax = MU_BIN_.at(xVar)[MU_BIN_.at(xVar).size()-1];
  histo.GetXaxis()->SetRangeUser(xMin , xMax);
  // Y-axis
  std::string yVar = ""; if (var=="Pt_Eta") { yVar = "Eta"; }
  std::string yLabel = objLbl; if (charge=="Plus") yLabel += "^{+}"; if (charge=="Minus") yLabel += "^{-}";
  if (yVar=="Eta") { yLabel += " #eta"; }
  if (yVar=="Pt" ) { yLabel += " p_{T} (GeV/c)"; }
  histo.GetYaxis()->CenterTitle(kFALSE);
  histo.GetYaxis()->SetTitleOffset(1.1);
  histo.GetYaxis()->SetTitleSize(0.045);
  histo.GetYaxis()->SetLabelSize(0.035);
  const double yMin = MU_BIN_.at(yVar)[0];
  const double yMax = MU_BIN_.at(yVar)[MU_BIN_.at(yVar).size()-1] + (MU_BIN_.at(yVar)[MU_BIN_.at(yVar).size()-1] - MU_BIN_.at(yVar)[0])*0.3;
  histo.GetYaxis()->SetLimits(yMin , yMax);
  histo.SetAxisRange(yMin , yMax , "Y");
  // Z-axis
  histo.GetZaxis()->SetLabelSize(0.022);
  auto palette = (TPaletteAxis*)histo.GetListOfFunctions()->FindObject("palette");
  if (palette==NULL) { std::cout << "[WARNING] Palette is NULL!" << std::endl; return; }
  palette->SetX1NDC(0.91);
  palette->SetX2NDC(0.93);
  palette->SetY1NDC(0.12);
  palette->SetY2NDC(0.92);
  // Set Axis Titles
  histo.SetTitle(Form(";%s;%s;", xLabel.c_str(), yLabel.c_str()));
};


void compareEff2D(const std::string& outDir, EffMap_t& effMap, const std::string comp)
{
  if (comp=="") { std::cout << "[ERROR] Comparison string is empty" << std::endl; return; }
  if (effMap.begin()->second.begin()->second.begin()->second.begin()->second.begin()->second.GetDimension() != 2) { std::cout << "[ERROR] compareEff2D only works for dimension == 2" << std::endl; return; }
  std::string compLabel; formatDecayLabel(compLabel, comp);
  // Draw all graphs
  for (auto& v : effMap) {
    for (auto& s : v.second) {
      if (s.first==comp) continue;
      for (auto& c : s.second) {
        if (c.first==comp) continue;
        for (auto& ch : c.second) {
          if (ch.first==comp) continue;
          for (auto& t : ch.second) {
            if (t.second.GetTotalHistogram() && t.second.GetTotalHistogram()->GetEntries()>0) {
              const std::string var    = v.first;
              const std::string sample = s.first;
              const std::string col    = c.first;
              const std::string charge = ch.first;
              const std::string type   = t.first;
              TEfficiency&      eff    = t.second;
              // Find the Reference 2D Eff
              std::string sample_Ref = sample, col_Ref = col, charge_Ref = charge;
              if (comp.find("MC_")!=std::string::npos) { sample_Ref = comp; }
              else if (comp=="Minus"||comp=="Plus") { charge_Ref = comp; }
              else if (comp=="PA") { col_Ref = comp; }
              else { std::cout << "[ERROR] Comparison string " << comp << " was not found!" << std::endl; return; }
              TEfficiency& ref = effMap.at(var).at(sample_Ref).at(col_Ref).at(charge_Ref).at(type);
              if (ref.GetTotalHistogram()==NULL) { std::cout << "[ERROR] Reference was not found!" << std::endl; return; }
              if (ref.GetTotalHistogram()->GetEntries()<10) continue;
              // Set the operation for the reference
              std::string opr = "-"; if (comp=="Minus"||comp=="Plus") { opr = "/"; }
              // Create Canvas
              TCanvas c("c", "c", 1000, 1000); c.cd();
              c.SetRightMargin(2.6);
              // Create the Text Info
              TLatex tex; tex.SetNDC(); tex.SetTextSize(0.028); float dy = 0;
              std::vector< std::string > textToPrint;
              std::string zLabel = type + " Efficiency";
              if (type == "Acceptance") { zLabel = "Acceptance"; }
              if (type == "Acceptance_Efficiency") { zLabel = "Acceptance x Efficiency"; }
              if (opr=="-") { textToPrint.push_back(zLabel + " Difference"); }
              if (opr=="/") { textToPrint.push_back(zLabel + " Ratio"); }
              std::string sampleLabel; formatDecayLabel(sampleLabel, sample);
              textToPrint.push_back(sampleLabel);
              if (comp.find("MC_")!=std::string::npos) { textToPrint.push_back("Ref: "+compLabel); }
              else { textToPrint.push_back("Ref: "+comp); }
              // Draw compison between histograms
              eff.Draw(); ref.Draw();
              gPad->Update();
              auto effHist = eff.GetPaintedHistogram();
              if (effHist==NULL) { std::cout << "[ERROR] Efficiency histogram is NULL!" << std::endl; return; }
              auto refHist = ref.GetPaintedHistogram();
              if (refHist==NULL) { std::cout << "[ERROR] Reference histogram is NULL!" << std::endl; return; }
              TH2D cmpHist(*((TH2D*)effHist)); cmpHist.Sumw2(kTRUE);
              if (opr=="-") { cmpHist.Add(refHist, -1); }
              if (opr=="/") { cmpHist.Divide(refHist);  }
              // Set Name
              std::string name = "";
              if (opr=="-") {
                if (comp.find("MC_")!=std::string::npos) { name = "dEff2D_MC_" + var +"_"+ sample +"_"+ col +"_"+ charge +"_"+ type; }
                else { name = "dEff2D_PA_" + var +"_"+ sample +"_"+ col +"_"+ charge +"_"+ type; }
              }
              if (opr=="/") { name = "rEff2D_MuChg_" + var +"_"+ sample +"_"+ col +"_"+ charge +"_"+ type; }
              cmpHist.SetName(name.c_str());
              // Add Min, Max and Mean value of the histogram
              double avg = 0. , err = 0. , min = 0. , max = 0.;
              getStats2D(avg, err, min, max, cmpHist);
              if (opr=="/") {
                textToPrint.push_back(Form("min: %.2f , max: %.2f", min, max));
                textToPrint.push_back(Form("mean: %.2f #pm %.2f", avg, err));
              }
              else {
                textToPrint.push_back(Form("min: %.2f %% , max: %.2f %%", (min*100.), (max*100.)));
                textToPrint.push_back(Form("mean: %.2f #pm %.2f %%", (avg*100.), (err*100.)));
              }
              // Draw the histogram
              cmpHist.Draw("colz");
              c.Modified(); c.Update();
              // Format the hisstogram
              formatComparisonEff2D(cmpHist, var, charge, type);
              const double refVal = ( (opr=="/") ? 1.0 : 0.0 );
              if (min==max) { min -= 0.01; max += 0.01; }
              if (min > refVal) { cmpHist.GetZaxis()->SetRangeUser( (std::floor(min*100.)/100.) , (std::ceil(max*100.)/100.) ); }
              else {
                const double range = std::max( std::abs(refVal - (std::floor(min*100.)/100.)) , std::abs(refVal - (std::ceil(max*100.)/100.)) );
                cmpHist.GetZaxis()->SetRangeUser( (refVal - range) , (refVal + range) );
              }
              c.Modified(); c.Update();
              // Draw the white box to hide under/overflow bins
              const double yMax = MU_BIN_2D_.at(var).second[MU_BIN_2D_.at(var).second.size()-1];
              TBox box(c.GetFrame()->GetX1(), yMax, c.GetFrame()->GetX2(), c.GetFrame()->GetY2());
              box.SetFillColor(kWhite); box.SetLineColor(kWhite);
              box.Draw("same");
              redrawBorder();
              c.Modified(); c.Update();
              // Draw the text
              for (const auto& s: textToPrint) { tex.DrawLatex(0.20, 0.87-dy, s.c_str()); dy+=0.035;}
              c.Modified(); c.Update();
              // set the CMS style
              int option = 114;
              if (col.find("pPb")!=std::string::npos) option = 112;
              if (col.find("Pbp")!=std::string::npos) option = 113;
              CMS_lumi(&c, option, 33, "");
              c.Modified(); c.Update();
              // Create Output Directory
              const std::string plotDir = outDir + "Comparison2D/" + comp+"/" + var+"/" + sample+"/" + col;
              makeDir(plotDir + "/png/");
              makeDir(plotDir + "/pdf/");
              makeDir(plotDir + "/root/");
              // Save Canvas
              c.SaveAs(( plotDir + "/png/" + cmpHist.GetName() + ".png" ).c_str());
              c.SaveAs(( plotDir + "/pdf/" + cmpHist.GetName() + ".pdf" ).c_str());
              c.SaveAs(( plotDir + "/root/" + cmpHist.GetName() + ".root" ).c_str());
              // Clean up memory
              c.Clear(); c.Close();
            }
          }
        }
      }
    }
  }
};
  

void mergeEff(EffMap_t& eff)
{
  // Step 1: Merge pPb and Pbp (inverted) -> PA
  if(std::find(COLL_.begin(), COLL_.end(), "PA") != COLL_.end()) {
    for (auto& v : eff) {
      for (auto& s : v.second) {
        if (s.second.count("PA")>0 && s.second.count("pPb")>0) {
          for (auto& ch : s.second.at("PA")) {
            for (auto& t : ch.second) {
              // Just add the pPb, no need to combine
              const TEfficiency& eff_pPb = s.second.at("pPb").at(ch.first).at(t.first);
              const TEfficiency& eff_Pbp = s.second.at("Pbp").at(ch.first).at(t.first);
              // Passed Histogram
              TH1D hPassed = *((TH1D*)eff_pPb.GetPassedHistogram());
              hPassed.Add(eff_pPb.GetPassedHistogram(), t.second.GetPassedHistogram(), eff_pPb.GetWeight(), eff_Pbp.GetWeight());
              t.second.SetPassedHistogram(hPassed, "f");
              // Total Histogram
              TH1D hTotal = *((TH1D*)eff_pPb.GetTotalHistogram());
              hTotal.Add(eff_pPb.GetTotalHistogram(), t.second.GetTotalHistogram(), eff_pPb.GetWeight(), eff_Pbp.GetWeight());
              t.second.SetTotalHistogram(hTotal, "f");
              // Set the statistics in case of weighted histograms
              if ( checkWeights(hPassed, hTotal) ) { t.second.SetStatisticOption(TEfficiency::kFNormal); }
            }
          }
        }
        else {
          std::cout << "[WARNING] Can't merge pPb and Pbp in " << s.first << " ! "<< std::endl;
          if (s.second.count("PA")>0) { s.second.erase("PA"); }
        }
      }
    }
  }
  // Step 2: Merge W Powheg samples
  const std::vector< std::string > sampleVec = { "MC_WToMuNu" , "MC_WToTauNu" };
  for (const auto& sample : sampleVec) {
    for (auto& v : eff) {
      if (v.second.count(sample+"_Minus")>0 && v.second.count(sample+"_Plus")>0) {
        for (auto& c : v.second.at(sample+"_Minus")) {
          for (auto& ch : c.second) {
            for (auto& t : ch.second) {
              std::string name = "";
              if (t.second.GetDimension()==1) { name = "eff1D_" + v.first +"_"+ sample +"_"+ c.first +"_"+ ch.first +"_"+ t.first; }
              if (t.second.GetDimension()==2) { name = "eff2D_" + v.first +"_"+ sample +"_"+ c.first +"_"+ ch.first +"_"+ t.first; }
              eff.at(v.first)[sample][c.first][ch.first][t.first] = eff.at(v.first).at(sample+"_"+ch.first).at(c.first).at(ch.first).at(t.first);
              eff.at(v.first).at(sample).at(c.first).at(ch.first).at(t.first).SetName(name.c_str());
            }
          }
        }
        eff.at(v.first).erase(sample+"_Minus");
        eff.at(v.first).erase(sample+"_Plus");
      }
    }
  }
};


void writeEff(TFile& file, const EffMap_t& eff, const std::string& mainDirName)
{
  if (eff.size()>0) {
    TDirectory* mainDir = file.mkdir(mainDirName.c_str());
    mainDir->cd();
    for (auto& v : eff) {
      TDirectory* varDir = mainDir->mkdir(v.first.c_str());
      varDir->cd();
      for (auto& s : v.second) {
        TDirectory* sampleDir = varDir->mkdir(s.first.c_str());
        sampleDir->cd();
        for (auto& c : s.second) {
          TDirectory* colDir = sampleDir->mkdir(c.first.c_str());
          colDir->cd();
          for (auto& ch : c.second) {
            for (auto& t : ch.second) {
              t.second.Write(t.second.GetName());
            }
          }
          sampleDir->cd();
        }
        varDir->cd();
      }
      mainDir->cd();
    }
  }
};


void saveEff(const std::string& outDir, const EffMap_t& eff1D, const EffMap_t& eff2D)
{
  // Step 0: Create the output file
  const std::string fileName = "efficiency.root";
  TFile file((outDir + fileName).c_str(), "RECREATE");
  file.cd();
  // Step 1: Store all the 1D Efficiency objects
  writeEff(file, eff1D, "Efficiency1D");
  // Step 2: Store all the 2D Efficiency objects
  writeEff(file, eff2D, "Efficiency2D");
  // Step 3: Write and Close the file
  file.Write();
  file.Close("R");
};


void setGlobalWeight(TH1DMap_t& h, const double& weight, const std::string& sample, const std::string& col)
{
  for (auto& v : h) {
    if (v.second.count(sample)>0 && v.second.at(sample).count(col)>0) {
      for (auto& ch : v.second.at(sample).at(col)) {
        for (auto& t : ch.second) {
          std::get<2>(t.second) = weight;
        }
      }
    }
  }
};


void setGlobalWeight(TH2DMap_t& h, const double& weight, const std::string& sample, const std::string& col)
{
  for (auto& v : h) {
    if (v.second.count(sample)>0 && v.second.at(sample).count(col)>0) {
      for (auto& ch : v.second.at(sample).at(col)) {
        for (auto& t : ch.second) {
          std::get<2>(t.second) = weight;
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
};


void formatLegendEntry(TLegendEntry& e)
{
  e.SetTextSize(0.028);
};


void formatDecayLabel(std::string& label, const std::string& inLabel)
{
  label = "";
  if (inLabel.find("MC_DYToMuMu")!=std::string::npos) {
    label = "Z/#gamma* #rightarrow #mu^{+} + #mu^{-}";
    if (inLabel.find("M_10_30")!=std::string::npos) { label = "Z/#gamma* #rightarrow #mu^{+} + #mu^{-} [10 < M < 30]"; }
    if (inLabel.find("M_30_")  !=std::string::npos) { label = "Z/#gamma* #rightarrow #mu^{+} + #mu^{-} [M > 30] ";      }
  }
  if (inLabel.find("MC_ZToMuMu")!=std::string::npos) {
    label = "Z #rightarrow #mu^{+} + #mu^{-}";
    if (inLabel.find("M_10_30")!=std::string::npos) { label = "Z #rightarrow #mu^{+} + #mu^{-} [10 < M < 30]"; }
    if (inLabel.find("M_30_")  !=std::string::npos) { label = "Z #rightarrow #mu^{+} + #mu^{-} [M > 30] ";      }
  }
  if (inLabel.find("MC_WToMu")!=std::string::npos) {
     label = "W #rightarrow #mu + #nu_{#mu}";
    if (inLabel.find("Plus") !=std::string::npos) { label = "W^{+} #rightarrow #mu^{+} + #nu_{#mu}"; }
    if (inLabel.find("Minus")!=std::string::npos) { label = "W^{-} #rightarrow #mu^{-} + #bar{#nu}_{#mu}"; }
  }
  if (inLabel.find("MC_WToTau")!=std::string::npos) {
     label = "W #rightarrow #tau #rightarrow #mu + #nu_{#mu} + #bar{#nu}_{#tau}";
    if (inLabel.find("Plus") !=std::string::npos) { label = "W^{+} #rightarrow #tau^{+} #rightarrow #mu^{+} + #nu_{#mu} + #bar{#nu}_{#tau}"; }
    if (inLabel.find("Minus")!=std::string::npos) { label = "W^{-} #rightarrow #tau^{-} #rightarrow #mu^{-} + #bar{#nu}_{#mu} + #nu_{#tau}"; }
  }
  if (inLabel == "MC_QCDToMu") { label = "QCD #rightarrow #mu"; }
  if (inLabel == "MC_TTall") { label = "t + #bar{t} #rightarrow All #rightarrow #mu"; }
};


void getStats2D(double& mean, double& rms, double& min, double& max, const TH2& hist)
{
  const uint nBins = (hist.GetNbinsX() * hist.GetNbinsY());
  // Get Mean Value
  mean = (hist.GetSumOfWeights() / nBins);
  // Get RMS, Min and Max Value
  double hMin = 999999999999., hMax = -999999999999.;
  for (int biny = 1; biny <= hist.GetNbinsY(); biny++) {
    for (int binx = 1; binx <= hist.GetNbinsX(); binx++) {
      const double value = hist.GetBinContent(hist.GetBin(binx, biny));
      rms += (value - mean)*(value - mean);
      if (value > hMax) { hMax = value; }
      if (value < hMin) { hMin = value; }
    }
  }
  rms = (rms / nBins);
  rms = std::sqrt(rms);
  min = hMin;
  max = hMax;
};


void fillErrorEff2D(TH2D& hist, const TEfficiency& eff, const std::string& var)
{
  if (eff.GetDimension() != 2) { std::cout << "[ERROR] fillErrorEff2D only works for dimension == 2" << std::endl; return; }
  const uint xbin = 3;
  const uint xIni = ( MU_BIN_2D_.at(var).first.size() - xbin );
  Double_t binX[xbin];
  for (uint i = 0; i < xbin; i++) { binX[i] = MU_BIN_2D_.at(var).first[i + xIni]; }
  Double_t binY[MU_BIN_2D_.at(var).second.size()];
  for (uint i = 0; i < MU_BIN_2D_.at(var).second.size(); i++) { binY[i] = MU_BIN_2D_.at(var).second[i]; }
  hist = TH2D("tmp","", (xbin-1), binX, (MU_BIN_2D_.at(var).second.size()-1), binY);
  for (int biny = 1; biny <= hist.GetNbinsY(); biny++) {
    for (int binx = 1; binx <= hist.GetNbinsX(); binx++) {
      const int effBinx = (binx + xIni);
      const int effBiny = biny;
      const int effBin = eff.GetGlobalBin(effBinx, effBiny);
      const double effValue = eff.GetEfficiency(effBin);
      const double effErrorLow = eff.GetEfficiencyErrorLow(effBin);
      const double effErrorUp = eff.GetEfficiencyErrorUp(effBin);
      const double effError = std::max(std::abs(effErrorLow), std::abs(effErrorUp));
      const double effRelError = ((effError*100.)/std::abs(effValue));
      hist.SetBinContent(binx, biny, effRelError);
    }
  }
};


void redrawBorder()
{
   gPad->Update();
   gPad->RedrawAxis();
   TLine l;
   l.DrawLine(gPad->GetUxmin(), gPad->GetUymax(), gPad->GetUxmax(), gPad->GetUymax());
   l.DrawLine(gPad->GetUxmax(), gPad->GetUymin(), gPad->GetUxmax(), gPad->GetUymax());
};


bool checkWeights(const TH1& pass , const TH1& total)
{
  if (pass.GetSumw2N() == 0 && total.GetSumw2N() == 0) return false;
 
  // check also that the total sum of weight and weight squares are consistent
  Double_t statpass[TH1::kNstat];
  Double_t stattotal[TH1::kNstat];
 
  pass.GetStats(statpass);
  total.GetStats(stattotal);
 
  double tolerance = (total.IsA() == TH1F::Class() ) ? 1.E-5 : 1.E-12; 
    
  //require: sum of weights == sum of weights^2
  if(!TMath::AreEqualRel(statpass[0],statpass[1],tolerance) ||
     !TMath::AreEqualRel(stattotal[0],stattotal[1],tolerance) ) {
    return true;
  }

  // histograms are not weighted 
  return false;  
};

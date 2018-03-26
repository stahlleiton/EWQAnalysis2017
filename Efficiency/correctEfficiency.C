#if !defined(__CINT__) || defined(__MAKECINT__)
// Auxiliary Headers
#include "../Utilities/HiMETTree.h"
#include "../Utilities/HiMuonTree.h"
#include "../Utilities/HiEvtTree.h"
#include "../Utilities/tnp_weight.h"
#include "../Utilities/HFweight.h"
#include "../Utilities/RoccoR.cc"
#include "../Utilities/EVENTUTILS.h"
// ROOT headers
#include "TROOT.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH1D.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TFrame.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TPaletteAxis.h"
#include "TVectorD.h"
#include "TMath.h"
#include "TRandom3.h"
// c++ headers
#include <dirent.h>
#include <iostream>
#include <map>
#include <vector>
#include <tuple>
#include <string>
// CMS headers
#include "../Utilities/CMS/tdrstyle.C"
#include "../Utilities/CMS/CMS_lumi.C"

#endif


// ------------------ TYPE -------------------------------
using TnPVec_t     =  std::map< std::string , std::vector< double > >;
using Unc1DVec_t   =  std::map< std::string , TVectorD >;
using Unc1DMap_t   =  std::map< std::string , std::map< std::string , std::map< std::string , std::map< std::string , std::map< std::string , Unc1DVec_t > > > > >;
using TH1DVec_t    =  std::map< std::string , std::vector< std::tuple< TH1D , TH1D , double > > >;
using TH1DMap_t    =  std::map< std::string , std::map< std::string , std::map< std::string , std::map< std::string , std::map< std::string , TH1DVec_t > > > > >;
using EffVec_t     =  std::map< std::string , std::vector< TEfficiency > >;
using EffMap_t     =  std::map< std::string , std::map< std::string , std::map< std::string , std::map< std::string , std::map< std::string , EffVec_t > > > > >;
using VarMap_t     =  std::map< std::string , std::map< std::string , std::map< std::string , std::map< std::string , std::vector< double > > > > >;
using BinMap_t     =  std::map< std::string , std::vector< double > >;
using BinMapMap_t  =  std::map< std::string , BinMap_t >;
using FileInfo_t   =  std::vector< std::pair< std::string , double > >;
using CorrMap_t    =  std::map< std::string , uint >;


// ------------------ FUNCTION -------------------------------
TnPVec_t getMCWeights        ( const std::unique_ptr<HiEvtTree>& evtTree , const CorrMap_t& corrType );
bool     getMCUncertainties  ( Unc1DVec_t& unc , const EffVec_t& eff );
bool     getMCUncertainties  ( Unc1DMap_t& unc , const EffMap_t& eff );
TnPVec_t getTnPScaleFactors  ( const double& pt, const double& eta , const CorrMap_t& corrType );
bool     getTnPUncertainties ( Unc1DVec_t& unc , const EffVec_t& eff );
bool     getTnPUncertainties ( Unc1DMap_t& unc , const EffMap_t& eff );
void     initEff1D           ( TH1DMap_t& h , const BinMapMap_t& binMap , const CorrMap_t& corrType );
bool     fillEff1D           ( TH1DVec_t& h , const bool& pass , const double& xVar , const TnPVec_t& sfTnP , const TnPVec_t& wMC , const double& evtWeight );
bool     fillEff1D           ( TH1DMap_t& h , const bool& pass , const std::string& type , const VarMap_t& var , const TnPVec_t& sfTnP , const TnPVec_t& wMC , const double& evtWeight , const BinMapMap_t& MU_BIN );
bool     loadEff1D           ( EffMap_t& eff, const TH1DMap_t& h );
void     mergeEff            ( EffMap_t& eff );
void     writeEff            ( TFile& file , const EffMap_t& eff , const Unc1DMap_t& unc , const std::string& mainDirName );
void     saveEff             ( const std::string& outDir , const EffMap_t& eff1D , const Unc1DMap_t& unc );
void     setGlobalWeight     ( TH1DMap_t& h , const double& weight , const std::string& sample , const std::string& col );
bool     checkWeights        ( const TH1& pass , const TH1& total );
const char* clStr            ( const std::string& in );


// ------------------ GLOBAL ------------------------------- 
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
const std::vector< std::string > effType = {"Total", "Acceptance"};
//
// Correction Categories
const CorrMap_t corrType_ = {
  { "NoCorr"           , 1   },
  { "TnP_Nominal"      , 1   },
  { "TnP_Stat_MuID"    , 100 },
  { "TnP_Stat_Trig"    , 2   },
  { "TnP_Stat_Iso"     , 100 },
  { "TnP_Syst_MuID"    , 2   },
  { "TnP_Syst_Trig"    , 2   },
  { "TnP_Syst_Iso"     , 2   },
  { "TnP_Syst_BinMuID" , 1   },
  { "TnP_Syst_BinIso"  , 1   },
  { "MC_Syst_PDF"      , 96  }, // EPPS16 + CT14
  { "MC_Syst_Scale"    , 6   }, // Renormalization + Factorization scale variations [0.5 , 2.0]
  { "MC_Syst_Alpha"    , 2   }  // CT14 Alpha_s variations
};
std::vector< double > absEtaTnP_ = { 0.0 , 1.2 , 2.1 , 2.4 };
std::vector< double > etaTnP_    = { -2.4 , -2.1 , -1.6 , -1.2 , -0.9 , -0.6 , -0.3 , 0.0 , 0.3 , 0.6 , 0.9 , 1.2 , 1.6 , 2.1 , 2.4 };
//
// Input Files for analysis
const std::string path_MC = "root://cms-xrd-global.cern.ch//store/group/phys_heavyions/anstahll/EWQAnalysis2017/pPb2016/8160GeV/MC/Embedded/Official";
const std::map< std::string , std::vector< std::pair< std::string , double > > > inputFileMap_ = {
  {"MC_WToMuNu_Plus_pPb"      , { { Form("%s/%s", path_MC.c_str(), "POWHEG/HiEWQForest_Embedded_Official_POWHEG_CT14_EPPS16_WToMuNu_Plus_pPb_8160GeV_20171003.root")  , POWHEG::XSec.at("WToMuNu_Plus").at("pPb")  } } },
  {"MC_WToMuNu_Minus_pPb"     , { { Form("%s/%s", path_MC.c_str(), "POWHEG/HiEWQForest_Embedded_Official_POWHEG_CT14_EPPS16_WToMuNu_Minus_pPb_8160GeV_20171003.root") , POWHEG::XSec.at("WToMuNu_Minus").at("pPb") } } },
  {"MC_WToMuNu_Plus_Pbp"      , { { Form("%s/%s", path_MC.c_str(), "POWHEG/HiEWQForest_Embedded_Official_POWHEG_CT14_EPPS16_WToMuNu_Plus_Pbp_8160GeV_20171003.root")  , POWHEG::XSec.at("WToMuNu_Plus").at("Pbp")  } } },
  {"MC_WToMuNu_Minus_Pbp"     , { { Form("%s/%s", path_MC.c_str(), "POWHEG/HiEWQForest_Embedded_Official_POWHEG_CT14_EPPS16_WToMuNu_Minus_Pbp_8160GeV_20171003.root") , POWHEG::XSec.at("WToMuNu_Minus").at("Pbp") } } }
};
std::map< std::string , std::vector< std::string > > sampleType_;


void correctEfficiency(const std::string workDirName = "NominalCM", const uint applyHFCorr = 1, const bool applyBosonPTCorr = false, const bool applyMuonPTCorr = false)
{
  //
  // Initialize the Kinematic Bin info
  BinMapMap_t  MU_BIN;
  if ( (workDirName == "Nominal") || (workDirName == "CutAndCount") ) {
    const BinMap_t  TMP = {
      { "Eta" , { -2.4 , -2.2 , -2.0 , -1.8 , -1.6 , -1.4 , -1.2 , -1.0 , -0.8 , -0.6 , -0.4 , -0.2 , 0.0 , 0.2 , 0.4 , 0.6 , 0.8 , 1.0 , 1.2 , 1.4 , 1.6 , 1.8 , 2.0 , 2.2 , 2.4 } }
    };
    for (const auto& v : TMP) { MU_BIN["PA"][v.first] = v.second; MU_BIN["pPb"][v.first] = v.second; MU_BIN["Pbp"][v.first] = v.second; }
  }
  else if ( (workDirName.find("NominalCM")!=std::string::npos) || (workDirName == "CutAndCountCM") ) {
    const BinMap_t  TMP_pPb = {
      { "EtaCM" , { -2.86 , -2.60 , -2.40 , -2.20 , -1.93 , -1.80 , -1.60 , -1.40 , -1.20 , -1.00 , -0.80 , -0.60 , -0.40 , -0.20 , 0.00 , 0.20 , 0.40 , 0.60 , 0.80 , 1.00 , 1.20 , 1.40 , 1.60 , 1.80, 1.93 } }
    };
    const BinMap_t  TMP_Pbp = {
      { "EtaCM" , { -1.93 , -1.80 , -1.60 , -1.40 , -1.20 , -1.00 , -0.80 , -0.60 , -0.40 , -0.20 , 0.00 , 0.20 , 0.40 , 0.60 , 0.80 , 1.00 , 1.20 , 1.40 , 1.60 , 1.80 , 1.93 , 2.20 , 2.40 , 2.60 , 2.86 } }
    };
    for (const auto& v : TMP_pPb) { MU_BIN["PA"][v.first] = v.second; MU_BIN["pPb"][v.first] = v.second; }
    for (const auto& v : TMP_Pbp) { MU_BIN["Pbp"][v.first] = v.second; }
  }
  else if (workDirName == "General") {
    const BinMap_t  TMP = {
      { "Eta" , { -2.4 , -2.2 , -2.0 , -1.8 , -1.6 , -1.4 , -1.2 , -1.0 , -0.8 , -0.6 , -0.4 , -0.2 , 0.0 , 0.2 , 0.4 , 0.6 , 0.8 , 1.0 , 1.2 , 1.4 , 1.6 , 1.8 , 2.0 , 2.2 , 2.4 } },
      { "Pt"  , { 25., 30., 35., 40., 45., 50, 80., 200. } }
    };
    for (const auto& v : TMP) { MU_BIN["PA"][v.first] = v.second; MU_BIN["pPb"][v.first] = v.second; MU_BIN["Pbp"][v.first] = v.second; }
  }
  else {
    std::cout << "[ERROR] WorkDirName " << workDirName << " has not been defined" << std::endl; return;
  }
  //
  // Define the working flags
  bool useEtaCM = false;      if ( (workDirName.find("NominalCM")!=std::string::npos) || (workDirName == "CutAndCountCM") ) { useEtaCM = true; }
  bool isCutAndCount = false; if ( (workDirName == "CutAndCount") || (workDirName == "CutAndCountCM") ) { isCutAndCount = true; }
  bool isMTCut = false;       if (workDirName.find("_MTCut")!=std::string::npos) { isMTCut = true; }
  //
  // Change the working directory
  const std::string CWD = getcwd(NULL, 0);
  const std::string mainDir = Form("%s/Output/", CWD.c_str());
  gSystem->mkdir(mainDir.c_str(), kTRUE);
  gSystem->ChangeDirectory(mainDir.c_str());
  //
  // Create list of samples
  sampleType_["sample"] = {};
  for (const auto& f : inputFileMap_) { std::string tmp = f.first; tmp = tmp.substr(0, tmp.find_last_of("_"));
    if (std::find(sampleType_.at("sample").begin(), sampleType_.at("sample").end(), tmp)==sampleType_.at("sample").end()) { sampleType_.at("sample").push_back(tmp); }
  }
  //
  // ------------------------------------------------------------------------------------------------------------------------
  //
  // Initialize the Correction Type Map
  CorrMap_t corrType;
  //
  corrType["NoCorr"] = 1;
  corrType["TnP_Nominal"] = 1;
  if (applyHFCorr>0) { corrType["HFCorr"] = 1; }
  //
  corrType["TnP_Syst_MuID"   ] = corrType_.at("TnP_Syst_MuID");
  corrType["TnP_Syst_Iso"    ] = corrType_.at("TnP_Syst_Iso");
  corrType["TnP_Syst_BinMuID"] = corrType_.at("TnP_Syst_BinMuID");
  corrType["TnP_Syst_BinIso" ] = corrType_.at("TnP_Syst_BinIso");
  for (uint iEta = 1; iEta < absEtaTnP_.size(); iEta++) {
    const std::string etaLbl = Form("p%.0f_p%.0f", absEtaTnP_[iEta-1]*10., absEtaTnP_[iEta]*10.);
    corrType[Form("TnP_Stat_MuID_%s" , etaLbl.c_str())] = corrType_.at("TnP_Stat_MuID");
    corrType[Form("TnP_Stat_Iso_%s"  , etaLbl.c_str())] = corrType_.at("TnP_Stat_Iso");
  }
  //
  corrType["TnP_Syst_Trig"] = corrType_.at("TnP_Syst_Trig");
  for (uint iEta = 1; iEta < etaTnP_.size(); iEta++) {
    const std::string etaLbl = Form("%s%.0f_%s%.0f", (etaTnP_[iEta-1]<0.?"m":"p") , std::abs(etaTnP_[iEta-1])*10. , (etaTnP_[iEta]<0.?"m":"p") , std::abs(etaTnP_[iEta])*10.);
    corrType[Form("TnP_Stat_Trig_%s" , etaLbl.c_str())] = corrType_.at("TnP_Stat_Trig");
  }
  //
  bool doMCWeight = false;
  for (const auto& cor : corrType_) { if (cor.first.find("MC_")!=std::string::npos) { corrType[cor.first] = cor.second; doMCWeight = true; } }
  //
  // Initialize the HF corrections
  std::unique_ptr<HFweight> corrHF;
  if (applyHFCorr>0) { corrHF = std::unique_ptr<HFweight>(new HFweight("/afs/cern.ch/work/e/echapon/public/DY_pA_2016/HFweight.root")); }
  //
  // Initialize the Rochester Corrector
  std::unique_ptr<RoccoR> RCCorr;
  std::unique_ptr<TRandom3> rnd;
  if (applyMuonPTCorr) {
    const std::string rcPath = "/home/llr/cms/stahl/ElectroWeakAnalysis/EWQAnalysis2017/Utilities/rcdata.2016.v3/";
    RCCorr.reset(new RoccoR(rcPath.c_str()));
    rnd.reset(new TRandom3());
  }
  //
  // ------------------------------------------------------------------------------------------------------------------------
  //
  // Declare the histograms for the efficiency
  TH1DMap_t h1D;   // Stores the total and passing histograms separately
  //
  // Initialize the efficiencies
  initEff1D(h1D , MU_BIN, corrType);
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
    metTree[sample] = std::unique_ptr<HiMETTree>(new HiMETTree());
    if (!metTree.at(sample)->GetTree(fileInfo)) return;
    if (metTree.at(sample)->GetEntries() != nentries.at(sample)) { std::cout << "[ERROR] Inconsistent number of entries between MET (" << 
        metTree.at(sample)->GetEntries() << ") and Muon Tree (" << muonTree.at(sample)->GetEntries() << ") !" << std::endl; return; }
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
    for (Long64_t jentry = 0; jentry < nentries.at(sample); jentry++) { //
      //
      // Get the entry in the trees
      if (muonTree.at(sample)->GetEntry(jentry)<0) { std::cout << "[ERROR] Muon Tree invalid entry!"  << std::endl; return; }
      if ( metTree.at(sample)->GetEntry(jentry)<0) { std::cout << "[ERROR] MET Tree invalid entry!"   << std::endl; return; }
      if ( evtTree.at(sample)->GetEntry(jentry)<0) { std::cout << "[ERROR] HiEVT Tree invalid entry!" << std::endl; return; }
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
      if (muonTree.at(sample)->Event_Run()    != metTree.at(sample)->Event_Run()   ) { std::cout << "[ERROR] MET Run does not agree!"     << std::endl; return; }
      if (muonTree.at(sample)->Event_Number() != metTree.at(sample)->Event_Number()) { std::cout << "[ERROR] MET Event does not agree!"   << std::endl; return; }
      if (muonTree.at(sample)->Event_Run()    != evtTree.at(sample)->run()         ) { std::cout << "[ERROR] HiEVT Run does not agree!"   << std::endl; return; }
      if (muonTree.at(sample)->Event_Number() != evtTree.at(sample)->evt()         ) { std::cout << "[ERROR] HiEVT Event does not agree!" << std::endl; return; }
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
      //const double mcWeight = ( ( lumi ) / muonTree.at(sample)->GetTreeEntries() );
      // Set the global weight only in the first event (i.e. once per sample)
      if (jentry==0) { setGlobalWeight(h1D, mcWeight, sampleType, col); }
      // Define the event weight (set to 1.0 by default)
      double evtWeight = 1.0;
      //
      // Determine the HF Weight
      if (applyHFCorr==1) {
        const double hf = ( (evtTree.at(sample)->hiHF()>=300.) ? 290. : evtTree.at(sample)->hiHF() ); // The histograms are made up to 300.
        if (hf<0.) { std::cout << "[ERROR] The hiHF is negative" << std::endl; return; }
        evtWeight *= corrHF->weight(hf , HFweight::HFside::both , false);
      }
      else if (applyHFCorr==2) {
        const double nTrks = ( (evtTree.at(sample)->hiNtracks()>=300.) ? 290. : evtTree.at(sample)->hiNtracks() ); // The histograms are made up to 300.
        if (nTrks<0.) { std::cout << "[ERROR] The nTracks is negative" << std::endl; return; }
        evtWeight *= corrHF->weight(nTrks, HFweight::HFside::track, false);
      }
      //
      // Get the MC PDF Weights
      TnPVec_t wMC; if (doMCWeight) { wMC = getMCWeights(evtTree.at(sample), corrType); }
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
        // Consider Generated Muons within the Pseudo-Rapidity acceptance of CMS
        if ( std::abs(muonTree.at(sample)->Gen_Muon_Mom()[iGenMu].Eta()) > 2.4 ) continue;
        // Apply Muon Corrections
        double evtW = evtWeight;
        // Determine the Boson PT Weight from Drell-Yan analysis
        if (applyBosonPTCorr) {
          // Determine the gen boson pT
          const auto momIdx = PA::getGenMom(iGenMu, sampleType, muonTree.at(sample));
          if (momIdx>=0) {
            const auto momPdg = std::abs(muonTree.at(sample)->Gen_Particle_PdgId()[momIdx]);
            if (momPdg!=24) { std::cout << "[ERROR] Mother pdg " << momPdg << " can not be used to correct the boson pT" << std::endl; return; }
            double bosonPT = muonTree.at(sample)->Gen_Particle_Mom()[momIdx].Pt();
            if (bosonPT<0.5) { bosonPT = 0.5; }
            evtW *= ( 1.0 / ( ( -0.37 * std::pow(bosonPT, -0.37) ) + 1.19 ) );
          }
          else { std::cout << "[ERROR] Mother of muon was not found!" << std::endl; return; }
        }
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
        genMuonVar[sampleType][col][chg]["Pt" ] = { mu_Gen_Pt  };
        if (useEtaCM) {
          const bool   ispPb = (col == "pPb");
          const double etaCM = PA::EtaLABtoCM( mu_Gen_Eta , ispPb );
          genMuonVar[sampleType][col][chg]["EtaCM"]    = { etaCM };
          genMuonVar[sampleType][col][chg]["Pt_EtaCM"] = { mu_Gen_Pt , etaCM };
        }
        else {
          genMuonVar[sampleType][col][chg]["Eta"]    = { mu_Gen_Eta };
          genMuonVar[sampleType][col][chg]["Pt_Eta"] = { mu_Gen_Pt , mu_Gen_Eta };
        }
        // Initialize the boolean flags
        bool passIdentification = false;
        bool passTrigger        = false;
        bool passIsolation      = false;
        // Initialize the Tag-And-Probe scale factos
        TnPVec_t sfTnP = {};
        //
        // Check that the generated muon is within the analysis kinematic range
        if (isGoodGenMuon) {
          // Find the reconstructed muon matched to gen
          const short iPFMu = muonTree.at(sample)->Gen_Muon_PF_Idx()[iGenMu];
          if (iPFMu >= 0) {
            //
            // PF Muon was matched to generated muon
            //
            // Extract the kinematic information of generated muon
            const double mu_PF_Pt  = muonTree.at(sample)->PF_Muon_Mom()[iPFMu].Pt();
            const double mu_PF_Eta = muonTree.at(sample)->PF_Muon_Mom()[iPFMu].Eta();
            // Determine the Tag-And-Probe scale factors
            sfTnP = getTnPScaleFactors(mu_PF_Pt, mu_PF_Eta, corrType);
            // Check the RECO Idx
            const short iRecoMu = muonTree.at(sample)->PF_Muon_Reco_Idx()[iPFMu];
            if (iRecoMu < 0) { std::cout << "[ERROR] Reco idx is negative" << std::endl; return; }
            // Determine the Muon PT (Rochester) scale factor
            double muSF = 1.0;
            if (applyMuonPTCorr) {
              const double mu_PF_Chg  = muonTree.at(sample)->PF_Muon_Charge()[iPFMu];
              const double mu_PF_Phi  = muonTree.at(sample)->PF_Muon_Mom()[iPFMu].Phi();
              const double mu_PF_NTrk = muonTree.at(sample)->Reco_Muon_InTrk_TrkLayers()[iRecoMu];
              if (mu_PF_NTrk>=0) { muSF = RCCorr->kScaleFromGenMC (mu_PF_Chg, mu_PF_Pt, mu_PF_Eta, mu_PF_Phi, mu_PF_NTrk, mu_Gen_Pt, rnd->Rndm(), 7, 6); }
              //else { std::cout << "[ERROR] Leading PF Muon Tracker Layers variable is negative!" << std::endl; return; }
            }
            // Check if the reconstructed muon pass muon ID and kinematic cuts
            passIdentification = PA::isGoodMuon(iPFMu , muonTree.at(sample), muSF);
            // Check if the reconstructed muon is matched to the trigger
            passTrigger = PA::isTriggerMatched(triggerIndex_, iPFMu, muonTree.at(sample));
            // Check if the reconstructed muon pass isolation cuts
            passIsolation = PA::isIsolatedMuon(iPFMu , muonTree.at(sample), muSF);
            // Check if we are using cut and count method
            if (isCutAndCount) {
              // Extract the MET
              const auto& MET = metTree.at(sample)->PF_MET_NoShift_Mom();
              // Recompute the Transver Mass based on the chosen MET
              const double mu_PF_Phi = muonTree.at(sample)->PF_Muon_Mom()[iPFMu].Phi();
              const double muMTVal = PA::getWTransverseMass(mu_PF_Pt, mu_PF_Phi, MET.Mod(), MET.Phi());
              // Define the MET and MT selection
              const bool passCutAndCountSel = ( (MET.Mod() >= 30.0) && (muMTVal >= 60.0) );
              // Include the Cut and Count selection in the Identification flag
              passIdentification = passIdentification && passCutAndCountSel;
            }
            // Check if we are using MT cut method
            if (isMTCut) {
              // Extract the MET
              const auto& MET = metTree.at(sample)->PF_MET_NoShift_Mom();
              // Recompute the Transver Mass based on the chosen MET
              const double mu_PF_Phi = muonTree.at(sample)->PF_Muon_Mom()[iPFMu].Phi();
              const double muMTVal = PA::getWTransverseMass(mu_PF_Pt, mu_PF_Phi, MET.Mod(), MET.Phi());
              // Define the MT selection
              const bool passMTSel = (muMTVal >= 60.0);
              // Include the MT cut selection in the Identification flag
              passIdentification = passIdentification && passMTSel;
            }
          }
          //
          // Total Efficiency (Based on Generated muons)
          //
          if (!fillEff1D(h1D, (passIdentification && passTrigger && passIsolation), "Total", genMuonVar, sfTnP, wMC, evtW, MU_BIN)) { return; }
        }
        //
        // Total Acceptance (Based on Generated muons)
        //
        if (!fillEff1D(h1D, isGoodGenMuon, "Acceptance", genMuonVar, sfTnP, wMC, evtW, MU_BIN)) { return; }
      }
    }
  }
  //
  // ------------------------------------------------------------------------------------------------------------------------
  //
  // Declare the efficiencies
  EffMap_t  eff1D; // Stores the efficiency
  //
  // Load the efficiencies with the histograms
  if (!loadEff1D(eff1D, h1D)) { return; };
  //
  // ------------------------------------------------------------------------------------------------------------------------
  //
  // Merge Efficiencies
  mergeEff(eff1D);
  //
  // ------------------------------------------------------------------------------------------------------------------------
  //
  // Declare the uncertainties container
  Unc1DMap_t  unc1D; // Stores the efficiency
  //
  // Calculate Uncertainties
  if (!getTnPUncertainties(unc1D, eff1D)) { return; };
  //
  if (doMCWeight) { if ( !getMCUncertainties(unc1D, eff1D)) { return; }; }
  //
  // ------------------------------------------------------------------------------------------------------------------------
  //
  // Define the output directory
  //
  // Store the Efficiencies
  //
  std::string outDir = mainDir + workDirName;
  if    (applyBosonPTCorr) { outDir += "_WithBosonPT"; }
  if    (applyMuonPTCorr ) { outDir += "_WithMuonPT";  }
  if      (applyHFCorr==1) { outDir += "_WithHF/";     }
  else if (applyHFCorr==2) { outDir += "_WithNTrack/"; }
  else { outDir += "/"; }
  gSystem->mkdir(outDir.c_str(), kTRUE);
  saveEff(outDir, eff1D, unc1D);
};


TnPVec_t getMCWeights(const std::unique_ptr<HiEvtTree>& evtTree, const CorrMap_t& corrType)
{
  TnPVec_t wMC;
  for (const auto& cor : corrType) {
    if (cor.first.find("MC_")==std::string::npos) continue;
    //
    wMC[cor.first].clear();
    //
    int    idx  = -1;
    double w_MC =  1.0;
    //
    for (uint i = 0; i < cor.second; i++) {
      if (cor.first=="MC_Syst_PDF") {
        if      (i < 40) { w_MC = evtTree->ttbar_w()[285 + i]; } // EPPS16 VARIATIONS
        else if (i < 96) { w_MC = evtTree->ttbar_w()[ 72 + i]; } // CT14 VARIATIONS
      }
      else if (cor.first=="MC_Syst_Scale") {
        // Renormalization+Factorization scale variations
        if      (i<4) { w_MC = evtTree->ttbar_w()[1 + i]; } // 1, 2, 3, 4
        else if (i<5) { w_MC = evtTree->ttbar_w()[2 + i]; } // 6
        else if (i<6) { w_MC = evtTree->ttbar_w()[3 + i]; } // 8
      }
      else if (cor.first=="MC_Syst_Alpha") {
        // CT14 alpha variations (0.117 - 0.119)
        w_MC = evtTree->ttbar_w()[168 + i];
      }
      else {
        std::cout << "[ERROR] The MC variation " << cor.first << " has not been defined" << std::endl;
      }
      //
      wMC.at(cor.first).push_back( w_MC );
    }
  }
  //
  return wMC;
};


TnPVec_t getTnPScaleFactors(const double& pt, const double& eta, const CorrMap_t& corrType)
{
  TnPVec_t sfTnP;
  for (const auto& cor : corrType) {
    //
    double sf_TnP  = 1.0;
    double sf_MuID = tnp_weight_muid_ppb( pt , eta , 0 );
    double sf_Trig = tnp_weight_trg_ppb (      eta , 0 );
    double sf_Iso  = tnp_weight_iso_ppb ( pt , eta , 0 );
    //
    sfTnP[cor.first].clear();
    for (uint i = 1; i <= cor.second; i++) {
      if ( (cor.first=="NoCorr") || (cor.first=="HFCorr") ) { sf_MuID = 1.0; sf_Trig = 1.0; sf_Iso = 1.0; }
      //
      else if (cor.first=="TnP_Stat_MuID"    ) { sf_MuID = tnp_weight_muid_ppb( pt , eta ,  i  ); }
      else if (cor.first=="TnP_Stat_Trig"    ) { sf_Trig = tnp_weight_trg_ppb (      eta ,  i  ); }
      else if (cor.first=="TnP_Stat_Iso"     ) { sf_Iso  = tnp_weight_iso_ppb ( pt , eta ,  i  ); }
      else if (cor.first=="TnP_Syst_MuID"    ) { sf_MuID = tnp_weight_muid_ppb( pt , eta , -i  ); }
      else if (cor.first=="TnP_Syst_Trig"    ) { sf_Trig = tnp_weight_trg_ppb (      eta , -i  ); }
      else if (cor.first=="TnP_Syst_Iso"     ) { sf_Iso  = tnp_weight_iso_ppb ( pt , eta , -i  ); }
      else if (cor.first=="TnP_Syst_BinMuID" ) { sf_MuID = tnp_weight_muid_ppb( pt , eta , -10 ); }
      else if (cor.first=="TnP_Syst_BinIso"  ) { sf_Iso  = tnp_weight_iso_ppb ( pt , eta , -10 ); }
      //
      else if (cor.first.find("MuID_")!=std::string::npos || cor.first.find("Iso_")!=std::string::npos) {
        for (uint iEta = 1; iEta < absEtaTnP_.size(); iEta++) {
          if (std::abs(eta) >= absEtaTnP_[iEta-1] && std::abs(eta) < absEtaTnP_[iEta]) {
            const std::string etaLbl = Form("p%.0f_p%.0f", absEtaTnP_[iEta-1]*10., absEtaTnP_[iEta]*10.);
            if      (cor.first==Form("TnP_Stat_MuID_%s"     , etaLbl.c_str())) { sf_MuID = tnp_weight_muid_ppb( pt , eta ,  i  ); }
            else if (cor.first==Form("TnP_Stat_Iso_%s"      , etaLbl.c_str())) { sf_Iso  = tnp_weight_iso_ppb ( pt , eta ,  i  ); }
            else if (cor.first==Form("TnP_Syst_MuID_%s"     , etaLbl.c_str())) { sf_MuID = tnp_weight_muid_ppb( pt , eta , -i  ); }
            else if (cor.first==Form("TnP_Syst_Iso_%s"      , etaLbl.c_str())) { sf_Iso  = tnp_weight_iso_ppb ( pt , eta , -i  ); }
            else if (cor.first==Form("TnP_Syst_BinMuID_%s"  , etaLbl.c_str())) { sf_MuID = tnp_weight_muid_ppb( pt , eta , -10 ); }
            else if (cor.first==Form("TnP_Syst_BinIso_%s"   , etaLbl.c_str())) { sf_Iso  = tnp_weight_iso_ppb ( pt , eta , -10 ); }
            break;
          }
        }
      }
      //
      else if (cor.first.find("Trig_")!=std::string::npos) {
        for (uint iEta = 1; iEta < etaTnP_.size(); iEta++) {
          if (eta >= etaTnP_[iEta-1] && eta < etaTnP_[iEta]) {
            const std::string etaLbl = Form("%s%.0f_%s%.0f", (etaTnP_[iEta-1]<0.?"m":"p") , std::abs(etaTnP_[iEta-1])*10. , (etaTnP_[iEta]<0.?"m":"p") , std::abs(etaTnP_[iEta])*10.);
            if      (cor.first==Form("TnP_Stat_Trig_%s" , etaLbl.c_str())) { sf_Trig = tnp_weight_trg_ppb( eta ,  i  ); }
            else if (cor.first==Form("TnP_Syst_Trig_%s" , etaLbl.c_str())) { sf_Trig = tnp_weight_trg_ppb( eta , -i  ); }
            break;
          }
        }
      }
      //
      sf_TnP = ( sf_MuID * sf_Trig * sf_Iso );
      //
      sfTnP.at(cor.first).push_back( sf_TnP );
    }
  }
  //
  return sfTnP;
};


bool getMCUncertainties(Unc1DVec_t& unc, const EffVec_t& eff)
{
  std::string label = "";
  if      (eff.count("TnP_Nominal")>0) { label = "TnP_Nominal"; }
  else if (eff.count("HFCorr"     )>0) { label = "HFCorr";      }
  else if (eff.count("NoCorr"     )>0) { label = "NoCorr";      }
  else { return true; }
  // Compute individual uncertainties
  const TEfficiency& nom = eff.at(label)[0];
  const uint nBin = nom.GetCopyTotalHisto()->GetNbinsX();
  for (const auto& co : eff) {
    if (co.first.find("MC_")==std::string::npos) continue;
    unc[co.first].ResizeTo(nBin);
    for (uint iBin = 1; iBin <= nBin; iBin++) {
      double uncVal = 0.0;
      if (co.first=="MC_Syst_PDF") {
        // Central Value
        const double ctVal = nom.GetEfficiency(iBin);
        // Variations (Use the offical EPPS16 approach)
        double errLo = 0.0 , errHi = 0.0;
        for(uint i = 0; i < (co.second.size()/2); i++) {
          const double miVal = co.second[(2*i)+0].GetEfficiency(iBin);
          const double plVal = co.second[(2*i)+1].GetEfficiency(iBin);
          errLo += std::pow( std::min( std::min( (plVal - ctVal) , (miVal - ctVal) ) , 0.0 ) , 2.0 );
          errHi += std::pow( std::max( std::max( (plVal - ctVal) , (miVal - ctVal) ) , 0.0 ) , 2.0 );
        }
        // Convert from 90% CL to 68% CL
        const double convFactor = TMath::ErfcInverse((1.-0.68))/TMath::ErfcInverse((1.-0.90));
        errLo = ( convFactor * std::sqrt( errLo ) );
        errHi = ( convFactor * std::sqrt( errHi ) );
        //
        uncVal = std::max( errLo , errHi );
      }
      else if (co.first=="MC_Syst_Scale") {
        // Central Value
        const double ctVal = nom.GetEfficiency(iBin);
        // Variations (Use the envelope approach)
        double errLo = 0.0 , errHi = 0.0;
        for(uint i = 0; i < co.second.size(); i++) {
          const double vrVal = co.second[i].GetEfficiency(iBin);
          errLo = std::min( std::min( (vrVal - ctVal) , errLo ) , 0.0 );
          errHi = std::max( std::max( (vrVal - ctVal) , errHi ) , 0.0 );
        }
        // We consider the scale uncertainty to be 68% CL, so no need to rescale
        errLo = std::abs( errLo );
        errHi = std::abs( errHi );
        //
        uncVal = std::max( errLo , errHi );
      }
      else if (co.first=="MC_Syst_Alpha") {
        if (co.second.size()>2) { std::cout << "[ERROR] Number of variations (" << co.second.size() << ") for MC_Syst_Alpha is larger than 2" << std::endl; return false; }
        // Variations (Use the PDF4LHC15 approach)
        const double miVal = co.second[0].GetEfficiency(iBin);
        const double plVal = co.second[1].GetEfficiency(iBin);
        double err = ( ( plVal - miVal ) / 2.0 );
        // Convert from deltaAlpha_s = 0.0010 to 0.0015 (68% CL)
        const double convFactor = (0.0015/0.0010);
        err = ( convFactor * std::abs( err ) );
        //
        uncVal = err;
      }
      else {
        std::cout << "[ERROR] MC systematic variation " << co.first << " has not been defined!" << std::endl; return false;
      }
      unc.at(co.first)[iBin-1] = uncVal;
    }
  }
  // Compute total systematic uncertainty
  unc["MC_Syst"].ResizeTo(nBin);
  for (uint iBin = 0; iBin < nBin; iBin++) {
    double uncVal = 0.0;
    for (const auto& co : eff) {
      if (co.first.find("MC_Syst")!=std::string::npos) {
        uncVal += ( unc.at(co.first)[iBin] * unc.at(co.first)[iBin] );
      }
    }
    unc.at("MC_Syst")[iBin] = std::sqrt( uncVal );
  }
  //
  return true;
};


bool getMCUncertainties(Unc1DMap_t& unc, const EffMap_t& eff)
{
  for (const auto& v : eff) {
    for (const auto& s : v.second) {
      for (const auto& c : s.second) {
        for (const auto& ch : c.second) {
          for (const auto& t : ch.second) {
            if (!getMCUncertainties(unc[v.first][s.first][c.first][ch.first][t.first], t.second)) { return false; };
          }
        }
      }
    }
  }
  return true;
};


bool getTnPUncertainties(Unc1DVec_t& unc, const EffVec_t& eff)
{
  if (eff.count("TnP_Nominal") == 0) { return true; }
  // Compute individual uncertainties
  const TEfficiency& nom = eff.at("TnP_Nominal")[0];
  const uint nBin = nom.GetCopyTotalHisto()->GetNbinsX();
  for (const auto& co : eff) {
    if (co.first=="TnP_Nominal" || co.first.find("TnP_")==std::string::npos) continue;
    if (eff.count(co.first+"_")>0) continue;
    unc[co.first].ResizeTo(nBin);
    for (uint iBin = 1; iBin <= nBin; iBin++) {
      double uncVal = 0.0;
      if (co.second.size() == 100) {
        double sum = 0.0;
        for(uint i = 0; i < co.second.size(); i++) {
          const double diff = ( co.second[i].GetEfficiency(iBin) - nom.GetEfficiency(iBin) );
          sum += ( diff * diff );
        }
        uncVal = std::sqrt( sum / 100.0 );
      }
      else if (co.second.size() == 2) {
        uncVal = std::max(std::abs(co.second[0].GetEfficiency(iBin) - nom.GetEfficiency(iBin)) , std::abs(co.second[1].GetEfficiency(iBin) - nom.GetEfficiency(iBin)));
      }
      else if (co.second.size() == 1) {
        uncVal = std::abs(co.second[0].GetEfficiency(iBin) - nom.GetEfficiency(iBin));
      }
      else {
        std::cout << "[ERROR] Number of variations " << co.second.size() << " is wrong!" << std::endl; return false;
      }
      unc.at(co.first)[iBin-1] = uncVal;
    }
  }
  // Compute total statistical uncertainty
  unc["TnP_Stat"].ResizeTo(nBin);
  for (uint iBin = 0; iBin < nBin; iBin++) {
    double uncVal = 0.0;
    for (const auto& co : eff) {
      if (co.first.find("TnP_Stat")!=std::string::npos) {
        uncVal += std::pow( unc.at(co.first)[iBin] , 2.0 );
      }
    }
    unc.at("TnP_Stat")[iBin] = std::sqrt( uncVal );
  }
  // Compute total systematic uncertainty
  unc["TnP_Syst"].ResizeTo(nBin);
  for (uint iBin = 0; iBin < nBin; iBin++) {
    double uncVal = 0.0;
    for (const auto& co : eff) {
      if (co.first.find("TnP_Syst")!=std::string::npos) {
        uncVal += std::pow( unc.at(co.first)[iBin] , 2.0 );
      }
    }
    uncVal += std::pow( ( 0.0034 * nom.GetEfficiency(iBin+1) ) , 2.0 ); // Impact of PileUp and Event Activity on Isolation, defined as relative uncertainty
    uncVal += std::pow( ( 0.0060 * nom.GetEfficiency(iBin+1) ) , 2.0 ); // STA Efficiency mismodelling, defined as relative uncertainty
    unc.at("TnP_Syst")[iBin] = std::sqrt( uncVal );
  }
  // Compute total TnP uncertainty
  unc["TnP_Tot"].ResizeTo(nBin);
  for (uint iBin = 0; iBin < nBin; iBin++) {
    double uncVal = 0.0;
    uncVal += std::pow( unc.at("TnP_Syst")[iBin] , 2.0 );
    uncVal += std::pow( unc.at("TnP_Stat")[iBin] , 2.0 );
    unc.at("TnP_Tot")[iBin] = std::sqrt( uncVal );
  }
  //
  return true;
};


bool getTnPUncertainties(Unc1DMap_t& unc, const EffMap_t& eff)
{
  for (const auto& v : eff) {
    for (const auto& s : v.second) {
      for (const auto& c : s.second) {
        for (const auto& ch : c.second) {
          for (const auto& t : ch.second) {
            if (!getTnPUncertainties(unc[v.first][s.first][c.first][ch.first][t.first], t.second)) { return false; };
          }
        }
      }
    }
  }
  return true;
};


void initEff1D(TH1DMap_t& h, const BinMapMap_t& binMap, const CorrMap_t& corrType)
{
  for (const auto& col : COLL_) {
    for (const auto& bins : binMap.at(col)) {
      Double_t bin[bins.second.size()];
      for (uint i = 0; i < bins.second.size(); i++) { bin[i] = bins.second[i]; }
      for (const auto& sample : sampleType_.at("sample")) {
        for (const auto& chg : CHG_) {
          for (const auto& type : effType) {
            for (const auto& cor : corrType) {
              if (type=="Acceptance" && (cor.first.find("TnP_")!=std::string::npos)) continue;
              for (uint i = 0; i < cor.second; i++) {
                h[bins.first][sample][col][chg][type][cor.first].push_back( std::make_tuple(TH1D() , TH1D() , 1.0) );
                std::get<0>(h.at(bins.first).at(sample).at(col).at(chg).at(type).at(cor.first)[i]) = TH1D("Passed", "Passed", (bins.second.size()-1), bin);
                std::get<1>(h.at(bins.first).at(sample).at(col).at(chg).at(type).at(cor.first)[i]) = TH1D("Total" , "Total" , (bins.second.size()-1), bin);
                std::get<2>(h.at(bins.first).at(sample).at(col).at(chg).at(type).at(cor.first)[i]) = 1.0;
                //
                const std::string name = "h1D_" + bins.first +"_"+ sample +"_"+ col +"_"+ chg +"_"+ type +"_"+ cor.first + ((cor.second>1) ? Form("_%d",(i+1)) : "");
                std::get<0>(h.at(bins.first).at(sample).at(col).at(chg).at(type).at(cor.first)[i]).SetName((name+"_Passed").c_str());
                std::get<1>(h.at(bins.first).at(sample).at(col).at(chg).at(type).at(cor.first)[i]).SetName((name+"_Total").c_str());
                std::get<0>(h.at(bins.first).at(sample).at(col).at(chg).at(type).at(cor.first)[i]).Sumw2();
                std::get<1>(h.at(bins.first).at(sample).at(col).at(chg).at(type).at(cor.first)[i]).Sumw2();
              }
            }
          }
        }
      }
    }
  }
};


bool fillEff1D(TH1DVec_t& h, const bool& pass, const double& xVar, const TnPVec_t& sfTnP, const TnPVec_t& wMC, const double& evtWeight, const std::string& type)
{
  for (auto& cor : h) {
    for (uint i = 0; i < cor.second.size(); i++) {
      //
      bool found = false;
      //
      double sf = 1.0;
      if (sfTnP.count(cor.first)>0 && sfTnP.at(cor.first).size()>i) { sf = sfTnP.at(cor.first)[i]; found = true; }
      else if (sfTnP.size()==0 || sfTnP.count(cor.first)==0) { sf = 1.0; }
      else { std::cout << "[ERROR] Correction " << cor.first << " has invalid number of entries: " << cor.second.size() << "  " << sfTnP.at(cor.first).size() << " !" << std::endl; return false; }
      if ((cor.first.find("TnP_")!=std::string::npos) && pass && sfTnP.size()==0) { std::cout << "[ERROR] TnP scale factor vector is empty!" << std::endl; return false; }
      //
      double w_MC = 1.0;
      if (wMC.count(cor.first)>0 && wMC.at(cor.first).size()>i) { w_MC = wMC.at(cor.first)[i]; found = true; }
      else if (wMC.size()==0 || wMC.count(cor.first)==0) { w_MC = 1.0; }
      else { std::cout << "[ERROR] MC PDF Weight " << cor.first << " has invalid number of entries: " << cor.second.size() << "  " << wMC.at(cor.first).size() << " !" << std::endl; return false; }
      //
      if (found==false && ( (cor.first.find("MC_")!=std::string::npos && wMC.size()>0) || (sfTnP.size()>0) )) {
        std::cout << "[ERROR] Correction " << cor.first << " was not found!" << std::endl; return false;
      }
      //
      // Fill the Pass histogram
      if      (cor.first=="NoCorr") { if (pass) { std::get<0>(cor.second[i]).Fill(xVar , 1.0                 ); } }
      else if (cor.first=="HFCorr") { if (pass) { std::get<0>(cor.second[i]).Fill(xVar , evtWeight           ); } }
      else if (type=="Acceptance" ) { if (pass) { std::get<0>(cor.second[i]).Fill(xVar , (evtWeight*w_MC)    ); } }
      else                          { if (pass) { std::get<0>(cor.second[i]).Fill(xVar , (evtWeight*sf*w_MC) ); } }
      // Fill the total histogram
      if      (cor.first=="NoCorr") { std::get<1>(cor.second[i]).Fill(xVar , 1.0              ); }
      else if (cor.first=="HFCorr") { std::get<1>(cor.second[i]).Fill(xVar , evtWeight        ); }
      else                          { std::get<1>(cor.second[i]).Fill(xVar , (evtWeight*w_MC) ); }
    }
  }
  return true;
};


bool fillEff1D(TH1DMap_t& h, const bool& pass, const std::string& type, const VarMap_t& var, const TnPVec_t& sfTnP, const TnPVec_t& wMC, const double& evtWeight, const BinMapMap_t& MU_BIN)
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
          if (c.first=="PA" && col=="Pbp") {
            if (v.first=="Eta") { xVar = -xVar; } // Invert in the LAB frame
            if (v.first=="EtaCM") {
              xVar = PA::EtaCMtoLAB(xVar, false); // Switch from CM -> LAB, in Pbp system
              xVar = -xVar;                       // Invert in the LAB frame, so from Pbp -> pPb
              xVar = PA::EtaLABtoCM(xVar, true);  // Switch from LAB -> CM, in pPb system
            }
          }
          if (xVar >= MU_BIN.at(c.first).at(v.first)[0] && xVar <= MU_BIN.at(c.first).at(v.first)[MU_BIN.at(c.first).at(v.first).size()-1]) { // Don't include values outside of range
            // Fill histograms
            if (!fillEff1D(ch.second.at(type), pass, xVar, sfTnP, wMC, evtWeight, type)) { return false; }
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
            for (const auto& co : t.second) {
              for (uint i = 0; i < co.second.size(); i++) {
                const TH1D&   hPassed = std::get<0>(co.second[i]);
                const TH1D&   hTotal  = std::get<1>(co.second[i]);
                const double& weight  = std::get<2>(co.second[i]);
                const std::string name = "eff1D_" + v.first +"_"+ s.first +"_"+ c.first +"_"+ ch.first +"_"+ t.first +"_"+ co.first + ((co.second.size()>1) ? Form("_%d",(i+1)) : "");
                if ( TEfficiency::CheckConsistency(hPassed, hTotal,"w") ) {
                  eff[v.first][s.first][c.first][ch.first][t.first][co.first].push_back( TEfficiency(hPassed , hTotal) );
                  eff.at(v.first).at(s.first).at(c.first).at(ch.first).at(t.first).at(co.first)[i].SetName(name.c_str());
                  // Set Global Weight
                  eff.at(v.first).at(s.first).at(c.first).at(ch.first).at(t.first).at(co.first)[i].SetWeight(weight);
                  if ( checkWeights(hPassed, hTotal) ) {
                    eff.at(v.first).at(s.first).at(c.first).at(ch.first).at(t.first).at(co.first)[i].SetStatisticOption(TEfficiency::kBJeffrey);
                  }
                }
                else {
                  std::cout << "[ERROR] 1D Histograms for " << name << " are not consistent!" << std::endl;
                  for(Int_t i = 0; i < (hPassed.GetNbinsX()+2); ++i) {
                    if(hPassed.GetBinContent(i) > hTotal.GetBinContent(i)) {
                      std::cout << "[ERROR] Bin " << i << " has Pass: " << hPassed.GetBinContent(i) << " Total: " << hTotal.GetBinContent(i) << std::endl;
                    }
                  }
                  eff[v.first][s.first][c.first][ch.first][t.first][co.first].push_back( TEfficiency() );
                  eff.at(v.first).at(s.first).at(c.first).at(ch.first).at(t.first).at(co.first)[i].SetName(name.c_str());
                  eff.at(v.first).at(s.first).at(c.first).at(ch.first).at(t.first).at(co.first)[i].SetTotalHistogram(hTotal, "f");
                  eff.at(v.first).at(s.first).at(c.first).at(ch.first).at(t.first).at(co.first)[i].SetPassedHistogram(hPassed, "f");
                  // Set Global Weight
                  eff.at(v.first).at(s.first).at(c.first).at(ch.first).at(t.first).at(co.first)[i].SetWeight(weight);
                  if ( checkWeights(hPassed, hTotal) ) {
                    eff.at(v.first).at(s.first).at(c.first).at(ch.first).at(t.first).at(co.first)[i].SetStatisticOption(TEfficiency::kBJeffrey);
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return true;
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
              for (auto& co : t.second) {
                for (uint i = 0; i < co.second.size(); i++) {
                  // Just add the pPb, no need to combine
                  const TEfficiency& eff_pPb = s.second.at("pPb").at(ch.first).at(t.first).at(co.first)[i];
                  const TEfficiency& eff_Pbp = s.second.at("Pbp").at(ch.first).at(t.first).at(co.first)[i]; // Only used for the weight
                  // Passed Histogram
                  TH1D hPassed = *((TH1D*)eff_pPb.GetPassedHistogram());
                  hPassed.Add(eff_pPb.GetPassedHistogram(), co.second[i].GetPassedHistogram(), eff_pPb.GetWeight(), eff_Pbp.GetWeight());
                  co.second[i].SetPassedHistogram(hPassed, "f");
                  // Total Histogram
                  TH1D hTotal = *((TH1D*)eff_pPb.GetTotalHistogram());
                  hTotal.Add(eff_pPb.GetTotalHistogram(), co.second[i].GetTotalHistogram(), eff_pPb.GetWeight(), eff_Pbp.GetWeight());
                  co.second[i].SetTotalHistogram(hTotal, "f");
                  // Set the statistics in case of weighted histograms
                  if ( checkWeights(hPassed, hTotal) ) { co.second[i].SetStatisticOption(TEfficiency::kBJeffrey); }
                }
              }
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
              for (auto& co : t.second) {
                for (uint i = 0; i < co.second.size(); i++) {
                  // Use the Plus W MC sample for the plus muon efficiency and the Minus W MC sample for the minus muon efficiency
                  const std::string name = "eff1D_" + v.first +"_"+ sample +"_"+ c.first +"_"+ ch.first +"_"+ t.first +"_"+ co.first + ((co.second.size()>1) ? Form("_%d",i) : "");
                  eff.at(v.first)[sample][c.first][ch.first][t.first][co.first].push_back(eff.at(v.first).at(sample+"_"+ch.first).at(c.first).at(ch.first).at(t.first).at(co.first)[i]);
                  eff.at(v.first).at(sample).at(c.first).at(ch.first).at(t.first).at(co.first)[i].SetName(name.c_str());
                }
              }
            }
          }
        }
        eff.at(v.first).erase(sample+"_Minus");
        eff.at(v.first).erase(sample+"_Plus");
      }
    }
  }
};


void writeEff(TFile& file, const EffMap_t& eff, const Unc1DMap_t& unc, const std::string& mainDirName)
{
  if (eff.size()>0) {
    TDirectory*  mainDir = file.mkdir(mainDirName.c_str());
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
              TDirectory* typeDir = colDir->GetDirectory(t.first.c_str());
              if (typeDir==NULL) { typeDir = colDir->mkdir(t.first.c_str()); }
              typeDir->cd();
              const std::map<std::string , TVectorD>& u = unc.at(v.first).at(s.first).at(c.first).at(ch.first).at(t.first);
              const std::string name = "unc1D_" + v.first +"_"+ s.first +"_"+ c.first +"_"+ ch.first +"_"+ t.first +"_";
              for (auto& co : t.second) {
                TDirectory* corDir = typeDir->GetDirectory(co.first.c_str());
                if (corDir==NULL) { corDir = typeDir->mkdir(co.first.c_str()); }
                corDir->cd();
                for (uint i = 0; i < co.second.size(); i++) {
                  co.second[i].Write(clStr(co.second[i].GetName()));
                }
                if (u.count(co.first)>0) {
                  if (co.first.find("TnP_S"  )!=std::string::npos) { u.at(co.first).Write((name + co.first).c_str()); }
                  if (co.first.find("MC_Syst")!=std::string::npos) { u.at(co.first).Write((name + co.first).c_str()); }
                }
                typeDir->cd();
              }
              if (u.count("TnP_Tot") > 0) { 
                u.at("TnP_Stat").Write((name + "TnP_Stat").c_str());
                u.at("TnP_Syst").Write((name + "TnP_Syst").c_str());
                u.at("TnP_Tot").Write((name + "TnP_Tot").c_str());
              }
              if (u.count("MC_Syst") > 0) { u.at("MC_Syst").Write((name + "MC_Syst").c_str()); }
              colDir->cd();
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


void saveEff(const std::string& outDir, const EffMap_t& eff1D, const Unc1DMap_t& unc1D)
{
  // Step 0: Create the output file
  const std::string fileName = "efficiencyTnP.root";
  TFile file((outDir + fileName).c_str(), "RECREATE");
  file.cd();
  // Step 1: Store all the 1D Efficiency objects
  writeEff(file, eff1D, unc1D, "TnPEfficiency1D");
  // Step 3: Write and Close the file
  file.Write();
  file.Close("R");
  std::cout << "[INFO] Information store in file " << outDir << fileName << std::endl;
};


void setGlobalWeight(TH1DMap_t& h, const double& weight, const std::string& sample, const std::string& col)
{
  for (auto& v : h) {
    if (v.second.count(sample)>0 && v.second.at(sample).count(col)>0) {
      for (auto& ch : v.second.at(sample).at(col)) {
        for (auto& t : ch.second) {
          for (auto& co : t.second) {
            for (auto& p : co.second) {
              std::get<2>(p) = weight;
            }
          }
        }
      }
    }
  }
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


const char* clStr(const std::string& in)
{
  std::string out = in;
  while (out.find("_copy")!=std::string::npos) { out.erase(out.find("_copy"), 5); }
  return out.c_str();
};

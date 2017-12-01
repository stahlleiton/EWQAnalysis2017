#if !defined(__CINT__) || defined(__MAKECINT__)
// Auxiliary Headers
#include "../../Utilities/HiMETTree.h"
#include "../../Utilities/HiMuonTree.h"
#include "../../Utilities/EVENTUTILS.h"
// ROOT headers
#include "TROOT.h"
#include "TSystem.h"
#include "TVectorD.h"
// c++ headers
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <string>

#endif


// ------------------ TYPE -------------------------------
using FilterMap_t  =  std::map< std::string , std::map< std::string , std::map< std::string , std::map< std::string , uint > > > >;


// ------------------ FUNCTION -------------------------------


// ------------------ GLOBAL ------------------------------- 
//
// Trigger Info
const ushort triggerIndex_ = PA::HLT_PAL3Mu12;
//
// Collision System
const std::vector< std::string > COLL_ = { "pPb" , "Pbp" , "PA" };
//
const std::vector< std::string > EVT_ = { "Region_DY" , "Region_QCD" , "Region_W" };
//
const std::vector< std::string > FILTERS_ = {
  "ALLEVENTS",
  "ALLFILTERS",
  "collisionEventSelectionPA",
  "collisionEventSelectionPA_rejectPU",
  "BadChargedCandidateFilter",
  "BadPFMuonFilter",
  "EcalDeadCellTriggerPrimitiveFilter",
  "HBHENoiseFilter",
  "HBHENoiseIsoFilter",
  "duplicateMuons",
  "eeBadScFilter",
  "globalTightHalo2016Filter",
  "goodVertices",
  "badMuons"
};
const std::vector< std::string > FTYPE_ = { "ONLY" , "WithPA" };
//
// Input Files for analysis
const std::string path_DATA = "root://cms-xrd-global.cern.ch//store/group/phys_heavyions/anstahll/EWQAnalysis2017/pPb2016/8160GeV/Data/";
const std::string path_MC   = "root://cms-xrd-global.cern.ch//store/group/phys_heavyions/anstahll/EWQAnalysis2017/pPb2016/8160GeV/MC/Embedded/Official";
const std::map< std::string , std::vector< std::string > > inputFileMap_ = {
  {"DATA" , { Form("%s/%s", path_DATA.c_str(), "HiEWQForest_PASingleMuon_pPb_Pbp_8160GeV_20171003.root") } },
  {"MC_TTall_pPb"    , { Form("%s/%s", path_MC.c_str(), "POWHEG/HiEWQForest_Embedded_Official_POWHEG_CT14_EPPS16_TTall_pPb_8160GeV_20171003.root") } },
  {"MC_TTall_Pbp"    , { Form("%s/%s", path_MC.c_str(), "POWHEG/HiEWQForest_Embedded_Official_POWHEG_CT14_EPPS16_TTall_Pbp_8160GeV_20171003.root") } },
  {"MC_DYToMuMu_M_10_30_pPb" , { Form("%s/%s", path_MC.c_str(), "POWHEG/HiEWQForest_Embedded_Official_POWHEG_CT14_EPPS16_DYtoMuMu_M_10_30_pPb_8160GeV_20171003.root") } },
  {"MC_DYToMuMu_M_30_Inf_pPb" , { Form("%s/%s", path_MC.c_str(), "POWHEG/HiEWQForest_Embedded_Official_POWHEG_CT14_EPPS16_DYtoMuMu_M_30_pPb_8160GeV_20171003.root") } },
  {"MC_DYToMuMu_M_10_30_Pbp" , { Form("%s/%s", path_MC.c_str(), "POWHEG/HiEWQForest_Embedded_Official_POWHEG_CT14_EPPS16_DYtoMuMu_M_10_30_Pbp_8160GeV_20171003.root") } },
  {"MC_DYToMuMu_M_30_Inf_Pbp" , { Form("%s/%s", path_MC.c_str(), "POWHEG/HiEWQForest_Embedded_Official_POWHEG_CT14_EPPS16_DYtoMuMu_M_30_Pbp_8160GeV_20171003.root") } },
  {"MC_WToTauNu_Plus_pPb" , { Form("%s/%s", path_MC.c_str(), "POWHEG/HiEWQForest_Embedded_Official_POWHEG_CT14_EPPS16_TAUOLA_WToTauNu_Plus_pPb_8160GeV_20171003.root") } },
  {"MC_WToTauNu_Minus_pPb" , { Form("%s/%s", path_MC.c_str(), "POWHEG/HiEWQForest_Embedded_Official_POWHEG_CT14_EPPS16_TAUOLA_WToTauNu_Minus_pPb_8160GeV_20171003.root") } },
  {"MC_WToTauNu_Plus_Pbp" , { Form("%s/%s", path_MC.c_str(), "POWHEG/HiEWQForest_Embedded_Official_POWHEG_CT14_EPPS16_TAUOLA_WToTauNu_Plus_Pbp_8160GeV_20171003.root") } },
  {"MC_WToTauNu_Minus_Pbp" , { Form("%s/%s", path_MC.c_str(), "POWHEG/HiEWQForest_Embedded_Official_POWHEG_CT14_EPPS16_TAUOLA_WToTauNu_Minus_Pbp_8160GeV_20171003.root") } },
  {"MC_WToMuNu_Plus_pPb"      , { Form("%s/%s", path_MC.c_str(), "POWHEG/HiEWQForest_Embedded_Official_POWHEG_CT14_EPPS16_WToMuNu_Plus_pPb_8160GeV_20171003.root")  } },
  {"MC_WToMuNu_Minus_pPb"     , { Form("%s/%s", path_MC.c_str(), "POWHEG/HiEWQForest_Embedded_Official_POWHEG_CT14_EPPS16_WToMuNu_Minus_pPb_8160GeV_20171003.root") } },
  {"MC_WToMuNu_Plus_Pbp"      , { Form("%s/%s", path_MC.c_str(), "POWHEG/HiEWQForest_Embedded_Official_POWHEG_CT14_EPPS16_WToMuNu_Plus_Pbp_8160GeV_20171003.root")  } },
  {"MC_WToMuNu_Minus_Pbp"     , { Form("%s/%s", path_MC.c_str(), "POWHEG/HiEWQForest_Embedded_Official_POWHEG_CT14_EPPS16_WToMuNu_Minus_Pbp_8160GeV_20171003.root") } }
};
std::map< std::string , std::vector< std::string > > sampleType_;


void checkMETFIlters(const bool print = true)
{
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
  // Initiliaze the Filter Map
  FilterMap_t Filter;
  for (const auto& c : COLL_) {
    for (const auto& e : EVT_) {
      for (const auto& f : FILTERS_) {
        for (const auto& t : FTYPE_) {
          for (const auto& s : sampleType_.at("sample")) {
            Filter[c][e][Form("%s_%s", f.c_str(), t.c_str())][s] = 0;
          }
        }
      }
    }
  };
  //
  // ------------------------------------------------------------------------------------------------------------------------
  //
  // Extract all the samples
  std::map< std::string , Long64_t > nentries;
  std::map< std::string , std::unique_ptr< HiMuonTree > > muonTree;
  std::map< std::string , std::unique_ptr< HiMETTree  > > metTree;
  for (const auto & inputFile : inputFileMap_) {
    const std::string sample = inputFile.first;
    const std::vector< std::string > fileInfo = inputFile.second;
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
    int treeIdx = -1;
    const bool isData = (sample=="DATA");
    const bool isMC   = (isData==false);
    std::cout << "[INFO] Starting to process " << nentries.at(sample) << " nentries" << std::endl;
    //
    for (Long64_t jentry = 0; jentry < nentries.at(sample); jentry++) { //
      //
      // Get the entry in the trees
      if (muonTree.at(sample)->GetEntry(jentry)<0    ) { std::cout << "[ERROR] Muon Tree invalid entry!"     << std::endl; return; }
      if ( metTree.at(sample)->GetEntry(jentry)<0    ) { std::cout << "[ERROR] MET Tree invalid entry!"      << std::endl; return; }
      //
      if (muonTree.at(sample)->Chain()->GetTreeNumber()!=treeIdx) {
        treeIdx = muonTree.at(sample)->Chain()->GetTreeNumber();
        std::cout << "[INFO] Processing Root File: " << inputFile.second[treeIdx] << std::endl;
      }
      //
      loadBar(jentry, nentries.at(sample));
      // 
      // Check that the different tree agrees well
      if (muonTree.at(sample)->Event_Run()    != metTree.at(sample)->Event_Run()   ) { std::cout << "[ERROR] MET Run does not agree!"     << std::endl; return; }
      if (muonTree.at(sample)->Event_Number() != metTree.at(sample)->Event_Number()) { std::cout << "[ERROR] MET Event does not agree!"   << std::endl; return; }
      //
      // Determine the collision system of the sample
      std::string evtCol = "";
      if (isData) {
        if (muonTree.at(sample)->Event_Run() >= 285410 && muonTree.at(sample)->Event_Run() <= 285951) evtCol = "Pbp"; // for Pbp
        if (muonTree.at(sample)->Event_Run() >= 285952 && muonTree.at(sample)->Event_Run() <= 286504) evtCol = "pPb"; // for pPb
      }
      else {
        if (sample.find("Pbp")!=std::string::npos) evtCol = "Pbp"; // for Pbp
        if (sample.find("pPb")!=std::string::npos) evtCol = "pPb"; // for pPb
      }
      if (evtCol=="") { std::cout << "[ERROR] Could not determine the collision system in the sample" << std::endl; return; }
      //
      // Determine the type of sample : i.e. MC_DYToMuMu
      std::string sampleType = sample;
      if (isMC) { sampleType = sampleType.substr(0, (sampleType.find(evtCol)-1)); }
      //
      // Event Based Information
      //
      //
      // Check Trigger Fired
      if (muonTree.at(sample)->Event_Trig_Fired()[triggerIndex_]==false) continue; // Trigger Event Fired
      //
      // Muon Based Information
      //
      // Find the Leading Pt Muon
      float maxPt = -99.; int leadMuPFIdx = -1;
      for (uint imu = 0; imu < muonTree.at(sample)->PF_Muon_Mom().size(); imu++) {
        if (muonTree.at(sample)->PF_Muon_Reco_Idx()[imu]==-1) { std::cout << "[ERROR] Reco idx is -1" << std::endl; return; }
        // Only consider Muons Matched to GEN in MC
        if (isMC && muonTree.at(sample)->PF_Muon_Gen_Idx()[imu]<0) continue;
        if (isMC && PA::checkGenMuon(muonTree.at(sample)->PF_Muon_Gen_Idx()[imu], sample, muonTree.at(sample)) == false) continue;
        // Only consider Tight Muons within the Pseudo-Rapidity acceptance of CMS
        if (PA::isTightMuon(imu, muonTree.at(sample))==false) continue;
        // Determine the highest pT muon
        if (maxPt < muonTree.at(sample)->PF_Muon_Mom()[imu].Pt()) { maxPt = muonTree.at(sample)->PF_Muon_Mom()[imu].Pt(); leadMuPFIdx = imu; }
      }
      if (leadMuPFIdx<=-1) continue;  // Only consider events with at least one good muon
      //
      // Apply PT Cut
      TLorentzVector muP4 = muonTree.at(sample)->PF_Muon_Mom()[leadMuPFIdx];
      if (muP4.Pt() < 25.) continue;  // Consider only leading muons with pT larger or equal than 25 GeV
      //
      // Apply Trigger Matching
      if (PA::isTriggerMatched(triggerIndex_, leadMuPFIdx, muonTree.at(sample))==false) continue;  // Only consider Muons Matched to Trigger
      //
      // Isolation and Charge of Leading Muon
      const float leadMuIso = muonTree.at(sample)->PF_Muon_IsoPFR03NoPUCorr()[leadMuPFIdx];
      const int   leadMuChg = int(muonTree.at(sample)->PF_Muon_Charge()[leadMuPFIdx]);
      //
      //
      // Classify the events
      //
      std::string eventType = "Other";
      //
      // Event Type: DrellYan->MuMu
      if (eventType=="Other") {
        if (PA::passDrellYanVeto(muonTree.at(sample)) == false) { eventType = "Region_DY"; }  // Found a Drell-Yan candidate
      }
      //
      // Event Type: QCD -> Muon + X
      //
      if ( eventType=="Other" ) {
        if (
            ( leadMuIso >= 0.4 ) && // Select Non-Isolated muons
            ( muP4.Pt() >  25. )    // Consider only leading muons with pT larger or equal than 25 GeV
            )
          {
            eventType = "Region_QCD";
          }
      }
      //
      // Event Type: W -> Muon + Neutrino
      //
      if ( eventType=="Other" ) {
        if (
            ( leadMuIso < 0.15 ) &&  // Select Isolated muons
            ( muP4.Pt() > 25.  )     // Consider only leading muons with pT larger or equal than 25 GeV
            )
          {
            eventType = "Region_W";
          }
      }
      if (eventType=="Other") continue;
      //
      // Apply Event Filters
      //
      for (const auto& c : COLL_) {
        if (c!="PA" && c!=evtCol) continue;
        for (const auto& e : EVT_) {
          if (e!=eventType) continue;
          for (const auto& t : FTYPE_) {
            for (const auto& s : sampleType_.at("sample")) {
              if (s!=sampleType) continue;
              //
              if (
                  (t=="ONLY") ||
                  (t=="WithPA" && metTree.at(sample)->Flag_collisionEventSelectionPA())
                  ) {
                Filter[c][eventType][Form("ALLEVENTS_%s", t.c_str())][s] += 1;
                if (metTree.at(sample)->Flag_collisionEventSelectionPA()         ) { Filter[c][eventType][Form("collisionEventSelectionPA_%s", t.c_str())         ][s] += 1; }
                if (metTree.at(sample)->Flag_collisionEventSelectionPA_rejectPU()) { Filter[c][eventType][Form("collisionEventSelectionPA_rejectPU_%s", t.c_str())][s] += 1; }
                if (metTree.at(sample)->Flag_goodVertices()                      ) { Filter[c][eventType][Form("goodVertices_%s", t.c_str())                      ][s] += 1; }
                if (metTree.at(sample)->Flag_globalTightHalo2016Filter()         ) { Filter[c][eventType][Form("globalTightHalo2016Filter_%s", t.c_str())         ][s] += 1; }
                if (metTree.at(sample)->Flag_HBHENoiseFilter()                   ) { Filter[c][eventType][Form("HBHENoiseFilter_%s", t.c_str())                   ][s] += 1; }
                if (metTree.at(sample)->Flag_HBHENoiseIsoFilter()                ) { Filter[c][eventType][Form("HBHENoiseIsoFilter_%s", t.c_str())                ][s] += 1; }
                if (metTree.at(sample)->Flag_EcalDeadCellTriggerPrimitiveFilter()) { Filter[c][eventType][Form("EcalDeadCellTriggerPrimitiveFilter_%s", t.c_str())][s] += 1; }
                if (metTree.at(sample)->Flag_BadPFMuonFilter()                   ) { Filter[c][eventType][Form("BadPFMuonFilter_%s", t.c_str())                   ][s] += 1; }
                if (metTree.at(sample)->Flag_BadChargedCandidateFilter()         ) { Filter[c][eventType][Form("BadChargedCandidateFilter_%s", t.c_str())         ][s] += 1; }
                if (metTree.at(sample)->Flag_eeBadScFilter()                     ) { Filter[c][eventType][Form("eeBadScFilter_%s", t.c_str())                     ][s] += 1; }
                if (metTree.at(sample)->Flag_badMuons()                          ) { Filter[c][eventType][Form("badMuons_%s", t.c_str())                          ][s] += 1; }
                if (metTree.at(sample)->Flag_duplicateMuons()                    ) { Filter[c][eventType][Form("duplicateMuons_%s", t.c_str())                    ][s] += 1; }

                if (
                    //(metTree.at(sample)->Flag_goodVertices()                      ) &&
                    (metTree.at(sample)->Flag_globalTightHalo2016Filter()         ) &&
                    (metTree.at(sample)->Flag_HBHENoiseFilter()                   ) &&
                    (metTree.at(sample)->Flag_HBHENoiseIsoFilter()                ) &&
                    (metTree.at(sample)->Flag_EcalDeadCellTriggerPrimitiveFilter()) &&
                    (metTree.at(sample)->Flag_BadPFMuonFilter()                   ) &&
                    (metTree.at(sample)->Flag_BadChargedCandidateFilter()         ) &&
                    (metTree.at(sample)->Flag_eeBadScFilter()                     )// &&
                    //(metTree.at(sample)->Flag_badMuons()                          ) &&
                    //(metTree.at(sample)->Flag_duplicateMuons()                    )
                    ){
                  Filter[c][eventType][Form("ALLFILTERS_%s", t.c_str())][s] += 1;
                }
              }
            }
          }
        }
      }
    }
  }
  //
  const std::string fileName = "test";
  std::ofstream file((mainDir + "/" + fileName + ".txt").c_str());
  if (print) {
    for (const auto& c : Filter) {
      for (const auto& e : c.second) {
        for (const auto& f : e.second) {
          for (const auto& s : f.second) {
            if (Filter.at(c.first).at(e.first).at("ALLEVENTS_ONLY").at(s.first)>0) {
              const uint all = ( (f.first.find("ONLY")!=std::string::npos) ? Filter.at(c.first).at(e.first).at("ALLEVENTS_ONLY").at(s.first) : Filter.at(c.first).at(e.first).at("ALLEVENTS_WithPA").at(s.first) );
              std::cout << c.first << "   " << e.first << "    "  << f.first << "    " << s.first << "  EVENTS NOT PASSED:  " << (all - s.second) << " of " << all << Form(" ( %.4f %% ) ", (double(all - s.second)*100./double(all))) << std::endl;
              file << c.first << "   " << e.first << "    "  << f.first << "    " << s.first << "  EVENTS NOT PASSED:  " << (all - s.second) << " of " << all << Form(" ( %.4f %% ) ", (double(all - s.second)*100./double(all))) << std::endl;
            }
          }
          std::cout << "   " << std::endl;
          file << "   " << std::endl;
        }
        std::cout << "   " << std::endl;
        file << "   " << std::endl;
      }
      std::cout << "   " << std::endl;
      file << "   " << std::endl;
    }
    std::cout << "   " << std::endl;
    file << "   " << std::endl;
  }
  file << std::endl; file << std::endl;
  //
}

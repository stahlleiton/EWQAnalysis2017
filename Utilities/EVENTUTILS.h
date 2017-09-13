//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////


#ifndef EVENTUTILS_H_
#define EVENTUTILS_H_

// Auxiliary Headers
#include "../Utilities/HiMETTree.h"
#include "../Utilities/HiMuonTree.h"
// c++ headers
#include <dirent.h>
#include <memory>
#include <iostream>
#include <map>
#include <vector>
#include <string>


namespace PA {
  // Trigger
  enum TRIGGERBIT {
    HLT_PAL2Mu12 = 0,
    HLT_PAL2Mu15 = 1,
    HLT_PAL3Mu3  = 2,
    HLT_PAL3Mu5  = 3,
    HLT_PAL3Mu7  = 4,
    HLT_PAL3Mu12 = 5,
    HLT_PAL3Mu15 = 6,
  };

  // Data Luminosities
  namespace LUMI {
    const double Data_pPb = 110.77;
    const double Data_Pbp = 62.59;
  }

};

// MC Related Information
//
namespace PYQUEN {
  std::map< std::string , std::map< std::string , double > > XSec = {
    { "WToMuNu"   , { { "pPb" , 1.159e4 * 208. * 1.e-3 } , { "Pbp" , 1.223e4 * 208. * 1.e-3 } } },
    { "DYToMuMu"  , { { "pPb" , 1.290e3 * 208. * 1.e-3 } , { "Pbp" , 1.342e3 * 208. * 1.e-3 } } },
    { "WToTauNu"  , { { "pPb" , 1.182e4 * 208. * 1.e-3 } , { "Pbp" , 1.213e4 * 208. * 1.e-3 } } }
  };
};
//
namespace POWHEG {
  std::map< std::string , std::map< std::string , double > > XSec = {
    {"WToMuNu_Plus"       , { { "pPb" , 1213.38 }  , { "Pbp" , 1214.06 } } },
    {"WToMuNu_Minus"      , { { "pPb" , 1082.24 }  , { "Pbp" , 1083.37 } } },
    {"DYToMuMu_M_30_Inf"  , { { "pPb" , 266.28  }  , { "Pbp" , 266.28  } } },
    {"DYToMuMu_M_10_30"   , { { "pPb" , 1182.24 }  , { "Pbp" , 1168.03 } } },
    {"WToTauNu_Plus"      , { { "pPb" , 1146.30 }  , { "Pbp" , 1147.44 } } },
    {"WToTauNu_Minus"     , { { "pPb" , 1026.32 }  , { "Pbp" , 1019.40 } } },
    {"TTall"              , { { "pPb" , 48.68   }  , { "Pbp" , 48.68   } } }
  };
};
//
namespace PYTHIA {
  std::map< std::string , std::map< std::string , double > > XSec = {  
    {"QCDToMu"  , { { "pPb" , (28919.33 * 0.78998) }  , { "Pbp" , 28927.43 * 0.81767 } } } // pPb: 3.67970e+05 * 3.77844e-4 * 208. 48.68  Pbp: 3.67966e+05 * 3.77954e-4 * 208.
  };
};


// Selections
namespace PA {
  //
  bool passEventFilter(const std::unique_ptr<HiMETTree>& metTree)
  {
    if (
        metTree->Flag_collisionEventSelectionPA() // pPb Collision Event Selection
      ) {
      return true; // All Event filters passed
    }
    //
    return false; // At least one event filter failed
  };
  //
  bool checkGenMuon(const ushort& iGenMu, const std::string& sample, const std::unique_ptr<HiMuonTree>& muonTree)
  {
    if (
        ( (sample.find("MC_DYToMuMu")!=std::string::npos) && (muonTree->findMuonMother(iGenMu, 23, 1).pdg==23 || muonTree->findMuonMother(iGenMu, 22, 1).pdg==22) ) || // Z/gamma* -> Muon + Muon
        ( (sample.find("MC_ZToMuMu") !=std::string::npos) && (muonTree->findMuonMother(iGenMu, 23, 1).pdg==23) ) || // Z -> Muon + Muon
        ( (sample.find("MC_WToMuNu") !=std::string::npos) && (muonTree->findMuonMother(iGenMu, 24, 1).pdg==24) ) || // W -> Muon + Neutrino
        ( (sample.find("MC_WToTauNu")!=std::string::npos) && (muonTree->findMuonMother(iGenMu, 24, 2).pdg==24) ) || // W -> Tau -> Muon + Neutrinos
        ( (sample.find("MC_QCDToMu") !=std::string::npos) ) || // QCD -> Muon
        ( (sample.find("MC_TTall")   !=std::string::npos) && (muonTree->findMuonMother(iGenMu, 6, 2).pdg==6 || muonTree->findMuonMother(iGenMu, 24, 1).pdg==24) ) // t -> W -> Muon + Neutrino
        ) {
      return true;
    }
    //
    return false;
  };
  //
  bool selectGenMuon(const ushort& iGenMu, const std::unique_ptr<HiMuonTree>& muonTree)
  {
    if (
        ( std::abs(muonTree->Gen_Muon_Mom()[iGenMu].Eta()) < 2.4 ) && // Consider Generated Muons within the Pseudo-Rapidity acceptance of CMS
        ( muonTree->Gen_Muon_Mom()[iGenMu].Pt() > 25.0           )    // Consider Generated Muons with pT > 25 GeV/c
        ) {
      return true;
    }
    //
    return false;
  };
  //
  bool isGoodMuon(const ushort& iPFMu, const std::unique_ptr<HiMuonTree>& muonTree)
  {
    const short iRecoMu = muonTree->PF_Muon_Reco_Idx()[iPFMu];
    if (iRecoMu < 0) { std::cout << "[ERROR] Reco idx is negative" << std::endl; return false; }
    if (
        ( std::abs(muonTree->PF_Muon_Mom()[iPFMu].Eta()) < 2.4 ) && // Consider Muons within the Pseudo-Rapidity acceptance of CMS
        ( muonTree->Reco_Muon_isTight()[iRecoMu] == true       ) && // Consider Tight Muons
        ( muonTree->PF_Muon_Mom()[iPFMu].Pt() > 25.0           )    // Consider Muons with pT > 25 GeV/c
        ) {
      return true;
    }
    //
    return false;
  };
  //
  bool isIsolatedMuon(const ushort& iPFMu, const std::unique_ptr<HiMuonTree>& muonTree)
  {
    const short iRecoMu = muonTree->PF_Muon_Reco_Idx()[iPFMu];
    if (iRecoMu < 0) { std::cout << "[ERROR] Reco idx is negative" << std::endl; return false; }
    if (
        ( muonTree->PF_Muon_IsoPFR03NoPUCorr()[iPFMu] < 0.15 ) // Consider Isolated Muons
        ) {
      return true;
    }
    //
    return false;
  };
  //
  bool isOfflineMuon(const ushort& iPFMu, const std::unique_ptr<HiMuonTree>& muonTree)
  {
    const short iRecoMu = muonTree->PF_Muon_Reco_Idx()[iPFMu];
    if (iRecoMu < 0) { std::cout << "[ERROR] Reco idx is negative" << std::endl; return false; }
    if (
        ( isGoodMuon(iPFMu, muonTree)     ) && // Consider Good Quality Muons
        ( isIsolatedMuon(iPFMu, muonTree) )    // Consider Isolated Muons
        ) {
      return true;
    }
    //
    return false;
  };
  //
  bool isTriggerMatched(const ushort& triggerIndex, const ushort& iPFMu, const std::unique_ptr<HiMuonTree>& muonTree)
  {
    const short iRecoMu = muonTree->PF_Muon_Reco_Idx()[iPFMu];
    if (
        ( muonTree->Event_Trig_Fired()[triggerIndex] == true       ) && // Consider events that fired the trigger HLT_PAL3Mu12_v1
        ( muonTree->Pat_Muon_Trig()[iRecoMu][triggerIndex] == true )    // Consider muons matched to the online muon that fired the trigger HLT_PAL3Mu12_v1
        ) {
      return true;
    }
    //
    return false;
  };
  //
  bool passDrellYanVeto(const std::unique_ptr<HiMuonTree>& muonTree)
  {
    for (ushort iPFMu1 = 0; iPFMu1 < muonTree->PF_Muon_Mom().size(); iPFMu1++) {
      for (ushort iPFMu2 = 0; iPFMu2 < muonTree->PF_Muon_Mom().size(); iPFMu2++) {
        if (
            ( isOfflineMuon(iPFMu1 , muonTree)            ) && // Consider muons passing offline selection
            ( isOfflineMuon(iPFMu2 , muonTree)            ) && // Consider muons passing offline selection
            ( muonTree->PF_Muon_Mom()[iPFMu1].Pt() > 15.0 ) && // Consider Muons with pT > 15 GeV
            ( muonTree->PF_Muon_Mom()[iPFMu1].Pt() > 15.0 ) && // Consider Muons with pT > 15 GeV
            ( muonTree->PF_Muon_Charge()[iPFMu1] != muonTree->PF_Muon_Charge()[iPFMu2] ) // Consider opposite-sign muons
            ) {
          return false; // Found possible Drell-Yan Candidate, so Drell-Yan Veto failed
        }
      }
    }
    //
    return true; // Drell-Yan veto succeeded
  };
  //
};


// Utiliy Functions

bool existDir(const std::string& dir)
{
  bool exist = false;
  void * dirp = gSystem->OpenDirectory(dir.c_str());
  if (dirp) { gSystem->FreeDirectory(dirp); exist = true; }
  return exist;
};

void makeDir(const std::string& dir)
{
  if (existDir(dir.c_str())==false){ 
    std::cout << "[INFO] DataSet directory: " << dir << " doesn't exist, will create it!" << std::endl;  
    gSystem->mkdir(dir.c_str(), kTRUE);
  }
};



#endif /* EVENTUTILS_H_ */

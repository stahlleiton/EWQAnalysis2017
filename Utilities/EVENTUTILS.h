//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////


#ifndef EVENTUTILS_H_
#define EVENTUTILS_H_

// Auxiliary Headers
#include "../Utilities/HiMETTree.h"
#include "../Utilities/HiMuonTree.h"
// ROOT Headers
#include "TLorentzVector.h"
// c++ headers
#include <sys/stat.h>
#include <dirent.h>
#include <stdio.h>
#include <memory>
#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <chrono>


// Utiliy Functions

bool existDir(const std::string& dir)
{
  bool exist = false;
  void * dirp = gSystem->OpenDirectory(dir.c_str());
  if (dirp) { gSystem->FreeDirectory(dirp); exist = true; }
  return exist;
};

bool existFile(const std::string& file)
{
  struct stat buffer;   
  return (stat (file.c_str(), &buffer) == 0);
};

void makeDir(const std::string& dir)
{
  if (existDir(dir.c_str())==false){ 
    std::cout << "[INFO] DataSet directory: " << dir << " doesn't exist, will create it!" << std::endl;  
    gSystem->mkdir(dir.c_str(), kTRUE);
  }
};

void roundValue( double& value , const uint& nDecimals )
{
  double tmp = value;
  tmp *= std::pow(10.0, nDecimals);
  tmp = std::round(tmp);
  tmp /= std::pow(10.0, nDecimals);
  value = tmp;
};

bool isEqual( const double inVal1 , const double inVal2 , const uint nDecimals )
{
  double val1 = inVal1; roundValue(val1, nDecimals);
  double val2 = inVal2; roundValue(val2, nDecimals);
  if (val1==val2) return true;
  return false;  
};

auto TIME_START_ = std::chrono::high_resolution_clock::now();
auto TIME_END_ = TIME_START_;
int  TIME_DIF_ = 0;
static inline void loadBar(const int iEvent, const int nEvents, const int r = 100, const int w = 100)
{
  // Only update r times.
  if ( iEvent == (nEvents-1) ) { std::cout << std::endl; }
  if ( (iEvent % ((nEvents/r) + 1)) != 0 ) return;
  // Calculuate the ratio of complete-to-incomplete.
  const float ratio = (iEvent / (float)nEvents);
  const int   c     = (ratio * w);
  // Get Time Difference
  TIME_END_   = std::chrono::high_resolution_clock::now();
  TIME_DIF_   = std::chrono::duration_cast<std::chrono::seconds>(TIME_END_ - TIME_START_).count();
  TIME_START_ = std::chrono::high_resolution_clock::now();
  // Show the percentage complete.
  const int sec  = int( double(TIME_DIF_) * (1.0-ratio) *100. );
  const int min  = (sec / 60);
  const int hour = (min / 60);
  printf("[INFO] %3d%% (%02d:%02d:%02d) [", (int)(ratio*100), hour, int(min%60), int(sec%60));
  // Show the load bar.
  for (int i = 0; i < c; i++) { std::cout << "="; }
  for (int i = c; i < w; i++) { std::cout << " "; }
  // ANSI Control codes to go back to the
  // previous line and clear it.
  std::cout << "]\r" << std::flush;
};


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
  };

};


// Variable Computation
namespace PA {
  //
  const double Muon_MASS_ = 0.1057;
  //
  double getWTransverseMass(const double& muPt, const double& muPhi, const double& nuPt, const double& nuPhi)
  {
    TLorentzVector pfMuonP4T = TLorentzVector(), METP4 = TLorentzVector();
    pfMuonP4T.SetPtEtaPhiM(muPt, 0.0, muPhi, Muon_MASS_);
    METP4.SetPtEtaPhiM( nuPt, 0.0, nuPhi, 0.0 );
    TLorentzVector muT = TLorentzVector( pfMuonP4T + METP4 );
    return muT.M();
  };
};

    
// Centre of Mass Frames
namespace PA {
  //
  double EtaLABtoCM(const double& etaLAB, const bool ispPb)
  {
    const double shift = ( ispPb ? 0.465 : -0.465 );
    double etaCM = etaLAB - shift;
    roundValue(etaCM, 4);
    return etaCM;
  };
  //
  double EtaCMtoLAB(const double& etaCM, const bool ispPb)
  {
    const double shift = ( ispPb ? 0.465 : -0.465 );
    double etaLAB = etaCM + shift;
    roundValue(etaLAB, 4);
    return etaLAB;
  };
  //
};

// MC Related Information
//
namespace PYQUEN {
  std::map< std::string , std::map< std::string , double > > XSec = {
    { "WToMuNu"            , { { "pPb" , 1.159e4 * 208. * 1.e-3 } , { "Pbp" , 1.223e4 * 208. * 1.e-3 } } },
    { "DYToMuMu_M_30_Inf"  , { { "pPb" , 1.290e3 * 208. * 1.e-3 } , { "Pbp" , 1.342e3 * 208. * 1.e-3 } } },
    { "WToTauNu"           , { { "pPb" , 1.182e4 * 208. * 1.e-3 } , { "Pbp" , 1.213e4 * 208. * 1.e-3 } } },
    { "TTall"              , { { "pPb" , 48.68                  } , { "Pbp" , 48.68                  } } }
  };
};  
//
namespace POWHEG {
  std::map< std::string , std::map< std::string , double > > XSec = {
    { "WToMuNu_Plus"       , { { "pPb" , 1213.38 }  , { "Pbp" , 1214.06 } } },
    { "WToMuNu_Minus"      , { { "pPb" , 1082.24 }  , { "Pbp" , 1083.37 } } },
    { "DYToMuMu_M_30_Inf"  , { { "pPb" , 266.28  }  , { "Pbp" , 266.28  } } },
    { "DYToMuMu_M_10_30"   , { { "pPb" , 1182.24 }  , { "Pbp" , 1168.03 } } },
    { "WToTauNu_Plus"      , { { "pPb" , 1146.30 }  , { "Pbp" , 1147.44 } } },
    { "WToTauNu_Minus"     , { { "pPb" , 1026.32 }  , { "Pbp" , 1019.40 } } },
    { "TTall"              , { { "pPb" , 48.68   }  , { "Pbp" , 48.68   } } }
  };
};
//
namespace PYTHIA {
  std::map< std::string , std::map< std::string , double > > XSec = {  
    { "QCDToMu"  , { { "pPb" , (28919.33 * 0.78998) }  , { "Pbp" , 28927.43 * 0.81767 } } } // pPb: 3.67970e+05 * 3.77844e-4 * 208. 48.68  Pbp: 3.67966e+05 * 3.77954e-4 * 208.
  };
};
//
// MC Tags
namespace PA {
  enum class MCType 
  {
    Invalid = 0,
      // For Pythia8
      PYTHIA_QCDToMu_pPb           =  10,
      PYTHIA_QCDToMu_Pbp           = -10,
      // For Pyquen
      PYQUEN_WToMuNu_pPb           =  20,
      PYQUEN_WToMuNu_Pbp           = -20,
      PYQUEN_DYToMuMu_M_30_Inf_pPb =  21,
      PYQUEN_DYToMuMu_M_30_Inf_Pbp = -21,
      PYQUEN_WToTauNu_pPb          =  22,
      PYQUEN_WToTauNu_Pbp          = -22,
      PYQUEN_TTall_pPb             =  23,
      PYQUEN_TTall_Pbp             = -23,
      // For Powheg
      POWHEG_WToMuNu_Plus_pPb      =  30,
      POWHEG_WToMuNu_Plus_Pbp      = -30,
      POWHEG_WToMuNu_Minus_pPb     =  31,
      POWHEG_WToMuNu_Minus_Pbp     = -31,
      POWHEG_DYToMuMu_M_30_Inf_pPb =  32,
      POWHEG_DYToMuMu_M_30_Inf_Pbp = -32,
      POWHEG_DYToMuMu_M_10_30_pPb  =  33,
      POWHEG_DYToMuMu_M_10_30_Pbp  = -33,
      POWHEG_WToTauNu_Plus_pPb     =  34,
      POWHEG_WToTauNu_Plus_Pbp     = -34,
      POWHEG_WToTauNu_Minus_pPb    =  35,
      POWHEG_WToTauNu_Minus_Pbp    = -35,
      POWHEG_TTall_pPb             =  36,
      POWHEG_TTall_Pbp             = -36
      };
  //
  const std::map< std::string , int > MCTypeDictionary = {
    { "Invalid"                      , int(MCType::Invalid)},
    // For Pythia8
    { "PYTHIA_QCDToMu_pPb"           , int(MCType::PYTHIA_QCDToMu_pPb)},
    { "PYTHIA_QCDToMu_Pbp"           , int(MCType::PYTHIA_QCDToMu_Pbp)},
    // For Pyquen
    { "PYQUEN_WToMuNu_pPb"           , int(MCType::PYQUEN_WToMuNu_pPb)},
    { "PYQUEN_WToMuNu_Pbp"           , int(MCType::PYQUEN_WToMuNu_Pbp)},
    { "PYQUEN_DYToMuMu_M_30_Inf_pPb" , int(MCType::PYQUEN_DYToMuMu_M_30_Inf_pPb)},
    { "PYQUEN_DYToMuMu_M_30_Inf_Pbp" , int(MCType::PYQUEN_DYToMuMu_M_30_Inf_Pbp)},
    { "PYQUEN_WToTauNu_pPb"          , int(MCType::PYQUEN_WToTauNu_pPb)},
    { "PYQUEN_WToTauNu_Pbp"          , int(MCType::PYQUEN_WToTauNu_Pbp)},
    { "PYQUEN_TTall_pPb"             , int(MCType::PYQUEN_TTall_pPb)},
    { "PYQUEN_TTall_Pbp"             , int(MCType::PYQUEN_TTall_Pbp)},
    // For Powheg
    { "POWHEG_WToMuNu_Plus_pPb"      , int(MCType::POWHEG_WToMuNu_Plus_pPb)},
    { "POWHEG_WToMuNu_Plus_Pbp"      , int(MCType::POWHEG_WToMuNu_Plus_Pbp)},
    { "POWHEG_WToMuNu_Minus_pPb"     , int(MCType::POWHEG_WToMuNu_Minus_pPb)},
    { "POWHEG_WToMuNu_Minus_Pbp"     , int(MCType::POWHEG_WToMuNu_Minus_Pbp)},
    { "POWHEG_DYToMuMu_M_30_Inf_pPb" , int(MCType::POWHEG_DYToMuMu_M_30_Inf_pPb)},
    { "POWHEG_DYToMuMu_M_30_Inf_Pbp" , int(MCType::POWHEG_DYToMuMu_M_30_Inf_Pbp)},
    { "POWHEG_DYToMuMu_M_10_30_pPb"  , int(MCType::POWHEG_DYToMuMu_M_10_30_pPb)},
    { "POWHEG_DYToMuMu_M_10_30_Pbp"  , int(MCType::POWHEG_DYToMuMu_M_10_30_Pbp)},
    { "POWHEG_WToTauNu_Plus_pPb"     , int(MCType::POWHEG_WToTauNu_Plus_pPb)},
    { "POWHEG_WToTauNu_Plus_Pbp"     , int(MCType::POWHEG_WToTauNu_Plus_Pbp)},
    { "POWHEG_WToTauNu_Minus_pPb"    , int(MCType::POWHEG_WToTauNu_Minus_pPb)},
    { "POWHEG_WToTauNu_Minus_Pbp"    , int(MCType::POWHEG_WToTauNu_Minus_Pbp)},
    { "POWHEG_TTall_pPb"             , int(MCType::POWHEG_TTall_pPb)},
    { "POWHEG_TTall_Pbp"             , int(MCType::POWHEG_TTall_Pbp)}
  };
  //
  std::string getMCTypeName(const int& MCTypeID)
    {
      for (const auto& MCT : MCTypeDictionary) { if (MCT.second == MCTypeID) { return MCT.first; } }
      return "Invalid";
    };
  //
  int getMCTypeID(const std::string& MCTypeName)
  {
    if (MCTypeDictionary.count(MCTypeName)>0) { return MCTypeDictionary.at(MCTypeName); }
    return int(MCType::Invalid);
  };
  //
  bool getCrossSection(double& xSection, const std::string& MCTypeName)
  {
    xSection = -1.0;
    std::string GENTAG  = MCTypeName; GENTAG.erase(GENTAG.find("_"), GENTAG.length());
    std::string NAMETAG = MCTypeName; NAMETAG = NAMETAG.substr(NAMETAG.find("_")+1); NAMETAG.erase(NAMETAG.find_last_of("_"), 5);
    std::string COLTAG  = MCTypeName; COLTAG = COLTAG.substr(COLTAG.find_last_of("_")+1);
    if (GENTAG == "POWHEG") {
      if (POWHEG::XSec.count(NAMETAG) == 0) { std::cout << "[ERROR] POWHEG Cross-section was not found for process: " << NAMETAG << std::endl; return false; }
      if (POWHEG::XSec.at(NAMETAG).count(COLTAG) == 0) { std::cout << "[ERROR] POWHEG Cross-section was not found for process: " << (NAMETAG+"_"+COLTAG) << std::endl; return false; }
      xSection = POWHEG::XSec.at(NAMETAG).at(COLTAG);
    }
    else if (GENTAG == "PYTHIA") {
      if (PYTHIA::XSec.count(NAMETAG) == 0) { std::cout << "[ERROR] PYTHIA Cross-section was not found for process: " << NAMETAG << std::endl; return false; }
      if (PYTHIA::XSec.at(NAMETAG).count(COLTAG) == 0) { std::cout << "[ERROR] PYTHIA Cross-section was not found for process: " << (NAMETAG+"_"+COLTAG) << std::endl; return false; }
      xSection = PYTHIA::XSec.at(NAMETAG).at(COLTAG);
    }
    else {
      std::cout << "[ERROR] Cross-section for enerator " << GENTAG << " has not been defined!" << std::endl; return false;
    }
    return true;
  };
  //
  bool getCrossSection(double& xSection, const int& MCTypeID)
  {
    return getCrossSection(xSection, getMCTypeName(MCTypeID));
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
        (sample.find("MC_TT")     == std::string::npos) && (sample.find("MC_DY")  == std::string::npos) &&
        (sample.find("MC_Z")      == std::string::npos) && (sample.find("MC_W")   == std::string::npos) &&
        (sample.find("MC_WToTau") == std::string::npos) && (sample.find("MC_QCD") == std::string::npos)
        )
      { std::cout << "[ERROR] GEN sample " << sample << " has an INVALID NAME!" << std::endl; return false; }
    if (
        ( (sample.find("MC_TT")    !=std::string::npos) && (muonTree->findMuonMother(iGenMu,  6, 2).pdg==6  || muonTree->findMuonMother(iGenMu, 24, 1).pdg==24) ) || // t -> W -> Muon + Neutrino
        ( (sample.find("MC_DY")    !=std::string::npos) && (muonTree->findMuonMother(iGenMu, 23, 1).pdg==23 || muonTree->findMuonMother(iGenMu, 22, 1).pdg==22) ) || // Z/gamma* -> Muon + Muon
        ( (sample.find("MC_Z")     !=std::string::npos) && (muonTree->findMuonMother(iGenMu, 23, 1).pdg==23) ) || // Z -> Muon + Muon
        ( (sample.find("MC_W")     !=std::string::npos) && (muonTree->findMuonMother(iGenMu, 24, 1).pdg==24) ) || // W -> Muon + Neutrino
        ( (sample.find("MC_WToTau")!=std::string::npos) && (muonTree->findMuonMother(iGenMu, 24, 2).pdg==24) ) || // W -> Tau -> Muon + Neutrinos
        ( (sample.find("MC_QCD")   !=std::string::npos) && (muonTree->MuonMother(iGenMu).pdg!=13) ) // QCD -> Muon (at least has one mother)
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
  bool isTightMuon(const ushort& iPFMu, const std::unique_ptr<HiMuonTree>& muonTree)
  {
    const short iRecoMu = muonTree->PF_Muon_Reco_Idx()[iPFMu];
    if (iRecoMu < 0) { std::cout << "[ERROR] Reco idx is negative" << std::endl; return false; }
    if (
        ( std::abs(muonTree->PF_Muon_Mom()[iPFMu].Eta()) < 2.4 ) && // Consider Muons within the Pseudo-Rapidity acceptance of CMS
        ( muonTree->Reco_Muon_isTight()[iRecoMu] == true       )    // Consider Tight Muons
        ) {
      return true;
    }
    //
    return false;
  };
  //
  bool isIsolatedMuon(const ushort& iPFMu, const std::unique_ptr<HiMuonTree>& muonTree)
  {
    if (
        ( muonTree->PF_Muon_IsoPFR03NoPUCorr()[iPFMu] < 0.15 ) // Consider Isolated Muons
        ) {
      return true;
    }
    //
    return false;
  };
  //
  bool isTightIsolatedMuon(const ushort& iPFMu, const std::unique_ptr<HiMuonTree>& muonTree)
  {
    if (
        ( isTightMuon(iPFMu, muonTree)    ) && // Consider Good Quality Muons
        ( isIsolatedMuon(iPFMu, muonTree) )    // Consider Isolated Muons
        ) {
      return true;
    }
    //
    return false;
  };
  //
  bool isGoodMuon(const ushort& iPFMu, const std::unique_ptr<HiMuonTree>& muonTree)
  {
    if (
        ( isTightMuon(iPFMu, muonTree)               ) && // Consider Tight Muons
        ( muonTree->PF_Muon_Mom()[iPFMu].Pt() > 25.0 )    // Consider Muons with pT > 25 GeV/c
        ) {
      return true;
    }
    //
    return false;
  };
  //
  bool isOfflineMuon(const ushort& iPFMu, const std::unique_ptr<HiMuonTree>& muonTree)
  {
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
            ( isTightIsolatedMuon(iPFMu1 , muonTree)          ) && // Consider muons passing Tight ID and isolation
            ( isTightIsolatedMuon(iPFMu2 , muonTree)          ) && // Consider muons passing Tight ID and isolation
            ( muonTree->PF_Muon_Mom()[iPFMu1].Pt() > 15.0     ) && // Consider Muons with pT > 15 GeV
            ( muonTree->PF_Muon_Mom()[iPFMu2].Pt() > 15.0     ) && // Consider Muons with pT > 15 GeV
            ( muonTree->PF_Muon_Charge()[iPFMu1] != muonTree->PF_Muon_Charge()[iPFMu2] ) // Consider opposite-sign muons
            ) {
          return false; // Found possible Drell-Yan Candidate, so Drell-Yan Veto failed
        }
      }
    }
    //
    return true; // Drell-Yan veto succeeded
  };
};



#endif /* EVENTUTILS_H_ */

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////


#ifndef EVENTUTILS_H_
#define EVENTUTILS_H_

// Auxiliary Headers
#include "../Utilities/HiMETTree.h"
#include "../Utilities/HiMuonTree.h"
// ROOT Headers
#include "TLorentzVector.h"
#include "TH1.h"
#include "TGraphAsymmErrors.h"
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

TH1D graphToHist( const TGraphAsymmErrors& gr , const int addOnly=0 )
{
  auto tmp = gr;
  // Remove all unwanted points
  if (addOnly!=0) {
    for (int i = 0; i < tmp.GetN(); ){
      double eta, xSec;
      tmp.GetPoint(i, eta, xSec);
      if      (addOnly<0 && eta>=0) { tmp.RemovePoint(i); }
      else if (addOnly>0 && eta<=0) { tmp.RemovePoint(i); }
      else { i++; }
    }
  }
  const uint nBin = tmp.GetN();
  double bins[nBin+1];
  for (uint iBin = 0; iBin <= nBin; iBin++) {
    double edge;
    if (iBin < nBin) {
      double x, y; tmp.GetPoint(iBin, x, y);
      double ex = tmp.GetErrorXlow(iBin);
      edge = x - ex;
    }
    else {
      double x, y; tmp.GetPoint((nBin-1), x, y);
      double ex = tmp.GetErrorXhigh(nBin-1);
      edge = x + ex;
    }
    bins[iBin] = edge;
  }
  TH1D h(tmp.GetName(), tmp.GetTitle(), nBin, bins);
  for (uint iBin = 0; iBin < nBin; iBin++) {
    double x, y; tmp.GetPoint(iBin, x, y);
    double ey = tmp.GetErrorY(iBin);
    h.SetBinContent(iBin+1, y);
    h.SetBinError(iBin+1, ey);
  }
  tmp.TAttLine::Copy(h);
  tmp.TAttFill::Copy(h);
  tmp.TAttMarker::Copy(h);
  return h;
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
    { "WToMuNu_Plus"        , { { "pPb" , 1213.38 }  , { "Pbp" , 1214.06 } } },
    { "WToMuNu_Minus"       , { { "pPb" , 1082.24 }  , { "Pbp" , 1083.37 } } },
    { "DYToMuMu_M_30_Inf"   , { { "pPb" , 266.28  }  , { "Pbp" , 266.28  } } },
    { "DYToMuMu_M_10_30"    , { { "pPb" , 1182.24 }  , { "Pbp" , 1168.03 } } },
    { "WToTauNu_Plus"       , { { "pPb" , 1146.30 }  , { "Pbp" , 1147.44 } } },
    { "WToTauNu_Minus"      , { { "pPb" , 1026.32 }  , { "Pbp" , 1019.40 } } },
    { "TTall"               , { { "pPb" , 48.68   }  , { "Pbp" , 48.68   } } },
    { "DYToTauTau_M_30_Inf" , { { "pPb" , 259.42  }  , { "Pbp" , 259.42  } } },
    { "DYToTauTau_M_10_30"  , { { "pPb" , 1143.74 }  , { "Pbp" , 1143.74 } } }
  };
};
//
namespace PYTHIA {
  std::map< std::string , std::map< std::string , double > > XSec = {  
    { "QCDToMu"  , { { "pPb" , (28919.33 * 0.6786) }  , { "Pbp" , 28927.43 * 0.6786 } } }, // pPb: 3.67970e+05 * 3.77844e-4 * 208. 48.68  Pbp: 3.67966e+05 * 3.77954e-4 * 208.
    { "WWall"    , { { "pPb" , 12.83 } , { "Pbp" , 12.83 } } },
    { "WZall"    , { { "pPb" , 4.63  } , { "Pbp" , 4.63  } } },
    { "ZZall"    , { { "pPb" , 2.09  } , { "Pbp" , 2.09  } } }
  };
};  
//
namespace EXTERN {
  std::map< std::string , std::map< std::string , double > > XSec = {
    { "TTall"    , { { "pPb" , 45.00   }  , { "Pbp" , 45.00   } } } // CMS Measurement: https://arxiv.org/pdf/1709.07411.pdf
  };
};
//
// MC Tags
namespace PA {
  enum class MCType 
  {
    Invalid = 0,
      // For Pythia8
      PYTHIA_QCDToMu_pPb             =  10,
      PYTHIA_QCDToMu_Pbp             = -10,
      PYTHIA_WWall_pPb               =  11,
      PYTHIA_WWall_Pbp               = -11,
      PYTHIA_WZall_pPb               =  12,
      PYTHIA_WZall_Pbp               = -12,
      PYTHIA_ZZall_pPb               =  13,
      PYTHIA_ZZall_Pbp               = -13,
      // For Pyquen
      PYQUEN_WToMuNu_pPb             =  20,
      PYQUEN_WToMuNu_Pbp             = -20,
      PYQUEN_DYToMuMu_M_30_Inf_pPb   =  21,
      PYQUEN_DYToMuMu_M_30_Inf_Pbp   = -21,
      PYQUEN_WToTauNu_pPb            =  22,
      PYQUEN_WToTauNu_Pbp            = -22,
      PYQUEN_TTall_pPb               =  23,
      PYQUEN_TTall_Pbp               = -23,
      // For Powheg
      POWHEG_WToMuNu_Plus_pPb        =  30,
      POWHEG_WToMuNu_Plus_Pbp        = -30,
      POWHEG_WToMuNu_Minus_pPb       =  31,
      POWHEG_WToMuNu_Minus_Pbp       = -31,
      POWHEG_DYToMuMu_M_30_Inf_pPb   =  32,
      POWHEG_DYToMuMu_M_30_Inf_Pbp   = -32,
      POWHEG_DYToMuMu_M_10_30_pPb    =  33,
      POWHEG_DYToMuMu_M_10_30_Pbp    = -33,
      POWHEG_WToTauNu_Plus_pPb       =  34,
      POWHEG_WToTauNu_Plus_Pbp       = -34,
      POWHEG_WToTauNu_Minus_pPb      =  35,
      POWHEG_WToTauNu_Minus_Pbp      = -35,
      POWHEG_TTall_pPb               =  36,
      POWHEG_TTall_Pbp               = -36,
      POWHEG_DYToTauTau_M_30_Inf_pPb =  37,
      POWHEG_DYToTauTau_M_30_Inf_Pbp = -37,
      POWHEG_DYToTauTau_M_10_30_pPb  =  38,
      POWHEG_DYToTauTau_M_10_30_Pbp  = -38
      };
  //
  // This Dictionary is used to assign a unique tag to each MC sample
  const std::map< std::string , int > MCTypeDictionary = {
    { "Invalid"                        , int(MCType::Invalid)},
    // For Pythia8
    { "PYTHIA_QCDToMu_pPb"             , int(MCType::PYTHIA_QCDToMu_pPb)},
    { "PYTHIA_QCDToMu_Pbp"             , int(MCType::PYTHIA_QCDToMu_Pbp)},
    { "PYTHIA_WWall_pPb"               , int(MCType::PYTHIA_WWall_pPb)},
    { "PYTHIA_WWall_Pbp"               , int(MCType::PYTHIA_WWall_Pbp)},
    { "PYTHIA_WZall_pPb"               , int(MCType::PYTHIA_WZall_pPb)},
    { "PYTHIA_WZall_Pbp"               , int(MCType::PYTHIA_WZall_Pbp)},
    { "PYTHIA_ZZall_pPb"               , int(MCType::PYTHIA_ZZall_pPb)},
    { "PYTHIA_ZZall_Pbp"               , int(MCType::PYTHIA_ZZall_Pbp)},
    // For Pyquen
    { "PYQUEN_WToMuNu_pPb"             , int(MCType::PYQUEN_WToMuNu_pPb)},
    { "PYQUEN_WToMuNu_Pbp"             , int(MCType::PYQUEN_WToMuNu_Pbp)},
    { "PYQUEN_DYToMuMu_M_30_Inf_pPb"   , int(MCType::PYQUEN_DYToMuMu_M_30_Inf_pPb)},
    { "PYQUEN_DYToMuMu_M_30_Inf_Pbp"   , int(MCType::PYQUEN_DYToMuMu_M_30_Inf_Pbp)},
    { "PYQUEN_WToTauNu_pPb"            , int(MCType::PYQUEN_WToTauNu_pPb)},
    { "PYQUEN_WToTauNu_Pbp"            , int(MCType::PYQUEN_WToTauNu_Pbp)},
    { "PYQUEN_TTall_pPb"               , int(MCType::PYQUEN_TTall_pPb)},
    { "PYQUEN_TTall_Pbp"               , int(MCType::PYQUEN_TTall_Pbp)},
    // For Powheg
    { "POWHEG_WToMuNu_Plus_pPb"        , int(MCType::POWHEG_WToMuNu_Plus_pPb)},
    { "POWHEG_WToMuNu_Plus_Pbp"        , int(MCType::POWHEG_WToMuNu_Plus_Pbp)},
    { "POWHEG_WToMuNu_Minus_pPb"       , int(MCType::POWHEG_WToMuNu_Minus_pPb)},
    { "POWHEG_WToMuNu_Minus_Pbp"       , int(MCType::POWHEG_WToMuNu_Minus_Pbp)},
    { "POWHEG_DYToMuMu_M_30_Inf_pPb"   , int(MCType::POWHEG_DYToMuMu_M_30_Inf_pPb)},
    { "POWHEG_DYToMuMu_M_30_Inf_Pbp"   , int(MCType::POWHEG_DYToMuMu_M_30_Inf_Pbp)},
    { "POWHEG_DYToMuMu_M_10_30_pPb"    , int(MCType::POWHEG_DYToMuMu_M_10_30_pPb)},
    { "POWHEG_DYToMuMu_M_10_30_Pbp"    , int(MCType::POWHEG_DYToMuMu_M_10_30_Pbp)},
    { "POWHEG_WToTauNu_Plus_pPb"       , int(MCType::POWHEG_WToTauNu_Plus_pPb)},
    { "POWHEG_WToTauNu_Plus_Pbp"       , int(MCType::POWHEG_WToTauNu_Plus_Pbp)},
    { "POWHEG_WToTauNu_Minus_pPb"      , int(MCType::POWHEG_WToTauNu_Minus_pPb)},
    { "POWHEG_WToTauNu_Minus_Pbp"      , int(MCType::POWHEG_WToTauNu_Minus_Pbp)},
    { "POWHEG_TTall_pPb"               , int(MCType::POWHEG_TTall_pPb)},
    { "POWHEG_TTall_Pbp"               , int(MCType::POWHEG_TTall_Pbp)},
    { "POWHEG_DYToTauTau_M_30_Inf_pPb" , int(MCType::POWHEG_DYToTauTau_M_30_Inf_pPb)},
    { "POWHEG_DYToTauTau_M_30_Inf_Pbp" , int(MCType::POWHEG_DYToTauTau_M_30_Inf_Pbp)},
    { "POWHEG_DYToTauTau_M_10_30_pPb"  , int(MCType::POWHEG_DYToTauTau_M_10_30_pPb)},
    { "POWHEG_DYToTauTau_M_10_30_Pbp"  , int(MCType::POWHEG_DYToTauTau_M_10_30_Pbp)}
  };
  //
  // This Dictionary is used to determine the cross-section used for each MC sample (based on their ID)
  const std::map< int , std::string > MCXSectionDictionary = {
    { int(MCType::Invalid)                        , "Invalid"},
    // For Pythia8
    { int(MCType::PYTHIA_QCDToMu_pPb)             , "PYTHIA_QCDToMu_pPb"},
    { int(MCType::PYTHIA_QCDToMu_Pbp)             , "PYTHIA_QCDToMu_Pbp"},
    { int(MCType::PYTHIA_WWall_pPb)               , "PYTHIA_WWall_pPb"},
    { int(MCType::PYTHIA_WWall_Pbp)               , "PYTHIA_WWall_Pbp"},
    { int(MCType::PYTHIA_WZall_pPb)               , "PYTHIA_WZall_pPb"},
    { int(MCType::PYTHIA_WZall_Pbp)               , "PYTHIA_WZall_Pbp"},
    { int(MCType::PYTHIA_ZZall_pPb)               , "PYTHIA_ZZall_pPb"},
    { int(MCType::PYTHIA_ZZall_Pbp)               , "PYTHIA_ZZall_Pbp"},
    // For Pyquen
    { int(MCType::PYQUEN_WToMuNu_pPb)             , "PYQUEN_WToMuNu_pPb"},
    { int(MCType::PYQUEN_WToMuNu_Pbp)             , "PYQUEN_WToMuNu_Pbp"},
    { int(MCType::PYQUEN_DYToMuMu_M_30_Inf_pPb)   , "PYQUEN_DYToMuMu_M_30_Inf_pPb"},
    { int(MCType::PYQUEN_DYToMuMu_M_30_Inf_Pbp)   , "PYQUEN_DYToMuMu_M_30_Inf_Pbp"},
    { int(MCType::PYQUEN_WToTauNu_pPb)            , "PYQUEN_WToTauNu_pPb"},
    { int(MCType::PYQUEN_WToTauNu_Pbp)            , "PYQUEN_WToTauNu_Pbp"},
    { int(MCType::PYQUEN_TTall_pPb)               , "PYQUEN_TTall_pPb"},
    { int(MCType::PYQUEN_TTall_Pbp)               , "PYQUEN_TTall_Pbp"},
    // For Powheg
    { int(MCType::POWHEG_WToMuNu_Plus_pPb)        , "POWHEG_WToMuNu_Plus_pPb"},
    { int(MCType::POWHEG_WToMuNu_Plus_Pbp)        , "POWHEG_WToMuNu_Plus_Pbp"},
    { int(MCType::POWHEG_WToMuNu_Minus_pPb)       , "POWHEG_WToMuNu_Minus_pPb"},
    { int(MCType::POWHEG_WToMuNu_Minus_Pbp)       , "POWHEG_WToMuNu_Minus_Pbp"},
    { int(MCType::POWHEG_DYToMuMu_M_30_Inf_pPb)   , "POWHEG_DYToMuMu_M_30_Inf_pPb"},
    { int(MCType::POWHEG_DYToMuMu_M_30_Inf_Pbp)   , "POWHEG_DYToMuMu_M_30_Inf_Pbp"},
    { int(MCType::POWHEG_DYToMuMu_M_10_30_pPb)    , "POWHEG_DYToMuMu_M_10_30_pPb"},
    { int(MCType::POWHEG_DYToMuMu_M_10_30_Pbp)    , "POWHEG_DYToMuMu_M_10_30_Pbp"},
    { int(MCType::POWHEG_WToTauNu_Plus_pPb)       , "POWHEG_WToTauNu_Plus_pPb"},
    { int(MCType::POWHEG_WToTauNu_Plus_Pbp)       , "POWHEG_WToTauNu_Plus_Pbp"},
    { int(MCType::POWHEG_WToTauNu_Minus_pPb)      , "POWHEG_WToTauNu_Minus_pPb"},
    { int(MCType::POWHEG_WToTauNu_Minus_Pbp)      , "POWHEG_WToTauNu_Minus_Pbp"},
    { int(MCType::POWHEG_TTall_pPb)               , "EXTERN_TTall_pPb"},
    { int(MCType::POWHEG_TTall_Pbp)               , "EXTERN_TTall_Pbp"},
    { int(MCType::POWHEG_DYToTauTau_M_30_Inf_pPb) , "POWHEG_DYToTauTau_M_30_Inf_pPb"},
    { int(MCType::POWHEG_DYToTauTau_M_30_Inf_Pbp) , "POWHEG_DYToTauTau_M_30_Inf_Pbp"},
    { int(MCType::POWHEG_DYToTauTau_M_10_30_pPb)  , "POWHEG_DYToTauTau_M_10_30_pPb"},
    { int(MCType::POWHEG_DYToTauTau_M_10_30_Pbp)  , "POWHEG_DYToTauTau_M_10_30_Pbp"}
  };
  //
  std::string getMCXSecName(const int& MCTypeID)
    {
      if (MCXSectionDictionary.count(MCTypeID)>0) { return MCXSectionDictionary.at(MCTypeID); }
      return "Invalid";
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
  bool getCrossSection(double& xSection, const std::string& MCXSecName)
  {
    xSection = -1.0;
    std::string GENTAG  = MCXSecName; GENTAG.erase(GENTAG.find("_"), GENTAG.length());
    std::string NAMETAG = MCXSecName; NAMETAG = NAMETAG.substr(NAMETAG.find("_")+1); NAMETAG.erase(NAMETAG.find_last_of("_"), 5);
    std::string COLTAG  = MCXSecName; COLTAG = COLTAG.substr(COLTAG.find_last_of("_")+1);
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
    else if (GENTAG == "PYQUEN") {
      if (PYQUEN::XSec.count(NAMETAG) == 0) { std::cout << "[ERROR] PYQUEN Cross-section was not found for process: " << NAMETAG << std::endl; return false; }
      if (PYQUEN::XSec.at(NAMETAG).count(COLTAG) == 0) { std::cout << "[ERROR] PYQUEN Cross-section was not found for process: " << (NAMETAG+"_"+COLTAG) << std::endl; return false; }
      xSection = PYQUEN::XSec.at(NAMETAG).at(COLTAG);
    }
    else if (GENTAG == "EXTERN") {
      if (EXTERN::XSec.count(NAMETAG) == 0) { std::cout << "[ERROR] External Cross-section was not found for process: " << NAMETAG << std::endl; return false; }
      if (EXTERN::XSec.at(NAMETAG).count(COLTAG) == 0) { std::cout << "[ERROR] External Cross-section was not found for process: " << (NAMETAG+"_"+COLTAG) << std::endl; return false; }
      xSection = EXTERN::XSec.at(NAMETAG).at(COLTAG);
    }
    else {
      std::cout << "[ERROR] Cross-section for generator " << GENTAG << " has not been defined!" << std::endl; return false;
    }
    return true;
  };
  //
  bool getCrossSection(double& xSection, const int& MCTypeID)
  {
    return getCrossSection(xSection, getMCXSecName(MCTypeID));
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
        (sample.find("MC_TT")     == std::string::npos) && (sample.find("MC_DY")  == std::string::npos) && (sample.find("MC_DYToTau") == std::string::npos) &&
        (sample.find("MC_ZZ")     == std::string::npos) && (sample.find("MC_WW")  == std::string::npos) && (sample.find("MC_WZ") == std::string::npos) &&
        (sample.find("MC_Z")      == std::string::npos) && (sample.find("MC_W")   == std::string::npos) &&
        (sample.find("MC_WToTau") == std::string::npos) && (sample.find("MC_QCD") == std::string::npos)
        )
      { std::cout << "[ERROR] GEN sample " << sample << " has an INVALID NAME!" << std::endl; return false; }
    //
    if (sample.find("MC_TT")!=std::string::npos) { // t -> W -> Muon + Neutrino or t -> W -> Tau -> Muon + Neutrino
      if (muonTree->findMuonMother(iGenMu, 24, 2).pdg==24 && muonTree->findMuonMother(iGenMu, 6, 3).pdg==6) { return true; }
    }
    else if (sample.find("MC_DYToTau")!=std::string::npos) { // Z/gamma* -> Tau -> Muon + Muon
      if (muonTree->MuonMother(iGenMu).pdg==15 && (muonTree->findMuonMother(iGenMu, 23, 2).pdg==23 || muonTree->findMuonMother(iGenMu, 22, 2).pdg==22)) { return true; }
    }
    else if (sample.find("MC_DY")!=std::string::npos) { // Z/gamma* -> Muon + Muon
      const auto& mom = muonTree->MuonMother(iGenMu);
      if (mom.pdg==23) { return true; }
      if (mom.pdg==22) { const auto& gMom = muonTree->Mother(mom.idx); if (gMom.idx!=mom.idx && (gMom.pdg==21 || gMom.pdg<10)) { return true; } }
    }
    else if (sample.find("MC_ZZ")!=std::string::npos) { // Z -> X -> Muon
      if (muonTree->MuonMother(iGenMu).pdg==23) { return true; }
    }
    else if (sample.find("MC_Z")!=std::string::npos) { // Z -> Muon + Muon
      if (muonTree->MuonMother(iGenMu).pdg==23) { return true; }
    }
    else if (sample.find("MC_WToTau")!=std::string::npos) { // W -> Tau -> Muon + Neutrinos
      if (muonTree->MuonMother(iGenMu).pdg==15 && muonTree->findMuonMother(iGenMu, 24, 2).pdg==24) { return true; }
    }
    else if (sample.find("MC_WZ")!=std::string::npos) { // Z -> X -> Muon or W -> X -> Muon
      const auto& mom = muonTree->MuonMother(iGenMu); if (mom.pdg==23 || mom.pdg==24) { return true; }
    }
    else if (sample.find("MC_WW")!=std::string::npos) { // W -> X -> Muon
      if (muonTree->MuonMother(iGenMu).pdg==24) { return true; }
    }
    else if (sample.find("MC_W")!=std::string::npos) { // W -> Muon + Neutrino
      if (muonTree->MuonMother(iGenMu).pdg==24) { return true; }
    }
    else if (sample.find("MC_QCD")!=std::string::npos) { // QCD -> Muon (at least has one mother)
      if (muonTree->MuonMother(iGenMu).pdg!=13) { return true; }
    }
    //
    return false;
  };
  //
  int getGenMom(const ushort& iGenMu, const std::string& sample, const std::unique_ptr<HiMuonTree>& muonTree)
  {
    if (
        (sample.find("MC_TT")     == std::string::npos) && (sample.find("MC_DY")  == std::string::npos) && (sample.find("MC_DYToTau") == std::string::npos) &&
        (sample.find("MC_ZZ")     == std::string::npos) && (sample.find("MC_WW")  == std::string::npos) && (sample.find("MC_WZ") == std::string::npos) &&
        (sample.find("MC_Z")      == std::string::npos) && (sample.find("MC_W")   == std::string::npos) &&
        (sample.find("MC_WToTau") == std::string::npos) && (sample.find("MC_QCD") == std::string::npos)
        )
      { std::cout << "[ERROR] GEN sample " << sample << " has an INVALID NAME!" << std::endl; return false; }
    //
    if (sample.find("MC_DYToTau")!=std::string::npos) { // Z/gamma* -> Tau -> Muon + Muon
      if (muonTree->MuonMother(iGenMu).pdg!=15) { return -1; }
      const auto& mom1 = muonTree->findMuonMother(iGenMu, 23, 2); if (mom1.pdg==23) { return mom1.idx; }
      const auto& mom2 = muonTree->findMuonMother(iGenMu, 22, 2); if (mom2.pdg==22) { return mom2.idx; }
    }
    else if (sample.find("MC_DY")!=std::string::npos) { // Z/gamma* -> Muon + Muon
      const auto& mom = muonTree->MuonMother(iGenMu);
      if (mom.pdg==23) { return mom.idx; }
      if (mom.pdg==22) { const auto& gMom = muonTree->Mother(mom.idx); if (gMom.idx!=mom.idx && (gMom.pdg==21 || gMom.pdg<10)) { return mom.idx; } }
    }
    else if (sample.find("MC_ZZ")!=std::string::npos) {
      const auto& mom = muonTree->MuonMother(iGenMu); if (mom.pdg==23) { return mom.idx; }
    }
    else if (sample.find("MC_Z")!=std::string::npos) { // Z -> X -> Muon
      const auto& mom = muonTree->MuonMother(iGenMu); if (mom.pdg==23) { return mom.idx; }
    }
    else if (sample.find("MC_WToTau")!=std::string::npos) { // W -> Tau -> Muon + Neutrinos
      if (muonTree->MuonMother(iGenMu).pdg!=15) { return -1; }
      const auto& mom = muonTree->findMuonMother(iGenMu, 24, 2); if (mom.pdg==24) { return mom.idx; }
    }
    else if (sample.find("MC_WZ")!=std::string::npos) { // Z -> X -> Muon or W -> X -> Muon
      const auto& mom = muonTree->MuonMother(iGenMu); if (mom.pdg==23 || mom.pdg==24) { return mom.idx; }
    }
    else if (sample.find("MC_WW")!=std::string::npos) { // W -> X -> Muon
      const auto& mom = muonTree->MuonMother(iGenMu); if (mom.pdg==24) { return mom.idx; }
    }
    else if (sample.find("MC_W")!=std::string::npos) { // W -> Muon + Neutrino
      const auto& mom = muonTree->MuonMother(iGenMu); if (mom.pdg==24) { return mom.idx; }
    }
    //
    return -1;
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
  bool isIsolatedMuon(const ushort& iPFMu, const std::unique_ptr<HiMuonTree>& muonTree, const double muSF=1.0)
  {
    const auto Muon_Iso = (muonTree->PF_Muon_IsoPFR03NoPUCorr()[iPFMu] / muSF);
    if (
        ( Muon_Iso < 0.15 ) // Consider Isolated Muons
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
  bool isGoodMuon(const ushort& iPFMu, const std::unique_ptr<HiMuonTree>& muonTree, const double muSF=1.0)
  {
    const auto Muon_Pt = (muonTree->PF_Muon_Mom()[iPFMu].Pt() * muSF);
    if (
        ( isTightMuon(iPFMu, muonTree) ) && // Consider Tight Muons
        ( Muon_Pt > 25.0               )    // Consider Muons with pT > 25 GeV/c
        ) {
      return true;
    }
    //
    return false;
  };
  //
  bool isOfflineMuon(const ushort& iPFMu, const std::unique_ptr<HiMuonTree>& muonTree, const double muSF=1.0)
  {
    if (
        ( isGoodMuon(iPFMu, muonTree, muSF)     ) && // Consider Good Quality Muons
        ( isIsolatedMuon(iPFMu, muonTree, muSF) )    // Consider Isolated Muons
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
  //
  bool hasZBoson(const std::unique_ptr<HiMuonTree>& muonTree)
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
          TLorentzVector diMuP4 = muonTree->PF_Muon_Mom()[iPFMu1] + muonTree->PF_Muon_Mom()[iPFMu2];
          if (diMuP4.M() > 80. && diMuP4.M() < 100.) {
            return true; // Found possible Z Boson Candidate
          }
        }
      }
    }
    //
    return false; // No Z Boson Candidate was found
  };
};



#endif /* EVENTUTILS_H_ */

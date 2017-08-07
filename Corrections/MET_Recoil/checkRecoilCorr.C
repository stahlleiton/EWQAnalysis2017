#if !defined(__CINT__) || defined(__MAKECINT__)
// Auxiliary Headers
#include "../../Utilities/HiMETTree.h"
#include "../../Utilities/HiMuonTree.h"
#include "../../Utilities/Histogram2.h"
#include "RecoilCorrector.C"
// ROOT Headers
#include <TFile.h>                    // file handle class
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TDirectory.h>
#include <TSystem.h>
// C++ Headers
#include <iostream>                   // standard I/O
#include <chrono>                     // for time measurement

#endif



void checkRecoilCorr(
                     const std::vector< std::string > METType = { "PF_RAW" , "PF_Type1" , "PF_NoHF_RAW" , "PF_NoHF_Type1" },
                     const std::vector< std::string > COLL = { "PA" }
                     )
{
  // Change the working directory
  const std::string CWD = getcwd(NULL, 0);
  const std::string mainDir = Form("%s/CheckRecoil/", CWD.c_str());
  gSystem->mkdir(Form("%s/CheckRecoil", CWD.c_str()), kTRUE);
  gSystem->ChangeDirectory(Form("%s/CheckRecoil", CWD.c_str()));

  auto t1_ALL = std::chrono::high_resolution_clock::now();
  auto t1 = t1_ALL;
  //
  // Define initial paramteres
  //
  // Kinematic Cuts for W
  const std::vector< double > W_MU_ETA = { -2.4 , -2.2 , -2.0 , -1.8 , -1.6 , -1.4 , -1.2 , -1.0 , -0.8 , -0.6 , -0.4 , -0.2 , 0.0 , 0.2 , 0.4 , 0.6 , 0.8 , 1.0 , 1.2 , 1.4 , 1.6 , 1.8 , 2.0 , 2.2 , 2.4 };
  const double W_MU_PT  = 25.0;
  const double W_MU_ISO = 0.15;
  // Kinematic Cuts for Z
  const std::vector< double > Z_RAP = { -2.4 , -2.1 , -1.6 , -1.2 , -0.8 , -0.3 , 0.3 , 0.8 , 1.2 , 1.6 , 2.1 , 2.4 };
  const std::vector< double > Z_M   = { 80.0 , 100.0 };
  const double Z_MU_PT  = 15.0;
  const double Z_MU_ISO = 0.15;
  // Kinematic Cuts for QCD
  const std::vector< double > QCD_MU_ETA = { -2.4 , -2.2 , -2.0 , -1.8 , -1.6 , -1.4 , -1.2 , -1.0 , -0.8 , -0.6 , -0.4 , -0.2 , 0.0 , 0.2 , 0.4 , 0.6 , 0.8 , 1.0 , 1.2 , 1.4 , 1.6 , 1.8 , 2.0 , 2.2 , 2.4 };
  const double QCD_MU_PT  = 25.0;
  const double QCD_MU_ISO = 0.4;
  //
  // Trigger Info
  const int triggerIndex = 5; // HLT_PAL3SingleMu12_v1
  // MET Info
  const std::map< std::string , std::string > METName = { { "PF_RAW" , "RAW" } , { "PF_Type1" , "TYPE1" } , { "PF_NoHF_RAW" , "NoHF_RAW" } , { "PF_NoHF_Type1" , "NoHF_TYPE1" } };
  const std::map< std::string , std::string > METLbl  = { { "PF_RAW" , "Raw" } , { "PF_Type1" , "Type1" } , { "PF_NoHF_RAW" , "NoHF_Raw" } , { "PF_NoHF_Type1" , "NoHF_Type1" } };
  //
  // Input Files for analysis
  //
  const std::map< std::string , std::vector< std::string > > fileName = {
    {"DATA"             , { "root://cms-xrd-global.cern.ch//store/group/phys_heavyions/anstahll/EWQAnalysis2017/pPb2016/8160GeV/Data/HiEWQForest_PASingleMuon_pPb_Pbp_8160GeV_20170717.root" } },
    {"MC_QCDToMu_pPb"   , { "root://cms-xrd-global.cern.ch//store/group/phys_heavyions/anstahll/EWQAnalysis2017/pPb2016/8160GeV/MC/Embedded/Official/HiEWQForest_Embedded_Official_QCDToMu_pPb_8160GeV_20170717.root"  } },
    {"MC_QCDToMu_Pbp"   , { "root://cms-xrd-global.cern.ch//store/group/phys_heavyions/anstahll/EWQAnalysis2017/pPb2016/8160GeV/MC/Embedded/Official/HiEWQForest_Embedded_Official_QCDToMu_Pbp_8160GeV_20170717.root"  } },
    {"MC_DYToMuMu_pPb"  , { "root://cms-xrd-global.cern.ch//store/group/phys_heavyions/anstahll/EWQAnalysis2017/pPb2016/8160GeV/MC/Embedded/Official/HiEWQForest_Embedded_Official_DYToMuMu_pPb_8160GeV_20170717.root" } },
    {"MC_DYToMuMu_Pbp"  , { "root://cms-xrd-global.cern.ch//store/group/phys_heavyions/anstahll/EWQAnalysis2017/pPb2016/8160GeV/MC/Embedded/Official/HiEWQForest_Embedded_Official_DYToMuMu_Pbp_8160GeV_20170717.root" } },
    {"MC_WToTauNu_pPb"  , { "root://cms-xrd-global.cern.ch//store/group/phys_heavyions/anstahll/EWQAnalysis2017/pPb2016/8160GeV/MC/Embedded/Official/HiEWQForest_Embedded_Official_WToTauNu_pPb_8160GeV_20170717.root" } },
    {"MC_WToTauNu_Pbp"  , { "root://cms-xrd-global.cern.ch//store/group/phys_heavyions/anstahll/EWQAnalysis2017/pPb2016/8160GeV/MC/Embedded/Official/HiEWQForest_Embedded_Official_WToTauNu_Pbp_8160GeV_20170717.root" } },
    {"MC_WToMuNu_pPb"   ,
     {
       "root://cms-xrd-global.cern.ch//store/group/phys_heavyions/anstahll/EWQAnalysis2017/pPb2016/8160GeV/MC/Embedded/Private/HiEWQForest_Embedded_PileUp_CT14_EPPS16_WToMuNu_Plus_pPb_8160GeV_20170717.root" ,
       "root://cms-xrd-global.cern.ch//store/group/phys_heavyions/anstahll/EWQAnalysis2017/pPb2016/8160GeV/MC/Embedded/Private/HiEWQForest_Embedded_PileUp_CT14_EPPS16_WToMuNu_Minus_pPb_8160GeV_20170717.root"
     } 
    },
    {"MC_WToMuNu_Pbp"  ,
     {
       "root://cms-xrd-global.cern.ch//store/group/phys_heavyions/anstahll/EWQAnalysis2017/pPb2016/8160GeV/MC/Embedded/Private/HiEWQForest_Embedded_PileUp_CT14_EPPS16_WToMuNu_Plus_Pbp_8160GeV_20170717.root" ,
       "root://cms-xrd-global.cern.ch//store/group/phys_heavyions/anstahll/EWQAnalysis2017/pPb2016/8160GeV/MC/Embedded/Private/HiEWQForest_Embedded_PileUp_CT14_EPPS16_WToMuNu_Minus_Pbp_8160GeV_20170717.root"
     }
    }
  };
  //
  std::map< std::string , std::vector< std::map< std::string , double > > > sampleInfo = {
    // Based on XSec Average
    {"MC_WToMuNu_pPb"  , 
     {
       { {"Entries" , 984000} , {"XSec" , (5833.54  * 1.000 * 208. * 1.e-3 * (1.0000))} , {"Lumi" , 110.77 } }, // For W plus
       { {"Entries" , 980000} , {"XSec" , (5203.09  * 1.000 * 208. * 1.e-3 * (1.0000))} , {"Lumi" , 110.77 } }  // For W minus
     }
    },
    {"MC_WToMuNu_Pbp" , 
     {
       { {"Entries" , 990000} , {"XSec" , (5833.54  * 1.000 * 208. * 1.e-3 * (1.0000))} , {"Lumi" ,  62.59 } }, // For W plus
       { {"Entries" , 976400} , {"XSec" , (5208.53  * 1.000 * 208. * 1.e-3 * (1.0000))} , {"Lumi" ,  62.59 } }  // For W minus
     }
    },
    {"MC_QCDToMu_pPb"   , { { {"Entries" , 967656} , {"XSec" , (3.67970e+08 * 3.77844e-4 * 208. * 1.e-3 * (0.7745*0.9896))} , {"Lumi" , 110.77 } } } },
    {"MC_QCDToMu_Pbp"   , { { {"Entries" , 968753} , {"XSec" , (3.67966e+08 * 3.77954e-4 * 208. * 1.e-3 * (0.8063*0.9839))} , {"Lumi" , 62.59  } } } },
    {"MC_DYToMuMu_pPb"  , { { {"Entries" , 978253} , {"XSec" , (1284.  * 1.000    * 208. * 1.e-3 * (1.0614*0.9968))}        , {"Lumi" , 110.77 } } } },
    {"MC_DYToMuMu_Pbp"  , { { {"Entries" , 988732} , {"XSec" , (1239.  * 1.000    * 208. * 1.e-3 * (1.1336*0.9964))}        , {"Lumi" , 62.59  } } } },
    {"MC_WToTauNu_pPb"  , { { {"Entries" , 992522} , {"XSec" , (12414. * 1.000    * 208. * 1.e-3 * (1.0000))}               , {"Lumi" , 110.77 } } } },
    {"MC_WToTauNu_Pbp"  , { { {"Entries" , 986860} , {"XSec" , (11256. * 1.000    * 208. * 1.e-3 * (1.0000))}               , {"Lumi" , 62.59  } } } }
  };
  //
  const std::map< std::string , std::vector< std::string > > sampleType = {
    { "sample" , { "DATA" , "MC_WToMuNu" , "MC_QCDToMu" , "MC_WToTauNu" , "MC_DYToMuMu" } },
    { "beam"   , { "pPb"  , "Pbp" } }
  };
  //
  // Create the sample labels
  std::vector< std::string > samples;
  for (const auto & sample : sampleType.at("sample")) {
    if (sample=="DATA" && (fileName.count(sample)>0)) { samples.push_back(sample); }
    else {
      for (auto & beam : sampleType.at("beam")) {
        std::string name = sample + "_" + beam;
        if ((fileName.count(name)>0)) { samples.push_back(name); }
      }
    }
  }
  //
  // Histograms
  //
  const std::map< std::string , std::vector< std::string > > histoType = {
    { "sample" , sampleType.at("sample") },
    { "beam"   , COLL }
  };
  // Create histogram labels
  std::vector< std::pair< std::string , std::string > > histLabel;
  for (const auto & sample : histoType.at("sample")) {
    for (const auto & beam : histoType.at("beam")) {
      histLabel.push_back(std::make_pair( sample , beam ));
    }
  }
  // 
  const std::map< std::string , std::vector< double > > objBin = {
    { "WToMu"    , W_MU_ETA   },
    { "ZToMuMu"  , Z_RAP      },
    { "QCDToMu"  , QCD_MU_ETA }
  };
  // For titles
  std::map< std::string , std::map< std::string , std::map< std::string , struct VarInfo > > > varInfo;
  for (const auto& obj : objBin) {
    for (const auto& met : METType) {
      // Eta Bins
      std::string label = "";
      const std::string etaLbl = ( obj.first=="ZToMuMu" ? "zRap" : "muEta" );
      const uint endIdx = obj.second.size()-1;
      for (uint i = 0; i < endIdx; i++) {
        label = Form("%s_NOCORR_%s_%.0f_%.0f_MET_Mod", METName.at(met).c_str(), etaLbl.c_str(), obj.second[i]*10., obj.second[i+1]*10.);
        varInfo[met][obj.first][label] = { Form("PF %s E^{Miss}_{T} Magnitud (GeV/c)", METLbl.at(met).c_str()) , 80 ,   0. , 160. , { false }  };
        label = Form("%s_SCALED_%s_%.0f_%.0f_MET_Mod", METName.at(met).c_str(), etaLbl.c_str(), obj.second[i]*10., obj.second[i+1]*10.);
        varInfo[met][obj.first][label] = { Form("PF %s E^{Miss}_{T} Magnitud (GeV/c)", METLbl.at(met).c_str()) , 80 ,   0. , 160. , { false }  };
        label = Form("%s_SMEARED_%s_%.0f_%.0f_MET_Mod", METName.at(met).c_str(), etaLbl.c_str(), obj.second[i]*10., obj.second[i+1]*10.);
        varInfo[met][obj.first][label] = { Form("PF %s E^{Miss}_{T} Magnitud (GeV/c)", METLbl.at(met).c_str()) , 80 ,   0. , 160. , { false }  };
      }
      // Inclusive Bin
      label = Form("%s_NOCORR_%s_%.0f_%.0f_MET_Mod", METName.at(met).c_str(), etaLbl.c_str(), obj.second[0]*10., obj.second[endIdx]*10.);
      varInfo[met][obj.first][label] = { Form("PF %s E^{Miss}_{T} Magnitud (GeV/c)", METLbl.at(met).c_str()) , 80 ,   0. , 160. , { false }  };
      label = Form("%s_SCALED_%s_%.0f_%.0f_MET_Mod", METName.at(met).c_str(), etaLbl.c_str(), obj.second[0]*10., obj.second[endIdx]*10.);
      varInfo[met][obj.first][label] = { Form("PF %s E^{Miss}_{T} Magnitud (GeV/c)", METLbl.at(met).c_str()) , 80 ,   0. , 160. , { false }  };
      label = Form("%s_SMEARED_%s_%.0f_%.0f_MET_Mod", METName.at(met).c_str(), etaLbl.c_str(), obj.second[0]*10., obj.second[endIdx]*10.);
      varInfo[met][obj.first][label] = { Form("PF %s E^{Miss}_{T} Magnitud (GeV/c)", METLbl.at(met).c_str()) , 80 ,   0. , 160. , { false }  };
    }
  }
  // For Plot Text
  std::map< std::string , struct TypeInfo > typeInfo = {
    { "WToMu"    , 
      { { "WToMuNu" , "WToTauNu" }   , {Form("p^{#mu}_{T}>%.0f GeV/c, Iso^{#mu}<%.2f", W_MU_PT, W_MU_ISO), "#mu Tight ID , Drell-Yan Veto"} }
    },
    { "QCDToMu"  , 
      { { "QCDToMu" }   ,{Form("p^{#mu}_{T}>%.0f GeV/c, Iso^{#mu}>%.2f", QCD_MU_PT, QCD_MU_ISO), "#mu Tight ID , Drell-Yan Veto"} }
    },
    { "DYToMuMu" , 
      { { "DYToMuMu" } , {Form("p^{#mu}_{T}>%.0f GeV/c, |#eta^{#mu}|<2.4, Iso^{#mu}<%.2f", Z_MU_PT, Z_MU_ISO), "#mu Tight ID, Opposite Sign", "|M^{#mu#mu} - 90| > 10 GeV/c^{2}"} }
    },
    { "ZToMuMu"  , 
      { { "DYToMuMu" } , {Form("p^{#mu}_{T}>%.0f GeV/c, |#eta^{#mu}|<2.4, Iso^{#mu}<%.2f", Z_MU_PT, Z_MU_ISO), "#mu Tight ID, Opposite Sign", "80 < M^{#mu^{+}#mu^{-}} < 100 GeV/c^{2}"} }
    }
  };
  //
  // Build the Histograms
  //
  std::map< std::string , Histogram2 > hist;
  // Add the different histograms
  for (const auto & h : histLabel) {
    const std::string sample = h.first + "_" + h.second;
    for (const auto & met : METType) {
      for (const auto & var : varInfo.at(met)) {
        hist[met].Book(sample, var.first, var.second);
      }
    }
  }
  //
  // Get the Recoil Corrections
  //
  RecoilCorrector recoilCorr = RecoilCorrector();
  for (const auto& met : METType) {
    const std::string fileName_MC   = Form("%s/FitRecoil/MC_DYToMuMu_PYQUEN/MET_%s/PA/Results/fits_RecoilPDF_%s_PA.root", CWD.c_str(), met.c_str(), met.c_str());
    const std::string fileName_DATA = Form("%s/FitRecoil/DATA/MET_%s/PA/Results/fits_RecoilPDF_%s_PA.root", CWD.c_str(), met.c_str(), met.c_str());
    recoilCorr.setInputFiles(met, fileName_MC, fileName_DATA);
  };
  //
  // ------------------------------------------------------------------------------------------------------------------------
  //
  // Extract all the samples
  std::map< std::string , Long64_t > nentries;
  std::map< std::string , std::unique_ptr< HiMuonTree > > muonTree;
  std::map< std::string , std::unique_ptr< HiMETTree  > > metTree , metNoHFTree;
  for (const auto & sample : samples) {
    //
    muonTree[sample] = std::unique_ptr<HiMuonTree>(new HiMuonTree());
    if (!muonTree.at(sample)->GetTree(fileName.at(sample))) return;
    nentries[sample] = muonTree.at(sample)->GetEntries();
    //
    metNoHFTree[sample] = std::unique_ptr<HiMETTree>(new HiMETTree());
    if (!metNoHFTree.at(sample)->GetTree(fileName.at(sample), "metAnaNoHF")) return;
    if (metNoHFTree.at(sample)->GetEntries() != nentries.at(sample)) { std::cout << "[ERROR] Inconsistent number of entries between MET NoHF (" << 
        metNoHFTree.at(sample)->GetEntries() << ") and Muon Tree (" << muonTree.at(sample)->GetEntries() << ") !" << std::endl; return; }
    //
    metTree[sample] = std::unique_ptr<HiMETTree>(new HiMETTree());
    if (!metTree.at(sample)->GetTree(fileName.at(sample), "metAna")) return;
    if (metTree.at(sample)->GetEntries() != nentries.at(sample)) { std::cout << "[ERROR] Inconsistent number of entries between MET (" << 
        metTree.at(sample)->GetEntries() << ") and Muon Tree (" << muonTree.at(sample)->GetEntries() << ") !" << std::endl; return; }
    //
    /*
    if (sample.find("MC_")!=std::string::npos) {
      Long64_t entry = 0;
      TChain* parent = muonTree.at(sample)->Tree();
      Long64_t* offsets = parent->GetTreeOffset();
      for(int itree = 0; itree < parent->GetNtrees(); itree++){
        entry += offsets[itree];
        parent->LoadTree(entry);
        sampleInfo[sample][itree].at("Entries") = parent->GetTree()->GetEntries();
      }
    }
    */
  }
  //
  // ------------------------------------------------------------------------------------------------------------------------
  //
  // Loop over the samples
  //
  for (auto & sample : samples) {
    //
    // Determine if the sample is Data
    bool isData = false; if (sample.find("DATA")!=std::string::npos) { isData = true; }
    //
    if (!isData) { recoilCorr.setInitialSetup(sample); }
    //
    for (Long64_t jentry = 0; jentry < nentries.at(sample); jentry++) {
      //
      // Get the entry in the trees
      if (metTree.at(sample)->GetEntry(jentry)<0) break;
      if (metNoHFTree.at(sample)->GetEntry(jentry)<0) break;
      if (muonTree.at(sample)->GetEntry(jentry)<0) break;
      const int itree = muonTree.at(sample)->Tree()->GetTreeNumber();
      //
      if (jentry%200000==0) std::cout << sample << " : " << jentry << "/" << nentries[sample] << std::endl;
      if (jentry%200000==0) { 
        std::cout << "[INFO] Processing time: " << std::chrono::duration_cast<std::chrono::seconds>( std::chrono::high_resolution_clock::now() - t1 ).count() << " sec" << std::endl;
        t1 = std::chrono::high_resolution_clock::now();
      }
      // 
      // Check that the different tree agrees well
      if (muonTree.at(sample)->Event_Run()    != metTree.at(sample)->Event_Run()        ) { std::cout << "[ERROR] MET Run does not agree!"  << std::endl; return; }
      if (muonTree.at(sample)->Event_Number() != metTree.at(sample)->Event_Number()     ) { std::cout << "[ERROR] MET Event does not agree!"      << std::endl; return; }
      if (muonTree.at(sample)->Event_Run()    != metNoHFTree.at(sample)->Event_Run()    ) { std::cout << "[ERROR] MET NoHF Run does not agree!"   << std::endl; return; }
      if (muonTree.at(sample)->Event_Number() != metNoHFTree.at(sample)->Event_Number() ) { std::cout << "[ERROR] MET NoHF Event does not agree!" << std::endl; return; }
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
      // Determine the type of sample
      std::string evtType = sample;
      if (!isData) { evtType = evtType.substr(0, (evtType.find(evtCol)-1)); }
      //
      // Get the Lumi re-weight for MC
      double w = ( isData ? 1.0 : ( (sampleInfo.at(sample)[itree].at("XSec")*sampleInfo.at(sample)[itree].at("Lumi")) / sampleInfo.at(sample)[itree].at("Entries") ) );
      //
      // Event Based Information
      //
      // Apply Event Filters
      if (metTree.at(sample)->Flag_collisionEventSelectionPA()==false) continue;    // PA Event Selection
      //
      // Check Trigger Fired
      if (muonTree.at(sample)->Event_Trig_Fired()[triggerIndex]==false) continue;
      //
      // Determine the MET
      std::map< std::string , TVector2 > MET;
      for (const auto& met : METType) {
        TVector2 evtMET = TVector2();
        if (met == "PF_RAW"        ) { evtMET = metTree.at(sample)->PF_MET_NoShift_Mom();        }
        if (met == "PF_Type1"      ) { evtMET = metTree.at(sample)->Type1_MET_NoShift_Mom();     }
        if (met == "PF_NoHF_RAW"   ) { evtMET = metNoHFTree.at(sample)->PF_MET_NoShift_Mom();    }
        if (met == "PF_NoHF_Type1" ) { evtMET = metNoHFTree.at(sample)->Type1_MET_NoShift_Mom(); }
        MET[met] = evtMET;
      }
      //
      // Muon Based Information
      //
      // Determine the Muon pT
      std::vector< TLorentzVector > mu_rawP4;
      for (uint imu = 0; imu < muonTree.at(sample)->PF_Muon_Mom().size(); imu++) {
        mu_rawP4.push_back( muonTree.at(sample)->PF_Muon_Mom()[imu] );
        // Proceed to invert the pseudo-rapidity of muons when collision is Pb-p and user wants PA
      }
      std::vector< TLorentzVector > mu_P4 = mu_rawP4;
      //
      // Determine the muon Isolation
      std::vector< double > mu_Iso;
      for (uint imu = 0; imu < mu_P4.size(); imu++) { mu_Iso.push_back( muonTree.at(sample)->PF_Muon_IsoPFR03NoPUCorr()[imu] ); }
      //
      /*
        if (applyMuonCorr) {
        for (uint imu = 0; imu < mu_P4.size(); imu++) {
        // Correct the muon pT
        corrMuonMom(mu_P4[imu], imu, muonTree.at(sample), !isData, rc);
        // Correct the muon isolation
        mu_Iso[imu] *= ( (mu_rawP4[imu].Pt() + 1E-12) / (mu_P4[imu].Pt() + 1E-12) );
        // Correct the MET
        for (auto& met : MET) { met.second += ( TVector2( mu_rawP4[imu].Px() , mu_rawP4[imu].Py() ) - TVector2( mu_P4[imu].Px() , mu_P4[imu].Py() ) ); }
        }
        }
      */
      //
      // Fix the muon PF - GEN matching
      std::vector<char> PF_Muon_Gen_Idx , Gen_Muon_PF_Idx;
      if (!isData) { muonTree.at(sample)->GetUniquePFGenMuonMatching(PF_Muon_Gen_Idx, Gen_Muon_PF_Idx, muonTree.at(sample)->PF_Muon_Gen_Idx()); }
      //
      // Find the Leading pT Muon
      float maxPt = -99.; int maxIdx = -1;
      for (uint imu = 0; imu < muonTree.at(sample)->PF_Muon_Mom().size(); imu++) {
        const short imuR = muonTree.at(sample)->PF_Muon_Reco_Idx()[imu];
        if (imuR == -1) { std::cout << "[ERROR] Reco idx is -1" << std::endl; return; }
        // Only consider Muons Matched to GEN in MC
        if (!isData && PF_Muon_Gen_Idx[imu]==-1) continue;
        if (sample.find("MC_WTo")!=std::string::npos      && muonTree.at(sample)->MuonMother(PF_Muon_Gen_Idx[imu]).pdg!=24 ) continue;
        if (sample.find("MC_DYToMuMu")!=std::string::npos && muonTree.at(sample)->MuonMother(PF_Muon_Gen_Idx[imu]).pdg!=23 ) continue;
        // Only consider Muons Matched to Trigger
        if (muonTree.at(sample)->Pat_Muon_Trig().at(imuR)[triggerIndex]==false) continue;
        // Only consider Tight Muons
        if (muonTree.at(sample)->Reco_Muon_isTight()[imuR]==false) continue;
          // Only consider Muons within the Pseudo-Rapidity acceptance of CMS
        if (abs(mu_P4[imu].Eta()) >= 2.4) continue;
        if (maxPt < mu_P4[imu].Pt()) { maxPt = mu_P4[imu].Pt(); maxIdx = imu; }
      }
      if (maxIdx==-1) continue;
      //
      // Correct the MET
      //
      // Determine the reference and boson pT vectors,
      if (!isData && !recoilCorr.getPtFromTree(maxIdx, muonTree.at(sample))) { return; }
      //
      std::map< std::string , TVector2 > MET_CORR_SCALING , MET_CORR_SMEARING;
      for (const auto& met : METType) {
        if (isData) {
          MET_CORR_SCALING[met] = MET[met];
          MET_CORR_SMEARING[met] = MET[met];
        }
        else {
          if (!recoilCorr.correctMET(MET_CORR_SMEARING[met] , MET[met], met, "Smearing")) { std::cout << "[ERROR] Failed to correct the MET using Smearing method!" << std::endl; return; }
          if (!recoilCorr.correctMET(MET_CORR_SCALING[met]  , MET[met], met, "Scaling" )) { std::cout << "[ERROR] Failed to correct the MET using Scaling method!"  << std::endl; return; }
        }
      }
      //
      // Classify the events
      //
      std::string eventType = "Other";
      //
      // Event Type: Drell-Yan -> Muon + Muon
      //
      TLorentzVector Z_P4;
      uint iZ = 0; double zM_diff = 99999999999999999.;
      if ( eventType=="Other" ) {
        for (uint iQQ = 0; iQQ < muonTree.at(sample)->PF_DiMuon_Charge().size(); iQQ++) {
          //
          ushort imu1 = muonTree.at(sample)->PF_DiMuon_Muon1_Idx()[iQQ];
          ushort imu2 = muonTree.at(sample)->PF_DiMuon_Muon2_Idx()[iQQ];
          short imu1R = muonTree.at(sample)->PF_Muon_Reco_Idx()[imu1];
          short imu2R = muonTree.at(sample)->PF_Muon_Reco_Idx()[imu2];
          //
          if (
              ( abs(mu_P4[imu1].Eta()) < 2.4 && abs(mu_P4[imu2].Eta()) < 2.4 ) &&   // Check that both muons are within the ETA acceptance
              ( mu_P4[imu1].Pt() > Z_MU_PT && mu_P4[imu2].Pt() > Z_MU_PT     ) &&   // Check that both muons have pt larger than 15
              ( mu_Iso[imu1] < Z_MU_ISO && mu_Iso[imu2] < Z_MU_ISO           ) &&   // Select isolated muons
              ( muonTree.at(sample)->Reco_Muon_isTight()[imu1R] && muonTree.at(sample)->Reco_Muon_isTight()[imu2R]) // Require both muons to pass Tight ID
              )
            {
              TLorentzVector diMu_P4 = mu_P4[imu1] + mu_P4[imu2];
              //
              if (diMu_P4.M() <= 30.) continue; // Select Dimuons with mass > 30 GeV to match with MC
              if (eventType=="Other") eventType = "DYToMuMu_ALL";
              if (
                  ( muonTree.at(sample)->PF_DiMuon_Charge()[iQQ] == 0      ) && // Select opposite sign dimuons
                  ( diMu_P4.M() > Z_M[0] && diMu_P4.M() < Z_M[1] )    // Select invariant mass within Z mass range
                  )
                {
                  // Find the best Z cadidate
                  if (std::abs(diMu_P4.M()-91.2) < zM_diff) { zM_diff = std::abs(diMu_P4.M()-91.2); iZ = iQQ; Z_P4 = diMu_P4; }
                  eventType = "ZToMuMu";
                }
            }
        }
      }
      //
      // Event Type: QCD -> Muon + X
      //
      if ( eventType=="Other" ) {
        if (
            ( mu_Iso[maxIdx]    >= QCD_MU_ISO ) && // Select Non-Isolated muons
            ( mu_P4[maxIdx].Pt() > QCD_MU_PT  )    // Consider only leading muons with pT larger or equal than 25 GeV
            )
          {
            eventType = "QCDToMu";
          }
      }
      //
      // Event Type: W -> Muon + Neutrino
      //
      if ( eventType=="Other" ) {
        if (
            ( mu_Iso[maxIdx]     < W_MU_ISO ) &&  // Select Isolated muons
            ( mu_P4[maxIdx].Pt() > W_MU_PT  )     // Consider only leading muons with pT larger or equal than 25 GeV
            )
          {
            eventType = "WToMu";
          }
      }
      //
      // Loop over each type of histogram
      for (const auto & h : histLabel) {
        const std::string name = h.first + "_" + h.second;
        //
        // Match the sample type with the one in the histogram
        if (h.first!=evtType) continue;
        //
        // Match the beam direction with the one in the histogram
        if (h.second!=evtCol && h.second!="PA") continue;
        //
        // Proceed to invert the pseudo-rapidity of muons when collision is Pb-p and user wants PA
        // NOTE: Since pZ = pT*sinh(eta) then inverting the sign of eta is the same as inverting the sign of pZ
        //
        std::vector< TLorentzVector > MU_P4 = mu_P4;
        if (h.second=="PA" && evtCol=="Pbp") { for (auto& p4 : MU_P4) {  p4.SetPz( -1.0*p4.Pz() ); } }
        //
        // Proceed to invert the pseudo-rapidity of dimuons when collision is Pb-p and user wants PA
        // NOTE: The sum of 2 muons with pZ sign-inverted is the same as inverting the pZ of the dimuon since pZ_DIMU = pZ_MU1 + pZ_MU2
        //
        TLorentzVector diMU_P4 = Z_P4;
        if (h.second=="PA" && evtCol=="Pbp") { diMU_P4.SetPz( -1.0*diMU_P4.Pz() ); }
        //
        // Fill the histograms
        //
        if ( eventType == "WToMu" || eventType == "ZToMuMu" || eventType == "QCDToMu" ) {
          for (const auto& met : METType) {
            std::map< std::string , std::map< std::string , std::pair < float , float > > > valueMap;
            // Eta bins
            std::string label = "";
            const std::string etaLbl = ( eventType=="ZToMuMu" ? "zRap" : "muEta" );
            const uint endIdx = objBin.at(eventType).size()-1;
            for (uint i = 0; i < endIdx; i++) {
              if ( 
                  ( eventType=="WToMu"   && ( MU_P4[maxIdx].Eta()>objBin.at("WToMu")[i]   && MU_P4[maxIdx].Eta()<=objBin.at("WToMu")[i+1]   ) ) ||
                  ( eventType=="QCDToMu" && ( MU_P4[maxIdx].Eta()>objBin.at("QCDToMu")[i] && MU_P4[maxIdx].Eta()<=objBin.at("QCDToMu")[i+1] ) ) ||
                  ( eventType=="ZToMuMu" && ( diMU_P4.Rapidity()>objBin.at("ZToMuMu")[i]  && diMU_P4.Rapidity()<=objBin.at("ZToMuMu")[i+1] ) )
                   )
                {
                  label = Form("%s_NOCORR_%s_%.0f_%.0f_MET_Mod", METName.at(met).c_str(), etaLbl.c_str(), objBin.at(eventType)[i]*10., objBin.at(eventType)[i+1]*10.);
                  valueMap[eventType][label] = std::make_pair( MET.at(met).Mod() , w );
                  label = Form("%s_SCALED_%s_%.0f_%.0f_MET_Mod", METName.at(met).c_str(), etaLbl.c_str(), objBin.at(eventType)[i]*10., objBin.at(eventType)[i+1]*10.);
                  valueMap[eventType][label] = std::make_pair( MET_CORR_SCALING.at(met).Mod() , w );
                  label = Form("%s_SMEARED_%s_%.0f_%.0f_MET_Mod", METName.at(met).c_str(), etaLbl.c_str(), objBin.at(eventType)[i]*10., objBin.at(eventType)[i+1]*10.);
                  valueMap[eventType][label] = std::make_pair( MET_CORR_SMEARING.at(met).Mod() , w );
                }
            }
            // Integrated bin
            if ( 
                ( eventType=="WToMu"   && ( MU_P4[maxIdx].Eta()>objBin.at("WToMu")[0]   && MU_P4[maxIdx].Eta()<=objBin.at("WToMu")[objBin.at("WToMu").size()-1]   ) ) ||
                ( eventType=="QCDToMu" && ( MU_P4[maxIdx].Eta()>objBin.at("QCDToMu")[0] && MU_P4[maxIdx].Eta()<=objBin.at("QCDToMu")[objBin.at("QCDToMu").size()-1] ) ) ||
                ( eventType=="ZToMuMu" && ( diMU_P4.Rapidity()>objBin.at("ZToMuMu")[0]  && diMU_P4.Rapidity()<=objBin.at("ZToMuMu")[objBin.at("ZToMuMu").size()-1] ) )
                 )
              {
                label = Form("%s_NOCORR_%s_%.0f_%.0f_MET_Mod", METName.at(met).c_str(), etaLbl.c_str(), objBin.at(eventType)[0]*10., objBin.at(eventType)[endIdx]*10.);
                valueMap[eventType][label] = std::make_pair( MET.at(met).Mod() , w );
                label = Form("%s_SCALED_%s_%.0f_%.0f_MET_Mod", METName.at(met).c_str(), etaLbl.c_str(), objBin.at(eventType)[0]*10., objBin.at(eventType)[endIdx]*10.);
                valueMap[eventType][label] = std::make_pair( MET_CORR_SCALING.at(met).Mod() , w );
                label = Form("%s_SMEARED_%s_%.0f_%.0f_MET_Mod", METName.at(met).c_str(), etaLbl.c_str(), objBin.at(eventType)[0]*10., objBin.at(eventType)[endIdx]*10.);
                valueMap[eventType][label] = std::make_pair( MET_CORR_SMEARING.at(met).Mod() , w );
              }
            for (const auto& v : valueMap) hist.at(met).Fill(name, v.first, v.second);
          }
        }
      } // End of loop of histograms
    } // End of loop over events
  } // End of loop over samples
  //
  // Draw the histograms
  for (const auto& met : METType) {
    gSystem->mkdir(Form("%s/CheckRecoil/%s", CWD.c_str(), met.c_str()), kTRUE);
    gSystem->ChangeDirectory(Form("%s/CheckRecoil/%s", CWD.c_str(), met.c_str()));
    hist.at(met).Draw("DATAvsMCStack", typeInfo);
    gSystem->ChangeDirectory(CWD.c_str());
  }
  //
  std::cout << "[INFO] Total processing time: " << std::chrono::duration_cast<std::chrono::minutes>( std::chrono::high_resolution_clock::now() - t1_ALL ).count() << " minutes" << std::endl;
}

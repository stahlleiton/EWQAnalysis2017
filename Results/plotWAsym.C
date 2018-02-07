// Auxiliary Headers
#include "Utilities/resultsUtils.h"
#include "resultsEWQ2tree.C"
// ROOT headers
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
// RooFit headers
// c++ headers
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <string>


/////////////////////
// OTHER FUNCTIONS //
/////////////////////

void extractInfo ( VarBinMap& inputVar , bool& useEtaCM , const std::string& workDirName , const std::string& colTag , const std::string& metTag ,
                   const std::string& dsTag , const std::string& thePoiName );
void getResult   ( BinPentaMap& var , BinSextaMapVec& systVar , VarBinMap& inputVar , const bool&  useEtaCM , const std::string& workDirName , const std::string& effType , const std::string& accType ,
		   const bool doSyst = true , const bool isNominal = true , const VarBinMap nomVar = {} , const uint systCorr = 0 );

/////////////////////


void plotWAsym(
               const std::string nominalWorkDirName = "NominalCM",
               const bool doSyst = true,
               const std::vector< std::string > collVec = { "PA" /*, "pPb" , "Pbp"  */},
               const std::string effType = "TnP",
               const std::string accType = ""
               )
{
  //
  // Define the Systematic varations
  // Variation   Method: (0 or 4): Fully Uncorrelated , 1: Eta Correlated , 2: Charge Correlated , 3: Fully Correlated
  WSDirMap workDirNames = {
    { "Nominal" , {
	{ "Nominal"  , { { nominalWorkDirName } , 0 } }
      }
    },
    // SIGNAL
    { "Signal" , {
	{ "METRange" , { { "NominalCM_METMax300" } , 0 } },
	{ "BinWidth" , { { "NominalCM_BinWidth1" , "NominalCM_BinWidth3" } , 3 } },
      }
    },
    // BACKGROUND: QCD
    { "QCD_Background" , {
	{ "QCD_ConstrainMean_Avg" , { { "SystematicCM_QCD_Constrain_Mean" } , 0 } },
	{ "QCD_MultiJet"          , { { "SystematicCM_QCD_MultiJet"       } , 0 } },
	{ "QCD_Iso0p08"           , { { "SystematicCM_QCD_Iso0p08"        } , 0 } },
      }
    },
    // BACKGROUND: ELECTROWEAK
    { "EWK_Background" , {
	{ "XSection_DrellYan" , { { "SystematicCM_XSECTION_DY_DOWN"     , "SystematicCM_XSECTION_DY_UP"     } , 3 } },
	{ "XSection_WToTau"   , { { "SystematicCM_XSECTION_WToTau_DOWN" , "SystematicCM_XSECTION_WToTau_UP" } , 3 } },
	{ "XSection_TTbar"    , { { "SystematicCM_XSECTION_TTbar_DOWN"  , "SystematicCM_XSECTION_TTbar_UP"  } , 3 } },
      }
    },
    // CORRECTION
    { "Recoil_Correction" , {
	{ "Recoil_Smearing" , { { "NominalCM_RecoilSmearing" } , 0 } },
	{ "Recoil_SystPtFunc" , { { "NominalCM_RecoilSystPtFunc" } , 0 } },
      }
    },
    // CORRECTION
    { "Event_Activity" , {
	{ "HF_NTrackCorr"  , { { "NominalCM_NTrackCorr" } , 0 } },
      }
    }
    // OTHER STUFF
    //{ "Other" , {
	//{ "Recoil_NoCorr" , { { "NominalCM_NoRecoilCorr" } , 3 } },
	//{ "NoCorr"        , { { "NominalCM_NoCorr"       } , 3 } },
	//{ "HF_Only"       , { { "NominalCM_HFCorrOnly"   } , 3 } },
	//{ "TnP_NoCorr"    , { { "NominalCM_NoTnPCorr"    } , 3 } },
	//{ "TnP_Only"    , { { "NominalCM_TnPCorrOnly"    } , 0 } },
	//{ "Recoil_Only" , { { "NominalCM_RecoilCorrOnly" } , 0 } },
	//{ "HF_NoCorr"   , { { "NominalCM_NoHFCorr"       } , 0 } },
    //}
    //}
  };
  workDirNames.at("Recoil_Correction")["Recoil_StatVar_MC"] = std::make_pair(std::vector<std::string>(), 0);
  for (uint i = 0; i < 100; i++) { workDirNames.at("Recoil_Correction").at("Recoil_StatVar_MC").first.push_back(Form("NominalCM_RecoilStatVar_MC/Variation_%d", i)); }
  workDirNames.at("Recoil_Correction")["Recoil_StatVar_DATA"] = std::make_pair(std::vector<std::string>(), 0);
  for (uint i = 0; i < 100; i++) { workDirNames.at("Recoil_Correction").at("Recoil_StatVar_DATA").first.push_back(Form("NominalCM_RecoilStatVar_DATA/Variation_%d", i)); }
  //
  const std::string metTag      = "METPF_RAW";
  const std::string dsTag       = "DATA";
  const std::string thePoiNames = "all";
  bool useEtaCM = false;
  //
  // Initialize the Result and Systematic variables
  BinSextaMap var;
  BinSeptaMapVec systVar;
  //
  // Get the Nominal Result
  VarBinMap inputVarNom;
  if (workDirNames.count("Nominal")>0) {
    const auto& wkDir = workDirNames.at("Nominal").at("Nominal").first[0];
    std::cout << "[INFO] Adding results for : " << wkDir << std::endl;
    for (const auto& colTag : collVec) { extractInfo( inputVarNom , useEtaCM , wkDir , colTag , metTag , dsTag , thePoiNames ); }
    getResult( var["Nominal"] , systVar["Efficiency"] , inputVarNom , useEtaCM , wkDir , effType , accType , doSyst );
  }
  if (var.count("Nominal")==0) { std::cout << "[ERROR] Nominal results are missing" << std::endl; return; }
  //
  // --------------------------------------------------------------------------------- //
  //
  // Print the Raw Yield tables
  const std::string CWD = getcwd(NULL, 0);
  const std::string outDir = CWD + "/Output/" + nominalWorkDirName+"/" + metTag+"/" + dsTag;
  if (!printYieldsTables(inputVarNom, outDir)) { return; }
  //
  // Print the Result tables
  if (!printResultsTables(var.at("Nominal"), inputVarNom, outDir)) { return; }
  //
  // --------------------------------------------------------------------------------- //
  //
  // Get the Systematic Results
  if (doSyst) {
    // Add all systematic uncertainty
    for (const auto& cat : workDirNames) {
      if (cat.first=="Nominal") continue;
      for (const auto& lbl : cat.second) {
	const auto& corrType = lbl.second.second;
	for (const auto& wkDir : lbl.second.first) {
	  std::cout << "[INFO] Adding results for : " << wkDir << std::endl;
	  VarBinMap inputVar;
	  for (const auto& colTag : collVec) { extractInfo( inputVar , useEtaCM , wkDir , colTag , metTag , dsTag , thePoiNames ); }
	  BinPentaMap tmpVar; BinSextaMapVec tmpSystVar;
	  getResult( tmpVar , tmpSystVar , inputVar , useEtaCM , wkDir , effType , accType , true , false , inputVarNom , corrType);
	  systVar[cat.first][lbl.first].push_back(tmpVar);
	}
      }
    }
    //
    // Compute all systematic uncertainties
    computeSystematic(var, systVar);
    //
    // --------------------------------------------------------------------------------- //
    // Print Covariance Matrix
    if (!printCovarianceMatrix(systVar, workDirNames, outDir, useEtaCM)) { return; }
    //
    // --------------------------------------------------------------------------------- //
    // Print systematic Tables
    //
    // Nominal Systematics
    if (!printSystematicTables(var.at("Nominal"), outDir , useEtaCM)) { return; }
    // Efficiency Systematics
    if (!printEffSystematicTables(systVar.at("Efficiency"), outDir , useEtaCM)) { return; }
    // Full Systematics
    if (!printFullSystematicTable(var, outDir, useEtaCM)) { return; }
  }
  //
  // Create the main plots
  GraphPentaMap graph;
  iniResultsGraph(graph, var);
  if (!fillResultsGraph(graph, var)) { return; }
  drawGraph(graph, outDir, useEtaCM, accType, effType);
  //
  // Create the systematic plots
  if (doSyst) {
    GraphPentaMap systGraph;
    for (const auto& sVar : systVar) { iniResultsGraph(systGraph, sVar.second); }
    for (const auto& sVar : systVar) { if (!fillResultsGraph(systGraph, sVar.second)) { return; } }
    drawGraph(systGraph, outDir, useEtaCM, accType, effType);
    drawSystematicGraph(systGraph, outDir, useEtaCM, accType, effType);
    drawSystematicGraph(graph, outDir, useEtaCM, accType, effType);
    drawCombineSystematicGraph(graph, outDir, useEtaCM, accType, effType);
  }
  //
  // Create the plots with predictions
  drawGraphWithTheory(graph, outDir, useEtaCM, accType, effType);
  //
  // Create the plots with 5 TeV results
  drawGraphWithpPb5TeV(graph, outDir, useEtaCM, accType, effType);
  //
};


void extractInfo(
                 VarBinMap& inputVar,
                 bool&  useEtaCM,
                 const std::string& workDirName,
                 const std::string& colTag,
                 const std::string& metTag,
                 const std::string& dsTag,
                 const std::string& thePoiNames
                 )
{
  //
  // --------------------------------------------------------------------------------- //
  //
  // Define the input file info
  const std::string CWD = getcwd(NULL, 0);
  const std::string inputDirPath = Form("%s/Tree/%s/%s/%s/W/%s", CWD.c_str(), workDirName.c_str(), metTag.c_str(), dsTag.c_str(), colTag.c_str());
  const std::string inputFileName = "tree_allvars.root";
  const std::string inputFilePath = Form("%s/%s", inputDirPath.c_str(), inputFileName.c_str());
  //
  // Open the input file
  std::unique_ptr<TFile> inputFile;
  if (existFile(inputFilePath)) { inputFile.reset(new TFile(inputFilePath.c_str(), "READ")); }
  if (!inputFile || inputFile->IsOpen()==false || inputFile->IsZombie()==true) {
    std::cout << "[WARNING] The input file " << inputFilePath << " was not found, will create it!" << std::endl; if (inputFile) { inputFile->Close(); }
    //
    std::string preCWD = CWD; preCWD.erase(preCWD.find_last_of("/"), 10);
    const std::string wsDirPath = Form("%s/Fitter/Output/%s/%s/%s/W/%s/result", preCWD.c_str(), workDirName.c_str(), metTag.c_str(), dsTag.c_str(), colTag.c_str());
    if (existDir(wsDirPath)==false) { std::cout << "[WARNING] Workspace directory " << wsDirPath << " was not found, will skip it!" << std::endl; return; }
    //
    if (!resultsEWQ2tree(workDirName, metTag, dsTag, colTag)) { return; };
    inputFile.reset(new TFile(inputFilePath.c_str(), "READ"));
    if (inputFile->IsOpen()==false || inputFile->IsZombie()==true) {
      std::cout << "[ERROR] The input file " << inputFilePath << " could not be re-created!" << std::endl; return;
    }
  }
  //
  // --------------------------------------------------------------------------------- //
  //
  // Define the tree info container
  TreeInfo info;
  // Initialize the tree info container
  iniResultsTreeInfo(info, thePoiNames);
  //
  // Extract the input tree
  TTree * tree = (TTree*) inputFile->Get("fitResults");
  if (tree==NULL) { std::cout << "[ERROR] The input tree fitResults was not found in " << inputFilePath << "" << std::endl; inputFile->Close(); return; }
  // Set the address of the tree branches
  setBranchAddress(*tree, info);
  //
  // Extract the input variables
  const uint nEntries = tree->GetEntries();
  for (uint i = 0; i < nEntries; i++) {
    tree->GetEntry(i);
    //
    // Determine the bin
    double etaMin = info.Var.at("VAR_Muon_Eta").at("Min");
    double etaMax = info.Var.at("VAR_Muon_Eta").at("Max");
    if ((etaMin == -99.) || (etaMax == -99.)) { std::cout << "[ERROR] The bin was not set properly!" << std::endl; return; }
    // Ignore inclusive bins
    if (std::abs(etaMax - etaMin) > 2.0) continue;
    //
    const std::string collSystem = *info.StrP.at("collSystem");
    const std::string charge     = *info.StrP.at("charge");
    //
    // If the bins where fitted in the Center of Mass, change from LAB to CM
    if (useEtaCM && !info.Flag.at("useEtaCM")) { std::cout << "[ERROR] The useEtaCM was true before but the current workspace is not CM, please don't mix them!" << std::endl; return; }
    if (info.Flag.at("useEtaCM")) {
      bool usepPb = false; if ( (collSystem == "PA") || (collSystem == "pPb") ) { usepPb = true; }
      etaMin = PA::EtaLABtoCM(std::max(-2.4, etaMin), usepPb);
      etaMax = PA::EtaLABtoCM(std::min(+2.4, etaMax), usepPb);
      useEtaCM = true;
    }
    // Round the bin boundaries to three decimals
    roundValue(etaMin, 3);
    roundValue(etaMax, 3);
    const anabin<0> bin(etaMin, etaMax);
    //
    for (auto& v : info.Var) {
      if (v.first.find("POI_")!=std::string::npos) {
        // Get the fit variables
        std::string name = v.first; name.erase(name.find("POI_"),4);
        inputVar[collSystem][bin][charge][name]["Val"] = v.second.at("Val");
        inputVar[collSystem][bin][charge][name]["Err_Stat_High"] = ( (v.second.at("ErrHi")>0) ? v.second.at("ErrHi") : v.second.at("Err") );
        inputVar[collSystem][bin][charge][name]["Err_Stat_Low" ] = ( (v.second.at("ErrLo")>0) ? v.second.at("ErrLo") : v.second.at("Err") );
        inputVar[collSystem][bin][charge][name]["Err_Syst_High"] = 0.0;
        inputVar[collSystem][bin][charge][name]["Err_Syst_Low" ] = 0.0;
      }
      else if (v.first.find("VAR_")!=std::string::npos) {
        // Get the dataset variables
        std::string name = v.first; name.erase(name.find("VAR_"),4);
        inputVar[collSystem][bin][charge][name]["Val"] = v.second.at("Val");
        inputVar[collSystem][bin][charge][name]["RMS"] = v.second.at("Err");
        inputVar[collSystem][bin][charge][name]["Min"] = v.second.at("Min");
        inputVar[collSystem][bin][charge][name]["Max"] = v.second.at("Max");
      }
      else {
        // Get the remaining variables
        const std::string name = v.first;
        inputVar[collSystem][bin][charge][name]["Val"] = v.second.at("Val");
      }
    }
    // Store if we want to use EtaCM
    inputVar[collSystem][bin][charge]["useEtaCM"]["Val"] = double(useEtaCM);
  }
  //
  //
  inputFile->Close();
  //
};


void getResult(
               BinPentaMap&       var,
               BinSextaMapVec&    systVar,
               VarBinMap&         inputVar,
               const bool&        useEtaCM,
               const std::string& workDirName,
               const std::string& effType,
               const std::string& accType,
               const bool         doSyst,
               const bool         isNominal,
               const VarBinMap    nomVar,
	       const uint         systCorr
               )
{
  //
  // --------------------------------------------------------------------------------- //
  //
  if (accType!="" || effType!="") {
    //
    // Initialize the Acceptance and Efficiency container
    BinPentaMap eff;
    iniAcceptanceAndEfficiency(eff, inputVar);
    //
    // Define the efficiency input file name
    const std::string CWD = getcwd(NULL, 0);
    std::string preCWD = CWD; preCWD.erase(preCWD.find_last_of("/"), 10);
    const bool useCM = (workDirName.find("CM")!=std::string::npos);
    std::string effWorkDirName = "";
    if (workDirName.find("CutAndCount")!=std::string::npos) { effWorkDirName = "CutAndCount"; }
    else { effWorkDirName = "Nominal"; }
    if (useCM) { effWorkDirName += "CM"; }
    uint useHFCorr = 1;
    if ( (workDirName.find("TnPCorrOnly")!=std::string::npos) || (workDirName.find("RecoilCorrOnly")!=std::string::npos) || (workDirName.find("NoHFCorr")!=std::string::npos) ) { useHFCorr = 0; }
    if (workDirName.find("NTrackCorr")!=std::string::npos) { useHFCorr = 2; }
    if      (useHFCorr==1) { effWorkDirName += "_WithHF";     }
    else if (useHFCorr==2) { effWorkDirName += "_WithNTrack"; }
    const std::string effDirPath  = Form("%s/Efficiency/Output/%s", preCWD.c_str(), effWorkDirName.c_str());
    const std::string effFileName = "efficiencyTnP.root";
    const std::string effFilePath = Form("%s/%s", effDirPath.c_str(), effFileName.c_str());
    //
    // Extract the Acceptance and Efficiency
    if (!getAcceptanceAndEfficiency(eff, effFilePath, useEtaCM, isNominal)) { return; }
    //
    // --------------------------------------------------------------------------------- //
    //
    // Proceed to correct the Raw Yields
    //
    if (!correctRawYields(inputVar, eff, accType, effType, isNominal)) { return; }
  }
  //
  // --------------------------------------------------------------------------------- //
  //
  // Compute the Charge Asymmetry
  //
  if (!computeChargeAsymmetry(var, systVar, inputVar, doSyst, nomVar, systCorr)) { return; }
  //
  // Compute the Forward-Backward ratio
  //
  if (!computeForwardBackwardRatio(var, systVar, inputVar, doSyst, nomVar, systCorr)) { return; }
  //
  // Compute the Cross-Section
  //
  if (!computeCrossSection(var, systVar, inputVar, doSyst, nomVar, systCorr)) { return; }
  //
};

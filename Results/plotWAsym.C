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

void plotWAsym(
               const std::string workDirName = "NominalCM",
               const std::vector< std::string > collVec = { "PA" },
               const std::string effType     = "TnP",
               const std::string accType     = "MC"
               )
{
  //
  const std::string metTag      = "METPF_RAW";
  const std::string dsTag       = "DATA";
  const std::string thePoiNames = "all";
  //
  for (const auto& colTag : collVec) {
    //
    // --------------------------------------------------------------------------------- //
    //
    // Define the input file info
    const std::string CWD = getcwd(NULL, 0);
    const std::string inputDirPath = Form("%s/Tree/%s/%s/%s/%s", CWD.c_str(), workDirName.c_str(), metTag.c_str(), dsTag.c_str(), colTag.c_str());
    const std::string inputFileName = "tree_allvars.root";
    const std::string inputFilePath = Form("%s/%s", inputDirPath.c_str(), inputFileName.c_str());
    // Open the input file
    TFile* inputFileP = new TFile(inputFilePath.c_str(), "READ");
    if (inputFileP->IsOpen()==false || inputFileP->IsZombie()==true) {
      std::cout << "[WARNING] The input file " << inputFilePath << " was not found, will create it!" << std::endl; inputFileP->Close(); delete inputFileP;
      if (!resultsEWQ2tree(workDirName, metTag, dsTag, colTag)) { return; };
      inputFileP = new TFile(inputFilePath.c_str(), "READ");
      if (inputFileP->IsOpen()==false || inputFileP->IsZombie()==true) {
        std::cout << "[ERROR] The input file " << inputFilePath << " could not be re-created!" << std::endl; return;
      }
    }
    TFile& inputFile = *inputFileP;
    //
    // --------------------------------------------------------------------------------- //
    //
    // Define the tree info container
    TreeInfo info;
    // Initialize the tree info container
    iniResultsTreeInfo(info, thePoiNames);
    //
    // Extract the input tree
    TTree * tree = (TTree*) inputFile.Get("fitResults");
    if (tree==NULL) { std::cout << "[ERROR] The input tree fitResults was not found in " << inputFilePath << "" << std::endl; inputFile.Close(); return; }
    // Set the address of the tree branches
    setBranchAddress(*tree, info);
    //
    // Extract the input variables
    VarBinMap inputVar;
    bool useEtaCM = false;
    //
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
      if (info.Flag.at("useEtaCM")) {
        bool usepPb = false; if ( (collSystem == "PA") || (collSystem == "pPb") ) { usepPb = true; }
        etaMin = PA::EtaLABtoCM(etaMin, usepPb);
        etaMax = PA::EtaLABtoCM(etaMax, usepPb);
        useEtaCM = true;
      }
      // Round the bin boundaries to two decimals
      roundValue(etaMin, 2);
      roundValue(etaMax, 2);
      anabin<0> bin(etaMin, etaMax);
      //
      // Get the fit variables
      for (auto& v : info.Var) {
        if (v.first.find("POI_")==std::string::npos) continue;
        std::string name = v.first; name.erase(name.find("POI_"),4);
        inputVar[collSystem][bin][charge][name]["Val"] = v.second.at("Val");
        inputVar[collSystem][bin][charge][name]["Err_Stat_High"] = v.second.at("Err");
        inputVar[collSystem][bin][charge][name]["Err_Stat_Low"]  = v.second.at("Err");
        inputVar[collSystem][bin][charge][name]["Err_Syst_High"] = 0.0;
        inputVar[collSystem][bin][charge][name]["Err_Syst_Low"]  = 0.0;
      }
      // Get the dataset variables
      for (auto& v : info.Var) {
        if (v.first.find("VAR_")==std::string::npos) continue;
        std::string name = v.first; name.erase(name.find("VAR_"),4);
        inputVar[collSystem][bin][charge][name]["Val"] = v.second.at("Val");
        inputVar[collSystem][bin][charge][name]["RMS"] = v.second.at("Err");
        inputVar[collSystem][bin][charge][name]["Min"] = v.second.at("Min");
        inputVar[collSystem][bin][charge][name]["Max"] = v.second.at("Max");
      }
    }
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
      std::string preCWD = CWD; preCWD.erase(preCWD.find_last_of("/"), 10);
      const std::string effDirPath  = Form("%s/Efficiency/TnPEfficiency/%s", preCWD.c_str(), workDirName.c_str());
      const std::string effFileName = "efficiencyTnP.root";
      const std::string effFilePath = Form("%s/%s", effDirPath.c_str(), effFileName.c_str());
      //
      // Extract the Acceptance and Efficiency
      if (!getAcceptanceAndEfficiency(eff, effFilePath, useEtaCM)) { return; }
      //
      // --------------------------------------------------------------------------------- //
      //
      // Proceed to correct the Raw Yields
      //
      if (!correctRawYields(inputVar, eff, accType, effType)) { return; }
    }
    //
    // --------------------------------------------------------------------------------- //
    //
    // Initialize the output variables
    BinPentaMap var;
    //
    const std::vector< std::string > varType = { "Var" , "Err_Stat_High" , "Err_Stat_Low" , "Err_Syst_High" , "Err_Syst_Low" };
    //
    // Compute the Charge Asymmetry
    //
    if (!computeChargeAsymmetry(var, inputVar)) { return; }
    //
    // Compute the Forward-Backward ratio
    //
    if (!computeForwardBackwardRatio(var, inputVar)) { return; }
    //
    // Compute the Cross-Section
    //
    if (!computeCrossSection(var, inputVar)) { return; }
    //
    // --------------------------------------------------------------------------------- //
    //
    // Initialize the Output Graphs
    GraphQuadMap graph;
    iniResultsGraph(graph, var);
    //
    // Fill the Output Graphs
    //
    if (!fillResultsGraph(graph, var)) { return; }
    //
    // --------------------------------------------------------------------------------- //
    //
    // Draw the Output Graphs
    drawGraph(graph, CWD, useEtaCM, accType, effType);
    //
    //
    inputFile.Close();
    //
    // Clean Memory
    delete inputFileP;
  }
}

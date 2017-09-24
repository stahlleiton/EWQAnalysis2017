#ifndef resultsEWQ2tree_C
#define resultsEWQ2tree_C

// Auxiliary Headers
#include "Utilities/resultsUtils.h"
// ROOT headers
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
// RooFit headers
#include "RooWorkspace.h"
#include "RooAbsData.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooFitResult.h"
#include "RooArgSet.h"
// c++ headers
#include <iostream>
#include <string>



bool resultsEWQ2tree(
                     const std::string workDirName = "NominalCM",
                     const std::string metTag      = "METPF_RAW",
                     const std::string dsTag       = "DATA",
                     const std::string colTag      = "PA",
                     const std::string thePoiNames = "all"
                     )
{
  //
  // --------------------------------------------------------------------------------- //
  //
  // Define the output file info
  const std::string CWD = getcwd(NULL, 0);
  const std::string outputDirPath = Form("%s/Tree/%s/%s/%s/%s", CWD.c_str(), workDirName.c_str(), metTag.c_str(), dsTag.c_str(), colTag.c_str());
  const std::string outputFileName = "tree_allvars.root";
  const std::string outputFilePath = Form("%s/%s", outputDirPath.c_str(), outputFileName.c_str());
  //
  // --------------------------------------------------------------------------------- //
  //
  // Define the tree info container
  TreeInfo info;
  // Initialize the tree info container
  iniResultsTreeInfo(info, thePoiNames);
  //
  // Initialize the tree
  TTree tree("fitResults", "Fit Results");
  //
  // Set the tree branches
  setBranches(tree, info);
  //
  // --------------------------------------------------------------------------------- //
  //
  // Get the list of input files
  //
  std::vector< std::string > inputFileNames;
  std::string preCWD = CWD; preCWD.erase(preCWD.find_last_of("/"), 10);
  const std::string inputDirPath = Form("%s/Fitter/Output/%s/%s/%s/%s/result", preCWD.c_str(), workDirName.c_str(), metTag.c_str(), dsTag.c_str(), colTag.c_str());
  if (!fileList(inputFileNames, inputDirPath)) { return false; };
  //
  // --------------------------------------------------------------------------------- //
  //
  // Loop over the input files
  //
  for (const auto& inputFileName : inputFileNames) {
    //
    std::cout << "Processing file: " << inputFileName << std::endl;
    //
    // Open input file
    const std::string inputFilePath = Form("%s/%s", inputDirPath.c_str(), inputFileName.c_str());
    TFile inputFile(inputFilePath.c_str(), "READ");
    //
    if (inputFile.IsOpen()==false || inputFile.IsZombie()==true) {
      std::cout << "[ERROR] The input file " << inputFilePath << " could not be created!" << std::endl;
    }
    //
    // Extract the Workspace
    RooWorkspace* ws = (RooWorkspace*) inputFile.Get("workspace");
    if (ws == NULL) { std::cout << "[Error] File: " << inputFilePath << " does not have the workspace!" << std::endl; inputFile.Close(); return false; }
    //
    // Extract the information from the workspace
    const std::string DSTAG = (ws->obj("DSTAG"))     ? ((TObjString*)ws->obj("DSTAG"))->GetString().Data()     : "";
    const std::string CHA   = (ws->obj("channel"))   ? ((TObjString*)ws->obj("channel"))->GetString().Data()   : "";
    const std::string COL   = (ws->obj("fitSystem")) ? ((TObjString*)ws->obj("fitSystem"))->GetString().Data() : "";
    const std::string CHG   = (ws->obj("fitCharge")) ? ((TObjString*)ws->obj("fitCharge"))->GetString().Data() : "";
    const std::string OBJ   = (ws->obj("fitObject")) ? ((TObjString*)ws->obj("fitObject"))->GetString().Data() : "";
    // Check the information
    if (DSTAG.find(dsTag)==std::string::npos) { std::cout << "[ERROR] Workspace DSTAG " << DSTAG << " is not consistent with input dsTag " << dsTag << std::endl; inputFile.Close(); return false; }
    if (COL != colTag ) { std::cout << "[ERROR] Workspace COL " << COL << " is not consistent with input colTag " << colTag << std::endl; inputFile.Close(); return false; }
    if (OBJ != "W"    ) { std::cout << "[ERROR] Only W fits are currently supported in result macros!"      << std::endl; return false; }
    if (CHA != "ToMu" ) { std::cout << "[ERROR] Only muon channel is currently supported in result macros!" << std::endl; return false; }
    //
    // Fill the information
    info.Str.at("collSystem") = COL;
    info.Str.at("charge"    ) = CHG;
    info.Str.at("fitObject" ) = OBJ;
    info.Str.at("channel"   ) = CHA;
    info.Str.at("metType"   ) = ( (ws->obj("METType")) ? ((TObjString*)ws->obj("METType"))->GetString().Data() : "" );
    //
    const std::string token  = ( CHG + "_" + COL );
    const std::string tag    = ( OBJ + CHA + token );
    const std::string dsName = ( "d" + CHG + "_" + DSTAG );
    //
    // Fill the Model Information
    const std::vector< std::string > objType = { "W" , "WToTau" , "DY" , "TTbar" , "QCD" };
    for (const auto& o : objType) {
      const std::string modelLabel = Form("Model_%s%s%s", o.c_str(), CHA.c_str(), token.c_str());
      info.Str.at("Model_"+o) = ( (ws->obj(modelLabel.c_str())) ? ((TObjString*)ws->obj(modelLabel.c_str()))->GetString().Data() : "" );
    }
    //
    // Fill the Cut Information
    info.Str.at("cutAndCount_W") = ( (ws->obj(Form("CutAndCount_%s", tag.c_str()))) ? ((TObjString*)ws->obj(Form("CutAndCount_%s", tag.c_str())))->GetString().Data() : "" );
    //
    // Fill the Flag Information
    bool useEtaCM = false;
    if ( (ws->var("useEtaCM") != NULL) && (ws->var("useEtaCM")->getVal() == 1.0) ) { useEtaCM = true; }
    info.Flag.at("useEtaCM") = useEtaCM;
    //
    // Fill the Dataset Variable Information
    RooAbsData* ds = (RooAbsData*) ws->data(Form("CutAndCount_%s", dsName.c_str()));
    if (ds == NULL) { ds = (RooAbsData*) ws->data(Form("%s", dsName.c_str())); }
    if (ds == NULL) { std::cout << "[ERROR] The Dataset " << dsName << " was not found in the workspace" << std::endl; inputFile.Close(); return false; }
    for (auto& v : info.Var) {
      if (v.first.find("VAR_")==std::string::npos) continue;
      std::string varName = v.first; varName.erase(varName.find("VAR_"), 4);
      //
      RooRealVar* var  = (RooRealVar*) ws->var(varName.c_str());
      RooRealVar* mean = ( (var && ds) ? (RooRealVar*) ds->meanVar(*var) : NULL );
      RooRealVar* rms  = ( (var && ds) ? (RooRealVar*) ds->rmsVar(*var) : NULL );
      //
      if (v.second.count("Min")) { v.second.at("Min") = var  ? var->getMin()   : -1.0; }
      if (v.second.count("Max")) { v.second.at("Max") = var  ? var->getMax()   : -1.0; }
      if (v.second.count("Val")) { v.second.at("Val") = mean ? mean->getVal()  : -1.0; }
      if (v.second.count("Err")) { v.second.at("Err") = rms  ? rms->getError() : -1.0; }
    }
    //
    // Get the Snapshots
    const RooArgSet *parIni = ws->getSnapshot("initialParameters");
    //
    // Fill the Fitted Variable Information
    for (auto& p : info.Var) {
      if (p.first.find("POI_")==std::string::npos) continue;
      std::string poiName = p.first; poiName.erase(poiName.find("POI_"), 4);
      poiName += token;
      //
      if (ws->var(poiName.c_str())!=NULL) {
        RooRealVar* poi = (RooRealVar*) ws->var(poiName.c_str());
        RooRealVar* poi_parIni = ( parIni ? (RooRealVar*) parIni->find(poiName.c_str()) : NULL );
        if (p.second.count("Min")) { p.second.at("Min") = poi ? poi->getMin()   : -1.0; }
        if (p.second.count("Max")) { p.second.at("Max") = poi ? poi->getMax()   : -1.0; }
        if (p.second.count("Val")) { p.second.at("Val") = poi ? poi->getVal()   : -1.0; }
        if (p.second.count("Err")) { p.second.at("Err") = poi ? poi->getError() : -1.0; }
        if (p.second.count("parIni_Val")) { p.second.at("parIni_Val") = poi_parIni ? poi_parIni->getMin() : -1.0; }
        if (p.second.count("parIni_Err")) { p.second.at("parIni_Err") = poi_parIni ? poi_parIni->getMax() : -1.0; }
      }
      //
      else if (ws->function(poiName.c_str())!=NULL) {
        const std::string pdfName = ( "pdfMET_Tot" + tag );
        RooFitResult* fitResult = (RooFitResult*) ws->obj(Form("fitResult_%s", pdfName.c_str()));
        RooFormulaVar* poi = (RooFormulaVar*) ws->function(poiName.c_str());
        if (p.second.count("Min")) { p.second.at("Min") = -1.0; }
        if (p.second.count("Max")) { p.second.at("Max") = -1.0; }
        if (p.second.count("Val")) { p.second.at("Val") = poi ? poi->getVal() : -1.0; }
        if (p.second.count("Err")) { p.second.at("Err") = (poi && fitResult) ? poi->getPropagatedError(*fitResult) : -1.0; }
        if (p.second.count("parIni_Val")) { p.second.at("parIni_Val") = (poi && parIni) ? poi->getVal(*parIni) : -1.0; }
        if (p.second.count("parIni_Err")) { p.second.at("parIni_Err") = -1.0; }
      }
    }
    //
    // Clean up the memory
    delete ws;
    inputFile.Close();
    //
    // Fill the tree
    tree.Fill();
  } // loop on the files
  //
  // Create the output file
  gSystem->mkdir(outputDirPath.c_str(), kTRUE);
  TFile outputFile(outputFilePath.c_str(), "RECREATE");
  if (outputFile.IsOpen()==false || outputFile.IsZombie()==true) {
    std::cout << "[ERROR] The output file " << outputFilePath << " could not be created!" << std::endl;
  }
  outputFile.cd();
  //
  // Write the output tree
  tree.Write();
  //
  // Write the output file
  outputFile.Write();
  // Close the output file
  outputFile.Close();
  //
  // return
  return true;
}


#endif // #ifndef resultsEWQ2tree_C

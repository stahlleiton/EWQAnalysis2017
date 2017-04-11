// -*- C++ -*-
//
// Package:    Fitter
//
/*
 Description: TTree to RooDataSet converter.
 Implementation:
 This program create RooDataSets from TTrees.
 */
// Original Author:  Andre Stahl,
//         Created:  Mar 09 19:08 CET 2017
//
//
#ifndef tree2DataSet_C
#define tree2DataSet_C


#include "./EWQ/EWQForest2DataSet.C"
#include "./Utilities/initClasses.h"
#include "./Utilities/EVENTUTILS.h"


bool checkFileInfo ( const StringVectorMap& FileInfo );
bool checkAnalysis ( const std::string& Analysis , GlobalInfo& info );


bool tree2DataSet(RooWorkspaceMap& Workspaces, const StringVectorMap& FileInfo, const std::string& Analysis, const bool& UpdateDS=false)
{
  if (!checkFileInfo(FileInfo)) return false;
  // Create Global Container with input information
  GlobalInfo info;
  info.Int["triggerIndex"] = PA::HLT_PAL3Mu12;
  info.Flag["applyWeight"] = false;
  info.Flag["updateDS"] = UpdateDS;
  info.Par["Analysis"] = Analysis;
  info.Par["VarType"] = FileInfo.at("VarType")[0];
  // Check if Analysis type is supported
  if (!checkAnalysis(Analysis, info)) return false;
  // Make RooDatasets
  if (info.Par.at("Tree_Type") == "EWQForest") { if (!EWQForest2DataSet(Workspaces, FileInfo, info)) return false; }
  return true;
};


bool checkAnalysis(const std::string& Analysis, GlobalInfo& info)
{
  if (Analysis=="WToMuNu") { 
    std::cout << "[INFO] Proceed to make datasets for W to Muon+Neutrino Analysis using EWQ Forest" << std::endl;
    info.Par["Tree_Type"] = "EWQForest";
    return true; 
  }
  std::cout << "[ERROR] The input Analysis: " << Analysis << " is not supported by this fitter!" << std::endl;
  return false;
};


bool checkFileInfo(const StringVectorMap& FileInfo)
{
  if ( FileInfo.at("OutputFileDir").size()==0  ) { std::cout << "[ERROR] OutputFileDir is empty!" << std::endl; return false;  }
  if ( FileInfo.at("InputFileNames").size()==0 ) { std::cout << "[ERROR] InputFileNames is empty!" << std::endl; return false; }
  if ( FileInfo.at("DSNames").size()==0        ) { std::cout << "[ERROR] DSNames is empty!" << std::endl; return false;        }
  return true;
};


#endif // #ifndef tree2DataSet_C

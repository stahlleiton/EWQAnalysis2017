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
#include "../../Utilities/EVENTUTILS.h"
#include "../../Utilities/tnp_weight.h"


bool checkFileInfo    ( const StringVectorMap_t& FileInfo );
bool checkAnalysis    ( const std::string& Analysis   , GlobalInfo& info );
bool createPADataset  ( RooWorkspaceMap_t& Workspaces , const std::string& sampleTag );
bool invertEtaAndFill ( RooDataSet& dsPA , const RooDataSet& dsPbp );


bool tree2DataSet(RooWorkspaceMap_t& Workspaces, const StringVectorMap_t& FileInfo, const std::string& VarType, const std::string& Analysis, const bool& UpdateDS=false)
{
  if (!checkFileInfo(FileInfo)) return false;
  // Create Global Container with input information
  GlobalInfo info;
  info.Int ["triggerIndex"] = PA::HLT_PAL3Mu12;
  info.Flag["applyCorr"]    = true;
  info.Flag["updateDS"]     = UpdateDS;
  info.Par ["Analysis"]     = Analysis;
  info.Par ["VarType"]      = VarType;
  // Check if Analysis type is supported
  if (!checkAnalysis(Analysis, info)) return false;
  // Make RooDatasets
  if (info.Par.at("Tree_Type") == "EWQForest") { if (!EWQForest2DataSet(Workspaces, FileInfo, info)) return false; }
  return true;
};


bool checkAnalysis(const std::string& Analysis, GlobalInfo& info)
{
  if (Analysis=="WToMuNu") { 
    std::cout << "[INFO] Proceed to make datasets for W to Muon + Neutrino Analysis using EWQ Forest" << std::endl;
    info.Par["Tree_Type"] = "EWQForest";
    return true; 
  }
  std::cout << "[ERROR] The input Analysis: " << Analysis << " is not supported by this fitter!" << std::endl;
  return false;
};


bool checkFileInfo(const StringVectorMap_t& FileInfo)
{
  if ( FileInfo.at("OutputFileDir").size()==0  ) { std::cout << "[ERROR] OutputFileDir is empty!" << std::endl; return false;  }
  if ( FileInfo.at("InputFileNames").size()==0 ) { std::cout << "[ERROR] InputFileNames is empty!" << std::endl; return false; }
  if ( FileInfo.at("DSNames").size()==0        ) { std::cout << "[ERROR] DSNames is empty!" << std::endl; return false;        }
  return true;
};


#endif // #ifndef tree2DataSet_C

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


bool checkFileInfo    ( const StringVectorMap& FileInfo );
bool checkAnalysis    ( const std::string& Analysis , GlobalInfo& info );
bool createPADataset  ( RooWorkspaceMap& Workspaces , const std::string& sampleTag );
bool invertEtaAndFill ( RooDataSet& dsPA , const RooDataSet& dsPbp );


bool tree2DataSet(RooWorkspaceMap& Workspaces, const StringVectorMap& FileInfo, const std::string& VarType, const std::string& Analysis, const bool& UpdateDS=false)
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


bool checkFileInfo(const StringVectorMap& FileInfo)
{
  if ( FileInfo.at("OutputFileDir").size()==0  ) { std::cout << "[ERROR] OutputFileDir is empty!" << std::endl; return false;  }
  if ( FileInfo.at("InputFileNames").size()==0 ) { std::cout << "[ERROR] InputFileNames is empty!" << std::endl; return false; }
  if ( FileInfo.at("DSNames").size()==0        ) { std::cout << "[ERROR] DSNames is empty!" << std::endl; return false;        }
  return true;
};


bool createPADataset(RooWorkspaceMap& Workspaces, const std::string& sampleTag)
{
  //
  const std::string sample_pPb = sampleTag + "_pPb";
  const std::string sample_Pbp = sampleTag + "_Pbp";
  const std::string sample_PA  = sampleTag + "_PA";
  //
  // Check input datasets
  if (Workspaces.count(sample_pPb)==0) { std::cout << "[ERROR] RooWorkspace for sample " << sample_pPb << " does not exist!" << std::endl; return false; }
  if (Workspaces.count(sample_Pbp)==0) { std::cout << "[ERROR] RooWorkspace for sample " << sample_Pbp << " does not exist!" << std::endl; return false; }
  //
  // Extract the RooDatasets
  auto dPl_pPb = (RooDataSet*)Workspaces.at(sample_pPb).data(("dPl_"+sample_pPb).c_str());
  if (dPl_pPb==NULL) { std::cout << "[ERROR] RooDataSet " << ("dPl_"+sample_pPb) << " does not exist!" << std::endl; return false; }
  auto dMi_pPb = (RooDataSet*)Workspaces.at(sample_pPb).data(("dMi_"+sample_pPb).c_str());
  if (dMi_pPb==NULL) { std::cout << "[ERROR] RooDataSet " << ("dMi_"+sample_pPb) << " does not exist!" << std::endl; return false; }
  //
  auto dPl_Pbp = (RooDataSet*)Workspaces.at(sample_Pbp).data(("dPl_"+sample_Pbp).c_str());
  if (dPl_Pbp==NULL) { std::cout << "[ERROR] RooDataSet " << ("dPl_"+sample_Pbp) << " does not exist!" << std::endl; return false; }
  auto dMi_Pbp = (RooDataSet*)Workspaces.at(sample_Pbp).data(("dMi_"+sample_Pbp).c_str());
  if (dMi_Pbp==NULL) { std::cout << "[ERROR] RooDataSet " << ("dMi_"+sample_Pbp) << " does not exist!" << std::endl; return false; }
  //
  std::cout << "[INFO] Creating combined PA RooDataSet for " << sampleTag << std::endl;
  //
  // Copy the pPb datasets
  auto dPl_PA = (RooDataSet*)dPl_pPb->Clone(("dPl_"+sample_PA).c_str());
  if (dPl_PA==NULL || dPl_PA->sumEntries()==0) { std::cout << "[ERROR] RooDataSet " << ("dPl_"+sample_PA) << " was not created!" << std::endl; return false; }
  auto dMi_PA = (RooDataSet*)dMi_pPb->Clone(("dMi_"+sample_PA).c_str());
  if (dMi_PA==NULL || dMi_PA->sumEntries()==0) { std::cout << "[ERROR] RooDataSet " << ("dMi_"+sample_PA) << " was not created!" << std::endl; return false; }
  //
  // Invert the eta of Pbp dataset and fill the PA dataset
  if (!invertEtaAndFill(*dPl_PA, *dPl_Pbp)) { return false; }
  if (!invertEtaAndFill(*dMi_PA, *dMi_Pbp)) { return false; }
  //
  // Check the consistency of the combined dataset
  if (dPl_PA->numEntries()!=(dPl_pPb->numEntries()+dPl_Pbp->numEntries())) { std::cout << "[ERROR] Number of entries for the combined Pl " << sampleTag << " dataset is inconsistent!" << std::endl; return false; }
  if (dMi_PA->numEntries()!=(dMi_pPb->numEntries()+dMi_Pbp->numEntries())) { std::cout << "[ERROR] Number of entries for the combined Mi " << sampleTag << " dataset is inconsistent!" << std::endl; return false; }
  if (std::abs(dPl_PA->sumEntries()-(dPl_pPb->sumEntries()+dPl_Pbp->sumEntries()))>0.05*(dPl_pPb->sumEntries()+dPl_Pbp->sumEntries())) { 
    std::cout << "[ERROR] Number of weighted entries for the combined Pl " << sampleTag << " dataset is inconsistent!" << std::endl; return false;
  }
  if (std::abs(dMi_PA->sumEntries()-(dMi_pPb->sumEntries()+dMi_Pbp->sumEntries()))>0.05*(dMi_pPb->sumEntries()+dMi_Pbp->sumEntries())) {
    std::cout << "[ERROR] Number of weighted entries for the combined Mi " << sampleTag << " dataset is inconsistent!" << std::endl; return false;
  }
  //
  // Import the new PA dataset
  Workspaces[sample_PA].import(*dPl_PA);
  Workspaces.at(sample_PA).import(*dMi_PA);
  Workspaces.at(sample_PA).import(*((TObjString*)Workspaces.at(sample_pPb).obj("METType")), "METType");
  //
  // Clean up the memory
  if (dPl_PA) delete dPl_PA;
  if (dMi_PA) delete dMi_PA;
  //
  // return
  return true;
};


bool invertEtaAndFill(RooDataSet& dsPA, const RooDataSet& dsPbp)
{
  for(int i = 0; i < dsPbp.numEntries(); i++){
    auto set = dsPbp.get(i);
    auto eta = (RooRealVar*)set->find("Muon_Eta");
    auto pt  = (RooRealVar*)set->find("Muon_Pt");
    if (eta==NULL) { std::cout << "[ERROR] RooRealVar Muon_Eta was not found!" << std::endl; return false; }
    // Correct the weight
    double weight = dsPbp.weight();
    if (std::string(dsPbp.GetName()).find("MC_")!=std::string::npos) {
      // Correct the TnP scale factor
      const double sf_MuID = tnp_weight_muid_ppb( pt->getVal() , eta->getVal() , 0 );
      const double sf_Trig = tnp_weight_trg_ppb (                eta->getVal() , 0 );
      const double sf_Iso  = tnp_weight_iso_ppb ( pt->getVal() , eta->getVal() , 0 );
      const double sf_MuID_Inv = tnp_weight_muid_ppb( pt->getVal() , -1.0*eta->getVal() , 0 );
      const double sf_Trig_Inv = tnp_weight_trg_ppb (                -1.0*eta->getVal() , 0 );
      const double sf_Iso_Inv  = tnp_weight_iso_ppb ( pt->getVal() , -1.0*eta->getVal() , 0 );
      weight = weight * ( (sf_MuID_Inv * sf_Trig_Inv * sf_Iso_Inv) / (sf_MuID * sf_Trig * sf_Iso) );
    }
    // Invert the pseudo-rapidity
    eta->setVal(-1.0*eta->getVal());
    // Add the new event
    dsPA.add(*set, weight);
  }
  return true;
};


#endif // #ifndef tree2DataSet_C

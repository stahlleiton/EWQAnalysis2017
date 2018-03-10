// -*- C++ -*-
//
// Package:    Fitter
//
/*
 Description: EWQ Forest to RooDataSet converter.
 Implementation:
 This program create RooDataSets from EWQ Forest.
 */
// Original Author:  Andre Stahl,
//         Created:  Mar 15 19:08 CET 2017
//
//

#ifndef EWQForest2DataSet_C
#define EWQForest2DataSet_C


#include "../../../Utilities/HiMETTree.h"
#include "../../../Utilities/HiMuonTree.h"
#include "../../../Utilities/HiEvtTree.h"
#include "../../../Utilities/tnp_weight.h"
#include "../../../Utilities/HFweight.h"
#include "../../../Corrections/MET_Recoil/RecoilCorrector.C"
#include "../Utilities/initClasses.h"


bool checkEWQDS          ( RooDataSet* DS                , const std::string& DSName         , const std::string& Analysis );
bool EWQForest_WToMuNu   ( RooWorkspaceMap_t& Workspaces , const StringVectorMap_t& FileInfo , const GlobalInfo&  info     );


bool EWQForest2DataSet(RooWorkspaceMap_t& Workspaces, const StringVectorMap_t& FileInfo, const GlobalInfo& info)
{
  std::string Analysis = info.Par.at("Analysis");
  if (Analysis=="WToMuNu") { if (!EWQForest_WToMuNu(Workspaces, FileInfo, info)) return false; }
  return true;
};

bool EWQForest_WToMuNu(RooWorkspaceMap_t& Workspaces, const StringVectorMap_t& FileInfo, const GlobalInfo& info)
{
  StringVector OutputFileNames;
  const std::string  chaDir = "Muon";
  const std::string  metTAG = "MET" + info.Par.at("VarType") + "_";
  const StringVector OutputFileDir  = FileInfo.at("OutputFileDir");
  const StringVector InputFileNames = FileInfo.at("InputFileNames");
  const StringVector DSNames        = FileInfo.at("DSNames");
  const bool isData = (DSNames[0].find("DATA")!=std::string::npos);
  const bool isMC   = (DSNames[0].find("MC")!=std::string::npos);
  for (auto& tag : DSNames) { 
    std::string o = (OutputFileDir[0] + chaDir + "/") + "DATASET_" + metTAG + tag + ".root"; 
    if (gSystem->AccessPathName(o.c_str())) { makeDir(OutputFileDir[1] + chaDir + "/"); o = (OutputFileDir[1] + chaDir + "/") + "DATASET_" + metTAG + tag + ".root"; }
    OutputFileNames.push_back(o);
  }
  // Extract Input Information
  const std::string  TYPE = info.Par.at("Analysis");
  const int triggerIndex  = info.Int.at("triggerIndex");
  // Create RooDataSets
  std::vector< std::unique_ptr<RooDataSet> > dataPl, dataMi, mcPl, mcMi, metPl, metMi;
  bool createDS = info.Flag.at("updateDS");
  // Check if RooDataSets exist and are not corrupt
  for (uint i=0; i<OutputFileNames.size(); i++) {
    if ( !gSystem->AccessPathName(OutputFileNames[i].c_str()) ) {
      std::cout << "[INFO] Loading RooDataSets from " << OutputFileNames[i] << std::endl;
      std::unique_ptr<TFile> DBFile = std::unique_ptr<TFile>(TFile::Open(OutputFileNames[i].c_str(),"READ"));
      if (!DBFile) { std::cout << "[ERROR] File: " << OutputFileNames[i] << " is corrupted!" << std::endl; return false; }
      dataPl.push_back( std::unique_ptr<RooDataSet>((RooDataSet*)DBFile->Get(Form("dPl_RAW_%s", DSNames[i].c_str()))) );
      dataMi.push_back( std::unique_ptr<RooDataSet>((RooDataSet*)DBFile->Get(Form("dMi_RAW_%s", DSNames[i].c_str()))) );
      if (isMC) { mcPl.push_back( std::unique_ptr<RooDataSet>((RooDataSet*)DBFile->Get(Form("mcPl_RAW_%s", DSNames[i].c_str()))) ); }
      if (isMC) { mcMi.push_back( std::unique_ptr<RooDataSet>((RooDataSet*)DBFile->Get(Form("mcMi_RAW_%s", DSNames[i].c_str()))) ); }
      metPl.push_back( std::unique_ptr<RooDataSet>((RooDataSet*)DBFile->Get(Form("metPl_RAW_%s", DSNames[i].c_str()))) );
      metMi.push_back( std::unique_ptr<RooDataSet>((RooDataSet*)DBFile->Get(Form("metMi_RAW_%s", DSNames[i].c_str()))) );
      if (checkEWQDS(dataPl[i].get(), DSNames[i], TYPE)==false) { createDS = true; }
      if (checkEWQDS(dataMi[i].get(), DSNames[i], TYPE)==false) { createDS = true; }
      DBFile->Close();
    }
    else { createDS = true; break; }
  }
  if (createDS) {
    ///// Input Forest
    //
    std::unique_ptr<HiMuonTree> muonTree = std::unique_ptr<HiMuonTree>(new HiMuonTree());
    if (!muonTree->GetTree(InputFileNames)) return false;
    Long64_t nentries = muonTree->GetEntries();
    std::unique_ptr<HiMETTree> metTree = std::unique_ptr<HiMETTree>(new HiMETTree());
    if (!metTree->GetTree(InputFileNames, "metAna")) return false;
    if (metTree->GetEntries() != nentries) { std::cout << "[ERROR] Inconsistent number of entries in metAna!" << std::endl; return false; }
    std::unique_ptr<HiMETTree> metNoHFTree = std::unique_ptr<HiMETTree>(new HiMETTree());
    if (!metNoHFTree->GetTree(InputFileNames, "metAnaNoHF")) return false;
    if (metNoHFTree->GetEntries() != nentries) { std::cout << "[ERROR] Inconsistent number of entries in metNoHFAna!" << std::endl; return false; }
    std::unique_ptr<HiEvtTree> evtTree = std::unique_ptr<HiEvtTree>(new HiEvtTree());
    if (!evtTree->GetTree(InputFileNames)) return false;
    if (evtTree->GetEntries() != nentries) { std::cout << "[ERROR] Inconsistent number of entries in evtTree!" << std::endl; return false; }
    ///// MC Generator Information
    RooRealVar mcNGen   = RooRealVar ( "NGen"      , "Number of Generated Events" ,   -1.0 ,   1000000000.0 , ""      );
    RooRealVar mcType   = RooRealVar ( "MC_Type"   , "MC Type Number"             , -100.0 ,          100.0 , ""      );
    RooRealVar bosonPt  = RooRealVar ( "Boson_Pt"  , "Boson p_{T}"                ,   -1.0 ,       100000.0 , "GeV/c" );
    RooRealVar bosonPhi = RooRealVar ( "Boson_Phi" , "Boson #phi"                 ,   -9.0 ,            9.0 , ""      );
    RooRealVar refPt    = RooRealVar ( "Ref_Pt"    , "Reference p_{T}"            ,   -1.0 ,       100000.0 , "GeV/c" );
    RooRealVar refPhi   = RooRealVar ( "Ref_Phi"   , "Reference #phi"             ,   -9.0 ,            9.0 , ""      );
    RooRealVar metPhi   = RooRealVar ( "MET_Phi"   , "#slash{E}_{T} #phi"         ,   -9.0 ,            9.0 , ""      );
    RooRealVar muPhi    = RooRealVar ( "Muon_Phi"  , "#mu #phi"                   ,   -9.0 ,            9.0 , ""      );
    RooRealVar hiHF     = RooRealVar ( "hiHF"      , "Total E_{HF}"               ,   -1.0 ,       100000.0 , "GeV/c" );
    RooRealVar nTracks  = RooRealVar ( "nTracks"   , "Number of Tracks"           ,   -1.0 ,   1000000000.0 , ""      );
    //
    RooArgSet  mcCols = RooArgSet(mcNGen, mcType, bosonPt, bosonPhi, refPt, refPhi);
    mcCols.add(metPhi); mcCols.add(muPhi); mcCols.add(hiHF); mcCols.add(nTracks);
    ///// MET Information
    RooRealVar metMag_RAW       = RooRealVar ( "MET_Mag_PF_RAW"        , "|#slash{E}_{T}|"     ,   -1.0 , 100000.0 , "GeV/c" );
    RooRealVar metPhi_RAW       = RooRealVar ( "MET_Phi_PF_RAW"        , "#slash{E}_{T} #phi"  ,   -9.0 ,      9.0 , ""      );
    RooRealVar metMag_TYPE1     = RooRealVar ( "MET_Mag_PF_Type1"      , "|#slash{E}_{T}|"     ,   -1.0 , 100000.0 , "GeV/c" );
    RooRealVar metPhi_TYPE1     = RooRealVar ( "MET_Phi_PF_Type1"      , "#slash{E}_{T} #phi"  ,   -9.0 ,      9.0 , ""      );
    RooRealVar metMag_NoHFRAW   = RooRealVar ( "MET_Mag_PF_NoHF_RAW"   , "|#slash{E}_{T}|"     ,   -1.0 , 100000.0 , "GeV/c" );
    RooRealVar metPhi_NoHFRAW   = RooRealVar ( "MET_Phi_PF_NoHF_RAW"   , "#slash{E}_{T} #phi"  ,   -9.0 ,      9.0 , ""      );
    RooRealVar metMag_NoHFTYPE1 = RooRealVar ( "MET_Mag_PF_NoHF_Type1" , "|#slash{E}_{T}|"     ,   -1.0 , 100000.0 , "GeV/c" );
    RooRealVar metPhi_NoHFTYPE1 = RooRealVar ( "MET_Phi_PF_NoHF_Type1" , "#slash{E}_{T} #phi"  ,   -9.0 ,      9.0 , ""      );
    //
    RooArgSet  metCols = RooArgSet(metMag_RAW, metPhi_RAW, metMag_TYPE1, metPhi_TYPE1);
    metCols.add(metMag_NoHFRAW);   metCols.add(metPhi_NoHFRAW);
    metCols.add(metMag_NoHFTYPE1); metCols.add(metPhi_NoHFTYPE1);
    ///// RooDataSet Variables
    RooRealVar   met    = RooRealVar ( "MET"        , "|#slash{E}_{T}|"   ,  -1.0 , 100000.0 ,  "GeV/c"     );
    RooRealVar   muPt   = RooRealVar ( "Muon_Pt"    , "#mu p_{T}"         ,  -1.0 , 100000.0 ,  "GeV/c"     );
    RooRealVar   muEta  = RooRealVar ( "Muon_Eta"   , "#mu #eta"          , -10.0 ,     10.0 ,  ""          );
    RooRealVar   muIso  = RooRealVar ( "Muon_Iso"   , "#mu Isolation"     ,  -1.0 , 100000.0 ,  ""          );
    RooRealVar   muMT   = RooRealVar ( "Muon_MT"    , "W Transverse Mass" ,  -1.0 , 100000.0 ,  "GeV/c^{2}" );
    RooRealVar   cent   = RooRealVar ( "Centrality" , "Centrality"        ,  -1.0 , 100000.0 ,  ""          );
    RooCategory  type   = RooCategory( "Event_Type" , "Event Type");
    type.defineType("Other", -1); type.defineType("DYToMuMu", 1); type.defineType("ZToMuMu", 2);
    RooArgSet cols = RooArgSet(met, muPt, muEta, muIso, muMT, cent);
    cols.add(type);
    ///// Initiliaze RooDataSets
    dataPl.clear(); dataMi.clear(); mcPl.clear(); mcMi.clear(); metPl.clear(); metMi.clear();
    for (uint i=0; i<DSNames.size(); i++) {
      std::cout << "[INFO] Creating " << "RooDataSet for " << DSNames[i] << std::endl;
      dataPl.push_back( std::unique_ptr<RooDataSet>(new RooDataSet(Form("dPl_RAW_%s", DSNames[i].c_str()), "dPl", cols)) );
      dataMi.push_back( std::unique_ptr<RooDataSet>(new RooDataSet(Form("dMi_RAW_%s", DSNames[i].c_str()), "dMi", cols)) );
      if (isMC) { mcPl.push_back( std::unique_ptr<RooDataSet>(new RooDataSet(Form("mcPl_RAW_%s", DSNames[i].c_str()), "mcPl", mcCols )) ); }
      if (isMC) { mcMi.push_back( std::unique_ptr<RooDataSet>(new RooDataSet(Form("mcMi_RAW_%s", DSNames[i].c_str()), "mcMi", mcCols )) ); }
      metPl.push_back( std::unique_ptr<RooDataSet>(new RooDataSet(Form("metPl_RAW_%s", DSNames[i].c_str()), "metPl", metCols)) );
      metMi.push_back( std::unique_ptr<RooDataSet>(new RooDataSet(Form("metMi_RAW_%s", DSNames[i].c_str()), "metMi", metCols)) );
    }
    ///// Iterate over the Input Forest
    int treeIdx = -1;
    std::cout << "[INFO] Starting to process " << nentries << " nentries" << std::endl;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
      //
      // Get the entry in the trees
      if (muonTree->GetEntry(jentry)<0) { std::cout << "[ERROR] Muon Tree invalid entry!"  << std::endl; return false; }
      if (metTree ->GetEntry(jentry)<0) { std::cout << "[ERROR] MET Tree invalid entry!"   << std::endl; return false; }
      if (evtTree ->GetEntry(jentry)<0) { std::cout << "[ERROR] Event Tree invalid entry!" << std::endl; return false; }
      if (metNoHFTree->GetEntry(jentry)<0) { std::cout << "[ERROR] MET NoHF Tree invalid entry!" << std::endl; return false; }
      // 
      // Check that the different tree agrees well
      if (muonTree->Event_Run()!=metTree->Event_Run()       ) { std::cout << "[ERROR] MET Run does not agree!"     << std::endl; return false; }
      if (muonTree->Event_Number()!=metTree->Event_Number() ) { std::cout << "[ERROR] MET Event does not agree!"   << std::endl; return false; }
      if (muonTree->Event_Run()!=evtTree->run()             ) { std::cout << "[ERROR] HiEVT Run does not agree!"   << std::endl; return false; }
      if (muonTree->Event_Number()!=evtTree->evt()          ) { std::cout << "[ERROR] HiEVT Event does not agree!" << std::endl; return false; }
      //
      if (muonTree->Chain()->GetTreeNumber()!=treeIdx) {
        treeIdx = muonTree->Chain()->GetTreeNumber();
        std::cout << "[INFO] Processing Root File: " << InputFileNames[treeIdx] << std::endl;
        // Get the MC Info
        if (isMC) {
          mcNGen.setVal ( muonTree->GetTreeEntries() );
          int mcID = PA::getMCTypeID(FileInfo.at("TreeTags")[treeIdx]);
          if (mcID==0) { std::cout << "[ERROR] MC Tag " << FileInfo.at("TreeTags")[treeIdx] << " is invalid!" << std::endl; return false; }
          mcType.setVal ( mcID );
        }
      }
      //
      loadBar(jentry, nentries);
      //
      // Determine the collision system of the sample
      std::string evtCol = "";
      if (isData) {
        if (muonTree->Event_Run() >= 285410 && muonTree->Event_Run() <= 285951) evtCol = "Pbp"; // for Pbp
        if (muonTree->Event_Run() >= 285952 && muonTree->Event_Run() <= 286504) evtCol = "pPb"; // for pPb
      }
      else {
        if (DSNames[0].find("Pbp")!=std::string::npos) evtCol = "Pbp"; // for Pbp
        if (DSNames[0].find("pPb")!=std::string::npos) evtCol = "pPb"; // for pPb
      }
      if (evtCol=="") { std::cout << "[ERROR] Could not determine the collision system in the sample" << std::endl; return false; }
      //
      // Event Based Information
      //
      // Apply Event Filters
      if (PA::passEventFilter(metTree)==false) continue; // PA Event Selection
      //
      // Check Trigger Fired
      if (muonTree->Event_Trig_Fired()[triggerIndex]==false) continue; // Trigger Event Fired
      //
      // Determine the MET
      TVector2 MET = TVector2();
      if      (info.Par.at("VarType") == "PF_RAW"       ) { MET = metTree->PF_MET_NoShift_Mom();        }
      else if (info.Par.at("VarType") == "PF_Type1"     ) { MET = metTree->Type1_MET_NoShift_Mom();     }
      else if (info.Par.at("VarType") == "PF_NoHF_RAW"  ) { MET = metNoHFTree->PF_MET_NoShift_Mom();    }
      else if (info.Par.at("VarType") == "PF_NoHF_Type1") { MET = metNoHFTree->Type1_MET_NoShift_Mom(); }
      else { std::cout << "[ERROR] MET Type " << info.Par.at("VarType") << " is not valid!" << std::endl; return false; }
      //
      // Muon Based Information
      //
      // Find the Leading Pt Muon
      float maxPt = -99.; int leadMuPFIdx = -1;
      for (uint imu = 0; imu < muonTree->PF_Muon_Mom().size(); imu++) {
        if (muonTree->PF_Muon_Reco_Idx()[imu]==-1) { std::cout << "[ERROR] Reco idx is -1" << std::endl; return false; }
        // Only consider Muons Matched to GEN in MC
        if (isMC && muonTree->PF_Muon_Gen_Idx()[imu]<0) continue;
        if (isMC && PA::checkGenMuon(muonTree->PF_Muon_Gen_Idx()[imu], DSNames[0], muonTree) == false) continue;
        // Only consider Tight Muons within the Pseudo-Rapidity acceptance of CMS
        if (PA::isTightMuon(imu, muonTree)==false) continue;
        // Determine the highest pT muon
        if (maxPt < muonTree->PF_Muon_Mom()[imu].Pt()) { maxPt = muonTree->PF_Muon_Mom()[imu].Pt(); leadMuPFIdx = imu; }
      }
      if (leadMuPFIdx<=-1) continue;  // Only consider events with at least one good muon
      //
      // Apply PT Cut
      TLorentzVector muP4 = muonTree->PF_Muon_Mom()[leadMuPFIdx];
      if (muP4.Pt() < 25.) continue;  // Consider only leading muons with pT larger or equal than 25 GeV
      //
      // Apply Trigger Matching
      if (PA::isTriggerMatched(triggerIndex, leadMuPFIdx, muonTree)==false) continue;  // Only consider Muons Matched to Trigger
      //
      //
      // Classify the events
      //
      std::string eventType = "Other";
      //
      // Event Type: DrellYan->MuMu
      if (eventType=="Other") {
        if (PA::passDrellYanVeto(muonTree) == false) {
          if (PA::hasZBoson(muonTree) == true) { eventType = "ZToMuMu"; }  // Found a Z->MuMu candidate
          else { eventType = "DYToMuMu"; } // Found a Drell-Yan candidate outside of Z mass region
        }
      }
      //
      // Isolation and Charge of Leading Muon
      const float leadMuIso = muonTree->PF_Muon_IsoPFR03NoPUCorr()[leadMuPFIdx];
      const int   leadMuChg = int(muonTree->PF_Muon_Charge()[leadMuPFIdx]);
      //
      // Recompute the Transver Mass based on the chosen MET
      const double muMTVal = PA::getWTransverseMass(muP4.Pt(), muP4.Phi(), MET.Mod(), MET.Phi());
      //
      //// Set the variables
      met.setVal    ( MET.Mod()  );
      muPt.setVal   ( muP4.Pt()  );
      muEta.setVal  ( muP4.Eta() );
      muIso.setVal  ( leadMuIso  );
      muMT.setVal   ( muMTVal    );
      cent.setVal   ( 0.0        );
      type.setLabel ( eventType.c_str() );
      //
      //// Set the MET information
      metMag_RAW.setVal       ( metTree->PF_MET_NoShift_Mom().Mod()        );
      metPhi_RAW.setVal       ( metTree->PF_MET_NoShift_Mom().Phi()        );
      metMag_TYPE1.setVal     ( metTree->Type1_MET_NoShift_Mom().Mod()     );
      metPhi_TYPE1.setVal     ( metTree->Type1_MET_NoShift_Mom().Phi()     );
      metMag_NoHFRAW.setVal   ( metNoHFTree->PF_MET_NoShift_Mom().Mod()    );
      metPhi_NoHFRAW.setVal   ( metNoHFTree->PF_MET_NoShift_Mom().Phi()    );
      metMag_NoHFTYPE1.setVal ( metNoHFTree->Type1_MET_NoShift_Mom().Mod() );
      metPhi_NoHFTYPE1.setVal ( metNoHFTree->Type1_MET_NoShift_Mom().Phi() );
      //
      // Get the Information needed for the Recoil Corrections
      if (isMC) {
        TVector2 ref_Pt , boson_Pt;
        if (!RecoilCorrector::getPtFromTree(ref_Pt, boson_Pt, leadMuPFIdx, muonTree, DSNames[0])) { return false; }
        bosonPt.setVal  ( boson_Pt.Mod() );
        bosonPhi.setVal ( boson_Pt.Phi() );
        refPt.setVal    ( ref_Pt.Mod()   );
        refPhi.setVal   ( ref_Pt.Phi()   );
        metPhi.setVal   ( MET.Phi()      );
        muPhi.setVal    ( muP4.Phi()     );
        hiHF.setVal     ( evtTree->hiHF());
        nTracks.setVal  ( evtTree->hiNtracks());
      }
      //
      //// Fill the RooDataSets
      for (uint i=0; i<DSNames.size(); i++) {
        if (DSNames[i].find(evtCol)!=std::string::npos) {
          if (leadMuChg > 0) { dataPl[i]->addFast(cols); }
          if (leadMuChg < 0) { dataMi[i]->addFast(cols); }
          if (isMC && (leadMuChg > 0)) { mcPl[i]->addFast(mcCols); }
          if (isMC && (leadMuChg < 0)) { mcMi[i]->addFast(mcCols); }
          if (leadMuChg > 0) { metPl[i]->addFast(metCols); }
          if (leadMuChg < 0) { metMi[i]->addFast(metCols); }
        }
      }
    }
    //// Save the RooDataSets
    for (uint i=0; i<DSNames.size(); i++) {
      std::unique_ptr<TFile> DBFile = std::unique_ptr<TFile>(TFile::Open(OutputFileNames[i].c_str(),"RECREATE"));
      DBFile->cd();
      dataPl[i]->Write(Form("dPl_RAW_%s", DSNames[i].c_str()));
      dataMi[i]->Write(Form("dMi_RAW_%s", DSNames[i].c_str()));
      if (isMC) { mcPl[i]->Write(Form("mcPl_RAW_%s", DSNames[i].c_str())); }
      if (isMC) { mcMi[i]->Write(Form("mcMi_RAW_%s", DSNames[i].c_str())); }
      metPl[i]->Write(Form("metPl_RAW_%s", DSNames[i].c_str()));
      metMi[i]->Write(Form("metMi_RAW_%s", DSNames[i].c_str()));
      DBFile->Write(); DBFile->Close();
    }
  }
  // Import datasets to the workspaces
  for (uint i=0; i<DSNames.size(); i++) {
    if(!dataPl[i]) { std::cout << "[ERROR] " << DSNames[i] << " plus dataset was not found" << std::endl;  return false; }
    if(!dataMi[i]) { std::cout << "[ERROR] " << DSNames[i] << " minus dataset was not found" << std::endl; return false; }
    if(dataPl[i]->numEntries()==0) { std::cout << "[WARNING] " << DSNames[i] << " plus dataset is empty!" << std::endl;  return false; }
    if(dataMi[i]->numEntries()==0) { std::cout << "[WARNING] " << DSNames[i] << " minus dataset is empty!" << std::endl; return false; }
    Workspaces[DSNames[i]].import(*dataPl[i]);
    Workspaces[DSNames[i]].import(*dataMi[i]);
    if (isMC) {
      if(mcPl[i] && mcPl[i]->numEntries()>0) Workspaces[DSNames[i]].import(*mcPl[i]);
      if(mcMi[i] && mcMi[i]->numEntries()>0) Workspaces[DSNames[i]].import(*mcMi[i]);
    }
    if (metPl[i]) Workspaces[DSNames[i]].import(*metPl[i]);
    if (metMi[i]) Workspaces[DSNames[i]].import(*metMi[i]);
    RooStringVar tmp; tmp.setVal(info.Par.at("VarType").c_str()); tmp.SetTitle("METType"); Workspaces[DSNames[i]].import(*((TObject*)&tmp), tmp.GetTitle());
  }
  dataPl.clear(); dataMi.clear(); mcPl.clear(); mcMi.clear(); metPl.clear(); metMi.clear();
  return true;
};


bool checkEWQDS(RooDataSet* DS, const string& DSName, const std::string& Analysis)
{
  if (DS->numEntries()==0 || DS->sumEntries()==0) { std::cout << "[WARNING] Original dataset: " << DS->GetName() << " is empty, will remake it!" << std::endl; return false; }
  const RooArgSet* row = DS->get();
  if (Analysis=="WToMuNu") {
    if (
        ( row->find("MET")        !=0 ) &&
        ( row->find("Muon_Eta")   !=0 ) &&
        ( row->find("Muon_Pt")    !=0 ) &&
        ( row->find("Muon_Iso")   !=0 ) &&
        ( row->find("Muon_MT")    !=0 ) &&
        ( row->find("Centrality") !=0 ) &&
        ( row->find("Event_Type") !=0 )
        ) 
      { return true; }
    else 
      { std::cout << "[WARNING] Original dataset: " << DS->GetName() << " is corrupted, will remake it!" << std::endl; }
  }
  return false;
};


bool setMETonDS(RooWorkspace& ws, const std::string& dDSName, const std::string& mcDSName, const std::string& metDSName, const std::string& metType)
{
  // Check if Data
  const bool isMC = (dDSName.find("MC_")!=std::string::npos);
  // Get the RooDatasets
  auto  dDS = (RooDataSet*)ws.data(dDSName.c_str());
  if ( dDS==NULL) { std::cout << "[ERROR] setMETonDS: Dataset " << dDSName << " was not found in the workspace!" << std::endl; return false; }
  auto metDS = (RooDataSet*)ws.data(metDSName.c_str());
  if (metDS==NULL) { std::cout << "[ERROR] setMETonDS: Dataset " << metDSName << " was not found in the workspace!" << std::endl; return false; }
  RooDataSet* mcDS = NULL; if (isMC) { mcDS = (RooDataSet*)ws.data(mcDSName.c_str()); }
  if (isMC && mcDS==NULL) { std::cout << "[ERROR] setMETonDS: Dataset " << mcDSName << " was not found in the workspace!" << std::endl; return false; }
  // Get the sample name
  std::string sample = dDS->GetName(); sample = sample.substr(sample.find((isMC ? "MC_" : "DATA_")));
  std::string chg    = dDS->GetName(); chg    = chg.substr(1, chg.find("_")-1);
  std::string col    = dDS->GetName(); col    = col.substr(col.find_last_of("_")+1);
  if ( (col!="pPb") && (col!="Pbp") ) { std::cout << "[ERROR] Collision system " << col << " is invalid!" << std::endl; return false; }
  if ( (chg!="Pl" ) && (chg!="Mi" ) ) { std::cout << "[ERROR] Charge "           << chg << " is invalid!" << std::endl; return false; }
  // Initialize the RooArgSets
  RooArgSet cols    = *dDS->get();
  RooArgSet metCols = *metDS->get();
  RooArgSet mcCols = (isMC ? *mcDS->get() : RooArgSet());
  // Initialize the new RooDataSets
  auto dWDS  = std::unique_ptr<RooDataSet>(new RooDataSet( Form("d%s_SET_%s" , chg.c_str(), sample.c_str()) , Form("d%s" , chg.c_str()) , cols) );
  std::unique_ptr<RooDataSet> mcWDS; if (isMC) { mcWDS.reset(new RooDataSet( Form("mc%s_SET_%s", chg.c_str(), sample.c_str()) , Form("mc%s", chg.c_str()) , mcCols)); }
  // Set the MET
  for (int i = 0; i < dDS->numEntries(); i++) {
    dDS->get(i);
    metDS->get(i);
    // Extract the User MET
    auto metMag = (RooRealVar*)metCols.find(Form("MET_Mag_%s", metType.c_str()));
    if (metMag==NULL) { std::cout << "[ERROR] setMETonDS: Variable MET Mag of " << metType << " was not found in " << metDSName << std::endl; return false; }
    auto metPhi = (RooRealVar*)metCols.find(Form("MET_Phi_%s", metType.c_str()));
    if (metPhi==NULL) { std::cout << "[ERROR] setMETonDS: Variable MET Phi of " << metType << " was not found in " << metDSName << std::endl; return false; }
    // Extract the Nominal MET
    auto metMag_Nom = (RooRealVar*)cols.find("MET");
    if (metMag_Nom==NULL) { std::cout << "[ERROR] setMETonDS: Variable MET Mag was not found in " << dDSName << std::endl; return false; }
    // Set the Nominal MET to the User MET
    metMag_Nom->setVal( metMag->getVal() );
    // Fill the new RooDataSets
    dWDS->addFast(cols);
    //
    if (isMC) {
      mcDS->get(i);
      auto metPhi_Nom = (RooRealVar*)mcCols.find("MET_Phi");
      if (metPhi_Nom==NULL) { std::cout << "[ERROR] setMETonDS: Variable MET Phi was not found in " << mcDSName << std::endl; return false; }
      metPhi_Nom->setVal( metPhi->getVal() );
      mcWDS->addFast(mcCols);
    }
  }
  // Import to RooWorkspace
  ws.import(*dWDS);
  if (isMC) { ws.import(*mcWDS); }
  // Return
  return true;
  //
};


bool setMET(RooWorkspaceMap_t& Workspaces, const std::string& metType)
{
  //
  for (auto& ws : Workspaces) {
    RooWorkspace& myws = ws.second;
    // Get the sample name
    const auto& dsList =  myws.allData();
    std::string sample = dsList.front()->GetName(); sample = sample.substr(sample.find("_")+1);
    const std::string METType = (myws.obj("METType")) ? ((RooStringVar*)myws.obj("METType"))->getVal() : "";
    if (metType == METType) { std::cout << "[INFO] User define MET for sample " << sample << " is the default (" << METType << "), skipping the setMET" << std::endl; continue; }
    std::cout << "[INFO] Setting MET of sample " << sample << " from " << METType << " to " << metType << std::endl;
    // Apply the Lumi re-weighting to positive muon dataset
    if (!setMETonDS(myws, ("dPl_"+sample), ("mcPl_"+sample), ("metPl_"+sample), metType)) { return false; }
    // Apply the Lumi re-weighting to negative muon dataset
    if (!setMETonDS(myws, ("dMi_"+sample), ("mcMi_"+sample), ("metMi_"+sample), metType)) { return false; }
  }
  //
  return true;
  //
};


bool applyLumiWeight(RooWorkspace& ws, GlobalInfo& info, const std::string dDSName, const std::string& mcDSName, const bool verbose = false)
{
  //
  if (dDSName.find("MC_")==std::string::npos) { std::cout << "[ERROR] applyLumiWeight: Dataset " << dDSName << " is not from MC!" << std::endl; return false; }
  // Get the RooDatasets
  auto  dDS = (RooDataSet*)ws.data(dDSName.c_str());
  if ( dDS==NULL) { std::cout << "[ERROR] applyLumiWeight: Dataset " << dDSName << " was not found in the workspace!" << std::endl; return false; }
  auto mcDS = (RooDataSet*)ws.data(mcDSName.c_str());
  if (mcDS==NULL) { std::cout << "[ERROR] applyLumiWeight: Dataset " << mcDSName << " was not found in the workspace!" << std::endl; return false; }
  // Get the sample name
  std::string sample = dDS->GetName(); sample = sample.substr(sample.find("MC_"));
  std::string chg    = dDS->GetName(); chg    = chg.substr(1, chg.find("_")-1);
  std::string col    = dDS->GetName(); col    = col.substr(col.find_last_of("_")+1);
  if ( (col!="pPb") && (col!="Pbp") ) { std::cout << "[ERROR] Collision system " << col << " is invalid!" << std::endl; return false; }
  if ( (chg!="Pl" ) && (chg!="Mi" ) ) { std::cout << "[ERROR] Charge "           << chg << " is invalid!" << std::endl; return false; }
  // Initialize the new RooRealVars
  RooRealVar weight = RooRealVar ( "Weight" , "Weight"        , -1.0 , 1000000000.0, ""        );
  RooRealVar mcXSec = RooRealVar ( "XSec"   , "Cross Section" , -1.0 , 1000000000.0, "nb"      );
  RooRealVar mcLumi = RooRealVar ( "Lumi"   , "Luminosity"    , -1.0 , 1000000000.0, "nb^{-1}" );
  // Initialize the RooArgSets
  RooArgSet cols   = *dDS->get();
  cols.add(weight);
  RooArgSet mcCols = *mcDS->get();
  mcCols.add(mcXSec);
  mcCols.add(mcLumi);
  // Initialize the new RooDataSets
  auto dWDS  = std::unique_ptr<RooDataSet>(new RooDataSet( Form("d%s_LUM_%s" , chg.c_str(), sample.c_str()) , Form("d%s" , chg.c_str()) , cols   , RooFit::WeightVar(weight)) );
  auto mcWDS = std::unique_ptr<RooDataSet>(new RooDataSet( Form("mc%s_LUM_%s", chg.c_str(), sample.c_str()) , Form("mc%s", chg.c_str()) , mcCols) );
  // Apply the Lumi re-weighting to positive muon dataset
  int mcID = -99; double xSection = -1.0 , luminosity = -1.0;
  for (int i = 0; i < dDS->numEntries(); i++) {
    dDS->get(i);
    mcDS->get(i);
    // Extract the Cross-Section
    auto mcType = (RooRealVar*)mcCols.find("MC_Type");
    if (mcType==NULL) { std::cout << "[ERROR] applyLumiWeight: Variable MC_Type was not found in " << mcDSName << std::endl; return false; }
    if (mcID != mcType->getVal()) { // Do every time there is a change of cross-section
      // Get the Luminosity
      luminosity = ( (col=="pPb") ? PA::LUMI::Data_pPb : PA::LUMI::Data_Pbp );
      if (info.Var.count("rLumi")>0) { luminosity *= info.Var.at("rLumi").at("Val"); }
      if (info.Var.count(("Luminosity_"+col).c_str())==0) { info.Var["Luminosity_"+col]["Val"] = luminosity; }
      // Get the Cross-Section
      const std::string mcXSecName = PA::getMCXSecName(int(mcType->getVal()));
      const std::string mcSample   = mcXSecName.substr(mcXSecName.find("_")+1);
      const std::string mcTag      = mcSample.substr(0, mcSample.find_last_of("_"));
      if (!PA::getCrossSection(xSection, mcXSecName)) { return false; }
      if (info.Var.count(("rXSection_"+mcTag).c_str())>0) { xSection *= info.Var.at(("rXSection_"+mcTag).c_str()).at("Val"); }
      mcID = mcType->getVal();
      if (info.Var.count(("XSection_"+mcSample).c_str())==0) { info.Var["XSection_"+mcSample]["Val"] = xSection; }
      std::string oStr = std::string("[INFO] ") + mcXSecName + " set to : Cross-section " + Form("%g", xSection) + " nb ";
      if (info.Var.count(("rXSection_"+mcTag).c_str())>0) {
        const double v = info.Var.at(("rXSection_"+mcTag).c_str()).at("Val");
        oStr += std::string("(") + Form("%s%.0f%%", ((v>=1.0) ? "+" : "-"), std::abs((v-1.0)*100.0)) + ") ";
      }
      oStr += std::string("and Luminosity ") + Form("%g", luminosity) + " /nb ";
      if (info.Var.count("rLumi")>0) {
        const double v = info.Var.at("rLumi").at("Val");
        oStr += std::string("(") + Form("%s%.0f%%", ((v>=1.0) ? "+" : "-"), std::abs((v-1.0)*100.0)) + ") ";
      }
      if (verbose) { std::cout << oStr << std::endl; }
    }
    mcLumi.setVal( luminosity );
    mcXSec.setVal( xSection   );
    // Compute the Luminosity MC Weight
    auto mcNGen = (RooRealVar*)mcCols.find("NGen");
    if (mcNGen==NULL) { std::cout << "[ERROR] applyLumiWeight: Variable NGen was not found in " << mcDSName << std::endl; return false; }
    weight.setVal( (mcXSec.getVal() * mcLumi.getVal()) / mcNGen->getVal() );
    // Fill the new RooDataSets
    if (weight.getVal() <= 0.0) { std::cout << "[ERROR] applyLumiWeight : Weight is negative ( " << weight.getVal() << " ) in " << mcDSName << std::endl; return false; }
    dWDS->addFast(cols, weight.getVal());
    mcWDS->addFast(mcCols);
  }
  // Import to RooWorkspace
  ws.import(*dWDS);
  ws.import(*mcWDS);
  // Return
  return true;
  //
};


bool reweightMCLumi(RooWorkspaceMap_t& Workspaces, GlobalInfo& info)
{
  //
  for (auto& ws : Workspaces) {
    if (ws.first.find("MC_")!=std::string::npos) {
      RooWorkspace& myws = ws.second;
      // Get the sample name
      std::string sample = "SET_" + ws.first;
      if (myws.data(Form("dPl_%s", sample.c_str()))==NULL) { sample = "RAW_" + ws.first; }
      if (myws.data(Form("dPl_%s", sample.c_str()))==NULL) { std::cout << "[ERROR] Sample " << sample << " was not found in the workspace!" << std::endl; return false; }
      std::cout << "[INFO] Re-weighting the Luminosity of " << sample << std::endl;
      // Apply the Lumi re-weighting to positive muon dataset
      if (!applyLumiWeight(myws, info, ("dPl_"+sample), ("mcPl_"+sample), true)) { return false; }
      // Apply the Lumi re-weighting to negative muon dataset
      const bool printXSec = (sample.find("MC_W_")!=std::string::npos);
      if (!applyLumiWeight(myws, info, ("dMi_"+sample), ("mcMi_"+sample), printXSec)) { return false; }
    }
  }
  //
  return true;
  //
};


bool readCorrectionDS(RooWorkspace& corrWS, const std::string sampleTag, const std::string recoilMethod, const bool applyHFCorr)
{
  //
  bool done = false;
  std::string sample = sampleTag; sample = sample.substr(sample.find("MC_"));
  const std::string corrDirName  = Form("%s/Correction/Recoil_%s", getcwd(NULL, 0), recoilMethod.c_str());
  const std::string corrFileName = Form("Corr_%s%s_Recoil_%s.root", sample.c_str(), (applyHFCorr?"_HFCorr":""), recoilMethod.c_str());
  const std::string corrFilePath = Form("%s/%s", corrDirName.c_str(), corrFileName.c_str());
  //
  if (existFile(corrFilePath)) {
    TFile corrFile(corrFilePath.c_str(), "READ");
    if (!corrFile.IsZombie() && corrFile.IsOpen() && corrFile.Get("workspace")) {
      std::cout << "[INFO] Corrected MET extracted from " << corrFilePath << std::endl;
      const RooWorkspace& ws = *(RooWorkspace*)corrFile.Get("workspace");
      if ( (ws.data(("mcPl_COR_"+sample).c_str()) && ws.data(("mcMi_COR_"+sample).c_str())) && (ws.data(("dPl_COR_"+sample).c_str()) && ws.data(("dMi_COR_"+sample).c_str())) ) {
        corrWS.import(*ws.data(("mcPl_COR_"+sample).c_str()));  corrWS.import(*ws.data(("mcMi_COR_"+sample).c_str()));
        corrWS.import(*ws.data(("dPl_COR_"+sample).c_str()));   corrWS.import(*ws.data(("dMi_COR_"+sample).c_str()));
        if (ws.data(("metPl_COR_"+sample).c_str()) && ws.data(("metMi_COR_"+sample).c_str())) {
          corrWS.import(*ws.data(("metPl_COR_"+sample).c_str()));  corrWS.import(*ws.data(("metMi_COR_"+sample).c_str()));
        }
        done = true;
      }
      else { std::cout << "[WARNING] Some correction datasets are missing, ignoring file " << corrFilePath << std::endl; }
    }
    corrFile.Close();
  }
  // Return
  return done;
};


void writeCorrectionDS(const RooWorkspace& ws, const std::string sampleTag, const std::string recoilMethod, const bool applyHFCorr)
{
  //
  std::string sample = sampleTag; sample = sample.substr(sample.find("MC_"));
  const std::string corrDirName  = Form("%s/Correction/Recoil_%s", getcwd(NULL, 0), recoilMethod.c_str());
  const std::string corrFileName = Form("Corr_%s%s_Recoil_%s.root", sample.c_str(), (applyHFCorr?"_HFCorr":""), recoilMethod.c_str());
  const std::string corrFilePath = Form("%s/%s", corrDirName.c_str(), corrFileName.c_str());
  if (!existDir(corrDirName)) { makeDir(corrDirName); }
  //
  if ( (ws.data(("mcPl_COR_"+sample).c_str()) && ws.data(("mcMi_COR_"+sample).c_str())) && (ws.data(("dPl_COR_"+sample).c_str()) && ws.data(("dMi_COR_"+sample).c_str())) ) {
    std::cout << "[INFO] Corrected MET stored in " << corrFilePath << std::endl;
    RooWorkspace corrWS;
    corrWS.import(*ws.data(("mcPl_COR_"+sample).c_str()));  corrWS.import(*ws.data(("mcMi_COR_"+sample).c_str()));
    corrWS.import(*ws.data(("dPl_COR_"+sample).c_str()));   corrWS.import(*ws.data(("dMi_COR_"+sample).c_str()));
    if (ws.data(("metPl_COR_"+sample).c_str()) && ws.data(("metMi_COR_"+sample).c_str())) {
      corrWS.import(*ws.data(("metPl_COR_"+sample).c_str()));  corrWS.import(*ws.data(("metMi_COR_"+sample).c_str()));
    }
    TFile corrFile(corrFilePath.c_str(), "RECREATE"); if (!corrFile.IsZombie() && corrFile.IsOpen()) { corrFile.cd(); corrWS.Write("workspace"); corrFile.Write(); }; corrFile.Close();
  }
  //
};


bool applyMCCorrection(RooWorkspace& ws, RooWorkspace& corrWS, const RooWorkspace& iniWS, const std::string dDSName, const std::string& mcDSName, 
                       const bool& applyTnPCorr, const bool& applyBosonPTCorr, HFweight* HFCorr, const std::string& HFMethod, RecoilCorrector& recoilCorr, const std::string& recoilMethod)
{
  //
  const bool applyHFCorr     = (HFCorr != NULL);
  const bool applyRecoilCorr = recoilCorr.isValid();
  if (!applyBosonPTCorr && !applyTnPCorr && !applyHFCorr && !applyRecoilCorr) { copyWorkspace(ws, iniWS, dDSName, false); return true; }
  if (dDSName.find("MC_")==std::string::npos) { std::cout << "[ERROR] applyMCCorrection: Dataset " << dDSName << " is not from MC!" << std::endl; return false; }
  // Get the RooDatasets
  auto  dDS = (RooDataSet*)iniWS.data(dDSName.c_str());
  if ( dDS==NULL) { std::cout << "[ERROR] applyMCCorrection: Dataset " << dDSName << " was not found in the workspace!" << std::endl; return false; }
  auto mcDS = (RooDataSet*)iniWS.data(mcDSName.c_str());
  if (mcDS==NULL) { std::cout << "[ERROR] applyMCCorrection: Dataset " << mcDSName << " was not found in the workspace!" << std::endl; return false; }
  // Get the sample name
  std::string sample = dDS->GetName(); sample = sample.substr(sample.find("MC_"));
  std::string chg    = dDS->GetName(); chg    = chg.substr(1, sample.find("_"));
  std::string col    = dDS->GetName(); col    = col.substr(col.find_last_of("_")+1);
  // Initialize the new RooRealVars
  RooRealVar mcRawMET = RooRealVar( "MET_RAW"     , "|#slash{E}_{T}|"          ,  -1.0 ,     100000.0 , "GeV/c"     );
  RooRealVar mcRawMT  = RooRealVar( "Muon_MT_RAW" , "W Transverse Mass"        ,  -1.0 ,     100000.0 , "GeV/c^{2}" );
  RooRealVar mcSFTnP  = RooRealVar( "SFTnP"       , "TagAndProbe Scale Factor" , -10.0 ,         10.0 , ""          );
  RooRealVar mcWEA    = RooRealVar( "wEA"         , "Event Activity Weight"    ,  -1.0 , 1000000000.0 , ""          );
  RooRealVar mcWBPT   = RooRealVar( "wBPT"        , "Boson pT Weight"          ,  -1.0 , 1000000000.0 , ""          );
  RooRealVar weight   = RooRealVar( "Weight"      , "Weight"                   ,  -1.0 , 1000000000.0 , ""          );
  // Initialize the RooArgSets
  RooArgSet cols   = *dDS->get();
  cols.add(weight);
  RooArgSet mcCols = *mcDS->get();
  if (applyBosonPTCorr ) { mcCols.add(mcWBPT);  }
  if (applyTnPCorr     ) { mcCols.add(mcSFTnP); }
  if (applyHFCorr      ) { mcCols.add(mcWEA);   }
  if (applyRecoilCorr  ) { mcCols.add(mcRawMET); mcCols.add(mcRawMT); }
  //
  // Initialize the new RooDataSets
  const std::string dCorrDSName  = Form("d%s_COR_%s" , chg.c_str(), sample.c_str());
  const std::string mcCorrDSName = Form("mc%s_COR_%s" , chg.c_str(), sample.c_str());
  auto dWDS  = std::unique_ptr<RooDataSet>(new RooDataSet( dCorrDSName.c_str()  , Form("d%s" , chg.c_str()) , cols   , RooFit::WeightVar(weight)) );
  auto mcWDS = std::unique_ptr<RooDataSet>(new RooDataSet( mcCorrDSName.c_str() , Form("mc%s", chg.c_str()) , mcCols) );
  // Extract the correction RooDataSet
  RooDataSet* corrDS   = NULL; if ( corrWS.data(dCorrDSName.c_str())  ) { corrDS   = (RooDataSet*)corrWS.data(dCorrDSName.c_str());  }
  RooDataSet* corrMCDS = NULL; if ( corrWS.data(mcCorrDSName.c_str()) ) { corrMCDS = (RooDataSet*)corrWS.data(mcCorrDSName.c_str()); }
  if ( (corrDS==NULL && corrMCDS!=NULL) || (corrDS!=NULL && corrMCDS==NULL) ) { std::cout << "[ERROR] One of the correction datasets is missing" << std::endl; return false; }
  if (corrDS!=NULL && corrMCDS!=NULL && corrDS->numEntries()!=corrMCDS->numEntries()) { std::cout << "[ERROR] Correction datasets have different number of entries" << std::endl; return false; }
  if (corrDS!=NULL && corrDS->numEntries()!=dDS->numEntries()) { std::cout << "[ERROR] Wrong number of entries corrDS (" << corrDS->numEntries() << ") and dDS (" << dDS->numEntries() << ")" << std::endl; return true; }
  // Apply the Lumi re-weighting to positive muon dataset
  int mcID = -99;
  for (int i = 0; i < dDS->numEntries(); i++) {
    dDS->get(i);
    mcDS->get(i);
    weight.setVal( dDS->weight() );
    //
    // Apply the Boson pT Re-Weighting
    //
    if (applyBosonPTCorr) {
      // Get the Boson pT
      auto bosonPt = (RooRealVar*)mcCols.find("Boson_Pt");
      if (bosonPt==NULL) { std::cout << "[ERROR] applyMCCorrection: Variable Boson_Pt was not found in " << dDSName << std::endl; return false; }
      double bosonPT = bosonPt->getVal(); if (bosonPT < 0.5) { bosonPT = 0.5; } // Weights derived down to 0.5 GeV
      // Compute the Boson pT weight
      const double weightBPT = ( 1.0 / ( ( -0.37 * std::pow(bosonPT, -0.37) ) + 1.19 ) );
      if (weightBPT==0.) { std::cout << "[ERROR] Weight is zero for boson pT " << bosonPT << " in " << dDSName << std::endl; return false; }
      if (weightBPT <0.) { std::cout << "[ERROR] Weight is negative for boson pT " << bosonPT << " in " << dDSName << std::endl; return false; }
      mcWBPT.setVal( weightBPT );
      // Set the Boson pT weight
      weight.setVal( weight.getVal() * mcWBPT.getVal() );
    }
    //
    // Apply the Tag And Probe corrections
    //
    if (applyTnPCorr) {
      // Get the Muon Kinematic Variables
      auto muPt = (RooRealVar*)cols.find("Muon_Pt");
      if (muPt==NULL) { std::cout << "[ERROR] applyMCCorrection: Variable Muon_Pt was not found in " << dDSName << std::endl; return false; }
      auto muEta = (RooRealVar*)cols.find("Muon_Eta");
      if (muEta==NULL) { std::cout << "[ERROR] applyMCCorrection: Variable Muon_Eta was not found in " << dDSName << std::endl; return false; }
      // Get the Nominal Tag and Probe Scale Factors
      const double sf_MuID = tnp_weight_muid_ppb( muPt->getVal() , muEta->getVal() , 0 );
      const double sf_Trig = tnp_weight_trg_ppb (                  muEta->getVal() , 0 );
      const double sf_Iso  = tnp_weight_iso_ppb ( muPt->getVal() , muEta->getVal() , 0 );
      mcSFTnP.setVal( sf_MuID * sf_Trig * sf_Iso );
      // Set the TnP Scale Factor weight
      weight.setVal( weight.getVal() * mcSFTnP.getVal() );
    }
    //
    // Apply Event Activity Re-Weighting
    //
    if (applyHFCorr) {
      double weightEA = 1.0;
      if (HFMethod=="HFBoth") {
        // Get the HF Energy
        auto hiHF = (RooRealVar*)mcCols.find("hiHF");
        if (hiHF==NULL) { std::cout << "[ERROR] applyMCCorrection: Variable HF was not found in " << dDSName << std::endl; return false; }
        // Get the Event Activity Weight
        weightEA = HFCorr->weight(hiHF->getVal(), HFweight::HFside::both, false);
        if (weightEA==0.) { std::cout << "[ERROR] Weight is zero for HF " << hiHF << std::endl; return false; }
      }
      else if (HFMethod=="NTracks") {
        // Get the Number of Tracks
        auto nTracks = (RooRealVar*)mcCols.find("nTracks");
        if (nTracks==NULL) { std::cout << "[ERROR] applyMCCorrection: Variable nTracks was not found in " << dDSName << std::endl; return false; }
        // Get the Event Activity Weight
        const double nTrks = ( (nTracks->getVal() >=300.) ? 290. : nTracks->getVal() ); // The histograms are made up to 300.
        weightEA = HFCorr->weight(nTrks, HFweight::HFside::track, false);
        if (weightEA==0.) { std::cout << "[ERROR] Weight is zero for nTracks " << nTrks << std::endl; return false; }
      }
      else { std::cout << "[ERROR] Event Activity correction method " << HFMethod << " has not been defined!" << std::endl; return false; }
      mcWEA.setVal( weightEA );
      // Set the TnP Scale Factor weight
      weight.setVal( weight.getVal() * mcWEA.getVal() );
    }
    //
    // Apply the MET Recoil corrections
    //
    if (applyRecoilCorr) {
      // Get the Muon Kinematic Variables
      auto muPt = (RooRealVar*)cols.find("Muon_Pt");
      if (muPt==NULL) { std::cout << "[ERROR] applyMCCorrection: Variable Muon_Pt was not found in " << dDSName << std::endl; return false; }
      auto muPhi = (RooRealVar*)mcCols.find("Muon_Phi");
      if (muPhi==NULL) { std::cout << "[ERROR] applyMCCorrection: Variable Muon_Phi was not found in " << dDSName << std::endl; return false; }
      // Get the MET vector
      auto met = (RooRealVar*)cols.find("MET");
      if (met==NULL) { std::cout << "[ERROR] applyMCCorrection: Variable MET was not found in " << dDSName << std::endl; return false; }
      mcRawMET.setVal(met->getVal());
      auto metPhi = (RooRealVar*)mcCols.find("MET_Phi");
      if (metPhi==NULL) { std::cout << "[ERROR] applyMCCorrection: Variable MET_Phi was not found in " << dDSName << std::endl; return false; }
      TVector2 MET_RAW; MET_RAW.SetMagPhi(met->getVal() , metPhi->getVal());
      // Get the Muon Transverse Mass
      auto muMT = (RooRealVar*)cols.find("Muon_MT");
      if (muMT==NULL) { std::cout << "[ERROR] applyMCCorrection: Variable Muon_MT was not found in " << dDSName << std::endl; return false; }
      mcRawMT.setVal(muMT->getVal());
      // Get the Reference pT vector
      auto refPt = (RooRealVar*)mcCols.find("Ref_Pt");
      if (refPt==NULL) { std::cout << "[ERROR] applyMCCorrection: Variable Ref_Pt was not found in " << dDSName << std::endl; return false; }
      auto refPhi = (RooRealVar*)mcCols.find("Ref_Phi");
      if (refPhi==NULL) { std::cout << "[ERROR] applyMCCorrection: Variable Ref_Phi was not found in " << dDSName << std::endl; return false; }
      TVector2 reference_pT; reference_pT.SetMagPhi(refPt->getVal() , refPhi->getVal());
      // Get the Boson pT vector
      auto bosonPt = (RooRealVar*)mcCols.find("Boson_Pt");
      if (bosonPt==NULL) { std::cout << "[ERROR] applyMCCorrection: Variable Boson_Pt was not found in " << dDSName << std::endl; return false; }
      auto bosonPhi = (RooRealVar*)mcCols.find("Boson_Phi");
      if (bosonPhi==NULL) { std::cout << "[ERROR] applyMCCorrection: Variable Boson_Phi was not found in " << dDSName << std::endl; return false; }
      TVector2 boson_pT; boson_pT.SetMagPhi(bosonPt->getVal() , bosonPhi->getVal());
      // Check correction dataset
      if (corrDS) {
        const RooArgSet& corrCols   = *corrDS->get(i);
        const RooArgSet& corrMCCols = *corrMCDS->get(i);
        if ( ( muPt->getVal()     != getVar(corrCols,   "Muon_Pt").getVal()     ) || ( muPhi->getVal()    != getVar(corrMCCols, "Muon_Phi").getVal()  ) ||
             ( met->getVal()      != getVar(corrMCCols, "MET_RAW").getVal()     ) || ( metPhi->getVal()   != getVar(corrMCCols, "MET_Phi").getVal()   ) ||
             ( refPt->getVal()    != getVar(corrMCCols, "Ref_Pt").getVal()      ) || ( refPhi->getVal()   != getVar(corrMCCols, "Ref_Phi").getVal()   ) ||
             ( bosonPt->getVal()  != getVar(corrMCCols, "Boson_Pt").getVal()    ) || ( bosonPhi->getVal() != getVar(corrMCCols, "Boson_Phi").getVal() ) ||
             ( muMT->getVal()     != getVar(corrMCCols, "Muon_MT_RAW").getVal() ) ) {
          std::cout << "[ERROR] Correction dataset " << mcCorrDSName << " is imcompatible with current dataset" << std::endl; return false;
        }
        else {
          met->setVal(getVar(corrCols, "MET").getVal());
          muMT->setVal(getVar(corrCols, "Muon_MT").getVal());
        }
      }
      else {
        // Set the Reference and Boson pT vectors
        recoilCorr.setPt(reference_pT, boson_pT);
        // Correct the MET
        TVector2 MET_CORR;
        std::string recoilMethodLbl = "";
        if (recoilMethod.find("Scaling" )!=std::string::npos) { recoilMethodLbl = "Scaling";  }
        if (recoilMethod.find("Smearing")!=std::string::npos) { recoilMethodLbl = "Smearing"; }
        if (!recoilCorr.correctMET(MET_CORR, MET_RAW, recoilMethodLbl)) { return false; }
        // Correct the Muon Transverse Mass
        const double muMTVal = PA::getWTransverseMass(muPt->getVal(), muPhi->getVal(), MET_CORR.Mod(), MET_CORR.Phi());
        // Set the corrected MET and Transverse Mass
        met->setVal(MET_CORR.Mod());
        muMT->setVal(muMTVal);
      }
    }
    // Fill the new RooDataSets
    if (weight.getVal() <= 0.0) { std::cout << "[ERROR] applyMCCorrection: Weight is negative ( " << weight.getVal() << " ) in " << dDSName << std::endl; return false; }
    dWDS->addFast(cols, weight.getVal());
    mcWDS->addFast(mcCols);
  }
  // Import to RooWorkspace
  ws.import(*dWDS);
  if (corrDS==NULL) {
    corrWS.import(*dWDS);
    corrWS.import(*mcWDS);
  }
  // Return
  return true;
};


bool correctMC(RooWorkspaceMap_t& Workspaces, const RooWorkspaceMap_t& iniWorkspaces, const GlobalInfo& info)
{
  //
  const bool applyBosonPTCorr = info.Flag.at("applyBosonPTCorr");
  const bool applyTnPCorr     = info.Flag.at("applyTnPCorr");
  const bool applyHFCorr      = info.Flag.at("applyHFCorr");
  const bool applyRecoilCorr  = info.Flag.at("applyRecoilCorr");
  //
  // Define the HF Weight
  std::unique_ptr<HFweight> HFCorr;
  std::string HFMethod;
  if (applyHFCorr) {
    HFMethod = info.Par.at("HFCorrMethod");
    HFCorr = std::unique_ptr<HFweight>(new HFweight("/afs/cern.ch/work/e/echapon/public/DY_pA_2016/HFweight.root"));
  }
  //
  // Define the MET Recoil Corrector
  RecoilCorrector recoilCorr(1);
  std::string recoilMethod;
  bool ignoreCorrDS = false;
  if (applyRecoilCorr) {
    recoilMethod = info.Par.at("RecoilCorrMethod");
    const std::string met = info.Par.at("RecoilMET");
    const std::string recoilPath = info.Par.at("RecoilPath");
    const std::string col = "PA";
    std::string fnc_MC = "doubleGauss" , fnc_DATA = "doubleGauss";
    if (recoilPath.find("BWGauss")!=std::string::npos) { fnc_MC   = "BWGauss" , fnc_DATA = "BWGauss"; }
    const std::string HFCorrLbl = (applyHFCorr ? "HFCorr" : "noHFCorr");
    std::string fileName_MC   = Form("%sMC_DYToMuMu_POWHEG/MET_%s/%s/%s/%s/Results/fits_RecoilPDF_%s_%s.root", recoilPath.c_str(), met.c_str(), col.c_str(), HFCorrLbl.c_str(), fnc_MC.c_str(), met.c_str(), col.c_str());
    std::string fileName_DATA = Form("%sDATA/MET_%s/%s/%s/Results/fits_RecoilPDF_%s_%s.root", recoilPath.c_str(), met.c_str(), col.c_str(), fnc_DATA.c_str(), met.c_str(), col.c_str());
    //
    const std::string recoilVarTyp = (info.Par.count("RecoilVarTyp")>0 ? info.Par.at("RecoilVarTyp") : "");
    const std::string recoilVarLbl = (info.Par.count("RecoilVarLbl")>0 ? info.Par.at("RecoilVarLbl") : "");
    if ((recoilVarTyp!="" && recoilVarLbl!="") || info.Flag.at("ignoreCorrDS")) { ignoreCorrDS = true; }
    if (recoilVarTyp.find("MC"  )!=std::string::npos) {
      fileName_MC = Form("%sMC_DYToMuMu_POWHEG/MET_%s/%s/%s/%s/Systematics/Sys_%s/fits_RecoilPDF_%s_%s.root", recoilPath.c_str(), met.c_str(), col.c_str(), HFCorrLbl.c_str(),
                         fnc_MC.c_str(), recoilVarLbl.c_str(), met.c_str(), col.c_str());
    }
    if (recoilVarTyp.find("DATA")!=std::string::npos) {
      fileName_DATA = Form("%sDATA/MET_%s/%s/%s/Systematics/Sys_%s/fits_RecoilPDF_%s_%s.root", recoilPath.c_str(), met.c_str(), col.c_str(),
                           fnc_DATA.c_str(), recoilVarLbl.c_str(), met.c_str(), col.c_str());
    }
    //
    if (!recoilCorr.setInputFiles(met, fileName_MC, fileName_DATA)) { return false; }
  }
  //
  for (const auto& ws : iniWorkspaces) {
    auto&  myws = Workspaces[ws.first];
    // Get the sample name
    std::string sample = "LUM_" + ws.first;
    if (ws.second.data(Form("dPl_%s", sample.c_str()))==NULL) { sample = "SET_" + ws.first; }
    if (ws.second.data(Form("dPl_%s", sample.c_str()))==NULL) { sample = "RAW_" + ws.first; }
    if (ws.second.data(Form("dPl_%s", sample.c_str()))==NULL) { std::cout << "[ERROR] Sample " << sample << " was not found in the workspace!" << std::endl; return false; }
    //
    if (ws.first.find("MC_")!=std::string::npos) {
      auto& iniWS = ws.second;
      const std::string corrString = std::string(applyBosonPTCorr ? " BosonPTCorr " : "") + std::string(applyTnPCorr ? " TnPCorr " : "") + std::string(applyHFCorr ? " HFCorr " : "")  + std::string(applyRecoilCorr ? " RecoilCorr " : "");
      std::cout << "[INFO] Applying corrections ( " << corrString << " ) to " << sample << std::endl;
      // Extract the corrections
      bool makeCorrFile = false; RooWorkspace corrWS;
      if (applyRecoilCorr && !ignoreCorrDS) { makeCorrFile = (!readCorrectionDS(corrWS, sample, recoilMethod, applyHFCorr)); }
      // Initialize the MET Recoil Corrector
      if (applyRecoilCorr) {
        if      (recoilMethod.find("OneGaussianMC"  )!=std::string::npos) { recoilCorr.setInitialSetup(sample, true, false, true); }
        else if (recoilMethod.find("OneGaussianDATA")!=std::string::npos) { recoilCorr.setInitialSetup(sample, false, true, true); }
        else if (recoilMethod.find("OneGaussian"    )!=std::string::npos) { recoilCorr.setInitialSetup(sample, true, true, true);  }
        else { recoilCorr.setInitialSetup(sample, false, false, true); }
      }
      // Apply the MC correction to positive muon dataset
      if (!applyMCCorrection(myws, corrWS, iniWS, ("dPl_"+sample), ("mcPl_"+sample), applyTnPCorr, applyBosonPTCorr, HFCorr.get(), HFMethod, recoilCorr, recoilMethod)) { return false; }
      // Apply the MC correction to negative muon dataset
      if (!applyMCCorrection(myws, corrWS, iniWS, ("dMi_"+sample), ("mcMi_"+sample), applyTnPCorr, applyBosonPTCorr, HFCorr.get(), HFMethod, recoilCorr, recoilMethod)) { return false; }
      // Store in the workspace the info regarding the corrections applied
      RooStringVar tmp; tmp.setVal(corrString.c_str()); tmp.SetTitle("CorrectionApplied"); myws.import(*((TObject*)&tmp), tmp.GetTitle());
      if (applyRecoilCorr) { tmp.setVal(recoilMethod.c_str()); tmp.SetTitle("RecoilMethod"); myws.import(*((TObject*)&tmp), tmp.GetTitle()); }
      copyWorkspace(myws, ws.second, "", false);
      // Save the corrections
      if (applyRecoilCorr && makeCorrFile && !ignoreCorrDS) {
	copyWorkspace(corrWS, myws, "", false);
	writeCorrectionDS(corrWS, sample, recoilMethod, applyHFCorr);
      }
    }
    else {
      copyWorkspace(myws, ws.second, Form("dPl_%s", sample.c_str()), false);
      copyWorkspace(myws, ws.second, Form("dMi_%s", sample.c_str()), false);
    }
  }
  //
  return true;
  //
};


bool invertEtaAndFill(RooDataSet& dsPA, const RooDataSet& dsPbp)
{
  for(int i = 0; i < dsPbp.numEntries(); i++){
    auto set = dsPbp.get(i);
    auto eta = (RooRealVar*)set->find("Muon_Eta");
    if (eta==NULL) { std::cout << "[ERROR] RooRealVar Muon_Eta was not found!" << std::endl; return false; }
    // Invert the pseudo-rapidity
    eta->setVal(-1.0*eta->getVal());
    // Add the new event
    if (dsPbp.weight() <= 0.0) { std::cout << "[ERROR] invertEtaAndFill: Weight is negative ( " << dsPbp.weight() << " )" << std::endl; return false; }
    dsPA.add(*set, dsPbp.weight());
  }
  return true;
};


bool createPADataset(RooWorkspaceMap_t& Workspaces, const std::string& sampleTag)
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
  std::string dsType = "COR";
  if (Workspaces.at(sample_pPb).data(("dPl_"+dsType+"_"+sample_pPb).c_str())==NULL) { dsType = "LUM"; }
  if (Workspaces.at(sample_pPb).data(("dPl_"+dsType+"_"+sample_pPb).c_str())==NULL) { dsType = "SET"; }
  if (Workspaces.at(sample_pPb).data(("dPl_"+dsType+"_"+sample_pPb).c_str())==NULL) { dsType = "RAW"; }
  if (Workspaces.at(sample_pPb).data(("dPl_"+dsType+"_"+sample_pPb).c_str())==NULL) { std::cout << "[ERROR] Sample " << sampleTag << " was not found!" << std::endl; return false; }
  const std::string dsSampleTag  = dsType + "_" + sampleTag;
  const std::string dsSample_pPb = dsSampleTag + "_pPb";
  const std::string dsSample_Pbp = dsSampleTag + "_Pbp";
  const std::string dsSample_PA  = dsSampleTag + "_PA";
  //
  auto dPl_pPb = (RooDataSet*)Workspaces.at(sample_pPb).data(("dPl_"+dsSample_pPb).c_str());
  if (dPl_pPb==NULL) { std::cout << "[ERROR] RooDataSet " << ("dPl_"+dsSample_pPb) << " does not exist!" << std::endl; return false; }
  auto dMi_pPb = (RooDataSet*)Workspaces.at(sample_pPb).data(("dMi_"+dsSample_pPb).c_str());
  if (dMi_pPb==NULL) { std::cout << "[ERROR] RooDataSet " << ("dMi_"+dsSample_pPb) << " does not exist!" << std::endl; return false; }
  auto  mcPl_pPb = (RooDataSet*)Workspaces.at(sample_pPb).data(("mcPl_"+dsSample_pPb).c_str());
  auto  mcMi_pPb = (RooDataSet*)Workspaces.at(sample_pPb).data(("mcMi_"+dsSample_pPb).c_str());
  auto metPl_pPb = (RooDataSet*)Workspaces.at(sample_pPb).data(("metPl_"+dsSample_pPb).c_str());
  auto metMi_pPb = (RooDataSet*)Workspaces.at(sample_pPb).data(("metMi_"+dsSample_pPb).c_str());
  //
  auto dPl_Pbp = (RooDataSet*)Workspaces.at(sample_Pbp).data(("dPl_"+dsSample_Pbp).c_str());
  if (dPl_Pbp==NULL) { std::cout << "[ERROR] RooDataSet " << ("dPl_"+dsSample_Pbp) << " does not exist!" << std::endl; return false; }
  auto dMi_Pbp = (RooDataSet*)Workspaces.at(sample_Pbp).data(("dMi_"+dsSample_Pbp).c_str());
  if (dMi_Pbp==NULL) { std::cout << "[ERROR] RooDataSet " << ("dMi_"+dsSample_Pbp) << " does not exist!" << std::endl; return false; }
  auto  mcPl_Pbp = (RooDataSet*)Workspaces.at(sample_Pbp).data(("mcPl_"+dsSample_Pbp).c_str());
  auto  mcMi_Pbp = (RooDataSet*)Workspaces.at(sample_Pbp).data(("mcMi_"+dsSample_Pbp).c_str());
  auto metPl_Pbp = (RooDataSet*)Workspaces.at(sample_Pbp).data(("metPl_"+dsSample_Pbp).c_str());
  auto metMi_Pbp = (RooDataSet*)Workspaces.at(sample_Pbp).data(("metMi_"+dsSample_Pbp).c_str());
  //
  std::cout << "[INFO] Creating combined PA RooDataSet for " << dsSampleTag << std::endl;
  //
  // Copy the pPb datasets
  auto dPl_PA = std::unique_ptr<RooDataSet>((RooDataSet*)dPl_pPb->Clone(("dPl_"+dsSample_PA).c_str()));
  if (dPl_PA==NULL || dPl_PA->sumEntries()==0) { std::cout << "[ERROR] RooDataSet " << ("dPl_"+dsSample_PA) << " was not created!" << std::endl; return false; }
  auto dMi_PA = std::unique_ptr<RooDataSet>((RooDataSet*)dMi_pPb->Clone(("dMi_"+dsSample_PA).c_str()));
  if (dMi_PA==NULL || dMi_PA->sumEntries()==0) { std::cout << "[ERROR] RooDataSet " << ("dMi_"+dsSample_PA) << " was not created!" << std::endl; return false; }
  //
  std::unique_ptr<RooDataSet>  mcPl_PA = NULL; if ( mcPl_pPb!=NULL) {  mcPl_PA = std::unique_ptr<RooDataSet>((RooDataSet*) mcPl_pPb->Clone(( "mcPl_"+dsSample_PA).c_str())); }
  std::unique_ptr<RooDataSet>  mcMi_PA = NULL; if ( mcMi_pPb!=NULL) {  mcMi_PA = std::unique_ptr<RooDataSet>((RooDataSet*) mcPl_pPb->Clone(( "mcMi_"+dsSample_PA).c_str())); }
  std::unique_ptr<RooDataSet> metPl_PA = NULL; if (metPl_pPb!=NULL) { metPl_PA = std::unique_ptr<RooDataSet>((RooDataSet*)metPl_pPb->Clone(("metPl_"+dsSample_PA).c_str())); }
  std::unique_ptr<RooDataSet> metMi_PA = NULL; if (metMi_pPb!=NULL) { metMi_PA = std::unique_ptr<RooDataSet>((RooDataSet*)metPl_pPb->Clone(("metMi_"+dsSample_PA).c_str())); }
  //
  // Invert the eta of Pbp dataset and fill the PA dataset
  if (!invertEtaAndFill(*dPl_PA, *dPl_Pbp)) { return false; }
  if (!invertEtaAndFill(*dMi_PA, *dMi_Pbp)) { return false; }
  if ( mcPl_Pbp!=NULL &&  mcPl_PA!=NULL) {  mcPl_PA->append( *mcPl_Pbp); }
  if ( mcMi_Pbp!=NULL &&  mcMi_PA!=NULL) {  mcMi_PA->append( *mcMi_Pbp); }
  if (metPl_Pbp!=NULL && metPl_PA!=NULL) { metPl_PA->append(*metPl_Pbp); }
  if (metMi_Pbp!=NULL && metMi_PA!=NULL) { metMi_PA->append(*metMi_Pbp); }
  //
  // Check the consistency of the combined dataset
  if (dPl_PA->numEntries()!=(dPl_pPb->numEntries()+dPl_Pbp->numEntries())) { std::cout << "[ERROR] Number of entries for the combined Pl " << dsSampleTag << " dataset is inconsistent!" << std::endl; return false; }
  if (dMi_PA->numEntries()!=(dMi_pPb->numEntries()+dMi_Pbp->numEntries())) { std::cout << "[ERROR] Number of entries for the combined Mi " << dsSampleTag << " dataset is inconsistent!" << std::endl; return false; }
  if (std::abs(dPl_PA->sumEntries()-(dPl_pPb->sumEntries()+dPl_Pbp->sumEntries()))>0.05*(dPl_pPb->sumEntries()+dPl_Pbp->sumEntries())) {
    std::cout << "[ERROR] Number of weighted entries for the combined Pl " << dsSampleTag << " dataset ( " 
              << dPl_PA->sumEntries() << " , " << dPl_pPb->sumEntries() << " , " << dPl_Pbp->sumEntries() << " ) is inconsistent!" << std::endl; return false;
  }
  if (std::abs(dMi_PA->sumEntries()-(dMi_pPb->sumEntries()+dMi_Pbp->sumEntries()))>0.05*(dMi_pPb->sumEntries()+dMi_Pbp->sumEntries())) {
    std::cout << "[ERROR] Number of weighted entries for the combined Mi " << dsSampleTag << " dataset ( "
              << dMi_PA->sumEntries() << " , " << dMi_pPb->sumEntries() << " , " << dMi_Pbp->sumEntries() << " ) is inconsistent!" << std::endl; return false;
  }
  //
  // Import the new PA dataset
  Workspaces[sample_PA].import(*dPl_PA);
  Workspaces.at(sample_PA).import(*dMi_PA);
  if ( mcPl_PA!=NULL) { Workspaces.at(sample_PA).import( *mcPl_PA); }
  if ( mcMi_PA!=NULL) { Workspaces.at(sample_PA).import( *mcMi_PA); }
  if (metPl_PA!=NULL) { Workspaces.at(sample_PA).import(*metPl_PA); }
  if (metMi_PA!=NULL) { Workspaces.at(sample_PA).import(*metMi_PA); }
  const std::vector< std::string > strLabels = { "METType" , "CorrectionApplied" , "RecoilMethod" };
  for (const auto& strL : strLabels) { if (Workspaces.at(sample_pPb).obj(strL.c_str())) { Workspaces.at(sample_PA).import(*Workspaces.at(sample_pPb).obj(strL.c_str()), strL.c_str()); } }
  // Import the Luminosity Ratio and xSections (used for systematics)
  if (sample_pPb.find("MC_")!=std::string::npos) {
    auto listVars = Workspaces.at(sample_pPb).allVars();
    std::unique_ptr<TIterator> parIt = std::unique_ptr<TIterator>(listVars.createIterator());
    for (auto it = (RooRealVar*)parIt->Next(); it!=NULL; it = (RooRealVar*)parIt->Next() ) {
      const std::string name = it->GetName();
      if ( (name=="rLumi") || (name.find("rXSection_")!=std::string::npos) ){
        if (Workspaces.at(sample_pPb).var(name.c_str())!=NULL) { Workspaces.at(sample_PA).import(*Workspaces.at(sample_pPb).var(name.c_str())); }
      }
    }
  }
  //
  // return
  return true;
};


#endif // #ifndef EWQForest2DataSet_C

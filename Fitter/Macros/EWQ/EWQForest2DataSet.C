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
#include "../../../Corrections/MET_Recoil/RecoilCorrector.C"
#include "../Utilities/initClasses.h"
#include <chrono>


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
  StringVector OutputFileDir   = FileInfo.at("OutputFileDir");
  StringVector InputFileNames  = FileInfo.at("InputFileNames");
  StringVector DSNames         = FileInfo.at("DSNames");
  const bool isData = (DSNames[0].find("DATA")!=std::string::npos);
  const bool isMC   = (DSNames[0].find("MC")!=std::string::npos);
  const bool useNoHFMET = (metTAG.find("NoHF")!=std::string::npos);
  for (auto& tag : DSNames) { 
    std::string o = (OutputFileDir[0] + chaDir + "/") + "DATASET_" + metTAG + tag + ".root"; 
    if (gSystem->AccessPathName(o.c_str())) { makeDir(OutputFileDir[1] + chaDir + "/"); o = (OutputFileDir[1] + chaDir + "/") + "DATASET_" + metTAG + tag + ".root"; }
    OutputFileNames.push_back(o);
  }
  // For computing time optimization
  auto t1_ALL = std::chrono::high_resolution_clock::now();
  auto t1 = t1_ALL;
  // Extract Input Information
  std::string  TYPE = info.Par.at("Analysis");
  int triggerIndex  = info.Int.at("triggerIndex");
  // Create RooDataSets
  std::vector< RooDataSet* > dataPl, dataMi, mcPl, mcMi;
  bool createDS = info.Flag.at("updateDS");
  // Check if RooDataSets exist and are not corrupt
  for (uint i=0; i<OutputFileNames.size(); i++) {
    if ( !gSystem->AccessPathName(OutputFileNames[i].c_str()) ) {
      std::cout << "[INFO] Loading RooDataSets from " << OutputFileNames[i] << std::endl;
      TFile *DBFile = TFile::Open(OutputFileNames[i].c_str(),"READ");
      if (!DBFile) { std::cout << "[ERROR] File: " << OutputFileNames[i] << " is corrupted!" << std::endl; return false; }
      dataPl.push_back( (RooDataSet*)DBFile->Get(Form("dPl_RAW_%s", DSNames[i].c_str())) );
      dataMi.push_back( (RooDataSet*)DBFile->Get(Form("dMi_RAW_%s", DSNames[i].c_str())) );
      if (isMC) { mcPl.push_back( (RooDataSet*)DBFile->Get(Form("mcPl_RAW_%s", DSNames[i].c_str())) ); }
      if (isMC) { mcMi.push_back( (RooDataSet*)DBFile->Get(Form("mcMi_RAW_%s", DSNames[i].c_str())) ); }
      if (checkEWQDS(dataPl[i], DSNames[i], TYPE)==false) { createDS = true; }
      if (checkEWQDS(dataMi[i], DSNames[i], TYPE)==false) { createDS = true; }
      DBFile->Close(); delete DBFile;
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
    if (useNoHFMET) { if (!metTree->GetTree(InputFileNames, "metAnaNoHF")) return false; }
    else { if (!metTree->GetTree(InputFileNames, "metAna")) return false; }
    if (metTree->GetEntries() != nentries) { std::cout << "[ERROR] Inconsistent number of entries!" << std::endl; return false; }
    std::unique_ptr<HiEvtTree> evtTree = std::unique_ptr<HiEvtTree>(new HiEvtTree());
    if (!evtTree->GetTree(InputFileNames)) return false;
    if (evtTree->GetEntries() != nentries) { std::cout << "[ERROR] Inconsistent number of entries!" << std::endl; return false; }
    ///// MC Generator Information
    RooRealVar mcNGen   = RooRealVar ( "NGen"      , "Number of Generated Events" ,   -1.0 ,   1000000000.0 , ""      );
    RooRealVar mcType   = RooRealVar ( "MC_Type"   , "MC Type Number"             , -100.0 ,          100.0 , ""      );
    RooRealVar bosonPt  = RooRealVar ( "Boson_Pt"  , "Boson p_{T}"                ,   -1.0 ,       100000.0 , "GeV/c" );
    RooRealVar bosonPhi = RooRealVar ( "Boson_Phi" , "Boson #phi"                 ,   -9.0 ,            9.0 , ""      );
    RooRealVar refPt    = RooRealVar ( "Ref_Pt"    , "Reference p_{T}"            ,   -1.0 ,       100000.0 , "GeV/c" );
    RooRealVar refPhi   = RooRealVar ( "Ref_Phi"   , "Reference #phi"             ,   -9.0 ,            9.0 , ""      );
    RooRealVar metPhi   = RooRealVar ( "MET_Phi"   , "#slash{E}_{T} #phi"         ,   -9.0 ,            9.0 , ""      );
    RooRealVar muPhi    = RooRealVar ( "Muon_Phi"  , "#mu #phi"                   ,   -9.0 ,            9.0 , ""      );
    RooArgSet  mcCols = RooArgSet(mcNGen, mcType, bosonPt, bosonPhi, refPt, refPhi);
    mcCols.add(metPhi); mcCols.add(muPhi);
    ///// RooDataSet Variables
    RooRealVar   met    = RooRealVar ( "MET"        , "|#slash{E}_{T}|"   ,  -1.0 , 100000.0 ,  "GeV/c"     );
    RooRealVar   muPt   = RooRealVar ( "Muon_Pt"    , "#mu p_{T}"         ,  -1.0 , 100000.0 ,  "GeV/c"     );
    RooRealVar   muEta  = RooRealVar ( "Muon_Eta"   , "#mu #eta"          , -10.0 ,     10.0 ,  ""          );
    RooRealVar   muIso  = RooRealVar ( "Muon_Iso"   , "#mu Isolation"     ,  -1.0 , 100000.0 ,  ""          );
    RooRealVar   muMT   = RooRealVar ( "Muon_MT"    , "W Transverse Mass" ,  -1.0 , 100000.0 ,  "GeV/c^{2}" );
    RooRealVar   cent   = RooRealVar ( "Centrality" , "Centrality"        ,  -1.0 , 100000.0 ,  ""          );
    RooCategory  type   = RooCategory( "Event_Type" , "Event Type");
    type.defineType("Other", -1); type.defineType("DYToMuMu", 1);
    RooArgSet cols = RooArgSet(met, muPt, muEta, muIso, muMT, cent);
    cols.add(type);
    ///// Initiliaze RooDataSets
    dataPl.clear(); dataMi.clear(); mcPl.clear(); mcMi.clear();
    for (uint i=0; i<DSNames.size(); i++) {
      std::cout << "[INFO] Creating " << "RooDataSet for " << DSNames[i] << std::endl;
      dataPl.push_back( new RooDataSet(Form("dPl_RAW_%s", DSNames[i].c_str()), "dPl", cols) );
      dataMi.push_back( new RooDataSet(Form("dMi_RAW_%s", DSNames[i].c_str()), "dMi", cols) );
      if (isMC) { mcPl.push_back( new RooDataSet(Form("mcPl_RAW_%s", DSNames[i].c_str()), "mcPl", mcCols ) ); }
      if (isMC) { mcMi.push_back( new RooDataSet(Form("mcMi_RAW_%s", DSNames[i].c_str()), "mcMi", mcCols ) ); }
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
      if (jentry%200000==0) std::cout << "[INFO] Processed : " << jentry << "/" << nentries << std::endl;
      if (jentry%200000==0) { 
        std::cout << "[INFO] Processing time: " << std::chrono::duration_cast<std::chrono::seconds>( std::chrono::high_resolution_clock::now() - t1 ).count() << " sec" << std::endl;
        t1 = std::chrono::high_resolution_clock::now();
      }
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
      if      (info.Par.at("VarType") == "PF_RAW"       ) { MET = metTree->PF_MET_NoShift_Mom();    }
      else if (info.Par.at("VarType") == "PF_Type1"     ) { MET = metTree->Type1_MET_NoShift_Mom(); }
      else if (info.Par.at("VarType") == "PF_NoHF_RAW"  ) { MET = metTree->PF_MET_NoShift_Mom();    }
      else if (info.Par.at("VarType") == "PF_NoHF_Type1") { MET = metTree->Type1_MET_NoShift_Mom(); }
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
        if (PA::passDrellYanVeto(muonTree) == false) { eventType = "DYToMuMu"; }  // Found a Drell-Yan candidate
      }
      if (eventType=="DYToMuMu") continue; // Remove Drell-Yan events
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
      }
      //
      //// Fill the RooDataSets
      for (uint i=0; i<DSNames.size(); i++) {
        if (DSNames[i].find(evtCol)!=std::string::npos) {
          if (leadMuChg > 0) { dataPl[i]->addFast(cols); }
          if (leadMuChg < 0) { dataMi[i]->addFast(cols); }
          if (isMC && (leadMuChg > 0)) { mcPl[i]->addFast(mcCols); }
          if (isMC && (leadMuChg < 0)) { mcMi[i]->addFast(mcCols); }
        }
      }
    }
    //// Save the RooDataSets
    for (uint i=0; i<DSNames.size(); i++) {
      TFile *DBFile = TFile::Open(OutputFileNames[i].c_str(),"RECREATE");
      DBFile->cd();
      dataPl[i]->Write(Form("dPl_RAW_%s", DSNames[i].c_str()));
      dataMi[i]->Write(Form("dMi_RAW_%s", DSNames[i].c_str()));
      if (isMC) { mcPl[i]->Write(Form("mcPl_RAW_%s", DSNames[i].c_str())); }
      if (isMC) { mcMi[i]->Write(Form("mcMi_RAW_%s", DSNames[i].c_str())); }
      DBFile->Write(); DBFile->Close(); delete DBFile;
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
      if (mcPl[i]) delete mcPl[i];
      if (mcMi[i]) delete mcMi[i];
    }
    TObjString tmp; tmp.SetString(info.Par.at("VarType").c_str()); Workspaces[DSNames[i]].import(*((TObject*)&tmp), "METType");
    // delete the local datasets
    delete dataPl[i];
    delete dataMi[i];
  }
  dataPl.clear(); dataMi.clear();
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


bool applyLumiWeight(RooWorkspace& ws, const std::string dDSName, const std::string& mcDSName, const bool verbose = false)
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
  std::string chg    = dDS->GetName(); chg    = chg.substr(1, sample.find("_"));
  std::string col    = dDS->GetName(); col    = col.substr(col.find_last_of("_")+1);
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
  auto dWDS  = new RooDataSet( Form("d%s_LUM_%s" , chg.c_str(), sample.c_str()) , Form("d%s" , chg.c_str()) , cols   , RooFit::WeightVar(weight) );
  auto mcWDS = new RooDataSet( Form("mc%s_LUM_%s", chg.c_str(), sample.c_str()) , Form("mc%s", chg.c_str()) , mcCols );
  // Apply the Lumi re-weighting to positive muon dataset
  int mcID = -99; double xSection = -1.0;
  for (int i = 0; i < dDS->numEntries(); i++) {
    dDS->get(i);
    mcDS->get(i);
    // Get the Luminosity
    mcLumi.setVal( (col=="pPb") ? PA::LUMI::Data_pPb : PA::LUMI::Data_Pbp );
    // Get the Cross-Section
    auto mcType = (RooRealVar*)mcCols.find("MC_Type");
    if (mcType==NULL) { std::cout << "[ERROR] applyLumiWeight: Variable MC_Type was not found in " << mcDSName << std::endl; return false; }
    if (mcID != mcType->getVal()) {
      if (!PA::getCrossSection(xSection, int(mcType->getVal()), verbose)) { return false; }
      mcID = mcType->getVal();
    }
    mcXSec.setVal( xSection );
    // Compute the Luminosity MC Weight
    auto mcNGen = (RooRealVar*)mcCols.find("NGen");
    if (mcNGen==NULL) { std::cout << "[ERROR] applyLumiWeight: Variable NGen was not found in " << mcDSName << std::endl; return false; }
    weight.setVal( (mcXSec.getVal() * mcLumi.getVal()) / mcNGen->getVal() );
    // Fill the new RooDataSets
    if (weight.getVal() <= 0.0) {
      std::cout << mcXSec.getVal() << "  " <<  mcLumi.getVal() << "  " << mcNGen->getVal() << std::endl;
      std::cout << "[ERROR] applyLumiWeight : Weight is negative ( " << weight.getVal() << " ) in " << mcDSName << std::endl; return false; }
    dWDS->addFast(cols, weight.getVal());
    mcWDS->addFast(mcCols);
  }
  // Import to RooWorkspace
  ws.import(*dWDS);
  ws.import(*mcWDS);
  // Clean the memory
  if (dWDS ) { delete dWDS;  }
  if (mcWDS) { delete mcWDS; }
  // Return
  return true;
  //
};


bool reweightMCLumi(RooWorkspaceMap_t& Workspaces)
{
  //
  for (auto& ws : Workspaces) {
    if (ws.first.find("MC_")!=std::string::npos) {
      RooWorkspace& myws = ws.second;
      // Get the sample name
      const auto& dsList =  myws.allData();
      std::string sample = dsList.front()->GetName(); sample = sample.substr(sample.find("_")+1);
      std::cout << "[INFO] Re-weighting the Luminosity of " << sample << std::endl;
      // Apply the Lumi re-weighting to positive muon dataset
      if (!applyLumiWeight(myws, ("dPl_"+sample), ("mcPl_"+sample), true)) { return false; }
      // Apply the Lumi re-weighting to negative muon dataset
      const bool printXSec = (sample.find("MC_W")!=std::string::npos);
      if (!applyLumiWeight(myws, ("dMi_"+sample), ("mcMi_"+sample), printXSec)) { return false; }
    }
  }
  //
  return true;
  //
};


bool applyMCCorrection(RooWorkspace& ws, const std::string dDSName, const std::string& mcDSName, const bool& applyTnPCorr, RecoilCorrector& recoilCorr, const std::string& recoilMethod)
{
  //
  const bool applyRecoilCorr = recoilCorr.isValid();
  if (!applyTnPCorr && !applyRecoilCorr) { return true; }
  if (dDSName.find("MC_")==std::string::npos) { std::cout << "[ERROR] applyMCCorrection: Dataset " << dDSName << " is not from MC!" << std::endl; return false; }
  // Get the RooDatasets
  auto  dDS = (RooDataSet*)ws.data(dDSName.c_str());
  if ( dDS==NULL) { std::cout << "[ERROR] applyMCCorrection: Dataset " << dDSName << " was not found in the workspace!" << std::endl; return false; }
  auto mcDS = (RooDataSet*)ws.data(mcDSName.c_str());
  if (mcDS==NULL) { std::cout << "[ERROR] applyMCCorrection: Dataset " << mcDSName << " was not found in the workspace!" << std::endl; return false; }
  // Get the sample name
  std::string sample = dDS->GetName(); sample = sample.substr(sample.find("MC_"));
  std::string chg    = dDS->GetName(); chg    = chg.substr(1, sample.find("_"));
  std::string col    = dDS->GetName(); col    = col.substr(col.find_last_of("_")+1);
  // Initialize the new RooRealVars
  RooRealVar mcRawMET = RooRealVar( "MET"     , "|#slash{E}_{T}|"          ,  -1.0 ,     100000.0 , "GeV/c"     );
  RooRealVar mcRawMT  = RooRealVar( "Muon_MT" , "W Transverse Mass"        ,  -1.0 ,     100000.0 , "GeV/c^{2}" );
  RooRealVar mcSFTnP  = RooRealVar( "SFTnP"   , "TagAndProbe Scale Factor" , -10.0 ,         10.0 , ""          );
  RooRealVar weight   = RooRealVar( "Weight"  , "Weight"                   ,  -1.0 , 1000000000.0 , ""          );
  // Initialize the RooArgSets
  RooArgSet cols   = *dDS->get();
  cols.add(weight);
  RooArgSet mcCols = *mcDS->get();
  if (applyTnPCorr) { mcCols.add(mcSFTnP); }
  if (applyRecoilCorr) { mcCols.add(mcRawMET); mcCols.add(mcRawMT); }
  //
  // Initialize the new RooDataSets
  auto dWDS  = new RooDataSet( Form("d%s_COR_%s" , chg.c_str(), sample.c_str()) , Form("d%s" , chg.c_str()) , cols   , RooFit::WeightVar(weight) );
  auto mcWDS = new RooDataSet( Form("mc%s_COR_%s", chg.c_str(), sample.c_str()) , Form("mc%s", chg.c_str()) , mcCols );
  // Apply the Lumi re-weighting to positive muon dataset
  int mcID = -99;
  for (int i = 0; i < dDS->numEntries(); i++) {
    dDS->get(i);
    mcDS->get(i);
    //
    // Apply the Tag And Probe corrections
    //
    weight.setVal( dDS->weight() );
    if (applyTnPCorr) {
      // Get the Muon Kinematic Variables
      auto muPt = (RooRealVar*)cols.find("Muon_Pt");
      if (muPt==NULL) { std::cout << "[ERROR] applyMCCorrection: Variable Muon_Pt was not found in " << dDSName << std::endl; return false; }
      auto muEta = (RooRealVar*)cols.find("Muon_Eta");
      if (muEta==NULL) { std::cout << "[ERROR] applyMCCorrection: Variable Muon_Eta was not found in " << dDSName << std::endl; return false; }
      // Get Nominal Tag and Probe Scale Factors
      const double sf_MuID = tnp_weight_muid_ppb( muPt->getVal() , muEta->getVal() , 0 );
      const double sf_Trig = tnp_weight_trg_ppb (                  muEta->getVal() , 0 );
      const double sf_Iso  = tnp_weight_iso_ppb ( muPt->getVal() , muEta->getVal() , 0 );
      mcSFTnP.setVal( sf_MuID * sf_Trig * sf_Iso );
      // Set the TnP Scale Factor weight
      weight.setVal( weight.getVal() * mcSFTnP.getVal() );
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
      // Set the Reference and Boson pT vectors
      recoilCorr.setPt(reference_pT, boson_pT);
      // Correct the MET
      TVector2 MET_CORR;
      if (!recoilCorr.correctMET(MET_CORR, MET_RAW, recoilMethod)) { return false; }
      // Correct the Muon Transverse Mass
      const double muMTVal = PA::getWTransverseMass(muPt->getVal(), muPhi->getVal(), MET_CORR.Mod(), MET_CORR.Phi());
      // Set the corrected MET and Transverse Mass
      met->setVal(MET_CORR.Mod());
      muMT->setVal(muMTVal);
    }
    // Fill the new RooDataSets
    if (weight.getVal() <= 0.0) { std::cout << "[ERROR] applyMCCorrection: Weight is negative ( " << weight.getVal() << " ) in " << dDSName << std::endl; return false; }
    dWDS->addFast(cols, weight.getVal());
    mcWDS->addFast(mcCols);
  }
  // Import to RooWorkspace
  ws.import(*dWDS);
  ws.import(*mcWDS);
  // Clean the memory
  if (dWDS ) { delete dWDS;  }
  if (mcWDS) { delete mcWDS; }
  // Return
  return true;
};


bool correctMC(RooWorkspaceMap_t& Workspaces, const GlobalInfo& info)
{
  //
  const bool applyTnPCorr    = info.Flag.at("applyTnPCorr");
  const bool applyRecoilCorr = info.Flag.at("applyRecoilCorr");
  if (!applyTnPCorr && !applyRecoilCorr) { std::cout << "[INFO] Corrections disabled, no MC corrections will be applied!" << std::endl; return true; }
  //
  // Define the MET Recoil Corrector
  RecoilCorrector recoilCorr = RecoilCorrector();
  std::string recoilMethod;
  if (applyRecoilCorr) {
    recoilMethod = info.Par.at("RecoilCorrMethod");
    const std::string met = info.Par.at("METType");
    std::string preCWD = getcwd(NULL, 0); preCWD.erase(preCWD.find_last_of("/"), 100);
    const std::string recoilDir     = Form("%s/Corrections/MET_Recoil", preCWD.c_str());
    const std::string fileName_MC   = Form("%s/FitRecoil/MC_DYToMuMu_PYQUEN/MET_%s/PA/Results/fits_RecoilPDF_%s_PA.root", recoilDir.c_str(), met.c_str(), met.c_str());
    const std::string fileName_DATA = Form("%s/FitRecoil/DATA/MET_%s/PA/Results/fits_RecoilPDF_%s_PA.root", recoilDir.c_str(), met.c_str(), met.c_str());
    if (!recoilCorr.setInputFiles(met, fileName_MC, fileName_DATA)) { return false; }
  }
  //
  for (auto& ws : Workspaces) {
    if (ws.first.find("MC_")!=std::string::npos) {
      RooWorkspace& myws = ws.second;
      // Get the sample name
      std::string sample = "LUM_" + ws.first;
      if (myws.data(Form("dPl_%s", sample.c_str()))==NULL) { sample = "RAW_" + ws.first; }
      std::cout << "[INFO] Applying corrections ( " << (applyTnPCorr ? "TnPCorr ," : "") << (applyRecoilCorr ? " RecoilCorr" : "") << " ) to " << sample << std::endl;
      // Initialize the MET Recoil Corrector
      if (applyRecoilCorr) { recoilCorr.setInitialSetup(sample); }
      // Apply the MC correction to positive muon dataset
      if (!applyMCCorrection(myws, ("dPl_"+sample), ("mcPl_"+sample), applyTnPCorr, recoilCorr, recoilMethod)) { return false; }
      // Apply the MC correction to negative muon dataset
      if (!applyMCCorrection(myws, ("dMi_"+sample), ("mcMi_"+sample), applyTnPCorr, recoilCorr, recoilMethod)) { return false; }
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
    auto pt  = (RooRealVar*)set->find("Muon_Pt");
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
  if (Workspaces.at(sample_pPb).data(("dPl_"+dsType+"_"+sample_pPb).c_str())==NULL) { dsType = "RAW"; }
  const std::string dsSampleTag  = dsType + "_" + sampleTag;
  const std::string dsSample_pPb = dsSampleTag + "_pPb";
  const std::string dsSample_Pbp = dsSampleTag + "_Pbp";
  const std::string dsSample_PA  = dsSampleTag + "_PA";
  //
  auto dPl_pPb = (RooDataSet*)Workspaces.at(sample_pPb).data(("dPl_"+dsSample_pPb).c_str());
  if (dPl_pPb==NULL) { std::cout << "[ERROR] RooDataSet " << ("dPl_"+dsSample_pPb) << " does not exist!" << std::endl; return false; }
  auto dMi_pPb = (RooDataSet*)Workspaces.at(sample_pPb).data(("dMi_"+dsSample_pPb).c_str());
  if (dMi_pPb==NULL) { std::cout << "[ERROR] RooDataSet " << ("dMi_"+dsSample_pPb) << " does not exist!" << std::endl; return false; }
  auto mcPl_pPb = (RooDataSet*)Workspaces.at(sample_pPb).data(("mcPl_"+dsSample_pPb).c_str());
  auto mcMi_pPb = (RooDataSet*)Workspaces.at(sample_pPb).data(("mcMi_"+dsSample_pPb).c_str());
  //
  auto dPl_Pbp = (RooDataSet*)Workspaces.at(sample_Pbp).data(("dPl_"+dsSample_Pbp).c_str());
  if (dPl_Pbp==NULL) { std::cout << "[ERROR] RooDataSet " << ("dPl_"+dsSample_Pbp) << " does not exist!" << std::endl; return false; }
  auto dMi_Pbp = (RooDataSet*)Workspaces.at(sample_Pbp).data(("dMi_"+dsSample_Pbp).c_str());
  if (dMi_Pbp==NULL) { std::cout << "[ERROR] RooDataSet " << ("dMi_"+dsSample_Pbp) << " does not exist!" << std::endl; return false; }
  auto mcPl_Pbp = (RooDataSet*)Workspaces.at(sample_Pbp).data(("mcPl_"+dsSample_Pbp).c_str());
  auto mcMi_Pbp = (RooDataSet*)Workspaces.at(sample_Pbp).data(("mcMi_"+dsSample_Pbp).c_str());
  //
  std::cout << "[INFO] Creating combined PA RooDataSet for " << dsSampleTag << std::endl;
  //
  // Copy the pPb datasets
  auto dPl_PA = (RooDataSet*)dPl_pPb->Clone(("dPl_"+dsSample_PA).c_str());
  if (dPl_PA==NULL || dPl_PA->sumEntries()==0) { std::cout << "[ERROR] RooDataSet " << ("dPl_"+dsSample_PA) << " was not created!" << std::endl; return false; }
  auto dMi_PA = (RooDataSet*)dMi_pPb->Clone(("dMi_"+dsSample_PA).c_str());
  if (dMi_PA==NULL || dMi_PA->sumEntries()==0) { std::cout << "[ERROR] RooDataSet " << ("dMi_"+dsSample_PA) << " was not created!" << std::endl; return false; }
  RooDataSet* mcPl_PA = NULL; if (mcPl_pPb!=NULL) { mcPl_PA = (RooDataSet*)mcPl_pPb->Clone(("mcPl_"+dsSample_PA).c_str()); }
  RooDataSet* mcMi_PA = NULL; if (mcMi_pPb!=NULL) { mcMi_PA = (RooDataSet*)mcPl_pPb->Clone(("mcMi_"+dsSample_PA).c_str()); }
  //
  // Invert the eta of Pbp dataset and fill the PA dataset
  if (!invertEtaAndFill(*dPl_PA, *dPl_Pbp)) { return false; }
  if (!invertEtaAndFill(*dMi_PA, *dMi_Pbp)) { return false; }
  if (mcPl_pPb!=NULL && mcPl_PA!=NULL) { mcPl_PA->append(*mcPl_Pbp); }
  if (mcMi_pPb!=NULL && mcMi_PA!=NULL) { mcMi_PA->append(*mcMi_Pbp); }
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
  if (mcPl_PA!=NULL) { Workspaces.at(sample_PA).import(*mcPl_PA); }
  if (mcMi_PA!=NULL) { Workspaces.at(sample_PA).import(*mcMi_PA); }
  Workspaces.at(sample_PA).import(*((TObjString*)Workspaces.at(sample_pPb).obj("METType")), "METType");
  //
  // Clean up the memory
  if (dPl_PA) delete dPl_PA;
  if (dMi_PA) delete dMi_PA;
  if (mcPl_PA) delete mcPl_PA;
  if (mcMi_PA) delete mcMi_PA;
  //
  // return
  return true;
};


#endif // #ifndef EWQForest2DataSet_C

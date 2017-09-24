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
#include "../../../Utilities/RecoilCorrector.C"
#include "../Utilities/initClasses.h"
#include <chrono>


bool checkEWQDS          ( RooDataSet* DS              , const std::string& DSName       , const std::string& Analysis );
bool EWQForest_WToMuNu   ( RooWorkspaceMap& Workspaces , const StringVectorMap& FileInfo , const GlobalInfo&  info     );


bool EWQForest2DataSet(RooWorkspaceMap& Workspaces, const StringVectorMap& FileInfo, const GlobalInfo& info)
{
  std::string Analysis = info.Par.at("Analysis");
  if (Analysis=="WToMuNu") { if (!EWQForest_WToMuNu(Workspaces, FileInfo, info)) return false; }
  return true;
};

bool EWQForest_WToMuNu(RooWorkspaceMap& Workspaces, const StringVectorMap& FileInfo, const GlobalInfo& info)
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
  bool applyCorr  = info.Flag.at("applyCorr");
  // Create RooDataSets
  std::vector< RooDataSet* > dataPl, dataMi, mcPl, mcMi;
  bool createDS = info.Flag.at("updateDS");
  // Check if RooDataSets exist and are not corrupt
  for (uint i=0; i<OutputFileNames.size(); i++) {
    if ( !gSystem->AccessPathName(OutputFileNames[i].c_str()) ) {
      std::cout << "[INFO] Loading RooDataSets from " << OutputFileNames[i] << std::endl;
      TFile *DBFile = TFile::Open(OutputFileNames[i].c_str(),"READ");
      if (!DBFile) { std::cout << "[ERROR] File: " << OutputFileNames[i] << " is corrupted!" << std::endl; return false; }
      dataPl.push_back( (RooDataSet*)DBFile->Get(Form("dPl_%s", DSNames[i].c_str())) );
      dataMi.push_back( (RooDataSet*)DBFile->Get(Form("dMi_%s", DSNames[i].c_str())) );
      if (isMC) { mcPl.push_back( (RooDataSet*)DBFile->Get(Form("mcPl_%s", DSNames[i].c_str())) ); }
      if (isMC) { mcMi.push_back( (RooDataSet*)DBFile->Get(Form("mcMi_%s", DSNames[i].c_str())) ); }
      if (checkEWQDS(dataPl[i], DSNames[i], TYPE)==false) { createDS = true; }
      if (checkEWQDS(dataMi[i], DSNames[i], TYPE)==false) { createDS = true; }
      DBFile->Close(); delete DBFile;
    }
    else { createDS = true; break; }
  }
  if (createDS) {
    ///// Input Forest
    //
    std::vector< std::pair< std::string , double > > inFileNames;
    for (uint i = 0; i < InputFileNames.size(); i++) {
      double xSection = 1.0; if (isMC && !PA::getCrossSection(xSection, FileInfo.at("XSectionTags")[i])) { return false; }
      inFileNames.push_back(std::make_pair(InputFileNames[i] , xSection));
    }
    // Get the Recoil Corrections
    //
    RecoilCorrector recoilCorr = RecoilCorrector();
    const std::string met = info.Par.at("VarType");
    const std::string recoilDir = "/home/llr/cms/stahl/ElectroWeakAnalysis/EWQAnalysis2017/Corrections/MET_Recoil";
    const std::string fileName_MC   = Form("%s/FitRecoil/MC_DYToMuMu_PYQUEN/MET_%s/PA/Results/fits_RecoilPDF_%s_PA.root", recoilDir.c_str(), met.c_str(), met.c_str());
    const std::string fileName_DATA = Form("%s/FitRecoil/DATA/MET_%s/PA/Results/fits_RecoilPDF_%s_PA.root", recoilDir.c_str(), met.c_str(), met.c_str());
    recoilCorr.setInputFiles(met, fileName_MC, fileName_DATA);
    if (isMC) { recoilCorr.setInitialSetup(DSNames[0]); }
    //
    std::unique_ptr<HiMuonTree> muonTree = std::unique_ptr<HiMuonTree>(new HiMuonTree());
    if (!muonTree->GetTree(inFileNames)) return false;
    Long64_t nentries = muonTree->GetEntries();
    std::unique_ptr<HiMETTree> metTree = std::unique_ptr<HiMETTree>(new HiMETTree());
    if (useNoHFMET) { if (!metTree->GetTree(inFileNames, "metAnaNoHF")) return false; }
    else { if (!metTree->GetTree(inFileNames, "metAna")) return false; }
    if (metTree->GetEntries() != nentries) { std::cout << "[ERROR] Inconsistent number of entries!" << std::endl; return false; }
    std::unique_ptr<HiEvtTree> evtTree = std::unique_ptr<HiEvtTree>(new HiEvtTree());
    if (!evtTree->GetTree(inFileNames)) return false;
    if (evtTree->GetEntries() != nentries) { std::cout << "[ERROR] Inconsistent number of entries!" << std::endl; return false; }
    ///// MC Generator Information
    RooRealVar mcNGen  = RooRealVar ( "NGen",  "Number of Generated Events", -1.0, 100000000000.0,  ""        );
    RooRealVar mcXSec  = RooRealVar ( "XSec",  "Cross Section"             , -1.0, 100000000000.0,  "nb"      );
    RooRealVar mcLumi  = RooRealVar ( "Lumi",  "Luminosity"                , -1.0, 100000000000.0,  "nb^{-1}" );
    RooRealVar mcSFTnP = RooRealVar ( "SFTnP", "TagAndProbe Scale Factor"  , -10.0,          10.0,  ""        );
    RooArgSet  mcCols  = RooArgSet  (mcNGen, mcXSec, mcLumi, mcSFTnP);
    ///// RooDataSet Variables
    RooRealVar   met    = RooRealVar ( "MET",         "|#slash{E}_{T}|",   -1.0, 100000.0,  "GeV/c"     );
    RooRealVar   muPt   = RooRealVar ( "Muon_Pt",     "#mu p_{T}",         -1.0, 100000.0,  "GeV/c"     );
    RooRealVar   muEta  = RooRealVar ( "Muon_Eta",    "#mu #eta",          -10., 10.,       ""          );
    RooRealVar   muIso  = RooRealVar ( "Muon_Iso",    "#mu Isolation",     -1.0, 100000.0,  ""          );
    RooRealVar   muMT   = RooRealVar ( "Muon_MT",     "W Transverse Mass", -1.0, 100000.0,  "GeV/c^{2}" );
    RooRealVar   cent   = RooRealVar ( "Centrality",  "Centrality",        -1.0, 100000.0,  ""          );
    RooRealVar   weight = RooRealVar ( "Weight",      "Weight",            -1.0, 10000000000.0,  ""     );
    RooCategory  type   = RooCategory( "Event_Type",  "Event Type");
    type.defineType("Other", -1); type.defineType("DYToMuMu", 1);
    RooArgSet cols = RooArgSet(met, muPt, muEta, muIso, muMT, cent, weight);
    cols.add(type);
    ///// Initiliaze RooDataSets
    dataPl.clear(); dataMi.clear(); mcPl.clear(); mcMi.clear();
    for (uint i=0; i<DSNames.size(); i++) {
      std::cout << "[INFO] Creating " << "RooDataSet for " << DSNames[i] << std::endl;
      dataPl.push_back( new RooDataSet(Form("dPl_%s", DSNames[i].c_str()), "dPl", cols, RooFit::WeightVar(weight)) );
      dataMi.push_back( new RooDataSet(Form("dMi_%s", DSNames[i].c_str()), "dMi", cols, RooFit::WeightVar(weight)) );
      if (isMC) { mcPl.push_back( new RooDataSet(Form("mcPl_%s", DSNames[i].c_str()), "mcPl", mcCols ) ); }
      if (isMC) { mcMi.push_back( new RooDataSet(Form("mcMi_%s", DSNames[i].c_str()), "mcMi", mcCols ) ); }
    }
    ///// Iterate over the Input Forest
    int treeIdx = -1;
    std::cout << "[INFO] Starting to process " << nentries << " nentries" << std::endl;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
      //
      // Get the entry in the trees
      if (muonTree->GetEntry(jentry)<0) break;
      if (metTree ->GetEntry(jentry)<0) break;
      if (evtTree ->GetEntry(jentry)<0) break;
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
      // Get the Lumi re-weight for MC
      double lumiWeight = 1.0;
      if (isMC) {
        mcNGen.setVal ( muonTree->GetTreeEntries()  );
        mcXSec.setVal ( muonTree->GetCrossSection() );
        mcLumi.setVal ( (evtCol=="pPb") ? PA::LUMI::Data_pPb : PA::LUMI::Data_Pbp );
        lumiWeight = ( (mcXSec.getVal() * mcLumi.getVal()) / mcNGen.getVal() );
      }
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
      // Get Nominal Tag and Probe Scale Factors
      double sf_TnP  = 1.0;
      if (isMC) {
        const double sf_MuID = tnp_weight_muid_ppb( muP4.Pt() , muP4.Eta() , 0 );
        const double sf_Trig = tnp_weight_trg_ppb (             muP4.Eta() , 0 );
        const double sf_Iso  = tnp_weight_iso_ppb ( muP4.Pt() , muP4.Eta() , 0 );
        sf_TnP  = ( sf_MuID * sf_Trig * sf_Iso );
        mcSFTnP.setVal( sf_TnP );
      }
      // Apply the MET Recoil corrections
      //
      // Event Corrections
      const double evtCorr = sf_TnP;
      //
      // Set Event Weight
      const double evtWeight = ( lumiWeight * (applyCorr ? evtCorr : 1.0) );
      //
      // Store isolation and charge of leading muon
      const float leadMuIso = muonTree->PF_Muon_IsoPFR03NoPUCorr()[leadMuPFIdx];
      const int   leadMuChg = int(muonTree->PF_Muon_Charge()[leadMuPFIdx]);
      if (isMC) { mcMuChg.setVal( leadMuChg ); }
      //
      // Recompute the Transver Mass based on the chosen MET
      TLorentzVector pfMuonP4T = TLorentzVector(), METP4 = TLorentzVector();
      pfMuonP4T.SetPtEtaPhiM(muP4.Pt(), 0.0, muP4.Phi(), muP4.M());
      METP4.SetPtEtaPhiM( MET.Mod(), 0.0, MET.Phi(), 0.0 );
      TLorentzVector muT = TLorentzVector( pfMuonP4T + METP4 );
      //
      // Classify the events
      //
      std::string eventType = "Other";
      //
      // Event Type: DrellYan->MuMu
      if (eventType=="Other") {
        if (PA::passDrellYanVeto(muonTree) == false) { eventType = "DYToMuMu"; }  // Found a Drell-Yan candidate
      }
      //
      //// Set the variables
      met.setVal    ( MET.Mod()  );
      muPt.setVal   ( muP4.Pt()  );
      muEta.setVal  ( muP4.Eta() );
      muIso.setVal  ( leadMuIso  );
      muMT.setVal   ( muT.M()    );
      cent.setVal   ( 0.0        );
      weight.setVal ( evtWeight  );
      type.setLabel ( eventType.c_str() );
      //
      //// Fill the RooDataSets
      for (uint i=0; i<DSNames.size(); i++) {
        if (DSNames[i].find(evtCol)!=std::string::npos) {
          if (leadMuChg > 0) { dataPl[i]->add(cols, weight.getVal()); }
          if (leadMuChg < 0) { dataMi[i]->add(cols, weight.getVal()); }
          if (isMC && (leadMuChg > 0)) { mcPl[i]->add(mcCols); }
          if (isMC && (leadMuChg < 0)) { mcMi[i]->add(mcCols); }
        }
      }
    }
    //// Save the RooDataSets
    for (uint i=0; i<DSNames.size(); i++) {
      TFile *DBFile = TFile::Open(OutputFileNames[i].c_str(),"RECREATE");
      DBFile->cd();
      dataPl[i]->Write(Form("dPl_%s", DSNames[i].c_str()));
      dataMi[i]->Write(Form("dMi_%s", DSNames[i].c_str()));
      if (isMC) { mcPl[i]->Write(Form("mcPl_%s", DSNames[i].c_str())); }
      if (isMC) { mcMi[i]->Write(Form("mcMi_%s", DSNames[i].c_str())); }
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


#endif // #ifndef EWQForest2DataSet_C

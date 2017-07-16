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


#include "../Utilities/HiMETTree.h"
#include "../Utilities/HiMuonTree.h"
#include "../Utilities/HiEvtTree.h"
#include "../Utilities/initClasses.h"


bool checkEWQDS          ( RooDataSet* DS              , const std::string& DSName       , const std::string& Analysis );
bool EWQForest_WToMuNu   ( RooWorkspaceMap& Workspaces    , const StringVectorMap& FileInfo , const GlobalInfo&  info  );


bool EWQForest2DataSet(RooWorkspaceMap& Workspaces, const StringVectorMap& FileInfo, const GlobalInfo& info)
{
  std::string Analysis = info.Par.at("Analysis");
  if (Analysis=="WToMuNu") { if (!EWQForest_WToMuNu(Workspaces, FileInfo, info)) return false; }
  return true;
};

bool EWQForest_WToMuNu(RooWorkspaceMap& Workspaces, const StringVectorMap& FileInfo, const GlobalInfo& info)
{
  StringVector OutputFileNames;
  std::string  chaDir = "Muon";
  std::string  metTAG = "MET" + info.Par.at("VarType") + "_";
  StringVector OutputFileDir   = FileInfo.at("OutputFileDir");
  StringVector InputFileNames  = FileInfo.at("InputFileNames");
  StringVector DSNames         = FileInfo.at("DSNames");
  bool isData = (DSNames[0].find("DATA")!=std::string::npos);
  bool useNoHFMET = (metTAG.find("NoHF")!=std::string::npos);
  for (auto& tag : DSNames) { 
    std::string o = (OutputFileDir[0] + chaDir + "/") + "DATASET_" + metTAG + tag + ".root"; 
    if (gSystem->AccessPathName(o.c_str())) { makeDir(OutputFileDir[1] + chaDir + "/"); o = (OutputFileDir[1] + chaDir + "/") + "DATASET_" + metTAG + tag + ".root"; }
    OutputFileNames.push_back(o);
  }
  // Extract Input Information
  std::string  TYPE = info.Par.at("Analysis");
  int triggerIndex  = info.Int.at("triggerIndex");
  bool applyWeight  = info.Flag.at("applyWeight");
  bool isMC = (DSNames[0].find("MC")!=std::string::npos);
  ///// Number of Generated Events
  std::vector< RooRealVar > NGen;
  // Create RooDataSets
  std::vector< RooDataSet* > dataPl, dataMi; RooDataSet* dataMC = 0;
  bool createDS = info.Flag.at("updateDS");
  // Check if RooDataSets exist and are not corrupt
  for (uint i=0; i<OutputFileNames.size(); i++) {
    if ( !gSystem->AccessPathName(OutputFileNames[i].c_str()) ) {
      cout << "[INFO] Loading RooDataSets from " << OutputFileNames[i] << endl;
      TFile *DBFile = TFile::Open(OutputFileNames[i].c_str(),"READ");
      if (!DBFile) { std::cout << "[ERROR] File: " << OutputFileNames[i] << " is corrupted!" << std::endl; return false; }
      dataPl.push_back( (RooDataSet*)DBFile->Get(Form("dPl_%s", DSNames[i].c_str())) );
      dataMi.push_back( (RooDataSet*)DBFile->Get(Form("dMi_%s", DSNames[i].c_str())) );
      if (checkEWQDS(dataPl[i], DSNames[i], TYPE)==false) { createDS = true; }
      if (checkEWQDS(dataMi[i], DSNames[i], TYPE)==false) { createDS = true; }
      if (isMC && !createDS) dataMC = (RooDataSet*)DBFile->Get(Form("dMC_%s", DSNames[i].c_str()));
      if (isMC && !createDS) NGen.push_back( *((RooRealVar*)DBFile->Get(Form("NGen_%s", DSNames[i].c_str()))) );
      DBFile->Close(); delete DBFile;
    }
    else { createDS = true; break; }
  }
  if (createDS) {
    ///// Input Forest
    std::unique_ptr<HiMuonTree> muonTree = std::unique_ptr<HiMuonTree>(new HiMuonTree());
    if (!muonTree->GetTree(InputFileNames[0])) return false;
    Long64_t nentries = muonTree->GetEntries();
    std::unique_ptr<HiMETTree> metTree = std::unique_ptr<HiMETTree>(new HiMETTree());
    if (useNoHFMET) { if (!metTree->GetTree(InputFileNames[0], 0, "metAnaNoHF")) return false; }
    else { if (!metTree->GetTree(InputFileNames[0], 0, "metAna")) return false; }
    if (metTree->GetEntries() != nentries) { std::cout << "[ERROR] Inconsistent number of entries!" << std::endl; return false; }
    muonTree->Tree()->AddFriend(metTree->Tree());
    std::unique_ptr<HiEvtTree> evtTree = std::unique_ptr<HiEvtTree>(new HiEvtTree());
    if (!evtTree->GetTree(InputFileNames[0])) return false;
    if (evtTree->GetEntries() != nentries) { std::cout << "[ERROR] Inconsistent number of entries!" << std::endl; return false; }
    muonTree->Tree()->AddFriend(evtTree->Tree());
    ///// RooDataSet Variables
    RooRealVar   met    = RooRealVar ( "MET",         "|#slash{E}_{T}|",   -1.0, 100000.0,  "GeV/c"     );
    RooRealVar   muPt   = RooRealVar ( "Muon_Pt",     "#mu p_{T}",         -1.0, 100000.0,  "GeV/c"     );
    RooRealVar   muEta  = RooRealVar ( "Muon_Eta",    "#mu #eta",          -10., 10.,       ""          );
    RooRealVar   muIso  = RooRealVar ( "Muon_Iso",    "#mu Isolation",     -1.0, 100000.0,  ""          );
    RooRealVar   muMT   = RooRealVar ( "Muon_MT",     "W Transverse Mass", -1.0, 100000.0,  "GeV/c^{2}" );
    RooRealVar   cent   = RooRealVar ( "Centrality",  "Centrality",        -1.0, 100000.0,  ""          );
    RooRealVar   weight = RooRealVar ( "Weight",      "Weight",            -1.0, 10000000000.0,  ""     );
    RooCategory  type   = RooCategory( "Event_Type",  "Event Type");
    type.defineType("Other", -1); type.defineType("DYToMuMu", 1); type.defineType("ZToMuMu", 2);
    RooArgSet cols = RooArgSet(met, muPt, muEta, muIso, muMT, cent, weight);
    cols.add(type);
    // For GEN MC
    RooRealVar muChg    = RooRealVar( "Muon_Chg", "#mu charge", -2., 2., "" );
    RooArgSet  genCols  = RooArgSet(muPt, muEta, muChg);
    ///// Initiliaze RooDataSets
    dataPl.clear(); dataMi.clear(); NGen.clear();
    for (uint i=0; i<DSNames.size(); i++) {
      cout << "[INFO] Creating " << "RooDataSet for " << DSNames[i] << endl;
      dataPl.push_back( new RooDataSet(Form("dPl_%s", DSNames[i].c_str()), "dPl", cols, RooFit::WeightVar(weight)) );
      dataMi.push_back( new RooDataSet(Form("dMi_%s", DSNames[i].c_str()), "dMi", cols, RooFit::WeightVar(weight)) );
      if (isMC) {
        NGen.push_back( RooRealVar(Form("NGen_%s", DSNames[i].c_str()), "Number of GEN Events", -1.0, 10000000000.0, "") );
        dataMC = new RooDataSet(Form("dMC_%s", DSNames[i].c_str()), "dMC", genCols);
      }
    }
    ///// For Lumi Count
    std::map<UShort_t, bool> LS;
    ///// Iterate over the Input Forest
    cout << "[INFO] Starting to process " << nentries << " nentries" << endl;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
      if (muonTree->GetEntry(jentry)<0) break;
      metTree->SetEntry(jentry); evtTree->SetEntry(jentry);
      if (muonTree->Event_Run()!=metTree->Event_Run()       ) { std::cout << "[ERROR] MET Run does not agree!"     << std::endl; return false; }
      if (muonTree->Event_Number()!=metTree->Event_Number() ) { std::cout << "[ERROR] MET Event does not agree!"   << std::endl; return false; }
      if (muonTree->Event_Run()!=evtTree->run()             ) { std::cout << "[ERROR] HiEVT Run does not agree!"   << std::endl; return false; }
      if (muonTree->Event_Number()!=evtTree->evt()          ) { std::cout << "[ERROR] HiEVT Event does not agree!" << std::endl; return false; }
      if (jentry%1000000==0) cout << "[INFO] " << jentry << "/" << nentries << endl;
      ///// For GEN MC
      if (isMC) {
        // Find the Leading Pt Gen Muon
        float maxPt = -99. , maxRPt = -99.;  int maxIdx = -1 , maxRIdx = -1;
        for (uint imu=0; imu<muonTree->Gen_Muon_Mom().size(); imu++) {
          TLorentzVector p4 = muonTree->Gen_Muon_Mom()[imu];
          if ((muonTree->Gen_Muon_Reco_Idx()[imu]!=-1) && (maxPt < p4.Pt())) { maxRPt = p4.Pt(); maxRIdx = imu; }
          if (maxPt < p4.Pt()) { maxPt = p4.Pt(); maxIdx = imu; }
        }
        if (maxRIdx!=-1) maxIdx = maxRIdx;
        if (maxIdx==-1) { std::cout << "[ERROR] MC Event does not contain any generated muon!" << std::endl; return false; }
        TLorentzVector muGenP4 = muonTree->Gen_Muon_Mom()[maxIdx];
        Int_t muGenChgV = int(muonTree->Gen_Muon_Charge()[maxIdx]);
        muPt.setVal  ( muGenP4.Pt()  );
        muEta.setVal ( muGenP4.Eta() );
        muChg.setVal ( muGenChgV     );
        dataMC->add(genCols);
        LS[muonTree->Event_Lumi()] = 1;
      }
      // Apply Event Filters
      if (metTree->Flag_collisionEventSelectionPA()==false) continue;           // PA Event Selection
      if (metTree->Flag_collisionEventSelectionPA_rejectPU()==false) continue;  // Reject PU events (NEW)
      if (metTree->Flag_goodVertices()==false) continue;                        // Primary Vertex Filter (JetMET Recommended)
      if (metTree->Flag_globalTightHalo2016Filter()==false) continue;           // beam halo filter (JetMET Recommended)
      if (metTree->Flag_HBHENoiseFilter()==false) continue;                     //  HBHE noise filter (JetMET Recommended)  
      if (metTree->Flag_HBHENoiseIsoFilter()==false) continue;                  // HBHEiso noise filter (JetMET Recommended)  
      if (metTree->Flag_EcalDeadCellTriggerPrimitiveFilter()==false) continue;  //  ECAL TP filter (JetMET Recommended)  
      if (metTree->Flag_BadPFMuonFilter()==false) continue;                     // Bad PF Muon Filter (JetMET Recommended)
      if (metTree->Flag_BadChargedCandidateFilter()==false) continue;           //  Bad Charged Hadron Filter (JetMET Recommended)
      if (!isMC && metTree->Flag_eeBadScFilter()==false) continue;              //  ee badSC noise filter (JetMET Recommended only for DATA)
      if (metTree->Flag_duplicateMuons()==false) continue;                      // Remove duplicate muon events (Geovanny's filter)
      if (metTree->Flag_badMuons()==false) continue;                            // Remove bad muon events (Geovanny's filter)
      // Check Trigger Fired
      if (muonTree->Event_Trig_Fired()[triggerIndex]==false) continue;          // Trigger Event Fired
      // Find the Leading Pt Muon
      float maxPt = -99.; int maxIdx = -1;
      float minPt = 999999999.; int minIdx = -1;
      for (uint imu=0; imu<muonTree->PF_Muon_Mom().size(); imu++) {
        ushort imuR = muonTree->PF_Muon_Reco_Idx()[imu];
        if (muonTree->Pat_Muon_Trig().at(imuR)[triggerIndex]==false) continue;  // Trigger Muon Matching
        if (muonTree->Reco_Muon_isTight()[imuR]==false) continue;  // Only consider Tight Muons
        TLorentzVector p4 = muonTree->PF_Muon_Mom()[imu];
        if (abs(p4.Eta())>=2.4) continue;                          // Only consider Muons within the ETA acceptance
        if (maxPt < p4.Pt()) { maxPt = p4.Pt(); maxIdx = imu; }
        if (minPt > p4.Pt()) { minPt = p4.Pt(); minIdx = imu; }
      }
      if (maxIdx==-1) continue;                                    // Only consider events with at least one good muon
      // Apply PT Cut
      TLorentzVector muP4 = muonTree->PF_Muon_Mom()[maxIdx];
      if (muP4.Pt() < 25.) continue;                              // Consider only leading muons with pT larger or equal than 25 GeV
      // Classify the events
      //
      std::string eventType = "Other";
      // Event Type: DYZ->MuMu
      if (eventType=="Other") {
        for (uint iQQ=0; iQQ<muonTree->PF_DiMuon_Charge().size(); iQQ++) {
          ushort idx1 = muonTree->PF_DiMuon_Muon1_Idx()[iQQ];
          ushort idx2 = muonTree->PF_DiMuon_Muon2_Idx()[iQQ];
          ushort idx1R = muonTree->PF_Muon_Reco_Idx()[idx1];
          ushort idx2R = muonTree->PF_Muon_Reco_Idx()[idx2];
          TLorentzVector p4_Mu1 = muonTree->PF_Muon_Mom()[idx1];
          TLorentzVector p4_Mu2 = muonTree->PF_Muon_Mom()[idx2];
          if (
              (abs(p4_Mu1.Eta())<2.4 && abs(p4_Mu2.Eta())<2.4 ) &&   // Check that both muons are within the ETA acceptance
              (p4_Mu1.Pt()>15. && p4_Mu2.Pt()>15.             ) &&   // Check that both muons have pt larger than 15
              (muonTree->Reco_Muon_isTight()[idx1R] && muonTree->Reco_Muon_isTight()[idx2R]) &&  // Require both muons to pass Tight ID
              (muonTree->PF_Muon_IsoPFR03NoPUCorr()[idx1]<0.15 && muonTree->PF_Muon_IsoPFR03NoPUCorr()[idx2]<0.15)  // Select isolated muons
               )
            {
              eventType = "DYToMuMu";
              if (
                  (muonTree->PF_DiMuon_Charge()[iQQ] == 0     ) &&   // Select opposite sign dimuons
                  (abs(muonTree->PF_DiMuon_Mom()[iQQ].M()-90.)<20.)      // Select invariant mass within [70. , 110.] GeV/c^2
                  )
                {
                  eventType = "ZToMuMu";
                }
              break;  // Found a Drell-Yan / Z candidate
            }
        }
      }
      // Choose which Type of MET we want to use
      TVector2 MET = TVector2();
      if (info.Par.at("VarType") == "PFRaw"    ) { MET = metTree->PF_MET_NoShift_Mom();    }
      if (info.Par.at("VarType") == "PFType1"  ) { MET = metTree->Type1_MET_NoShift_Mom(); }
      if (info.Par.at("VarType") == "NoHFRaw"  ) { MET = metTree->PF_MET_NoShift_Mom();    }
      if (info.Par.at("VarType") == "NoHFType1") { MET = metTree->Type1_MET_NoShift_Mom(); }
      // Recompute the Transver Mass based on the chosen MET
      TLorentzVector pfMuonP4T = TLorentzVector(), METP4 = TLorentzVector();
      pfMuonP4T.SetPtEtaPhiM(muP4.Pt(), 0.0, muP4.Phi(), muP4.M());
      METP4.SetPtEtaPhiM( MET.Mod(), 0.0, MET.Phi(), 0.0 );
      TLorentzVector muT = TLorentzVector( pfMuonP4T + METP4 );

      // Get the PV-Z correction
      double w_zPV = 1.0;
      if (isMC) { 
        TF1 fWeight = TF1("fWeight","gaus(0)/(gaus(3))", -30., 30.);
        fWeight.SetParameters(0.0207, 1.5839, 4.8070, 0.0176, 1.5073, 5.6747);
        w_zPV = fWeight.Eval(muonTree->Event_PriVtx_Pos().Z());
      }
      double w_Evt = w_zPV;
      
      //// Set the variables
      uint runNumber = muonTree->Event_Run();
      Float_t        muIsoV = muonTree->PF_Muon_IsoPFR03NoPUCorr()[maxIdx];
      Int_t          muChgV = int(muonTree->PF_Muon_Charge()[maxIdx]);
      met.setVal    ( MET.Mod()  );
      muPt.setVal   ( muP4.Pt()  );
      muEta.setVal  ( muP4.Eta() );
      muIso.setVal  ( muIsoV     );
      muMT.setVal   ( muT.M()    );
      cent.setVal   ( 0.0        );
      weight.setVal ( w_Evt      );
      type.setLabel ( eventType.c_str() );
      //// Fill the RooDataSets
      for (uint i=0; i<DSNames.size(); i++) {
        bool save = false;
        if (isMC) { 
          save = true; 
        }
        else {
          if ( (DSNames[i].find("pPb")!=std::string::npos) && (runNumber>=285952 && runNumber<=286504) ) save = true;
          if ( (DSNames[i].find("Pbp")!=std::string::npos) && (runNumber>=285410 && runNumber<=285951) ) save = true;
          if ( (DSNames[i].find("PA")!=std::string::npos)  && (runNumber>=285410 && runNumber<=286504) ) save = true;
        }
        if (save && (muChgV == +1)) dataPl[i]->add(cols, ((applyWeight) ? weight.getVal() : 1.0));
        if (save && (muChgV == -1)) dataMi[i]->add(cols, ((applyWeight) ? weight.getVal() : 1.0));
      }
    }
    //// Save the RooDataSets
    for (uint i=0; i<DSNames.size(); i++) {
      TFile *DBFile = TFile::Open(OutputFileNames[i].c_str(),"RECREATE");
      DBFile->cd();
      dataPl[i]->Write(Form("dPl_%s", DSNames[i].c_str()));
      dataMi[i]->Write(Form("dMi_%s", DSNames[i].c_str()));
      if (isMC) {
        Long64_t genEntries = nentries;
        double diff = (abs((LS.size()*600000) - (nentries/3.779e-4))/(nentries/3.779e-4));
        if (DSNames[i].find("QCD")!=std::string::npos) { if (diff<0.01) { genEntries = LS.size()*600000; } else { genEntries = (nentries/3.779e-4); } }
        NGen[i].setVal(genEntries); NGen[i].Write(Form("NGen_%s", DSNames[i].c_str()));
      }
      DBFile->Write(); DBFile->Close(); delete DBFile;
    }
  }
  // Import datasets to the workspaces
  for (uint i=0; i<DSNames.size(); i++) {
    if(!dataPl[i]) { cout << "[ERROR] " << DSNames[i] << " plus dataset was not found" << endl;  return false; }
    if(!dataMi[i]) { cout << "[ERROR] " << DSNames[i] << " minus dataset was not found" << endl; return false; }
    if(dataPl[i]->numEntries()==0) { cout << "[WARNING] " << DSNames[i] << " plus dataset is empty!" << endl;  return false; }
    if(dataMi[i]->numEntries()==0) { cout << "[WARNING] " << DSNames[i] << " minus dataset is empty!" << endl; return false; }
    Workspaces[DSNames[i]].import(*dataPl[i]);
    Workspaces[DSNames[i]].import(*dataMi[i]);
    if (isMC) {
      Workspaces[DSNames[i]].import(NGen[i]);
      if(dataMC && dataMC->numEntries()>0) Workspaces[DSNames[i]].import(*dataMC);
      if (dataMC) delete dataMC;
    }
    TObjString tmp; tmp.SetString(info.Par.at("VarType").c_str()); Workspaces[DSNames[i]].import(*((TObject*)&tmp), "METType");
    // delete the local datasets
    delete dataPl[i];
    delete dataMi[i];
  }
  dataPl.clear(); dataMi.clear(); NGen.clear();
  return true;
};

bool checkEWQDS(RooDataSet* DS, const string& DSName, const std::string& Analysis)
{
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
      { cout << "[WARNING] Original dataset: " << DS->GetName() << " is corrupted, will remake it!" << endl; }
  }
  return false;
};


#endif // #ifndef EWQForest2DataSet_C

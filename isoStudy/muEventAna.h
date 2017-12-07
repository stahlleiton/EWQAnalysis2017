//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Feb 22 19:01:46 2017 by ROOT version 6.08/04
// from TTree Muon_Event/
// found on file: root://eoscms//eos/cms/store/group/phys_heavyions/anstahll/EWQAnalysis2017/pPb2016/8160GeV/Data/HiEWQForest_PASingleMuon_pPb_Pbp_8160GeV_20170213.root
//////////////////////////////////////////////////////////

#ifndef muEventAna_h
#define muEventAna_h

// include the classes for the two other trees
#include "muPFAna.C"
#include "muRecoAna.C"
#include "METPFAna.C"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>

#include <map>
#include <iostream>

using namespace std;

// Header file for the classes stored in the TTree if any.
#include "TVector3.h"
#include "vector"
#include "vector"

class muEventAna {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // friend trees
   muPFAna   *fFriend1;
   muRecoAna *fFriend2;
   METPFAna  *fFriend3;

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          Event_Run;
   UShort_t        Event_Lumi;
   UInt_t          Event_Bx;
   ULong64_t       Event_Orbit;
   ULong64_t       Event_Number;
   UChar_t         Event_nPV;
   TVector3        *Event_PriVtx_Pos;
   TVector3        *Event_PriVtx_Err;
   vector<bool>    *Event_Trig_Fired;
   vector<int>     *Event_Trig_Presc;

   // trigger indices (from https://github.com/stahlleiton/cmssw/blob/EWQAnalysis2017/HeavyIonsAnalysis/ElectroWeakAnalysis/python/muonAnalyzer_cfi.py#L13):
   // triggerPathNames   = cms.vstring(
   //       "HLT_PAL2Mu12_v1",
   //       "HLT_PAL2Mu15_v1",
   //       "HLT_PAL3Mu3_v1",
   //       "HLT_PAL3Mu5_v3",
   //       "HLT_PAL3Mu7_v1",
   //       "HLT_PAL3Mu12_v1",
   //       "HLT_PAL3Mu15_v1"
   //       ),

   // List of branches
   TBranch        *b_Event_Run;   //!
   TBranch        *b_Event_Lumi;   //!
   TBranch        *b_Event_Bx;   //!
   TBranch        *b_Event_Orbit;   //!
   TBranch        *b_Event_Number;   //!
   TBranch        *b_Event_nPV;   //!
   TBranch        *b_Event_PriVtx_Pos;   //!
   TBranch        *b_Event_PriVtx_Err;   //!
   TBranch        *b_Event_Trig_Fired;   //!
   TBranch        *b_Event_Trig_Presc;   //!

   muEventAna(TTree *treeEvent=0,TTree *treePF=0,TTree *treeReco=0, TTree *treePFMET=0);
   muEventAna(TFile *filein);
   virtual ~muEventAna();
   virtual Int_t    Cut(TString var, int iqq, float cutval, bool isPFidx);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *treeEvent, TTree *treePF, TTree *treeReco, TTree *treeMETPF);
   virtual void     Loop(const char* filename);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   bool IsAccept(double pt, double eta);
   float Get(TString varname, int idx);
   float myPFIso(int idx, float deltaR);
  
};

#endif

#ifdef muEventAna_cxx
muEventAna::muEventAna(TFile *filein) {
   TTree *te = (TTree*) filein->Get("muonAna/Muon_Event");
   TTree *tp = (TTree*) filein->Get("muonAna/Muon_PF");
   TTree *tr = (TTree*) filein->Get("muonAna/Muon_Reco");
   TTree *tm = (TTree*) filein->Get("metAna/MET_PF");
   Init(te,tp,tr,tm);
}
muEventAna::muEventAna(TTree *treeEvent, TTree *treePF, TTree *treeReco, TTree *treeMETPF) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (treeEvent == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("root://eoscms//eos/cms/store/group/phys_heavyions/anstahll/EWQAnalysis2017/pPb2016/8160GeV/Data/HiEWQForest_PASingleMuon_pPb_Pbp_8160GeV_20170213.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("root://eoscms//eos/cms/store/group/phys_heavyions/anstahll/EWQAnalysis2017/pPb2016/8160GeV/Data/HiEWQForest_PASingleMuon_pPb_Pbp_8160GeV_20170213.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("root://eoscms//eos/cms/store/group/phys_heavyions/anstahll/EWQAnalysis2017/pPb2016/8160GeV/Data/HiEWQForest_PASingleMuon_pPb_Pbp_8160GeV_20170213.root:/muonAna");
      dir->GetObject("Muon_Event",treeEvent);

   }
   Init(treeEvent, treePF, treeReco, treeMETPF);
}

muEventAna::~muEventAna()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();

   if (fFriend1) delete fFriend1;
   if (fFriend2) delete fFriend2;
   if (fFriend3) delete fFriend3;
}

Int_t muEventAna::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t muEventAna::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void muEventAna::Init(TTree *treeEvent, TTree *treePF, TTree *treeReco, TTree *treeMETPF)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Event_PriVtx_Pos = 0;
   Event_PriVtx_Err = 0;
   Event_Trig_Fired = 0;
   Event_Trig_Presc = 0;
   // Set branch addresses and branch pointers
   if (!treeEvent) return;
   fChain = treeEvent;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Event_Run", &Event_Run, &b_Event_Run);
   fChain->SetBranchAddress("Event_Lumi", &Event_Lumi, &b_Event_Lumi);
   fChain->SetBranchAddress("Event_Bx", &Event_Bx, &b_Event_Bx);
   fChain->SetBranchAddress("Event_Orbit", &Event_Orbit, &b_Event_Orbit);
   fChain->SetBranchAddress("Event_Number", &Event_Number, &b_Event_Number);
   fChain->SetBranchAddress("Event_nPV", &Event_nPV, &b_Event_nPV);
   fChain->SetBranchAddress("Event_PriVtx_Pos", &Event_PriVtx_Pos, &b_Event_PriVtx_Pos);
   fChain->SetBranchAddress("Event_PriVtx_Err", &Event_PriVtx_Err, &b_Event_PriVtx_Err);
   fChain->SetBranchAddress("Event_Trig_Fired", &Event_Trig_Fired, &b_Event_Trig_Fired);
   fChain->SetBranchAddress("Event_Trig_Presc", &Event_Trig_Presc, &b_Event_Trig_Presc);
   Notify();
   
   fFriend1 = new muPFAna(treePF);
   fFriend2 = new muRecoAna(treeReco);
   fFriend3 = new METPFAna(treeMETPF);
}

Bool_t muEventAna::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void muEventAna::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t muEventAna::Cut(TString var, int iqq, float cutval, bool isPFidx)
{
   // int idx1=-1, idx2=-1;
   // int pfidx1=-1, pfidx2=-1;
   // if (isPFidx) {
   //    pfidx1 = PF_DiMuon_Muon1_Idx[iqq];
   //    pfidx2 = PF_DiMuon_Muon2_Idx[iqq];
   //    idx1 = PF_Muon_Reco_Idx[pfidx1];
   //    idx2 = PF_Muon_Reco_Idx[pfidx2];
   // } else {
   //    idx1 = Reco_DiMuon_Muon1_Idx[iqq];
   //    idx2 = Reco_DiMuon_Muon2_Idx[iqq];
   //    pfidx1 = Reco_Muon_PF_Idx[idx1];
   //    pfidx2 = Reco_Muon_PF_Idx[idx2];
   // }

   // float mu1=0, mu2=0;
   // if (var == "PF_IsoPFR03") {
   //    if (pfidx1<0 || pfidx2<0) return false;
   //    mu1 = PF_Muon_IsoPFR03[pfidx1];
   //    mu2 = PF_Muon_IsoPFR03[pfidx2];
   //    return (mu1>=cutval && mu2>=cutval);
   // }
   // if (var == "PF_IsoPFR03NoPUCorr") {
   //    if (pfidx1<0 || pfidx2<0) return false;
   //    mu1 = PF_Muon_IsoPFR03NoPUCorr[pfidx1];
   //    mu2 = PF_Muon_IsoPFR03NoPUCorr[pfidx2];
   //    return (mu1>=cutval && mu2>=cutval);
   // }
   // if (var == "PF_IsoPFR04") {
   //    if (pfidx1<0 || pfidx2<0) return false;
   //    mu1 = PF_Muon_IsoPFR04[pfidx1];
   //    mu2 = PF_Muon_IsoPFR04[pfidx2];
   //    return (mu1>=cutval && mu2>=cutval);
   // }
   // if (var == "PF_IsoPFR04NoPUCorr") {
   //    if (pfidx1<0 || pfidx2<0) return false;
   //    mu1 = PF_Muon_IsoPFR04NoPUCorr[pfidx1];
   //    mu2 = PF_Muon_IsoPFR04NoPUCorr[pfidx2];
   //    return (mu1>=cutval && mu2>=cutval);
   // }
   // if (var == "PF_EM_Chg_sumR03Pt") {
   //    if (pfidx1<0 || pfidx2<0) return false;
   //    mu1 = PF_Muon_EM_Chg_sumR03Pt[pfidx1];
   //    mu2 = PF_Muon_EM_Chg_sumR03Pt[pfidx2];
   //    return (mu1>=cutval && mu2>=cutval);
   // }
   // if (var == "PF_EM_Chg_sumR04Pt") {
   //    if (pfidx1<0 || pfidx2<0) return false;
   //    mu1 = PF_Muon_EM_Chg_sumR04Pt[pfidx1];
   //    mu2 = PF_Muon_EM_Chg_sumR04Pt[pfidx2];
   //    return (mu1>=cutval && mu2>=cutval);
   // }
   // if (var == "PF_EM_Neu_sumR03Et") {
   //    if (pfidx1<0 || pfidx2<0) return false;
   //    mu1 = PF_Muon_EM_Neu_sumR03Et[pfidx1];
   //    mu2 = PF_Muon_EM_Neu_sumR03Et[pfidx2];
   //    return (mu1>=cutval && mu2>=cutval);
   // }
   // if (var == "PF_EM_Neu_sumR04Et") {
   //    if (pfidx1<0 || pfidx2<0) return false;
   //    mu1 = PF_Muon_EM_Neu_sumR04Et[pfidx1];
   //    mu2 = PF_Muon_EM_Neu_sumR04Et[pfidx2];
   //    return (mu1>=cutval && mu2>=cutval);
   // }
   // if (var == "PF_Had_Chg_sumR03Pt") {
   //    if (pfidx1<0 || pfidx2<0) return false;
   //    mu1 = PF_Muon_Had_Chg_sumR03Pt[pfidx1];
   //    mu2 = PF_Muon_Had_Chg_sumR03Pt[pfidx2];
   //    return (mu1>=cutval && mu2>=cutval);
   // }
   // if (var == "PF_Had_Chg_sumR04Pt") {
   //    if (pfidx1<0 || pfidx2<0) return false;
   //    mu1 = PF_Muon_Had_Chg_sumR04Pt[pfidx1];
   //    mu2 = PF_Muon_Had_Chg_sumR04Pt[pfidx2];
   //    return (mu1>=cutval && mu2>=cutval);
   // }
   // if (var == "PF_Had_Neu_sumR03Et") {
   //    if (pfidx1<0 || pfidx2<0) return false;
   //    mu1 = PF_Muon_Had_Neu_sumR03Et[pfidx1];
   //    mu2 = PF_Muon_Had_Neu_sumR03Et[pfidx2];
   //    return (mu1>=cutval && mu2>=cutval);
   // }
   // if (var == "PF_Had_Neu_sumR04Et") {
   //    if (pfidx1<0 || pfidx2<0) return false;
   //    mu1 = PF_Muon_Had_Neu_sumR04Et[pfidx1];
   //    mu2 = PF_Muon_Had_Neu_sumR04Et[pfidx2];
   //    return (mu1>=cutval && mu2>=cutval);
   // }
   // if (var == "PF_Had_PU_sumR03Pt") {
   //    if (pfidx1<0 || pfidx2<0) return false;
   //    mu1 = PF_Muon_Had_PU_sumR03Pt[pfidx1];
   //    mu2 = PF_Muon_Had_PU_sumR03Pt[pfidx2];
   //    return (mu1>=cutval && mu2>=cutval);
   // }
   // if (var == "PF_Had_PU_sumR04Pt") {
   //    if (pfidx1<0 || pfidx2<0) return false;
   //    mu1 = PF_Muon_Had_PU_sumR04Pt[pfidx1];
   //    mu2 = PF_Muon_Had_PU_sumR04Pt[pfidx2];
   //    return (mu1>=cutval && mu2>=cutval);
   // }
   // if (var == "PF_VtxProb") {
   //    if (!isPFidx) return false;
   //    mu1 = PF_DiMuon_VtxProb[iqq];
   //    return (mu1>cutval);
   // }

   // // Reco variables
   // if (var == "Reco_isPF") {
   //    if (idx1<0 || idx2<0) return false;
   //    mu1 = Reco_Muon_isPF[idx1];
   //    mu2 = Reco_Muon_isPF[idx2];
   //    return (mu1>=cutval && mu2>=cutval);
   // }
   // if (var == "Reco_isGlobal") {
   //    if (idx1<0 || idx2<0) return false;
   //    mu1 = Reco_Muon_isGlobal[idx1];
   //    mu2 = Reco_Muon_isGlobal[idx2];
   //    return (mu1>=cutval && mu2>=cutval);
   // }
   // if (var == "Reco_isTracker") {
   //    if (idx1<0 || idx2<0) return false;
   //    mu1 = Reco_Muon_isTracker[idx1];
   //    mu2 = Reco_Muon_isTracker[idx2];
   //    return (mu1>=cutval && mu2>=cutval);
   // }
   // if (var == "Reco_isStandAlone") {
   //    if (idx1<0 || idx2<0) return false;
   //    mu1 = Reco_Muon_isStandAlone[idx1];
   //    mu2 = Reco_Muon_isStandAlone[idx2];
   //    return (mu1>=cutval && mu2>=cutval);
   // }
   // if (var == "Reco_isLoose") {
   //    if (idx1<0 || idx2<0) return false;
   //    mu1 = Reco_Muon_isLoose[idx1];
   //    mu2 = Reco_Muon_isLoose[idx2];
   //    return (mu1>=cutval && mu2>=cutval);
   // }
   // if (var == "Reco_isMedium") {
   //    if (idx1<0 || idx2<0) return false;
   //    mu1 = Reco_Muon_isMedium[idx1];
   //    mu2 = Reco_Muon_isMedium[idx2];
   //    return (mu1>=cutval && mu2>=cutval);
   // }
   // if (var == "Reco_isHighPt") {
   //    if (idx1<0 || idx2<0) return false;
   //    mu1 = Reco_Muon_isHighPt[idx1];
   //    mu2 = Reco_Muon_isHighPt[idx2];
   //    return (mu1>=cutval && mu2>=cutval);
   // }
   // if (var == "Reco_isSoft") {
   //    if (idx1<0 || idx2<0) return false;
   //    mu1 = Reco_Muon_isSoft[idx1];
   //    mu2 = Reco_Muon_isSoft[idx2];
   //    return (mu1>=cutval && mu2>=cutval);
   // }
   // if (var == "Reco_isTight") {
   //    if (idx1<0 || idx2<0) return false;
   //    mu1 = Reco_Muon_isTight[idx1];
   //    mu2 = Reco_Muon_isTight[idx2];
   //    return (mu1>=cutval && mu2>=cutval);
   // }
   // if (var == "Reco_isArbitrated") {
   //    if (idx1<0 || idx2<0) return false;
   //    mu1 = Reco_Muon_isArbitrated[idx1];
   //    mu2 = Reco_Muon_isArbitrated[idx2];
   //    return (mu1>=cutval && mu2>=cutval);
   // }
   // if (var == "Reco_TrackerArbitrated") {
   //    if (idx1<0 || idx2<0) return false;
   //    mu1 = Reco_Muon_TrackerArbitrated[idx1];
   //    mu2 = Reco_Muon_TrackerArbitrated[idx2];
   //    return (mu1>=cutval && mu2>=cutval);
   // }

   return false;
}
bool muEventAna::IsAccept(double pt, double eta)
{
  // // Acceptance cuts 2015 V2
  // 
  // return ( (fabs(eta)<1.2 && pt>=3.5) ||
  //         (1.2<=fabs(eta) && fabs(eta)<2.1 && pt>=5.77-1.89*fabs(eta)) ||
  //         (2.1<=fabs(eta) && fabs(eta)<2.4 && pt>=1.8) );

  // high pt acceptance
  return (fabs(eta)<2.4 && pt>15);
}
float muEventAna::Get(TString varname, int idx) {
   if (idx<0) return -1e99;
   if (varname == "PF_Candidate_isPU") return fFriend1->PF_Candidate_isPU->at(idx);
   if (varname == "PF_Candidate_Id") return fFriend1->PF_Candidate_Id->at(idx);
   if (varname == "PF_Candidate_Eta") return fFriend1->PF_Candidate_Eta->at(idx);
   if (varname == "PF_Candidate_Phi") return fFriend1->PF_Candidate_Phi->at(idx);
   if (varname == "PF_Candidate_Pt") return fFriend1->PF_Candidate_Pt->at(idx);
   if (varname == "PF_Muon_Charge") return fFriend1->PF_Muon_Charge->at(idx);
   if (varname == "PF_Muon_Reco_Idx") return fFriend1->PF_Muon_Reco_Idx->at(idx);
   if (varname == "PF_Muon_IsoPFR03") return fFriend1->PF_Muon_IsoPFR03->at(idx);
   if (varname == "PF_Muon_IsoPFR03NoPUCorr") return fFriend1->PF_Muon_IsoPFR03NoPUCorr->at(idx);
   if (varname == "PF_Muon_IsoPFR04") return fFriend1->PF_Muon_IsoPFR04->at(idx);
   if (varname == "PF_Muon_IsoPFR04NoPUCorr") return fFriend1->PF_Muon_IsoPFR04NoPUCorr->at(idx);
   if (varname == "PF_Muon_EM_Chg_sumR03Pt") return fFriend1->PF_Muon_EM_Chg_sumR03Pt->at(idx);
   if (varname == "PF_Muon_EM_Chg_sumR04Pt") return fFriend1->PF_Muon_EM_Chg_sumR04Pt->at(idx);
   if (varname == "PF_Muon_EM_Neu_sumR03Et") return fFriend1->PF_Muon_EM_Neu_sumR03Et->at(idx);
   if (varname == "PF_Muon_EM_Neu_sumR04Et") return fFriend1->PF_Muon_EM_Neu_sumR04Et->at(idx);
   if (varname == "PF_Muon_Had_Chg_sumR03Pt") return fFriend1->PF_Muon_Had_Chg_sumR03Pt->at(idx);
   if (varname == "PF_Muon_Had_Chg_sumR04Pt") return fFriend1->PF_Muon_Had_Chg_sumR04Pt->at(idx);
   if (varname == "PF_Muon_Had_Neu_sumR03Et") return fFriend1->PF_Muon_Had_Neu_sumR03Et->at(idx);
   if (varname == "PF_Muon_Had_Neu_sumR04Et") return fFriend1->PF_Muon_Had_Neu_sumR04Et->at(idx);
   if (varname == "PF_Muon_Had_PU_sumR03Pt") return fFriend1->PF_Muon_Had_PU_sumR03Pt->at(idx);
   if (varname == "PF_Muon_Had_PU_sumR04Pt") return fFriend1->PF_Muon_Had_PU_sumR04Pt->at(idx);
   if (varname == "PF_DiMuon_Charge") return fFriend1->PF_DiMuon_Charge->at(idx);
   if (varname == "PF_DiMuon_Muon1_Idx") return fFriend1->PF_DiMuon_Muon1_Idx->at(idx);
   if (varname == "PF_DiMuon_Muon2_Idx") return fFriend1->PF_DiMuon_Muon2_Idx->at(idx);
   if (varname == "PF_DiMuon_VtxProb") return fFriend1->PF_DiMuon_VtxProb->at(idx);
   if (varname == "PF_DiMuon_DCA") return fFriend1->PF_DiMuon_DCA->at(idx);
   if (varname == "PF_DiMuon_MassErr") return fFriend1->PF_DiMuon_MassErr->at(idx);

   if (varname == "Reco_Muon_Charge") return fFriend2->Reco_Muon_Charge->at(idx);
   if (varname == "Reco_Muon_PF_Idx") return fFriend2->Reco_Muon_PF_Idx->at(idx);
   if (varname == "Reco_Muon_isPF") return fFriend2->Reco_Muon_isPF->at(idx);
   if (varname == "Reco_Muon_isGlobal") return fFriend2->Reco_Muon_isGlobal->at(idx);
   if (varname == "Reco_Muon_isTracker") return fFriend2->Reco_Muon_isTracker->at(idx);
   if (varname == "Reco_Muon_isStandAlone") return fFriend2->Reco_Muon_isStandAlone->at(idx);
   if (varname == "Reco_Muon_isLoose") return fFriend2->Reco_Muon_isLoose->at(idx);
   if (varname == "Reco_Muon_isMedium") return fFriend2->Reco_Muon_isMedium->at(idx);
   if (varname == "Reco_Muon_isHighPt") return fFriend2->Reco_Muon_isHighPt->at(idx);
   if (varname == "Reco_Muon_isSoft") return fFriend2->Reco_Muon_isSoft->at(idx);
   if (varname == "Reco_Muon_isTight") return fFriend2->Reco_Muon_isTight->at(idx);
   if (varname == "Reco_Muon_isArbitrated") return fFriend2->Reco_Muon_isArbitrated->at(idx);
   if (varname == "Reco_Muon_TrackerArbitrated") return fFriend2->Reco_Muon_TrackerArbitrated->at(idx);
   if (varname == "Reco_Muon_GlobalPromptTight") return fFriend2->Reco_Muon_GlobalPromptTight->at(idx);
   if (varname == "Reco_Muon_TMLastStationLoose") return fFriend2->Reco_Muon_TMLastStationLoose->at(idx);
   if (varname == "Reco_Muon_TMLastStationTight") return fFriend2->Reco_Muon_TMLastStationTight->at(idx);
   if (varname == "Reco_Muon_TM2DCompatibilityLoose") return fFriend2->Reco_Muon_TM2DCompatibilityLoose->at(idx);
   if (varname == "Reco_Muon_TM2DCompatibilityTight") return fFriend2->Reco_Muon_TM2DCompatibilityTight->at(idx);
   if (varname == "Reco_Muon_TMOneStationLoose") return fFriend2->Reco_Muon_TMOneStationLoose->at(idx);
   if (varname == "Reco_Muon_TMOneStationTight") return fFriend2->Reco_Muon_TMOneStationTight->at(idx);
   if (varname == "Reco_Muon_GMTkChiCompatibility") return fFriend2->Reco_Muon_GMTkChiCompatibility->at(idx);
   if (varname == "Reco_Muon_GMStaChiCompatibility") return fFriend2->Reco_Muon_GMStaChiCompatibility->at(idx);
   if (varname == "Reco_Muon_GMTkKinkTight") return fFriend2->Reco_Muon_GMTkKinkTight->at(idx);
   if (varname == "Reco_Muon_TMLastStationAngLoose") return fFriend2->Reco_Muon_TMLastStationAngLoose->at(idx);
   if (varname == "Reco_Muon_TMLastStationAngTight") return fFriend2->Reco_Muon_TMLastStationAngTight->at(idx);
   if (varname == "Reco_Muon_TMOneStationAngLoose") return fFriend2->Reco_Muon_TMOneStationAngLoose->at(idx);
   if (varname == "Reco_Muon_TMOneStationAngTight") return fFriend2->Reco_Muon_TMOneStationAngTight->at(idx);
   if (varname == "Reco_Muon_MatchedStations") return fFriend2->Reco_Muon_MatchedStations->at(idx);
   if (varname == "Reco_Muon_Matches") return fFriend2->Reco_Muon_Matches->at(idx);
   if (varname == "Reco_Muon_SegmentComp") return fFriend2->Reco_Muon_SegmentComp->at(idx);
   if (varname == "Reco_Muon_Chi2Pos") return fFriend2->Reco_Muon_Chi2Pos->at(idx);
   if (varname == "Reco_Muon_TrkKink") return fFriend2->Reco_Muon_TrkKink->at(idx);
   if (varname == "Reco_Muon_InTrk_PtErr") return fFriend2->Reco_Muon_InTrk_PtErr->at(idx);
   if (varname == "Reco_Muon_InTrk_isHighPurity") return fFriend2->Reco_Muon_InTrk_isHighPurity->at(idx);
   if (varname == "Reco_Muon_InTrk_ValidHits") return fFriend2->Reco_Muon_InTrk_ValidHits->at(idx);
   if (varname == "Reco_Muon_InTrk_LostHits") return fFriend2->Reco_Muon_InTrk_LostHits->at(idx);
   if (varname == "Reco_Muon_InTrk_ValidPixHits") return fFriend2->Reco_Muon_InTrk_ValidPixHits->at(idx);
   if (varname == "Reco_Muon_InTrk_TrkLayers") return fFriend2->Reco_Muon_InTrk_TrkLayers->at(idx);
   if (varname == "Reco_Muon_InTrk_PixLayers") return fFriend2->Reco_Muon_InTrk_PixLayers->at(idx);
   if (varname == "Reco_Muon_InTrk_dXY") return fFriend2->Reco_Muon_InTrk_dXY->at(idx);
   if (varname == "Reco_Muon_InTrk_dXYErr") return fFriend2->Reco_Muon_InTrk_dXYErr->at(idx);
   if (varname == "Reco_Muon_InTrk_dZ") return fFriend2->Reco_Muon_InTrk_dZ->at(idx);
   if (varname == "Reco_Muon_InTrk_dZErr") return fFriend2->Reco_Muon_InTrk_dZErr->at(idx);
   if (varname == "Reco_Muon_InTrk_ValFrac") return fFriend2->Reco_Muon_InTrk_ValFrac->at(idx);
   if (varname == "Reco_Muon_InTrk_NormChi2") return fFriend2->Reco_Muon_InTrk_NormChi2->at(idx);
   if (varname == "Reco_Muon_GlbTrk_PtErr") return fFriend2->Reco_Muon_GlbTrk_PtErr->at(idx);
   if (varname == "Reco_Muon_GlbTrk_ValidMuonHits") return fFriend2->Reco_Muon_GlbTrk_ValidMuonHits->at(idx);
   if (varname == "Reco_Muon_GlbTrk_NormChi2") return fFriend2->Reco_Muon_GlbTrk_NormChi2->at(idx);
   if (varname == "Reco_Muon_BestTrk_Type") return fFriend2->Reco_Muon_BestTrk_Type->at(idx);
   if (varname == "Reco_Muon_BestTrk_PtErr") return fFriend2->Reco_Muon_BestTrk_PtErr->at(idx);
   if (varname == "Reco_Muon_BestTrk_dXY") return fFriend2->Reco_Muon_BestTrk_dXY->at(idx);
   if (varname == "Reco_Muon_BestTrk_dXYErr") return fFriend2->Reco_Muon_BestTrk_dXYErr->at(idx);
   if (varname == "Reco_Muon_BestTrk_dZ") return fFriend2->Reco_Muon_BestTrk_dZ->at(idx);
   if (varname == "Reco_Muon_BestTrk_dZErr") return fFriend2->Reco_Muon_BestTrk_dZErr->at(idx);
   if (varname == "Reco_Muon_IsoPFR03") return fFriend2->Reco_Muon_IsoPFR03->at(idx);
   if (varname == "Reco_Muon_IsoPFR03NoPUCorr") return fFriend2->Reco_Muon_IsoPFR03NoPUCorr->at(idx);
   if (varname == "Reco_Muon_IsoPFR04") return fFriend2->Reco_Muon_IsoPFR04->at(idx);
   if (varname == "Reco_Muon_IsoPFR04NoPUCorr") return fFriend2->Reco_Muon_IsoPFR04NoPUCorr->at(idx);
   if (varname == "Reco_Muon_EM_Chg_sumR03Pt") return fFriend2->Reco_Muon_EM_Chg_sumR03Pt->at(idx);
   if (varname == "Reco_Muon_EM_Chg_sumR04Pt") return fFriend2->Reco_Muon_EM_Chg_sumR04Pt->at(idx);
   if (varname == "Reco_Muon_EM_Neu_sumR03Et") return fFriend2->Reco_Muon_EM_Neu_sumR03Et->at(idx);
   if (varname == "Reco_Muon_EM_Neu_sumR04Et") return fFriend2->Reco_Muon_EM_Neu_sumR04Et->at(idx);
   if (varname == "Reco_Muon_Had_Chg_sumR03Pt") return fFriend2->Reco_Muon_Had_Chg_sumR03Pt->at(idx);
   if (varname == "Reco_Muon_Had_Chg_sumR04Pt") return fFriend2->Reco_Muon_Had_Chg_sumR04Pt->at(idx);
   if (varname == "Reco_Muon_Had_Neu_sumR03Et") return fFriend2->Reco_Muon_Had_Neu_sumR03Et->at(idx);
   if (varname == "Reco_Muon_Had_Neu_sumR04Et") return fFriend2->Reco_Muon_Had_Neu_sumR04Et->at(idx);
   if (varname == "Reco_Muon_Had_PU_sumR03Pt") return fFriend2->Reco_Muon_Had_PU_sumR03Pt->at(idx);
   if (varname == "Reco_Muon_Had_PU_sumR04Pt") return fFriend2->Reco_Muon_Had_PU_sumR04Pt->at(idx);
   if (varname == "Reco_Muon_IsoR03") return fFriend2->Reco_Muon_IsoR03->at(idx);
   if (varname == "Reco_Muon_IsoR05") return fFriend2->Reco_Muon_IsoR05->at(idx);
   if (varname == "Reco_Muon_Trk_sumR03Pt") return fFriend2->Reco_Muon_Trk_sumR03Pt->at(idx);
   if (varname == "Reco_Muon_Trk_sumR05Pt") return fFriend2->Reco_Muon_Trk_sumR05Pt->at(idx);
   if (varname == "Reco_DiMuon_Charge") return fFriend2->Reco_DiMuon_Charge->at(idx);
   if (varname == "Reco_DiMuon_Muon1_Idx") return fFriend2->Reco_DiMuon_Muon1_Idx->at(idx);
   if (varname == "Reco_DiMuon_Muon2_Idx") return fFriend2->Reco_DiMuon_Muon2_Idx->at(idx);
   if (varname == "Reco_DiMuon_isCowBoy") return fFriend2->Reco_DiMuon_isCowBoy->at(idx);
   if (varname == "Reco_DiMuon_VtxProb") return fFriend2->Reco_DiMuon_VtxProb->at(idx);
   if (varname == "Reco_DiMuon_DCA") return fFriend2->Reco_DiMuon_DCA->at(idx);
   if (varname == "Reco_DiMuon_MassErr") return fFriend2->Reco_DiMuon_MassErr->at(idx);
   
   // custom made isolation
   if (varname == "PF_Muon_myIsoPFR010") return myPFIso(idx,0.1);
   if (varname == "PF_Muon_myIsoPFR015") return myPFIso(idx,0.15);
   if (varname == "PF_Muon_myIsoPFR020") return myPFIso(idx,0.2);
   if (varname == "PF_Muon_myIsoPFR025") return myPFIso(idx,0.25);
   if (varname == "PF_Muon_myIsoPFR030") return myPFIso(idx,0.3);
   if (varname == "PF_Muon_myIsoPFR035") return myPFIso(idx,0.35);
   if (varname == "PF_Muon_myIsoPFR040") return myPFIso(idx,0.4);
   if (varname == "PF_Muon_myIsoPFR045") return myPFIso(idx,0.45);

   return -1e99;
}

float muEventAna::myPFIso(int idx, float deltaR) {
   TLorentzVector *tlvpfmu = (TLorentzVector*) fFriend1->PF_Muon_Mom->At(idx);
   double ptmu = tlvpfmu->Pt();
   double val=0;
   for (unsigned int i=0; i<fFriend1->PF_Candidate_isPU->size(); i++) {
      if (fFriend1->PF_Candidate_isPU->at(i)) continue;
      int id = fFriend1->PF_Candidate_Id->at(i);
      if (!(id==1 || id==4 || id==5)) continue; // keep only h, gamma, h0
      TLorentzVector tlvpf; tlvpf.SetPtEtaPhiM(fFriend1->PF_Candidate_Pt->at(i), fFriend1->PF_Candidate_Eta->at(i), fFriend1->PF_Candidate_Phi->at(i), 0);
      double dr = tlvpfmu->DeltaR(tlvpf);
      if (dr<0.0001 || dr>deltaR) continue;
      val += tlvpf.Pt();
   }
   return val/ptmu;
}

#endif // #ifdef muEventAna_cxx

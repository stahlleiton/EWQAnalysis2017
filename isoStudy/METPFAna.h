//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Feb 24 14:13:35 2017 by ROOT version 6.08/04
// from TTree MET_PF/
// found on file: root://eoscms//eos/cms/store/group/phys_heavyions/anstahll/EWQAnalysis2017/pPb2016/8160GeV/Data/HiEWQForest_PASingleMuon_pPb_Pbp_8160GeV_20170213.root
//////////////////////////////////////////////////////////

#ifndef METPFAna_h
#define METPFAna_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TVector2.h"

class METPFAna {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   TVector2        *PF_MET_Mom;
   Float_t         PF_MET_Muon_Et;
   Float_t         PF_MET_Muon_EtFrac;
   Float_t         PF_MET_EM_Chg_Et;
   Float_t         PF_MET_EM_Chg_EtFrac;
   Float_t         PF_MET_EM_Neu_Et;
   Float_t         PF_MET_EM_Neu_EtFrac;
   Float_t         PF_MET_EM_HF_Et;
   Float_t         PF_MET_EM_HF_EtFrac;
   Float_t         PF_MET_Had_Chg_Et;
   Float_t         PF_MET_Had_Chg_EtFrac;
   Float_t         PF_MET_Had_Neu_Et;
   Float_t         PF_MET_Had_Neu_EtFrac;
   Float_t         PF_MET_Had_HF_Et;
   Float_t         PF_MET_Had_HF_EtFrac;
   TVector2        *PF_MET_NoShift_Mom;
   Float_t         PF_MET_NoShift_sumEt;
   TVector2        *PF_MET_ElectronEnDown_Mom;
   Float_t         PF_MET_ElectronEnDown_sumEt;
   TVector2        *PF_MET_ElectronEnUp_Mom;
   Float_t         PF_MET_ElectronEnUp_sumEt;
   TVector2        *PF_MET_JetEnDown_Mom;
   Float_t         PF_MET_JetEnDown_sumEt;
   TVector2        *PF_MET_JetEnUp_Mom;
   Float_t         PF_MET_JetEnUp_sumEt;
   TVector2        *PF_MET_JetResDown_Mom;
   Float_t         PF_MET_JetResDown_sumEt;
   TVector2        *PF_MET_JetResUp_Mom;
   Float_t         PF_MET_JetResUp_sumEt;
   TVector2        *PF_MET_MuonEnDown_Mom;
   Float_t         PF_MET_MuonEnDown_sumEt;
   TVector2        *PF_MET_MuonEnUp_Mom;
   Float_t         PF_MET_MuonEnUp_sumEt;
   TVector2        *PF_MET_PhotonEnDown_Mom;
   Float_t         PF_MET_PhotonEnDown_sumEt;
   TVector2        *PF_MET_PhotonEnUp_Mom;
   Float_t         PF_MET_PhotonEnUp_sumEt;
   TVector2        *PF_MET_TauEnDown_Mom;
   Float_t         PF_MET_TauEnDown_sumEt;
   TVector2        *PF_MET_TauEnUp_Mom;
   Float_t         PF_MET_TauEnUp_sumEt;
   TVector2        *PF_MET_UnclusEnDown_Mom;
   Float_t         PF_MET_UnclusEnDown_sumEt;
   TVector2        *PF_MET_UnclusEnUp_Mom;
   Float_t         PF_MET_UnclusEnUp_sumEt;

   // List of branches
   TBranch        *b_PF_MET_Mom;   //!
   TBranch        *b_PF_MET_Muon_Et;   //!
   TBranch        *b_PF_MET_Muon_EtFrac;   //!
   TBranch        *b_PF_MET_EM_Chg_Et;   //!
   TBranch        *b_PF_MET_EM_Chg_EtFrac;   //!
   TBranch        *b_PF_MET_EM_Neu_Et;   //!
   TBranch        *b_PF_MET_EM_Neu_EtFrac;   //!
   TBranch        *b_PF_MET_EM_HF_Et;   //!
   TBranch        *b_PF_MET_EM_HF_EtFrac;   //!
   TBranch        *b_PF_MET_Had_Chg_Et;   //!
   TBranch        *b_PF_MET_Had_Chg_EtFrac;   //!
   TBranch        *b_PF_MET_Had_Neu_Et;   //!
   TBranch        *b_PF_MET_Had_Neu_EtFrac;   //!
   TBranch        *b_PF_MET_Had_HF_Et;   //!
   TBranch        *b_PF_MET_Had_HF_EtFrac;   //!
   TBranch        *b_PF_MET_NoShift_Mom;   //!
   TBranch        *b_PF_MET_NoShift_sumEt;   //!
   TBranch        *b_PF_MET_ElectronEnDown_Mom;   //!
   TBranch        *b_PF_MET_ElectronEnDown_sumEt;   //!
   TBranch        *b_PF_MET_ElectronEnUp_Mom;   //!
   TBranch        *b_PF_MET_ElectronEnUp_sumEt;   //!
   TBranch        *b_PF_MET_JetEnDown_Mom;   //!
   TBranch        *b_PF_MET_JetEnDown_sumEt;   //!
   TBranch        *b_PF_MET_JetEnUp_Mom;   //!
   TBranch        *b_PF_MET_JetEnUp_sumEt;   //!
   TBranch        *b_PF_MET_JetResDown_Mom;   //!
   TBranch        *b_PF_MET_JetResDown_sumEt;   //!
   TBranch        *b_PF_MET_JetResUp_Mom;   //!
   TBranch        *b_PF_MET_JetResUp_sumEt;   //!
   TBranch        *b_PF_MET_MuonEnDown_Mom;   //!
   TBranch        *b_PF_MET_MuonEnDown_sumEt;   //!
   TBranch        *b_PF_MET_MuonEnUp_Mom;   //!
   TBranch        *b_PF_MET_MuonEnUp_sumEt;   //!
   TBranch        *b_PF_MET_PhotonEnDown_Mom;   //!
   TBranch        *b_PF_MET_PhotonEnDown_sumEt;   //!
   TBranch        *b_PF_MET_PhotonEnUp_Mom;   //!
   TBranch        *b_PF_MET_PhotonEnUp_sumEt;   //!
   TBranch        *b_PF_MET_TauEnDown_Mom;   //!
   TBranch        *b_PF_MET_TauEnDown_sumEt;   //!
   TBranch        *b_PF_MET_TauEnUp_Mom;   //!
   TBranch        *b_PF_MET_TauEnUp_sumEt;   //!
   TBranch        *b_PF_MET_UnclusEnDown_Mom;   //!
   TBranch        *b_PF_MET_UnclusEnDown_sumEt;   //!
   TBranch        *b_PF_MET_UnclusEnUp_Mom;   //!
   TBranch        *b_PF_MET_UnclusEnUp_sumEt;   //!

   METPFAna(TTree *tree=0);
   virtual ~METPFAna();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef METPFAna_cxx
METPFAna::METPFAna(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("root://eoscms//eos/cms/store/group/phys_heavyions/anstahll/EWQAnalysis2017/pPb2016/8160GeV/Data/HiEWQForest_PASingleMuon_pPb_Pbp_8160GeV_20170213.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("root://eoscms//eos/cms/store/group/phys_heavyions/anstahll/EWQAnalysis2017/pPb2016/8160GeV/Data/HiEWQForest_PASingleMuon_pPb_Pbp_8160GeV_20170213.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("root://eoscms//eos/cms/store/group/phys_heavyions/anstahll/EWQAnalysis2017/pPb2016/8160GeV/Data/HiEWQForest_PASingleMuon_pPb_Pbp_8160GeV_20170213.root:/metAna");
      dir->GetObject("MET_PF",tree);

   }
   Init(tree);
}

METPFAna::~METPFAna()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t METPFAna::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t METPFAna::LoadTree(Long64_t entry)
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

void METPFAna::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   PF_MET_Mom = 0;
   PF_MET_NoShift_Mom = 0;
   PF_MET_ElectronEnDown_Mom = 0;
   PF_MET_ElectronEnUp_Mom = 0;
   PF_MET_JetEnDown_Mom = 0;
   PF_MET_JetEnUp_Mom = 0;
   PF_MET_JetResDown_Mom = 0;
   PF_MET_JetResUp_Mom = 0;
   PF_MET_MuonEnDown_Mom = 0;
   PF_MET_MuonEnUp_Mom = 0;
   PF_MET_PhotonEnDown_Mom = 0;
   PF_MET_PhotonEnUp_Mom = 0;
   PF_MET_TauEnDown_Mom = 0;
   PF_MET_TauEnUp_Mom = 0;
   PF_MET_UnclusEnDown_Mom = 0;
   PF_MET_UnclusEnUp_Mom = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("PF_MET_Mom", &PF_MET_Mom, &b_PF_MET_Mom);
   fChain->SetBranchAddress("PF_MET_Muon_Et", &PF_MET_Muon_Et, &b_PF_MET_Muon_Et);
   fChain->SetBranchAddress("PF_MET_Muon_EtFrac", &PF_MET_Muon_EtFrac, &b_PF_MET_Muon_EtFrac);
   fChain->SetBranchAddress("PF_MET_EM_Chg_Et", &PF_MET_EM_Chg_Et, &b_PF_MET_EM_Chg_Et);
   fChain->SetBranchAddress("PF_MET_EM_Chg_EtFrac", &PF_MET_EM_Chg_EtFrac, &b_PF_MET_EM_Chg_EtFrac);
   fChain->SetBranchAddress("PF_MET_EM_Neu_Et", &PF_MET_EM_Neu_Et, &b_PF_MET_EM_Neu_Et);
   fChain->SetBranchAddress("PF_MET_EM_Neu_EtFrac", &PF_MET_EM_Neu_EtFrac, &b_PF_MET_EM_Neu_EtFrac);
   fChain->SetBranchAddress("PF_MET_EM_HF_Et", &PF_MET_EM_HF_Et, &b_PF_MET_EM_HF_Et);
   fChain->SetBranchAddress("PF_MET_EM_HF_EtFrac", &PF_MET_EM_HF_EtFrac, &b_PF_MET_EM_HF_EtFrac);
   fChain->SetBranchAddress("PF_MET_Had_Chg_Et", &PF_MET_Had_Chg_Et, &b_PF_MET_Had_Chg_Et);
   fChain->SetBranchAddress("PF_MET_Had_Chg_EtFrac", &PF_MET_Had_Chg_EtFrac, &b_PF_MET_Had_Chg_EtFrac);
   fChain->SetBranchAddress("PF_MET_Had_Neu_Et", &PF_MET_Had_Neu_Et, &b_PF_MET_Had_Neu_Et);
   fChain->SetBranchAddress("PF_MET_Had_Neu_EtFrac", &PF_MET_Had_Neu_EtFrac, &b_PF_MET_Had_Neu_EtFrac);
   fChain->SetBranchAddress("PF_MET_Had_HF_Et", &PF_MET_Had_HF_Et, &b_PF_MET_Had_HF_Et);
   fChain->SetBranchAddress("PF_MET_Had_HF_EtFrac", &PF_MET_Had_HF_EtFrac, &b_PF_MET_Had_HF_EtFrac);
   fChain->SetBranchAddress("PF_MET_NoShift_Mom", &PF_MET_NoShift_Mom, &b_PF_MET_NoShift_Mom);
   fChain->SetBranchAddress("PF_MET_NoShift_sumEt", &PF_MET_NoShift_sumEt, &b_PF_MET_NoShift_sumEt);
   fChain->SetBranchAddress("PF_MET_ElectronEnDown_Mom", &PF_MET_ElectronEnDown_Mom, &b_PF_MET_ElectronEnDown_Mom);
   fChain->SetBranchAddress("PF_MET_ElectronEnDown_sumEt", &PF_MET_ElectronEnDown_sumEt, &b_PF_MET_ElectronEnDown_sumEt);
   fChain->SetBranchAddress("PF_MET_ElectronEnUp_Mom", &PF_MET_ElectronEnUp_Mom, &b_PF_MET_ElectronEnUp_Mom);
   fChain->SetBranchAddress("PF_MET_ElectronEnUp_sumEt", &PF_MET_ElectronEnUp_sumEt, &b_PF_MET_ElectronEnUp_sumEt);
   fChain->SetBranchAddress("PF_MET_JetEnDown_Mom", &PF_MET_JetEnDown_Mom, &b_PF_MET_JetEnDown_Mom);
   fChain->SetBranchAddress("PF_MET_JetEnDown_sumEt", &PF_MET_JetEnDown_sumEt, &b_PF_MET_JetEnDown_sumEt);
   fChain->SetBranchAddress("PF_MET_JetEnUp_Mom", &PF_MET_JetEnUp_Mom, &b_PF_MET_JetEnUp_Mom);
   fChain->SetBranchAddress("PF_MET_JetEnUp_sumEt", &PF_MET_JetEnUp_sumEt, &b_PF_MET_JetEnUp_sumEt);
   fChain->SetBranchAddress("PF_MET_JetResDown_Mom", &PF_MET_JetResDown_Mom, &b_PF_MET_JetResDown_Mom);
   fChain->SetBranchAddress("PF_MET_JetResDown_sumEt", &PF_MET_JetResDown_sumEt, &b_PF_MET_JetResDown_sumEt);
   fChain->SetBranchAddress("PF_MET_JetResUp_Mom", &PF_MET_JetResUp_Mom, &b_PF_MET_JetResUp_Mom);
   fChain->SetBranchAddress("PF_MET_JetResUp_sumEt", &PF_MET_JetResUp_sumEt, &b_PF_MET_JetResUp_sumEt);
   fChain->SetBranchAddress("PF_MET_MuonEnDown_Mom", &PF_MET_MuonEnDown_Mom, &b_PF_MET_MuonEnDown_Mom);
   fChain->SetBranchAddress("PF_MET_MuonEnDown_sumEt", &PF_MET_MuonEnDown_sumEt, &b_PF_MET_MuonEnDown_sumEt);
   fChain->SetBranchAddress("PF_MET_MuonEnUp_Mom", &PF_MET_MuonEnUp_Mom, &b_PF_MET_MuonEnUp_Mom);
   fChain->SetBranchAddress("PF_MET_MuonEnUp_sumEt", &PF_MET_MuonEnUp_sumEt, &b_PF_MET_MuonEnUp_sumEt);
   fChain->SetBranchAddress("PF_MET_PhotonEnDown_Mom", &PF_MET_PhotonEnDown_Mom, &b_PF_MET_PhotonEnDown_Mom);
   fChain->SetBranchAddress("PF_MET_PhotonEnDown_sumEt", &PF_MET_PhotonEnDown_sumEt, &b_PF_MET_PhotonEnDown_sumEt);
   fChain->SetBranchAddress("PF_MET_PhotonEnUp_Mom", &PF_MET_PhotonEnUp_Mom, &b_PF_MET_PhotonEnUp_Mom);
   fChain->SetBranchAddress("PF_MET_PhotonEnUp_sumEt", &PF_MET_PhotonEnUp_sumEt, &b_PF_MET_PhotonEnUp_sumEt);
   fChain->SetBranchAddress("PF_MET_TauEnDown_Mom", &PF_MET_TauEnDown_Mom, &b_PF_MET_TauEnDown_Mom);
   fChain->SetBranchAddress("PF_MET_TauEnDown_sumEt", &PF_MET_TauEnDown_sumEt, &b_PF_MET_TauEnDown_sumEt);
   fChain->SetBranchAddress("PF_MET_TauEnUp_Mom", &PF_MET_TauEnUp_Mom, &b_PF_MET_TauEnUp_Mom);
   fChain->SetBranchAddress("PF_MET_TauEnUp_sumEt", &PF_MET_TauEnUp_sumEt, &b_PF_MET_TauEnUp_sumEt);
   fChain->SetBranchAddress("PF_MET_UnclusEnDown_Mom", &PF_MET_UnclusEnDown_Mom, &b_PF_MET_UnclusEnDown_Mom);
   fChain->SetBranchAddress("PF_MET_UnclusEnDown_sumEt", &PF_MET_UnclusEnDown_sumEt, &b_PF_MET_UnclusEnDown_sumEt);
   fChain->SetBranchAddress("PF_MET_UnclusEnUp_Mom", &PF_MET_UnclusEnUp_Mom, &b_PF_MET_UnclusEnUp_Mom);
   fChain->SetBranchAddress("PF_MET_UnclusEnUp_sumEt", &PF_MET_UnclusEnUp_sumEt, &b_PF_MET_UnclusEnUp_sumEt);
   Notify();
}

Bool_t METPFAna::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void METPFAna::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t METPFAna::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef METPFAna_cxx

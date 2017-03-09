//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Feb 22 19:02:06 2017 by ROOT version 6.08/04
// from TTree Muon_PF/
// found on file: root://eoscms//eos/cms/store/group/phys_heavyions/anstahll/EWQAnalysis2017/pPb2016/8160GeV/Data/HiEWQForest_PASingleMuon_pPb_Pbp_8160GeV_20170213.root
//////////////////////////////////////////////////////////

#ifndef muPFAna_h
#define muPFAna_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "TClonesArray.h"
#include "vector"
#include "vector"
#include "vector"
#include "TVector2.h"
#include "TLorentzVector.h"

class muPFAna {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   vector<bool>    *PF_Candidate_isPU;
   vector<unsigned char> *PF_Candidate_Id;
   vector<float>   *PF_Candidate_Eta;
   vector<float>   *PF_Candidate_Phi;
   vector<float>   *PF_Candidate_Pt;
   TClonesArray    *PF_Muon_Mom;
   vector<char>    *PF_Muon_Charge;
   vector<short>   *PF_Muon_Reco_Idx;
   vector<float>   *PF_Muon_IsoPFR03;
   vector<float>   *PF_Muon_IsoPFR03NoPUCorr;
   vector<float>   *PF_Muon_IsoPFR04;
   vector<float>   *PF_Muon_IsoPFR04NoPUCorr;
   vector<float>   *PF_Muon_EM_Chg_sumR03Pt;
   vector<float>   *PF_Muon_EM_Chg_sumR04Pt;
   vector<float>   *PF_Muon_EM_Neu_sumR03Et;
   vector<float>   *PF_Muon_EM_Neu_sumR04Et;
   vector<float>   *PF_Muon_Had_Chg_sumR03Pt;
   vector<float>   *PF_Muon_Had_Chg_sumR04Pt;
   vector<float>   *PF_Muon_Had_Neu_sumR03Et;
   vector<float>   *PF_Muon_Had_Neu_sumR04Et;
   vector<float>   *PF_Muon_Had_PU_sumR03Pt;
   vector<float>   *PF_Muon_Had_PU_sumR04Pt;
   TClonesArray    *PF_DiMuon_Mom;
   vector<char>    *PF_DiMuon_Charge;
   vector<unsigned short> *PF_DiMuon_Muon1_Idx;
   vector<unsigned short> *PF_DiMuon_Muon2_Idx;
   TClonesArray    *PF_DiMuon_Vertex;
   vector<float>   *PF_DiMuon_VtxProb;
   vector<float>   *PF_DiMuon_DCA;
   vector<float>   *PF_DiMuon_MassErr;
   TVector2        *PF_MET_Mom;
   TClonesArray    *PF_MuonMET_TransMom;

   // List of branches
   TBranch        *b_PF_Candidate_isPU;   //!
   TBranch        *b_PF_Candidate_Id;   //!
   TBranch        *b_PF_Candidate_Eta;   //!
   TBranch        *b_PF_Candidate_Phi;   //!
   TBranch        *b_PF_Candidate_Pt;   //!
   TBranch        *b_PF_Muon_Mom;   //!
   TBranch        *b_PF_Muon_Charge;   //!
   TBranch        *b_PF_Muon_Reco_Idx;   //!
   TBranch        *b_PF_Muon_IsoPFR03;   //!
   TBranch        *b_PF_Muon_IsoPFR03NoPUCorr;   //!
   TBranch        *b_PF_Muon_IsoPFR04;   //!
   TBranch        *b_PF_Muon_IsoPFR04NoPUCorr;   //!
   TBranch        *b_PF_Muon_EM_Chg_sumR03Pt;   //!
   TBranch        *b_PF_Muon_EM_Chg_sumR04Pt;   //!
   TBranch        *b_PF_Muon_EM_Neu_sumR03Et;   //!
   TBranch        *b_PF_Muon_EM_Neu_sumR04Et;   //!
   TBranch        *b_PF_Muon_Had_Chg_sumR03Pt;   //!
   TBranch        *b_PF_Muon_Had_Chg_sumR04Pt;   //!
   TBranch        *b_PF_Muon_Had_Neu_sumR03Et;   //!
   TBranch        *b_PF_Muon_Had_Neu_sumR04Et;   //!
   TBranch        *b_PF_Muon_Had_PU_sumR03Pt;   //!
   TBranch        *b_PF_Muon_Had_PU_sumR04Pt;   //!
   TBranch        *b_PF_DiMuon_Mom;   //!
   TBranch        *b_PF_DiMuon_Charge;   //!
   TBranch        *b_PF_DiMuon_Muon1_Idx;   //!
   TBranch        *b_PF_DiMuon_Muon2_Idx;   //!
   TBranch        *b_PF_DiMuon_Vertex;   //!
   TBranch        *b_PF_DiMuon_VtxProb;   //!
   TBranch        *b_PF_DiMuon_DCA;   //!
   TBranch        *b_PF_DiMuon_MassErr;   //!
   TBranch        *b_PF_MET_Mom;   //!
   TBranch        *b_PF_MuonMET_TransMom;   //!

   muPFAna(TTree *tree=0);
   virtual ~muPFAna();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   float Get(TString varname, int idx);
};

#endif

#ifdef muPFAna_cxx
muPFAna::muPFAna(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("root://eoscms//eos/cms/store/group/phys_heavyions/anstahll/EWQAnalysis2017/pPb2016/8160GeV/Data/HiEWQForest_PASingleMuon_pPb_Pbp_8160GeV_20170213.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("root://eoscms//eos/cms/store/group/phys_heavyions/anstahll/EWQAnalysis2017/pPb2016/8160GeV/Data/HiEWQForest_PASingleMuon_pPb_Pbp_8160GeV_20170213.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("root://eoscms//eos/cms/store/group/phys_heavyions/anstahll/EWQAnalysis2017/pPb2016/8160GeV/Data/HiEWQForest_PASingleMuon_pPb_Pbp_8160GeV_20170213.root:/muonAna");
      dir->GetObject("Muon_PF",tree);

   }
   Init(tree);
}

muPFAna::~muPFAna()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t muPFAna::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t muPFAna::LoadTree(Long64_t entry)
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

void muPFAna::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   PF_Candidate_isPU = 0;
   PF_Candidate_Id = 0;
   PF_Candidate_Eta = 0;
   PF_Candidate_Phi = 0;
   PF_Candidate_Pt = 0;
   PF_Muon_Mom = 0;
   PF_Muon_Charge = 0;
   PF_Muon_Reco_Idx = 0;
   PF_Muon_IsoPFR03 = 0;
   PF_Muon_IsoPFR03NoPUCorr = 0;
   PF_Muon_IsoPFR04 = 0;
   PF_Muon_IsoPFR04NoPUCorr = 0;
   PF_Muon_EM_Chg_sumR03Pt = 0;
   PF_Muon_EM_Chg_sumR04Pt = 0;
   PF_Muon_EM_Neu_sumR03Et = 0;
   PF_Muon_EM_Neu_sumR04Et = 0;
   PF_Muon_Had_Chg_sumR03Pt = 0;
   PF_Muon_Had_Chg_sumR04Pt = 0;
   PF_Muon_Had_Neu_sumR03Et = 0;
   PF_Muon_Had_Neu_sumR04Et = 0;
   PF_Muon_Had_PU_sumR03Pt = 0;
   PF_Muon_Had_PU_sumR04Pt = 0;
   PF_DiMuon_Mom = 0;
   PF_DiMuon_Charge = 0;
   PF_DiMuon_Muon1_Idx = 0;
   PF_DiMuon_Muon2_Idx = 0;
   PF_DiMuon_Vertex = 0;
   PF_DiMuon_VtxProb = 0;
   PF_DiMuon_DCA = 0;
   PF_DiMuon_MassErr = 0;
   PF_MET_Mom = 0;
   PF_MuonMET_TransMom = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("PF_Candidate_isPU", &PF_Candidate_isPU, &b_PF_Candidate_isPU);
   fChain->SetBranchAddress("PF_Candidate_Id", &PF_Candidate_Id, &b_PF_Candidate_Id);
   fChain->SetBranchAddress("PF_Candidate_Eta", &PF_Candidate_Eta, &b_PF_Candidate_Eta);
   fChain->SetBranchAddress("PF_Candidate_Phi", &PF_Candidate_Phi, &b_PF_Candidate_Phi);
   fChain->SetBranchAddress("PF_Candidate_Pt", &PF_Candidate_Pt, &b_PF_Candidate_Pt);
   fChain->SetBranchAddress("PF_Muon_Mom", &PF_Muon_Mom, &b_PF_Muon_Mom);
   fChain->SetBranchAddress("PF_Muon_Charge", &PF_Muon_Charge, &b_PF_Muon_Charge);
   fChain->SetBranchAddress("PF_Muon_Reco_Idx", &PF_Muon_Reco_Idx, &b_PF_Muon_Reco_Idx);
   fChain->SetBranchAddress("PF_Muon_IsoPFR03", &PF_Muon_IsoPFR03, &b_PF_Muon_IsoPFR03);
   fChain->SetBranchAddress("PF_Muon_IsoPFR03NoPUCorr", &PF_Muon_IsoPFR03NoPUCorr, &b_PF_Muon_IsoPFR03NoPUCorr);
   fChain->SetBranchAddress("PF_Muon_IsoPFR04", &PF_Muon_IsoPFR04, &b_PF_Muon_IsoPFR04);
   fChain->SetBranchAddress("PF_Muon_IsoPFR04NoPUCorr", &PF_Muon_IsoPFR04NoPUCorr, &b_PF_Muon_IsoPFR04NoPUCorr);
   fChain->SetBranchAddress("PF_Muon_EM_Chg_sumR03Pt", &PF_Muon_EM_Chg_sumR03Pt, &b_PF_Muon_EM_Chg_sumR03Pt);
   fChain->SetBranchAddress("PF_Muon_EM_Chg_sumR04Pt", &PF_Muon_EM_Chg_sumR04Pt, &b_PF_Muon_EM_Chg_sumR04Pt);
   fChain->SetBranchAddress("PF_Muon_EM_Neu_sumR03Et", &PF_Muon_EM_Neu_sumR03Et, &b_PF_Muon_EM_Neu_sumR03Et);
   fChain->SetBranchAddress("PF_Muon_EM_Neu_sumR04Et", &PF_Muon_EM_Neu_sumR04Et, &b_PF_Muon_EM_Neu_sumR04Et);
   fChain->SetBranchAddress("PF_Muon_Had_Chg_sumR03Pt", &PF_Muon_Had_Chg_sumR03Pt, &b_PF_Muon_Had_Chg_sumR03Pt);
   fChain->SetBranchAddress("PF_Muon_Had_Chg_sumR04Pt", &PF_Muon_Had_Chg_sumR04Pt, &b_PF_Muon_Had_Chg_sumR04Pt);
   fChain->SetBranchAddress("PF_Muon_Had_Neu_sumR03Et", &PF_Muon_Had_Neu_sumR03Et, &b_PF_Muon_Had_Neu_sumR03Et);
   fChain->SetBranchAddress("PF_Muon_Had_Neu_sumR04Et", &PF_Muon_Had_Neu_sumR04Et, &b_PF_Muon_Had_Neu_sumR04Et);
   fChain->SetBranchAddress("PF_Muon_Had_PU_sumR03Pt", &PF_Muon_Had_PU_sumR03Pt, &b_PF_Muon_Had_PU_sumR03Pt);
   fChain->SetBranchAddress("PF_Muon_Had_PU_sumR04Pt", &PF_Muon_Had_PU_sumR04Pt, &b_PF_Muon_Had_PU_sumR04Pt);
   fChain->SetBranchAddress("PF_DiMuon_Mom", &PF_DiMuon_Mom, &b_PF_DiMuon_Mom);
   fChain->SetBranchAddress("PF_DiMuon_Charge", &PF_DiMuon_Charge, &b_PF_DiMuon_Charge);
   fChain->SetBranchAddress("PF_DiMuon_Muon1_Idx", &PF_DiMuon_Muon1_Idx, &b_PF_DiMuon_Muon1_Idx);
   fChain->SetBranchAddress("PF_DiMuon_Muon2_Idx", &PF_DiMuon_Muon2_Idx, &b_PF_DiMuon_Muon2_Idx);
   fChain->SetBranchAddress("PF_DiMuon_Vertex", &PF_DiMuon_Vertex, &b_PF_DiMuon_Vertex);
   fChain->SetBranchAddress("PF_DiMuon_VtxProb", &PF_DiMuon_VtxProb, &b_PF_DiMuon_VtxProb);
   fChain->SetBranchAddress("PF_DiMuon_DCA", &PF_DiMuon_DCA, &b_PF_DiMuon_DCA);
   fChain->SetBranchAddress("PF_DiMuon_MassErr", &PF_DiMuon_MassErr, &b_PF_DiMuon_MassErr);
   fChain->SetBranchAddress("PF_MET_Mom", &PF_MET_Mom, &b_PF_MET_Mom);
   fChain->SetBranchAddress("PF_MuonMET_TransMom", &PF_MuonMET_TransMom, &b_PF_MuonMET_TransMom);
   Notify();
}

Bool_t muPFAna::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void muPFAna::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t muPFAna::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
float muPFAna::Get(TString varname, int idx) {
   if (varname == "PF_Candidate_isPU") return PF_Candidate_isPU->at(idx);
   if (varname == "PF_Candidate_Id") return PF_Candidate_Id->at(idx);
   if (varname == "PF_Candidate_Eta") return PF_Candidate_Eta->at(idx);
   if (varname == "PF_Candidate_Phi") return PF_Candidate_Phi->at(idx);
   if (varname == "PF_Candidate_Pt") return PF_Candidate_Pt->at(idx);
   if (varname == "PF_Muon_Charge") return PF_Muon_Charge->at(idx);
   if (varname == "PF_Muon_Reco_Idx") return PF_Muon_Reco_Idx->at(idx);
   if (varname == "PF_Muon_IsoPFR03") return PF_Muon_IsoPFR03->at(idx);
   if (varname == "PF_Muon_IsoPFR03NoPUCorr") return PF_Muon_IsoPFR03NoPUCorr->at(idx);
   if (varname == "PF_Muon_IsoPFR04") return PF_Muon_IsoPFR04->at(idx);
   if (varname == "PF_Muon_IsoPFR04NoPUCorr") return PF_Muon_IsoPFR04NoPUCorr->at(idx);
   if (varname == "PF_Muon_EM_Chg_sumR03Pt") return PF_Muon_EM_Chg_sumR03Pt->at(idx);
   if (varname == "PF_Muon_EM_Chg_sumR04Pt") return PF_Muon_EM_Chg_sumR04Pt->at(idx);
   if (varname == "PF_Muon_EM_Neu_sumR03Et") return PF_Muon_EM_Neu_sumR03Et->at(idx);
   if (varname == "PF_Muon_EM_Neu_sumR04Et") return PF_Muon_EM_Neu_sumR04Et->at(idx);
   if (varname == "PF_Muon_Had_Chg_sumR03Pt") return PF_Muon_Had_Chg_sumR03Pt->at(idx);
   if (varname == "PF_Muon_Had_Chg_sumR04Pt") return PF_Muon_Had_Chg_sumR04Pt->at(idx);
   if (varname == "PF_Muon_Had_Neu_sumR03Et") return PF_Muon_Had_Neu_sumR03Et->at(idx);
   if (varname == "PF_Muon_Had_Neu_sumR04Et") return PF_Muon_Had_Neu_sumR04Et->at(idx);
   if (varname == "PF_Muon_Had_PU_sumR03Pt") return PF_Muon_Had_PU_sumR03Pt->at(idx);
   if (varname == "PF_Muon_Had_PU_sumR04Pt") return PF_Muon_Had_PU_sumR04Pt->at(idx);
   if (varname == "PF_DiMuon_Charge") return PF_DiMuon_Charge->at(idx);
   if (varname == "PF_DiMuon_Muon1_Idx") return PF_DiMuon_Muon1_Idx->at(idx);
   if (varname == "PF_DiMuon_Muon2_Idx") return PF_DiMuon_Muon2_Idx->at(idx);
   if (varname == "PF_DiMuon_VtxProb") return PF_DiMuon_VtxProb->at(idx);
   if (varname == "PF_DiMuon_DCA") return PF_DiMuon_DCA->at(idx);
   if (varname == "PF_DiMuon_MassErr") return PF_DiMuon_MassErr->at(idx);
   return -1e99;
}
#endif // #ifdef muPFAna_cxx

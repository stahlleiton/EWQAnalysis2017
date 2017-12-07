//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Dec  5 15:49:23 2017 by ROOT version 6.06/00
// from TTree Muon_Reco/
// found on file: root://cms-xrd-global.cern.ch//store/group/phys_heavyions/anstahll/EWQAnalysis2017/pPb2016/8160GeV/Data/HiEWQForest_PASingleMuon_pPb_Pbp_8160GeV_20171003.root
//////////////////////////////////////////////////////////

#ifndef muRecoAna_h
#define muRecoAna_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TClonesArray.h"
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"

class muRecoAna {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UChar_t         Reco_Muon_N;
   TClonesArray    *Reco_Muon_Mom;
   vector<char>    *Reco_Muon_Charge;
   vector<char>    *Reco_Muon_PF_Idx;
   vector<vector<unsigned char> > *Pat_Muon_Trig;
   vector<float>   *Pat_Muon_dB;
   vector<float>   *Pat_Muon_dBErr;
   vector<bool>    *Reco_Muon_isPF;
   vector<bool>    *Reco_Muon_isGlobal;
   vector<bool>    *Reco_Muon_isTracker;
   vector<bool>    *Reco_Muon_isStandAlone;
   vector<bool>    *Reco_Muon_isLoose;
   vector<bool>    *Reco_Muon_isMedium;
   vector<bool>    *Reco_Muon_isHighPt;
   vector<bool>    *Reco_Muon_isSoft;
   vector<bool>    *Reco_Muon_isTight;
   vector<bool>    *Reco_Muon_isArbitrated;
   vector<bool>    *Reco_Muon_TrackerArbitrated;
   vector<bool>    *Reco_Muon_GlobalPromptTight;
   vector<bool>    *Reco_Muon_TMLastStationLoose;
   vector<bool>    *Reco_Muon_TMLastStationTight;
   vector<bool>    *Reco_Muon_TM2DCompatibilityLoose;
   vector<bool>    *Reco_Muon_TM2DCompatibilityTight;
   vector<bool>    *Reco_Muon_TMOneStationLoose;
   vector<bool>    *Reco_Muon_TMOneStationTight;
   vector<bool>    *Reco_Muon_GMTkChiCompatibility;
   vector<bool>    *Reco_Muon_GMStaChiCompatibility;
   vector<bool>    *Reco_Muon_GMTkKinkTight;
   vector<bool>    *Reco_Muon_TMLastStationAngLoose;
   vector<bool>    *Reco_Muon_TMLastStationAngTight;
   vector<bool>    *Reco_Muon_TMOneStationAngLoose;
   vector<bool>    *Reco_Muon_TMOneStationAngTight;
   vector<short>   *Reco_Muon_MatchedStations;
   vector<short>   *Reco_Muon_Matches;
   vector<float>   *Reco_Muon_SegmentComp;
   vector<float>   *Reco_Muon_Chi2Pos;
   vector<float>   *Reco_Muon_TrkKink;
   TClonesArray    *Reco_Muon_InTrk_Mom;
   vector<float>   *Reco_Muon_InTrk_PtErr;
   vector<bool>    *Reco_Muon_InTrk_isHighPurity;
   vector<short>   *Reco_Muon_InTrk_ValidHits;
   vector<short>   *Reco_Muon_InTrk_LostHits;
   vector<short>   *Reco_Muon_InTrk_ValidPixHits;
   vector<short>   *Reco_Muon_InTrk_TrkLayers;
   vector<short>   *Reco_Muon_InTrk_PixLayers;
   vector<float>   *Reco_Muon_InTrk_dXY;
   vector<float>   *Reco_Muon_InTrk_dXYErr;
   vector<float>   *Reco_Muon_InTrk_dZ;
   vector<float>   *Reco_Muon_InTrk_dZErr;
   vector<float>   *Reco_Muon_InTrk_ValFrac;
   vector<float>   *Reco_Muon_InTrk_NormChi2;
   TClonesArray    *Reco_Muon_GlbTrk_Mom;
   vector<float>   *Reco_Muon_GlbTrk_PtErr;
   vector<short>   *Reco_Muon_GlbTrk_ValidMuonHits;
   vector<float>   *Reco_Muon_GlbTrk_NormChi2;
   vector<char>    *Reco_Muon_BestTrk_Type;
   TClonesArray    *Reco_Muon_BestTrk_Mom;
   TClonesArray    *Reco_Muon_BestTrk_Vertex;
   vector<float>   *Reco_Muon_BestTrk_PtErr;
   vector<float>   *Reco_Muon_BestTrk_dXY;
   vector<float>   *Reco_Muon_BestTrk_dXYErr;
   vector<float>   *Reco_Muon_BestTrk_dZ;
   vector<float>   *Reco_Muon_BestTrk_dZErr;
   vector<float>   *Reco_Muon_IsoPFR03;
   vector<float>   *Reco_Muon_IsoPFR03NoPUCorr;
   vector<float>   *Reco_Muon_IsoPFR04;
   vector<float>   *Reco_Muon_IsoPFR04NoPUCorr;
   vector<float>   *Reco_Muon_EM_Chg_sumR03Pt;
   vector<float>   *Reco_Muon_EM_Chg_sumR04Pt;
   vector<float>   *Reco_Muon_EM_Neu_sumR03Et;
   vector<float>   *Reco_Muon_EM_Neu_sumR04Et;
   vector<float>   *Reco_Muon_Had_Chg_sumR03Pt;
   vector<float>   *Reco_Muon_Had_Chg_sumR04Pt;
   vector<float>   *Reco_Muon_Had_Neu_sumR03Et;
   vector<float>   *Reco_Muon_Had_Neu_sumR04Et;
   vector<float>   *Reco_Muon_Had_PU_sumR03Pt;
   vector<float>   *Reco_Muon_Had_PU_sumR04Pt;
   vector<float>   *Reco_Muon_IsoR03;
   vector<float>   *Reco_Muon_IsoR05;
   vector<float>   *Reco_Muon_Trk_sumR03Pt;
   vector<float>   *Reco_Muon_Trk_sumR05Pt;
   UShort_t        Reco_DiMuon_N;
   TClonesArray    *Reco_DiMuon_Mom;
   vector<char>    *Reco_DiMuon_Charge;
   vector<unsigned char> *Reco_DiMuon_Muon1_Idx;
   vector<unsigned char> *Reco_DiMuon_Muon2_Idx;
   vector<bool>    *Reco_DiMuon_isCowBoy;
   TClonesArray    *Reco_DiMuon_Vertex;
   vector<float>   *Reco_DiMuon_VtxProb;
   vector<float>   *Reco_DiMuon_DCA;
   vector<float>   *Reco_DiMuon_CTau;
   vector<float>   *Reco_DiMuon_CTauErr;
   vector<float>   *Reco_DiMuon_CosAlpha;
   vector<float>   *Reco_DiMuon_MassErr;

   // List of branches
   TBranch        *b_Reco_Muon_N;   //!
   TBranch        *b_Reco_Muon_Mom;   //!
   TBranch        *b_Reco_Muon_Charge;   //!
   TBranch        *b_Reco_Muon_PF_Idx;   //!
   TBranch        *b_Pat_Muon_Trig;   //!
   TBranch        *b_Pat_Muon_dB;   //!
   TBranch        *b_Pat_Muon_dBErr;   //!
   TBranch        *b_Reco_Muon_isPF;   //!
   TBranch        *b_Reco_Muon_isGlobal;   //!
   TBranch        *b_Reco_Muon_isTracker;   //!
   TBranch        *b_Reco_Muon_isStandAlone;   //!
   TBranch        *b_Reco_Muon_isLoose;   //!
   TBranch        *b_Reco_Muon_isMedium;   //!
   TBranch        *b_Reco_Muon_isHighPt;   //!
   TBranch        *b_Reco_Muon_isSoft;   //!
   TBranch        *b_Reco_Muon_isTight;   //!
   TBranch        *b_Reco_Muon_isArbitrated;   //!
   TBranch        *b_Reco_Muon_TrackerArbitrated;   //!
   TBranch        *b_Reco_Muon_GlobalPromptTight;   //!
   TBranch        *b_Reco_Muon_TMLastStationLoose;   //!
   TBranch        *b_Reco_Muon_TMLastStationTight;   //!
   TBranch        *b_Reco_Muon_TM2DCompatibilityLoose;   //!
   TBranch        *b_Reco_Muon_TM2DCompatibilityTight;   //!
   TBranch        *b_Reco_Muon_TMOneStationLoose;   //!
   TBranch        *b_Reco_Muon_TMOneStationTight;   //!
   TBranch        *b_Reco_Muon_GMTkChiCompatibility;   //!
   TBranch        *b_Reco_Muon_GMStaChiCompatibility;   //!
   TBranch        *b_Reco_Muon_GMTkKinkTight;   //!
   TBranch        *b_Reco_Muon_TMLastStationAngLoose;   //!
   TBranch        *b_Reco_Muon_TMLastStationAngTight;   //!
   TBranch        *b_Reco_Muon_TMOneStationAngLoose;   //!
   TBranch        *b_Reco_Muon_TMOneStationAngTight;   //!
   TBranch        *b_Reco_Muon_MatchedStations;   //!
   TBranch        *b_Reco_Muon_Matches;   //!
   TBranch        *b_Reco_Muon_SegmentComp;   //!
   TBranch        *b_Reco_Muon_Chi2Pos;   //!
   TBranch        *b_Reco_Muon_TrkKink;   //!
   TBranch        *b_Reco_Muon_InTrk_Mom;   //!
   TBranch        *b_Reco_Muon_InTrk_PtErr;   //!
   TBranch        *b_Reco_Muon_InTrk_isHighPurity;   //!
   TBranch        *b_Reco_Muon_InTrk_ValidHits;   //!
   TBranch        *b_Reco_Muon_InTrk_LostHits;   //!
   TBranch        *b_Reco_Muon_InTrk_ValidPixHits;   //!
   TBranch        *b_Reco_Muon_InTrk_TrkLayers;   //!
   TBranch        *b_Reco_Muon_InTrk_PixLayers;   //!
   TBranch        *b_Reco_Muon_InTrk_dXY;   //!
   TBranch        *b_Reco_Muon_InTrk_dXYErr;   //!
   TBranch        *b_Reco_Muon_InTrk_dZ;   //!
   TBranch        *b_Reco_Muon_InTrk_dZErr;   //!
   TBranch        *b_Reco_Muon_InTrk_ValFrac;   //!
   TBranch        *b_Reco_Muon_InTrk_NormChi2;   //!
   TBranch        *b_Reco_Muon_GlbTrk_Mom;   //!
   TBranch        *b_Reco_Muon_GlbTrk_PtErr;   //!
   TBranch        *b_Reco_Muon_GlbTrk_ValidMuonHits;   //!
   TBranch        *b_Reco_Muon_GlbTrk_NormChi2;   //!
   TBranch        *b_Reco_Muon_BestTrk_Type;   //!
   TBranch        *b_Reco_Muon_BestTrk_Mom;   //!
   TBranch        *b_Reco_Muon_BestTrk_Vertex;   //!
   TBranch        *b_Reco_Muon_BestTrk_PtErr;   //!
   TBranch        *b_Reco_Muon_BestTrk_dXY;   //!
   TBranch        *b_Reco_Muon_BestTrk_dXYErr;   //!
   TBranch        *b_Reco_Muon_BestTrk_dZ;   //!
   TBranch        *b_Reco_Muon_BestTrk_dZErr;   //!
   TBranch        *b_Reco_Muon_IsoPFR03;   //!
   TBranch        *b_Reco_Muon_IsoPFR03NoPUCorr;   //!
   TBranch        *b_Reco_Muon_IsoPFR04;   //!
   TBranch        *b_Reco_Muon_IsoPFR04NoPUCorr;   //!
   TBranch        *b_Reco_Muon_EM_Chg_sumR03Pt;   //!
   TBranch        *b_Reco_Muon_EM_Chg_sumR04Pt;   //!
   TBranch        *b_Reco_Muon_EM_Neu_sumR03Et;   //!
   TBranch        *b_Reco_Muon_EM_Neu_sumR04Et;   //!
   TBranch        *b_Reco_Muon_Had_Chg_sumR03Pt;   //!
   TBranch        *b_Reco_Muon_Had_Chg_sumR04Pt;   //!
   TBranch        *b_Reco_Muon_Had_Neu_sumR03Et;   //!
   TBranch        *b_Reco_Muon_Had_Neu_sumR04Et;   //!
   TBranch        *b_Reco_Muon_Had_PU_sumR03Pt;   //!
   TBranch        *b_Reco_Muon_Had_PU_sumR04Pt;   //!
   TBranch        *b_Reco_Muon_IsoR03;   //!
   TBranch        *b_Reco_Muon_IsoR05;   //!
   TBranch        *b_Reco_Muon_Trk_sumR03Pt;   //!
   TBranch        *b_Reco_Muon_Trk_sumR05Pt;   //!
   TBranch        *b_Reco_DiMuon_N;   //!
   TBranch        *b_Reco_DiMuon_Mom;   //!
   TBranch        *b_Reco_DiMuon_Charge;   //!
   TBranch        *b_Reco_DiMuon_Muon1_Idx;   //!
   TBranch        *b_Reco_DiMuon_Muon2_Idx;   //!
   TBranch        *b_Reco_DiMuon_isCowBoy;   //!
   TBranch        *b_Reco_DiMuon_Vertex;   //!
   TBranch        *b_Reco_DiMuon_VtxProb;   //!
   TBranch        *b_Reco_DiMuon_DCA;   //!
   TBranch        *b_Reco_DiMuon_CTau;   //!
   TBranch        *b_Reco_DiMuon_CTauErr;   //!
   TBranch        *b_Reco_DiMuon_CosAlpha;   //!
   TBranch        *b_Reco_DiMuon_MassErr;   //!

   muRecoAna(TTree *tree=0);
   virtual ~muRecoAna();
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

#ifdef muRecoAna_cxx
muRecoAna::muRecoAna(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("root://cms-xrd-global.cern.ch//store/group/phys_heavyions/anstahll/EWQAnalysis2017/pPb2016/8160GeV/Data/HiEWQForest_PASingleMuon_pPb_Pbp_8160GeV_20171003.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("root://cms-xrd-global.cern.ch//store/group/phys_heavyions/anstahll/EWQAnalysis2017/pPb2016/8160GeV/Data/HiEWQForest_PASingleMuon_pPb_Pbp_8160GeV_20171003.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("root://cms-xrd-global.cern.ch//store/group/phys_heavyions/anstahll/EWQAnalysis2017/pPb2016/8160GeV/Data/HiEWQForest_PASingleMuon_pPb_Pbp_8160GeV_20171003.root:/muonAna");
      dir->GetObject("Muon_Reco",tree);

   }
   Init(tree);
}

muRecoAna::~muRecoAna()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t muRecoAna::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t muRecoAna::LoadTree(Long64_t entry)
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

void muRecoAna::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Reco_Muon_Mom = 0;
   Reco_Muon_Charge = 0;
   Reco_Muon_PF_Idx = 0;
   Pat_Muon_Trig = 0;
   Pat_Muon_dB = 0;
   Pat_Muon_dBErr = 0;
   Reco_Muon_isPF = 0;
   Reco_Muon_isGlobal = 0;
   Reco_Muon_isTracker = 0;
   Reco_Muon_isStandAlone = 0;
   Reco_Muon_isLoose = 0;
   Reco_Muon_isMedium = 0;
   Reco_Muon_isHighPt = 0;
   Reco_Muon_isSoft = 0;
   Reco_Muon_isTight = 0;
   Reco_Muon_isArbitrated = 0;
   Reco_Muon_TrackerArbitrated = 0;
   Reco_Muon_GlobalPromptTight = 0;
   Reco_Muon_TMLastStationLoose = 0;
   Reco_Muon_TMLastStationTight = 0;
   Reco_Muon_TM2DCompatibilityLoose = 0;
   Reco_Muon_TM2DCompatibilityTight = 0;
   Reco_Muon_TMOneStationLoose = 0;
   Reco_Muon_TMOneStationTight = 0;
   Reco_Muon_GMTkChiCompatibility = 0;
   Reco_Muon_GMStaChiCompatibility = 0;
   Reco_Muon_GMTkKinkTight = 0;
   Reco_Muon_TMLastStationAngLoose = 0;
   Reco_Muon_TMLastStationAngTight = 0;
   Reco_Muon_TMOneStationAngLoose = 0;
   Reco_Muon_TMOneStationAngTight = 0;
   Reco_Muon_MatchedStations = 0;
   Reco_Muon_Matches = 0;
   Reco_Muon_SegmentComp = 0;
   Reco_Muon_Chi2Pos = 0;
   Reco_Muon_TrkKink = 0;
   Reco_Muon_InTrk_Mom = 0;
   Reco_Muon_InTrk_PtErr = 0;
   Reco_Muon_InTrk_isHighPurity = 0;
   Reco_Muon_InTrk_ValidHits = 0;
   Reco_Muon_InTrk_LostHits = 0;
   Reco_Muon_InTrk_ValidPixHits = 0;
   Reco_Muon_InTrk_TrkLayers = 0;
   Reco_Muon_InTrk_PixLayers = 0;
   Reco_Muon_InTrk_dXY = 0;
   Reco_Muon_InTrk_dXYErr = 0;
   Reco_Muon_InTrk_dZ = 0;
   Reco_Muon_InTrk_dZErr = 0;
   Reco_Muon_InTrk_ValFrac = 0;
   Reco_Muon_InTrk_NormChi2 = 0;
   Reco_Muon_GlbTrk_Mom = 0;
   Reco_Muon_GlbTrk_PtErr = 0;
   Reco_Muon_GlbTrk_ValidMuonHits = 0;
   Reco_Muon_GlbTrk_NormChi2 = 0;
   Reco_Muon_BestTrk_Type = 0;
   Reco_Muon_BestTrk_Mom = 0;
   Reco_Muon_BestTrk_Vertex = 0;
   Reco_Muon_BestTrk_PtErr = 0;
   Reco_Muon_BestTrk_dXY = 0;
   Reco_Muon_BestTrk_dXYErr = 0;
   Reco_Muon_BestTrk_dZ = 0;
   Reco_Muon_BestTrk_dZErr = 0;
   Reco_Muon_IsoPFR03 = 0;
   Reco_Muon_IsoPFR03NoPUCorr = 0;
   Reco_Muon_IsoPFR04 = 0;
   Reco_Muon_IsoPFR04NoPUCorr = 0;
   Reco_Muon_EM_Chg_sumR03Pt = 0;
   Reco_Muon_EM_Chg_sumR04Pt = 0;
   Reco_Muon_EM_Neu_sumR03Et = 0;
   Reco_Muon_EM_Neu_sumR04Et = 0;
   Reco_Muon_Had_Chg_sumR03Pt = 0;
   Reco_Muon_Had_Chg_sumR04Pt = 0;
   Reco_Muon_Had_Neu_sumR03Et = 0;
   Reco_Muon_Had_Neu_sumR04Et = 0;
   Reco_Muon_Had_PU_sumR03Pt = 0;
   Reco_Muon_Had_PU_sumR04Pt = 0;
   Reco_Muon_IsoR03 = 0;
   Reco_Muon_IsoR05 = 0;
   Reco_Muon_Trk_sumR03Pt = 0;
   Reco_Muon_Trk_sumR05Pt = 0;
   Reco_DiMuon_Mom = 0;
   Reco_DiMuon_Charge = 0;
   Reco_DiMuon_Muon1_Idx = 0;
   Reco_DiMuon_Muon2_Idx = 0;
   Reco_DiMuon_isCowBoy = 0;
   Reco_DiMuon_Vertex = 0;
   Reco_DiMuon_VtxProb = 0;
   Reco_DiMuon_DCA = 0;
   Reco_DiMuon_CTau = 0;
   Reco_DiMuon_CTauErr = 0;
   Reco_DiMuon_CosAlpha = 0;
   Reco_DiMuon_MassErr = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Reco_Muon_N", &Reco_Muon_N, &b_Reco_Muon_N);
   fChain->SetBranchAddress("Reco_Muon_Mom", &Reco_Muon_Mom, &b_Reco_Muon_Mom);
   fChain->SetBranchAddress("Reco_Muon_Charge", &Reco_Muon_Charge, &b_Reco_Muon_Charge);
   fChain->SetBranchAddress("Reco_Muon_PF_Idx", &Reco_Muon_PF_Idx, &b_Reco_Muon_PF_Idx);
   fChain->SetBranchAddress("Pat_Muon_Trig", &Pat_Muon_Trig, &b_Pat_Muon_Trig);
   fChain->SetBranchAddress("Pat_Muon_dB", &Pat_Muon_dB, &b_Pat_Muon_dB);
   fChain->SetBranchAddress("Pat_Muon_dBErr", &Pat_Muon_dBErr, &b_Pat_Muon_dBErr);
   fChain->SetBranchAddress("Reco_Muon_isPF", &Reco_Muon_isPF, &b_Reco_Muon_isPF);
   fChain->SetBranchAddress("Reco_Muon_isGlobal", &Reco_Muon_isGlobal, &b_Reco_Muon_isGlobal);
   fChain->SetBranchAddress("Reco_Muon_isTracker", &Reco_Muon_isTracker, &b_Reco_Muon_isTracker);
   fChain->SetBranchAddress("Reco_Muon_isStandAlone", &Reco_Muon_isStandAlone, &b_Reco_Muon_isStandAlone);
   fChain->SetBranchAddress("Reco_Muon_isLoose", &Reco_Muon_isLoose, &b_Reco_Muon_isLoose);
   fChain->SetBranchAddress("Reco_Muon_isMedium", &Reco_Muon_isMedium, &b_Reco_Muon_isMedium);
   fChain->SetBranchAddress("Reco_Muon_isHighPt", &Reco_Muon_isHighPt, &b_Reco_Muon_isHighPt);
   fChain->SetBranchAddress("Reco_Muon_isSoft", &Reco_Muon_isSoft, &b_Reco_Muon_isSoft);
   fChain->SetBranchAddress("Reco_Muon_isTight", &Reco_Muon_isTight, &b_Reco_Muon_isTight);
   fChain->SetBranchAddress("Reco_Muon_isArbitrated", &Reco_Muon_isArbitrated, &b_Reco_Muon_isArbitrated);
   fChain->SetBranchAddress("Reco_Muon_TrackerArbitrated", &Reco_Muon_TrackerArbitrated, &b_Reco_Muon_TrackerArbitrated);
   fChain->SetBranchAddress("Reco_Muon_GlobalPromptTight", &Reco_Muon_GlobalPromptTight, &b_Reco_Muon_GlobalPromptTight);
   fChain->SetBranchAddress("Reco_Muon_TMLastStationLoose", &Reco_Muon_TMLastStationLoose, &b_Reco_Muon_TMLastStationLoose);
   fChain->SetBranchAddress("Reco_Muon_TMLastStationTight", &Reco_Muon_TMLastStationTight, &b_Reco_Muon_TMLastStationTight);
   fChain->SetBranchAddress("Reco_Muon_TM2DCompatibilityLoose", &Reco_Muon_TM2DCompatibilityLoose, &b_Reco_Muon_TM2DCompatibilityLoose);
   fChain->SetBranchAddress("Reco_Muon_TM2DCompatibilityTight", &Reco_Muon_TM2DCompatibilityTight, &b_Reco_Muon_TM2DCompatibilityTight);
   fChain->SetBranchAddress("Reco_Muon_TMOneStationLoose", &Reco_Muon_TMOneStationLoose, &b_Reco_Muon_TMOneStationLoose);
   fChain->SetBranchAddress("Reco_Muon_TMOneStationTight", &Reco_Muon_TMOneStationTight, &b_Reco_Muon_TMOneStationTight);
   fChain->SetBranchAddress("Reco_Muon_GMTkChiCompatibility", &Reco_Muon_GMTkChiCompatibility, &b_Reco_Muon_GMTkChiCompatibility);
   fChain->SetBranchAddress("Reco_Muon_GMStaChiCompatibility", &Reco_Muon_GMStaChiCompatibility, &b_Reco_Muon_GMStaChiCompatibility);
   fChain->SetBranchAddress("Reco_Muon_GMTkKinkTight", &Reco_Muon_GMTkKinkTight, &b_Reco_Muon_GMTkKinkTight);
   fChain->SetBranchAddress("Reco_Muon_TMLastStationAngLoose", &Reco_Muon_TMLastStationAngLoose, &b_Reco_Muon_TMLastStationAngLoose);
   fChain->SetBranchAddress("Reco_Muon_TMLastStationAngTight", &Reco_Muon_TMLastStationAngTight, &b_Reco_Muon_TMLastStationAngTight);
   fChain->SetBranchAddress("Reco_Muon_TMOneStationAngLoose", &Reco_Muon_TMOneStationAngLoose, &b_Reco_Muon_TMOneStationAngLoose);
   fChain->SetBranchAddress("Reco_Muon_TMOneStationAngTight", &Reco_Muon_TMOneStationAngTight, &b_Reco_Muon_TMOneStationAngTight);
   fChain->SetBranchAddress("Reco_Muon_MatchedStations", &Reco_Muon_MatchedStations, &b_Reco_Muon_MatchedStations);
   fChain->SetBranchAddress("Reco_Muon_Matches", &Reco_Muon_Matches, &b_Reco_Muon_Matches);
   fChain->SetBranchAddress("Reco_Muon_SegmentComp", &Reco_Muon_SegmentComp, &b_Reco_Muon_SegmentComp);
   fChain->SetBranchAddress("Reco_Muon_Chi2Pos", &Reco_Muon_Chi2Pos, &b_Reco_Muon_Chi2Pos);
   fChain->SetBranchAddress("Reco_Muon_TrkKink", &Reco_Muon_TrkKink, &b_Reco_Muon_TrkKink);
   fChain->SetBranchAddress("Reco_Muon_InTrk_Mom", &Reco_Muon_InTrk_Mom, &b_Reco_Muon_InTrk_Mom);
   fChain->SetBranchAddress("Reco_Muon_InTrk_PtErr", &Reco_Muon_InTrk_PtErr, &b_Reco_Muon_InTrk_PtErr);
   fChain->SetBranchAddress("Reco_Muon_InTrk_isHighPurity", &Reco_Muon_InTrk_isHighPurity, &b_Reco_Muon_InTrk_isHighPurity);
   fChain->SetBranchAddress("Reco_Muon_InTrk_ValidHits", &Reco_Muon_InTrk_ValidHits, &b_Reco_Muon_InTrk_ValidHits);
   fChain->SetBranchAddress("Reco_Muon_InTrk_LostHits", &Reco_Muon_InTrk_LostHits, &b_Reco_Muon_InTrk_LostHits);
   fChain->SetBranchAddress("Reco_Muon_InTrk_ValidPixHits", &Reco_Muon_InTrk_ValidPixHits, &b_Reco_Muon_InTrk_ValidPixHits);
   fChain->SetBranchAddress("Reco_Muon_InTrk_TrkLayers", &Reco_Muon_InTrk_TrkLayers, &b_Reco_Muon_InTrk_TrkLayers);
   fChain->SetBranchAddress("Reco_Muon_InTrk_PixLayers", &Reco_Muon_InTrk_PixLayers, &b_Reco_Muon_InTrk_PixLayers);
   fChain->SetBranchAddress("Reco_Muon_InTrk_dXY", &Reco_Muon_InTrk_dXY, &b_Reco_Muon_InTrk_dXY);
   fChain->SetBranchAddress("Reco_Muon_InTrk_dXYErr", &Reco_Muon_InTrk_dXYErr, &b_Reco_Muon_InTrk_dXYErr);
   fChain->SetBranchAddress("Reco_Muon_InTrk_dZ", &Reco_Muon_InTrk_dZ, &b_Reco_Muon_InTrk_dZ);
   fChain->SetBranchAddress("Reco_Muon_InTrk_dZErr", &Reco_Muon_InTrk_dZErr, &b_Reco_Muon_InTrk_dZErr);
   fChain->SetBranchAddress("Reco_Muon_InTrk_ValFrac", &Reco_Muon_InTrk_ValFrac, &b_Reco_Muon_InTrk_ValFrac);
   fChain->SetBranchAddress("Reco_Muon_InTrk_NormChi2", &Reco_Muon_InTrk_NormChi2, &b_Reco_Muon_InTrk_NormChi2);
   fChain->SetBranchAddress("Reco_Muon_GlbTrk_Mom", &Reco_Muon_GlbTrk_Mom, &b_Reco_Muon_GlbTrk_Mom);
   fChain->SetBranchAddress("Reco_Muon_GlbTrk_PtErr", &Reco_Muon_GlbTrk_PtErr, &b_Reco_Muon_GlbTrk_PtErr);
   fChain->SetBranchAddress("Reco_Muon_GlbTrk_ValidMuonHits", &Reco_Muon_GlbTrk_ValidMuonHits, &b_Reco_Muon_GlbTrk_ValidMuonHits);
   fChain->SetBranchAddress("Reco_Muon_GlbTrk_NormChi2", &Reco_Muon_GlbTrk_NormChi2, &b_Reco_Muon_GlbTrk_NormChi2);
   fChain->SetBranchAddress("Reco_Muon_BestTrk_Type", &Reco_Muon_BestTrk_Type, &b_Reco_Muon_BestTrk_Type);
   fChain->SetBranchAddress("Reco_Muon_BestTrk_Mom", &Reco_Muon_BestTrk_Mom, &b_Reco_Muon_BestTrk_Mom);
   fChain->SetBranchAddress("Reco_Muon_BestTrk_Vertex", &Reco_Muon_BestTrk_Vertex, &b_Reco_Muon_BestTrk_Vertex);
   fChain->SetBranchAddress("Reco_Muon_BestTrk_PtErr", &Reco_Muon_BestTrk_PtErr, &b_Reco_Muon_BestTrk_PtErr);
   fChain->SetBranchAddress("Reco_Muon_BestTrk_dXY", &Reco_Muon_BestTrk_dXY, &b_Reco_Muon_BestTrk_dXY);
   fChain->SetBranchAddress("Reco_Muon_BestTrk_dXYErr", &Reco_Muon_BestTrk_dXYErr, &b_Reco_Muon_BestTrk_dXYErr);
   fChain->SetBranchAddress("Reco_Muon_BestTrk_dZ", &Reco_Muon_BestTrk_dZ, &b_Reco_Muon_BestTrk_dZ);
   fChain->SetBranchAddress("Reco_Muon_BestTrk_dZErr", &Reco_Muon_BestTrk_dZErr, &b_Reco_Muon_BestTrk_dZErr);
   fChain->SetBranchAddress("Reco_Muon_IsoPFR03", &Reco_Muon_IsoPFR03, &b_Reco_Muon_IsoPFR03);
   fChain->SetBranchAddress("Reco_Muon_IsoPFR03NoPUCorr", &Reco_Muon_IsoPFR03NoPUCorr, &b_Reco_Muon_IsoPFR03NoPUCorr);
   fChain->SetBranchAddress("Reco_Muon_IsoPFR04", &Reco_Muon_IsoPFR04, &b_Reco_Muon_IsoPFR04);
   fChain->SetBranchAddress("Reco_Muon_IsoPFR04NoPUCorr", &Reco_Muon_IsoPFR04NoPUCorr, &b_Reco_Muon_IsoPFR04NoPUCorr);
   fChain->SetBranchAddress("Reco_Muon_EM_Chg_sumR03Pt", &Reco_Muon_EM_Chg_sumR03Pt, &b_Reco_Muon_EM_Chg_sumR03Pt);
   fChain->SetBranchAddress("Reco_Muon_EM_Chg_sumR04Pt", &Reco_Muon_EM_Chg_sumR04Pt, &b_Reco_Muon_EM_Chg_sumR04Pt);
   fChain->SetBranchAddress("Reco_Muon_EM_Neu_sumR03Et", &Reco_Muon_EM_Neu_sumR03Et, &b_Reco_Muon_EM_Neu_sumR03Et);
   fChain->SetBranchAddress("Reco_Muon_EM_Neu_sumR04Et", &Reco_Muon_EM_Neu_sumR04Et, &b_Reco_Muon_EM_Neu_sumR04Et);
   fChain->SetBranchAddress("Reco_Muon_Had_Chg_sumR03Pt", &Reco_Muon_Had_Chg_sumR03Pt, &b_Reco_Muon_Had_Chg_sumR03Pt);
   fChain->SetBranchAddress("Reco_Muon_Had_Chg_sumR04Pt", &Reco_Muon_Had_Chg_sumR04Pt, &b_Reco_Muon_Had_Chg_sumR04Pt);
   fChain->SetBranchAddress("Reco_Muon_Had_Neu_sumR03Et", &Reco_Muon_Had_Neu_sumR03Et, &b_Reco_Muon_Had_Neu_sumR03Et);
   fChain->SetBranchAddress("Reco_Muon_Had_Neu_sumR04Et", &Reco_Muon_Had_Neu_sumR04Et, &b_Reco_Muon_Had_Neu_sumR04Et);
   fChain->SetBranchAddress("Reco_Muon_Had_PU_sumR03Pt", &Reco_Muon_Had_PU_sumR03Pt, &b_Reco_Muon_Had_PU_sumR03Pt);
   fChain->SetBranchAddress("Reco_Muon_Had_PU_sumR04Pt", &Reco_Muon_Had_PU_sumR04Pt, &b_Reco_Muon_Had_PU_sumR04Pt);
   fChain->SetBranchAddress("Reco_Muon_IsoR03", &Reco_Muon_IsoR03, &b_Reco_Muon_IsoR03);
   fChain->SetBranchAddress("Reco_Muon_IsoR05", &Reco_Muon_IsoR05, &b_Reco_Muon_IsoR05);
   fChain->SetBranchAddress("Reco_Muon_Trk_sumR03Pt", &Reco_Muon_Trk_sumR03Pt, &b_Reco_Muon_Trk_sumR03Pt);
   fChain->SetBranchAddress("Reco_Muon_Trk_sumR05Pt", &Reco_Muon_Trk_sumR05Pt, &b_Reco_Muon_Trk_sumR05Pt);
   fChain->SetBranchAddress("Reco_DiMuon_N", &Reco_DiMuon_N, &b_Reco_DiMuon_N);
   fChain->SetBranchAddress("Reco_DiMuon_Mom", &Reco_DiMuon_Mom, &b_Reco_DiMuon_Mom);
   fChain->SetBranchAddress("Reco_DiMuon_Charge", &Reco_DiMuon_Charge, &b_Reco_DiMuon_Charge);
   fChain->SetBranchAddress("Reco_DiMuon_Muon1_Idx", &Reco_DiMuon_Muon1_Idx, &b_Reco_DiMuon_Muon1_Idx);
   fChain->SetBranchAddress("Reco_DiMuon_Muon2_Idx", &Reco_DiMuon_Muon2_Idx, &b_Reco_DiMuon_Muon2_Idx);
   fChain->SetBranchAddress("Reco_DiMuon_isCowBoy", &Reco_DiMuon_isCowBoy, &b_Reco_DiMuon_isCowBoy);
   fChain->SetBranchAddress("Reco_DiMuon_Vertex", &Reco_DiMuon_Vertex, &b_Reco_DiMuon_Vertex);
   fChain->SetBranchAddress("Reco_DiMuon_VtxProb", &Reco_DiMuon_VtxProb, &b_Reco_DiMuon_VtxProb);
   fChain->SetBranchAddress("Reco_DiMuon_DCA", &Reco_DiMuon_DCA, &b_Reco_DiMuon_DCA);
   fChain->SetBranchAddress("Reco_DiMuon_CTau", &Reco_DiMuon_CTau, &b_Reco_DiMuon_CTau);
   fChain->SetBranchAddress("Reco_DiMuon_CTauErr", &Reco_DiMuon_CTauErr, &b_Reco_DiMuon_CTauErr);
   fChain->SetBranchAddress("Reco_DiMuon_CosAlpha", &Reco_DiMuon_CosAlpha, &b_Reco_DiMuon_CosAlpha);
   fChain->SetBranchAddress("Reco_DiMuon_MassErr", &Reco_DiMuon_MassErr, &b_Reco_DiMuon_MassErr);
   Notify();
}

Bool_t muRecoAna::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void muRecoAna::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t muRecoAna::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

float muRecoAna::Get(TString varname, int idx) {
   if (varname == "Reco_Muon_Charge") return Reco_Muon_Charge->at(idx);
   if (varname == "Reco_Muon_PF_Idx") return Reco_Muon_PF_Idx->at(idx);
   if (varname == "Reco_Muon_isPF") return Reco_Muon_isPF->at(idx);
   if (varname == "Reco_Muon_isGlobal") return Reco_Muon_isGlobal->at(idx);
   if (varname == "Reco_Muon_isTracker") return Reco_Muon_isTracker->at(idx);
   if (varname == "Reco_Muon_isStandAlone") return Reco_Muon_isStandAlone->at(idx);
   if (varname == "Reco_Muon_isLoose") return Reco_Muon_isLoose->at(idx);
   if (varname == "Reco_Muon_isMedium") return Reco_Muon_isMedium->at(idx);
   if (varname == "Reco_Muon_isHighPt") return Reco_Muon_isHighPt->at(idx);
   if (varname == "Reco_Muon_isSoft") return Reco_Muon_isSoft->at(idx);
   if (varname == "Reco_Muon_isTight") return Reco_Muon_isTight->at(idx);
   if (varname == "Reco_Muon_isArbitrated") return Reco_Muon_isArbitrated->at(idx);
   if (varname == "Reco_Muon_TrackerArbitrated") return Reco_Muon_TrackerArbitrated->at(idx);
   if (varname == "Reco_Muon_GlobalPromptTight") return Reco_Muon_GlobalPromptTight->at(idx);
   if (varname == "Reco_Muon_TMLastStationLoose") return Reco_Muon_TMLastStationLoose->at(idx);
   if (varname == "Reco_Muon_TMLastStationTight") return Reco_Muon_TMLastStationTight->at(idx);
   if (varname == "Reco_Muon_TM2DCompatibilityLoose") return Reco_Muon_TM2DCompatibilityLoose->at(idx);
   if (varname == "Reco_Muon_TM2DCompatibilityTight") return Reco_Muon_TM2DCompatibilityTight->at(idx);
   if (varname == "Reco_Muon_TMOneStationLoose") return Reco_Muon_TMOneStationLoose->at(idx);
   if (varname == "Reco_Muon_TMOneStationTight") return Reco_Muon_TMOneStationTight->at(idx);
   if (varname == "Reco_Muon_GMTkChiCompatibility") return Reco_Muon_GMTkChiCompatibility->at(idx);
   if (varname == "Reco_Muon_GMStaChiCompatibility") return Reco_Muon_GMStaChiCompatibility->at(idx);
   if (varname == "Reco_Muon_GMTkKinkTight") return Reco_Muon_GMTkKinkTight->at(idx);
   if (varname == "Reco_Muon_TMLastStationAngLoose") return Reco_Muon_TMLastStationAngLoose->at(idx);
   if (varname == "Reco_Muon_TMLastStationAngTight") return Reco_Muon_TMLastStationAngTight->at(idx);
   if (varname == "Reco_Muon_TMOneStationAngLoose") return Reco_Muon_TMOneStationAngLoose->at(idx);
   if (varname == "Reco_Muon_TMOneStationAngTight") return Reco_Muon_TMOneStationAngTight->at(idx);
   if (varname == "Reco_Muon_MatchedStations") return Reco_Muon_MatchedStations->at(idx);
   if (varname == "Reco_Muon_Matches") return Reco_Muon_Matches->at(idx);
   if (varname == "Reco_Muon_SegmentComp") return Reco_Muon_SegmentComp->at(idx);
   if (varname == "Reco_Muon_Chi2Pos") return Reco_Muon_Chi2Pos->at(idx);
   if (varname == "Reco_Muon_TrkKink") return Reco_Muon_TrkKink->at(idx);
   if (varname == "Reco_Muon_InTrk_PtErr") return Reco_Muon_InTrk_PtErr->at(idx);
   if (varname == "Reco_Muon_InTrk_isHighPurity") return Reco_Muon_InTrk_isHighPurity->at(idx);
   if (varname == "Reco_Muon_InTrk_ValidHits") return Reco_Muon_InTrk_ValidHits->at(idx);
   if (varname == "Reco_Muon_InTrk_LostHits") return Reco_Muon_InTrk_LostHits->at(idx);
   if (varname == "Reco_Muon_InTrk_ValidPixHits") return Reco_Muon_InTrk_ValidPixHits->at(idx);
   if (varname == "Reco_Muon_InTrk_TrkLayers") return Reco_Muon_InTrk_TrkLayers->at(idx);
   if (varname == "Reco_Muon_InTrk_PixLayers") return Reco_Muon_InTrk_PixLayers->at(idx);
   if (varname == "Reco_Muon_InTrk_dXY") return Reco_Muon_InTrk_dXY->at(idx);
   if (varname == "Reco_Muon_InTrk_dXYErr") return Reco_Muon_InTrk_dXYErr->at(idx);
   if (varname == "Reco_Muon_InTrk_dZ") return Reco_Muon_InTrk_dZ->at(idx);
   if (varname == "Reco_Muon_InTrk_dZErr") return Reco_Muon_InTrk_dZErr->at(idx);
   if (varname == "Reco_Muon_InTrk_ValFrac") return Reco_Muon_InTrk_ValFrac->at(idx);
   if (varname == "Reco_Muon_InTrk_NormChi2") return Reco_Muon_InTrk_NormChi2->at(idx);
   if (varname == "Reco_Muon_GlbTrk_PtErr") return Reco_Muon_GlbTrk_PtErr->at(idx);
   if (varname == "Reco_Muon_GlbTrk_ValidMuonHits") return Reco_Muon_GlbTrk_ValidMuonHits->at(idx);
   if (varname == "Reco_Muon_GlbTrk_NormChi2") return Reco_Muon_GlbTrk_NormChi2->at(idx);
   if (varname == "Reco_Muon_BestTrk_Type") return Reco_Muon_BestTrk_Type->at(idx);
   if (varname == "Reco_Muon_BestTrk_PtErr") return Reco_Muon_BestTrk_PtErr->at(idx);
   if (varname == "Reco_Muon_BestTrk_dXY") return Reco_Muon_BestTrk_dXY->at(idx);
   if (varname == "Reco_Muon_BestTrk_dXYErr") return Reco_Muon_BestTrk_dXYErr->at(idx);
   if (varname == "Reco_Muon_BestTrk_dZ") return Reco_Muon_BestTrk_dZ->at(idx);
   if (varname == "Reco_Muon_BestTrk_dZErr") return Reco_Muon_BestTrk_dZErr->at(idx);
   if (varname == "Reco_Muon_IsoPFR03") return Reco_Muon_IsoPFR03->at(idx);
   if (varname == "Reco_Muon_IsoPFR03NoPUCorr") return Reco_Muon_IsoPFR03NoPUCorr->at(idx);
   if (varname == "Reco_Muon_IsoPFR04") return Reco_Muon_IsoPFR04->at(idx);
   if (varname == "Reco_Muon_IsoPFR04NoPUCorr") return Reco_Muon_IsoPFR04NoPUCorr->at(idx);
   if (varname == "Reco_Muon_EM_Chg_sumR03Pt") return Reco_Muon_EM_Chg_sumR03Pt->at(idx);
   if (varname == "Reco_Muon_EM_Chg_sumR04Pt") return Reco_Muon_EM_Chg_sumR04Pt->at(idx);
   if (varname == "Reco_Muon_EM_Neu_sumR03Et") return Reco_Muon_EM_Neu_sumR03Et->at(idx);
   if (varname == "Reco_Muon_EM_Neu_sumR04Et") return Reco_Muon_EM_Neu_sumR04Et->at(idx);
   if (varname == "Reco_Muon_Had_Chg_sumR03Pt") return Reco_Muon_Had_Chg_sumR03Pt->at(idx);
   if (varname == "Reco_Muon_Had_Chg_sumR04Pt") return Reco_Muon_Had_Chg_sumR04Pt->at(idx);
   if (varname == "Reco_Muon_Had_Neu_sumR03Et") return Reco_Muon_Had_Neu_sumR03Et->at(idx);
   if (varname == "Reco_Muon_Had_Neu_sumR04Et") return Reco_Muon_Had_Neu_sumR04Et->at(idx);
   if (varname == "Reco_Muon_Had_PU_sumR03Pt") return Reco_Muon_Had_PU_sumR03Pt->at(idx);
   if (varname == "Reco_Muon_Had_PU_sumR04Pt") return Reco_Muon_Had_PU_sumR04Pt->at(idx);
   if (varname == "Reco_Muon_IsoR03") return Reco_Muon_IsoR03->at(idx);
   if (varname == "Reco_Muon_IsoR05") return Reco_Muon_IsoR05->at(idx);
   if (varname == "Reco_Muon_Trk_sumR03Pt") return Reco_Muon_Trk_sumR03Pt->at(idx);
   if (varname == "Reco_Muon_Trk_sumR05Pt") return Reco_Muon_Trk_sumR05Pt->at(idx);
   if (varname == "Reco_DiMuon_Charge") return Reco_DiMuon_Charge->at(idx);
   if (varname == "Reco_DiMuon_Muon1_Idx") return Reco_DiMuon_Muon1_Idx->at(idx);
   if (varname == "Reco_DiMuon_Muon2_Idx") return Reco_DiMuon_Muon2_Idx->at(idx);
   if (varname == "Reco_DiMuon_isCowBoy") return Reco_DiMuon_isCowBoy->at(idx);
   if (varname == "Reco_DiMuon_VtxProb") return Reco_DiMuon_VtxProb->at(idx);
   if (varname == "Reco_DiMuon_DCA") return Reco_DiMuon_DCA->at(idx);
   if (varname == "Reco_DiMuon_MassErr") return Reco_DiMuon_MassErr->at(idx);
   return -1e99;
}
#endif // #ifdef muRecoAna_cxx

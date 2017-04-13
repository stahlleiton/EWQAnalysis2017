#include <iostream>
#include <vector>

#include <TSystem.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1D.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLine.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>

#include "HiMuonTree.h"

using namespace std;


class MuIDTest {
  public:
    MuIDTest(const bool _isMC, const string _inputFile, const string _outputFile, const double _muPtCut);
    ~MuIDTest();
    void init();
    void variableTest();

  private:
    bool isValidZmumu(const TLorentzVector *di, const TLorentzVector *mu1, const TLorentzVector *mu2, const double muptcut);
    void getEfficiency();
    void setHistStyle(TH1 &h);
    
    bool isMC;
    string inputFile, outputFile;
    double muPtCut;

    TFile *fout;

    // Check tested entries:
    // All dimuons
    TH1D hPt, hMass;
    TH1D hDxy, hDz;
    // [0] all Z, [1] Z with TightMuonID passed, [2] efficiency
    // [3] Z with MediumMuonID passed, [4] efficiency
    TH1D h_Pt[5];
    TH1D h_Mass[5];

    // Check muon ID variables
    // [0]: TightMuonID selection criteria - denominator (no cuts applied)
    // [1]: TightMuonID selection criteria - all but this one excluded
    // [2]: [1]/[0]
    // [3]: TightMuonID selection criteria - all excluded except this one
    // [4]: [3]/[0]
    TH1D h_isTight[5]; // [0]: no cuts, [1]: TightMuonID applied
    TH1D h_isGlobal[5];
    TH1D h_isPF[5];
    TH1D h_GlbTrk_NormChi2[5];
    TH1D h_GlbTrk_ValidMuonHits[5];
    TH1D h_MatchedStations[5];
    TH1D h_BestTrk_dXY[5];
    TH1D h_BestTrk_dZ[5];
    TH1D h_InTrk_ValidPixHits[5];
    TH1D h_InTrk_TrkLayers[5];

    // [0]: MediumMuonID selection criteria - denominator (no cuts applied except isPF)
    // [1]: MediumMuonID selection criteria - all but this one excluded
    // [2]: [1]/[0]
    // [3]: MediumMuonID selection criteria - all excluded except this one
    // [4]: [3]/[0]
    TH1D h_isMedium[5];
    TH1D h_isLoose[5];
    TH1D h_isGlobal_Medium[5];
    TH1D h_GlbTrk_NormChi2_Medium[5];
    TH1D h_Chi2Pos[5];
    TH1D h_TrkKink[5];
    TH1D h_SegmentComp[5];
    int num_SegmentComp[5];
    TH1D h_InTrk_ValFrac[5];
};

MuIDTest::MuIDTest(const bool _isMC, const string _inputFile, const string _outputFile, const double _muPtCut) {
  isMC = _isMC;
  inputFile = _inputFile;
  outputFile = _outputFile;
  muPtCut = _muPtCut;
  cout << "MuIDTest::MuIDTest() " << isMC << " " << inputFile << " " << outputFile << " " << muPtCut << endl;
}

MuIDTest::~MuIDTest() {
}

bool MuIDTest::isValidZmumu(const TLorentzVector *dimu, const TLorentzVector *mu1, const TLorentzVector *mu2, const double muptcut) {
  bool result = false;
  if (dimu->M()>=80 && dimu->M()<100 && TMath::Abs(dimu->Rapidity())<2.4 &&
      TMath::Abs(mu1->Eta())<2.4 && TMath::Abs(mu2->Eta())<2.4 &&
      TMath::Abs(mu1->Pt())>=muptcut && TMath::Abs(mu2->Pt())>=muptcut) {
    result = true;
  }
  return result;
}

void MuIDTest::setHistStyle(TH1 &h) {
  if (isMC) {
    h.SetLineColor(kRed);
    h.SetFillColor(kRed);
    h.SetFillStyle(3445);
    h.SetLineWidth(1);
  } else {
    h.SetLineColor(kGreen+1);
    h.SetFillColor(kGreen);
    h.SetFillStyle(0);
    h.SetLineWidth(2);
  }
}

void MuIDTest::init() {
  cout << "MuIDTest::init() " << endl;
  fout = new TFile(outputFile.c_str(),"recreate");
  TH1::SetDefaultSumw2();

  hPt = TH1D("hPt",";#mu^{+}#mu^{-} p_{T} [GeV/c];Counts / (0.2 GeV/c)",500,0,100);
  hMass = TH1D("hMass",";#mu^{+}#mu^{-} mass [GeV/c^{2}];Counts / (0.25 GeV/c^{2})",440,0,110);
  hDxy = TH1D("hDxy",";#mu d_{XY} [cm];Counts / (0.001 cm)",1500,-3,3);
  hDz = TH1D("hDz",";#mu d_{Z} [cm];Counts / (0.001 cm)",1500,-3,3);
  setHistStyle(hPt);
  setHistStyle(hMass);
  setHistStyle(hDxy);
  setHistStyle(hDz);
  for (int i=0; i<5; i++) {
    h_Pt[i] = TH1D(Form("h_Pt_%d",i),";#mu^{+}#mu^{-} p_{T} [GeV/c];",100,0,100);
    h_Mass[i] = TH1D(Form("h_Mass_%d",i),";#mu^{+}#mu^{-} mass [GeV/c^{2}];",80,70,110);
    setHistStyle(h_Pt[i]);
    setHistStyle(h_Mass[i]);
  }
  for (int i=0; i<5; i++) {
    h_isTight[i] = TH1D(Form("h_isTight_%d",i),";isTight;",2,0,2);
    h_isGlobal[i] = TH1D(Form("h_isGlobal_%d",i),";isGlobal;",2,0,2);
    h_isPF[i] = TH1D(Form("h_isPF_%d",i),";isPF;",2,0,2);
    h_GlbTrk_NormChi2[i] = TH1D(Form("h_GlbTrk_NormChi2_%d",i),";GlbTrk_NormChi2;",100,0,25);
    h_GlbTrk_ValidMuonHits[i] = TH1D(Form("h_GlbTrk_ValidMuonHits_%d",i),";GlbTrk_ValidMuonHits;",60,0,60);
    h_MatchedStations[i] = TH1D(Form("h_MatchedStations_%d",i),";MatchedStations;",7,0,7);
    h_BestTrk_dXY[i] = TH1D(Form("h_BestTrk_dXY_%d",i),";BestTrk_dXY;",1500,-1,1);
    h_BestTrk_dZ[i] = TH1D(Form("h_BestTrk_dZ_%d",i),";BestTrk_dZ;",1500,-1,1);
    h_InTrk_ValidPixHits[i] = TH1D(Form("h_InTrk_ValidPixHits_%d",i),";InTrk_ValidPixHits;",10,0,10);
    h_InTrk_TrkLayers[i] = TH1D(Form("h_InTrk_TrkLayers_%d",i),";InTrk_TrkLayers;",20,0,20);
    setHistStyle(h_isTight[i]);
    setHistStyle(h_isGlobal[i]);
    setHistStyle(h_isPF[i]);
    setHistStyle(h_GlbTrk_NormChi2[i]);
    setHistStyle(h_GlbTrk_ValidMuonHits[i]);
    setHistStyle(h_MatchedStations[i]);
    setHistStyle(h_BestTrk_dXY[i]);
    setHistStyle(h_BestTrk_dZ[i]);
    setHistStyle(h_InTrk_ValidPixHits[i]);
    setHistStyle(h_InTrk_TrkLayers[i]);

    h_isMedium[i] = TH1D(Form("h_isMedium_%d",i),";isMedium;",2,0,2);
    h_isLoose[i] = TH1D(Form("h_isLoose_%d",i),";isLoose;",2,0,2);
    h_isGlobal_Medium[i] = TH1D(Form("h_isGlobal_Medium_%d",i),";isGlobal_Medium;",2,0,2);
    h_GlbTrk_NormChi2_Medium[i] = TH1D(Form("h_GlbTrk_NormChi2_Medium_%d",i),";GlbTrk_NormChi2_Medium;",100,0,25);
    h_Chi2Pos[i] = TH1D(Form("h_Chi2Pos_%d",i),";Chi2Pos;",50,0,50);
    h_TrkKink[i] = TH1D(Form("h_TrkKink_%d",i),";TrkKink;",50,0,50);
    h_SegmentComp[i] = TH1D(Form("h_SegmentComp_%d",i),";;",100,0,1);
    h_InTrk_ValFrac[i] = TH1D(Form("h_InTrk_ValFrac_%d",i),";InTrk_ValFrac;",100,0,2);
    setHistStyle(h_isMedium[i]);
    setHistStyle(h_isLoose[i]);
    setHistStyle(h_isGlobal_Medium[i]);
    setHistStyle(h_GlbTrk_NormChi2_Medium[i]);
    setHistStyle(h_Chi2Pos[i]);
    setHistStyle(h_TrkKink[i]);
    setHistStyle(h_SegmentComp[i]);
    setHistStyle(h_InTrk_ValFrac[i]);

    num_SegmentComp[i]=0;
  }
}

void MuIDTest::getEfficiency() {
  h_Pt[2].Divide(&h_Pt[1],&h_Pt[0]);
  h_Mass[2].Divide(&h_Mass[1],&h_Mass[0]);
  h_Pt[4].Divide(&h_Pt[3],&h_Pt[0]);
  h_Mass[4].Divide(&h_Mass[3],&h_Mass[0]);
  
  h_isTight[2].Divide(&h_isTight[1],&h_isTight[0]);
  h_isGlobal[2].Divide(&h_isGlobal[1],&h_isGlobal[0]);
  h_isPF[2].Divide(&h_isPF[1],&h_isPF[0]);
  h_GlbTrk_NormChi2[2].Divide(&h_GlbTrk_NormChi2[1],&h_GlbTrk_NormChi2[0]);
  h_GlbTrk_ValidMuonHits[2].Divide(&h_GlbTrk_ValidMuonHits[1],&h_GlbTrk_ValidMuonHits[0]);
  h_MatchedStations[2].Divide(&h_MatchedStations[1],&h_MatchedStations[0]);
  h_BestTrk_dXY[2].Divide(&h_BestTrk_dXY[1],&h_BestTrk_dXY[0]);
  h_BestTrk_dZ[2].Divide(&h_BestTrk_dZ[1],&h_BestTrk_dZ[0]);
  h_InTrk_ValidPixHits[2].Divide(&h_InTrk_ValidPixHits[1],&h_InTrk_ValidPixHits[0]);
  h_InTrk_TrkLayers[2].Divide(&h_InTrk_TrkLayers[1],&h_InTrk_TrkLayers[0]);
 
  h_isTight[4].Divide(&h_isTight[3],&h_isTight[0]);
  h_isGlobal[4].Divide(&h_isGlobal[3],&h_isGlobal[0]);
  h_isPF[4].Divide(&h_isPF[3],&h_isPF[0]);
  h_GlbTrk_NormChi2[4].Divide(&h_GlbTrk_NormChi2[3],&h_GlbTrk_NormChi2[0]);
  h_GlbTrk_ValidMuonHits[4].Divide(&h_GlbTrk_ValidMuonHits[3],&h_GlbTrk_ValidMuonHits[0]);
  h_MatchedStations[4].Divide(&h_MatchedStations[3],&h_MatchedStations[0]);
  h_BestTrk_dXY[4].Divide(&h_BestTrk_dXY[3],&h_BestTrk_dXY[0]);
  h_BestTrk_dZ[4].Divide(&h_BestTrk_dZ[3],&h_BestTrk_dZ[0]);
  h_InTrk_ValidPixHits[4].Divide(&h_InTrk_ValidPixHits[3],&h_InTrk_ValidPixHits[0]);
  h_InTrk_TrkLayers[4].Divide(&h_InTrk_TrkLayers[3],&h_InTrk_TrkLayers[0]);

  h_isMedium[2].Divide(&h_isMedium[1],&h_isMedium[0]);
  h_isLoose[2].Divide(&h_isLoose[1],&h_isLoose[0]);
  h_isGlobal_Medium[2].Divide(&h_isGlobal_Medium[1],&h_isGlobal_Medium[0]);
  h_GlbTrk_NormChi2_Medium[2].Divide(&h_GlbTrk_NormChi2_Medium[1],&h_GlbTrk_NormChi2_Medium[0]);
  h_Chi2Pos[2].Divide(&h_Chi2Pos[1],&h_Chi2Pos[0]);
  h_TrkKink[2].Divide(&h_TrkKink[1],&h_TrkKink[0]);
  h_SegmentComp[2].Divide(&h_SegmentComp[1],&h_SegmentComp[0]);
  num_SegmentComp[2] = num_SegmentComp[1]/num_SegmentComp[0];
  h_InTrk_ValFrac[2].Divide(&h_InTrk_ValFrac[1],&h_InTrk_ValFrac[0]);

  h_isMedium[4].Divide(&h_isMedium[3],&h_isMedium[0]);
  h_isLoose[4].Divide(&h_isLoose[3],&h_isLoose[0]);
  h_isGlobal_Medium[4].Divide(&h_isGlobal_Medium[3],&h_isGlobal_Medium[0]);
  h_GlbTrk_NormChi2_Medium[4].Divide(&h_GlbTrk_NormChi2_Medium[3],&h_GlbTrk_NormChi2_Medium[0]);
  h_Chi2Pos[4].Divide(&h_Chi2Pos[3],&h_Chi2Pos[0]);
  h_TrkKink[4].Divide(&h_TrkKink[3],&h_TrkKink[0]);
  h_SegmentComp[4].Divide(&h_SegmentComp[3],&h_SegmentComp[0]);
  num_SegmentComp[4] = num_SegmentComp[3]/num_SegmentComp[0];
  h_InTrk_ValFrac[4].Divide(&h_InTrk_ValFrac[3],&h_InTrk_ValFrac[0]);
  cout << "MuIDTest::getEfficiency() done" << endl;
}

void MuIDTest::variableTest() {
  cout << "MuIDTest::variableTest() " << endl;
  // Load input file and trees
  HiMuonTree tree = HiMuonTree();
  tree.GetTree(inputFile.c_str());
  cout << "MuIDTest::variableTest() after HiMuonTree GetTree " << endl;

  Long64_t nentries = tree.GetEntries();
  
  for (Long64_t evt=0; evt<nentries; evt++) {
    if (tree.GetEntry(evt)<0) break;
    
    // Pick up single mu event for dxy, dz test
    TClonesArray Reco_Muon_Mom = tree.Reco_Muon_Mom();
    unsigned int nmu = Reco_Muon_Mom.GetEntries();
    vector<float> v_BestTrk_dXY = tree.Reco_Muon_BestTrk_dXY();
    vector<float> v_BestTrk_dZ = tree.Reco_Muon_BestTrk_dZ();

    for (unsigned int imu=0; imu<nmu; imu++) {
      hDxy.Fill(v_BestTrk_dXY.at(imu));
      hDz.Fill(v_BestTrk_dZ.at(imu));
    }
    
    // Pick up Z->mumu event
    TClonesArray Reco_DiMuon_Mom = tree.Reco_DiMuon_Mom();
    vector<char> Reco_DiMuon_Charge = tree.Reco_DiMuon_Charge();
    vector<float> Reco_DiMuon_VtxProb = tree.Reco_DiMuon_VtxProb();
    unsigned int ndimu = Reco_DiMuon_Mom.GetEntries();
    
    for (unsigned int idimu=0; idimu<ndimu; idimu++) {
      TLorentzVector *dimu = (TLorentzVector *)Reco_DiMuon_Mom.At(idimu);
      hMass.Fill(dimu->M());
      hPt.Fill(dimu->Pt());

      // Track down who are Z's daughters
      vector<unsigned short> Muon1_Idx = tree.Reco_DiMuon_Muon1_Idx();
      vector<unsigned short> Muon2_Idx = tree.Reco_DiMuon_Muon2_Idx();
      unsigned short muIdx1 = Muon1_Idx.at(idimu);
      unsigned short muIdx2 = Muon2_Idx.at(idimu);
      TLorentzVector *mu1 = (TLorentzVector *)Reco_Muon_Mom.At(muIdx1);
      TLorentzVector *mu2 = (TLorentzVector *)Reco_Muon_Mom.At(muIdx2);
      
      // Z->mumu event is selected within 80-100 GeV/c2 mass range
      // + daughter muons within certain eta & pT ragions
      if (isValidZmumu(dimu, mu1, mu2, muPtCut) && Reco_DiMuon_Charge.at(idimu)==0 && Reco_DiMuon_VtxProb.at(idimu)>0.01) {
        h_Mass[0].Fill(dimu->M());
        h_Pt[0].Fill(dimu->Pt());
        
        // Consider only PF muons for study
        vector<bool>  v_isPF = tree.Reco_Muon_isPF();
        if (!v_isPF.at(muIdx1) && !v_isPF.at(muIdx2)) continue;

        vector<bool>  v_isTight = tree.Reco_Muon_isTight();
        // Load conditions of isTightMuon()
        vector<bool>  v_isGlobal = tree.Reco_Muon_isGlobal();
        vector<float> v_GlbTrk_NormChi2 = tree.Reco_Muon_GlbTrk_NormChi2();
        vector<short> v_GlbTrk_ValidMuonHits = tree.Reco_Muon_GlbTrk_ValidMuonHits();
        vector<short> v_MatchedStations = tree.Reco_Muon_MatchedStations();
//        vector<float> v_BestTrk_dXY = tree.Reco_Muon_BestTrk_dXY();
//        vector<float> v_BestTrk_dZ = tree.Reco_Muon_BestTrk_dZ();
        vector<short> v_InTrk_ValidPixHits = tree.Reco_Muon_InTrk_ValidPixHits();
        vector<short> v_InTrk_TrkLayers = tree.Reco_Muon_InTrk_TrkLayers();
        
        vector<bool>  v_isMedium = tree.Reco_Muon_isMedium();
        // Load conditions of isMediumMuon()
        vector<bool>  v_isLoose = tree.Reco_Muon_isLoose();
        vector<float> v_Chi2Pos = tree.Reco_Muon_Chi2Pos();
        vector<float> v_TrkKink = tree.Reco_Muon_TrkKink();
        vector<float> v_SegmentComp = tree.Reco_Muon_SegmentComp();
        vector<float> v_InTrk_ValFrac = tree.Reco_Muon_InTrk_ValFrac();

        // Test isTight cuts
        bool isTight1 = v_isTight.at(muIdx1);
        bool isTight2 = v_isTight.at(muIdx2);
        // Pick up the best muon as a tag muon
        unsigned short muIdx = 0;
        if ( (isTight1 && !isTight2) || (isTight1 && isTight2 && mu1->Pt() > mu2->Pt()) )
          muIdx = muIdx2; //muIdx1 is the best quality -> probe is muIdx2
        else if ( (!isTight1 && isTight2) || (isTight1 && isTight2 && mu1->Pt() < mu2->Pt()) )
          muIdx = muIdx2; //muIdx2 is the best quality -> probe is muIdx1
        else //all muons aren't good
          muIdx = 9999;
        
        if (muIdx!=9999) { // Don't test if non of muons in the pair are good
          bool isTight    = v_isTight.at(muIdx);
          bool isGlobal   = v_isGlobal.at(muIdx);
          bool isPF       = v_isPF.at(muIdx);
          float GlbTrk_NormChi2      = v_GlbTrk_NormChi2.at(muIdx);
          short GlbTrk_ValidMuonHits = v_GlbTrk_ValidMuonHits.at(muIdx);
          short MatchedStations      = v_MatchedStations.at(muIdx);
          float BestTrk_dXY          = v_BestTrk_dXY.at(muIdx);
          float BestTrk_dZ           = v_BestTrk_dZ.at(muIdx);
          short InTrk_ValidPixHits   = v_InTrk_ValidPixHits.at(muIdx);
          short InTrk_TrkLayers      = v_InTrk_TrkLayers.at(muIdx);
          
          // Fill denominator histograms (no cuts applied except isPF)
          h_isTight[0].Fill(isTight);
          h_isGlobal[0].Fill(isGlobal);
          h_isPF[0].Fill(isPF);
          h_GlbTrk_NormChi2[0].Fill(GlbTrk_NormChi2);
          h_GlbTrk_ValidMuonHits[0].Fill(GlbTrk_ValidMuonHits);
          h_MatchedStations[0].Fill(MatchedStations);
          h_BestTrk_dXY[0].Fill(BestTrk_dXY);
          h_BestTrk_dZ[0].Fill(BestTrk_dZ);
          h_InTrk_ValidPixHits[0].Fill(InTrk_ValidPixHits);
          h_InTrk_TrkLayers[0].Fill(InTrk_TrkLayers);

          if (isTight1 && isTight2) {
            h_Mass[1].Fill(dimu->M());
            h_Pt[1].Fill(dimu->Pt());
          }

          bool pass_isGlobal=false, pass_isPF=false;
          bool pass_GlbTrk_NormChi2=false, pass_GlbTrk_ValidMuonHits=false;
          bool pass_MatchedStations=false, pass_BestTrk_dXY=false, pass_BestTrk_dZ=false;
          bool pass_InTrk_ValidPixHits=false, pass_InTrk_TrkLayers=false;

          if (isGlobal) pass_isGlobal = true;
          if (isPF) pass_isPF = true;
          if (GlbTrk_NormChi2<10) pass_GlbTrk_NormChi2 = true;
          if (GlbTrk_ValidMuonHits>0) pass_GlbTrk_ValidMuonHits = true;
          if (MatchedStations>1) pass_MatchedStations = true;
          if (TMath::Abs(BestTrk_dXY)<0.2) pass_BestTrk_dXY = true;
          if (TMath::Abs(BestTrk_dZ)<0.5) pass_BestTrk_dZ = true;
          if (InTrk_ValidPixHits>0) pass_InTrk_ValidPixHits = true;
          if (InTrk_TrkLayers>5) pass_InTrk_TrkLayers = true;

          // Fill numerator histograms (TightMuonID cuts applied except that cut)
          if (isGlobal && isPF && GlbTrk_NormChi2 && GlbTrk_ValidMuonHits && MatchedStations && BestTrk_dXY && BestTrk_dZ && InTrk_ValidPixHits && InTrk_TrkLayers) 
            h_isTight[1].Fill(isTight);
          
          if (isPF && GlbTrk_NormChi2 && GlbTrk_ValidMuonHits && MatchedStations && BestTrk_dXY && BestTrk_dZ && InTrk_ValidPixHits && InTrk_TrkLayers) 
            h_isGlobal[1].Fill(isGlobal);
          
          if (isGlobal && GlbTrk_NormChi2 && GlbTrk_ValidMuonHits && MatchedStations && BestTrk_dXY && BestTrk_dZ && InTrk_ValidPixHits && InTrk_TrkLayers) 
            h_isPF[1].Fill(isPF);
          
          if (isGlobal && isPF && GlbTrk_ValidMuonHits && MatchedStations && BestTrk_dXY && BestTrk_dZ && InTrk_ValidPixHits && InTrk_TrkLayers) 
            h_GlbTrk_NormChi2[1].Fill(GlbTrk_NormChi2);
          
          if (isGlobal && isPF && GlbTrk_NormChi2 && MatchedStations && BestTrk_dXY && BestTrk_dZ && InTrk_ValidPixHits && InTrk_TrkLayers) 
            h_GlbTrk_ValidMuonHits[1].Fill(GlbTrk_ValidMuonHits);
          
          if (isGlobal && isPF && GlbTrk_NormChi2 && GlbTrk_ValidMuonHits && BestTrk_dXY && BestTrk_dZ && InTrk_ValidPixHits && InTrk_TrkLayers) 
            h_MatchedStations[1].Fill(MatchedStations);
          
          if (isGlobal && isPF && GlbTrk_NormChi2 && GlbTrk_ValidMuonHits && MatchedStations && BestTrk_dZ && InTrk_ValidPixHits && InTrk_TrkLayers) 
            h_BestTrk_dXY[1].Fill(BestTrk_dXY);
          
          if (isGlobal && isPF && GlbTrk_NormChi2 && GlbTrk_ValidMuonHits && MatchedStations && BestTrk_dXY && InTrk_ValidPixHits && InTrk_TrkLayers) 
            h_BestTrk_dZ[1].Fill(BestTrk_dZ);
          
          if (isGlobal && isPF && GlbTrk_NormChi2 && GlbTrk_ValidMuonHits && MatchedStations && BestTrk_dXY && BestTrk_dZ && InTrk_TrkLayers) 
            h_InTrk_ValidPixHits[1].Fill(InTrk_ValidPixHits);
          
          if (isGlobal && isPF && GlbTrk_NormChi2 && GlbTrk_ValidMuonHits && MatchedStations && BestTrk_dXY && BestTrk_dZ && InTrk_ValidPixHits) 
            h_InTrk_TrkLayers[1].Fill(InTrk_TrkLayers);
          
          // Fill numerator histograms (TightMuonID cuts NOT applied except that cut)
          if (isGlobal && isPF && GlbTrk_NormChi2 && GlbTrk_ValidMuonHits && MatchedStations && BestTrk_dXY && BestTrk_dZ && InTrk_ValidPixHits && InTrk_TrkLayers) 
            h_isTight[3].Fill(isTight);
          
          if (isGlobal) h_isGlobal[3].Fill(isGlobal);
          
          if (isPF) h_isPF[3].Fill(isPF);
          
          if (GlbTrk_NormChi2) h_GlbTrk_NormChi2[3].Fill(GlbTrk_NormChi2);
          
          if (GlbTrk_ValidMuonHits) h_GlbTrk_ValidMuonHits[3].Fill(GlbTrk_ValidMuonHits);
          
          if (MatchedStations) h_MatchedStations[3].Fill(MatchedStations);
          
          if (BestTrk_dXY) h_BestTrk_dXY[3].Fill(BestTrk_dXY);
          
          if (BestTrk_dZ) h_BestTrk_dZ[3].Fill(BestTrk_dZ);
          
          if (InTrk_ValidPixHits) h_InTrk_ValidPixHits[3].Fill(InTrk_ValidPixHits);
          
          if (InTrk_TrkLayers) h_InTrk_TrkLayers[3].Fill(InTrk_TrkLayers);
        } // End of test isTight muons

        // Test isMedium cuts
        bool isMedium1        = v_isMedium.at(muIdx1);
        bool isMedium2        = v_isMedium.at(muIdx2);
        if ( (isMedium1 && !isMedium2) || (isMedium1 && isMedium2 && mu1->Pt() > mu2->Pt()) )
          muIdx = muIdx2; //muIdx1 is the best quality -> probe is muIdx2
        else if ( (!isMedium1 && isMedium2) || (isMedium1 && isMedium2 && mu1->Pt() < mu2->Pt()) )
          muIdx = muIdx2; //muIdx2 is the best quality -> probe is muIdx1
        else //all muons aren't good
          muIdx = 9999;
        
        if (muIdx!=9999) {
          bool isMedium       = v_isMedium.at(muIdx);
          bool isLoose        = v_isLoose.at(muIdx);
          bool isGlobal       = v_isGlobal.at(muIdx);
          float GlbTrk_NormChi2      = v_GlbTrk_NormChi2.at(muIdx);
          float Chi2Pos       = v_Chi2Pos.at(muIdx);
          float TrkKink       = v_TrkKink.at(muIdx);
          float SegmentComp   = v_SegmentComp.at(muIdx);
          float InTrk_ValFrac = v_InTrk_ValFrac.at(muIdx);

          // Fill denominator histograms (no cuts applied except isPF)
          h_isMedium[0].Fill(isMedium);
          h_isLoose[0].Fill(isLoose);
          h_isGlobal_Medium[0].Fill(isGlobal);
          h_GlbTrk_NormChi2_Medium[0].Fill(GlbTrk_NormChi2);
          h_Chi2Pos[0].Fill(Chi2Pos);
          h_TrkKink[0].Fill(TrkKink);
          h_SegmentComp[0].Fill(SegmentComp);
          num_SegmentComp[0]++;
          h_InTrk_ValFrac[0].Fill(InTrk_ValFrac);

          // Fill numerator histograms (MediumMuonID cuts applied except that cut)
          if (isMedium1 && isMedium2) {
            h_Mass[3].Fill(dimu->M());
            h_Pt[3].Fill(dimu->Pt());
          }

          bool pass_isMedium=false, pass_isLoose=false, pass_isGlobal=false;
          bool pass_GlbTrk_NormChi2=false, pass_Chi2Pos=false, pass_TrkKink=false;
          bool pass_InTrk_ValFrac=false;

          if (isMedium) pass_isMedium=true;
          if (isGlobal) pass_isGlobal=true;
          if (isLoose) pass_isLoose=true;
          if (GlbTrk_NormChi2<3) pass_GlbTrk_NormChi2=true;
          if (Chi2Pos<12) pass_Chi2Pos=true;
          if (TrkKink<20) pass_TrkKink=true;
          if (InTrk_ValFrac>0.8) pass_InTrk_ValFrac=true;
          
          //IsMedium conditions are defined as follow
          //bool goodGlob = (isGlobal && (GlbTrk_NormChi2<3) && (Chi2Pos<12) && (TrkKink<20); 
          //bool isMedium = (isLoose && (InTrk_ValFrac>0.8) && (SegmentComp > (goodGlob ? 0.303 : 0.451)));

          if (isLoose && pass_InTrk_ValFrac && ( SegmentComp > ((pass_isGlobal && pass_GlbTrk_NormChi2 && pass_Chi2Pos && pass_TrkKink) ? 0.303 : 0.451) ) )
            h_isMedium[1].Fill(isMedium);
          
          if (pass_InTrk_ValFrac && ( SegmentComp > ((pass_isGlobal && pass_GlbTrk_NormChi2 && pass_Chi2Pos && pass_TrkKink) ? 0.303 : 0.451) ) )
            h_isLoose[1].Fill(isLoose);

          if (isLoose && pass_InTrk_ValFrac && ( SegmentComp > ((pass_GlbTrk_NormChi2 && pass_Chi2Pos && pass_TrkKink) ? 0.303 : 0.451) ) )
            h_isGlobal_Medium[1].Fill(isGlobal);
          
          if (isLoose && pass_InTrk_ValFrac && ( SegmentComp > ((pass_isGlobal && pass_Chi2Pos && pass_TrkKink) ? 0.303 : 0.451) ) )
            h_GlbTrk_NormChi2_Medium[1].Fill(GlbTrk_NormChi2);
          
          if (isLoose && pass_InTrk_ValFrac && ( SegmentComp > ((pass_isGlobal && pass_GlbTrk_NormChi2 && pass_TrkKink) ? 0.303 : 0.451) ) )
            h_Chi2Pos[1].Fill(Chi2Pos);
          
          if (isLoose && pass_InTrk_ValFrac && ( SegmentComp > ((pass_isGlobal && pass_GlbTrk_NormChi2 && pass_Chi2Pos) ? 0.303 : 0.451) ) )
            h_TrkKink[1].Fill(TrkKink);
          
          if (isLoose && pass_InTrk_ValFrac) {
            h_SegmentComp[1].Fill(SegmentComp);
            num_SegmentComp[1]++;
          }

          if (isLoose && ( SegmentComp > ((pass_isGlobal && pass_GlbTrk_NormChi2 && pass_Chi2Pos && pass_TrkKink) ? 0.303 : 0.451) ) )
            h_InTrk_ValFrac[1].Fill(InTrk_ValFrac);

          // Fill numerator histograms (MediumMuonID cuts NOT applied except that cut)
          if (isMedium) h_isMedium[3].Fill(isMedium);
          
          if (isLoose) h_isLoose[3].Fill(isLoose);

          if (isGlobal) h_isGlobal_Medium[3].Fill(isGlobal);
          
          if (GlbTrk_NormChi2<3) h_GlbTrk_NormChi2_Medium[3].Fill(GlbTrk_NormChi2);
          
          if (Chi2Pos<12) h_Chi2Pos[3].Fill(Chi2Pos);
          
          if (TrkKink<20) h_TrkKink[3].Fill(TrkKink);
          
          if ( SegmentComp > ((pass_isGlobal && pass_GlbTrk_NormChi2 && pass_Chi2Pos && pass_TrkKink) ? 0.303 : 0.451) ) {
            h_SegmentComp[3].Fill(SegmentComp);
            num_SegmentComp[3]++;
          }

          if (InTrk_ValFrac>0.8) h_InTrk_ValFrac[3].Fill(InTrk_ValFrac);

        } // End of test isMedium

      } // end of isValidZmumu() condition
    } // end of dimuon loop
  } // end of evt loop

  getEfficiency();

  // Store all histograms
  fout->cd();
  hPt.Write();
  hMass.Write();
  hDxy.Write();
  hDz.Write();
  for (int i=0; i<5; i++) {
    h_Pt[i].Write();
    h_Mass[i].Write();
  }
  for (int i=0; i<5; i++) {
    h_isTight[i].Write();
    h_isGlobal[i].Write();
    h_isPF[i].Write();
    h_GlbTrk_NormChi2[i].Write();
    h_GlbTrk_ValidMuonHits[i].Write();
    h_MatchedStations[i].Write();
    h_BestTrk_dXY[i].Write();
    h_BestTrk_dZ[i].Write();
    h_InTrk_ValidPixHits[i].Write();
    h_InTrk_TrkLayers[i].Write();

    h_isMedium[i].Write();
    h_isLoose[i].Write();
    h_isGlobal_Medium[i].Write();
    h_GlbTrk_NormChi2_Medium[i].Write();
    h_Chi2Pos[i].Write();
    h_TrkKink[i].Write();
    h_SegmentComp[i].Write();
    h_InTrk_ValFrac[i].Write();
  }
  fout->Close();
  cout << "MuIDTest::variableTest() after saving histograms" << endl;

  cout << "\t\tnum_SegmentComp ";
  for (int i=0; i<5; i++) {
    cout << num_SegmentComp[i] << " ";
  }
  cout << endl;

}







//////////////////////////////////////////////////
////////////// S U B R O U T I N E ///////////////
//////////////////////////////////////////////////
void setHistoError(TH1 *h) {
  const int nbins = h->GetNbinsX();
  for (int i=0; i<=nbins; i++) {
    h->SetBinError(i,0);
  }
}


void MuIDVariableTest_Sub(const double muPtCut) {
  string inputFileMC="root://eoscms//eos/cms/store/group/phys_heavyions/anstahll/EWQAnalysis2017/pPb2016/8160GeV/MC/HiEWQForest_DYToMuMu_pPb_8160GeV_20170323.root";
  string outputFileMC=Form("outMC_muPtCut_%.0f.root",muPtCut);
  string inputFileData="root://eoscms//eos/cms/store/group/phys_heavyions/anstahll/EWQAnalysis2017/pPb2016/8160GeV/Data/HiEWQForest_PASingleMuon_pPb_Pbp_8160GeV_20170323.root";
  string outputFileData=Form("outData_muPtCut_%.0f.root",muPtCut);
  
  MuIDTest mctest(true, inputFileMC, outputFileMC, muPtCut);
  mctest.init();
  mctest.variableTest();

  MuIDTest datatest(false, inputFileData, outputFileData, muPtCut);
  datatest.init();
  datatest.variableTest();
}

void drawHisto_detail(TFile *fMC, TFile *fData, const double muPtCut, const string plot, const double xvalue1=-99, const double xvalue2=-99) {
  TCanvas canv("canv","canv",600,600);
  TPad pad1("pad1","pad1",0,0.3,1,1);
  pad1.SetTopMargin(0.04);
  pad1.SetBottomMargin(0);
  pad1.SetRightMargin(0.04);
  pad1.Draw();
  TPad pad2("pad2","pad2",0,0,1,0.3);
  pad2.SetTopMargin(0);
  pad2.SetBottomMargin(0.24);
  pad2.SetRightMargin(0.04);
  pad2.Draw();

  TLine line;
  line.SetLineWidth(2);
  line.SetLineStyle(7);
  line.SetLineColor(kGray+3);

  TLatex lat;
  lat.SetTextSize(0.04);
  lat.SetNDC(kTRUE);

  TLegend leg(0.67,0.8,0.95,0.95);
  leg.SetFillColor(0);
  leg.SetFillStyle(4000);
  leg.SetBorderSize(0);
  leg.SetMargin(0.15);


  // MC
  // All cuts except this cut(1), NO cuts except this cut(3)
  TH1D *hMC0 = (TH1D*)fMC->Get(Form("h_%s_0",plot.c_str()));
  TH1D *hMC1 = (TH1D*)fMC->Get(Form("h_%s_1",plot.c_str()));
  TH1D *hMC3 = (TH1D*)fMC->Get(Form("h_%s_3",plot.c_str()));

  TH1D *hRatio1 = (TH1D*)fMC->Get(Form("h_%s_2",plot.c_str()));
  TH1D *hRatio3 = (TH1D*)fMC->Get(Form("h_%s_4",plot.c_str()));
  setHistoError(hRatio1);
  setHistoError(hRatio3);

  double xaxismin, xaxismax; // for line drawing

  if (plot.find("dXY")!=std::string::npos ||
      plot.find("Dxy")!=std::string::npos) {
    hMC0->GetXaxis()->SetRangeUser(-0.15,0.15);
    hMC1->GetXaxis()->SetRangeUser(-0.15,0.15);
    hMC3->GetXaxis()->SetRangeUser(-0.15,0.15);
    hRatio1->GetXaxis()->SetRangeUser(-0.15,0.15);
    hRatio3->GetXaxis()->SetRangeUser(-0.15,0.15);
    pad1.SetLogy(1);
    xaxismin = -0.15;
    xaxismax = 0.15;
  } else if (plot.find("dZ")!=std::string::npos ||
             plot.find("Dz")!=std::string::npos) {
    hMC0->GetXaxis()->SetRangeUser(-0.2,0.2);
    hMC1->GetXaxis()->SetRangeUser(-0.2,0.2);
    hMC3->GetXaxis()->SetRangeUser(-0.2,0.2);
    hRatio1->GetXaxis()->SetRangeUser(-0.2,0.2);
    hRatio3->GetXaxis()->SetRangeUser(-0.2,0.2);
    pad1.SetLogy(1);
    xaxismin = -0.2;
    xaxismax = 0.2;
  } else if (plot.find("GlbTrk_NormChi2")!=std::string::npos) {
    pad1.SetLogy(1);
    xaxismin = hMC0->GetXaxis()->GetBinLowEdge(1);
    xaxismax = hMC0->GetXaxis()->GetBinLowEdge( hMC0->GetNbinsX() ) + hMC0->GetXaxis()->GetBinWidth( hMC0->GetNbinsX() );
  } else {
    pad1.SetLogy(0);
    xaxismin = hMC0->GetXaxis()->GetBinLowEdge(1);
    xaxismax = hMC0->GetXaxis()->GetBinLowEdge( hMC0->GetNbinsX() ) + hMC0->GetXaxis()->GetBinWidth( hMC0->GetNbinsX() );
  }

  hMC0->SetLineColor(kRed+1);
  hMC0->SetFillColor(kRed+1);
  hMC0->SetFillStyle(3445);
  hMC0->SetLineWidth(1);
  hMC1->SetLineColor(kGreen+1);
  hMC1->SetFillColor(kGreen);
  hMC1->SetLineWidth(2);
  hMC1->SetFillStyle(0);
  hMC3->SetMarkerStyle(kOpenCircle);
  hMC3->SetLineColor(kBlue+1);
  hMC3->SetMarkerColor(kBlue+1);
  hMC3->SetLineWidth(2);
  hMC3->SetFillStyle(0);

  hRatio1->SetMarkerStyle(kOpenCircle);
  hRatio1->SetMarkerColor(kGreen+1);
  hRatio1->SetLineColor(kGreen+1);
  hRatio1->SetFillStyle(0);
  hRatio3->SetMarkerStyle(kOpenCircle);
  hRatio3->SetMarkerColor(kBlue+1);
  hRatio3->SetLineColor(kBlue+1);
  hRatio3->SetFillStyle(0);

  leg.AddEntry(hMC0,"All Z events", "lf");
  if (plot.find("Mass")!=std::string::npos ||
      plot.find("Pt")!=std::string::npos) {
    leg.AddEntry(hMC1,"TightMuonID", "l");
    leg.AddEntry(hMC3,"MediumMuonID", "l");
  } else {
    leg.AddEntry(hMC1,"All cuts except this cut", "lf");
  }
//  leg.AddEntry(hMC3,"Only this cut ", "lf");

  pad1.cd();
  hMC0->Draw("HIST");
  hMC1->Draw("HIST same");
  if (plot.find("Mass")!=std::string::npos ||
      plot.find("Pt")!=std::string::npos) {
    hMC3->Draw("HIST same");
  }
  leg.Draw();
  lat.DrawLatex(0.65,0.75,"Z #rightarrow #mu#mu MC");
  double min = hMC0->GetMinimum();
  double max = hMC0->GetMaximum();
  if (xvalue1!=-99) line.DrawLine(xvalue1,min,xvalue1,max);
  if (xvalue2!=-99) line.DrawLine(xvalue2,min,xvalue2,max);

  pad2.cd();
  hRatio1->GetXaxis()->SetTitleSize(0.1);
  hRatio1->GetYaxis()->SetTitleSize(0.1);
  hRatio1->GetXaxis()->SetLabelSize(0.1);
  hRatio1->GetYaxis()->SetLabelSize(0.1);
  hRatio1->GetYaxis()->SetTitleOffset(0.4);
  hRatio1->GetYaxis()->SetTitle("Efficiency");
  hRatio1->GetYaxis()->SetRangeUser(0.,1.2);
  hRatio1->Draw("p");
  if (plot.find("Mass")!=std::string::npos ||
      plot.find("Pt")!=std::string::npos) {
    hRatio3->Draw("p same");
  }
  line.DrawLine(xaxismin,1,xaxismax,1);
  if (xvalue1!=-99) line.DrawLine(xvalue1,0,xvalue1,1.2);
  if (xvalue2!=-99) line.DrawLine(xvalue2,0,xvalue2,1.2);

  string dir = "figs";
  gSystem->mkdir(dir.c_str(),kTRUE);
  canv.SaveAs(Form("%s/%s_MC_muPtCut%.0f.png",dir.c_str(),plot.c_str(),muPtCut));
  canv.SaveAs(Form("%s/%s_MC_muPtCut%.0f.pdf",dir.c_str(),plot.c_str(),muPtCut));

  canv.Clear();
  pad1.Draw();
  pad2.Draw();

  // Data
  TH1D *hData0 = (TH1D*)fData->Get(Form("h_%s_0",plot.c_str()));
  TH1D *hData1 = (TH1D*)fData->Get(Form("h_%s_1",plot.c_str()));
  TH1D *hData3 = (TH1D*)fData->Get(Form("h_%s_3",plot.c_str()));

  hRatio1 = (TH1D*)fData->Get(Form("h_%s_2",plot.c_str()));
  hRatio3 = (TH1D*)fData->Get(Form("h_%s_4",plot.c_str()));
  setHistoError(hRatio1);
  setHistoError(hRatio3);
 
  if (plot.find("dXY")!=std::string::npos ||
      plot.find("Dxy")!=std::string::npos) {
    hData0->GetXaxis()->SetRangeUser(-0.15,0.15);
    hData1->GetXaxis()->SetRangeUser(-0.15,0.15);
    hData3->GetXaxis()->SetRangeUser(-0.15,0.15);
    hRatio1->GetXaxis()->SetRangeUser(-0.15,0.15);
    hRatio3->GetXaxis()->SetRangeUser(-0.15,0.15);
    pad1.SetLogy(1);
    xaxismin = -0.15;
    xaxismax = 0.15;
  } else if (plot.find("dZ")!=std::string::npos ||
             plot.find("Dz")!=std::string::npos) {
    hData0->GetXaxis()->SetRangeUser(-0.2,0.2);
    hData1->GetXaxis()->SetRangeUser(-0.2,0.2);
    hData3->GetXaxis()->SetRangeUser(-0.2,0.2);
    hRatio1->GetXaxis()->SetRangeUser(-0.2,0.2);
    hRatio3->GetXaxis()->SetRangeUser(-0.2,0.2);
    pad1.SetLogy(1);
    xaxismin = -0.2;
    xaxismax = 0.2;
  } else if (plot.find("GlbTrk_NormChi2")!=std::string::npos) {
    pad1.SetLogy(1);
    xaxismin = hData0->GetXaxis()->GetBinLowEdge(1);
    xaxismax = hData0->GetXaxis()->GetBinLowEdge( hData0->GetNbinsX() ) + hData0->GetXaxis()->GetBinWidth( hData0->GetNbinsX() );
  } else {
    pad1.SetLogy(0);
    xaxismin = hData0->GetXaxis()->GetBinLowEdge(1);
    xaxismax = hData0->GetXaxis()->GetBinLowEdge( hData0->GetNbinsX() ) + hData0->GetXaxis()->GetBinWidth( hData0->GetNbinsX() );
  }
  
  hData0->SetLineColor(kRed+1);
  hData0->SetFillColor(kRed+1);
  hData0->SetFillStyle(3445);
  hData0->SetLineWidth(1);
  hData1->SetLineColor(kGreen+1);
  hData1->SetFillColor(kGreen);
  hData1->SetLineWidth(2);
  hData1->SetFillStyle(0);
  hData3->SetMarkerStyle(kOpenCircle);
  hData3->SetLineColor(kBlue+1);
  hData3->SetMarkerColor(kBlue+1);
  hData3->SetFillStyle(0);

  hRatio1->SetMarkerStyle(kOpenCircle);
  hRatio1->SetMarkerColor(kGreen+1);
  hRatio1->SetLineColor(kGreen+1);
  hRatio1->SetFillStyle(0);
  hRatio3->SetMarkerStyle(kOpenCircle);
  hRatio3->SetMarkerColor(kBlue+1);
  hRatio3->SetLineColor(kBlue+1);
  hRatio3->SetFillStyle(0);

  pad1.cd();
  hData0->Draw("HIST");
  hData1->Draw("HIST same");
  if (plot.find("Mass")!=std::string::npos ||
      plot.find("Pt")!=std::string::npos) {
    hData3->Draw("HIST same");
  }
  leg.Draw();
  lat.DrawLatex(0.65,0.75,"Z #rightarrow #mu#mu Data");
  min = hData0->GetMinimum();
  max = hData0->GetMaximum();
  if (xvalue1!=-99) line.DrawLine(xvalue1,min,xvalue1,max);
  if (xvalue2!=-99) line.DrawLine(xvalue2,min,xvalue2,max);

  pad2.cd();
  hRatio1->GetXaxis()->SetTitleSize(0.1);
  hRatio1->GetYaxis()->SetTitleSize(0.1);
  hRatio1->GetXaxis()->SetLabelSize(0.1);
  hRatio1->GetYaxis()->SetLabelSize(0.1);
  hRatio1->GetYaxis()->SetTitleOffset(0.4);
  hRatio1->GetYaxis()->SetTitle("Efficiency");
  hRatio1->GetYaxis()->SetRangeUser(0.,1.2);
  hRatio1->Draw("p");
  if (plot.find("Mass")!=std::string::npos ||
      plot.find("Pt")!=std::string::npos) {
    hRatio3->Draw("p same");
  }
  line.DrawLine(xaxismin,1,xaxismax,1);
  if (xvalue1!=-99) line.DrawLine(xvalue1,0,xvalue1,1.2);
  if (xvalue2!=-99) line.DrawLine(xvalue2,0,xvalue2,1.2);

  canv.SaveAs(Form("%s/%s_DATA_muPtCut%.0f.png",dir.c_str(),plot.c_str(),muPtCut));
  canv.SaveAs(Form("%s/%s_DATA_muPtCut%.0f.pdf",dir.c_str(),plot.c_str(),muPtCut));
 
  //// Get integrated efficiency!
  int nbins = hMC1->GetNbinsX();
  int xbin1 = (xvalue1!=-99) ? hMC1->FindBin(xvalue1) : 1;
  int xbin2 = (xvalue2!=-99) ? hMC1->FindBin(xvalue2) : nbins;
  cout << "\t\t" << plot << " type2 MC ";
  cout << "\t\t" << nbins << "  " << xbin1 << "  " << xbin2 << endl;
  cout << "\t\t" << hMC0->Integral() << "  "
       << hMC1->Integral(1,xbin1) << "  " << hMC1->Integral(xbin1,nbins) << "  " << hMC1->Integral(xbin1,xbin2) << endl;
  cout << "\t\t\t" << hMC1->Integral(1,xbin1-1) << "  " << hMC1->Integral(xbin1+1,nbins) << endl;
  cout << "\t\t" << plot << " type2 DATA ";
  cout << "\t\t" << nbins << "  " << xbin1 << "  " << xbin2 << endl;
  cout << "\t\t" << hData0->Integral() << "  "
       << hData1->Integral(1,xbin1) << "  " << hData1->Integral(xbin1,nbins) << "  " << hData1->Integral(xbin1,xbin2) << endl;
  cout << "\t\t\t" << hData1->Integral(1,xbin1-1) << "  " << hData1->Integral(xbin1+1,nbins) << endl;

}

void drawHisto_simple(TFile *fMC, TFile *fData, const double muPtCut, const string plot) {
  TCanvas canv("canv","canv",600,600);
  TPad pad1("pad1","pad1",0,0.3,1,1);
  pad1.SetTopMargin(0.04);
  pad1.SetBottomMargin(0);
  pad1.SetRightMargin(0.04);
  pad1.Draw();
  TPad pad2("pad2","pad2",0,0,1,0.3);
  pad2.SetTopMargin(0);
  pad2.SetBottomMargin(0.24);
  pad2.SetRightMargin(0.04);
  pad2.Draw();

  TLine line;
  line.SetLineWidth(2);
  line.SetLineStyle(7);
  line.SetLineColor(kGray+3);

  TLegend leg(0.67,0.85,0.95,0.95);
  leg.SetFillColor(0);
  leg.SetFillStyle(4000);
  leg.SetBorderSize(0);
  leg.SetMargin(0.15);

  // All dimuons
  string hist;
  hist = "h" + plot;
  TH1D *hMC = (TH1D*)fMC->Get(hist.c_str());
  TH1D *hData = (TH1D*)fData->Get(hist.c_str());
  TH1D *hRatio = (TH1D*)hMC->Clone("hRatio_2");
  setHistoError(hRatio);

  if (plot.find("dXY")!=std::string::npos ||
      plot.find("Dxy")!=std::string::npos) {
    pad1.SetLogy(1);
  }
  else if (plot.find("dZ")!=std::string::npos ||
           plot.find("Dz")!=std::string::npos) {
    pad1.SetLogy(1);
  }
  else pad1.SetLogy(0);

  leg.AddEntry(hMC,"Z #rightarrow #mu#mu MC", "lf");
  leg.AddEntry(hData,"Z #rightarrow #mu#mu data", "lf");

  pad1.cd();
  hMC->DrawNormalized("HIST");
  hData->DrawNormalized("HIST same");
  leg.Draw();

  pad2.cd();
  hRatio->SetLineColor(kBlack);
  hRatio->SetMarkerSize(1);
  hRatio->SetMarkerStyle(kOpenCircle);
  hRatio->GetXaxis()->SetTitleSize(0.1);
  hRatio->GetYaxis()->SetTitleSize(0.1);
  hRatio->GetXaxis()->SetLabelSize(0.1);
  hRatio->GetYaxis()->SetLabelSize(0.1);
  hRatio->Divide(hMC,hData);
  hRatio->GetYaxis()->SetRangeUser(0.5,1.5);
  hRatio->Draw("p");
  double xaxismin = hMC->GetXaxis()->GetBinLowEdge(1);
  double xaxismax = hMC->GetXaxis()->GetBinLowEdge( hMC->GetNbinsX() ) + hMC->GetXaxis()->GetBinWidth( hMC->GetNbinsX() );
  line.DrawLine(xaxismin,1,xaxismax,1);

  string dir = "figs";
  gSystem->mkdir(dir.c_str(),kTRUE);
  canv.SaveAs(Form("%s/%s_simple_muPtCut%.0f.png",dir.c_str(),plot.c_str(),muPtCut));
  canv.SaveAs(Form("%s/%s_simple_muPtCut%.0f.pdf",dir.c_str(),plot.c_str(),muPtCut));
}

void drawHisto_eff(TFile *fMC, TFile *fData, const double muPtCut, const string plot, const double xvalue1=-99, const double xvalue2=-99) {
  TCanvas canv("canv","canv",600,600);
  TPad pad1("pad1","pad1",0,0.3,1,1);
  pad1.SetTopMargin(0.04);
  pad1.SetBottomMargin(0);
  pad1.SetRightMargin(0.04);
  pad1.Draw();
  TPad pad2("pad2","pad2",0,0,1,0.3);
  pad2.SetTopMargin(0);
  pad2.SetBottomMargin(0.24);
  pad2.SetRightMargin(0.04);
  pad2.Draw();

  TLine line;
  line.SetLineWidth(2);
  line.SetLineStyle(7);
  line.SetLineColor(kGray+3);

  TLegend leg(0.67,0.85,0.95,0.95);
  leg.SetFillColor(0);
  leg.SetFillStyle(4000);
  leg.SetBorderSize(0);
  leg.SetMargin(0.15);

  // All cuts except this cut
  string hist;
  hist = "h_" + plot + "_2";
  TH1D *hMC = (TH1D*)fMC->Get(hist.c_str());
  TH1D *hData = (TH1D*)fData->Get(hist.c_str());
  TH1D *hRatio = (TH1D*)hMC->Clone("hRatio_2");

  leg.AddEntry(hMC,"Z #rightarrow #mu#mu MC", "lf");
  leg.AddEntry(hData,"Z #rightarrow #mu#mu data", "lf");

  if (plot.find("Mass")!=std::string::npos ||
      plot.find("Pt")!=std::string::npos) {
    hMC->GetYaxis()->SetRangeUser(0.8,1.2);
  }
  else hMC->GetYaxis()->SetRangeUser(0,1.2);

  pad1.cd();
  hMC->GetYaxis()->SetTitle("Efficiency");
  hMC->Draw("HIST");
  hData->Draw("HIST same");
  if (xvalue1!=-99) line.DrawLine(xvalue1,0,xvalue1,1.2);
  if (xvalue2!=-99) line.DrawLine(xvalue2,0,xvalue2,1.2);
  leg.Draw();

  pad2.cd();
  hRatio->SetLineColor(kBlack);
  hRatio->SetMarkerSize(1.1);
  hRatio->SetMarkerStyle(kOpenCircle);
  hRatio->GetXaxis()->SetTitleSize(0.1);
  hRatio->GetYaxis()->SetTitleSize(0.1);
  hRatio->GetXaxis()->SetLabelSize(0.1);
  hRatio->GetYaxis()->SetLabelSize(0.1);
  hRatio->GetYaxis()->SetTitleOffset(0.4);
  hRatio->Divide(hMC,hData);
  hRatio->GetYaxis()->SetTitle("Eff MC/Eff Data");
  hRatio->GetYaxis()->SetRangeUser(0.5,1.5);
  setHistoError(hRatio);
  hRatio->Draw("p");
  double xaxismin = hMC->GetXaxis()->GetBinLowEdge(1);
  double xaxismax = hMC->GetXaxis()->GetBinLowEdge( hMC->GetNbinsX() ) + hMC->GetXaxis()->GetBinWidth( hMC->GetNbinsX() );
  line.DrawLine(xaxismin,1,xaxismax,1);
  if (xvalue1!=-99) line.DrawLine(xvalue1,0.5,xvalue1,1.5);
  if (xvalue2!=-99) line.DrawLine(xvalue2,0.5,xvalue2,1.5);

  string dir = "figs";
  gSystem->mkdir(dir.c_str(),kTRUE);
  canv.SaveAs(Form("%s/%s_type2_muPtCut%.0f.png",dir.c_str(),plot.c_str(),muPtCut));
  canv.SaveAs(Form("%s/%s_type2_muPtCut%.0f.pdf",dir.c_str(),plot.c_str(),muPtCut));

  canv.Clear();
  pad1.Draw();
  pad2.Draw();

  // NO cuts except this cut
  hist = "h_" + plot + "_4";
  hMC = (TH1D*)fMC->Get(hist.c_str());
  hData = (TH1D*)fData->Get(hist.c_str());

  if (plot.find("Mass")!=std::string::npos ||
      plot.find("Pt")!=std::string::npos) {
    hMC->GetYaxis()->SetRangeUser(0.8,1.2);
  }
  else hMC->GetYaxis()->SetRangeUser(0,1.2);

  pad1.cd();
  hMC->GetYaxis()->SetTitle("Efficiency");
  hMC->Draw("HIST");
  hData->Draw("HIST same");
  if (xvalue1!=-99) line.DrawLine(xvalue1,0,xvalue1,1.2);
  if (xvalue2!=-99) line.DrawLine(xvalue2,0,xvalue2,1.2);
  leg.Draw();

  pad2.cd();
  hRatio->SetLineColor(kBlack);
  hRatio->SetMarkerSize(1.1);
  hRatio->SetMarkerStyle(kOpenCircle);
  hRatio->GetXaxis()->SetTitleSize(0.1);
  hRatio->GetYaxis()->SetTitleSize(0.1);
  hRatio->GetXaxis()->SetLabelSize(0.1);
  hRatio->GetYaxis()->SetLabelSize(0.1);
  hRatio->GetYaxis()->SetTitleOffset(0.4);
  hRatio->Divide(hMC,hData);
  hRatio->GetYaxis()->SetTitle("Eff MC/Eff Data");
  hRatio->GetYaxis()->SetRangeUser(0.5,1.5);
  setHistoError(hRatio);
  hRatio->Draw("p");
  xaxismin = hMC->GetXaxis()->GetBinLowEdge(1);
  xaxismax = hMC->GetXaxis()->GetBinLowEdge( hMC->GetNbinsX() ) + hMC->GetXaxis()->GetBinWidth( hMC->GetNbinsX() );
  line.DrawLine(xaxismin,1,xaxismax,1);
  if (xvalue1!=-99) line.DrawLine(xvalue1,0.5,xvalue1,1.5);
  if (xvalue2!=-99) line.DrawLine(xvalue2,0.5,xvalue2,1.5);

  canv.SaveAs(Form("%s/%s_type4_muPtCut%.0f.png",dir.c_str(),plot.c_str(),muPtCut));
  canv.SaveAs(Form("%s/%s_type4_muPtCut%.0f.pdf",dir.c_str(),plot.c_str(),muPtCut));

}



//////////////////////////////////////////////////
//////////////////// M A I N /////////////////////
//////////////////////////////////////////////////
void drawHisto(const double muPtCut) {
  gStyle->SetOptStat(0);

  string FileMC=Form("outMC_muPtCut_%.0f.root",muPtCut);
  string FileData=Form("outData_muPtCut_%.0f.root",muPtCut);
  
  TFile *fMC = new TFile(FileMC.c_str());
  TFile *fData = new TFile(FileData.c_str());
  if (!fMC->IsOpen() || !fData->IsOpen()) return;

  string plot;

  // all muons
  drawHisto_simple(fMC, fData, muPtCut, "Mass");
  drawHisto_simple(fMC, fData, muPtCut, "Pt");
  drawHisto_simple(fMC, fData, muPtCut, "Dxy");
  drawHisto_simple(fMC, fData, muPtCut, "Dz");



  // One tight muon
  drawHisto_eff(fMC, fData, muPtCut, "Mass");
  drawHisto_eff(fMC, fData, muPtCut, "Pt");
  drawHisto_eff(fMC, fData, muPtCut, "isTight", 1);
  drawHisto_eff(fMC, fData, muPtCut, "isGlobal", 1);
  drawHisto_eff(fMC, fData, muPtCut, "isPF", 1);
  drawHisto_eff(fMC, fData, muPtCut, "GlbTrk_NormChi2", 10);
  drawHisto_eff(fMC, fData, muPtCut, "GlbTrk_ValidMuonHits", 1); // cut is 0 ! 1 is for line drawing only
  drawHisto_eff(fMC, fData, muPtCut, "MatchedStations", 2); // cut is 1 ! 2 is for line drawing only
  drawHisto_eff(fMC, fData, muPtCut, "BestTrk_dXY", 0.2, -0.2);
  drawHisto_eff(fMC, fData, muPtCut, "BestTrk_dZ", 0.5, -0.5);
  drawHisto_eff(fMC, fData, muPtCut, "InTrk_ValidPixHits", 1); // cut is 0 ! 1 is for line drawing only
  drawHisto_eff(fMC, fData, muPtCut, "InTrk_TrkLayers", 6); // cut is 5 ! 6 is for line drawing only

  //bool goodGlob = (isGlobal && (GlbTrk_NormChi2<3) && (Chi2Pos<12) && (TrkKink<20); 
  //bool isMedium = (isLoose && (InTrk_ValFrac>0.8) && (SegmentComp > (goodGlob ? 0.303 : 0.451)));
  // One medium muon
  drawHisto_eff(fMC, fData, muPtCut, "isMedium", 1);
  drawHisto_eff(fMC, fData, muPtCut, "isLoose", 1);
  drawHisto_eff(fMC, fData, muPtCut, "isGlobal_Medium", 1);
  drawHisto_eff(fMC, fData, muPtCut, "GlbTrk_NormChi2_Medium", 3);
  drawHisto_eff(fMC, fData, muPtCut, "Chi2Pos", 12);
  drawHisto_eff(fMC, fData, muPtCut, "TrkKink", 20);
  drawHisto_eff(fMC, fData, muPtCut, "SegmentComp", 0.303, 0.451);
  drawHisto_eff(fMC, fData, muPtCut, "InTrk_ValFrac", 0.8);

 

  // One tight muon
  drawHisto_detail(fMC, fData, muPtCut, "Mass");
  drawHisto_detail(fMC, fData, muPtCut, "Pt");
  drawHisto_detail(fMC, fData, muPtCut, "isTight", 1);
  drawHisto_detail(fMC, fData, muPtCut, "isGlobal", 1);
  drawHisto_detail(fMC, fData, muPtCut, "isPF", 1);
  drawHisto_detail(fMC, fData, muPtCut, "GlbTrk_NormChi2", 10);
  drawHisto_detail(fMC, fData, muPtCut, "GlbTrk_ValidMuonHits", 1); // cut is 0 ! 1 is for line drawing only
  drawHisto_detail(fMC, fData, muPtCut, "MatchedStations", 2); // cut is 1 ! 2 is for line drawing only
  drawHisto_detail(fMC, fData, muPtCut, "BestTrk_dXY", 0.2, -0.2);
  drawHisto_detail(fMC, fData, muPtCut, "BestTrk_dZ", 0.5, -0.5);
  drawHisto_detail(fMC, fData, muPtCut, "InTrk_ValidPixHits", 1); // cut is 0 ! 1 is for line drawing only
  drawHisto_detail(fMC, fData, muPtCut, "InTrk_TrkLayers", 6); // cut is 5 ! 6 is for line drawing only

  // One medium muon
  drawHisto_detail(fMC, fData, muPtCut, "isMedium", 1);
  drawHisto_detail(fMC, fData, muPtCut, "isLoose", 1);
  drawHisto_detail(fMC, fData, muPtCut, "isGlobal_Medium", 1);
  drawHisto_detail(fMC, fData, muPtCut, "GlbTrk_NormChi2_Medium", 3);
  drawHisto_detail(fMC, fData, muPtCut, "Chi2Pos", 12);
  drawHisto_detail(fMC, fData, muPtCut, "TrkKink", 20);
  drawHisto_detail(fMC, fData, muPtCut, "SegmentComp", 0.303, 0.451);
  drawHisto_detail(fMC, fData, muPtCut, "InTrk_ValFrac", 0.8);

}


void MuIDVariableTest() {
  double muPtCut = 15;
  MuIDVariableTest_Sub(muPtCut);
  drawHisto(muPtCut);
  
//  muPtCut = 20;
//  MuIDVariableTest_Sub(muPtCut);
//  drawHisto(muPtCut);
//  
//  muPtCut = 25;
//  MuIDVariableTest_Sub(muPtCut);
//  drawHisto(muPtCut);

}




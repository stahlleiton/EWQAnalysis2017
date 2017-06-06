#include "GetAccEff.h"
#include "HiMuonTree.h"
#include "HiMETTree.h"



//////////////////////////////////////
////////// M A I N ///////////////////
//////////////////////////////////////
void GetAccEff(
    const string beamDir="pPb",
    const int CombineAccEff=0
    ) {

  // Input/Output files
  string inputFile;
  if (beamDir.compare("pPb")==0) inputFile="root://eoscms//eos/cms/store/group/phys_heavyions/anstahll/EWQAnalysis2017/pPb2016/8160GeV/MC/Embedded/HiEWQForest_Embedded_WToMuNu_pPb_8160GeV_20170518.root";
  else if (beamDir.compare("Pbp")==0) inputFile="root://eoscms//eos/cms/store/group/phys_heavyions/anstahll/EWQAnalysis2017/pPb2016/8160GeV/MC/Embedded/HiEWQForest_Embedded_WToMuNu_Pbp_8160GeV_20170518.root";
  string outputFile=Form("AccEff_%s.root",beamDir.c_str());

  // Define efficiency histograms
  //[0]: mu+
  //[1]: mu-
  // CombineAccEff == 0: acc*eff
  // CombineAccEff == 1: eff, 2: 5 TeV W analysis setting eff
  TH1D hden_eta[2], hnum_eta[2];
  TH1D hden_pt[2][nbins_eta], hnum_pt[2][nbins_eta];
  TGraphAsymmErrors heff_eta[2];
  TGraphAsymmErrors heff_pt[2][nbins_eta];

  // CombineAccEff == 1: acc, 2: 5 TeV W analysis setting acc
  TH1D hden_acc_eta[2], hnum_acc_eta[2];
  TH1D hden_acc_pt[2][nbins_eta], hnum_acc_pt[2][nbins_eta];
  TGraphAsymmErrors heff_acc_eta[2];
  TGraphAsymmErrors heff_acc_pt[2][nbins_eta];

  string suffix = beamDir;
  for (int i=0; i<2; i++) {
    if (i==0) suffix = beamDir + "_muPlus";
    if (i==1) suffix = beamDir + "_muMinus";
    hden_eta[i] = TH1D(Form("hden_eta_%s",suffix.c_str()),";#eta_{lab};",nbins_eta,bins_eta);
    hnum_eta[i] = TH1D(Form("hnum_eta_%s",suffix.c_str()),";#eta_{lab};",nbins_eta,bins_eta);
    heff_eta[i] = TGraphAsymmErrors();
    heff_eta[i].SetName(Form("heff_eta_%s",suffix.c_str()));
    heff_eta[i].GetXaxis()->SetTitle("#eta_{lab}");
    heff_eta[i].GetYaxis()->SetTitle("Efficiency");
    if (CombineAccEff>0) {
      hden_acc_eta[i] = TH1D(Form("hden_acc_eta_%s",suffix.c_str()),";#eta_{lab};",nbins_eta,bins_eta);
      hnum_acc_eta[i] = TH1D(Form("hnum_acc_eta_%s",suffix.c_str()),";#eta_{lab};",nbins_eta,bins_eta);
      heff_acc_eta[i] = TGraphAsymmErrors();
      heff_acc_eta[i].SetName(Form("heff_acc_eta_%s",suffix.c_str()));
      heff_acc_eta[i].GetXaxis()->SetTitle("#eta_{lab}");
      heff_acc_eta[i].GetYaxis()->SetTitle("Efficiency");
    }
    for (int j=0; j<nbins_eta; j++) {
      hden_pt[i][j] = TH1D(Form("hden_pt_%.0f-%.0f_%s",bins_eta[j]*10,bins_eta[j+1]*10,suffix.c_str()),";p_{T};",nbins_pt,bins_pt);
      hnum_pt[i][j] = TH1D(Form("hnum_pt_%.0f-%.0f_%s",bins_eta[j]*10,bins_eta[j+1]*10,suffix.c_str()),";p_{T};",nbins_pt,bins_pt);
      heff_pt[i][j] = TGraphAsymmErrors();
      heff_pt[i][j].SetName(Form("heff_pt_%.0f-%.0f_%s",bins_eta[j]*10,bins_eta[j+1]*10,suffix.c_str()));
      heff_pt[i][j].GetXaxis()->SetTitle("p_{T}");
      heff_pt[i][j].GetYaxis()->SetTitle("Efficiency");
      if (CombineAccEff>0) {
        hden_acc_pt[i][j] = TH1D(Form("hden_acc_pt_%.0f-%.0f_%s",bins_eta[j]*10,bins_eta[j+1]*10,suffix.c_str()),";p_{T};",nbins_pt,bins_pt);
        hnum_acc_pt[i][j] = TH1D(Form("hnum_acc_pt_%.0f-%.0f_%s",bins_eta[j]*10,bins_eta[j+1]*10,suffix.c_str()),";p_{T};",nbins_pt,bins_pt);
        heff_acc_pt[i][j] = TGraphAsymmErrors();
        heff_acc_pt[i][j].SetName(Form("heff_acc_pt_%.0f-%.0f_%s",bins_eta[j]*10,bins_eta[j+1]*10,suffix.c_str()));
        heff_acc_pt[i][j].GetXaxis()->SetTitle("p_{T}");
        heff_acc_pt[i][j].GetYaxis()->SetTitle("Efficiency");
      }
    }
    
  }

  // Read TTree, check default selection conditions
  HiMuonTree mutree = HiMuonTree();
  mutree.GetTree(inputFile.c_str());
  cout << "GetAccEff:: after HiMuonTree GetTree " << endl;
  Long64_t nentries = mutree.GetEntries();

  HiMETTree mettree = HiMETTree();
  mettree.GetTree(inputFile.c_str());
  cout << "GetAccEff:: after HiMETTree GetTree " << endl;

  // Loop over events
  vector<Long64_t> skip_event;
  for (Long64_t evt=0; evt<nentries; evt++) {
    if (mutree.GetEntry(evt)<0) break;

    /////// ###### GEN objects will be proceed!
    // Load gen particles, check if it is W
    vector< pair<int, unsigned short> > v_GenMother; //gen mu idx, its mother's idx
    vector<int> v_Gen_Particle_PdgId = mutree.Gen_Particle_PdgId();
    vector<unsigned char> v_Gen_Particle_Status = mutree.Gen_Particle_Status();
    vector < vector < UShort_t > > vv_Gen_Particle_Mother_Idx = mutree.Gen_Particle_Mother_Idx();
    vector<char> v_Gen_Muon_Reco_Idx = mutree.Gen_Muon_Reco_Idx();
    vector<unsigned short> v_Gen_Muon_Particle_Idx = mutree.Gen_Muon_Particle_Idx();
    vector<TLorentzVector> Gen_Muon_Mom = mutree.Gen_Muon_Mom();
    vector<char> Gen_Muon_Charge = mutree.Gen_Muon_Charge();

    cout << endl;
    for (auto v_mu=v_Gen_Muon_Particle_Idx.begin(); v_mu!=v_Gen_Muon_Particle_Idx.end(); ++v_mu) {
      int ipar = *v_mu;
      int pdgid = v_Gen_Particle_PdgId.at(ipar);
      unsigned char status = v_Gen_Particle_Status.at(ipar);
      vector <UShort_t> v_Gen_Particle_Mother_Idx = vv_Gen_Particle_Mother_Idx.at(ipar);
      if (pdgid!=24) {
      cout << "\tin the main, the pre-mother vector size "
           << " " << evt << " " << ipar << ": " << pdgid << " " << static_cast<unsigned short>(status) << " "
           << vv_Gen_Particle_Mother_Idx.at(ipar).size() << "==" << v_Gen_Particle_Mother_Idx.size() << "/" << vv_Gen_Particle_Mother_Idx.size() << endl;
      }
      // Is this W? Trackdown to the mother 
      // gen muon is anti-particle, it has opposite sign to gen W, but same charge in reco muon
      unsigned short mother = trackDownMothers(v_Gen_Particle_PdgId, v_Gen_Particle_Status, vv_Gen_Particle_Mother_Idx, ipar, static_cast<int>(-24*(pdgid/13)), 3);

      if (mother==9999) { // no mother information available.
        continue;
      } else if (mother==8888) { // mother's idx is idx of itself, skip this muon and move to next muon
        continue;
      }
      
      // this gen muon is going to be used in reco too
      int igenmu = distance(v_Gen_Muon_Particle_Idx.begin(),v_mu);
      pair<int, unsigned short> pair = make_pair(igenmu, mother);
      v_GenMother.push_back(pair);
    
      cout << "W found " << evt << " gen mu " << igenmu
           << ": mother's idx " << mother << " and muon's particle idx " << *v_mu << "\n\t"
           << v_Gen_Particle_PdgId.at(*v_mu) << "/" << static_cast<unsigned short>(v_Gen_Particle_Status.at(*v_mu)) << " "
           << v_Gen_Particle_PdgId.at(mother) << "/" << static_cast<unsigned short>(v_Gen_Particle_Status.at(mother)) << " "
           << endl;

    } // end of W mother muon & genmu matching

    // who's the highest pT gen muon?
    cout << "v_GenMother: " << v_GenMother.size() << endl;

    int igenmu_keep = -1;
    double highest_pt = -1;
    for (auto v_genidx=v_GenMother.begin(); v_genidx!=v_GenMother.end(); ++v_genidx) {
      pair<int, unsigned short> pair = *v_genidx;
      int igenmu = pair.first;
      TLorentzVector genmu = Gen_Muon_Mom.at(igenmu);
      cout << "igenmu: " << igenmu << ", " << genmu.Eta() << ", " << genmu.Pt() << endl;
      if (genmu.Pt()>highest_pt) {
        highest_pt = genmu.Pt();
        igenmu_keep = distance(v_GenMother.begin(),v_genidx);
      }
    }

    // Fill denominator histograms
    if (igenmu_keep!=-1) {
      pair<int, unsigned short> igenmu = v_GenMother.at(igenmu_keep);
      TLorentzVector mu = Gen_Muon_Mom.at(igenmu.first);
      char charge = Gen_Muon_Charge.at(igenmu.first);

      if (charge>0) {
        if (CombineAccEff==0) {
          hden_eta[0].Fill(mu.Eta());
        } else { 
          if ( (CombineAccEff==1 && TMath::Abs(mu.Eta())<bins_eta[nbins_eta]) ||
               (CombineAccEff==2 && TMath::Abs(mu.Eta())<bins_eta[nbins_eta] && mu.Pt()>25) 
             ) {
            hden_eta[0].Fill(mu.Eta());
            hnum_acc_eta[0].Fill(mu.Eta());
          }
          hden_acc_eta[0].Fill(mu.Eta());
        }
        for (int j=0; j<nbins_eta; j++) {
          if (bins_eta[j]<=mu.Eta() && bins_eta[j+1]>mu.Eta()) {
            if (CombineAccEff==0) {
              hden_pt[0][j].Fill(mu.Pt());
            } else {
              if ( (CombineAccEff==1 && TMath::Abs(mu.Eta())<bins_eta[nbins_eta]) ||
                   (CombineAccEff==2 && TMath::Abs(mu.Eta())<bins_eta[nbins_eta] && mu.Pt()>25) 
                 ) {
                hden_pt[0][j].Fill(mu.Pt());
                hnum_acc_pt[0][j].Fill(mu.Pt());
              }
              hden_acc_pt[0][j].Fill(mu.Pt());
            }
          }
        }
      } else {
        if (CombineAccEff==0) {
          hden_eta[1].Fill(mu.Eta());
        } else {
          if ( (CombineAccEff==1 && TMath::Abs(mu.Eta())<bins_eta[nbins_eta]) ||
               (CombineAccEff==2 && TMath::Abs(mu.Eta())<bins_eta[nbins_eta] && mu.Pt()>25) 
             ) {
            hden_eta[1].Fill(mu.Eta());
            hnum_acc_eta[1].Fill(mu.Eta());
          }
          hden_acc_eta[1].Fill(mu.Eta());
        }
        for (int j=0; j<nbins_eta; j++) {
          if (bins_eta[j]<=mu.Eta() && bins_eta[j+1]>mu.Eta()) {
            if (CombineAccEff==0) {
              hden_pt[1][j].Fill(mu.Pt());
            } else {
              if ( (CombineAccEff==1 && TMath::Abs(mu.Eta())<bins_eta[nbins_eta]) ||
                   (CombineAccEff==2 && TMath::Abs(mu.Eta())<bins_eta[nbins_eta] && mu.Pt()>25) 
                 ) {
                hden_pt[1][j].Fill(mu.Pt());
                hnum_acc_pt[1][j].Fill(mu.Pt());
              }
              hden_acc_pt[1][j].Fill(mu.Pt());
            }
          }
        }
      } // fill denominator depending on its charge
    } // end of filling up denominator with the highest pT muon



    /////// ###### RECO objects will be proceed!
    
    // Is this event fine at gen level? if not, skip the event
    auto skip_evt_idx = find(skip_event.begin(),skip_event.end(),evt);
    if (skip_evt_idx != skip_event.end()) continue;
    
    if (mettree.GetEntry(evt)<0) {
      cout << "Cannot retrieve MET tree, break" << endl;
      break;
    }
    
    // Select collision events without pile-up rejection
    bool Flag_collisionEventSelectionPA = mettree.Flag_collisionEventSelectionPA();
    // MET filters for 2017 Moriond
    bool Flag_goodVertices = mettree.Flag_goodVertices();
    bool Flag_CSCTightHaloFilter = mettree.Flag_CSCTightHaloFilter();
    bool Flag_HBHENoiseFilter = mettree.Flag_HBHENoiseFilter();
    bool Flag_HBHENoiseIsoFilter = mettree.Flag_HBHENoiseIsoFilter();
    bool Flag_EcalDeadCellTriggerPrimitiveFilter = mettree.Flag_EcalDeadCellTriggerPrimitiveFilter();
    bool Flag_chargedHadronTrackResolutionFilter = mettree.Flag_chargedHadronTrackResolutionFilter();
    bool Flag_muonBadTrackFilter = mettree.Flag_muonBadTrackFilter();
    bool Flag_eeBadScFilter = mettree.Flag_eeBadScFilter();
    // Fake muon rejections, can be changed
    bool Flag_duplicateMuons = mettree.Flag_duplicateMuons();
    bool Flag_badMuons = mettree.Flag_badMuons();
    
    if (Flag_collisionEventSelectionPA != 1 ||
        Flag_goodVertices != 1 ||
        Flag_CSCTightHaloFilter != 1 ||
        Flag_HBHENoiseFilter != 1 ||
        Flag_HBHENoiseIsoFilter != 1 ||
        Flag_EcalDeadCellTriggerPrimitiveFilter != 1 || // MC doesn't work
        Flag_chargedHadronTrackResolutionFilter != 1 ||
        Flag_muonBadTrackFilter != 1 ||
//        Flag_eeBadScFilter != 1 || // only for data
        Flag_duplicateMuons != 1 ||
        Flag_badMuons != 1 ) continue;

    // RECO muons to be checked
    vector<TLorentzVector> Reco_Muon_Mom = mutree.Reco_Muon_Mom();
    vector<char> v_Reco_Muon_Gen_Idx = mutree.Reco_Muon_Gen_Idx();
    vector<char> v_Reco_Muon_Charge = mutree.Reco_Muon_Charge();
    vector<bool> v_Reco_Muon_isTight = mutree.Reco_Muon_isTight();
    vector<bool> v_Reco_Muon_isMedium = mutree.Reco_Muon_isMedium();
    vector<char>  v_Gen_Muon_PF_Idx = mutree.Gen_Muon_PF_Idx();
    vector<float> v_PF_Muon_IsoPFR03NoPUCorr = mutree.PF_Muon_IsoPFR03NoPUCorr();

    // gen-reco matching with v_GenMother
    if (igenmu_keep!=-1) {
      pair<int, unsigned short> p_igenmu = v_GenMother.at(igenmu_keep);
      int igenmu = p_igenmu.first;

      // use reco muons only when there's matched gen muon
      cout << "RECO" << endl;

      char irecmu = v_Gen_Muon_Reco_Idx.at(igenmu);
      cout << "irecmu: " << static_cast<short>(irecmu) << endl;
      if (irecmu<0) {
        cout << "W eff: skip recmu cause a reco muon has to be found" << endl;
        continue; // no matched gen-reco muon pair, skip reco step
      }

      // Check if this muon passes quality cuts
      if (!v_Reco_Muon_isTight.at(irecmu)) continue;
      
      // Check if trigger is fired (Currently, trigger isn't working in sample)
      vector<bool> v_Event_Trig_Fired = mutree.Event_Trig_Fired();
      cout << "after v_Event_Trig_Fired " << endl;
      vector< vector< UChar_t > > vv_Pat_Muon_Trig = mutree.Pat_Muon_Trig();
      cout << "after vv_Pat_Muon_Trig " << endl;
      vector< UChar_t > v_Pat_Muon_Trig = vv_Pat_Muon_Trig.at(irecmu);
      cout << "after v_Pat_Muon_Trig " << endl;

      //5th trigger is L3Mu12
      if (!v_Event_Trig_Fired.at(5) || !v_Pat_Muon_Trig.at(5)) continue;
     
      char ipfmu = v_Gen_Muon_PF_Idx.at(igenmu);
      
      TLorentzVector recmu = Reco_Muon_Mom.at(irecmu);
      TLorentzVector genmu = Gen_Muon_Mom.at(igenmu);
      cout << "At numerator, " << static_cast<short>(ipfmu) << " , genmu: " << igenmu << " " << genmu.Eta() << " " << genmu.Pt() << endl;
      // Check if this muon pass kinematic selections
      if (recmu.Pt()>=25 && TMath::Abs(recmu.Eta())<bins_eta[nbins_eta] && v_PF_Muon_IsoPFR03NoPUCorr.at(ipfmu)<0.15) {
        char charge = v_Reco_Muon_Charge.at(irecmu);
        cout << "recmu charge " << static_cast<short>(charge) << endl;
        if (charge>0) {
          hnum_eta[0].Fill(genmu.Eta());
          for (int j=0; j<nbins_eta; j++) {
            if (bins_eta[j]<=recmu.Eta() && bins_eta[j+1]>recmu.Eta()) {
              hnum_pt[0][j].Fill(genmu.Pt());
            }
          }
        } else {
          hnum_eta[1].Fill(genmu.Eta());
          for (int j=0; j<nbins_eta; j++) {
            if (bins_eta[j]<=recmu.Eta() && bins_eta[j+1]>recmu.Eta()) {
              hnum_pt[1][j].Fill(genmu.Pt());
            }
          }
        }
      } // end of numerator histograms filled up
    } // end of gen-reco matching with v_GenMother
  } // end of evt loop


  // Calculate efficiency
  for (int i=0; i<2; i++) {
    CheckUnderFlow(hnum_eta[i],hden_eta[i]);
    heff_eta[i].Divide(&hnum_eta[i],&hden_eta[i]);
    if (CombineAccEff>0) {
      CheckUnderFlow(hnum_acc_eta[i],hden_acc_eta[i]);
      heff_acc_eta[i].Divide(&hnum_acc_eta[i],&hden_acc_eta[i]);
    }
    for (int j=0; j<nbins_eta; j++) {
      CheckUnderFlow(hnum_pt[i][j],hden_pt[i][j]);
      heff_pt[i][j].Divide(&hnum_pt[i][j],&hden_pt[i][j]);
      if (CombineAccEff>0) {
        CheckUnderFlow(hnum_acc_pt[i][j],hden_acc_pt[i][j]);
        heff_acc_pt[i][j].Divide(&hnum_acc_pt[i][j],&hden_acc_pt[i][j]);
      }
    }
  }


  // Write histograms into output file
  TFile foutput(outputFile.c_str(),"recreate");
  foutput.cd();
  for (int i=0; i<2; i++) {
    hden_eta[i].Write();
    hnum_eta[i].Write();
    heff_eta[i].Write();
    if (CombineAccEff>0) {
      hden_acc_eta[i].Write();
      hnum_acc_eta[i].Write();
      heff_acc_eta[i].Write();
    }
    for (int j=0; j<nbins_eta; j++) {
      hden_pt[i][j].Write();
      hnum_pt[i][j].Write();
      heff_pt[i][j].Write();
      if (CombineAccEff>0) {
        hden_acc_pt[i][j].Write();
        hnum_acc_pt[i][j].Write();
        heff_acc_pt[i][j].Write();
      }
    }
  }
  foutput.Close();


}

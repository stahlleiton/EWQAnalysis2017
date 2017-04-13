#include "GetAccEff.h"
#include "HiMuonTree.h"
#include "HiMETTree.h"



//////////////////////////////////////
////////// M A I N ///////////////////
//////////////////////////////////////
void GetAccEff_DYZ(
    const string muIDType="Tight",
    const string beamDir="pPb"
    ) {

  // Input/Output files
  string inputFile;
  if (beamDir.compare("pPb")==0) inputFile="root://eoscms//eos/cms/store/group/phys_heavyions/anstahll/EWQAnalysis2017/pPb2016/8160GeV/MC/HiEWQForest_DYToMuMu_pPb_8160GeV_20170323.root";
  else if (beamDir.compare("Pbp")==0) inputFile="root://eoscms//eos/cms/store/group/phys_heavyions/anstahll/EWQAnalysis2017/pPb2016/8160GeV/MC/HiEWQForest_DYToMuMu_Pbp_8160GeV_20170323.root";
  string outputFile=Form("AccEff_%s_%s.root",muIDType.c_str(),beamDir.c_str());

  // Define efficiency histograms
  // CombineAccEff == true:  acc*eff
  // num[0]: normal
  // num[1]: veto eff
  TH1D hden_eta, hnum_eta[2];
  TGraphAsymmErrors heff_eta[2];

  string suffix_ = beamDir + "_" + muIDType;
  string suffix = suffix_;
  hden_eta = TH1D(Form("hden_eta_%s",suffix.c_str()),";#eta_{lab};",nbins_eta,bins_eta);
  for (int i=0; i<2; i++) {
    if (i==0) suffix = suffix_ + "_normal";
    if (i==1) suffix = suffix_ + "_veto";
    hnum_eta[i] = TH1D(Form("hnum_eta_%s",suffix.c_str()),";#eta_{lab};",nbins_eta,bins_eta);
    heff_eta[i] = TGraphAsymmErrors();
    heff_eta[i].SetName(Form("heff_eta_%s",suffix.c_str()));
    heff_eta[i].GetXaxis()->SetTitle("#eta_{lab}");
    heff_eta[i].GetYaxis()->SetTitle("Efficiency");
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
  for (Long64_t evt=763; evt<nentries; evt++) {
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

    cout << endl;
    for (auto v_mu=v_Gen_Muon_Particle_Idx.begin(); v_mu!=v_Gen_Muon_Particle_Idx.end(); ++v_mu) {
      int ipar = *v_mu;
      int pdgid = v_Gen_Particle_PdgId.at(ipar);
      unsigned char status = v_Gen_Particle_Status.at(ipar);
      vector <UShort_t> v_Gen_Particle_Mother_Idx = vv_Gen_Particle_Mother_Idx.at(ipar);
      if (pdgid!=22) {
      cout << "\tin the main, the pre-mother vector size "
           << " " << evt << " " << ipar << ": " << pdgid << " " << static_cast<unsigned short>(status) << " "
           << vv_Gen_Particle_Mother_Idx.at(ipar).size() << "==" << v_Gen_Particle_Mother_Idx.size() << "/" << vv_Gen_Particle_Mother_Idx.size() << endl;
      }
      // Is this DY or Z? Trackdown to the mother
      unsigned short mother = trackDownMothers(v_Gen_Particle_PdgId, v_Gen_Particle_Status, vv_Gen_Particle_Mother_Idx, ipar, 23, 3);

      if (mother==9999) { // no mother information available. check photon instead
        mother = trackDownMothers(v_Gen_Particle_PdgId, v_Gen_Particle_Status, vv_Gen_Particle_Mother_Idx, ipar, 22, 3);
        if (mother==9999) { // no mother information available for photon too. skip this muon and move to next muon
          continue;
        }
      } else if (mother==8888) { // mother's idx is idx of itself, skip this muon and move to next muon
        continue;
      }
      
      // this gen muon is going to be used in reco too
      int igenmu = distance(v_Gen_Muon_Particle_Idx.begin(),v_mu);
      pair<int, unsigned short> pair = make_pair(igenmu, mother);
      v_GenMother.push_back(pair);
    
      cout << "DY,Z found " << evt << " gen mu " << igenmu
           << ": mother's idx " << mother << " and muon's particle idx " << *v_mu << "\n\t"
           << v_Gen_Particle_PdgId.at(*v_mu) << "/" << static_cast<unsigned short>(v_Gen_Particle_Status.at(*v_mu)) << " "
           << v_Gen_Particle_PdgId.at(mother) << "/" << static_cast<unsigned short>(v_Gen_Particle_Status.at(mother)) << " "
           << endl;

    } // end of DY, Z mother muon & genmu matching

    // check if gen muons have same mother
    int momid_check = -1;
    int igenmu_keep[2] = {-1,-1};
    for (auto v_genidx=v_GenMother.begin(); (v_GenMother.size()>=2) && v_genidx!=v_GenMother.end(); ++v_genidx) {
      
      pair<int, unsigned short> pair = *v_genidx;
      int igenmu = pair.first;
      int momid = pair.second;

      cout << "momid check " << momid_check << " " << momid << endl;
      if (momid_check==-1) {
        momid_check = momid;
        igenmu_keep[0] = distance(v_GenMother.begin(),v_genidx);
      } else if (momid_check==momid) { // 2 muons have same mother
        igenmu_keep[1] = distance(v_GenMother.begin(),v_genidx);
        break; // fill up histograms once per 1 event
      }
    } // end of 2 muons have same mother

    // Fill denominator with the highest pT gen muon
    if (igenmu_keep[0] != -1 && igenmu_keep[1] != -1) {
      pair<int, unsigned short> igenmu1 = v_GenMother.at(igenmu_keep[0]);
      pair<int, unsigned short> igenmu2 = v_GenMother.at(igenmu_keep[1]);
      TLorentzVector mu1 = Gen_Muon_Mom.at(igenmu1.first);
      TLorentzVector mu2 = Gen_Muon_Mom.at(igenmu2.first);
      if (mu1.Pt() > mu2.Pt()) hden_eta.Fill(mu1.Eta());
      else hden_eta.Fill(mu2.Eta());
    } else continue; // skip the event if this doesn't have 2 muons from a same mother




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
    
        cout<<Flag_collisionEventSelectionPA<<" "
        <<Flag_goodVertices<<" "
        <<Flag_CSCTightHaloFilter<<" "
        <<Flag_HBHENoiseFilter<<" "
        <<Flag_HBHENoiseIsoFilter<<" "
        <<Flag_EcalDeadCellTriggerPrimitiveFilter<<" "
        <<Flag_chargedHadronTrackResolutionFilter<<" "
        <<Flag_muonBadTrackFilter<<" "
        <<Flag_duplicateMuons<<" "
        <<Flag_badMuons<<endl;

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
    vector<bool> v_Reco_Muon_isTight = mutree.Reco_Muon_isTight();
    vector<char>  v_Reco_Muon_PF_Idx = mutree.Reco_Muon_PF_Idx();
    vector<float> v_PF_Muon_IsoPFR03NoPUCorr = mutree.PF_Muon_IsoPFR03NoPUCorr();

    // gen-reco matching with v_GenMother
    if (igenmu_keep[0] != -1 && igenmu_keep[1] != -1) {
      pair<int, unsigned short> p_igenmu1 = v_GenMother.at(igenmu_keep[0]);
      pair<int, unsigned short> p_igenmu2 = v_GenMother.at(igenmu_keep[1]);
      int igenmu1 = p_igenmu1.first;
      int igenmu2 = p_igenmu2.first;

      // use reco muons only when there's matched gen muon
      cout << "RECO : v_GenMother" << endl;

      char irecmu1 = v_Gen_Muon_Reco_Idx.at(igenmu1);
      char irecmu2 = v_Gen_Muon_Reco_Idx.at(igenmu2);
      cout << "irecmu: " << static_cast<short>(irecmu1) << " " << static_cast<short>(irecmu2) << endl;
      if (irecmu1<0 || irecmu2<0) {
        cout << "DY,Z eff: skip recmu cause both of reco muons have to be found" << endl;
        continue; // no matched gen-reco muon pair, skip reco step
      }

      char ipfmu1 = v_Reco_Muon_PF_Idx.at(irecmu1);
      char ipfmu2 = v_Reco_Muon_PF_Idx.at(irecmu2);
      cout << "ipfmu: " << static_cast<short>(ipfmu1) << " " << static_cast<short>(ipfmu2) << endl;
      if (ipfmu1<0 || ipfmu2<0) {
        cout << "DY,Z eff: skip pfmu cause both of pf muons have to be found" << endl;
        continue; // no matched gen-reco muon pair, skip reco step
      }

      // Check if trigger is fired
//      vector<bool> v_Event_Trig_Fired = mutree.Event_Trig_Fired();
//      cout << "after v_Event_Trig_Fired " << endl;
//      vector< vector< UChar_t > > vv_Pat_Muon_Trig = mutree.Pat_Muon_Trig();
//      cout << "after vv_Pat_Muon_Trig " << endl;
//      vector< UChar_t > v_Pat_Muon_Trig1 = vv_Pat_Muon_Trig.at(irecmu1);
//      vector< UChar_t > v_Pat_Muon_Trig2 = vv_Pat_Muon_Trig.at(irecmu2);
//      cout << "after v_Pat_Muon_Trig " << endl;
//
//      //5th trigger is L3Mu12
//      if (!v_Event_Trig_Fired.at(5) || !v_Pat_Muon_Trig1.at(5) || !v_Pat_Muon_Trig2.at(5)) continue;
     
      // Check if this muon passes quality cuts, kinematic selections
      TLorentzVector recmu1 = Reco_Muon_Mom.at(irecmu1);
      TLorentzVector recmu2 = Reco_Muon_Mom.at(irecmu2);
      TLorentzVector genmu1 = Gen_Muon_Mom.at(igenmu1);
      TLorentzVector genmu2 = Gen_Muon_Mom.at(igenmu2);
      cout << "At numerator, genmu1: " << igenmu1 << " " << genmu1.Eta() << " " << genmu1.Pt() << endl;
      cout << "At numerator, genmu2: " << igenmu2 << " " << genmu2.Eta() << " " << genmu2.Pt() << endl;

      // check the pT,eta of the highest pT muon only, not both of them 
      double highest_pt  = (recmu1.Pt()>recmu2.Pt()) ? recmu1.Pt()  : recmu2.Pt() ;
      double highest_eta = (recmu1.Pt()>recmu2.Pt()) ? recmu1.Eta() : recmu2.Eta() ;

      if ( (TMath::Abs(recmu1.Eta())<bins_eta[nbins_eta] && TMath::Abs(recmu2.Eta())<bins_eta[nbins_eta]) &&
           (v_PF_Muon_IsoPFR03NoPUCorr.at(ipfmu1)<0.15 && v_PF_Muon_IsoPFR03NoPUCorr.at(ipfmu2)<0.15) &&
           (recmu1.Pt()>=15 && recmu2.Pt()>=15 && highest_pt>=25) &&
           (v_Reco_Muon_isTight.at(irecmu1) && v_Reco_Muon_isTight.at(irecmu2))
         )
      { // events that will be vetoed
        cout << "pass analysis cut" << endl;
        hnum_eta[0].Fill(highest_eta);
      } else { // events that will NOT be vetoed
        cout << "DONT pass analysis cut" << endl;
        hnum_eta[1].Fill(highest_eta);
      }

    } // end of gen-reco matching with v_GenMother
  } // end of evt loop



  // Calculate efficiency
  CheckUnderFlow(hnum_eta[0],hden_eta);
  heff_eta[0].Divide(&hnum_eta[0],&hden_eta);
  CheckUnderFlow(hnum_eta[1],hden_eta);
  heff_eta[1].Divide(&hnum_eta[1],&hden_eta);


  // Write histograms into output file
  TFile foutput(outputFile.c_str(),"recreate");
  foutput.cd();
  hden_eta.Write();
  hnum_eta[0].Write();
  hnum_eta[1].Write();
  heff_eta[0].Write();
  heff_eta[1].Write();
  foutput.Close();


}

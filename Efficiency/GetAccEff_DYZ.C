#include "GetAccEff.h"
#include "HiMuonTree.h"
#include "HiMETTree.h"



//////////////////////////////////////
////////// M A I N ///////////////////
//////////////////////////////////////
void GetAccEff_DYZ(
    const string beamDir="pPb"
    ) {

  // Input/Output files
  string inputFile;
  if (beamDir.compare("pPb")==0) inputFile="root://eoscms//eos/cms/store/group/phys_heavyions/anstahll/EWQAnalysis2017/pPb2016/8160GeV/MC/Embedded/HiEWQForest_Embedded_DYToMuMu_pPb_8160GeV_20170518.root";
  else if (beamDir.compare("Pbp")==0) inputFile="root://eoscms//eos/cms/store/group/phys_heavyions/anstahll/EWQAnalysis2017/pPb2016/8160GeV/MC/Embedded/HiEWQForest_Embedded_DYToMuMu_Pbp_8160GeV_20170518.root";
  string outputFile=Form("DYZveto_%s.root",beamDir.c_str());

  // Define efficiency histograms
  // den[0]: all generated candidates
  // den[1]: all reconstructed candidates before applying veto cuts
  // num[0]: all reconstructed candidates before applying veto cuts
  // num[1]: all reconstructed candidates after veto cuts applied
  // num[2]: all reconstructed candidates after veto cuts applied and failed
  TH1D hden_eta[2], hnum_eta[3];
  TGraphAsymmErrors heff_eta[3];

  string suffix = beamDir;
  for (int i=0; i<3; i++) {
    if (i<2) hden_eta[i] = TH1D(Form("hden_eta_%s_%d",suffix.c_str(),i),";#eta_{lab};",nbins_eta,bins_eta);
    hnum_eta[i] = TH1D(Form("hnum_eta_%s_%d",suffix.c_str(),i),";#eta_{lab};",nbins_eta,bins_eta);
    heff_eta[i] = TGraphAsymmErrors();
    heff_eta[i].SetName(Form("heff_eta_%s_%d",suffix.c_str(),i));
    heff_eta[i].GetXaxis()->SetTitle("#eta_{lab}");
    heff_eta[i].GetYaxis()->SetTitle("");
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
  for (Long64_t evt=0; evt<50000; evt++) {
    if (mutree.GetEntry(evt)<0) break;

    /////// ###### GEN objects will be proceed!
    // Load gen particles, check if it is W
    vector< pair<int, unsigned short> > v_GenMother; //gen mu idx, its mother's idx
    vector<int> v_Gen_Particle_PdgId = mutree.Gen_Particle_PdgId();
    vector<unsigned char> v_Gen_Particle_Status = mutree.Gen_Particle_Status();
    vector < vector < UShort_t > > vv_Gen_Particle_Mother_Idx = mutree.Gen_Particle_Mother_Idx();
    vector<char> v_Gen_Muon_Reco_Idx = mutree.Gen_Muon_Reco_Idx();
    vector<char> v_Gen_Muon_PF_Idx = mutree.Gen_Muon_PF_Idx();
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
    char highest_pt_genmu = -1;
    if (igenmu_keep[0] != -1 && igenmu_keep[1] != -1) {
      pair<int, unsigned short> igenmu1 = v_GenMother.at(igenmu_keep[0]);
      pair<int, unsigned short> igenmu2 = v_GenMother.at(igenmu_keep[1]);
      TLorentzVector mu1 = Gen_Muon_Mom.at(igenmu1.first);
      TLorentzVector mu2 = Gen_Muon_Mom.at(igenmu2.first);
      highest_pt_genmu = (mu1.Pt()>mu2.Pt()) ? igenmu1.first : igenmu2.first ;
      
      if (mu1.Pt() > mu2.Pt()) hden_eta[0].Fill(mu1.Eta());
      else hden_eta[0].Fill(mu2.Eta());
    
    } else continue; // skip the event if this doesn't have 2 muons from a same mother


    /////// ###### RECO objects will be proceed!
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

    // Check if trigger is fired
    vector<bool> v_Event_Trig_Fired = mutree.Event_Trig_Fired();
    if (v_Event_Trig_Fired.size()<=0) continue;
    cout << "after v_Event_Trig_Fired " << endl;
    vector< vector< UChar_t > > vv_Pat_Muon_Trig = mutree.Pat_Muon_Trig();
    if (vv_Pat_Muon_Trig.size()<=0) continue;
    cout << "after vv_Pat_Muon_Trig " << endl;
    
    // RECO muons to be checked
    vector<TLorentzVector> Reco_Muon_Mom = mutree.Reco_Muon_Mom();
    vector<char> v_Reco_Muon_Gen_Idx = mutree.Reco_Muon_Gen_Idx();
    vector<bool> v_Reco_Muon_isTight = mutree.Reco_Muon_isTight();
    vector<char>  v_Reco_Muon_PF_Idx = mutree.Reco_Muon_PF_Idx();
    vector<float> v_PF_Muon_IsoPFR03NoPUCorr = mutree.PF_Muon_IsoPFR03NoPUCorr();

    // Check the highest pT reco muon and fill den[1] and num[0]
    char irecmu0 = v_Gen_Muon_Reco_Idx.at(highest_pt_genmu);
    char ipfmu0 = v_Gen_Muon_PF_Idx.at(highest_pt_genmu);
    cout << "irecmu0, ipfmu0: " << static_cast<short>(irecmu0) << " " << static_cast<short>(ipfmu0) << endl;
 
    // gen-reco matching with v_GenMother
    pair<int, unsigned short> p_igenmu1 = v_GenMother.at(igenmu_keep[0]);
    pair<int, unsigned short> p_igenmu2 = v_GenMother.at(igenmu_keep[1]);
    int igenmu1 = p_igenmu1.first;
    int igenmu2 = p_igenmu2.first;
    char irecmu1 = v_Gen_Muon_Reco_Idx.at(igenmu1);
    char irecmu2 = v_Gen_Muon_Reco_Idx.at(igenmu2);
    char ipfmu1 = v_Gen_Muon_PF_Idx.at(igenmu1);
    char ipfmu2 = v_Gen_Muon_PF_Idx.at(igenmu2);
    cout << "irecmu: " << static_cast<short>(irecmu1) << " " << static_cast<short>(irecmu2) << endl;
    cout << "ipfmu: " << static_cast<short>(ipfmu1) << " " << static_cast<short>(ipfmu2) << endl;
    
    if ((irecmu1<0 && irecmu2>=0) && highest_pt_genmu == igenmu1) {
      irecmu0 = irecmu2;
      ipfmu0 = ipfmu2;
    } else if ((irecmu1>=0 && irecmu2<0) && highest_pt_genmu == igenmu2) {
      irecmu0 = irecmu1;
      ipfmu0 = ipfmu1;
    }
    if (irecmu0<0) continue;

    vector< UChar_t > v_Pat_Muon_Trig0 = vv_Pat_Muon_Trig.at(irecmu0);
    //5th trigger is L3Mu12
    if (!v_Pat_Muon_Trig0.at(5)) continue;
    cout << "after v_Pat_Muon_Trig " << endl;

    TLorentzVector recmu0= Reco_Muon_Mom.at(irecmu0);
    if ( TMath::Abs(recmu0.Eta())>bins_eta[nbins_eta] || recmu0.Pt()<25 ||
         !v_Reco_Muon_isTight.at(irecmu0) ||
         (v_Reco_Muon_isTight.at(irecmu0) && v_PF_Muon_IsoPFR03NoPUCorr.at(ipfmu0)>0.15)
       ) continue;

    TLorentzVector genmu0 = Gen_Muon_Mom.at(highest_pt_genmu);
    // all reconstructed candidates before applying veto cuts
    hnum_eta[0].Fill(genmu0.Eta());
    hden_eta[1].Fill(genmu0.Eta());

    // Check if both muons pass veto contidions, kinematic selections
    if (irecmu1<0 || irecmu2<0) { // when only 1 muon survived, fill up num[2]
      hnum_eta[2].Fill(genmu0.Eta());
      continue;
    }
    
    bool passed_event = false;
    
    TLorentzVector recmu1 = Reco_Muon_Mom.at(irecmu1);
    TLorentzVector recmu2 = Reco_Muon_Mom.at(irecmu2);
    cout << "recmu1: " << recmu1.Eta() << " " << recmu1.Pt() << endl;
    cout << "recmu2: " << recmu2.Eta() << " " << recmu2.Pt() << endl;
  
    if ( (TMath::Abs(recmu1.Eta())<bins_eta[nbins_eta] && TMath::Abs(recmu2.Eta())<bins_eta[nbins_eta]) &&
         (recmu1.Pt()>=15 && recmu2.Pt()>=15) &&
         (v_Reco_Muon_isTight.at(irecmu1) && v_Reco_Muon_isTight.at(irecmu2))
       )
    { 
      if (v_PF_Muon_IsoPFR03NoPUCorr.at(ipfmu1)<0.15 && v_PF_Muon_IsoPFR03NoPUCorr.at(ipfmu2)<0.15) {
        // events that will be vetoed
        cout << "pass veto cuts" << endl;
        passed_event = true;
      }
    } 

    if (passed_event) hnum_eta[1].Fill(genmu0.Eta());
    else hnum_eta[2].Fill(genmu0.Eta());

  } // end of evt loop



  // Calculate efficiency
  for (int i=0; i<2; i++) {
    CheckUnderFlow(hnum_eta[i],hden_eta[i]);
    heff_eta[i].Divide(&hnum_eta[i],&hden_eta[i]);
  }
  CheckUnderFlow(hnum_eta[2],hden_eta[1]);
  heff_eta[2].Divide(&hnum_eta[2],&hden_eta[1]);

  // Write histograms into output file
  TFile foutput(outputFile.c_str(),"recreate");
  foutput.cd();
  for (int i=0; i<3; i++) {
    if (i<2) hden_eta[i].Write();
    hnum_eta[i].Write();
    heff_eta[i].Write();
  }
  foutput.Close();


}

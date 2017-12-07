#define muEventAna_cxx
#include "muEventAna.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include <TTree.h>

// settings
const int itrig = 5; // HLT_PAL3Mu12

void muEventAna::Loop(const char* filename)
{
   //   In a ROOT session, you can do:
   //      root> .L muEventAna.C
   //      root> muEventAna t
   //      root> t.GetEntry(12); // Fill t data members with entry number 12
   //      root> t.Show();       // Show values of entry 12
   //      root> t.Show(16);     // Read and show values of entry 16
   //      root> t.Loop();       // Loop on all entries
   //

   //     This is the loop skeleton where:
   //    jentry is the global entry number in the chain
   //    ientry is the entry number in the current Tree
   //  Note that the argument to GetEntry must be:
   //    jentry for TChain::GetEntry
   //    ientry for TTree::GetEntry and TBranch::GetEntry
   //
   //       To read only selected branches, Insert statements like:
   // METHOD1:
   //    fChain->SetBranchStatus("*",0);  // disable all branches
   //    fChain->SetBranchStatus("branchname",1);  // activate branchname
   // METHOD2: replace line
   //    fChain->GetEntry(jentry);       //read all branches
   //by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   TH1::SetDefaultSumw2(kTRUE);

   Long64_t nentries = fChain->GetEntriesFast();

   TFile *f = new TFile(filename,"RECREATE");

   double mass_sig_min = 80;
   double mass_sig_max = 100;

   const int nvars = 45;
   TString varnames[nvars] = {
      "PF_Muon_IsoPFR03","PF_Muon_IsoPFR03NoPUCorr","PF_Muon_IsoPFR04","PF_Muon_IsoPFR04NoPUCorr","PF_Muon_EM_Chg_sumR03Pt",
      "PF_Muon_EM_Chg_sumR04Pt","PF_Muon_EM_Neu_sumR03Pt","PF_Muon_EM_Neu_sumR04Pt","PF_Muon_Had_Chg_sumR03Pt","PF_Muon_Had_Chg_sumR04Pt",
      "PF_Muon_Had_Neu_sumR03Pt","PF_Muon_Had_Neu_sumR04Pt","PF_Muon_Had_PU_sumR03Pt","PF_Muon_Had_PU_sumR04Pt","PF_DiMuon_VtxProb",
      "PF_Muon_IsoPFR03_reliso","PF_Muon_IsoPFR03NoPUCorr_reliso","PF_Muon_IsoPFR04_reliso","PF_Muon_IsoPFR04NoPUCorr_reliso","PF_Muon_EM_Chg_sumR03Pt_reliso",
      "PF_Muon_EM_Chg_sumR04Pt_reliso","PF_Muon_EM_Neu_sumR03Pt_reliso","PF_Muon_EM_Neu_sumR04Pt_reliso","PF_Muon_Had_Chg_sumR03Pt_reliso","PF_Muon_Had_Chg_sumR04Pt_reliso",
      "PF_Muon_Had_Neu_sumR03Pt_reliso","PF_Muon_Had_Neu_sumR04Pt_reliso","PF_Muon_Had_PU_sumR03Pt_reliso","PF_Muon_Had_PU_sumR04Pt_reliso", "Reco_Muon_IsoR03_reliso",
      "Reco_Muon_IsoR05", "Reco_Muon_Trk_sumR03Pt", "Reco_Muon_Trk_sumR05Pt", "Reco_Muon_IsoR03", "Reco_Muon_IsoR05_reliso",
      "Reco_Muon_Trk_sumR03Pt_reliso", "Reco_Muon_Trk_sumR05Pt_reliso", "PF_Muon_myIsoPFR010", "PF_Muon_myIsoPFR015", "PF_Muon_myIsoPFR020",
      "PF_Muon_myIsoPFR025", "PF_Muon_myIsoPFR030", "PF_Muon_myIsoPFR035", "PF_Muon_myIsoPFR040", "PF_Muon_myIsoPFR045"
   };
   int nbins[nvars] = {100,100,100,100,100,
      100,100,100,100,100,
      100,100,100,100,100,
      100,100,100,100,100,
      100,100,100,100,100,
      100,100,100,100,100,
      100,100,100,100,100,
      100,100,100,100,100,
      100,100,100,100,100
   };
   double varmin[nvars] = {0,0,0,0,0,
      0,0,0,0,0,
      0,0,0,0,0,
      0,0,0,0,0,
      0,0,0,0,0,
      0,0,0,0,0,
      0,0,0,0,0,
      0,0,0,0,0,
      0,0,0,0,0
   };
   double varmax[nvars] = {1,1,1,1,20,
      20,20,20,100,100,
      100,100,100,100,1,
      1,1,1,1,1,
      1,1,1,1,1,
      1,1,1,1,20,
      20,20,20,1,1,
      1,1,1,1,1,
      1,1,1,1,1
   };

   TH1F *hist_OS[nvars] = {0};
   TH1F *hist_SS[nvars] = {0};
   TH1F *hist_bkg[nvars] = {0};

   for (int i=0; i<nvars; i++) {
      TString hname = "hist_" + varnames[i] + "_OS";
      hist_OS[i] = new TH1F(hname,varnames[i]+";"+varnames[i]+";Entries",nbins[i],varmin[i],varmax[i]);
      hname = "hist_" + varnames[i] + "_SS";
      hist_SS[i] = new TH1F(hname,varnames[i]+";"+varnames[i]+";Entries",nbins[i],varmin[i],varmax[i]);
      hname = "hist_" + varnames[i] + "_bkg";
      hist_bkg[i] = new TH1F(hname,varnames[i]+";"+varnames[i]+";Entries",nbins[i],varmin[i],varmax[i]);
   }

   // global plots
   TH1F *hm_OS = new TH1F("hm_OS","mass (OS)",200,0,200);
   TH1F *hm_SS = new TH1F("hm_SS","mass (SS)",200,0,200);
   TH1F *hprescale = new TH1F("hprescale","prescale",100,0,100);

   // disable branches
   fChain->SetBranchStatus("*",0);
   fChain->SetBranchStatus("Event_Trig*",1);
   fChain->SetBranchStatus("Event_nPV",1);
   fFriend1->fChain->SetBranchStatus("*",0);
   fFriend1->fChain->SetBranchStatus("*Mom",1);
   fFriend1->fChain->SetBranchStatus("*Charge",1);
   fFriend1->fChain->SetBranchStatus("*Idx",1);
   fFriend1->fChain->SetBranchStatus("PF_Candidate*",1);
   fFriend2->fChain->SetBranchStatus("*",0);
   fFriend2->fChain->SetBranchStatus("*Mom",1);
   fFriend2->fChain->SetBranchStatus("*Charge",1);
   fFriend2->fChain->SetBranchStatus("*Idx",1);
   fFriend2->fChain->SetBranchStatus("*Tight",1);
   fFriend3->fChain->SetBranchStatus("*",0);
   fFriend3->fChain->SetBranchStatus("PF_MET_Mom",1);
   for (int ivar=0; ivar<nvars; ivar++) {
      fChain->SetBranchStatus(varnames[ivar],1);
      fFriend1->fChain->SetBranchStatus(varnames[ivar],1);
      fFriend2->fChain->SetBranchStatus(varnames[ivar],1);
   }

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;

      // trigger selection
      b_Event_Trig_Fired->GetEntry(jentry);
      if (!Event_Trig_Fired->at(itrig)) continue;

      // nPV selection
      b_Event_nPV->GetEntry(jentry);
      if ((jentry%100000)==0) cout << (int) Event_nPV << " " << jentry << endl;
      //if (Event_nPV != 1) continue;

      nb = fChain->GetEntry(jentry);   nbytes += nb;
      fFriend1->GetEntry(jentry);
      fFriend2->GetEntry(jentry);
      fFriend3->GetEntry(jentry);
      // if (Cut(ientry) < 0) continue;
      if ( jentry % 100000 == 0 ) std::cout << jentry << "/" << nentries << std::endl;
      // if (jentry>100000) break;

      double weight = Event_Trig_Presc->at(itrig);
      hprescale->Fill(weight);

      // dimuons (OS and SS)
      for (unsigned int i=0; i<fFriend1->PF_DiMuon_Charge->size(); i++) {
         // selection
         // acceptance
         int pfidx1 = fFriend1->PF_DiMuon_Muon1_Idx->at(i);
         int pfidx2 = fFriend1->PF_DiMuon_Muon2_Idx->at(i);
         TLorentzVector *tlvpf1 = (TLorentzVector*) fFriend1->PF_Muon_Mom->At(pfidx1);
         TLorentzVector *tlvpf2 = (TLorentzVector*) fFriend1->PF_Muon_Mom->At(pfidx2);
         if (!(IsAccept(tlvpf1->Pt(), tlvpf1->Eta()) && IsAccept(tlvpf2->Pt(), tlvpf2->Eta()))) continue;
         
         // ID
         int idx1 = fFriend1->PF_Muon_Reco_Idx->at(pfidx1);
         int idx2 = fFriend1->PF_Muon_Reco_Idx->at(pfidx2);
         if (!(fFriend2->Reco_Muon_isTight->at(idx1) && fFriend2->Reco_Muon_isTight->at(idx2))) continue;
         
         // charge
         bool isOS = (fFriend1->PF_DiMuon_Charge->at(i) == 0);

         TLorentzVector *tlvmumu = (TLorentzVector*) fFriend1->PF_DiMuon_Mom->At(i);
         double mass = tlvmumu->M();
         if (isOS) hm_OS->Fill(mass,weight);
         else hm_SS->Fill(mass,weight);

         if (mass<mass_sig_min || mass>mass_sig_max) continue;

         // ok, this is a good dimuon. Now fill histos.
         for (int ivar=0; ivar<nvars; ivar++) {
            double var=0; int i1=-1,i2=-1;
            if (varnames[ivar].Contains("PF_DiMuon")) i1=i;
            if (varnames[ivar].Contains("Reco_Muon")) {i1=idx1; i2=idx2;}
            if (varnames[ivar].Contains("PF_Muon")) {i1=pfidx1; i2=pfidx2;}

            var = Get(varnames[ivar].ReplaceAll("_reliso",""),i1);
            if (varnames[ivar].Contains("_reliso")) var = var / tlvpf1->Pt();
            if (isOS) hist_OS[ivar]->Fill(var,weight);
            else hist_SS[ivar]->Fill(var,weight);
            if (i2>=0) {
               var = Get(varnames[ivar].ReplaceAll("_reliso",""),i2);
               if (varnames[ivar].Contains("_reliso")) var = var / tlvpf2->Pt();
               if (isOS) hist_OS[ivar]->Fill(var,weight);
               else hist_SS[ivar]->Fill(var,weight);
            }
         } // var loop

         // we found one good dimuon... skip the rest of the event
         break;
      } // PF dimuon loop

      // if there is only one muon and it is low MET: consider it for bkg
      if (fFriend1->PF_Muon_Charge->size()==1 && fFriend3->PF_MET_Mom->Mod()<5) {
         // selection
         // acceptance
         TLorentzVector *tlvpf1 = (TLorentzVector*) fFriend1->PF_Muon_Mom->At(0);
         if (!IsAccept(tlvpf1->Pt(), tlvpf1->Eta())) continue;
         
         // ID
         int idx1 = fFriend1->PF_Muon_Reco_Idx->at(0);
         if (!(fFriend2->Reco_Muon_isTight->at(idx1))) continue;

         // ok, this is a good dimuon. Now fill histos.
         for (int ivar=0; ivar<nvars; ivar++) {
            int i1=0;
            if (varnames[ivar].Contains("PF_DiMuon")) continue;
            if (varnames[ivar].Contains("Reco_Muon")) {i1=idx1;}

            double var = Get(varnames[ivar].ReplaceAll("_reliso",""),i1);
            if (varnames[ivar].Contains("_reliso")) var = var / tlvpf1->Pt();
            hist_bkg[ivar]->Fill(var,weight);
         } // var loop
      } // if single muon && MET<5
   } // event loop

   f->Write();
   f->Close();
   if (f) delete f;

}


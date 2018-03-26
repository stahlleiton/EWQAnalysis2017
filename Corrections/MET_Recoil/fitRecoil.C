//================================================================================================
//
// Perform fits to recoil against Z->mumu events
//
//  * Outputs a ROOT file of functions parametrizing the distribution of recoil components
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "../../Utilities/HiMuonTree.h"
#include "../../Utilities/HiMETTree.h"
#include "../../Utilities/HiEvtTree.h"
#include "../../Utilities/EVENTUTILS.h"
#include "../../Utilities/HFweight.h"
#include "../../Fitter/Macros/Utilities/initClasses.h"
#include "../../Utilities/CMS/tdrstyle.C"
#include "../../Utilities/CMS/CMS_lumi.C"
#include <iostream>                   // standard I/O
#include <fstream>                    // standard I/O
#include <TFile.h>                    // file handle class
#include <TTree.h>                    // class to access ntuples
#include <TF1.h>                      // 1D function
#include <TFitResult.h>               // class to handle fit results
#include <TGraphAsymmErrors.h>        // graph class
#include <TLorentzVector.h>           // 4-vector class
#include <TH1D.h>                     // plots
#include <TH2D.h>                     // plots
#include <TCanvas.h>                  // canvas
#include <TSystem.h>
#include <TPaveText.h>

#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "RooGaussian.h"
#include "RooDataSet.h"
#include "RooAbsPdf.h"
#include "RooHist.h"
#include "RooAddPdf.h"
#include "RooHistPdf.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooDataHist.h"
#include "RooFormulaVar.h"
#endif

//=== FUNCTION DECLARATIONS ======================================================================================


// print the chi2 result
bool printChi2(TPad& pad, const RooWorkspace& ws, const RooPlot& frame, const string& varLabel, const RooDataSet& dataset, const string& pdfLabel);

// Axis and Text Utility functions
void updateYAxisRange ( RooPlot& frame, const RooWorkspace& myws, std::string varLabel, 
                        const RooDataSet& dataset, const bool& logScale );
void updateYAxisRange ( std::unique_ptr<TGraphAsymmErrors>& graph );

std::string formatText(const std::string& text);


// generate web page
void makeHTML(const std::string outDir,
              const std::map< std::string , std::unique_ptr<TGraphAsymmErrors> >& u1Graph,
              const std::map< std::string , std::unique_ptr<TGraphAsymmErrors> >& u2Graph,
              const std::string uparName, const std::string uprpName, const uint nbins);


//--------------------------------------------------------------------------------------------------
// perform fit of recoil component
bool performFit(
                RooWorkspace& ws,
                const std::vector< std::unique_ptr<TH1D> >& hv,
                const std::vector< double >& ptBins,
                const std::string uName,
                const std::string rapRange,
                const uint model,
                const std::string dsName,
                const std::string col,
                const bool isData,
                TCanvas& c,
                std::map< std::string , std::vector< std::vector< double > > >& varArr,
                const std::string outputDir
                );


//=== MAIN MACRO ================================================================================================= 

void fitRecoil(
               const bool isData = false,
               uint pfumodel = 2, // u1/2 model (1 => single Gaussian, 2 => double Gaussian, 3 => triple Gaussian, 4 => Breit-Wigner plus Gaussian)
               const int applyHFCorr = 0, // Only used for MC, in data it is ignored. 1: hiHF , 2: NTracksv
               const bool applyBosonPTCorr = false, // Only used for MC, in data it is ignored.
               const std::string yRange = "full", // Z rapidity (full => [-2.4,2.4], mid => [-1.6,1.6], fwd => |y| < 1.6)
               const std::vector< std::string > metType = { "PF_RAW"/*,"PF_RAW_JetEnDown","PF_RAW_JetEnUp","PF_RAW_JetResDown","PF_RAW_JetResUp","PF_RAW_MuonEnDown","PF_RAW_MuonEnUp" , "PF_Type1" , "PF_NoHF_RAW", "PF_NoHF_Type1" */},
               const std::vector< std::string > COLL    = { "PA" , "Pbp" , "pPb" },
               const bool remakeDS = false
               )
{
  //
  const std::string uparName   = "u1";
  const std::string uprpName   = "u2";
  //
  if (pfumodel>4 || pfumodel<1){
    std::cout << "[ERROR] The supported models are: 1 => single Gaussian, 2 => double Gaussian, 3 => triple Gaussian" << std::endl;
    return;
  }
  //
  if (yRange!="full" && yRange!="mid" && yRange!="fwd"){
    std::cout << "[ERROR] The supported Z rapidity ranges are: full => [-2.4,2.4], mid => [-1.6,1.6], fwd => |y| < 1.6" << std::endl;
    return;
  }
  //
  // Change the working directory
  const std::string CWD = getcwd(NULL, 0);
  const std::string mainDir = Form("%s/FitRecoil/", CWD.c_str());
  gSystem->mkdir(Form("%s/FitRecoil", CWD.c_str()), kTRUE);
  gSystem->ChangeDirectory(Form("%s/FitRecoil", CWD.c_str()));
  

  // RooFit, you are annoying, so just shut up! XD
  RooMsgService::instance().getStream(0).removeTopic(RooFit::Caching);  
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Caching);
  RooMsgService::instance().getStream(0).removeTopic(RooFit::Plotting);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Plotting);
  RooMsgService::instance().getStream(0).removeTopic(RooFit::Integration);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Integration);
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;
  
  std::vector< std::string > fileName;
  std::string dsLabel;
  uint pfu1model = pfumodel; // u1 model (1 => single Gaussian, 2 => double Gaussian, 3 => triple Gaussian, 4 => Breit-Wigner plus Gaussian)
  uint pfu2model = pfumodel; // u2 model (1 => single Gaussian, 2 => double Gaussian, 3 => triple Gaussian, 4 => Breit-Wigner plus Gaussian)
  if (isData) {
//    fileName.push_back("/store/group/phys_heavyions/anstahll/EWQAnalysis2017/pPb2016/8160GeV/Data/HiEWQForest_PASingleMuon_pPb_Pbp_8160GeV_20171003.root");
    fileName.push_back("/home/llr/cms/blanco/Analysis/WAnalysis/DATASETS/DATA/HiEWQForest_PASingleMuon_pPb_Pbp_8160GeV_20171003.root");
    dsLabel = "DATA";
//    pfu1model = 1;
    //pfu2model = 1;
  }
  else {
//    fileName.push_back("/store/group/phys_heavyions/anstahll/EWQAnalysis2017/pPb2016/8160GeV/MC/Embedded/Official/POWHEG/HiEWQForest_Embedded_Official_POWHEG_CT14_EPPS16_DYtoMuMu_M_30_pPb_8160GeV_20171003.root");
//    fileName.push_back("/store/group/phys_heavyions/anstahll/EWQAnalysis2017/pPb2016/8160GeV/MC/Embedded/Official/POWHEG/HiEWQForest_Embedded_Official_POWHEG_CT14_EPPS16_DYtoMuMu_M_30_Pbp_8160GeV_20171003.root");
    fileName.push_back("/home/llr/cms/blanco/Analysis/WAnalysis/DATASETS/MC/HiEWQForest_Embedded_Official_POWHEG_CT14_EPPS16_DYtoMuMu_M_30_pPb_8160GeV_20171003.root");
    fileName.push_back("/home/llr/cms/blanco/Analysis/WAnalysis/DATASETS/MC/HiEWQForest_Embedded_Official_POWHEG_CT14_EPPS16_DYtoMuMu_M_30_Pbp_8160GeV_20171003.root");
    dsLabel = "MC_DYToMuMu_POWHEG";
//    pfu1model = 2;
    //pfu2model = 2;
  }
 
  // Define Boson pT Binning
  std::vector< double > ptBins;
  std::string rapRange;
  if (yRange=="full") {
    ptBins = { 0., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 19., 21., 23., 25., 27., 30., 35., 40., 45., 50, 60., 70., 90., 140.};
    rapRange = "[-2.4,2.4]";
  }
  else {
    ptBins = { 0., 3., 5., 7., 9., 12., 15., 20., 30., 40., 60., 140.};
    if (yRange=="mid") rapRange = "[-1.6,1.6]";
    if (yRange=="fwd") rapRange = "< 1.6";
  }
  const uint nbins = (ptBins.size()-1);

  // Define the Boson Kinematic Cuts
  const Double_t MASS_LOW  = 80;
  const Double_t MASS_HIGH = 110;
  const Double_t PT_CUT    = 15.;
  const Double_t ETA_CUT   = 2.4;

  // Define Trigger Event Selection
  const int triggerIndex = PA::HLT_PAL3Mu12;

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  

  // Initialize the histograms
  std::map< std::string , std::map< std::string , std::vector< std::unique_ptr<TH1D> > > > hPFu1v , hPFu2v;
  for (const auto& met : metType) { for (const auto& col : COLL) { hPFu1v[met][col].resize(nbins); hPFu2v[met][col].resize(nbins); } }
  //
  // Initialize the RooDataSets
  std::map< std::string , std::map< std::string , RooWorkspace > > myws_u1 , myws_u2;
  //
  // Check if DataSet have already been made
  //
  bool makeDS = false;
  //
  if (remakeDS==false) {
    // Check the directory
    const std::string dsDir = Form("%s/DataSet/", CWD.c_str());
    const std::string inputDSDir = dsDir + dsLabel +"/";
    if (existDir(inputDSDir)) {
      // Open the input file
      for (const auto& met : metType) {
        for (const auto& col : COLL) {
          std::string label = "";
          if (!isData) {
            if (applyHFCorr==1) { label += "_HFCorr";  } else if (applyHFCorr==2) { label += "_NTrack"; } else { label += ""; }
            if (applyBosonPTCorr) { label += "_BosonPT"; }
          }
          const std::string inDSname = (inputDSDir + Form("DATASET_%s_%s%s.root", met.c_str(), col.c_str(), label.c_str()));
          if (existFile(inDSname)) {
            bool prodNotFound = false;
            std::cout << "[INFO] Extracting DataSets from " << inDSname << endl;
            auto inDSfile = std::unique_ptr<TFile>(new TFile(inDSname.c_str(), "READ"));
            if (inDSfile!=NULL && inDSfile->IsOpen() && !inDSfile->IsZombie()) {
              if (inDSfile->Get("workspace_u1")) { myws_u1[met][col].import(*((RooWorkspace*)inDSfile->Get("workspace_u1"))->data(Form("dPF%s_%s", uparName.c_str(), dsLabel.c_str()))); } else { prodNotFound = true; }
              if (inDSfile->Get("workspace_u2")) { myws_u2[met][col].import(*((RooWorkspace*)inDSfile->Get("workspace_u2"))->data(Form("dPF%s_%s", uprpName.c_str(), dsLabel.c_str()))); } else { prodNotFound = true; }
              for (uint ibin = 0; ibin < nbins; ibin++) {
                if (inDSfile->Get(Form("hPFu1_%i", ibin))) {
                  hPFu1v.at(met).at(col).at(ibin) = std::unique_ptr<TH1D>((TH1D*)inDSfile->Get(Form("hPFu1_%i", ibin))); hPFu1v.at(met).at(col).at(ibin)->SetDirectory(0);
                } else { prodNotFound = true; }
                if (inDSfile->Get(Form("hPFu2_%i", ibin))) {
                  hPFu2v.at(met).at(col).at(ibin).reset((TH1D*)inDSfile->Get(Form("hPFu2_%i", ibin))); hPFu2v.at(met).at(col).at(ibin)->SetDirectory(0);
                } else { prodNotFound = true; }
              }
            }
            else { std::cout << "[INFO] Input DataSet file " << inDSname << " fail to open, will make the datasets" << std::endl; makeDS = true; }
            if (prodNotFound) { std::cout << "[WARNING] Input DataSet file is corrupt, will remake it!" << std::endl; makeDS = true; }
            if (inDSfile!=NULL) { inDSfile->Close(); }
          }
          else { std::cout << "[INFO] Input DataSet file " << inDSname << " was not found, will make the datasets" << std::endl; makeDS = true; }
          if (makeDS) break;
        }
        if (makeDS) break;
      }
    }
    else { std::cout << "[INFO] Input DataSet Directory " << inputDSDir << " was not found, will make the datasets" << std::endl; makeDS = true; }
  }
  //
  // Make the datasets
  //
  if (makeDS) {
    //
    for (const auto& met : metType) {
      for (const auto& col : COLL) {
        for(uint ibin = 0; ibin < nbins; ibin++) {
          // For u1
          hPFu1v.at(met).at(col).at(ibin).reset(new TH1D(Form("hPFu1_%s_%s_%i", met.c_str(), col.c_str(), ibin), "", 200, -200-ptBins[ibin], 200-ptBins[ibin]));
          hPFu1v.at(met).at(col).at(ibin)->Sumw2();
          // For u2
          hPFu2v.at(met).at(col).at(ibin).reset(new TH1D(Form("hPFu2_%s_%s_%i", met.c_str(), col.c_str(), ibin), "", 200, -200, 200));
          hPFu2v.at(met).at(col).at(ibin)->Sumw2();
        }
      }
    }
    //
    for (const auto& met : metType) {
      for (const auto& col : COLL) {
        //
        RooRealVar ptVar  = RooRealVar("pt", "pt", -1000., 1000., "GeV/c");
        RooRealVar weight = RooRealVar("weight", "weight", -1.0, 1000000000.0, "");
        RooRealVar u1Var  = RooRealVar(uparName.c_str(), uparName.c_str(), -1000., 1000., "GeV/c");
        std::unique_ptr<RooDataSet> ds;
        if (!isData && (applyHFCorr>0 || applyBosonPTCorr)) {
          ds.reset(new RooDataSet(Form("dPF%s_%s", uparName.c_str(), dsLabel.c_str()), "", RooArgSet(u1Var, ptVar, weight), RooFit::WeightVar(weight)));
        }
        else { ds.reset(new RooDataSet(Form("dPF%s_%s", uparName.c_str(), dsLabel.c_str()), "", RooArgSet(u1Var, ptVar))); }
        myws_u1[met][col].import(*ds);
        //
        RooRealVar u2Var  = RooRealVar(uprpName.c_str(), uprpName.c_str(), -1000., 1000., "GeV/c");
        if (!isData && (applyHFCorr>0 || applyBosonPTCorr)) {
          ds.reset(new RooDataSet(Form("dPF%s_%s", uprpName.c_str(), dsLabel.c_str()), "", RooArgSet(u2Var, ptVar, weight), RooFit::WeightVar(weight)));
        }
        else { ds.reset(new RooDataSet(Form("dPF%s_%s", uprpName.c_str(), dsLabel.c_str()), "", RooArgSet(u2Var, ptVar))); }
        myws_u2[met][col].import(*ds);
      }
    }
    //
    // Define the HF Weight
    std::unique_ptr<HFweight> HFCorr;
    if (!isData && applyHFCorr>0) { HFCorr = std::unique_ptr<HFweight>(new HFweight("/afs/cern.ch/work/e/echapon/public/DY_pA_2016/HFweight.root")); }
    //
    // Proceed to loop over each input file
    for(const auto& file : fileName) {
      std::cout << "[INFO] Processing " << file << endl;
      //
      std::unique_ptr<HiMuonTree> muonTree = std::unique_ptr<HiMuonTree>(new HiMuonTree());
      if (!muonTree->GetTree(file)) return;
      Long64_t nentries = muonTree->GetEntries();
      //
      std::unique_ptr<HiEvtTree> evtTree = std::unique_ptr<HiEvtTree>(new HiEvtTree());
      if (!evtTree->GetTree(file)) return;
      if (evtTree->GetEntries() != nentries) { std::cout << "[ERROR] Inconsistent number of entries!" << std::endl; return; }
      //
      std::unique_ptr<HiMETTree> metTree = std::unique_ptr<HiMETTree>(new HiMETTree());
      if (!metTree->GetTree(file, "metAna")) return;
      if (metTree->GetEntries() != nentries) { std::cout << "[ERROR] Inconsistent number of entries!" << std::endl; return; }
      //
      std::unique_ptr<HiMETTree> metNoHFTree = std::unique_ptr<HiMETTree>(new HiMETTree());
      if (!metNoHFTree->GetTree(file, "metAnaNoHF")) return;
      if (metNoHFTree->GetEntries() != nentries) { std::cout << "[ERROR] Inconsistent number of entries!" << std::endl; return; }
      //
      // Loop over events
      //
      for (Long64_t jentry=0; jentry<nentries; jentry++) {
        //
        // Get the entry in the trees
        if (   muonTree->GetEntry(jentry)<0) { std::cout << "[ERROR] Muon Tree invalid entry!"     << std::endl; return; }
        if (    evtTree->GetEntry(jentry)<0) { std::cout << "[ERROR] Event Tree invalid entry!"    << std::endl; return; }
        if (    metTree->GetEntry(jentry)<0) { std::cout << "[ERROR] MET Tree invalid entry!"      << std::endl; return; }
        if (metNoHFTree->GetEntry(jentry)<0) { std::cout << "[ERROR] MET NoHF Tree invalid entry!" << std::endl; return; }
        // 
        // Check that the different tree agrees well
        if (muonTree->Event_Run()!=metTree->Event_Run()          ) { std::cout << "[ERROR] MET Run does not agree!"        << std::endl; return; }
        if (muonTree->Event_Number()!=metTree->Event_Number()    ) { std::cout << "[ERROR] MET Event does not agree!"      << std::endl; return; }
        if (muonTree->Event_Run()!=metNoHFTree->Event_Run()      ) { std::cout << "[ERROR] MET NoHF Run does not agree!"   << std::endl; return; }
        if (muonTree->Event_Number()!=metNoHFTree->Event_Number()) { std::cout << "[ERROR] MET NoHF Event does not agree!" << std::endl; return; }
        if (muonTree->Event_Run()!=evtTree->run()                ) { std::cout << "[ERROR] HiEVT Run does not agree!"      << std::endl; return; }
        if (muonTree->Event_Number()!=evtTree->evt()             ) { std::cout << "[ERROR] HiEVT Event does not agree!"    << std::endl; return; }
        //
        loadBar(jentry, nentries);
        //
        // Determine the collision system of the sample
        std::string evtCol = "";
        if (isData) {
          if (muonTree->Event_Run() >= 285410 && muonTree->Event_Run() <= 285951) evtCol = "Pbp"; // for Pbp
          if (muonTree->Event_Run() >= 285952 && muonTree->Event_Run() <= 286504) evtCol = "pPb"; // for pPb
        }
        else {
          if (file.find("_Pbp_")!=std::string::npos) evtCol = "Pbp"; // for Pbp
          if (file.find("_pPb_")!=std::string::npos) evtCol = "pPb"; // for pPb
        }
        if (evtCol=="") { std::cout << "[ERROR] Could not determine the collision system in the sample" << std::endl; return; }
        //
        // Event Based Information
        //
        // Apply Event Filters
        if (PA::passEventFilter(metTree)==false) continue; // PA Event Selection
        //
        // Check Trigger Fired
        if (muonTree->Event_Trig_Fired()[triggerIndex]==false) continue; // Trigger Event Fired
        //
        // Get the HF Weight
        double weightEvt = 1.0;
        if (!isData && (HFCorr!=NULL)) {
          if (applyHFCorr==1) {
            const double hf = (evtTree->hiHF()>=300. ? 290. : evtTree->hiHF());
            weightEvt = HFCorr->weight(hf, HFweight::HFside::both, false);
          }
          else if (applyHFCorr==2) {
            const double nTrks = (evtTree->hiNtracks()>=300. ? 290. : evtTree->hiNtracks());
            weightEvt = HFCorr->weight(nTrks, HFweight::HFside::track, false);
          }
        }
        //
        // Loop over dimuons
        //
        for (uint idimu = 0; idimu < muonTree->PF_DiMuon_Mom().size(); idimu++) {
          // Do some simple selection on the events being used
          const ushort iPFMu1 = muonTree->PF_DiMuon_Muon1_Idx()[idimu];
          const ushort iPFMu2 = muonTree->PF_DiMuon_Muon2_Idx()[idimu];
          //
          // Only consider Muons Matched to GEN in MC
          if (!isData && muonTree->PF_Muon_Gen_Idx()[iPFMu1]<0) continue;
          if (!isData && muonTree->PF_Muon_Gen_Idx()[iPFMu2]<0) continue;
          // Check that the mother is a Z boson in MC and both share the same
          if (!isData) {
            const auto mom_1 = muonTree->MuonMother(muonTree->PF_Muon_Gen_Idx()[iPFMu1]);
            const auto mom_2 = muonTree->MuonMother(muonTree->PF_Muon_Gen_Idx()[iPFMu2]);
            if ( (mom_1.pdg != 23) || (mom_2.pdg != 23) || (mom_1.idx != mom_2.idx) ) continue;
          }
          //
          double weight = weightEvt;
          // Apply the Boson pT correction
          if (!isData && applyBosonPTCorr) {
            // Determine the gen boson pT
            const auto momIdx = PA::getGenMom(muonTree->PF_Muon_Gen_Idx()[iPFMu1], "MC_ZToMuMu", muonTree);
            if (momIdx>=0) {
              const auto momPdg = std::abs(muonTree->Gen_Particle_PdgId()[momIdx]);
              if (momPdg!=23) { std::cout << "[ERROR] Mother pdg " << momPdg << " can not be used to correct the boson pT" << std::endl; return; }
              double bosonPT = muonTree->Gen_Particle_Mom()[momIdx].Pt();
              if (bosonPT<0.5) { bosonPT = 0.5; }
              weight *= ( 1.0 / ( ( -0.37 * std::pow(bosonPT, -0.37) ) + 1.19 ) );
            }
            else { std::cout << "[ERROR] Mother of muon was not found!" << std::endl; return; }
          }
          // Consider muons passing Tight ID and isolation
          if (PA::isTightIsolatedMuon(iPFMu1 , muonTree)==false) continue;
          if (PA::isTightIsolatedMuon(iPFMu2 , muonTree)==false) continue;
          // Double Muon Charge Cut
          if(muonTree->PF_DiMuon_Charge()[idimu] != 0) continue; // Require opposite sign
          //
          // Get the momentum of each muon
          const TLorentzVector lep1 = muonTree->PF_Muon_Mom()[iPFMu1];
          const TLorentzVector lep2 = muonTree->PF_Muon_Mom()[iPFMu2];
          //
          // Kinematic Cuts
          //
          // Get the momentum of the dimuon
          const TLorentzVector dilep = lep1 + lep2;
          if(dilep.M() < MASS_LOW || dilep.M() > MASS_HIGH) continue; // Require to be within the Z mass window
          if(yRange=="mid" && (dilep.Rapidity()<-1.6 || dilep.Rapidity()>1.6)) continue; // Require dimuon to be within the defined rapidity window
          if(yRange=="fwd" && abs(dilep.Rapidity())<1.6) continue; // Require dimuon to be within the defined rapidity window
          if(lep1.Pt()        < PT_CUT  || lep2.Pt()        < PT_CUT)  continue;
          if(fabs(lep1.Eta()) > ETA_CUT || fabs(lep2.Eta()) > ETA_CUT) continue;
          //
          TVector2 dilep2D = TVector2(); dilep2D.SetMagPhi(dilep.Pt(),dilep.Phi()); // pT vector of dimuon
          //
          // Find the pt bin index
          int ipt = -1;
          for(uint ibin = 0; ibin < nbins; ibin++) {
            if(dilep.Pt() > ptBins[ibin] && dilep.Pt() <= ptBins[ibin+1]) ipt = ibin;
          }
          if(ipt<0) continue;
          //
          // Loop over met type
          //
          for (const auto& met : metType) {
            //
            // Extract the MET
            TVector2 met2D;
            if (met=="PF_RAW"       ) met2D = metTree->PF_MET_NoShift_Mom();
            if (met=="PF_Type1"     ) met2D = metTree->Type1_MET_NoShift_Mom();
            if (met=="PF_NoHF_RAW"  ) met2D = metNoHFTree->PF_MET_NoShift_Mom();
            if (met=="PF_NoHF_Type1") met2D = metNoHFTree->Type1_MET_NoShift_Mom();
            if (met=="PF_RAW_JetEnDown") met2D = metTree->PF_MET_JetEnDown_Mom();
            if (met=="PF_RAW_JetEnUp") met2D = metTree->PF_MET_JetEnUp_Mom();
            if (met=="PF_RAW_JetResDown") met2D = metTree->PF_MET_JetResDown_Mom();
            if (met=="PF_RAW_JetResUp") met2D = metTree->PF_MET_JetResUp_Mom();
            if (met=="PF_RAW_MuonEnDown") met2D = metTree->PF_MET_MuonEnDown_Mom();
            if (met=="PF_RAW_MuonEnUp") met2D = metTree->PF_MET_MuonEnUp_Mom();
            //
            // Compute Recoil
            TVector2 U = -1.*(dilep2D + met2D); // Recoil Vector
            // u1 : Recoil Component Parallel to the Dimuon pT Vector
            double  u1 = ( U * dilep2D ) / dilep2D.Mod();
            // u2 : Recoil Component Perpendicular to the Dimuon pT Vector
            double  u2 = ( ( U.Px() * dilep2D.Py() ) - ( U.Py() * dilep2D.Px() ) ) / dilep2D.Mod();
            //
            for (const auto& col : COLL) {
              if (col==evtCol || col=="PA") {
                // Fill the histograms with the recoil (u1 or u2) in each pT bin
                hPFu1v.at(met).at(col).at(ipt)->Fill(u1, weight);
                hPFu2v.at(met).at(col).at(ipt)->Fill(u2, weight);
                // Fill the RooWorkspaces with the recoil
                //
                RooDataSet* u1DS = ((RooDataSet*)myws_u1.at(met).at(col).data(Form("dPF%s_%s", uparName.c_str(), dsLabel.c_str())));
                RooDataSet* u2DS = ((RooDataSet*)myws_u2.at(met).at(col).data(Form("dPF%s_%s", uprpName.c_str(), dsLabel.c_str())));
                //
                const auto& u1Set = *u1DS->get();
                const auto& u2Set = *u2DS->get();
                //
                ((RooRealVar*)u1Set.find("pt"))->setVal(dilep.Pt());
                ((RooRealVar*)u2Set.find("pt"))->setVal(dilep.Pt());
                ((RooRealVar*)u1Set.find(uparName.c_str()))->setVal(u1);
                ((RooRealVar*)u2Set.find(uprpName.c_str()))->setVal(u2);
                //
                u1DS->addFast(u1Set, weight);
                u2DS->addFast(u2Set, weight);
              }
            }
          }
        }
      }
    }
    // Save the workspaces
    for (const auto& met : metType) {
      for (const auto& col : COLL) {
        // Create output directory
        const std::string dsDir = Form("%s/DataSet/", CWD.c_str());
        const std::string outputDir = dsDir + dsLabel +"/";
        if (!existDir(outputDir)) { makeDir(outputDir); }
        // Create output file
        std::string label = "";
        if (!isData) {
          if (applyHFCorr==1  ) { label += "_HFCorr";  } else if (applyHFCorr==2) { label += "_NTrack"; }
          if (applyBosonPTCorr) { label += "_BosonPT"; }
        }
        const std::string outfname = (outputDir + Form("DATASET_%s_%s%s.root", met.c_str(), col.c_str(), label.c_str()));
        auto outfile = std::unique_ptr<TFile>(new TFile(outfname.c_str(), "RECREATE"));
        if (outfile!=NULL && outfile->IsOpen() && !outfile->IsZombie()) {
          myws_u1.at(met).at(col).Write("workspace_u1");
          myws_u2.at(met).at(col).Write("workspace_u2");
          for (uint ibin = 0; ibin < nbins; ibin++) {
            if (hPFu1v.at(met).at(col).at(ibin)) { hPFu1v.at(met).at(col).at(ibin)->Write(Form("hPFu1_%i", ibin)); }
            if (hPFu2v.at(met).at(col).at(ibin)) { hPFu2v.at(met).at(col).at(ibin)->Write(Form("hPFu2_%i", ibin)); }
          }
        }
        if (outfile!=NULL) { outfile->Close(); }
      }
    }
  }
 
  // ------- Arrays and graphs to store fit results -----------
  std::map< std::string , std::map< std::string , std::map< std::string , std::vector< std::vector< double > > > > > u1VarArr;
  std::map< std::string , std::map< std::string , std::map< std::string , std::vector< std::vector< double > > > > > u2VarArr;
 
  // Initialize Canvas
  setTDRStyle();
  auto c = std::unique_ptr<TCanvas>(new TCanvas("c", "c", 800, 800));

  for (const auto& met : metType) {
    for (const auto& col : COLL) {
      std::cout << "[INFO] Working with MET " << met << " and coll: " << col  << std::endl;
      // Create output directories
      std::string labelDir = "";
      if (!isData) {
        if (applyHFCorr==1  ) { labelDir += "HFCorr";  } else if (applyHFCorr==2) { labelDir += "NTrack"; } else { labelDir += "noHFCorr"; }
        if (applyBosonPTCorr) { labelDir += "_BosonPT"; }
      }
      if (labelDir!="") { labelDir += "/"; }
      const std::string outputDir = mainDir + dsLabel +"/"+ ("MET_"+met) +"/"+ col +"/" + labelDir.c_str() + (pfu1model==1?"singleGauss/":(pfu1model==2?"doubleGauss/":(pfu1model==3?"doubleGauss/":"BWGauss/")));
      if (!existDir(outputDir)) { makeDir(outputDir); }
      // Do fits on u1
      std::cout << "[INFO] Proceed to perform the fit for u1" << std::endl;
      performFit(
                 myws_u1.at(met).at(col),
                 hPFu1v.at(met).at(col),
                 ptBins,
                 uparName,
                 rapRange,
                 pfu1model,
                 Form("dPF%s_%s", uparName.c_str(), dsLabel.c_str()),
                 col,
                 isData,
                 *c,
                 u1VarArr[met][col],
                 outputDir
                 );
      std::cout << "[INFO] Proceed to perform the fit for u2" << std::endl;
      // Do fits on u2
      performFit(
                 myws_u2.at(met).at(col),
                 hPFu2v.at(met).at(col),
                 ptBins,
                 uprpName,
                 rapRange,
                 pfu2model,
                 Form("dPF%s_%s", uprpName.c_str(), dsLabel.c_str()),
                 col,
                 isData,
                 *c,
                 u2VarArr[met][col],
                 outputDir
                 );
    }
  }

  //--------------------------------------------------------------------------------------------------------------
  // Make plots
  //==============================================================================================================  

  for (const auto& met : metType) {
    for (const auto& col : COLL) {
      // Create output directories
      std::string labelDir = "";
      if (!isData) {
        if (applyHFCorr==1) { labelDir += "HFCorr";  } else if (applyHFCorr==2) { labelDir += "NTrack"; } else { labelDir += "noHFCorr"; }
        if (applyBosonPTCorr) { labelDir += "_BosonPT"; }
      }
      if (labelDir!="") { labelDir += "/"; }
      const std::string outputDir = mainDir + dsLabel +"/"+ ("MET_"+met) +"/"+ col +"/" + labelDir.c_str()+ (pfu1model==1?"singleGauss/":(pfu1model==2?"doubleGauss/":(pfu1model==3?"doubleGauss/":"BWGauss/"))) +"Fits/";
      gSystem->mkdir(TString(outputDir + uparName + "/" + "png/"), true);
      gSystem->mkdir(TString(outputDir + uparName + "/" + "pdf/"), true);
      gSystem->mkdir(TString(outputDir + uparName + "/" + "root/"), true);
      gSystem->mkdir(TString(outputDir + uprpName + "/" + "png/"), true);
      gSystem->mkdir(TString(outputDir + uprpName + "/" + "pdf/"), true);
      gSystem->mkdir(TString(outputDir + uprpName + "/" + "root/"), true);
      // Define output graphs
      std::map< std::string , std::unique_ptr<TGraphAsymmErrors> > u1Graph;
      std::map< std::string , std::unique_ptr<TGraphAsymmErrors> > u2Graph;
      // Plotting u1 vs. dilepton pT
      for (const auto& var : u1VarArr.at(met).at(col)) {
        //
        std::vector<double> xval , xeLo , xeHi , yval , yeLo , yeHi;
        for (uint ibin = 0; ibin < nbins; ibin++) {
          if ( (var.second[ibin][1]>0) || (var.second[ibin][2]>0) ) {
            yval.push_back(var.second[ibin][0]);
            yeLo.push_back(var.second[ibin][1]);
            yeHi.push_back(var.second[ibin][2]);
            xval.push_back(var.second[ibin][3]);
            xeLo.push_back(var.second[ibin][4]);
            xeHi.push_back(var.second[ibin][5]);
          }
        }
        //
        u1Graph[var.first] = std::unique_ptr<TGraphAsymmErrors>(new TGraphAsymmErrors(yval.size(), &xval[0], &yval[0], &xeLo[0], &xeHi[0], &yeLo[0], &yeHi[0]));
        u1Graph[var.first]->SetName(Form("grPF%s%s", uparName.c_str(), var.first.c_str()));
        u1Graph[var.first]->SetTitle("");
        u1Graph[var.first]->SetMarkerColor(kBlack);
        u1Graph[var.first]->SetMarkerStyle(kOpenCircle);
        u1Graph[var.first]->Draw();
        u1Graph[var.first]->GetXaxis()->SetTitle("q_{T}(ll) [GeV/c]");
        std::string varLbl;
        if ( (var.first.find("rsigma")!=std::string::npos) || (var.first.find("dmean")!=std::string::npos) || (var.first.find("frac")!=std::string::npos) ) {
          varLbl = Form("%s(u_{1})", formatText(var.first).c_str());
        }
        else { varLbl = Form("%s(u_{1}) [GeV/c]", formatText(var.first).c_str()); }
        u1Graph[var.first]->GetYaxis()->SetTitle(varLbl.c_str());
        u1Graph[var.first]->GetXaxis()->SetTitleSize(0.045);
        u1Graph[var.first]->GetXaxis()->SetLabelSize(0.040);
        u1Graph[var.first]->GetXaxis()->SetTitleFont(42);
        u1Graph[var.first]->GetXaxis()->SetTitleOffset(1.28);
        u1Graph[var.first]->GetXaxis()->SetRangeUser(0.0, 200.);
        u1Graph[var.first]->GetYaxis()->SetLabelSize(0.045);
        u1Graph[var.first]->GetYaxis()->SetTitleSize(0.040);
        u1Graph[var.first]->GetYaxis()->SetTitleOffset(1.7);
        u1Graph[var.first]->GetYaxis()->SetTitleFont(42);
        updateYAxisRange(u1Graph[var.first]);
        u1Graph[var.first]->Draw();
        c->SetLogx(0); c->Update();
        int lumiId = 0;
        if (isData ) { if (col=="pPb") { lumiId = 109; } else if (col=="Pbp") { lumiId = 110; } else if (col=="PA") { lumiId = 111; } }
        if (!isData) { if (col=="pPb") { lumiId = 112; } else if (col=="Pbp") { lumiId = 113; } else if (col=="PA") { lumiId = 114; } }
        std::unique_ptr<TPaveText> tb = std::unique_ptr<TPaveText>(new TPaveText(0.17,0.65,0.45,0.85,"NDC"));
        tb->SetTextColor(kBlack);
        tb->SetFillStyle(0);
        tb->SetBorderSize(0);
        tb->AddText(Form("y %s%s",((rapRange.find("<")!=std::string::npos)?"":"#in "),rapRange.c_str()));
        tb->Draw("same");
        CMS_lumi(c.get(), lumiId, 33, "");
        c->SaveAs( (outputDir + uparName + "/" + "png/"  + Form("pf%s%s.png" , uparName.c_str(), var.first.c_str())).c_str() );
        c->SaveAs( (outputDir + uparName + "/" + "pdf/"  + Form("pf%s%s.pdf" , uparName.c_str(), var.first.c_str())).c_str() );
        c->SaveAs( (outputDir + uparName + "/" + "root/" + Form("pf%s%s.root", uparName.c_str(), var.first.c_str())).c_str() );
        c->Clear();
      }
      // Plotting u2 vs. dilepton pT
      for (const auto& var : u2VarArr.at(met).at(col)) {
        //
        std::vector<double> xval , xeLo , xeHi , yval , yeLo , yeHi;
        for (uint ibin = 0; ibin < nbins; ibin++) {
          if ( (var.second[ibin][1]>0) || (var.second[ibin][2]>0) ) {
            yval.push_back(var.second[ibin][0]);
            yeLo.push_back(var.second[ibin][1]);
            yeHi.push_back(var.second[ibin][2]);
            xval.push_back(var.second[ibin][3]);
            xeLo.push_back(var.second[ibin][4]);
            xeHi.push_back(var.second[ibin][5]);
          }
        }
        //
        u2Graph[var.first] = std::unique_ptr<TGraphAsymmErrors>(new TGraphAsymmErrors(yval.size(), &xval[0], &yval[0], &xeLo[0], &xeHi[0], &yeLo[0], &yeHi[0]));
        u2Graph[var.first]->SetName(Form("grPF%s%s", uprpName.c_str(), var.first.c_str()));
        u2Graph[var.first]->SetTitle("");
        u2Graph[var.first]->Draw();
        u2Graph[var.first]->SetMarkerColor(kBlack);
        u2Graph[var.first]->SetMarkerStyle(kOpenCircle);
        u2Graph[var.first]->GetXaxis()->SetTitle("q_{T}(ll) [GeV/c]");
        std::string varLbl;
        if ( (var.first.find("rsigma")!=std::string::npos) || (var.first.find("dmean")!=std::string::npos) || (var.first.find("frac")!=std::string::npos) ) {
          varLbl = Form("%s(u_{2})", formatText(var.first).c_str());
        }
        else { varLbl = Form("%s(u_{2}) [GeV/c]", formatText(var.first).c_str()); }
        u2Graph[var.first]->GetYaxis()->SetTitle(varLbl.c_str());
        u2Graph[var.first]->GetXaxis()->SetTitleSize(0.045);
        u2Graph[var.first]->GetXaxis()->SetLabelSize(0.040);
        u2Graph[var.first]->GetXaxis()->SetTitleFont(42);
        u2Graph[var.first]->GetXaxis()->SetTitleOffset(1.28);
        u2Graph[var.first]->GetXaxis()->SetRangeUser(0.0, 200.);
        u2Graph[var.first]->GetYaxis()->SetLabelSize(0.045);
        u2Graph[var.first]->GetYaxis()->SetTitleSize(0.040);
        u2Graph[var.first]->GetYaxis()->SetTitleOffset(1.7);
        u2Graph[var.first]->GetYaxis()->SetTitleFont(42);
        u2Graph[var.first]->Draw();
        updateYAxisRange(u2Graph[var.first]);
        c->SetLogx(0); c->Update();
        if (var.first.find("mean")!=std::string::npos) {
          u2Graph[var.first]->GetYaxis()->SetRangeUser(-5.0, 5.0);
        }
        int lumiId = 0;
        if (isData ) { if (col=="pPb") { lumiId = 109; } else if (col=="Pbp") { lumiId = 110; } else if (col=="PA") { lumiId = 111; } }
        if (!isData) { if (col=="pPb") { lumiId = 112; } else if (col=="Pbp") { lumiId = 113; } else if (col=="PA") { lumiId = 114; } }
        CMS_lumi(c.get(), lumiId, 33, "");
        std::unique_ptr<TPaveText> tb = std::unique_ptr<TPaveText>(new TPaveText(0.17,0.65,0.45,0.85,"NDC"));
        tb->SetTextColor(kBlack);
        tb->SetFillStyle(0);
        tb->SetBorderSize(0);
        tb->AddText(Form("y %s%s",((rapRange.find("<")!=std::string::npos)?"":"#in "),rapRange.c_str()));
        tb->Draw("same");
        c->SaveAs( (outputDir + uprpName + "/" + "png/"  + Form("pf%s%s.png" , uprpName.c_str(), var.first.c_str())).c_str() );
        c->SaveAs( (outputDir + uprpName + "/" + "pdf/"  + Form("pf%s%s.pdf" , uprpName.c_str(), var.first.c_str())).c_str() );
        c->SaveAs( (outputDir + uprpName + "/" + "root/" + Form("pf%s%s.root", uprpName.c_str(), var.first.c_str())).c_str() );
        c->Clear();
      }
      // clean up

      //--------------------------------------------------------------------------------------------------------------
      // Output
      //==============================================================================================================

      const std::string outfname = (outputDir + Form("plots_RecoilPDF_%s_%s.root", met.c_str(), col.c_str()));
      auto outfile = std::unique_ptr<TFile>(new TFile(outfname.c_str(), "RECREATE"));
      for (const auto& graph : u1Graph) { if (graph.second) graph.second->Write(); }
      for (const auto& graph : u2Graph) { if (graph.second) graph.second->Write(); }
      outfile->Close();

      const std::string htmlDir = mainDir + dsLabel +"/"+ ("MET_"+met) +"/"+ col + "/" + labelDir.c_str() + (pfu1model==1?"singleGauss/":(pfu1model==2?"doubleGauss/":(pfu1model==3?"doubleGauss/":"BWGauss/")));
      makeHTML(htmlDir, u1Graph, u2Graph, uparName, uprpName, nbins);
  
      cout << "  <> Output saved in " << outputDir << endl;
    }
  }

  // clean up 
  c->Close();
}
  
//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
bool performFit(
                RooWorkspace& ws,
                const vector< std::unique_ptr<TH1D> >& hv,
                const std::vector< double >& ptBins,
                const std::string uName,
                const std::string rapRange,
                const uint model,
                const std::string dsName,
                const std::string col,
                const bool isData,
                TCanvas& c,
                std::map< std::string , std::vector< std::vector< double > > >& varArr,
                const std::string outputDir
                ) 
{
  //
  // Check input settings
  //
  if (hv.size()==0) { std::cout << "[ERROR] hv vector is empty"; return false; }
  //if (hv[0]->GetEntries()==0) { std::cout << "[ERROR] hv[0] has zero entries"; return false; }
  if ((ptBins.size()-1)!=hv.size()) { std::cout << "[ERROR] hv and ptBins don't have the same number of entries"; return false; }
  if (model>4 || model<1) { std::cout << "[ERROR] Wrong input model!"; return false; }
  if (ws.var(uName.c_str())==NULL) {  std::cout << "[ERROR] Recoil parameter was not found in the workspace!"; return false; }
  //
  // Create output directories
  //
  const std::string outDir = outputDir + "Plots/" + uName + "/";
  gSystem->mkdir( (outDir + "png/").c_str(), true);
  gSystem->mkdir( (outDir + "pdf/").c_str(), true);
  //
  // Set up fit parameters
  //
  std::cout << "[INFO] Creating the models and parameters for the fit" << std::endl;
  //
  // Width Value
  if (model>=1) { ws.factory( Form("sigma1[%.6f, %.6f, %.6f]" , hv[0]->GetRMS(), 0.001*(hv[0]->GetRMS()), 1000.0*(hv[0]->GetRMS())) ); }
  if (model>=2) {
    ws.factory( Form("rsigma2[%.6f, %.6f, %.6f]" , isData?2.0:0.55, isData?0.1:-0.5, isData?5.0:1.0) );
    ws.factory( "RooFormulaVar::sigma2( '@0*@1' , {sigma1, rsigma2})" );
  }
  if (model==3) {
    ws.factory( Form("rsigma3[%.6f, %.6f, %.6f]" , 2.0, 0.1, 5.0) );
    ws.factory( "RooFormulaVar::sigma3( '@0*@1' , {sigma2, rsigma3})" );
  }
  //
  // Mean Value
  if (model>=1) { ws.factory( Form("mean1[%.6f, %.6f, %.6f]" , hv[0]->GetMean(), hv[0]->GetXaxis()->GetXmin(), hv[0]->GetXaxis()->GetXmax()) ); }
  if (model>=2) {
    ws.factory( Form("dmean2[%.6f, %.6f, %.6f]", 0.0, -6.0, 6.0) );
    if (isData || uName=="u2") ws.var("dmean2")->setConstant(true);
    ws.factory( "RooFormulaVar::mean2( '@0 + @1*@2' , {mean1, sigma1, dmean2})" );
  }
  if (model==3) {
    ws.factory( Form("dmean3[%.6f, %.6f, %.6f]", 0.0, -0.25, 0.4) );
    ws.var("dmean3")->setConstant(true);
    ws.factory( "RooFormulaVar::mean3( '@0 + @1*@2' , {mean2, sigma2, dmean3})" );
  }
  //
  // Fractions
  if (model==2) { ws.factory( Form("frac[%.6f, %.6f, %.6f]" , isData?0.70:0.45, 0.0, 1.0) ); ws.var("frac")->setConstant(true); }
  if (model==3) { ws.factory( Form("frac2[%.6f, %.6f, %.6f]", 0.70, 0.0, 1.0) ); ws.var("frac2")->setConstant(true);}
  if (model==4) { ws.factory( Form("frac[%.6f, %.6f, %.6f]" , isData?0.75:0.85, 0.0, 1.0) ); ws.var("frac")->setConstant(true); }

  //
  // Signal Model Functions
  ws.factory( Form("nsig[%.6f, %.6f, %.6f]" , hv[0]->Integral(), 0., 2.0*(hv[0]->Integral())) );
  if (model>=1) { ws.factory( Form("Gaussian::gauss1(%s, mean1, sigma1)", uName.c_str()) ); }
  if (model==2) { ws.factory( Form("Gaussian::gauss2(%s, mean2, sigma2)", uName.c_str()) ); }
  if (model==3) { ws.factory( Form("Gaussian::gauss3(%s, mean3, sigma3)", uName.c_str()) ); }
  if (model==4) { ws.factory( Form("RooBreitWigner::BW(%s, mean2, sigma2)", uName.c_str()) ); }
  //
  // Signal Model using recursive sum
  if (model==1) { ws.factory( "SUM::sig( gauss1 )" ); }
  if (model==2) { ws.factory( "RSUM::sig( frac*gauss1 , gauss2 )" ); }
  if (model==3) { ws.factory( "RSUM::sig( frac*gauss1 , frac2*gauss2 , gauss3 )" ); }
  if (model==4) { ws.factory( "RSUM::sig( frac*gauss1 , BW )" ); }
  //
  // Define formula for overall mean (mean)
  if (model==2 || model==4) { ws.factory( "RooFormulaVar::mean(   '((@0)*@1 + (1.0-@0)*@2)' , {frac , mean1, mean2 } )" ); }
  if (model==3) { ws.factory( "RooFormulaVar::mean23( '((@0)*@1 + (1.0-@0)*@2)' , {frac2, mean2, mean3 } )" ); }
  if (model==3) { ws.factory( "RooFormulaVar::mean(   '((@0)*@1 + (1.0-@0)*@2)' , {frac , mean1, mean23} )" ); }
  //
  // Define formula for overall width (sigma)
  //
  //    Average sigma
  if (model==2 || model==4) { ws.factory( "RooFormulaVar::sigma(   '((@0)*@1 + (1.0-@0)*@2)' , {frac , sigma1, sigma2 } )" ); }
  if (model==3) { ws.factory( "RooFormulaVar::sigma23( '((@0)*@1 + (1.0-@0)*@2)' , {frac2, sigma2, sigma3 } )" ); }
  if (model==3) { ws.factory( "RooFormulaVar::sigma(   '((@0)*@1 + (1.0-@0)*@2)' , {frac , sigma1, sigma23} )" ); }
  //
  //    Approximate sigma
  if (model>=2) { ws.factory( "RooFormulaVar::varR1( '@0*@0 + @1*@1 - @2*@2' , {sigma1, mean1, mean} )" ); }
  if (model>=2) { ws.factory( "RooFormulaVar::varR2( '@0*@0 + @1*@1 - @2*@2' , {sigma2, mean2, mean} )" ); }
  if (model==3) { ws.factory( "RooFormulaVar::varR3( '@0*@0 + @1*@1 - @2*@2' , {sigma3, mean3, mean} )" ); }
  if (model==3) { ws.factory( "RooFormulaVar::varR23( '( (@0)*@1 + (1.0-@0)*@2 )' , {frac2, varR2, varR3} )" ); }
  if (model==2 || model==4) { ws.factory( "RooFormulaVar::sigmaG( 'TMath::Sqrt( (@0)*@1 + (1.0-@0)*@2 )' , {frac, varR1, varR2} )" ); }
  if (model==3) { ws.factory( "RooFormulaVar::sigmaG( 'TMath::Sqrt( (@0)*@1 + (1.0-@0)*@2 )' , {frac, varR1, varR23} )" ); }
  //
  // Construct Fit Model
  RooAddPdf  modelpdf("model", "model", RooArgList(*ws.pdf("sig")), RooArgList(*ws.var("nsig")));
  ws.import(modelpdf);
  //
  // Initialize the output var arrays
  varArr.clear();
  RooArgSet listVar = ws.allVars();
  auto parIt = std::unique_ptr<TIterator>(listVar.createIterator());
  for (RooRealVar* it = (RooRealVar*)parIt->Next(); it!=NULL; it = (RooRealVar*)parIt->Next() ) {
    if (it->isConstant()) continue;
    std::string name = it->GetName(); if (name==uName || name=="pt" || name=="nsig" || name=="chi2" || name=="ndof") continue;
    for (uint ibin = 0; ibin < (ptBins.size()-1); ibin++) { varArr[name].push_back( { 0.0 , -1.0 , 0.0 , -1.0 } ); }
  }
  RooArgSet listFunc = ws.allFunctions();
  auto parFunIt = std::unique_ptr<TIterator>(listFunc.createIterator());
  for (RooRealVar* it = (RooRealVar*)parFunIt->Next(); it!=NULL; it = (RooRealVar*)parFunIt->Next() ) {
    if (std::string(it->GetName()).find("recursive")!=std::string::npos) continue;
    if (std::string(it->GetName()).find("mean")==std::string::npos && std::string(it->GetName()).find("sigma")==std::string::npos) continue;
    for (uint ibin = 0; ibin < (ptBins.size()-1); ibin++) { varArr[it->GetName()].push_back( { 0.0 , -1.0 , 0.0 , -1.0 } ); }
  }
  //
  // Loop over each bin
  //
  for(uint ibin = 0; ibin < (ptBins.size()-1); ibin++) {
    if (hv[ibin]->GetEntries()==0) { std::cout << "[WARNING] pT bin : " << ibin << " has zero entries!" << std::endl;  continue; }
    //
    // Define the dataset to fit
    //
    string strCut = Form("(%.8f < pt && pt <= %.8f)", ptBins[ibin], ptBins[ibin+1]);
    cout << "[INFO] Using local RooDataSet with cuts: " << strCut << endl;
    std::unique_ptr<RooDataSet> dataset = std::unique_ptr<RooDataSet>(((RooDataSet*)ws.data(dsName.c_str())->reduce(strCut.c_str())));
    if (dataset->sumEntries()==0){ std::cout << "[ERROR] No events passed the cuts!" << endl; return false; }
    //
    // Define the fit range
    //
    ws.var(uName.c_str())->setRange("FitWindow", hv[ibin]->GetXaxis()->GetXmin(), hv[ibin]->GetXaxis()->GetXmax());
    ws.var(uName.c_str())->setMin( hv[ibin]->GetXaxis()->GetXmin() );
    ws.var(uName.c_str())->setMax( hv[ibin]->GetXaxis()->GetXmax() );
    //
    // Update fit parameters
    //
    // Check if we can use the previous result
    bool usePrevResult = false;
    if ( ibin>0 && (varArr.count("mean1")>0) && (varArr.at("mean1").size()>(ibin-1)) ) { usePrevResult = true; }
    //
    // Mean Value
    ws.var("mean1")->setVal( hv[ibin]->GetMean() );
    ws.var("mean1")->setMin( hv[ibin]->GetXaxis()->GetXmin() );
    ws.var("mean1")->setMax( hv[ibin]->GetXaxis()->GetXmax() );
    //
    // Width Value
    if (usePrevResult) { ws.var("sigma1")->setVal( varArr.at("sigma1")[ibin-1][0] ); }
    //ws.var("sigma1")->setVal( hv[ibin]->GetRMS() );
    ws.var("sigma1")->setMin(0.1*(hv[ibin]->GetRMS()) );
    ws.var("sigma1")->setMax(10.0*(hv[ibin]->GetRMS()) );
    //
    // Yield Value
    ws.var("nsig")->setVal( hv[ibin]->Integral() );
    ws.var("nsig")->setMax( 2.0*(hv[ibin]->Integral()) );
    if (usePrevResult) {
      if (ws.var("dmean2") && varArr.count("dmean2")>0) { ws.var("dmean2")->setVal( -varArr.at("dmean2")[ibin-1][0] ); }
    }
    else {
      if (ws.var("dmean2")) { ws.var("dmean2")->setVal(0.0); }
      if (ws.var("dmean3")) { ws.var("dmean3")->setVal(0.0); }
    }
    if (ws.var("rsigma2")) { ws.var("rsigma2")->setVal(isData?2.0:0.60); }
    if (ws.var("rsigma3")) { ws.var("rsigma3")->setVal(2.0); }
    if (!isData && uName=="u1" && ( ptBins[ibin]>89.9) && model!=4) {ws.var("frac")->setVal(0.70); }
    //
    // Perform fit
    //
    auto fitResult = std::unique_ptr<RooFitResult>(ws.pdf("model")->fitTo(*dataset, RooFit::Extended(kTRUE), RooFit::Range("FitWindow"),
                                                                          RooFit::Minos(!dataset->isWeighted()), RooFit::NumCPU(32), RooFit::Save(), RooFit::PrintLevel(-1)));
    if (!fitResult || fitResult->status()!=0) { std::cout << "[ERROR] Fit failed for pt bin : " << ibin << std::endl; }
    if (!fitResult) { std::cout << "[ERROR] Fit results empty!" << std::endl; return false; }
    fitResult->Print("v");
    //
    // Store the fit results
    //
    const double xVal = dataset->mean(*ws.var("pt"));
    auto  rmsVar = std::unique_ptr<RooRealVar>(dataset->rmsVar(*ws.var("pt")));
    const double xErrLo = rmsVar->getVal();
    const double xErrHi = xErrLo;
    //
    std::map< std::string , bool > isAtLimit;
    std::cout << "[INFO] Extracting all variables from the workspace" << std::endl;
    RooArgSet listVar = ws.allVars(); // Needed to avoid segmentation fault
    auto parIt = std::unique_ptr<TIterator>(listVar.createIterator());
    for (RooRealVar* it = (RooRealVar*)parIt->Next(); it!=NULL; it = (RooRealVar*)parIt->Next() ) {
      if (it->isConstant()) continue;
      std::string name = it->GetName(); if (name==uName || name=="pt" || name=="nsig" || name=="chi2" || name=="ndof") continue;
      isAtLimit[name] = isParAtLimit(*ws.var(name.c_str()));
      double yVal = ws.var(name.c_str())->getVal(); if (name.find("mean")!=std::string::npos) { yVal *= -1.0; }
      const double yErrLo = getErrorLo(*ws.var(name.c_str()));
      const double yErrHi = getErrorHi(*ws.var(name.c_str()));
      varArr[name][ibin] = { yVal , yErrLo , yErrHi , xVal , xErrLo , xErrHi };
    }
    std::cout << "[INFO] Extracting all functions from the workspace" << std::endl;
    RooArgSet listFunc = ws.allFunctions(); // Needed to avoid segmentation fault
    auto parFunIt = std::unique_ptr<TIterator>(listFunc.createIterator());
    for (RooRealVar* it = (RooRealVar*)parFunIt->Next(); it!=NULL; it = (RooRealVar*)parFunIt->Next() ) {
      if (std::string(it->GetName()).find("recursive")!=std::string::npos) continue;
      if (std::string(it->GetName()).find("mean")==std::string::npos && std::string(it->GetName()).find("sigma")==std::string::npos) continue;
      double yVal = ws.function(it->GetName())->getVal(); if (std::string(it->GetName()).find("mean")!=std::string::npos) { yVal *= -1.0; }
      const double yErrLo = ws.function(it->GetName())->getPropagatedError(*fitResult);
      const double yErrHi = yErrLo;
      varArr[it->GetName()][ibin] = { yVal , yErrLo , yErrHi , xVal , xErrLo , xErrHi };
    }
    //
    // Plot the fit results
    //
    std::cout << "[INFO] Plotting the fit results" << std::endl;
    auto frame = std::unique_ptr<RooPlot>(ws.var(uName.c_str())->frame(RooFit::Bins(hv[ibin]->GetNbinsX())));
    dataset->plotOn(frame.get(), RooFit::MarkerStyle(kFullCircle), RooFit::MarkerSize(0.8), RooFit::DrawOption("ZP"));
    //ws.pdf("model")->plotOn(frame.get(), RooFit::VisualizeError(*fitResult, 1, 0), RooFit::FillStyle(3225), RooFit::LineColor(kOrange), RooFit::FillColor(kOrange),
    //                        RooFit::Normalization(dataset->sumEntries(), RooAbsReal::NumEvent), RooFit::Range("FitWindow"), RooFit::NormRange("FitWindow"));
    if(model>=2) ws.pdf("model")->plotOn(frame.get(), RooFit::Normalization(dataset->sumEntries(), RooAbsReal::NumEvent), RooFit::Range("FitWindow"), RooFit::NormRange("FitWindow"), 
                                         RooFit::Components("gauss1"), RooFit::LineStyle(kDashed), RooFit::LineColor(kRed));
    if(model==2) ws.pdf("model")->plotOn(frame.get(), RooFit::Normalization(dataset->sumEntries(), RooAbsReal::NumEvent), RooFit::Range("FitWindow"), RooFit::NormRange("FitWindow"),
                                         RooFit::Components("gauss2"), RooFit::LineStyle(kDashed), RooFit::LineColor(kCyan+2));
    if(model==3) ws.pdf("model")->plotOn(frame.get(), RooFit::Normalization(dataset->sumEntries(), RooAbsReal::NumEvent), RooFit::Range("FitWindow"), RooFit::NormRange("FitWindow"),
                                         RooFit::Components("gauss3"), RooFit::LineStyle(kDashed), RooFit::LineColor(kGreen+2));
    if(model==4) ws.pdf("model")->plotOn(frame.get(), RooFit::Normalization(dataset->sumEntries(), RooAbsReal::NumEvent), RooFit::Range("FitWindow"), RooFit::NormRange("FitWindow"),
                                         RooFit::Components("BW"), RooFit::LineStyle(kDashed), RooFit::LineColor(kGreen+2));
    ws.pdf("model")->plotOn(frame.get(), RooFit::Normalization(dataset->sumEntries(), RooAbsReal::NumEvent), RooFit::Range("FitWindow"), RooFit::NormRange("FitWindow"), RooFit::LineColor(kBlue));
    //
    // Extract the goodness of fit using Chi2
    auto frameTMP = std::unique_ptr<RooPlot>((RooPlot*)frame->Clone("TMP"));
    RooHist* hPull = new RooHist(2.0); // !!!DONT USE UNIQUE POINTER!!!!, 2 represents the binWidth
    if (!makePullHist(*hPull, *frameTMP, "", "", true)) { return false; }
    hPull->SetName("hPull");
    auto frame2 = std::unique_ptr<RooPlot>(ws.var(uName.c_str())->frame(RooFit::Title("Pull Distribution"), RooFit::Bins(hv[ibin]->GetNbinsX())));
    frame2->addPlotable(hPull, "EP");
    //
    // Create the text labels
    std::cout << "[INFO] Creating the text labels for the plots" << std::endl;
    std::string pname = "";
    const std::string xlabel   = Form("PF %c_{%c} [GeV]", uName[0], uName[1]);
    const std::string ylabel   = Form("Events / %.1f GeV/c", hv[ibin]->GetBinWidth(1));
    const std::string binlabel = Form("%.0f < q_{T} < %.0f ; y %s%s", ptBins[ibin], ptBins[ibin+1],((rapRange.find("<")!=std::string::npos)?"":"#in ") ,rapRange.c_str());
    const std::string nsigtext = Form("N_{evts} = %.0f #pm %.0f", ws.var("nsig")->getVal(), ws.var("nsig")->getError());
    std::map< std::string , std::string > varText;
    for (const auto& var : varArr) {
      if (isEqual(var.second[ibin][1], var.second[ibin][2], 2)) {
        varText[var.first] = Form("%s = %.2f #pm %.2f" , formatText(var.first).c_str(), var.second[ibin][0] , var.second[ibin][1]);
      }
      else {
        varText[var.first] = Form("%s = %.2f + %.2f - %.2f" , formatText(var.first).c_str(), var.second[ibin][0] , var.second[ibin][2] , var.second[ibin][1]);
      }
      if (isAtLimit[var.first]) { varText.at(var.first) += " (!)"; }
    }
    std::unique_ptr<TPaveText> tb;
    if (model==1) { tb = std::unique_ptr<TPaveText>(new TPaveText(0.17,0.65,0.45,0.85,"NDC")); }
    if (model==2) { tb = std::unique_ptr<TPaveText>(new TPaveText(0.17,0.50,0.45,0.85,"NDC")); }
    if (model==3) { tb = std::unique_ptr<TPaveText>(new TPaveText(0.17,0.35,0.45,0.85,"NDC")); }
    if (model==4) { tb = std::unique_ptr<TPaveText>(new TPaveText(0.17,0.20,0.45,0.85,"NDC")); }
    tb->SetTextColor(kBlack);
    tb->SetFillStyle(0);
    tb->SetBorderSize(0);
    tb->AddText(binlabel.c_str());
    tb->AddText(nsigtext.c_str());
    for (const auto& var : varText) { if (var.first.find("nsig" )!=std::string::npos) tb->AddText(var.second.c_str()); }
    if (varText.count("mean")>0 ) { tb->AddText(varText["mean"].c_str());  }
    if (varText.count("mean1")>0) { tb->AddText(varText["mean1"].c_str()); }
    for (const auto& var : varText) { if (var.first.find("dmean" )!=std::string::npos) tb->AddText(var.second.c_str()); }
    if (varText.count("sigma")>0 )  { tb->AddText(varText["sigma"].c_str());   }
    if (varText.count("sigmaG")>0)  { tb->AddText(varText["sigmaG"].c_str());  }
    if (varText.count("sigma1")>0 ) { tb->AddText(varText["sigma1"].c_str());  }
    for (const auto& var : varText) { if (var.first.find("rsigma")!=std::string::npos) tb->AddText(var.second.c_str()); }
    for (const auto& var : varText) { if (var.first.find("frac" )!=std::string::npos) tb->AddText(var.second.c_str()); }
    //
    // Format the frames
    std::cout << "[INFO] Format the frames" << std::endl;
    //  Main Frame
    frame->SetTitle("");
    frame->GetXaxis()->SetTitle("");
    frame->GetXaxis()->SetTitleSize(0.05);
    frame->GetXaxis()->SetTitleFont(42);
    frame->GetXaxis()->SetTitleOffset(3);
    frame->GetXaxis()->SetLabelOffset(3);
    //frame->GetXaxis()->SetTitle(ylabel.c_str());
    frame->GetYaxis()->SetLabelSize(0.044);
    frame->GetYaxis()->SetTitleSize(0.044);
    frame->GetYaxis()->SetTitleOffset(1.7);
    frame->GetYaxis()->SetTitleFont(42);
    //  Pull Frame
    frame2->SetTitle("");
    frame2->GetYaxis()->CenterTitle(kTRUE);
    frame2->GetYaxis()->SetTitleOffset(0.4);
    frame2->GetYaxis()->SetTitleSize(0.16);
    frame2->GetYaxis()->SetLabelSize(0.1);
    frame2->GetYaxis()->SetTitle("Pull");
    frame2->GetXaxis()->SetTitleOffset(1);
    frame2->GetXaxis()->SetTitleSize(0.16);
    frame2->GetXaxis()->SetLabelSize(0.14);
    frame2->GetXaxis()->SetTitle(xlabel.c_str());
    frame2->GetYaxis()->SetRangeUser(-7.0, 7.0);
    // Define the plotting pads
    TPad *pad1  = new TPad("pad1", "", 0, 0.23, 1, 1);  // Unique Pointer does produce Segmentation Fault, so don't use it
    TPad *pad2  = new TPad("pad2", "", 0, 0, 1, 0.228); // Unique Pointer does produce Segmentation Fault, so don't use it
    auto  pline = std::unique_ptr<TLine>(new TLine(ws.var(uName.c_str())->getMin(), 0.0,  ws.var(uName.c_str())->getMax(), 0.0));
    // Format the pads
    c.cd();
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.4);
    pad2->SetFillStyle(4000); 
    pad2->SetFrameFillStyle(4000);
    pad1->SetBottomMargin(0.015);
    //
    // Draw the frames
    std::cout << "[INFO] Drawing the frames" << std::endl;
    //
    // Main Frame
    pad1->Draw();
    pad1->cd(); 
    frame->Draw();
    tb->Draw("same");
    //
    // Apply CMS style to pad
    std::cout << "[INFO] Setting the CMS style on the plot" << std::endl;
    int lumiId = 0;
    if (isData ) { if (col=="pPb") { lumiId = 109; } else if (col=="Pbp") { lumiId = 110; } else if (col=="PA") { lumiId = 111; } }
    if (!isData) { if (col=="pPb") { lumiId = 112; } else if (col=="Pbp") { lumiId = 113; } else if (col=="PA") { lumiId = 114; } }
    CMS_lumi(pad1, lumiId, 33, "");
    gStyle->SetTitleFontSize(0.05);
    pad1->Update();
    //
    // Pull Frame
    c.cd();
    pad2->Draw();
    pad2->cd();
    frame2->Draw();
    //
    if(!printChi2(*pad2, ws, *frameTMP, uName, *dataset, "model")) { return false; };
    pline->Draw("same");
    pad2->Update();
    //
    // Save the pads
    //
    std::cout << "[INFO] Saving the plots with linear scale" << std::endl;
    updateYAxisRange(*frame, ws, uName, *dataset, false);
    pad1->SetLogy(false);
    pad1->SetLogx(false);
    pad1->Update();
    pname = Form("pf%sfit_%i", uName.c_str(), ibin);
    c.SaveAs((outDir + "png/" + pname + ".png").c_str());
    c.SaveAs((outDir + "pdf/" + pname + ".pdf").c_str());
    //
    std::cout << "[INFO] Saving the plots with logarithmic scale" << std::endl;
    updateYAxisRange(*frame, ws, uName, *dataset, true);
    pad1->SetLogy(true);
    pad1->SetLogx(false);
    pad1->Update();
    c.Update();
    pname = Form("pf%sfitlog_%i", uName.c_str(), ibin);
    c.SaveAs((outDir + "png/" + pname + ".png").c_str());
    c.SaveAs((outDir + "pdf/" + pname + ".pdf").c_str());
    //
    std::cout << "[INFO] All done with this plot" << std::endl;
    c.Clear();
  }
  // Clean up
  // return
  return true;
};


//--------------------------------------------------------------------------------------------------
void makeHTML(
              const std::string outDir,
              const std::map< std::string , std::unique_ptr<TGraphAsymmErrors> >& u1Graph,
              const std::map< std::string , std::unique_ptr<TGraphAsymmErrors> >& u2Graph,
              const std::string uparName   = "u1",
              const std::string uprpName   = "u2",
              const uint nbins = 10
              )
{

  ofstream htmlfile;
  std::string htmlfname;
  htmlfname = Form("%s/plots.html", outDir.c_str());
  std::cout << "[INFO] Proceed to create HLTM file: " << htmlfname << std::endl;
  htmlfile.open(htmlfname);
  if (!htmlfile.is_open()) { std::cout << "[ERROR] Was not able to create file: " << htmlfname << std::endl; return; }
  htmlfile << "<!DOCTYPE html" << endl;
  htmlfile << "    PUBLIC \"-//W3C//DTD HTML 3.2//EN\">" << endl;
  htmlfile << "<html>" << endl;
  htmlfile << "<head><title>Recoil Fits</title></head>" << endl;
  htmlfile << "<body bgcolor=\"EEEEEE\">" << endl;

  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl;
  uint counter = 0;
  for (const auto& gr : u1Graph) {
    if (counter%5 == 0) { htmlfile << "<tr>" << endl; }
    htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"Fits/" << uparName << "/png/pf" << uparName << gr.first << ".png\"><img src=\"Fits/" << uparName << "/png/pf" << uparName << gr.first << ".png\" alt=\"Fits/" << uparName << "/png/pf" << uparName << gr.first << ".png\" width=\"100%\"></a></td>" << endl;
    if (counter%5 == 4) { htmlfile << "</tr>" << endl; }
    counter++;
  }
  if (u1Graph.count("sigma2")==0) {
    htmlfile << "<td width=\"20%\"></td>" << endl;
    htmlfile << "<td width=\"20%\"></td>" << endl;
    htmlfile << "<td width=\"20%\"></td>" << endl;
  }
  htmlfile << "</tr>" << endl;
  
  htmlfile << "</table>" << endl;
  htmlfile << "PF " << uparName << " fits:";
  htmlfile << " <a target=\"_blank\" href=\"pf" << uparName << "fits.html\">[linear scale]</a>";
  htmlfile << " <a target=\"_blank\" href=\"pf" << uparName << "fitslog.html\">[log scale]</a>" << endl;
  htmlfile << "<hr />" << endl;

  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl;
  htmlfile << "<tr>" << endl;
  counter = 0;
  for (const auto& gr : u2Graph) {
    if (counter%5 == 0) { htmlfile << "<tr>" << endl; }
    htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"Fits/" << uprpName << "/png/pf" << uprpName << gr.first << ".png\"><img src=\"Fits/" << uprpName << "/png/pf" << uprpName << gr.first << ".png\" alt=\"Fits/" << uprpName << "/png/pf" << uprpName << gr.first << ".png\" width=\"100%\"></a></td>" << endl;
    if (counter%5 == 4) { htmlfile << "</tr>" << endl; }
    counter++;
  }
  if (u2Graph.count("sigma2")==0) {
    htmlfile << "<td width=\"20%\"></td>" << endl;
    htmlfile << "<td width=\"20%\"></td>" << endl;
    htmlfile << "<td width=\"20%\"></td>" << endl;
  }
  htmlfile << "</tr>" << endl;
  htmlfile << "</table>" << endl;
  htmlfile << "PF " << uprpName << " fits:";
  htmlfile << " <a target=\"_blank\" href=\"pf" << uprpName << "fits.html\">[linear scale]</a>";
  htmlfile << " <a target=\"_blank\" href=\"pf" << uprpName << "fitslog.html\">[log scale]</a>" << endl;
  htmlfile << "<hr />" << endl;
  
  htmlfile << "</body>" << endl;
  htmlfile << "</html>" << endl;
  htmlfile.close();
  
  uint ibin=0;
  
  //
  // PF u1 fits page
  //
  htmlfname = Form("%s/pf%sfits.html", outDir.c_str(), uparName.c_str());
  std::cout << "[INFO] Proceed to create HLTM file: " << htmlfname << std::endl;
  htmlfile.open(htmlfname);
  if (!htmlfile.is_open()) { std::cout << "[ERROR] Was not able to create file: " << htmlfname << std::endl; return; }
  htmlfile << "<!DOCTYPE html" << endl;
  htmlfile << "    PUBLIC \"-//W3C//DTD HTML 3.2//EN\">" << endl;
  htmlfile << "<html>" << endl;
  htmlfile << "<body bgcolor=\"EEEEEE\">" << endl;
  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl;   
  for(ibin=0; ibin<nbins; ibin++) {
    if(ibin%5 == 0) htmlfile << "<tr>" << endl;
    htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"Plots/" << uparName << "/png/pf" << uparName << "fit_" << ibin << ".png\"><img src=\"Plots/" << uparName << "/png/pf" << uparName << "fit_" << ibin << ".png\" alt=\"Plots/" << uparName << "/png/pf" << uparName << "fit_" << ibin << ".png\" width=\"100%\"></a></td>" << endl;
    if(ibin%5 == 4) htmlfile << "</tr>" << endl;
  }
  if(ibin%5 != 4) {
    while(ibin%5 != 4) {
      htmlfile << "<td width=\"20%\"></td>" << endl;
      ibin++;
    }
    htmlfile << "</tr>" << endl;
  }
  htmlfile << "</table>" << endl;
  htmlfile << "</body>" << endl;
  htmlfile << "</html>" << endl;
  htmlfile.close();  

  htmlfname = Form("%s/pf%sfitslog.html", outDir.c_str(), uparName.c_str());
  std::cout << "[INFO] Proceed to create HLTM file: " << htmlfname << std::endl;
  htmlfile.open(htmlfname);
  if (!htmlfile.is_open()) { std::cout << "[ERROR] Was not able to create file: " << htmlfname << std::endl; return; }
  htmlfile << "<!DOCTYPE html" << endl;
  htmlfile << "    PUBLIC \"-//W3C//DTD HTML 3.2//EN\">" << endl;
  htmlfile << "<html>" << endl;
  htmlfile << "<body bgcolor=\"EEEEEE\">" << endl;
  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl;   
  for(ibin=0; ibin<nbins; ibin++) {
    if(ibin%5 == 0) htmlfile << "<tr>" << endl;
    htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"Plots/" << uparName << "/png/pf" << uparName << "fitlog_" << ibin << ".png\"><img src=\"Plots/" << uparName << "/png/pf" << uparName << "fitlog_" << ibin << ".png\"alt=\"Plots/" << uparName << "/png/pf" << uparName << "fitlog_" << ibin << ".png\" width=\"100%\"></a></td>" << endl;
    if(ibin%5 == 4) htmlfile << "</tr>" << endl;
  }
  if(ibin%5 != 4) {
    while(ibin%5 != 4) {
      htmlfile << "<td width=\"20%\"></td>" << endl;
      ibin++;
    }
    htmlfile << "</tr>" << endl;
  }
  htmlfile << "</table>" << endl;
  htmlfile << "</body>" << endl;
  htmlfile << "</html>" << endl;
  htmlfile.close(); 
  
  //
  // PF " << uprpName << " fits page
  //
  htmlfname = Form("%s/pf%sfits.html", outDir.c_str(), uprpName.c_str());
  std::cout << "[INFO] Proceed to create HLTM file: " << htmlfname << std::endl;
  htmlfile.open(htmlfname);
  if (!htmlfile.is_open()) { std::cout << "[ERROR] Was not able to create file: " << htmlfname << std::endl; return; }
  htmlfile << "<!DOCTYPE html" << endl;
  htmlfile << "    PUBLIC \"-//W3C//DTD HTML 3.2//EN\">" << endl;
  htmlfile << "<html>" << endl;
  htmlfile << "<body bgcolor=\"EEEEEE\">" << endl;
  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl; 
  for(ibin=0; ibin<nbins; ibin++) {
    if(ibin%5 == 0) htmlfile << "<tr>" << endl;
    htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"Plots/" << uprpName << "/png/pf" << uprpName << "fit_" << ibin << ".png\"><img src=\"Plots/" << uprpName << "/png/pf" << uprpName << "fit_" << ibin << ".png\"alt=\"Plots/" << uprpName << "/png/pf" << uprpName << "fit_" << ibin << ".png\" width=\"100%\"></a></td>" << endl;
    if(ibin%5 == 4) htmlfile << "</tr>" << endl;
  }
  if(ibin%5 != 4) {
    while(ibin%5 != 4) {
      htmlfile << "<td width=\"20%\"></td>" << endl;
      ibin++;
    }
    htmlfile << "</tr>" << endl;
  }
  htmlfile << "</table>" << endl;
  htmlfile << "</body>" << endl;
  htmlfile << "</html>" << endl;
  htmlfile.close(); 

  htmlfname = Form("%s/pf%sfitslog.html", outDir.c_str(), uprpName.c_str());
  std::cout << "[INFO] Proceed to create HLTM file: " << htmlfname << std::endl;
  htmlfile.open(htmlfname);
  if (!htmlfile.is_open()) { std::cout << "[ERROR] Was not able to create file: " << htmlfname << std::endl; return; }
  htmlfile << "<!DOCTYPE html" << endl;
  htmlfile << "    PUBLIC \"-//W3C//DTD HTML 3.2//EN\">" << endl;
  htmlfile << "<html>" << endl;
  htmlfile << "<body bgcolor=\"EEEEEE\">" << endl;
  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl; 
  for(ibin=0; ibin<nbins; ibin++) {
    if(ibin%5 == 0) htmlfile << "<tr>" << endl;
    htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"Plots/" << uprpName << "/png/pf" << uprpName << "fitlog_" << ibin << ".png\"><img src=\"Plots/" << uprpName << "/png/pf" << uprpName << "fitlog_" << ibin << ".png\"alt=\"Plots/" << uprpName << "/png/pf" << uprpName << "fitlog_" << ibin << ".png\" width=\"100%\"></a></td>" << endl;
    if(ibin%5 == 4) htmlfile << "</tr>" << endl;
  }
  if(ibin%5 != 4) {
    while(ibin%5 != 4) {
      htmlfile << "<td width=\"20%\"></td>" << endl;
      ibin++;
    }
    htmlfile << "</tr>" << endl;
  }
  htmlfile << "</table>" << endl;
  htmlfile << "</body>" << endl;
  htmlfile << "</html>" << endl;
  htmlfile.close(); 
};


bool printChi2(TPad& pad, const RooWorkspace& ws, const RooPlot& frame, const string& varLabel, const RooDataSet& dataset, const string& pdfLabel)
{
  //
  if (ws.pdf(pdfLabel.c_str())==NULL) { std::cout << "[ERROR] PDF "     << pdfLabel  << " was not found!" << std::endl; return false; }
  //
  pad.cd();
  //
  TH1D hData, hFit;
  if(!rooPlotToTH1(hData, hFit, frame)) { return false; }
  //
  auto parList = std::unique_ptr<RooArgSet>(ws.pdf(pdfLabel.c_str())->getParameters(dataset));
  uint nFitPar = parList->selectByAttrib("Constant", kFALSE)->getSize();
  RooHist hPull(hData.GetBinWidth(1));
  if (!makePullHist(hPull, frame, "", "", true)) { return false; }
  double* ypulls = hPull.GetY();
  uint nFullBins = 0; double chi2=0;
  for (int i = 0; i < hData.GetNbinsX(); i++) {
    if ( (hData.GetBinError(i+1) != 0.0) && (hData.GetBinContent(i+1) != 0.0) ) {
      chi2 += ypulls[i]*ypulls[i];
      nFullBins++;
    }
  }
  const int ndof = (nFullBins - nFitPar);
  std::cout << "[INFO] Using Standard method gives Chi2/NDoF " << chi2 << " / " << ndof << std::endl;
  //
  TLatex t = TLatex(); t.SetNDC(); t.SetTextSize(0.12);
  t.DrawLatex(0.70, 0.85, Form("#chi^{2}/ndof = %.0f / %d ", chi2, ndof));
  //
  return true;
};


void updateYAxisRange(RooPlot& frame, const RooWorkspace& myws, std::string varLabel, const RooDataSet& dataset, const bool& logScale)
{
// Find maximum and minimum points of Plot to rescale Y axis
  auto h = std::unique_ptr<TH1>(dataset.createHistogram("hist", *myws.var(varLabel.c_str()), RooFit::Binning(frame.GetNbinsX(),frame.GetXaxis()->GetXmin(),frame.GetXaxis()->GetXmax())));
  Double_t yMax = h->GetBinContent(h->GetMaximumBin());
  Double_t yMin = 1e99; for (int i=1; i<=h->GetNbinsX(); i++) if (h->GetBinContent(i)>0) yMin = min(yMin, h->GetBinContent(i));
  yMax = yMax + sqrt(yMax);
  yMin = yMin - sqrt(yMin);
  // Set the up and down of y axis
  double fMax = 0.7 , fMin = 0.0;
  double yDown , yUp;
  if (logScale) {
    yDown = 0.1;
    yUp = yMax*100.;
  }
  else {
    yDown = std::floor((fMax*yMin - fMin*yMax)/(fMax - fMin));
    yUp   = std::ceil(yDown + (yMax-yMin)/(fMax-fMin));
  }
  // Update the y range
  frame.GetYaxis()->SetRangeUser(yDown, yUp);
};


void updateYAxisRange(std::unique_ptr<TGraphAsymmErrors>& graph)
{
  if (graph==NULL) return;
  // Find the max and min of graph
  double yMin = 9999999999. , yMax = -9999999999.;
  for (int i = 0; i < graph->GetN(); i++) {
    double x, y, yL, yH;
    graph->GetPoint(i, x, y);
    yL = y - graph->GetErrorY(i);
    yH = y + graph->GetErrorY(i);
    if (yH > yMax) { yMax = yH; }
    if (yL < yMin) { yMin = yL; }
  }
  // Set the up and down of y axis
  double fMax = 0.7 , fMin = 0.1;
  double yDown , yUp;
  yDown = (fMax*yMin - fMin*yMax)/(fMax - fMin);
  yUp   = yDown + (yMax-yMin)/(fMax-fMin);
  // Update the y range
  graph->GetYaxis()->SetRangeUser(yDown, yUp);
};


std::string formatText(const std::string& text)
{
  std::string out = text;
  if (out=="mean")  { out = "#minus#mu";   return out; }
  if (out=="sigma") { out = "#sigma"; return out; }
  if (out.find("mean")!=std::string::npos)  { out.replace(out.find("mean")  , std::string("mean").length()  , "#minus#mu_{"); out += "}"; }
  if (out.find("sigma")!=std::string::npos) { out.replace(out.find("sigma") , std::string("sigma").length() , "#sigma_{");    out += "}"; }
  if (out.find("frac")!=std::string::npos)  { out.replace(out.find("frac")  , std::string("frac").length()  , "f_{");         out += "}"; }
  return out;
};

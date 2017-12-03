#if !defined(__CINT__) || defined(__MAKECINT__)
// Auxiliary Headers
#include "../../Utilities/HiMuonTree.h"
#include "../../Utilities/EVENTUTILS.h"
#include "../../Results/Utilities/resultsUtils.h"
// ROOT headers
#include "TROOT.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
// c++ headers
#include <dirent.h>
#include <iostream>
#include <map>
#include <vector>
#include <tuple>
#include <string>

#endif


// ------------------ TYPE -------------------------------
using TH1DMap_t    =  std::map< std::string , std::map< std::string , std::map< std::string , std::map< std::string , std::map< std::string , TH1D > > > > >;
using BinMap_t     =  std::map< std::string , std::vector< double > >;
using BinMapMap_t  =  std::map< std::string , BinMap_t >;
using FileInfo_t   =  std::vector< std::pair< std::string , double > >;


// ------------------ FUNCTION -------------------------------
void     initHist1D  ( TH1DMap_t& h , const BinMapMap_t& binMap );
void     setInputVar ( std::map<std::string , VarBinMap>& inputVar, const TH1DMap_t& h , const BinMapMap_t& binMap );


// ------------------ GLOBAL ------------------------------- 
//
//
// Collision System
const std::vector< std::string > COLL_ = { "pPb" , "Pbp" , "PA" };
//
// Muon Charge
const std::vector< std::string > CHG_  = { "Pl" , "Mi" };
//
const std::vector< std::string > VAR_  = { "FSR" };
const std::vector< std::string > SUB_  = { "PRE" , "POST" };
//
// Input Files for analysis
const std::string path_MC = "root://cms-xrd-global.cern.ch//store/group/phys_heavyions/anstahll/EWQAnalysis2017/pPb2016/8160GeV/MC/Embedded/Official";
const std::map< std::string , std::vector< std::pair< std::string , double > > > inputFileMap_ = {
  {"MC_WToMuNu_Plus_pPb"      , { { Form("%s/%s", path_MC.c_str(), "POWHEG/HiEWQForest_Embedded_Official_POWHEG_CT14_EPPS16_WToMuNu_Plus_pPb_8160GeV_20170813.root")  , POWHEG::XSec.at("WToMuNu_Plus").at("pPb")  } } },
  {"MC_WToMuNu_Minus_pPb"     , { { Form("%s/%s", path_MC.c_str(), "POWHEG/HiEWQForest_Embedded_Official_POWHEG_CT14_EPPS16_WToMuNu_Minus_pPb_8160GeV_20170813.root") , POWHEG::XSec.at("WToMuNu_Minus").at("pPb") } } },
  {"MC_WToMuNu_Plus_Pbp"      , { { Form("%s/%s", path_MC.c_str(), "POWHEG/HiEWQForest_Embedded_Official_POWHEG_CT14_EPPS16_WToMuNu_Plus_Pbp_8160GeV_20170813.root")  , POWHEG::XSec.at("WToMuNu_Plus").at("Pbp")  } } },
  {"MC_WToMuNu_Minus_Pbp"     , { { Form("%s/%s", path_MC.c_str(), "POWHEG/HiEWQForest_Embedded_Official_POWHEG_CT14_EPPS16_WToMuNu_Minus_Pbp_8160GeV_20170813.root") , POWHEG::XSec.at("WToMuNu_Minus").at("Pbp") } } }
};
std::map< std::string , std::vector< std::string > > sampleType_;


void correctFSR(const bool useEtaCM = true)
{
  //
  // Initialize the Kinematic Bin info
  BinMapMap_t  MU_BIN;
  if (useEtaCM == false) {
    const BinMap_t  TMP = {
      { "Eta" , { -2.4 , -2.2 , -2.0 , -1.8 , -1.6 , -1.4 , -1.2 , -1.0 , -0.8 , -0.6 , -0.4 , -0.2 , 0.0 , 0.2 , 0.4 , 0.6 , 0.8 , 1.0 , 1.2 , 1.4 , 1.6 , 1.8 , 2.0 , 2.2 , 2.4 } }
    };
    for (const auto& v : TMP) { MU_BIN["PA"][v.first] = v.second; MU_BIN["pPb"][v.first] = v.second; MU_BIN["Pbp"][v.first] = v.second; }
  }
  else {
    const BinMap_t  TMP_pPb = {
      { "EtaCM" , { -2.86 , -2.60 , -2.40 , -2.20 , -2.00 , -1.80 , -1.60 , -1.40 , -1.20 , -1.00 , -0.80 , -0.60 , -0.40 , -0.20 , 0.00 , 0.20 , 0.40 , 0.60 , 0.80 , 1.00 , 1.20 , 1.40 , 1.60 , 1.80, 1.93 } }
    };
    const BinMap_t  TMP_Pbp = {
      { "EtaCM" , { -1.93 , -1.80 , -1.60 , -1.40 , -1.20 , -1.00 , -0.80 , -0.60 , -0.40 , -0.20 , 0.00 , 0.20 , 0.40 , 0.60 , 0.80 , 1.00 , 1.20 , 1.40 , 1.60 , 1.80 , 2.00 , 2.20 , 2.40 , 2.60 , 2.86 } }
    };
    for (const auto& v : TMP_pPb) { MU_BIN["PA"][v.first] = v.second; MU_BIN["pPb"][v.first] = v.second; }
    for (const auto& v : TMP_Pbp) { MU_BIN["Pbp"][v.first] = v.second; }
  }
  //
  // Change the working directory
  const std::string CWD = getcwd(NULL, 0);
  const std::string mainDir = Form("%s/FSRCorrection/", CWD.c_str());
  gSystem->mkdir(mainDir.c_str(), kTRUE);
  gSystem->ChangeDirectory(mainDir.c_str());
  //
  // Create list of samples
  sampleType_["sample"] = {};
  for (const auto& f : inputFileMap_) { std::string tmp = f.first; tmp = tmp.substr(0, tmp.find_last_of("_"));
    if (std::find(sampleType_.at("sample").begin(), sampleType_.at("sample").end(), tmp)==sampleType_.at("sample").end()) { sampleType_.at("sample").push_back(tmp); }
  }
  //
  // ------------------------------------------------------------------------------------------------------------------------
  //
  // Declare the histograms for the efficiency
  TH1DMap_t h1D;   // Stores the total and passing histograms separately
  //
  // Initialize the efficiencies
  initHist1D(h1D , MU_BIN);
  //
  // ------------------------------------------------------------------------------------------------------------------------
  //
  // Extract all the samples
  std::map< std::string , Long64_t > nentries;
  std::map< std::string , std::unique_ptr< HiMuonTree > > muonTree;
  for (const auto & inputFile : inputFileMap_) {
    const std::string sample = inputFile.first;
    std::vector< std::string >  fileInfo;
    for (const auto& f : inputFile.second) { fileInfo.push_back(f.first); }
    for (const auto& f : fileInfo) { std::cout << "[INFO] Loading File: " << f << std::endl; }
    //
    muonTree[sample] = std::unique_ptr<HiMuonTree>(new HiMuonTree());
    if (!muonTree.at(sample)->GetTree(fileInfo)) return;
    nentries[sample] = muonTree.at(sample)->GetEntries();
  }
  //
  // ------------------------------------------------------------------------------------------------------------------------
  //
  // Loop over the samples
  //
  for (auto & inputFile : inputFileMap_) {
    const std::string sample = inputFile.first;
    // Loop over the events
    int treeIdx = -1;
    double crossSection = -1.0;
    std::cout << "[INFO] Starting to process " << nentries.at(sample) << " nentries" << std::endl;
    for (Long64_t jentry = 0; jentry < nentries.at(sample); jentry++) { //
      //
      // Get the entry in the trees
      if (muonTree.at(sample)->GetEntry(jentry)<0) { std::cout << "[ERROR] Muon Tree invalid entry!"  << std::endl; return; }
      //
      if (muonTree.at(sample)->Chain()->GetTreeNumber()!=treeIdx) {
        treeIdx = muonTree.at(sample)->Chain()->GetTreeNumber();
        std::cout << "[INFO] Processing Root File: " << inputFile.second[treeIdx].first << std::endl;
        // Get the MC Cross-Section
        crossSection = inputFile.second[treeIdx].second;
      }
      //
      loadBar(jentry, nentries.at(sample));
      // 
      // Determine the collision system of the sample
      std::string col = "";
      if (sample.find("Pbp")!=std::string::npos) col = "Pbp"; // for Pbp
      if (sample.find("pPb")!=std::string::npos) col = "pPb"; // for pPb
      if (col=="") { std::cout << "[ERROR] Could not determine the collision system in the sample" << std::endl; return; }
      //
      // Determine the type of sample : i.e. MC_DYToMuMu
      std::string sampleType = sample;
      sampleType = sampleType.substr(0, (sampleType.find(col)-1));
      //
      // Get the Lumi re-weight for MC (global weight)
      const double lumi = ( (col=="pPb") ? PA::LUMI::Data_pPb : PA::LUMI::Data_Pbp );
      const double mcWeight = ( ( crossSection * lumi ) / muonTree.at(sample)->GetTreeEntries() );
      //
      // Check Muon Conditions
      //
      // Loop over the generated final muons
      for (ushort iGenMu = 0; iGenMu < muonTree.at(sample)->Gen_Muon_Mom().size(); iGenMu++) {
        // Check the content of the MC
        const auto mom = muonTree.at(sample)->findMuonMother(iGenMu, 24, 1);
        if (mom.pdg!=24) continue;
        //
        // Find the index of the PRE FSR muon
        int iGenMu_PRE = -1;
        for (uint iDau = 0; iDau < muonTree.at(sample)->Gen_Particle_Daughter_Idx()[mom.idx].size(); iDau++) {
          const auto dauIdx = muonTree.at(sample)->Gen_Particle_Daughter_Idx()[mom.idx][iDau];
          if (std::abs(muonTree.at(sample)->Gen_Particle_PdgId()[dauIdx])==13) { iGenMu_PRE = dauIdx; break; }
        }
        if (iGenMu_PRE == -1) { std::cout << "[ERROR] The PRE FSR muon was not found" << std::endl; return; }
        //if (iGenMu_PRE == muonTree.at(sample)->Gen_Muon_Particle_Idx()[iGenMu]) { std::cout << "[ERROR] The PRE and POST FSR muons have the same index" << std::endl; return; }
        //
        // Extract Kinematic Info from PRE-FSR muon
        const double mu_PRE_FSR_Gen_Pt  = muonTree.at(sample)->Gen_Particle_Mom()[iGenMu_PRE].Pt();
        const double mu_PRE_FSR_Gen_Eta = muonTree.at(sample)->Gen_Particle_Mom()[iGenMu_PRE].Eta();
        const double mu_PRE_FSR_Gen_EtaCM = PA::EtaLABtoCM( mu_PRE_FSR_Gen_Eta , (col == "pPb") );
        //
        // Extract Kinematic Info from POST-FSR muon
        const double mu_POST_FSR_Gen_Pt  = muonTree.at(sample)->Gen_Muon_Mom()[iGenMu].Pt();
        const double mu_POST_FSR_Gen_Eta = muonTree.at(sample)->Gen_Muon_Mom()[iGenMu].Eta();
        const double mu_POST_FSR_Gen_EtaCM = PA::EtaLABtoCM( mu_POST_FSR_Gen_Eta , (col == "pPb") );
        //
        const std::string var = ( useEtaCM ? "EtaCM" : "Eta" );
        //
        // Determine the charge of the generated muon
        const short charge = muonTree.at(sample)->Gen_Muon_Charge()[iGenMu];
        const std::string chg = ( (charge < 0) ? "Mi" : "Pl" );
        //
        // Fill Histograms
        //
        for (const auto& coll : COLL_) {
          if (coll!="PA" && coll!=col) continue;
          if ( (mu_PRE_FSR_Gen_Pt > 25.) && (std::abs(mu_PRE_FSR_Gen_Eta) < 2.4) ) {
            double value = ( useEtaCM ? mu_PRE_FSR_Gen_EtaCM : mu_PRE_FSR_Gen_Eta );
            if (coll=="PA" && col=="Pbp") {
              if (useEtaCM) {
                value = PA::EtaCMtoLAB(value, false); // Switch from CM -> LAB, in Pbp system
                value = -value;                       // Invert in the LAB frame, so from Pbp -> pPb
                value = PA::EtaLABtoCM(value, true);  // Switch from LAB -> CM, in pPb system
              }
              else {
                value = -value;                       // Invert in the LAB frame, so from Pbp -> pPb
              }
            }
            h1D.at(var).at(coll).at(chg).at("FSR").at("PRE").Fill(value , mcWeight);
          }
          //
          if ( (mu_POST_FSR_Gen_Pt > 25.) && (std::abs(mu_POST_FSR_Gen_Eta) < 2.4) ) {
            double value = ( useEtaCM ? mu_POST_FSR_Gen_EtaCM : mu_POST_FSR_Gen_Eta );
            if (coll=="PA" && col=="Pbp") {
              if (useEtaCM) {
                value = PA::EtaCMtoLAB(value, false); // Switch from CM -> LAB, in Pbp system
                value = -value;                       // Invert in the LAB frame, so from Pbp -> pPb
                value = PA::EtaLABtoCM(value, true);  // Switch from LAB -> CM, in pPb system
              }
              else {
                value = -value;                       // Invert in the LAB frame, so from Pbp -> pPb
              }
            }
            h1D.at(var).at(coll).at(chg).at("FSR").at("POST").Fill(value , mcWeight);
          }
        }
        //
      }
    }
  }
  //
  // ------------------------------------------------------------------------------------------------------------------------
  //
  // Extract the input variables
  std::map< std::string , VarBinMap > inputVar;
  setInputVar(inputVar, h1D, MU_BIN);
  //
  BinSextaMap var;
  //
  // --------------------------------------------------------------------------------- //
  //
  for(const auto& v : inputVar) {
    BinPentaMap tmpVar;
    //
    // Compute the Charge Asymmetry
    //
    if (!computeChargeAsymmetry(tmpVar, v.second)) { return; }
    //
    // Compute the Forward-Backward ratio
    //
    if (!computeForwardBackwardRatio(tmpVar, v.second)) { return; }
    //
    // Compute the Cross-Section
    //
    if (!computeCrossSection(tmpVar, v.second)) { return; }
    //
    const std::string name = ( (v.first=="FSR_POST") ? "Nominal" : v.first );
    var[name].push_back(tmpVar);
    //
  }
  //
  // --------------------------------------------------------------------------------- //
  //
  // Initialize the Output Graphs
  GraphPentaMap graph;
  iniResultsGraph(graph, var);
  //
  // --------------------------------------------------------------------------------- //
  //
  // Fill the Output Graphs
  //
  if (!fillResultsGraph(graph, var)) { return; }
  //
  // --------------------------------------------------------------------------------- //
  //
  // Draw the Output Graphs
  const std::string outDir = CWD + "/Output";
  drawGraph(graph, outDir, useEtaCM, "", "");
  //
  // Draw combined systematic graph
  drawCombineSystematicGraph(graph, outDir, useEtaCM, "", "");
  //
};


void initHist1D(TH1DMap_t& h, const BinMapMap_t& binMap)
{
  //
  for (const auto& col : COLL_) {
    for (const auto& bins : binMap.at(col)) {
      Double_t bin[bins.second.size()];
      for (uint i = 0; i < bins.second.size(); i++) { bin[i] = bins.second[i]; }
      for (const auto& chg : CHG_) {
        for (const auto& type : VAR_) {
          for (const auto& sub : SUB_) {
            h[bins.first][col][chg][type][sub] = TH1D();
            h.at(bins.first).at(col).at(chg).at(type).at(sub) = TH1D(sub.c_str(), sub.c_str(), (bins.second.size()-1), bin);
            //
            const std::string name = "h1D_" + bins.first +"_"+ col +"_"+ chg +"_"+ type +"_"+ sub;
            h.at(bins.first).at(col).at(chg).at(type).at(sub).SetName(name.c_str());
            h.at(bins.first).at(col).at(chg).at(type).at(sub).Sumw2();
          }
        }
      }
    }
  }
};


void setInputVar(std::map<std::string , VarBinMap>& inputVar, const TH1DMap_t& h, const BinMapMap_t& binMap)
{
  for (const auto& col : COLL_) {
    for (const auto& bins : binMap.at(col)) {
      for (uint iBin = 1; iBin < bins.second.size(); iBin++) {
        for (const auto& ch : h.at(bins.first).at(col)) {
          for (const auto& t : ch.second) {
            for (const auto& sub : t.second) {
              //
              const std::string collSystem = col;
              const std::string charge     = ch.first;
              const anabin<0>   bin(bins.second[iBin-1], bins.second[iBin]);
              //
              const auto& hist = sub.second;
              // Compute the center value of the bin
              const double binVal = ( ( bins.second[iBin-1] + bins.second[iBin] ) / 2.0 );
              const int iHistBin = hist.GetXaxis()->FindBin(binVal);
              //
              const std::string name = Form("%s_%s", t.first.c_str(), sub.first.c_str());
              //
              inputVar[name][collSystem][bin][charge]["N_WToMu"]["Val"] = hist.GetBinContent(iHistBin);
              inputVar[name][collSystem][bin][charge]["N_WToMu"]["Err_Stat_High"] = hist.GetBinError(iHistBin);
              inputVar[name][collSystem][bin][charge]["N_WToMu"]["Err_Stat_Low" ] = hist.GetBinError(iHistBin);
              inputVar[name][collSystem][bin][charge]["N_WToMu"]["Err_Syst_High"] = 0.0;
              inputVar[name][collSystem][bin][charge]["N_WToMu"]["Err_Syst_Low" ] = 0.0;
              //
              if (collSystem == "pPb") { inputVar[name][collSystem][bin][charge]["Luminosity"]["Val"] = PA::LUMI::Data_pPb; }
              if (collSystem == "Pbp") { inputVar[name][collSystem][bin][charge]["Luminosity"]["Val"] = PA::LUMI::Data_Pbp; }
              if (collSystem == "PA" ) { inputVar[name][collSystem][bin][charge]["Luminosity"]["Val"] = ( PA::LUMI::Data_pPb + PA::LUMI::Data_Pbp ); }
            }
          }
        }
      }
    }
  }
};
              

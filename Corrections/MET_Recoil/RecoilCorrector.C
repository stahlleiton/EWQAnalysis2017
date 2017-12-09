#if !defined(__CINT__) || defined(__MAKECINT__)
// Auxiliary Headers
#include "RecoilCorrector.h"
#include "Math/RootFinder.h"
#include "Math/Functor.h"
#include "Math/QuantFuncMathCore.h"
// ROOT Headers
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>                    // file handle class
#include "TMath.h"                    // Mathematical tool class
#include <TFitResult.h>               // class to handle fit results
#include <TGraphErrors.h>             // graph class
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TKey.h>
#include <TClass.h>
// C++ Headers
#include <iostream>                   // standard I/O
#include <fstream>                    // standard I/O

#endif


//-----------------------------------------------------------------------------------------------------------------------------------------
RecoilCorrector::RecoilCorrector(const int seed)
{
  // Initialize all local parameters
  this->clear();
  // Initialize the random generator
  this->fRandom_.reset(new TRandom3(seed));
};


//-----------------------------------------------------------------------------------------------------------------------------------------
RecoilCorrector::~RecoilCorrector()
{
  // Clean up memory
  this->clear();
};


//-----------------------------------------------------------------------------------------------------------------------------------------
bool RecoilCorrector::isValid()
{
  if ( (u1Fits_MC_.size()==0) || (u2Fits_MC_.size()==0) || (u1Fits_DATA_.size()==0) || (u2Fits_DATA_.size()==0) ) { return false; }
  if ( (std::abs(reference_pT_.Mod())<0.000001) && (std::abs(reference_pT_.Phi())<0.000001) ) { return false; }
  if ( (std::abs(boson_pT_.Mod())<0.000001) && (std::abs(boson_pT_.Phi())<0.000001) ) { return false; }
  if ( fRandom_==NULL ) { return false; }
  return true;
};


//-----------------------------------------------------------------------------------------------------------------------------------------
bool RecoilCorrector::setInputFiles(
                                    const std::string& METType,
                                    const std::string& mcFileName,
                                    const std::string& dataFileName
                                    )
{
  //
  // Extract the fits from Z MC/DATA samples
  //
  if (!getRecoilFits(mcFileName   , u1Fits_MC_[METType]   , u2Fits_MC_[METType]  )) { return false; }
  if (!getRecoilFits(dataFileName , u1Fits_DATA_[METType] , u2Fits_DATA_[METType])) { return false; }
  this->METType_ =  METType;
  //
  // Return back
  //
  return true;
};

//-----------------------------------------------------------------------------------------------------------------------------------------
void RecoilCorrector::setInitialSetup(
                                      const std::string& sample,
                                      const bool useOneGaussian_MC,
                                      const bool useOneGaussian_DATA,
                                      const bool corrMean
                                      )
{
  //
  // Setup all initial settings
  //
  // Initialize all other variables
  this->useOneGaussian_MC_   = useOneGaussian_MC;
  this->useOneGaussian_DATA_ = useOneGaussian_DATA;
  this->corrMean_            = corrMean;
  if (sample.find("QCD")!=std::string::npos) { this->corrMean_ = false; }
};


//-----------------------------------------------------------------------------------------------------------------------------------------
void RecoilCorrector::clear()
{
  //
  // Clear all global variables
  //
  this->reference_pT_        = TVector2();
  this->boson_pT_            = TVector2();
  this->useOneGaussian_MC_   = false;
  this->useOneGaussian_DATA_ = false;
  this->corrMean_            = true;
  this->METType_             = "";
  //
  this->u1Fits_MC_.clear();
  this->u2Fits_MC_.clear();
  this->u1Fits_DATA_.clear();
  this->u2Fits_DATA_.clear();
};


//-----------------------------------------------------------------------------------------------------------------------------------------
bool RecoilCorrector::getRecoilFits(
                                    const std::string& fileName,
                                    std::map< std::string , std::unique_ptr<TF1> >& u1Fit,
                                    std::map< std::string , std::unique_ptr<TF1> >& u2Fit
                                   )
{
  // Initialize output variables
  u1Fit.clear(); u2Fit.clear();
  //
  // Open the input file
  //
  std::cout << "[INFO] Reading file: " << fileName << std::endl;
  auto file = std::unique_ptr<TFile>(TFile::Open(fileName.c_str(), "READ"));
  if (!file || file->IsZombie()) { std::cout << "[ERROR] Failed to open file: " << fileName << std::endl; return false; }
  //
  // Loop over all TF1 objects inside the file
  //
  TIter next(file->GetListOfKeys());
  for (TKey* key = (TKey*)next(); key!=NULL; key = (TKey*)next() ) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TF1")) continue;
    TF1 *f = (TF1*)key->ReadObj();
    std::string name = f->GetName();
    std::cout << "[INFO] Importing TF1: " << name << std::endl;
    if ( (name.find("PFu1")!=std::string::npos) ) { u1Fit[name.substr(name.find("u1")+2)].reset(f); }
    if ( (name.find("PFu2")!=std::string::npos) ) { u2Fit[name.substr(name.find("u2")+2)].reset(f); }
  }
  //
  // Clean up the memory
  //
  file->Close("R");
  //
  // Return back
  //
  return true;
};


//-----------------------------------------------------------------------------------------------------------------------------------------
void RecoilCorrector::evalRecoilFits(
                                     std::map< std::string , double >& uPar,
                                     const double& pT,
                                     const std::map< std::string , std::unique_ptr<TF1> >& uFit,
                                     const bool corrMean
                                     )
{
  // Initialize output variables
  uPar.clear();
  //
  // Evaluate the Recoil Fit functions
  //
  for (const auto& fit : uFit) {
    if (fit.first.find("mean")!=std::string::npos) {
      uPar[fit.first] = -1.0*fit.second->Eval(pT); // Mean values are defined with a negative in the functions
    }
    else { uPar[fit.first] = fit.second->Eval(pT); }
  }
  //
  // Compute the other parameters
  //
  const int model = findModel(uPar);
  //
  if (model==1) {
    uPar["mean"]   = uPar.at("mean1");
    uPar["sigmaG"] = uPar.at("sigma1");
    uPar["sigma"]  = uPar.at("sigma1");
  }    
  if (model==2) {
    double f = -1;
    if      (uPar.at("sigma1") != uPar.at("sigma2")) { f = ( (uPar.at("sigma") - uPar.at("sigma2")) / (uPar.at("sigma1") - uPar.at("sigma2")) ); }
    else if (uPar.at("mean1")  != uPar.at("mean2") ) { f = ( (uPar.at("mean")  - uPar.at("mean2") ) / (uPar.at("mean1")  - uPar.at("mean2") ) ); }
    if (f > 1.0) { f = 1.0; } else if (f < 0.0) { f = 0.0; }
    uPar["f"] = f;
    /*
    uPar["mean"] = uPar.at("f")*uPar.at("mean1") + (1.0 - uPar.at("f"))*uPar.at("mean2");
    double var1 = uPar.at("sigma1")*uPar.at("sigma1") + uPar.at("mean1")*uPar.at("mean1") - uPar.at("mean")*uPar.at("mean");
    double var2 = uPar.at("sigma2")*uPar.at("sigma2") + uPar.at("mean2")*uPar.at("mean2") - uPar.at("mean")*uPar.at("mean");
    uPar["sigmaG"] = TMath::Sqrt( uPar.at("f")*var1 + (1.0 - uPar.at("f"))*var2 );
    uPar["sigma"] = uPar.at("f")*uPar.at("sigma1") + (1.0 - uPar.at("f"))*uPar.at("sigma2");
    */
  }
  if (model==3) {
    double f2 = -1;
    if      (uPar.at("sigma2") != uPar.at("sigma3")) { f2 = ( (uPar.at("sigma23") - uPar.at("sigma3")) / (uPar.at("sigma2") - uPar.at("sigma3")) ); }
    else if (uPar.at("mean2")  != uPar.at("mean3") ) { f2 = ( (uPar.at("mean23")  - uPar.at("mean3") ) / (uPar.at("mean2")  - uPar.at("mean3") ) ); }
    if (f2 > 1.0) { f2 = 1.0; } else if (f2 < 0.0) { f2 = 0.0; }
    double f = -1;
    if      (uPar.at("sigma1") != uPar.at("sigma23")) { f = ( (uPar.at("sigma") - uPar.at("sigma23")) / (uPar.at("sigma1") - uPar.at("sigma23")) ); }
    else if (uPar.at("mean1")  != uPar.at("mean23") ) { f = ( (uPar.at("mean")  - uPar.at("mean23") ) / (uPar.at("mean1")  - uPar.at("mean23") ) ); }
    if (f > 1.0) { f = 1.0; } else if (f < 0.0) { f = 0.0; }
    uPar["f2"] = f2;
    uPar["f"]  = f;
    /*
    uPar["mean23"] = uPar.at("f2")*uPar.at("mean2") + (1.0 - uPar.at("f2"))*uPar.at("mean3");
    uPar["mean"]   = uPar.at("f")*uPar.at("mean1")  + (1.0 - uPar.at("f"))*uPar.at("mean23");
    uPar["sigma23"] = uPar.at("f2")*uPar.at("sigma2") + (1.0 - uPar.at("f2"))*uPar.at("sigma3");
    double var1 = uPar.at("sigma1")*uPar.at("sigma1") + uPar.at("mean1")*uPar.at("mean1") - uPar.at("mean")*uPar.at("mean");
    double var2 = uPar.at("sigma2")*uPar.at("sigma2") + uPar.at("mean2")*uPar.at("mean2") - uPar.at("mean")*uPar.at("mean");
    double var3 = uPar.at("sigma3")*uPar.at("sigma3") + uPar.at("mean3")*uPar.at("mean3") - uPar.at("mean")*uPar.at("mean");
    uPar["sigmaG"] = TMath::Sqrt( uPar.at("f")*var1 + (1.0 - uPar.at("f"))*( uPar.at("f2")*var2 + (1.0 - uPar.at("f2"))*var3 ) );
    uPar["sigma"] = uPar.at("f")*uPar.at("sigma1") + (1.0 - uPar.at("f"))*uPar.at("sigma23");
    */
  }
};


//-----------------------------------------------------------------------------------------------------------------------------------------
int RecoilCorrector::findModel( const std::map< std::string , double >& uPar )
{
  int model = 0;
  if ( (model == 0) && (uPar.count("sigma1")>0 && uPar.count("mean1")>0) ) { model = 1; }
  if ( (model == 1) && (uPar.count("sigma2")>0 && uPar.count("mean2")>0) ) { model = 2; }
  if ( (model == 2) && (uPar.count("sigma3")>0 && uPar.count("mean3")>0) ) { model = 3; }
  if ( (model >= 2) && (uPar.count("sigma")==0   || uPar.count("mean")==0  ) ) { std::cout << "[ERROR] Parameters missing for CDF using model " << model << std::endl; return -1; }
  if ( (model >= 3) && (uPar.count("sigma23")==0 || uPar.count("mean23")==0) ) { std::cout << "[ERROR] Parameters missing for CDF using model " << model << std::endl; return -1; }
  return model;
};


//-----------------------------------------------------------------------------------------------------------------------------------------
TF1* RecoilCorrector::getFunctionCDF(
                                     const std::map< std::string , double >& uPar,
                                     const double& CDF,
                                     const bool useOneGaussian
                                     )
{
  // Initialize output variables
  TF1* fcnCDF = NULL;
  //
  // Determine which CDF model to use
  //
  const int model = findModel(uPar);
  if (model<1) { return fcnCDF; }
  //
  // Build the CDF model
  //
  if (model==1 || useOneGaussian) {
    if (uPar.at("sigma") == 0.) { std::cout << "[ERROR] Sigma of model 1 is zero!" << std::endl; return NULL; }
    fcnCDF = new TF1("fcnCDF", "( ROOT::Math::gaussian_cdf(x, [0], [1]) ) - [2]", -1000., 1000.);
    fcnCDF->SetParameters(uPar.at("sigma") , uPar.at("mean") , CDF);
  }
  else {
    if (model==2) {
      if ( (uPar.at("sigma1") == 0.) || (uPar.at("sigma2") == 0.) ) { std::cout << "[ERROR] Sigmas of model 2 are zero!" << std::endl; return NULL; }
      std::string fcnText = "([0])*( ROOT::Math::gaussian_cdf(x, [1], [2]) ) + ";
      fcnText +=      "(1.0 - [0])*( ROOT::Math::gaussian_cdf(x, [3], [4]) ) - [5]";
      fcnCDF = new TF1("fcnCDF", fcnText.c_str(), -1000., 1000.);
      fcnCDF->SetParameters(uPar.at("f") , uPar.at("sigma1") , uPar.at("mean1") , uPar.at("sigma2") , uPar.at("mean2") , CDF);
    }
    if (model==3) {
      if ( (uPar.at("sigma1") == 0.) || (uPar.at("sigma2") == 0.) || (uPar.at("sigma3") == 0.) ) { std::cout << "[ERROR] Sigmas of model 3 are zero!" << std::endl; return NULL; }
      std::string fcnText = "([0])*( ROOT::Math::gaussian_cdf(x, [2], [3]) ) + ";
      fcnText +=      "(1.0 - [0])*( ([1])*( ROOT::Math::gaussian_cdf(x, [4], [5]) ) + ";
      fcnText +=              "(1.0 - [1])*( ROOT::Math::gaussian_cdf(x, [6], [7]) ) - [8]";
      fcnCDF = new TF1("fcnCDF", fcnText.c_str(), -1000., 1000.);
      fcnCDF->SetParameters(uPar.at("f") , uPar.at("f2") , uPar.at("sigma1") , uPar.at("mean1") , uPar.at("sigma2") , uPar.at("mean2") , uPar.at("sigma3") , uPar.at("mean3") , CDF);
    }
  }
  //
  // Return back
  //
  return fcnCDF;
};


//-----------------------------------------------------------------------------------------------------------------------------------------
bool RecoilCorrector::getPtFromTree(
                                    TVector2& reference_pT,
                                    TVector2& boson_pT,
                                    const uint& muonIdx,
                                    const std::unique_ptr<HiMuonTree>& muonTree,
                                    const std::string& sample
                                    )
{
  //
  if (sample=="") { std::cout << "[ERROR] Sample is not defined!" << std::endl; return false; }
  //
  // Initialize output variables
  reference_pT = TVector2(); boson_pT = TVector2();
  //
  // Determine the reference and the boson pT vectors
  //
  // Case: DrellYan -> Muon + Muon
  if (sample.find("MC_DY")!=std::string::npos) {
    // Find the number of RECO muons matched to GEN muons from Z/gamma boson decay
    std::vector< std::vector< uint > > idx;
    // Get the mom of the leading muon
    const int genMuIdx = muonTree->PF_Muon_Gen_Idx()[muonIdx];
    if ( genMuIdx > -1 ) {
      const auto& muMom = muonTree->MuonMother(genMuIdx);
      if      ( muMom.pdg == 23 ) { idx.push_back({ muonIdx , uint(genMuIdx) , muMom.idx }); }
      else if ( muMom.pdg == 22 ) { idx.push_back({ muonIdx , uint(genMuIdx) , muMom.idx }); }
      else { std::cout << "[ERROR] Leading Muon in DrellYan MC comes from invalid mother " << muMom.pdg << " !" << std::endl; return false; }
      for (uint iMu = 0; iMu < muonTree->PF_Muon_Mom().size(); iMu++) {
        const int& iGenMu = muonTree->PF_Muon_Gen_Idx()[iMu];
        if ( (iGenMu != genMuIdx) && (iGenMu > -1) ) {
          const auto& mom = muonTree->MuonMother(iGenMu);
          if ( mom.idx == muMom.idx ) { idx.push_back({ iMu , uint(iGenMu) , mom.idx }); break; }
        }
      }
    }
    if ( idx.size() == 0 ) { std::cout << "[ERROR] Found 0 RECO muons matched to GEN muons from Z decay!" << std::endl; return false; }
    // Found 2 RECO muons from DrellYan decay
    if ( idx.size() == 2 ) {
      if ( (idx[0][0] != muonIdx) && (idx[1][0] != muonIdx) ) { std::cout << "[ERROR] Input muon index is inconsistent with the RECO muon from DrellYan->MuMu decay!" << std::endl; return false; }
      if ( idx[0][2] != idx[1][2] ) { std::cout << "[ERROR] The 2 muons don't share the same mother (" << idx[0][2] << " , " << idx[1][2] << ") !" << std::endl; return false; }
      const TLorentzVector& DY_P4 = muonTree->PF_Muon_Mom()[idx[0][0]] + muonTree->PF_Muon_Mom()[idx[1][0]];
      // Reference is DY pT
      reference_pT.SetMagPhi( DY_P4.Pt(), DY_P4.Phi() );
      // Boson pT is the sum of both RECO Muons pT
      boson_pT.SetMagPhi( DY_P4.Pt(), DY_P4.Phi() );
    }
    // Found 1 RECO muon from Drell Yan decay
    if ( idx.size() == 1 ) {
      if ( idx[0][0] != muonIdx ) { std::cout << "[ERROR] Input muon index is inconsistent with the RECO muon from DrellYan->Mu decay!" << std::endl; return false; }
      // Find the GEN muon that was not matched to a RECO muon
      std::vector< uint > iGenMissingMuon;
      for (uint iGenMu = 0; iGenMu < muonTree->Gen_Muon_Mom().size(); iGenMu++) {
        if (iGenMu != idx[0][1]) {
          const auto& mom = muonTree->MuonMother(iGenMu);
          if (mom.idx == idx[0][2]) { iGenMissingMuon.push_back(iGenMu); }
        }
      }
      if ( iGenMissingMuon.size() != 1 ) { std::cout << "[ERROR] Inconsistent number of GEN muons (" << iGenMissingMuon.size() << ") from Z decay!" << std::endl; return false; }
      // Reference is RECO Muon pT
      reference_pT.SetMagPhi( muonTree->PF_Muon_Mom()[idx[0][0]].Pt() , muonTree->PF_Muon_Mom()[idx[0][0]].Phi() );
      // Boson pT is the sum of RECO Muon pT and missing GEN Muon pT
      boson_pT.SetMagPhi( muonTree->Gen_Muon_Mom()[iGenMissingMuon[0]].Pt() , muonTree->Gen_Muon_Mom()[iGenMissingMuon[0]].Phi() );
      boson_pT += reference_pT;
    }
  }
  //
  // Case: W -> Tau -> Muon + Neutrinos
  else if (sample.find("MC_WToTau")!=std::string::npos) {
    // Find the number of RECO muons matched to GEN muons from W->Tau boson decay
    std::vector< uint > idx;
    const int& iGenMu = muonTree->PF_Muon_Gen_Idx()[muonIdx];
    if ( iGenMu > -1 ) {
      const auto& mom = muonTree->findMuonMother(iGenMu, 24, 2);
      if ( mom.pdg == 24 ) { idx = { muonIdx , uint(iGenMu) , mom.idx }; }
    }
    if ( idx.size() == 0 ) { std::cout << "[ERROR] Found 0 RECO muons matched to GEN muon from W->Tau decay!" << std::endl; return false; }
    // Found 1 RECO muon from W->Tau decay
    // Reference is RECO Muon pT
    reference_pT.SetMagPhi( muonTree->PF_Muon_Mom()[idx[0]].Pt() , muonTree->PF_Muon_Mom()[idx[0]].Phi() );
    // Boson pT is the GEN W boson pT (mother of muon)
    boson_pT.SetMagPhi( muonTree->Gen_Particle_Mom()[idx[2]].Pt() , muonTree->Gen_Particle_Mom()[idx[2]].Phi() );
  }
  //
  // Case: W -> Muon + Neutrino
  else if (sample.find("MC_W")!=std::string::npos) {
    // Find the number of RECO muons matched to GEN muons from W boson decay
    std::vector< uint > idx;
    const int& iGenMu = muonTree->PF_Muon_Gen_Idx()[muonIdx];
    if ( iGenMu > -1 ) {
      const auto& mom = muonTree->MuonMother(iGenMu);
      if ( mom.pdg == 24 ) { idx = { muonIdx , uint(iGenMu) , mom.idx }; }
    }
    if ( idx.size() == 0 ) { std::cout << "[ERROR] Found 0 RECO muons matched to GEN muon from W decay!" << std::endl; return false; }
    // Found 1 RECO muon from W decay
    // Find the GEN neutrino
    int iGenNeutrino = -1;
    for (uint i = 0; i < muonTree->Gen_Particle_PdgId().size(); i++) {
      if ( (muonTree->Gen_Particle_Status()[i]==1) && (abs(muonTree->Gen_Particle_PdgId()[i])==14) ) {
        const auto& mom = muonTree->Mother(i);
        if (mom.idx == idx[2]) { iGenNeutrino = i; break; }
      }
    }
    if ( iGenNeutrino == -1 ) { std::cout << "[ERROR] No GEN neutrinos were found from W decay!" << std::endl; return false; }
    // Reference is RECO Muon pT
    reference_pT.SetMagPhi( muonTree->PF_Muon_Mom()[idx[0]].Pt() , muonTree->PF_Muon_Mom()[idx[0]].Phi() );
    // Boson pT is the sum of RECO Muon pT and missing GEN Neutrino pT
    boson_pT.SetMagPhi( muonTree->Gen_Particle_Mom()[iGenNeutrino].Pt() , muonTree->Gen_Particle_Mom()[iGenNeutrino].Phi() );
    boson_pT += reference_pT;
  }
  //
  // Case: t -> W -> Muon + Neutrinos
  else if (sample.find("MC_TT")!=std::string::npos) {
    // Find the number of RECO muons matched to GEN muons from W->Muon boson decay
    std::vector< uint > idx;
    const int& iGenMu = muonTree->PF_Muon_Gen_Idx()[muonIdx];
    if ( iGenMu > -1 ) {
      const auto& mom = muonTree->MuonMother(iGenMu);
      if ( mom.pdg == 24 ) { idx = { muonIdx , uint(iGenMu) , mom.idx }; }
    }
    if ( idx.size() == 0 ) { std::cout << "[ERROR] Found 0 RECO muons matched to GEN muon from t->W->Mu decay!" << std::endl; return false; }
    // Found 1 RECO muon from W->Muon decay
    // Reference is RECO Muon pT
    reference_pT.SetMagPhi( muonTree->PF_Muon_Mom()[idx[0]].Pt() , muonTree->PF_Muon_Mom()[idx[0]].Phi() );
    // Boson pT is the GEN W boson pT (mother of muon)
    boson_pT.SetMagPhi( muonTree->Gen_Particle_Mom()[idx[2]].Pt() , muonTree->Gen_Particle_Mom()[idx[2]].Phi() );
  }
  //
  // Case: QCD -> Muon + X
  else if (sample.find("MC_QCD")!=std::string::npos) {
    // Reference is RECO Muon pT
    reference_pT.SetMagPhi( muonTree->PF_Muon_Mom()[muonIdx].Pt() , muonTree->PF_Muon_Mom()[muonIdx].Phi() );
    // Boson pT is also the RECO Muon pT
    boson_pT = reference_pT;
  }
  //
  else {
    std::cout << "[ERROR] Sample : " << sample << " has not been implemented in the Recoil Corrector Class!" << std::endl; return false;
  }
  //
  // Return back
  //
  return true;
};


//-----------------------------------------------------------------------------------------------------------------------------------------
void RecoilCorrector::computeMET(
                                 TVector2& MET,
                                 const double& u1Val,
                                 const double& u2Val,
                                 const TVector2& reference_pT,
                                 const TVector2& boson_pT
                                 )
{
  // Initialize output variables
  MET = TVector2();
  //
  // Compute the MET vector based on the input Recoil components
  //
  // MET X Component
  const double MET_X = ( -1.0 * reference_pT.Mod() * TMath::Cos( reference_pT.Phi() ) ) - ( u1Val * TMath::Cos( boson_pT.Phi() ) ) + ( u2Val * TMath::Sin( boson_pT.Phi() ) );
  // MET Y Component
  const double MET_Y = ( -1.0 * reference_pT.Mod() * TMath::Sin( reference_pT.Phi() ) ) - ( u1Val * TMath::Sin( boson_pT.Phi() ) ) - ( u2Val * TMath::Cos( boson_pT.Phi() ) );
  // MET Vector
  MET.Set(MET_X, MET_Y);
};


//-----------------------------------------------------------------------------------------------------------------------------------------
void RecoilCorrector::computeRecoil(
                                    double& u1Val,
                                    double& u2Val,
                                    const TVector2& reference_pT,
                                    const TVector2& boson_pT,
                                    const TVector2& MET
                                    )
{
  // Initialize output variables
  u1Val = -9999.0; u1Val = -9999.0;
  //
  // Compute Recoil components
  //
  const TVector2 U = -1.0 * ( reference_pT + MET );
  // Compute Recoil Parallel component with respect to Boson pT
  u1Val = ( U * boson_pT ) / boson_pT.Mod();
  // Compute Recoil Perpendicular component with respect to Boson pT
  u2Val = ( ( U.Px() * boson_pT.Py() ) - ( U.Py() * boson_pT.Px() ) ) / boson_pT.Mod();
};


//-----------------------------------------------------------------------------------------------------------------------------------------
bool RecoilCorrector::computeRecoilCDF(
                                       double& CDF,
                                       const double& uVal,
                                       const std::map< std::string , double >& uPar,
                                       const bool useOneGaussian
                                       )
{
  // Initialize output variables
  CDF = -1.0;
  //
  // Evaluate the model
  //
  auto fcnCDF = std::unique_ptr<TF1>(getFunctionCDF(uPar, useOneGaussian));
  if (fcnCDF==NULL) { return false; }
  CDF = fcnCDF->Eval(uVal);
  if ( (CDF < 0.0) || (CDF > 1.0) ) { std::cout << "[ERROR] Computed CDF (" << CDF << ") is out of range!" << std::endl; return false; }
  //
  // Return back
  //
  return true;
};


//-----------------------------------------------------------------------------------------------------------------------------------------
bool RecoilCorrector::computeRecoilInvCDF(
                                          double& uVal,
                                          const double& CDF,
                                          const std::map< std::string , double >& uPar,
                                          const bool useOneGaussian
                                          )
{
  //
  if ( (CDF < 0.0) || (CDF > 1.0) ) { std::cout << "[ERROR] Input CDF (" << CDF << ") is out of range!" << std::endl; return false; }
  //
  // Initialize output variables
  uVal = -9999.0;
  //
  // Determine which CDF model to use
  //
  const int model = findModel(uPar);
  if (model<1) { return false; }
  //
  // Determine the inverse CDF model
  //
  if (model>=1) {
    // Use the ROOT::Math Gaussian CDF inverse (also used as initial guess for multiple gaussian models)
    uVal = uPar.at("mean") + ROOT::Math::gaussian_quantile( CDF , uPar.at("sigma") );
  }
  if (!useOneGaussian && model>1) {
    // Build CDF function
    auto fcnCDF = std::unique_ptr<TF1>(getFunctionCDF(uPar, useOneGaussian, CDF));
    if (fcnCDF==NULL) { return false; }
    // Initiliaze the RootFinder
    ROOT::Math::RootFinder rf(ROOT::Math::RootFinder::kGSL_STEFFENSON);
    // Assign the CDF function
    ROOT::Math::GradFunctor1D wf(*fcnCDF);
    // Use the 1 Gaussian inverse CDF as initial guess value
    rf.SetFunction(wf, uVal);
    // Proceed to find the root value
    rf.Solve();
    if (rf.Status() != 0) {
      std::cout << "[WARNING] Computing the inverse CDF for model " << model << " failed with status " << rf.Status() << std::endl;
      uVal = uPar.at("mean") + ROOT::Math::gaussian_quantile( CDF , uPar.at("sigma") );
    }
    else { double uVal = rf.Root(); }
    // Check result
    if (uVal <= -1000. || uVal >= 1000.) { std::cout << "[ERROR] The inverse CDF for model " << model << " was not found!" << std::endl; return false; }
  }
  //
  // Return back
  //
  return true;
};


//-----------------------------------------------------------------------------------------------------------------------------------------
bool RecoilCorrector::applyScalingRecoilCorrection(
                                                   double& uVal_CORR,
                                                   const double& uVal_RAW,
                                                   const std::map< std::string , double >& uPar_MC,
                                                   const std::map< std::string , double >& uPar_DATA,
                                                   const bool useOneGaussian_MC,
                                                   const bool useOneGaussian_DATA
                                                   )
{
  //
  if (uPar_MC.at("sigma") == 0.) { std::cout << "[ERROR] MC sigma used in Scaling method is zero!" << std::endl; return false; }
  //
  // Initialize output variables
  uVal_CORR = -9999.0;
  //
  const int model_MC = findModel(uPar_MC);
  const int model_DATA = findModel(uPar_DATA);
  //
  // Apply scaling method to correct the MC Recoil
  //
  // Simple case: Both MC and Data are modeled with 1 Guassian
  if ( ( useOneGaussian_MC || (model_MC == 1) )  && ( useOneGaussian_DATA || (model_DATA == 1) ) ) {
    const double z_MC = ( uVal_RAW - uPar_MC.at("mean") ) / uPar_MC.at("sigma");
    uVal_CORR = uPar_DATA.at("mean") + (z_MC * uPar_DATA.at("sigma"));
  }
  // General case: Either the MC or the Data are modeled with more than 1 Gaussian
  else {
    double CDF_MC = -1.0;
    if (!computeRecoilCDF(CDF_MC, uVal_RAW, uPar_MC, useOneGaussian_MC)) { return false; }
    if ( (CDF_MC > 0.0) && (CDF_MC < 1.0) ) {
      if (!computeRecoilInvCDF(uVal_CORR, CDF_MC, uPar_DATA, useOneGaussian_DATA)) { return false; }
    }
    else {
      std::cout << "[WARNING] Can't compute inverse CDF since CDF_MC is " << CDF_MC << " , will go back to simple case!" << std::endl;
      const double z_MC = ( uVal_RAW - uPar_MC.at("mean") ) / uPar_MC.at("sigma");
      uVal_CORR = uPar_DATA.at("mean") + (z_MC * uPar_DATA.at("sigma"));
    } 
  }
  //
  // Return back
  //
  return true;
};


//-----------------------------------------------------------------------------------------------------------------------------------------
bool RecoilCorrector::applySmearingRecoilCorrection(
                                                    double& uVal_CORR,
                                                    const double& uVal_RAW,
                                                    const std::map< std::string , double >& uPar_MC,
                                                    const std::map< std::string , double >& uPar_DATA
                                                    )
{
  // Initialize output variables
  uVal_CORR = -9999.0;
  //
  // Determine which CDF model to use
  //
  const int model_MC = findModel(uPar_MC);
  if (model_MC<1) { return false; }
  const int model_DATA = findModel(uPar_DATA);
  if (model_DATA<1) { return false; }
  //
  // Apply smearing method to correct the MC Recoil
  //
  // Compute the mean value ( u - u_MC + u_DATA )
  double mean = ( uVal_RAW - uPar_MC.at("mean") + uPar_DATA.at("mean") );
  // Compute the width  sqrt( u_DATA^2 - u_MC^2 ) 
  double sigma = ( TMath::Power(uPar_DATA.at("sigma"), 2.0) - TMath::Power(uPar_MC.at("sigma"), 2.0) );
  if (sigma<0.0) { std::cout << "[Warning] Variance in MC (" <<  uPar_MC.at("sigma") << ") is larger than in DATA (" << uPar_DATA.at("sigma") << ") !" << std::endl; }
  if (sigma<0.0) { std::cout << "[Info] Setting sigma for smearing algorithm to zero!" << std::endl; sigma = 0.0; }
  sigma = TMath::Sqrt(sigma); // Take the square root of the variance to compute the gauusian width
  // Sample a random number based on a gaussian distribution
  if (sigma==0) { uVal_CORR = mean; }
  else { uVal_CORR = fRandom_->Gaus(mean, sigma); }
  //
  // Return back
  //
  return true;
};


//-----------------------------------------------------------------------------------------------------------------------------------------
bool RecoilCorrector::correctMCRecoil(
                                      double& uVal_CORR,
                                      const double& uVal_RAW,
                                      const std::map< std::string , double >& uPar_MC,
                                      const std::map< std::string , double >& uPar_DATA,
                                      const std::string& method,
                                      const bool useOneGaussian_MC,
                                      const bool useOneGaussian_DATA
                                      )
{
  // Initialize output variables
  uVal_CORR = -9999.0;
  //
  // Correct the MC Recoil
  //
  if (method == "Scaling") {
    if (!applyScalingRecoilCorrection(uVal_CORR, uVal_RAW, uPar_MC, uPar_DATA, useOneGaussian_MC, useOneGaussian_DATA)) { return false; };
  }
  if (method == "Smearing") {
    if (!applySmearingRecoilCorrection(uVal_CORR, uVal_RAW, uPar_MC, uPar_DATA)) { return false; }
  }
  //
  // Return back
  //
  return true;
};


//-----------------------------------------------------------------------------------------------------------------------------------------
bool RecoilCorrector::correctMCRecoil(
                                      double& uVal_CORR,
                                      const double& uVal_RAW,
                                      const double& pT,
                                      const std::map< std::string , std::unique_ptr<TF1> >& uFit_MC,
                                      const std::map< std::string , std::unique_ptr<TF1> >& uFit_DATA,
                                      const std::string& method,
                                      const bool useOneGaussian_MC,
                                      const bool useOneGaussian_DATA,
                                      const bool corrMean
                                      )
{
  // Initialize output variables
  uVal_CORR = -9999.0;
  //
  // Evaluate the Recoil Fit functions
  //
  // For MC
  std::map< std::string , double > uPar_MC;
  evalRecoilFits(uPar_MC, pT, uFit_MC, corrMean);
  // For DATA
  std::map< std::string , double > uPar_DATA;
  evalRecoilFits(uPar_DATA, pT, uFit_DATA, corrMean);
  //
  if (corrMean==false) { for(const auto& u : uPar_MC) { if (u.first.find("mean")!=std::string::npos) { uPar_DATA[u.first] = u.second; } } }
  //
  // Correct the MC Recoil
  //
  if (!correctMCRecoil(uVal_CORR, uVal_RAW, uPar_MC, uPar_DATA, method, useOneGaussian_MC, useOneGaussian_DATA)) { return false; }
  //
  // Return back
  //
  return true;
};


//-----------------------------------------------------------------------------------------------------------------------------------------
bool RecoilCorrector::correctMET(
                                 TVector2& MET_CORR,
                                 const TVector2& MET_RAW,
                                 const std::string& method
                                 )
{
  //
  if ( (reference_pT_.Px() == 0.0) && (reference_pT_.Py() == 0.0) ) { std::cout << "[ERROR] The reference pT has not been defined!" << std::endl; return false; }
  if ( (    boson_pT_.Px() == 0.0) && (    boson_pT_.Py() == 0.0) ) { std::cout << "[ERROR] The boson pT has not been defined!"     << std::endl; return false; }
  if ( u1Fits_MC_.count(METType_) == 0 ) { std::cout << "[ERROR] MET type " << METType_ << " has not been defined!" << std::endl; return false; }
  //
  // Initialize output variables
  MET_CORR = TVector2();
  //
  // Compute the Recoil components
  //
  double u1_RAW, u2_RAW;
  computeRecoil(u1_RAW, u2_RAW, reference_pT_, boson_pT_, MET_RAW);
  //
  // Correct the Recoil components
  //
  // For parallel recoil
  double u1_CORR;
  if(!correctMCRecoil(u1_CORR, u1_RAW, boson_pT_.Mod(), u1Fits_MC_.at(METType_), u1Fits_DATA_.at(METType_), method, useOneGaussian_MC_, useOneGaussian_DATA_, corrMean_)) { return false; }
  // For perpendicular recoil
  double u2_CORR;
  if(!correctMCRecoil(u2_CORR, u2_RAW, boson_pT_.Mod(), u2Fits_MC_.at(METType_), u2Fits_DATA_.at(METType_), method, useOneGaussian_MC_, useOneGaussian_DATA_, corrMean_)) { return false; }
  //
  // Compute the MET based on the corrected recoil components
  //  
  computeMET(MET_CORR, u1_CORR, u2_CORR, reference_pT_, boson_pT_);
  //
  // Return back
  //
  return true;
};

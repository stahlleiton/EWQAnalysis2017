#ifndef RecoilCorrections_H
#define RecoilCorrections_H

// Auxiliary Headers
#include "../../Utilities/HiMuonTree.h"
// ROOT Headers
#include <TRandom3.h>
#include <TF1.h>                      // 1D function
// C++ Headers
#include <memory>
#include <string>
#include <map>
#include <vector>
#include <ctime>                      // time

//
// ** apply recoil corrections **
// 
// usage: 
//    double met=rawMetValue;
//    double metphi=rawMetPhiValue;
//    RecoilCorrector corrector;
//    corrector->Correct(met,metphi,GenZPt,GenZPhi,leptonPt,leptonPhi);
//  
//
// where leptonPt, leptonPhi are dilepton kinematics for z->ll and single lepton kinematics for w->lnu
//

class RecoilCorrector
{
  
 public:

  RecoilCorrector(const int seed = time(NULL));
  ~RecoilCorrector();

  //
  // Public Functions
  //
  void  clear();
  bool  isValid();
  bool  setInputFiles       ( const std::string& METType , const std::string& mcFileName , const std::string& dataFileName );
  void  setInitialSetup     ( const std::string& sample , const bool useOneGaussian_MC = false, const bool useOneGaussian_DATA = false , const bool corrMean = true );
  bool  correctMET          ( TVector2& MET_CORR , const TVector2& MET_RAW , const std::string& method );
  void  setPt               ( const TVector2& reference_pT , const TVector2& boson_pT ) { this->reference_pT_ = reference_pT; this->boson_pT_ = boson_pT; };

  static bool getPtFromTree ( TVector2& reference_pT , TVector2& boson_pT , const uint& muonIdx , const std::unique_ptr<HiMuonTree>& muonTree , const std::string& sample );


 protected:

  //
  // Protected Functions
  //
  bool  getRecoilFits       ( const std::string& fileName , std::map< std::string , std::unique_ptr<TF1> >& u1Fit , std::map< std::string , std::unique_ptr<TF1> >& u2Fit );
  void  evalRecoilFits      ( std::map< std::string , double >& uPar , const double& pT , const std::map< std::string , std::unique_ptr<TF1> >& uFit , const bool setMeanToZero = false );
  int   findModel           ( const std::map< std::string , double >& uPar );
  TF1*  getFunctionCDF      ( const std::map< std::string , double >& uPar , const double& CDF , const bool useOneGaussian = false );
  //
  void  computeMET          ( TVector2& MET , const double& u1Val , const double& u2Val , const TVector2& reference_pT , const TVector2& boson_pT );
  void  computeRecoil       ( double& u1Val , double& u2Val , const TVector2& reference_pT , const TVector2& boson_pT , const TVector2& MET );
  bool  computeRecoilCDF    ( double& CDF , const double& uVal , const std::map< std::string , double >& uPar , const bool useOneGaussian = false );
  bool  computeRecoilInvCDF ( double& uVal , const double& CDF , const std::map< std::string , double >& uPar , const bool useOneGaussian = false );
  //
  bool  applyScalingRecoilCorrection  ( double& uVal_CORR ,  const double& uVal_RAW , const std::map< std::string , double >& uPar_MC , 
                                        const std::map< std::string , double >& uPar_DATA , const bool useOneGaussian_MC = false , const bool useOneGaussian_DATA = false );
  bool  applySmearingRecoilCorrection ( double& uVal_CORR , const double& uVal_RAW , const std::map< std::string , double >& uPar_MC , const std::map< std::string , double >& uPar_DATA );
  //
  bool  correctMCRecoil     ( double& uVal_CORR , const double& uVal_RAW , const std::map< std::string , double >& uPar_MC , const std::map< std::string , double >& uPar_DATA , 
                              const std::string& method , const bool useOneGaussian_MC = false , const bool useOneGaussian_DATA = false );
  bool  correctMCRecoil     ( double& uVal_CORR , const double& uVal_RAW , const double& pT , const std::map< std::string , std::unique_ptr<TF1> >& uFit_MC , const std::map< std::string , std::unique_ptr<TF1> >& uFit_DATA , 
                              const std::string& method , const bool useOneGaussian_MC = false , const bool useOneGaussian_DATA = false , const bool corrMean = true , const bool setMeanToZero = false );


  //
  // Global Variables
  //
  std::unique_ptr<TRandom3> fRandom_;
  //
  bool useOneGaussian_MC_;
  bool useOneGaussian_DATA_;
  bool corrMean_;
  //
  std::map< std::string , std::map< std::string , std::unique_ptr<TF1> > > u1Fits_MC_;
  std::map< std::string , std::map< std::string , std::unique_ptr<TF1> > > u1Fits_DATA_;
  std::map< std::string , std::map< std::string , std::unique_ptr<TF1> > > u2Fits_MC_;
  std::map< std::string , std::map< std::string , std::unique_ptr<TF1> > > u2Fits_DATA_;
  //
  TVector2  reference_pT_;
  TVector2  boson_pT_;
  //
  std::string METType_;

};



#endif

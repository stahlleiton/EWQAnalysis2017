#ifndef Utilities_resultUtils_h
#define Utilities_resultUtils_h
// Auxiliary Headers
#include "bin.h"
#include "../../Utilities/EVENTUTILS.h"
// ROOT headers
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraphAsymmErrors.h"
#include "TEfficiency.h"
#include "TVectorD.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TFrame.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TPaletteAxis.h"
#include "TMatrixDSym.h"
// RooFit headers
// c++ headers
#include <dirent.h>
#include <iostream>
#include <fstream>
#include <map>
#include <array>
// CMS headers
#include "../../Utilities/CMS/tdrstyle.C"
#include "../../Utilities/CMS/CMS_lumi.C"
// MCFM result
#include "MCFM.h"
#include "MCFM_v2.h"
#include "MCFM_REW.h"
#include "HIN13007.h"


const double LUMIUNC_ = 0.035; //0.035;


// ------------------ TYPE -------------------------------
typedef std::map< std::string , std::string       > StringMap;
typedef std::map< std::string , std::string*      > StringPMap;
typedef std::map< std::string , bool              > BoolMap;
typedef std::map< std::string , TGraphAsymmErrors > GraphMap;
typedef std::map< std::string , GraphMap          > GraphMapMap;
typedef std::map< std::string , GraphMapMap       > GraphTriMap;
typedef std::map< std::string , GraphTriMap       > GraphQuadMap;
typedef std::map< std::string , GraphQuadMap      > GraphPentaMap;
typedef std::map< anabin<0>   , double            > BinMap;
typedef std::map< std::string , BinMap            > BinMapMap;
typedef std::map< std::string , BinMapMap         > BinTriMap;
typedef std::map< std::string , BinTriMap         > BinQuadMap;
typedef std::map< std::string , BinQuadMap        > BinPentaMap;
typedef std::map< std::string , double            > DoubleMap;
typedef std::map< std::string , DoubleMap         > DoubleMapMap;
typedef std::map< std::string , DoubleMapMap      > DoubleTriMap;
typedef std::map< anabin<0>   , DoubleTriMap      > DoubleBinTriMap;
typedef std::map< std::string , DoubleBinTriMap   > VarBinMap;
typedef std::pair< VarBinMap  , uint              > VarBinPairMap;
typedef std::vector< BinPentaMap                  > BinPentaMapVec;
typedef std::map< std::string , BinPentaMapVec    > BinSextaMapVec;
typedef std::map< std::string , BinSextaMapVec    > BinSeptaMapVec;
typedef std::map< std::string , BinPentaMap       > BinSextaMap;
typedef std::map< std::string , BinSextaMap       > BinSeptaMap;
typedef std::map< std::string , std::vector< double > > DoubleVecMap;
typedef std::map< std::string , std::pair< uint , uint > > CorrMap;
typedef std::map< std::string , TMatrixDSym       > CovMatrixMap;
typedef std::map< std::string , CovMatrixMap      > CovMatrixDiMap;
typedef std::map< std::string , std::map< std::string , std::pair< std::vector< std::string > , uint > > > WSDirMap;


// Tree Info Structure (wrapper to carry information around)
typedef struct TreeInfo {
  DoubleMapMap     Var;
  StringMap        Str;
  StringPMap       StrP;
  BoolMap          Flag;
  void             Clear() { this->Var.clear(); this->Str.clear(); this->StrP.clear(); this->Flag.clear(); }
  TreeInfo() {}
  TreeInfo(const TreeInfo &ref, bool keep = false) {
    this->Copy(ref, keep);
  }
  ~TreeInfo() {
    for (auto& p : this->StrP) { if (p.second!=NULL) { delete p.second; } }
    this->Clear();
  }
  void Copy(const DoubleMapMap &ref, bool keep = true) {
    if (!keep) this->Var.clear();
    for (const auto& var : ref) {
      for (const auto& ele : var.second) {
        this->Var[var.first][ele.first] = ele.second;
      }
    }
  }
  void Copy(const StringMap &ref, bool keep = true) {
    if (!keep) this->Str.clear();
    for (const auto& str : ref) {
      this->Str[str.first] = str.second;
    }
  }
  void Copy(const BoolMap &ref, bool keep = true) {
    if (!keep) this->Flag.clear();
    for (const auto& flag : ref) {
      this->Flag[flag.first] = flag.second;
    }
  }
  void Copy(const TreeInfo &ref, bool keep = true) {
    this->Copy(ref.Var, keep);
    this->Copy(ref.Str, keep);
    this->Copy(ref.Flag, keep);
  }
  bool operator == (const DoubleMapMap &ref) const 
  {
    if (ref.size() != this->Var.size()) return false;
    for (const auto& var : this->Var) {
      if (ref.count(var.first)==0) return false;
      for (const auto& ele : var.second) {
        if (ref.at(var.first).count(ele.first)==0) return false;
        if (ele.second != ref.at(var.first).at(ele.first)) return false;
      }
    }
    return true;
  }
  bool operator == (const StringMap &ref) const 
  {
    if (ref.size() != this->Str.size()) return false;
    for (auto& str : this->Str) {
      if (ref.count(str.first)==0) return false;
      if (str.second != ref.at(str.first)) return false;
    }
    return true;
  }
  bool operator == (const BoolMap &ref) const 
  {
    if (ref.size() != this->Flag.size()) return false;
    for (auto& flag : this->Flag) {
      if (ref.count(flag.first)==0) return false;
      if (flag.second != ref.at(flag.first)) return false;
    }
    return true;
  }
  bool operator != ( const DoubleMapMap  &ref ) const { if (*this == ref) { return false; } return true; }
  bool operator != ( const StringMap     &ref ) const { if (*this == ref) { return false; } return true; }
  bool operator != ( const BoolMap       &ref ) const { if (*this == ref) { return false; } return true; }
  bool operator == (const TreeInfo  &ref) const 
  {
    if ( *this != ref.Var  ) return false;
    if ( *this != ref.Str  ) return false;
    if ( *this != ref.Flag ) return false;
    return true;
  }
} TreeInfo;


// ------------------ FUNCTION -------------------------------

bool fileList(std::vector< std::string >& fileNames, const std::string& dirPath)
{
  // Open the directory
  DIR * dpdf = opendir(dirPath.c_str());
  // Search for all the files inside the directory
  if (dpdf != NULL){
    struct dirent *epdf;
    while ((epdf = readdir(dpdf))){
      if (strcmp(epdf->d_name,".")!=0 && strcmp(epdf->d_name,"..")!=0 ) {
        //std::cout << "[INFO] Adding file: " << epdf->d_name << std::endl;
        fileNames.push_back(epdf->d_name);
      }
    }
  } else {
    std::cout << "[ERROR] Working directory ( " << dirPath << " ) was not found!" << endl; return false;
  }
  return true;
};


void iniResultsTreeInfo( TreeInfo& info , const std::string& thePoiNames )
{
  //
  // Define the names of the Dataset RooRealVars
  const std::vector< std::string > dsVarNames = { "MET" , "Muon_Pt" , "Muon_Eta" , "Muon_Iso" , "Muon_MT" , "Centrality" };
  const std::vector< std::string > varType    = { "Min" , "Max" , "Val" , "Err" };
  //
  for (const auto& v : dsVarNames) {
    for (const auto& t : varType) {
      info.Var["VAR_"+v][t] = -99.0;
    }
  }
  //
  // Define the POI names
  std::string thePoiNamesStr = thePoiNames;
  if (thePoiNamesStr == "all") {
    thePoiNamesStr  = "N_WToMu,N_WToTauToMu,N_DYToMu,N_DYToTauToMu,N_QCDToMu,N_TTbarToMu";
    thePoiNamesStr += ",Alpha_QCDToMu,Beta_QCDToMu,x0_QCDToMu,Sigma0_QCDToMu,Sigma1_QCDToMu,Sigma2_QCDToMu";
  }
  thePoiNamesStr = thePoiNamesStr + ","; // Just to include the last name XD
  std::vector< std::string > thePoiNamesVec;
  while (thePoiNamesStr.find(",")!=std::string::npos) {
    thePoiNamesVec.push_back( thePoiNamesStr.substr(0, thePoiNamesStr.find(",")) );
    thePoiNamesStr = thePoiNamesStr.substr(thePoiNamesStr.find(",")+1);
  }
  const std::vector< std::string > poiVarType = { "Min" , "Max" , "Val" , "Err" , "ErrLo" , "ErrHi" , "parIni_Val" , "parIni_Err" };
  //
  for (const auto& p : thePoiNamesVec) {
    for (const auto& t : poiVarType) {
      info.Var["POI_"+p][t] = -99.0;
    }
  }
  //
  // Define the names of the remaining variables
  info.Var["Luminosity"]["Val"] = -1.0;
  info.Var["N_DS_Entries"]["Val"] = -1.0;
  info.Var["N_FIT_Entries"]["Val"] = -1.0;
  info.Var["TEST_FIT"]["Chi2"] = -1.0;
  info.Var["TEST_FIT"]["NDoF"] = -1.0;
  info.Var["TEST_FIT"]["Val"] = -1.0;
  info.Var["TEST_FIT_KS"]["KS"] = -1.0;
  info.Var["TEST_FIT_KS"]["Val"] = -1.0;
  info.Var["TEST_FIT_AD"]["AD"] = -1.0;
  info.Var["TEST_FIT_AD"]["Val"] = -1.0;
  info.Var["TEST_FIT_BCChi2"]["Chi2"] = -1.0;
  info.Var["TEST_FIT_BCChi2"]["NDoF"] = -1.0;
  info.Var["TEST_FIT_BCChi2"]["Val"] = -1.0;
  info.Var["TEST_FIT_PChi2"]["Chi2"] = -1.0;
  info.Var["TEST_FIT_PChi2"]["NDoF"] = -1.0;
  info.Var["TEST_FIT_PChi2"]["Val"] = -1.0;
  //
  // Initialize the model name containers
  const std::vector< std::string > objType = { "W" , "WToTau" , "DY" , "DYToTau" , "TTbar" , "QCD" };
  //
  for (const auto& o : objType) {
    info.Str["Model_"+o] = "";
  }
  //
  // Initialize the cut value containers
  info.Str["cutAndCount_W"] = "";
  //
  // Initiliaze the collision and charge information
  info.Str["collSystem"] = "";
  info.Str["charge"] = "";
  info.Str["fitObject"] = "";
  info.Str["channel"] = "";
  info.Str["metType"] = "";
  //
  // Boolean flags
  info.Flag["useEtaCM"] = false;
  //
};


void setBranches( TTree& tree , TreeInfo& info )
{
  //
  // Create Tree branches for variables
  for (auto& v : info.Var) {
    for (auto& p : v.second) {
      tree.Branch(Form("%s_%s", v.first.c_str(), p.first.c_str()) , &(p.second) , Form("%s_%s/D", v.first.c_str(), p.first.c_str()));
    }
  }
  //
  // Create Tree branches for flags
  for (auto& f : info.Flag) {
    tree.Branch(Form("%s", f.first.c_str()) , &(f.second) , Form("%s/O", f.first.c_str()));
  }
  //
  // Create Tree branches for char arrays
  for (auto& c : info.Str) {
    tree.Branch(Form("%s", c.first.c_str()) , &(c.second));
  }
  //
};


void setBranchAddress( TTree& tree , TreeInfo& info )
{
  //
  // Create Tree branches for variables
  for (auto& v : info.Var) {
    for (auto& p : v.second) {
      if (tree.GetBranch(Form("%s_%s", v.first.c_str(), p.first.c_str()))) {
        tree.SetBranchAddress(Form("%s_%s", v.first.c_str(), p.first.c_str()) , &(p.second));
      }
    }
  }
  //
  // Create Tree branches for flags
  for (auto& f : info.Flag) {
    if (tree.GetBranch(Form("%s", f.first.c_str()))) {
      tree.SetBranchAddress(Form("%s", f.first.c_str()) , &(f.second));
    }
  }
  //
  // Create Tree branches for char arrays
  for (auto& c : info.Str) {
    if (tree.GetBranch(Form("%s", c.first.c_str()))) {
      tree.SetBranchAddress(Form("%s", c.first.c_str()) , &(info.StrP[c.first]));
    }
  }
  //
};


void iniAcceptanceAndEfficiency( BinPentaMap& eff , const VarBinMap& inputVar )
{
  //
  std::vector< std::string > lbl = { "Val" , "Err_Stat_High" , "Err_Stat_Low" , "Err_Syst_High" , "Err_Syst_Low" , "Err_Tot_High" , "Err_Tot_Low" };
  //
  for (const auto& c : inputVar) {
    for (const auto& b : c.second) {
      for (const auto& ch : b.second) {
        // For MC Acceptance
        for (const auto& l : lbl) {
          eff[c.first][ch.first]["Acceptance_MC"][l][b.first] = ( (l=="Val") ? 1.0 : 0.0 );
          // For MC Efficiency
          eff[c.first][ch.first]["Efficiency_MC"][l][b.first] = ( (l=="Val") ? 1.0 : 0.0 );
          // For Corrected Efficiency
          eff[c.first][ch.first]["Efficiency_TnP"][l][b.first] = ( (l=="Val") ? 1.0 : 0.0 );
        }
      }
    }
  }
};


double sumErrors( const double& ErrorA , const double& ErrorB )
{
  const double errorA2 = std::pow( ErrorA , 2.0 );
  const double errorB2 = std::pow( ErrorB , 2.0 );
  return ( std::sqrt( errorA2 + errorB2 ) );
};
// Propagation Method: 0: Fully Uncorrelated
// Variation   Method: 4: Fully Uncorrelated , 5: Fully Correlated
// Mixed       Method: 1: Eta Correlated , 2: Charge Correlated , 3: Fully Correlated
CorrMap effTnPType_ = {
  { "TnP_Stat_Iso"    , { 100 , 1 } } , { "TnP_Stat_MuID"    , { 100 , 1 } } , { "TnP_Stat_Trig" , { 2 , 1 } } ,
  { "TnP_Syst_BinIso" , {   1 , 2 } } , { "TnP_Syst_BinMuID" , {   1 , 2 } } , { "TnP_Syst_Iso"  , { 2 , 2 } } , { "TnP_Syst_MuID" , { 2 , 2 } } , { "TnP_Syst_Trig" , { 2 , 2 } } ,
  { "TnP_Syst_PU"     , {   2 , 2 } } , { "TnP_Syst_STA"     , {   2 , 2 } } ,
  { "MC_Syst_PDF"     , {  96 , 5 } } , { "MC_Syst_Scale"    , {   6 , 5 } } , { "MC_Syst_Alpha" , { 2 , 5 } }
};
std::vector< double > absEtaTnP_ = { 0.0 , 1.2 , 2.1 , 2.4 };
std::vector< double > etaTnP_    = { -2.4 , -2.1 , -1.6 , -1.2 , -0.9 , -0.6 , -0.3 , 0.0 , 0.3 , 0.6 , 0.9 , 1.2 , 1.6 , 2.1 , 2.4 };

void fillTnPType( CorrMap& effTnPType , const std::vector< double >& absEtaTnP , const std::vector< double >& etaTnP )
{
  // Initialize the Correction Type Map
  CorrMap effTnP;
  effTnP["TnP_Syst_STA"]     = effTnPType.at("TnP_Syst_STA");
  effTnP["TnP_Syst_PU" ]     = effTnPType.at("TnP_Syst_PU");
  effTnP["TnP_Syst_MuID"]    = effTnPType.at("TnP_Syst_MuID");
  effTnP["TnP_Syst_Iso"]     = effTnPType.at("TnP_Syst_Iso");
  effTnP["TnP_Syst_Trig"]    = effTnPType.at("TnP_Syst_Trig");
  effTnP["TnP_Syst_BinMuID"] = effTnPType.at("TnP_Syst_BinMuID");
  effTnP["TnP_Syst_BinIso"]  = effTnPType.at("TnP_Syst_BinIso");
  //
  for (uint iEta = 1; iEta < absEtaTnP.size(); iEta++) {
    const std::string etaLbl = Form("%s%.0f_%s%.0f", (absEtaTnP[iEta-1]<0.?"m":"p") , std::abs(absEtaTnP[iEta-1])*10. , (absEtaTnP[iEta]<0.?"m":"p") , std::abs(absEtaTnP[iEta])*10.);
    effTnP[Form("TnP_Stat_MuID_%s" , etaLbl.c_str())] = effTnPType.at("TnP_Stat_MuID");
    effTnP[Form("TnP_Stat_Iso_%s"  , etaLbl.c_str())] = effTnPType.at("TnP_Stat_Iso");
  }
  //
  for (uint iEta = 1; iEta < etaTnP.size(); iEta++) {
    const std::string etaLbl = Form("%s%.0f_%s%.0f", (etaTnP[iEta-1]<0.?"m":"p") , std::abs(etaTnP[iEta-1])*10. , (etaTnP[iEta]<0.?"m":"p") , std::abs(etaTnP[iEta])*10.);
    effTnP[Form("TnP_Stat_Trig_%s" , etaLbl.c_str())] = effTnPType.at("TnP_Stat_Trig");
  }
  //
  effTnP["MC_Syst_PDF"]   = effTnPType.at("MC_Syst_PDF");
  effTnP["MC_Syst_Scale"] = effTnPType.at("MC_Syst_Scale");
  effTnP["MC_Syst_Alpha"] = effTnPType.at("MC_Syst_Alpha");
  //
  effTnPType.clear();
  effTnPType = effTnP;
  //
};


void computeTnPUncertainty( double& Error_High, double& Error_Low , const std::vector<double>& variation , const double& nominal , const std::string& type )
{
  double Err_High = 0.0 , Err_Low = 0.0;
  if (variation.size() == 1) {
    Err_High = std::abs(variation[0] - nominal);
    Err_Low  = Err_High;
  }
  else if (type.find("TnP_")!=std::string::npos) {
    if (variation.size() > 2) {
      double sum = 0.0;
      for(uint i = 0; i < variation.size(); i++) {
        const double diff = ( variation[i] - nominal );
        sum += ( diff * diff );
      }
      Err_High = std::sqrt( sum / variation.size() );
      Err_Low  = Err_High;
    }
    else if (variation.size() == 2) {
      Err_High = std::max( std::max( (variation[0] - nominal) , (variation[1] - nominal) ) , 0.0 );
      Err_Low  = std::max( std::max( (nominal - variation[0]) , (nominal - variation[1]) ) , 0.0 );
    }
  }
  else if (type.find("MC_Syst_PDF")!=std::string::npos) {
    // Variations (Use the offical EPPS16 approach)
    for(uint i = 0; i < (variation.size()/2); i++) {
      const double ctVal = nominal;
      const double miVal = variation[(2*i)+0];
      const double plVal = variation[(2*i)+1];
      Err_High += std::pow( std::max( std::max( (plVal - ctVal) , (miVal - ctVal) ) , 0.0 ) , 2.0 );
      Err_Low  += std::pow( std::min( std::min( (plVal - ctVal) , (miVal - ctVal) ) , 0.0 ) , 2.0 );
    }
    // Convert from 90% CL to 68% CL
    const double convFactor = TMath::ErfcInverse((1.00-0.68))/TMath::ErfcInverse((1.00-0.90));
    Err_High = ( convFactor * std::sqrt( Err_High ) );
    Err_Low  = ( convFactor * std::sqrt( Err_Low  ) );
  }
  else if (type.find("MC_Syst_Scale")!=std::string::npos) {
    // Variations (Use the envelope approach)
    for(uint i = 0; i < variation.size(); i++) {
      const double ctVal = nominal;
      const double vrVal = variation[i];
      Err_High = std::max( std::max( (vrVal - ctVal) , Err_High ) , 0.0 );
      Err_Low  = std::min( std::min( (vrVal - ctVal) , Err_Low  ) , 0.0 );
    }
    Err_High = std::abs( Err_High );
    Err_Low  = std::abs( Err_Low  );
  }
  else if (type.find("MC_Syst_Alpha")!=std::string::npos) {
    // Variations (Use the PDF4LHC15 approach)
    const double miVal = variation[0];
    const double plVal = variation[1];
    double err = ( ( plVal - miVal ) / 2.0 );
    // Convert from deltaAlpha_s = 0.001 to 0.0015 (68% CL)
    const double convFactor = (0.0015/0.0010);
    Err_High = ( convFactor * std::abs( err ) );
    Err_Low  = Err_High;
  }
  else {
    std::cout << "[ERROR] The systematic uncertainty " << type << " has not been implemented!" << std::endl; return;
  }
  // Let's symmetryze the errors, more conservative
  if (Err_High > Err_Low) { Err_Low = Err_High; } else { Err_High = Err_Low; }
  //
  Error_High = sumErrors( Error_High , Err_High );
  Error_Low  = sumErrors( Error_Low  , Err_Low  );
};


void getEffContent( double& val , double& err_High , double& err_Low , const TEfficiency& eff , const double& binVal )
{
  const int iBin = eff.FindFixBin(binVal);
  val      = eff.GetEfficiency(iBin);
  err_High = sumErrors( err_High , eff.GetEfficiencyErrorUp(iBin)  );
  err_Low  = sumErrors( err_Low  , eff.GetEfficiencyErrorLow(iBin) );
};


void getUncContent( double& err_High , double& err_Low , const TVectorD& unc , const TEfficiency& eff , const double& binVal )
{
  const int idx = ( eff.FindFixBin(binVal) - 1 ); // index is bin number - 1
  const double unc_High = unc[idx];
  const double unc_Low  = unc[idx];
  err_High = sumErrors( err_High , unc_High );
  err_Low  = sumErrors( err_Low  , unc_Low  );
};


bool getAcceptanceAndEfficiency( BinPentaMap& effMap , const std::string& inputFilePath , const bool useEtaCM = true , const bool isNominal = true )
{
  //
  // Update the TnP Type container later used to derive the efficiency variations
  if (isNominal) { fillTnPType(effTnPType_, absEtaTnP_, etaTnP_); }
  //
  // Open the input file
  TFile inputFile(inputFilePath.c_str(), "READ");
  if (inputFile.IsOpen()==false || inputFile.IsZombie()==true) { std::cout << "[ERROR] The input efficiency file " << inputFilePath << " was not found!" << std::endl; return false; }
  inputFile.cd();
  //
  std::string var = "Eta"; if (useEtaCM) { var = "EtaCM"; }
  const std::string sample = "MC_WToMuNu";
  const std::string corr = ( (inputFilePath.find("_With")!=std::string::npos) ? "HFCorr" : "NoCorr" );
  //
  for (const auto& c : effMap) {
    for (const auto& ch : c.second) {
      if (ch.second.count("Acceptance_MC")==0 || ch.second.at("Acceptance_MC").count("Val")==0) { std::cout << "[ERROR] The efficiency container is not valid" << std::endl; return false; }
      auto& eff = effMap.at(c.first).at(ch.first);
      //
      const std::string col = c.first;
      std::string charge = ""; if (ch.first=="Pl") { charge = "Plus"; } if (ch.first=="Mi") { charge = "Minus"; }
      //
      std::map< std::string , std::string > effFileDir;
      effFileDir["Acceptance_MC" ] = Form("TnPEfficiency1D/%s/%s/%s/Acceptance/%s"     , var.c_str(), sample.c_str(), col.c_str(), corr.c_str());
      effFileDir["Efficiency_MC" ] = Form("TnPEfficiency1D/%s/%s/%s/Total/%s"          , var.c_str(), sample.c_str(), col.c_str(), corr.c_str());
      effFileDir["Efficiency_TnP"] = Form("TnPEfficiency1D/%s/%s/%s/Total/TnP_Nominal" , var.c_str(), sample.c_str(), col.c_str());
      //
      std::map< std::string , std::vector< std::string > > effFileName;
      effFileName["Acceptance_MC" ].push_back( Form("eff1D_%s_%s_%s_%s_Acceptance_%s"     , var.c_str(), sample.c_str(), col.c_str(), charge.c_str(), corr.c_str()) );
      effFileName["Efficiency_MC" ].push_back( Form("eff1D_%s_%s_%s_%s_Total_%s"          , var.c_str(), sample.c_str(), col.c_str(), charge.c_str(), corr.c_str()) );
      effFileName["Efficiency_TnP"].push_back( Form("eff1D_%s_%s_%s_%s_Total_TnP_Nominal" , var.c_str(), sample.c_str(), col.c_str(), charge.c_str()) );
      //
      std::map< std::string , std::string > uncTnPDir , uncTnPName;
      if (isNominal) {
        uncTnPDir ["Efficiency_TnP_Stat"] = Form("TnPEfficiency1D/%s/%s/%s/Total"   , var.c_str(), sample.c_str(), col.c_str());
        uncTnPDir ["Efficiency_TnP_Syst"] = Form("TnPEfficiency1D/%s/%s/%s/Total"   , var.c_str(), sample.c_str(), col.c_str());
        uncTnPName["Efficiency_TnP_Stat"] = Form("unc1D_%s_%s_%s_%s_Total_TnP_Stat" , var.c_str(), sample.c_str(), col.c_str(), charge.c_str());
        uncTnPName["Efficiency_TnP_Syst"] = Form("unc1D_%s_%s_%s_%s_Total_TnP_Syst" , var.c_str(), sample.c_str(), col.c_str(), charge.c_str());
        uncTnPDir ["Efficiency_MC_Syst" ] = Form("TnPEfficiency1D/%s/%s/%s/Total"   , var.c_str(), sample.c_str(), col.c_str());
        uncTnPName["Efficiency_MC_Syst" ] = Form("unc1D_%s_%s_%s_%s_Total_MC_Syst"  , var.c_str(), sample.c_str(), col.c_str(), charge.c_str());
      }
      //
      // Only extract the components for the uncertainties if we are in the nominal case (systematics only care about the value)
      if (isNominal) {
        // Used to derive the TnP and PDF variations at observable level
        for (const auto& effT : effTnPType_) {
          if ( (effT.first.find("_STA")!=std::string::npos) || (effT.first.find("_PU")!=std::string::npos) ) continue;
          const std::string name = Form("Efficiency_%s", effT.first.c_str());
          effFileDir[name] = Form("TnPEfficiency1D/%s/%s/%s/Total/%s" , var.c_str(), sample.c_str(), col.c_str(), effT.first.c_str());
          for (uint i = 0; i < effT.second.first; i++) {
            std::string effName = Form("eff1D_%s_%s_%s_%s_Total_%s" , var.c_str(), sample.c_str(), col.c_str(), charge.c_str(), effT.first.c_str());
            if (effT.second.first>1) { effName = Form("%s_%d" , effName.c_str(), i); }
            effFileName[name].push_back( effName );
          }
        }
        // Used to propagate the TnP and PDF uncertainties to the observable
        for (const auto& effT : effTnPType_) {
          if ( (effT.first.find("_STA")!=std::string::npos) || (effT.first.find("_PU")!=std::string::npos) ) continue;
          const std::string name = Form("Efficiency_%s", effT.first.c_str());
          uncTnPDir[name]  = Form("TnPEfficiency1D/%s/%s/%s/Total/%s" , var.c_str(), sample.c_str(), col.c_str(), effT.first.c_str());
          uncTnPName[name] = Form("unc1D_%s_%s_%s_%s_Total_%s" , var.c_str(), sample.c_str(), col.c_str(), charge.c_str(), effT.first.c_str());
        }
      }
      //
      // Procced to extract the efficiencies
      std::map< std::string , TEfficiency* > effFileObj;
      for (const auto& effT : effFileName) {
        for (uint i = 0; i < effT.second.size(); i++) {
          const std::string name = effT.first + ( (effT.second.size()>1) ? Form("_%d", i) : "");
          effFileObj[name] = (TEfficiency*) inputFile.Get(Form("%s/%s", effFileDir.at(effT.first).c_str(), effT.second[i].c_str()));
          if (effFileObj.at(name)==NULL) {
            if (name.find("MC_Syst")!=std::string::npos) {
              std::cout << "[INFO] " << effT.second[i] << " located in " << effFileDir.at(effT.first) << " was not found in file " << inputFilePath << ". Will ignore it!" << std::endl;
              effFileObj.erase(name);
            }
            else { std::cout << "[ERROR] " << effT.second[i] << " located in " << effFileDir.at(effT.first) << " was not found in file " << inputFilePath << std::endl; return false; }
          }
        }
      }
      // Proceed to extract the TnP and PDF uncertainties
      std::map< std::string , TVectorD* > uncTnPObj;
      for (const auto& uncT : uncTnPName) {
        uncTnPObj[uncT.first] = (TVectorD*) inputFile.Get(Form("%s/%s", uncTnPDir.at(uncT.first).c_str(), uncT.second.c_str()));
        if (uncTnPObj.at(uncT.first)==NULL) {
          if (uncT.first.find("MC_Syst")!=std::string::npos) {
            std::cout << "[INFO] " << uncT.second << " located in " << uncTnPDir.at(uncT.first) << " was not found in file " << inputFilePath << ". Will ignore it!" << std::endl;
            uncTnPObj.erase(uncT.first);
          }
          else { std::cout << "[ERROR] " << uncT.second << " located in " << uncTnPDir.at(uncT.first) << " was not found in file " << inputFilePath << std::endl; return false; }
        }
      }
      //
      for (const auto& b : ch.second.at("Acceptance_MC").at("Val")) {
        // Compute the center value of the bin
        const double etaVal = ( ( b.first.etabin().high() + b.first.etabin().low() ) / 2.0 );
        // Fill the Efficiency Objects
        for (const auto& effT : effFileObj) {
          getEffContent(eff[effT.first]["Val"][b.first] , eff[effT.first]["Err_Stat_High"][b.first] , eff[effT.first]["Err_Stat_Low"][b.first] , *(effT.second) , etaVal);
          eff[effT.first]["Err_Syst_High"][b.first] = 0.0; eff[effT.first]["Err_Syst_Low"][b.first] = 0.0;
          if ( (isNominal==false) || (effT.first.find("TnP_")!=std::string::npos) || (effT.first.find("MC_Syst")!=std::string::npos) ){
            eff[effT.first]["Err_Stat_High"][b.first] = 0.0; eff[effT.first]["Err_Stat_Low"][b.first] = 0.0;
          }
        }
        //
        // Only compute uncertainties for Nominal case , for systematics we only care about their value
        if (isNominal) {
          // Set the effiency values for the PileUp and Event Activity
          std::vector< std::string > effLbl = { "PU_0" , "PU_1" , "STA_0" , "STA_1" , "PU_PROP" , "STA_PROP" };
          std::vector< std::string > parLbl = { "Val" , "Err_Syst_High" , "Err_Syst_Low" , "Err_Stat_High" , "Err_Stat_Low" };
          for (const auto& e : effLbl) {
            for (const auto& p : parLbl) {
              const std::string name = ("Efficiency_TnP_Syst_" + e);
              //
              double uncVal = 0.;
              if (e.find("PU" )!=std::string::npos) { uncVal = 0.0034; } // Impact of PileUp and Event Activity on Isolation +- 0.0034 (Relative uncertainty)
              if (e.find("STA")!=std::string::npos) { uncVal = 0.0060; } // STA Efficiency mismodelling +-0.0060 (Relative uncertainty)
              //
              if ( (e.find("_PROP")!=std::string::npos) && (p.find("Err_Syst")!=std::string::npos) ) {
                eff[name][p][b.first] = uncVal*eff.at("Efficiency_TnP").at("Val").at(b.first);
              }
              else if ( (e.find("_0")!=std::string::npos) && (p=="Val") ) {
                eff[name][p][b.first] = eff.at("Efficiency_TnP").at("Val").at(b.first)*(1.0 + uncVal);
              }
              else if ( (e.find("_1")!=std::string::npos) && (p=="Val") ) {
                eff[name][p][b.first] = eff.at("Efficiency_TnP").at("Val").at(b.first)*(1.0 - uncVal);
              }
              else if (p=="Val") {
                eff[name][p][b.first] = eff.at("Efficiency_TnP").at("Val").at(b.first);
              }
              else {
                eff[name][p][b.first] = 0.0;
              }
            }
          }
          // Fill TnP and PDF Uncertainties
          for (const auto& uncT : uncTnPObj) {
            const std::string name = Form("%s_PROP", uncT.first.c_str());
            getUncContent(eff[name]["Err_Syst_High"][b.first] , eff[name]["Err_Syst_Low"][b.first] , *(uncT.second) , *(effFileObj.at("Efficiency_TnP")) , etaVal);
            eff[name]["Err_Stat_High"][b.first] = 0.0; eff[name]["Err_Stat_Low"][b.first] = 0.0;
            eff[name]["Val"][b.first] = eff.at("Efficiency_TnP").at("Val").at(b.first);
            // Add the TnP Stat and Syst, and the MC PDF Syst to the corrected efficiency systematic error
            if (uncT.first=="Efficiency_TnP_Stat" || uncT.first=="Efficiency_TnP_Syst" || uncT.first=="Efficiency_MC_Syst") {
              getUncContent(eff.at("Efficiency_TnP").at("Err_Syst_High").at(b.first) , eff.at("Efficiency_TnP").at("Err_Syst_Low").at(b.first) , *(uncT.second) , *(effFileObj.at("Efficiency_TnP")) , etaVal);
            }
          }
        }
        //
        // Fill the Total Uncertainties
        for (auto& e : eff) {
          e.second["Err_Tot_Low" ][b.first] = sumErrors(e.second["Err_Stat_Low" ][b.first], e.second["Err_Syst_Low" ][b.first]);
          e.second["Err_Tot_High"][b.first] = sumErrors(e.second["Err_Stat_High"][b.first], e.second["Err_Syst_High"][b.first]);
        }
      }
    }
  }
  // Close Input File
  inputFile.Close();
  // Return
  return true;
};

  
double getCorrectedYieldValue( const double& N_Raw , const double& Acceptance , const double& Efficiency )
{
  return ( N_Raw / ( Acceptance * Efficiency ) );
};

  
double getCorrectedYieldError( const double& N_Raw , const double& Acceptance , const double& Efficiency , const double& Err_N_Raw , const double& Err_Acceptance , const double& Err_Efficiency )
{
  const double common    = getCorrectedYieldValue( N_Raw , Acceptance , Efficiency );
  const double error2N   = std::pow( ( Err_N_Raw / N_Raw ) , 2.0 );
  const double error2Eff = std::pow( ( Err_Efficiency / Efficiency ) , 2.0 );
  const double error2Acc = std::pow( ( Err_Acceptance / Acceptance ) , 2.0 );
  return ( std::abs(common) * std::sqrt( error2N + error2Eff + error2Acc ) );
};


bool correctRawYields( VarBinMap& inputVar , const BinPentaMap& effMap , const std::string accType = "MC" , const std::string effType = "TnP" , const bool isNominal = true )
{
  //
  const std::string accName = ( (accType!="") ? Form("Acceptance_%s", accType.c_str()) : "" );
  const std::string effName = ( (effType!="") ? Form("Efficiency_%s", effType.c_str()) : "" );
  //
  for (const auto& c : inputVar) {
    for (const auto& b : c.second) {
      for (const auto& ch : b.second) {
        if (effMap.count(c.first)==0) { std::cout << "[ERROR] Efficiency container does not have the collision " << c.first << std::endl; return false; }
        if (effMap.at(c.first).count(ch.first)==0) { std::cout << "[ERROR] Efficiency container does not have the charge " << ch.first << std::endl; return false; }
        // Check that eveythin is fine
        const auto& eff = effMap.at(c.first).at(ch.first);
        if (accName!="" && eff.count(accName)==0) { std::cout << "[ERROR] Efficiency container does not have the variable " << accName << std::endl; return false; }
        if (effName!="" && eff.count(effName)==0) { std::cout << "[ERROR] Efficiency container does not have the variable " << effName << std::endl; return false; }
        if ( (accName!="" && eff.at(accName).at("Val").count(b.first)==0) || (effName!="" && eff.at(effName).at("Val").count(b.first)==0) ) {
          std::cout << "[ERROR] Efficiency container does not have the bin [" << b.first.etabin().high() << " , " << b.first.etabin().low() << "] " << std::endl; return false;
        }
        // Extract the Acceptance and Efficiency
        std::map< std::string , double > Acceptance , Efficiency;
        if (accName!="") { for (const auto& t : eff.at(accName)) { Acceptance[t.first] = t.second.at(b.first); } }
        else {
          Acceptance["Val"] = 1.0; Acceptance["Err_Tot_Low"] = 0.0; Acceptance["Err_Tot_High"] = 0.0;
          for (const auto& t : ch.second.at("N_WToMu")) { if (t.first=="Val") continue; Acceptance[t.first] = 0.0; }
        }
        if (effName!="") { for (const auto& t : eff.at(effName)) { Efficiency[t.first] = t.second.at(b.first); } }
        else {
          Efficiency["Val"] = 1.0; Efficiency["Err_Tot_Low"] = 0.0; Efficiency["Err_Tot_High"] = 0.0;
          for (const auto& t : ch.second.at("N_WToMu")) { if (t.first=="Val") continue; Efficiency[t.first] = 0.0; }
        }
        //
        // Create a backup
        auto& N_Co = inputVar.at(c.first).at(b.first).at(ch.first);
        for (const auto& t : N_Co.at("N_WToMu")) { N_Co["N_WToMu_RAW"][t.first] = t.second; }
        if (accName=="Acceptance_MC" ) { for (const auto& t : eff.at(accName)) { N_Co[accName][t.first] = t.second.at(b.first); } }
        if (effName=="Efficiency_MC" || effName=="Efficiency_TnP") { for (const auto& t : eff.at(effName)) { N_Co[effName][t.first] = t.second.at(b.first); } }
        //
        auto&       N_Corr = N_Co.at("N_WToMu");
        const auto& N_Raw  = N_Co.at("N_WToMu_RAW");
        //
        // Fill with the corrected values;
        N_Corr.at("Val") = getCorrectedYieldValue(N_Raw.at("Val") , Acceptance.at("Val") , Efficiency.at("Val"));
        N_Corr.at("Err_Stat_Low" ) = getCorrectedYieldError(N_Raw.at("Val") , Acceptance.at("Val") , Efficiency.at("Val") , N_Raw.at("Err_Stat_Low" ) , 0.0 , 0.0);
        N_Corr.at("Err_Stat_High") = getCorrectedYieldError(N_Raw.at("Val") , Acceptance.at("Val") , Efficiency.at("Val") , N_Raw.at("Err_Stat_High") , 0.0 , 0.0);
        N_Corr.at("Err_Syst_Low" ) = getCorrectedYieldError(N_Raw.at("Val") , Acceptance.at("Val") , Efficiency.at("Val") , N_Raw.at("Err_Syst_Low" ) , Acceptance.at("Err_Tot_Low" ) , Efficiency.at("Err_Tot_Low" ));
        N_Corr.at("Err_Syst_High") = getCorrectedYieldError(N_Raw.at("Val") , Acceptance.at("Val") , Efficiency.at("Val") , N_Raw.at("Err_Syst_High") , Acceptance.at("Err_Tot_High") , Efficiency.at("Err_Tot_High"));
        if (isNominal == false) { for (auto& t : N_Corr) { if (t.first!="Val") { t.second = 0.0; } } }
        //
        if (isNominal) {
          for (const auto& effT : eff) {
            if (
                (accType == "MC"  && effT.first=="Acceptance_MC") ||
                (effType == "MC"  && effT.first=="Efficiency_MC") ||
                (effType == "TnP" && effT.first.find("Efficiency")!=std::string::npos)
              ) {
              auto& N_Corr = inputVar.at(c.first).at(b.first).at(ch.first)[Form("N_WToMu_%s", effT.first.c_str())];
              N_Corr["Val"] = getCorrectedYieldValue(N_Raw.at("Val") , Acceptance.at("Val") , eff.at(effT.first).at("Val").at(b.first));
              N_Corr["Err_Stat_Low" ] = getCorrectedYieldError(N_Raw.at("Val") , Acceptance.at("Val") , eff.at(effT.first).at("Val").at(b.first) ,
                                                               N_Raw.at("Err_Syst_Low" ) , Acceptance.at("Err_Stat_Low" ) , eff.at(effT.first).at("Err_Stat_Low").at(b.first));
              N_Corr["Err_Stat_High"] = getCorrectedYieldError(N_Raw.at("Val") , Acceptance.at("Val") , eff.at(effT.first).at("Val").at(b.first) ,
                                                               N_Raw.at("Err_Syst_High") , Acceptance.at("Err_Stat_High") , eff.at(effT.first).at("Err_Stat_High").at(b.first));
              N_Corr["Err_Syst_Low" ] = getCorrectedYieldError(N_Raw.at("Val") , Acceptance.at("Val") , eff.at(effT.first).at("Val").at(b.first) ,
                                                               N_Raw.at("Err_Syst_Low" ) , Acceptance.at("Err_Syst_Low" ) , eff.at(effT.first).at("Err_Syst_Low").at(b.first));
              N_Corr["Err_Syst_High"] = getCorrectedYieldError(N_Raw.at("Val") , Acceptance.at("Val") , eff.at(effT.first).at("Val").at(b.first) ,
                                                               N_Raw.at("Err_Syst_High") , Acceptance.at("Err_Syst_High") , eff.at(effT.first).at("Err_Syst_High").at(b.first));
            }
          }
        }
      }
    }
  }
  return true;
};


double getChargeAsymmetryValue( const double& N_Plus , const double& N_Minus )
{
  const double diff = ( N_Plus - N_Minus );
  const double sum  = ( N_Plus + N_Minus );
  return ( diff / sum );
};

  
double getChargeAsymmetryError( const double& N_Plus , const double& N_Minus , const double& Err_Plus , const double& Err_Minus )
{
  const double sum      = ( N_Plus + N_Minus );
  const double common   = ( (2.0 * N_Minus * N_Plus) / (sum * sum) );
  const double error2Pl = std::pow( ( Err_Plus  / N_Plus  ) , 2.0 );
  const double error2Mi = std::pow( ( Err_Minus / N_Minus ) , 2.0 );
  return ( std::abs(common) * std::sqrt( error2Pl + error2Mi ) );
};


bool computeChargeAsymmetry( BinPentaMap& var , BinSextaMapVec& systVar , const VarBinMap& inputVar , const bool& doSyst , const VarBinMap& nomVar , const uint& corrSyst )
{
  //
  for (const auto& c : inputVar) {
    for (const auto& b : c.second) {
      // Check that everything is fine
      if (b.second.count("Pl")==0) { std::cout << "[ERROR] Plus charge in " << c.first << " is missing in bin ["  << b.first.etabin().low() << " , " << b.first.etabin().high() << "]" << std::endl; return false; }
      if (b.second.count("Mi")==0) { std::cout << "[ERROR] Minus charge in " << c.first << " is missing in bin [" << b.first.etabin().low() << " , " << b.first.etabin().high() << "]" << std::endl; return false; }
      if (b.second.at("Mi").count("N_WToMu")==0 || b.second.at("Pl").count("N_WToMu")==0) {
        std::cout << "[ERROR] N_WToMu variable is missing in bin [" << b.first.etabin().high() << " , " << b.first.etabin().low() << "]" << std::endl; return false;
      }
      //
      auto&       oVar    = var[c.first][""]["Charge_Asymmetry"];
      const auto& iVar_Pl = b.second.at("Pl").at("N_WToMu");
      const auto& iVar_Mi = b.second.at("Mi").at("N_WToMu");
      //
      // Compute charge asymmetry value
      oVar["Val"][b.first] = getChargeAsymmetryValue(iVar_Pl.at("Val"), iVar_Mi.at("Val"));
      //
      // Compute charge asymmetry statistical error
      const std::vector< std::string > statT = { "Err_Stat_High" , "Err_Stat_Low" };
      for (const auto& t : statT) { oVar[t][b.first] = getChargeAsymmetryError(iVar_Pl.at("Val"), iVar_Mi.at("Val"), iVar_Pl.at(t), iVar_Mi.at(t)); }
      //
      // Compute the charge asymmetry systematic error (in this case from efficiency TnP)
      const std::vector< std::string > systT = { "Err_Syst_High" , "Err_Syst_Low" };
      // Initialize the Systematic Uncertainties
      for (const auto& t : systT) { oVar[t][b.first] = 0.0; }
      //
      if (doSyst && nomVar.size()==0) {
        //
        std::string effType = ""; if (b.second.at("Pl").count("N_WToMu_Efficiency_TnP")>0) { effType = "TnP"; } else if (b.second.at("Pl").count("N_WToMu_Efficiency_MC")>0) { effType = "MC"; }
        //
        // If we are not using TnP efficieny, simply progate the uncertainties
        if (effType!="TnP") {
          for (const auto& t : systT) { oVar[t][b.first] = getChargeAsymmetryError(iVar_Pl.at("Val"), iVar_Mi.at("Val"), iVar_Pl.at(t), iVar_Mi.at(t)); }
        }
        else {
          // Add the Statistical Error of the Efficiency
          const auto& iMCVar_Pl = b.second.at("Pl").at(Form("N_WToMu_Efficiency_%s", effType.c_str()));
          const auto& iMCVar_Mi = b.second.at("Mi").at(Form("N_WToMu_Efficiency_%s", effType.c_str()));
          //
          if (systVar["MC_Statistics"].size()<1) { systVar.at("MC_Statistics").resize(1); }
          auto& sMCVar = systVar.at("MC_Statistics")[0][c.first][""]["Charge_Asymmetry"];
          sMCVar["Val"][b.first] = oVar.at("Val").at(b.first);
          for (uint i=0; i<systT.size(); i++) {
            const double error = getChargeAsymmetryError(iMCVar_Pl.at("Val"), iMCVar_Mi.at("Val"), iMCVar_Pl.at(statT[i]), iMCVar_Mi.at(statT[i]));
            oVar.at(systT[i]).at(b.first) = sumErrors( oVar.at(systT[i]).at(b.first) , error );
            sMCVar[systT[i]][b.first] = error;
          }
          // Add the other Systematic Errors of the Efficiency
          DoubleVecMap variation;
          for (const auto& effT : effTnPType_) {
            //
            // Case: Apply Propagation Method
            if (effT.second.second<2) {
              const std::string vL = Form("N_WToMu_Efficiency_%s_PROP", effT.first.c_str());
              if (b.second.at("Mi").count(vL)==0 || b.second.at("Pl").count(vL)==0) {
                std::cout << "[ERROR] " << vL << " variable is missing in bin [" << b.first.etabin().high() << " , " << b.first.etabin().low() << "]" << std::endl; return false;
              }
              if (systVar[effT.first].size()<1) { systVar.at(effT.first).resize(1); }
              auto& sVar = systVar.at(effT.first)[0][c.first][""]["Charge_Asymmetry"];
              sVar["Val"][b.first] = oVar.at("Val").at(b.first);
              for (const auto& t : systT) {
                // Case: Charge UnCorrelated
                if (effT.second.second==0 || effT.second.second==1) {
                  const double error = getChargeAsymmetryError(b.second.at("Pl").at(vL).at("Val"), b.second.at("Mi").at(vL).at("Val"),
                                                               b.second.at("Pl").at(vL).at(t)    , b.second.at("Mi").at(vL).at(t));
                  oVar.at(t).at(b.first) = sumErrors( oVar.at(t).at(b.first) , error );
                  sVar[t][b.first] = error;
                }
              }
            }
            // Case: Apply Variation Method
            if (effT.second.second>=2) {
              if (systVar[effT.first].size()<effT.second.first) { systVar.at(effT.first).resize(effT.second.first); }
              for (uint i = 0; i < effT.second.first; i++) {
                const std::string vL = std::string("N_WToMu_Efficiency_") + effT.first + ( (effT.second.first>1) ? Form("_%d", i) : "" );
                if (b.second.at("Mi").count(vL)==0 || b.second.at("Pl").count(vL)==0) {
                  std::cout << "[ERROR] " << vL << " variable is missing in bin [" << b.first.etabin().high() << " , " << b.first.etabin().low() << "]" << std::endl; return false;
                }
                auto& sVar = systVar.at(effT.first)[i][c.first][""]["Charge_Asymmetry"];
                for (const auto& t : systT) { sVar[t][b.first] = 0.0; }
                sVar["Val"][b.first] = oVar.at("Val").at(b.first);
                // Case: Charge Correlated
                if (effT.second.second==2 || effT.second.second==3 || effT.second.second==5) {
                  const double vVar = getChargeAsymmetryValue(b.second.at("Pl").at(vL).at("Val"), b.second.at("Mi").at(vL).at("Val"));
                  variation[effT.first].push_back( vVar );
                  for (const auto& t : systT) { sVar[t][b.first] = std::abs(vVar - oVar.at("Val").at(b.first)); }
                  sVar["Val"][b.first] = vVar;
                }
                // Case: Charge Uncorrelated
                if (effT.second.second==4) {
                  const double plVar = getChargeAsymmetryValue(b.second.at("Pl").at(vL).at("Val") ,                  iVar_Mi.at("Val"));
                  const double miVar = getChargeAsymmetryValue(                 iVar_Pl.at("Val") , b.second.at("Mi").at(vL).at("Val"));
                  variation[(effT.first+"_Pl")].push_back( plVar );
                  variation[(effT.first+"_Mi")].push_back( miVar );
                  for (const auto& t : systT) { systVar.at(effT.first)[i][c.first][""]["Charge_Asymmetry_Pl"][t][b.first] = std::abs(plVar - oVar.at("Val").at(b.first)); }
                  for (const auto& t : systT) { systVar.at(effT.first)[i][c.first][""]["Charge_Asymmetry_Mi"][t][b.first] = std::abs(miVar - oVar.at("Val").at(b.first)); }
                }
              }
            }
          }
          for (const auto& vEr : variation) {
            computeTnPUncertainty(oVar.at("Err_Syst_High").at(b.first), oVar.at("Err_Syst_Low").at(b.first), variation.at(vEr.first), oVar.at("Val").at(b.first), vEr.first);
          }
        }
      }
      else if (doSyst && nomVar.size()>0) {
	//
	const auto& nVar_Pl = nomVar.at(c.first).at(b.first).at("Pl").at("N_WToMu");
	const auto& nVar_Mi = nomVar.at(c.first).at(b.first).at("Mi").at("N_WToMu");
	//
	// Add the Systematic Error
	//
	// Case: Apply Variation Method
        //
	const double nomVal = getChargeAsymmetryValue(nVar_Pl.at("Val"), nVar_Mi.at("Val"));
	oVar.at("Val").at(b.first) = nomVal;
        //
        // Case: Charge Correlated
        if (corrSyst==2 || corrSyst==3 || corrSyst==5) {
          const double varVal = getChargeAsymmetryValue(iVar_Pl.at("Val"), iVar_Mi.at("Val"));
          for (const auto& t : systT) { oVar.at(t).at(b.first) = std::abs(varVal - nomVal); }
          oVar.at("Val").at(b.first) = varVal;
        }
        // Case: Charge Uncorrelated
        if (corrSyst==0 || corrSyst==1 || corrSyst==4) {
          const double plVar = getChargeAsymmetryValue(iVar_Pl.at("Val"), nVar_Mi.at("Val"));
          const double miVar = getChargeAsymmetryValue(nVar_Pl.at("Val"), iVar_Mi.at("Val"));
          for (const auto& t : systT) { var.at(c.first).at("")["Charge_Asymmetry_Pl"][t][b.first] = std::abs(plVar - nomVal); }
          for (const auto& t : systT) { var.at(c.first).at("")["Charge_Asymmetry_Mi"][t][b.first] = std::abs(miVar - nomVal); }
        }
      }
    }
  }
  return true;
};

  
double getForwardBackwardRatioValue( const double& N_Forward , const double& N_Backward )
{
  return ( N_Forward / N_Backward );
};

  
double getForwardBackwardRatioError( const double& N_Forward , const double& N_Backward , const double& Err_Forward , const double& Err_Backward )
{
  const double common   = getForwardBackwardRatioValue( N_Forward , N_Backward );
  const double error2Fw = std::pow( ( Err_Forward  / N_Forward  ) , 2.0 );
  const double error2Bw = std::pow( ( Err_Backward / N_Backward ) , 2.0 );
  return ( std::abs(common) * std::sqrt( error2Fw + error2Bw ) );
};

  
double getForwardBackwardRatioError( const double& N_Forward_Plus   , const double& N_Backward_Plus   , const double& N_Forward_Minus   , const double& N_Backward_Minus ,
                                     const double& Err_Forward_Plus , const double& Err_Backward_Plus , const double& Err_Forward_Minus , const double& Err_Backward_Minus )
{
  const double common     = getForwardBackwardRatioValue( ( N_Forward_Plus + N_Forward_Minus ) , ( N_Backward_Plus + N_Backward_Minus ) );
  const double error2FwPl = std::pow( ( Err_Forward_Plus   / ( N_Forward_Plus  + N_Forward_Minus  ) ) , 2.0 );
  const double error2FwMi = std::pow( ( Err_Forward_Minus  / ( N_Forward_Plus  + N_Forward_Minus  ) ) , 2.0 );
  const double error2BwPl = std::pow( ( Err_Backward_Plus  / ( N_Backward_Plus + N_Backward_Minus ) ) , 2.0 );
  const double error2BwMi = std::pow( ( Err_Backward_Minus / ( N_Backward_Plus + N_Backward_Minus ) ) , 2.0 );
  return ( std::abs(common) * std::sqrt( error2FwPl + error2FwMi + error2BwPl + error2BwMi ) );
};


bool computeForwardBackwardRatio( BinPentaMap& var ,BinSextaMapVec& systVar , const VarBinMap& inputVar , const bool& doSyst , const VarBinMap& nomVar , const uint& corrSyst )
{
  //
  for (const auto& c : inputVar) {
    for (const auto& b : c.second) {
      if (b.first.etabin().low() < 0.0 ) continue; // Ignore the backward bins when looping
      // Build the backward bin
      const auto& binFw = b.first;
      const auto  binBw = anabin<0>(-1.0*binFw.etabin().high() , ( (binFw.etabin().low() == 0.0) ? 0.0 : -1.0*binFw.etabin().low() ));
      // Check that everything is fine
      if (c.second.count(binBw)==0) { /*std::cout << "[WARNING] Backward bin [" << binBw.etabin().low() << " , " << binBw.etabin().high() << "] is not defined" << std::endl;*/ continue; }
      if (c.second.count(binFw)==0) { std::cout << "[ERROR] Forward bin [" << binBw.etabin().low() << " , " << binBw.etabin().high() << "] is not defined" << std::endl; return false; }
      if (c.second.at(binFw).count("Pl")==0) { std::cout << "[ERROR] Plus charge in " << c.first << " is missing in bin [" << binFw.etabin().low() << " , " << binFw.etabin().high() << "]" << std::endl; return false; }
      if (c.second.at(binFw).count("Mi")==0) { std::cout << "[ERROR] Minus charge in " << c.first << " is missing in bin [" << binFw.etabin().low() << " , " << binFw.etabin().high() << "]" << std::endl; return false; }
      if (c.second.at(binBw).count("Pl")==0) { std::cout << "[ERROR] Plus charge in " << c.first << " is missing in bin [" << binBw.etabin().low() << " , " << binBw.etabin().high() << "]" << std::endl; return false; }
      if (c.second.at(binBw).count("Mi")==0) { std::cout << "[ERROR] Minus charge in " << c.first << " is missing in bin [" << binBw.etabin().low() << " , " << binBw.etabin().high() << "]" << std::endl; return false; }
      if (c.second.at(binFw).at("Mi").count("N_WToMu")==0 || c.second.at(binFw).at("Pl").count("N_WToMu")==0) {
        std::cout << "[ERROR] N_WToMu variable is missing in bin [" << binFw.etabin().low() << " , " << binFw.etabin().high() << "]" << std::endl; return false;
      }
      if (c.second.at(binBw).at("Mi").count("N_WToMu")==0 || c.second.at(binBw).at("Pl").count("N_WToMu")==0) {
        std::cout << "[ERROR] N_WToMu variable is missing in bin [" << binBw.etabin().low() << " , " << binBw.etabin().high() << "]" << std::endl; return false;
      }
      //
      auto&       oVar_Inc  = var[c.first][""]["ForwardBackward_Ratio"];
      auto&       oVar_Pl   = var[c.first]["Pl"]["ForwardBackward_Ratio"];
      auto&       oVar_Mi   = var[c.first]["Mi"]["ForwardBackward_Ratio"];
      const auto& iVar_FwPl = c.second.at(binFw).at("Pl").at("N_WToMu");
      const auto& iVar_FwMi = c.second.at(binFw).at("Mi").at("N_WToMu");
      const auto& iVar_BwPl = c.second.at(binBw).at("Pl").at("N_WToMu");
      const auto& iVar_BwMi = c.second.at(binBw).at("Mi").at("N_WToMu");
      //
      // Compute forward-backward ratio value
      oVar_Inc["Val"][b.first] = getForwardBackwardRatioValue( ( iVar_FwPl.at("Val") + iVar_FwMi.at("Val") ) , ( iVar_BwPl.at("Val") + iVar_BwMi.at("Val") ) );
      oVar_Pl["Val"][b.first]  = getForwardBackwardRatioValue( iVar_FwPl.at("Val") , iVar_BwPl.at("Val") );
      oVar_Mi["Val"][b.first]  = getForwardBackwardRatioValue( iVar_FwMi.at("Val") , iVar_BwMi.at("Val") );
      //
      // Compute forward-backward ratio statistical error
      const std::vector< std::string > statT = { "Err_Stat_High" , "Err_Stat_Low" };
      for (const auto& t : statT) {
        oVar_Inc[t][b.first] = getForwardBackwardRatioError(
                                                            iVar_FwPl.at("Val"), iVar_BwPl.at("Val"), iVar_FwMi.at("Val"), iVar_BwMi.at("Val"), 
                                                            iVar_FwPl.at(t)    , iVar_BwPl.at(t)    , iVar_FwMi.at(t)    , iVar_BwMi.at(t)
                                                            );
      }
      for (const auto& t : statT) { oVar_Pl[t][b.first] = getForwardBackwardRatioError(iVar_FwPl.at("Val"), iVar_BwPl.at("Val"), iVar_FwPl.at(t), iVar_BwPl.at(t)); }
      for (const auto& t : statT) { oVar_Mi[t][b.first] = getForwardBackwardRatioError(iVar_FwMi.at("Val"), iVar_BwMi.at("Val"), iVar_FwMi.at(t), iVar_BwMi.at(t)); }
      //
      // Compute the forward-backward ratio systematic error (in this case from efficiency TnP)
      const std::vector< std::string > systT = { "Err_Syst_High" , "Err_Syst_Low" };
      // Initialize the Systematic Uncertainties
      for (const auto& t : systT) { oVar_Inc[t][b.first] = 0.0; }
      for (const auto& t : systT) { oVar_Pl[t][b.first]  = 0.0; }
      for (const auto& t : systT) { oVar_Mi[t][b.first]  = 0.0; }
      //
      if (doSyst && nomVar.size()==0) {
        //
        std::string effType = ""; if (b.second.at("Pl").count("N_WToMu_Efficiency_TnP")>0) { effType = "TnP"; } else if (b.second.at("Pl").count("N_WToMu_Efficiency_MC")>0) { effType = "MC"; }
        //
        // If we are not using TnP efficieny, simply progate the uncertainties
        if (effType!="TnP") {
          for (const auto& t : systT) {
            oVar_Inc[t][b.first] = getForwardBackwardRatioError(
                                                                iVar_FwPl.at("Val"), iVar_BwPl.at("Val"), iVar_FwMi.at("Val"), iVar_BwMi.at("Val"), 
                                                                iVar_FwPl.at(t)    , iVar_BwPl.at(t)    , iVar_FwMi.at(t)    , iVar_BwMi.at(t)
                                                                );
          }
          for (const auto& t : systT) { oVar_Pl[t][b.first] = getForwardBackwardRatioError(iVar_FwPl.at("Val"), iVar_BwPl.at("Val"), iVar_FwPl.at(t), iVar_BwPl.at(t)); }
          for (const auto& t : systT) { oVar_Mi[t][b.first] = getForwardBackwardRatioError(iVar_FwMi.at("Val"), iVar_BwMi.at("Val"), iVar_FwMi.at(t), iVar_BwMi.at(t)); }
        }
        else {
          // Add the Statistical Error of the Efficiency
          const auto& iMCVar_FwPl = c.second.at(binFw).at("Pl").at(Form("N_WToMu_Efficiency_%s", effType.c_str()));
          const auto& iMCVar_FwMi = c.second.at(binFw).at("Mi").at(Form("N_WToMu_Efficiency_%s", effType.c_str()));
          const auto& iMCVar_BwPl = c.second.at(binBw).at("Pl").at(Form("N_WToMu_Efficiency_%s", effType.c_str()));
          const auto& iMCVar_BwMi = c.second.at(binBw).at("Mi").at(Form("N_WToMu_Efficiency_%s", effType.c_str()));
          if (systVar["MC_Statistics"].size()<1) { systVar.at("MC_Statistics").resize(1); }
          // Plus
          auto& sMCVar_Pl = systVar.at("MC_Statistics")[0][c.first]["Pl"]["ForwardBackward_Ratio"];
          sMCVar_Pl["Val"][b.first] = oVar_Pl.at("Val").at(b.first);
          for (uint i=0; i<systT.size(); i++) {
            const double error = getForwardBackwardRatioError(iMCVar_FwPl.at("Val"), iMCVar_BwPl.at("Val"), iMCVar_FwPl.at(statT[i]), iMCVar_BwPl.at(statT[i]));
            oVar_Pl.at(systT[i]).at(b.first) =  sumErrors( oVar_Pl.at(systT[i]).at(b.first) , error );
            sMCVar_Pl[systT[i]][b.first] = error;
          }
          // Minus
          auto& sMCVar_Mi = systVar.at("MC_Statistics")[0][c.first]["Mi"]["ForwardBackward_Ratio"];
          sMCVar_Mi["Val"][b.first] = oVar_Mi.at("Val").at(b.first);
          for (uint i=0; i<systT.size(); i++) {
            const double error = getForwardBackwardRatioError(iMCVar_FwMi.at("Val"), iMCVar_BwMi.at("Val"), iMCVar_FwMi.at(statT[i]), iMCVar_BwMi.at(statT[i]));
            oVar_Mi.at(systT[i]).at(b.first) =  sumErrors( oVar_Mi.at(systT[i]).at(b.first) , error );
            sMCVar_Mi[systT[i]][b.first] = error;
          }
          // Inclusive
          auto& sMCVar_Inc = systVar.at("MC_Statistics")[0][c.first][""]["ForwardBackward_Ratio"];
          sMCVar_Inc["Val"][b.first] = oVar_Inc.at("Val").at(b.first);
          for (uint i=0; i<systT.size(); i++) {
            const double error = getForwardBackwardRatioError( iMCVar_FwPl.at("Val")    , iMCVar_BwPl.at("Val")    , iMCVar_FwMi.at("Val")    , iMCVar_BwMi.at("Val"), 
                                                               iMCVar_FwPl.at(statT[i]) , iMCVar_BwPl.at(statT[i]) , iMCVar_FwMi.at(statT[i]) , iMCVar_BwMi.at(statT[i]) );
            oVar_Inc.at(systT[i]).at(b.first) = sumErrors( oVar_Inc.at(systT[i]).at(b.first) , error );
            sMCVar_Inc[systT[i]][b.first] = error;
          }
          // Add the other Systematic Errors of the Efficiency
          DoubleVecMap variation_Inc , variation_Pl , variation_Mi;
          for (const auto& effT : effTnPType_) {
            //
            // Charged Forward-Backward Ratios
            //
            // Case: Apply Propagation Method
            if (effT.second.second==0 || effT.second.second==2) {
              const std::string vL = Form("N_WToMu_Efficiency_%s_PROP", effT.first.c_str());
              if (c.second.at(binFw).at("Pl").count(vL)==0 || c.second.at(binBw).at("Pl").count(vL)==0 || c.second.at(binFw).at("Mi").count(vL)==0 || c.second.at(binBw).at("Mi").count(vL)==0) {
                std::cout << "[ERROR] " << vL << " variable is missing in bin [" << b.first.etabin().high() << " , " << b.first.etabin().low() << "]" << std::endl; return false;
              }
              const auto& iEVar_FwPl = c.second.at(binFw).at("Pl").at(vL);
              const auto& iEVar_FwMi = c.second.at(binFw).at("Mi").at(vL);
              const auto& iEVar_BwPl = c.second.at(binBw).at("Pl").at(vL);
              const auto& iEVar_BwMi = c.second.at(binBw).at("Mi").at(vL);
              if (systVar[effT.first].size()<1) { systVar.at(effT.first).resize(1); }
              // Plus
              auto& sVar_Pl = systVar.at(effT.first)[0][c.first]["Pl"]["ForwardBackward_Ratio"];
              for (const auto& t : systT) { sVar_Pl[t][b.first] = 0.0; }
              sVar_Pl["Val"][b.first] = oVar_Pl.at("Val").at(b.first);
              // Case: Eta UnCorrelated
              if (effT.second.second==0 || effT.second.second==2) {
                for (const auto& t : systT) {
                  const double error = getForwardBackwardRatioError(iEVar_FwPl.at("Val"), iEVar_BwPl.at("Val"), iEVar_FwPl.at(t), iEVar_BwPl.at(t));
                  oVar_Pl.at(t).at(b.first) = sumErrors( oVar_Pl.at(t).at(b.first) , error );
                  sVar_Pl[t][b.first] = error;
                }
              }
              // Minus
              auto& sVar_Mi = systVar.at(effT.first)[0][c.first]["Mi"]["ForwardBackward_Ratio"];
              for (const auto& t : systT) { sVar_Mi[t][b.first] = 0.0; }
              sVar_Mi["Val"][b.first] = oVar_Mi.at("Val").at(b.first);
              // Case: Eta UnCorrelated
              if (effT.second.second==0 || effT.second.second==2) {
                for (const auto& t : systT) {
                  const double error = getForwardBackwardRatioError(iEVar_FwMi.at("Val"), iEVar_BwMi.at("Val"), iEVar_FwMi.at(t), iEVar_BwMi.at(t));
                  oVar_Mi.at(t).at(b.first) = sumErrors( oVar_Mi.at(t).at(b.first) , error );
                  sVar_Mi[t][b.first] = error;
                }
              }
            }
            // Case: Apply Variation Method
            if (effT.second.second==1 || effT.second.second==3 || effT.second.second>=4) {
              if (systVar[effT.first].size()<effT.second.first) { systVar.at(effT.first).resize(effT.second.first); }
              for (uint i = 0; i < effT.second.first; i++) {
                const std::string vL = std::string("N_WToMu_Efficiency_") + effT.first + ( (effT.second.first>1) ? Form("_%d", i) : "" );
                if (c.second.at(binFw).at("Pl").count(vL)==0 || c.second.at(binBw).at("Pl").count(vL)==0 || c.second.at(binFw).at("Mi").count(vL)==0 || c.second.at(binBw).at("Mi").count(vL)==0) {
                  std::cout << "[ERROR] " << vL << " variable is missing in bin [" << b.first.etabin().high() << " , " << b.first.etabin().low() << "]" << std::endl; return false;
                }
                const auto& iEVar_FwPl = c.second.at(binFw).at("Pl").at(vL);
                const auto& iEVar_FwMi = c.second.at(binFw).at("Mi").at(vL);
                const auto& iEVar_BwPl = c.second.at(binBw).at("Pl").at(vL);
                const auto& iEVar_BwMi = c.second.at(binBw).at("Mi").at(vL);
                // Plus
                auto& sVar_Pl = systVar.at(effT.first)[i][c.first]["Pl"]["ForwardBackward_Ratio"];
                for (const auto& t : systT) { sVar_Pl[t][b.first] = 0.0; }
                sVar_Pl["Val"][b.first] = oVar_Pl.at("Val").at(b.first);
                // Case: Eta Correlated
                if (effT.second.second==1 || effT.second.second==3 || effT.second.second==5) {
                  const double vVar_Pl = getForwardBackwardRatioValue( iEVar_FwPl.at("Val") , iEVar_BwPl.at("Val") );
                  variation_Pl[effT.first].push_back( vVar_Pl );
                  for (const auto& t : systT) { sVar_Pl[t][b.first] = std::abs(vVar_Pl - oVar_Pl.at("Val").at(b.first)); }
                  sVar_Pl["Val"][b.first] = vVar_Pl;
                }
                // Case: Eta UnCorrelated
                else if (effT.second.second==4) {
                  const double vVar_Fw_Pl = getForwardBackwardRatioValue( iEVar_FwPl.at("Val") , iVar_BwPl.at("Val") );
                  const double vVar_Bw_Pl = getForwardBackwardRatioValue( iVar_FwPl.at("Val") , iEVar_BwPl.at("Val") );
                  variation_Pl[(effT.first+"_Fw")].push_back( vVar_Fw_Pl );
                  variation_Pl[(effT.first+"_Bw")].push_back( vVar_Bw_Pl );
                  for (const auto& t : systT) { systVar.at(effT.first)[i][c.first]["Pl"]["ForwardBackward_Ratio_Fw"][t][b.first] = std::abs(vVar_Fw_Pl - oVar_Pl.at("Val").at(b.first)); }
                  for (const auto& t : systT) { systVar.at(effT.first)[i][c.first]["Pl"]["ForwardBackward_Ratio_Bw"][t][b.first] = std::abs(vVar_Bw_Pl - oVar_Pl.at("Val").at(b.first)); }
                }
                // Minus
                auto& sVar_Mi = systVar.at(effT.first)[i][c.first]["Mi"]["ForwardBackward_Ratio"];
                for (const auto& t : systT) { sVar_Mi[t][b.first] = 0.0; }
                sVar_Mi["Val"][b.first] = oVar_Mi.at("Val").at(b.first);
                // Case: Eta Correlated
                if (effT.second.second==1 || effT.second.second==3 || effT.second.second==5) {
                  const double vVar_Mi = getForwardBackwardRatioValue( iEVar_FwMi.at("Val") , iEVar_BwMi.at("Val") );
                  variation_Mi[effT.first].push_back( vVar_Mi );
                  for (const auto& t : systT) { sVar_Mi[t][b.first] = std::abs(vVar_Mi - oVar_Mi.at("Val").at(b.first)); }
                  sVar_Mi["Val"][b.first] = vVar_Mi;
                }
                // Case: Eta UnCorrelated
                else if (effT.second.second==4) {
                  const double vVar_Fw_Mi = getForwardBackwardRatioValue( iEVar_FwMi.at("Val") , iVar_BwMi.at("Val") );
                  const double vVar_Bw_Mi = getForwardBackwardRatioValue( iVar_FwMi.at("Val") , iEVar_BwMi.at("Val") );
                  variation_Mi[(effT.first+"_Fw")].push_back( vVar_Fw_Mi );
                  variation_Mi[(effT.first+"_Bw")].push_back( vVar_Bw_Mi );
                  for (const auto& t : systT) { systVar.at(effT.first)[i][c.first]["Mi"]["ForwardBackward_Ratio_Fw"][t][b.first] = std::abs(vVar_Fw_Mi - oVar_Mi.at("Val").at(b.first)); }
                  for (const auto& t : systT) { systVar.at(effT.first)[i][c.first]["Mi"]["ForwardBackward_Ratio_Bw"][t][b.first] = std::abs(vVar_Bw_Mi - oVar_Mi.at("Val").at(b.first)); }
                }
              }
            }
            //
            // Inclusive Forward-Backward Ratio
            //
            // Case: Apply Propagation Method
            if (effT.second.second==0) {
              const std::string vL = Form("N_WToMu_Efficiency_%s_PROP", effT.first.c_str());
              if (c.second.at(binFw).at("Pl").count(vL)==0 || c.second.at(binBw).at("Pl").count(vL)==0 || c.second.at(binFw).at("Mi").count(vL)==0 || c.second.at(binBw).at("Mi").count(vL)==0) {
                std::cout << "[ERROR] " << vL << " variable is missing in bin [" << b.first.etabin().high() << " , " << b.first.etabin().low() << "]" << std::endl; return false;
              }
              const auto& iEVar_FwPl = c.second.at(binFw).at("Pl").at(vL);
              const auto& iEVar_FwMi = c.second.at(binFw).at("Mi").at(vL);
              const auto& iEVar_BwPl = c.second.at(binBw).at("Pl").at(vL);
              const auto& iEVar_BwMi = c.second.at(binBw).at("Mi").at(vL);
              if (systVar[effT.first].size()<1) { systVar.at(effT.first).resize(1); }
              // Plus
              auto& sVar_Inc = systVar.at(effT.first)[0][c.first][""]["ForwardBackward_Ratio"];
              sVar_Inc["Val"][b.first] = oVar_Inc.at("Val").at(b.first);
              for (const auto& t : systT) {
                // Case: Eta UnCorrelated
                if (effT.second.second==0) {
                  const double error = getForwardBackwardRatioError( iEVar_FwPl.at("Val") , iEVar_BwPl.at("Val") , iEVar_FwMi.at("Val") , iEVar_BwMi.at("Val"), 
                                                                     iEVar_FwPl.at(t)     , iEVar_BwPl.at(t)     , iEVar_FwMi.at(t)     , iEVar_BwMi.at(t) );
                  oVar_Inc.at(t).at(b.first) = sumErrors( oVar_Inc.at(t).at(b.first) , error );
                  sVar_Inc[t][b.first] = error;
                }
              }
            }
            // Case: Apply Variation Method
            if (effT.second.second>0) {
              if (systVar[effT.first].size()<effT.second.first) { systVar.at(effT.first).resize(effT.second.first); }
              for (uint i = 0; i < effT.second.first; i++) {
                const std::string vL = std::string("N_WToMu_Efficiency_") + effT.first + ( (effT.second.first>1) ? Form("_%d", i) : "" );
                if (c.second.at(binFw).at("Pl").count(vL)==0 || c.second.at(binBw).at("Pl").count(vL)==0 || c.second.at(binFw).at("Mi").count(vL)==0 || c.second.at(binBw).at("Mi").count(vL)==0) {
                  std::cout << "[ERROR] " << vL << " variable is missing in bin [" << b.first.etabin().high() << " , " << b.first.etabin().low() << "]" << std::endl; return false;
                }
                const auto& iEVar_FwPl = c.second.at(binFw).at("Pl").at(vL);
                const auto& iEVar_FwMi = c.second.at(binFw).at("Mi").at(vL);
                const auto& iEVar_BwPl = c.second.at(binBw).at("Pl").at(vL);
                const auto& iEVar_BwMi = c.second.at(binBw).at("Mi").at(vL);
                // Inclusive
                auto& sVar_Inc = systVar.at(effT.first)[i][c.first][""]["ForwardBackward_Ratio"];
                for (const auto& t : systT) { sVar_Inc[t][b.first] = 0.0; }
                sVar_Inc["Val"][b.first] = oVar_Inc.at("Val").at(b.first);
                // Case: Fully Correlated
                if (effT.second.second==3 || effT.second.second==5) {
                  const double vVar_Inc = getForwardBackwardRatioValue( ( iEVar_FwPl.at("Val") + iEVar_FwMi.at("Val") ) , ( iEVar_BwPl.at("Val") + iEVar_BwMi.at("Val") ) );
                  variation_Inc[effT.first].push_back( vVar_Inc );
                  for (const auto& t : systT) { sVar_Inc[t][b.first] = std::abs(vVar_Inc - oVar_Inc.at("Val").at(b.first)); }
                  sVar_Inc["Val"][b.first] = vVar_Inc;
                }
                // Case: Charge Correlated and Eta UnCorrelated
                if (effT.second.second==2) {
                  const double vVar_Fw = getForwardBackwardRatioValue( ( iEVar_FwPl.at("Val") + iEVar_FwMi.at("Val") ) , (  iVar_BwPl.at("Val") +  iVar_BwMi.at("Val") ) );
                  const double vVar_Bw = getForwardBackwardRatioValue( (  iVar_FwPl.at("Val") +  iVar_FwMi.at("Val") ) , ( iEVar_BwPl.at("Val") + iEVar_BwMi.at("Val") ) );
                  variation_Inc[(effT.first+"_Fw")].push_back( vVar_Fw );
                  variation_Inc[(effT.first+"_Bw")].push_back( vVar_Bw );
                  for (const auto& t : systT) { systVar.at(effT.first)[i][c.first][""]["ForwardBackward_Ratio_Fw"][t][b.first] = std::abs(vVar_Fw - oVar_Inc.at("Val").at(b.first)); }
                  for (const auto& t : systT) { systVar.at(effT.first)[i][c.first][""]["ForwardBackward_Ratio_Bw"][t][b.first] = std::abs(vVar_Bw - oVar_Inc.at("Val").at(b.first)); }
                }
                // Case: Charge UnCorrelated and Eta Correlated
                if (effT.second.second==1) {
                  const double vVar_Pl = getForwardBackwardRatioValue( ( iEVar_FwPl.at("Val") +  iVar_FwMi.at("Val") ) , ( iEVar_BwPl.at("Val") +  iVar_BwMi.at("Val") ) );
                  const double vVar_Mi = getForwardBackwardRatioValue( (  iVar_FwPl.at("Val") + iEVar_FwMi.at("Val") ) , (  iVar_BwPl.at("Val") + iEVar_BwMi.at("Val") ) );
                  variation_Inc[(effT.first+"_Pl")].push_back( vVar_Pl );
                  variation_Inc[(effT.first+"_Mi")].push_back( vVar_Mi );
                  for (const auto& t : systT) { systVar.at(effT.first)[i][c.first][""]["ForwardBackward_Ratio_Pl"][t][b.first] = std::abs(vVar_Pl - oVar_Inc.at("Val").at(b.first)); }
                  for (const auto& t : systT) { systVar.at(effT.first)[i][c.first][""]["ForwardBackward_Ratio_Mi"][t][b.first] = std::abs(vVar_Mi - oVar_Inc.at("Val").at(b.first)); }
                }
                // Case: Fully UnCorrelated
                if (effT.second.second==4) {
                  const double vVar_Fw_Pl = getForwardBackwardRatioValue( ( iEVar_FwPl.at("Val") +  iVar_FwMi.at("Val") ) , (  iVar_BwPl.at("Val") +  iVar_BwMi.at("Val") ) );
                  const double vVar_Bw_Pl = getForwardBackwardRatioValue( (  iVar_FwPl.at("Val") +  iVar_FwMi.at("Val") ) , ( iEVar_BwPl.at("Val") +  iVar_BwMi.at("Val") ) );
                  const double vVar_Fw_Mi = getForwardBackwardRatioValue( (  iVar_FwPl.at("Val") + iEVar_FwMi.at("Val") ) , (  iVar_BwPl.at("Val") +  iVar_BwMi.at("Val") ) );
                  const double vVar_Bw_Mi = getForwardBackwardRatioValue( (  iVar_FwPl.at("Val") +  iVar_FwMi.at("Val") ) , (  iVar_BwPl.at("Val") + iEVar_BwMi.at("Val") ) );
                  variation_Inc[(effT.first+"_Fw_Pl")].push_back( vVar_Fw_Pl );
                  variation_Inc[(effT.first+"_Bw_Pl")].push_back( vVar_Bw_Pl );
                  variation_Inc[(effT.first+"_Fw_Mi")].push_back( vVar_Fw_Mi );
                  variation_Inc[(effT.first+"_Bw_Mi")].push_back( vVar_Bw_Mi );
                  for (const auto& t : systT) { systVar.at(effT.first)[i][c.first][""]["ForwardBackward_Ratio_Fw_Pl"][t][b.first] = std::abs(vVar_Fw_Pl - oVar_Inc.at("Val").at(b.first)); }
                  for (const auto& t : systT) { systVar.at(effT.first)[i][c.first][""]["ForwardBackward_Ratio_Bw_Pl"][t][b.first] = std::abs(vVar_Bw_Pl - oVar_Inc.at("Val").at(b.first)); }
                  for (const auto& t : systT) { systVar.at(effT.first)[i][c.first][""]["ForwardBackward_Ratio_Fw_Mi"][t][b.first] = std::abs(vVar_Fw_Mi - oVar_Inc.at("Val").at(b.first)); }
                  for (const auto& t : systT) { systVar.at(effT.first)[i][c.first][""]["ForwardBackward_Ratio_Bw_Mi"][t][b.first] = std::abs(vVar_Bw_Mi - oVar_Inc.at("Val").at(b.first)); }
                }
              }
            }
          }
          for (const auto& vEr : variation_Inc) {
            computeTnPUncertainty(oVar_Inc.at("Err_Syst_High").at(b.first), oVar_Inc.at("Err_Syst_Low").at(b.first), variation_Inc.at(vEr.first), oVar_Inc.at("Val").at(b.first), vEr.first);
          }
          for (const auto& vEr : variation_Pl) {
            computeTnPUncertainty(oVar_Pl.at("Err_Syst_High").at(b.first), oVar_Pl.at("Err_Syst_Low").at(b.first), variation_Pl.at(vEr.first), oVar_Pl.at("Val").at(b.first), vEr.first);
            computeTnPUncertainty(oVar_Mi.at("Err_Syst_High").at(b.first), oVar_Mi.at("Err_Syst_Low").at(b.first), variation_Mi.at(vEr.first), oVar_Mi.at("Val").at(b.first), vEr.first);
          }
        }
      }
      else if (doSyst && nomVar.size()>0) {
	//
	const auto& nVar_FwPl = nomVar.at(c.first).at(binFw).at("Pl").at("N_WToMu");
	const auto& nVar_FwMi = nomVar.at(c.first).at(binFw).at("Mi").at("N_WToMu");
	const auto& nVar_BwPl = nomVar.at(c.first).at(binBw).at("Pl").at("N_WToMu");
	const auto& nVar_BwMi = nomVar.at(c.first).at(binBw).at("Mi").at("N_WToMu");
	//
	// Add the Systematic Error
	//
	// Charged Forward-Backward Ratios
	//
	// Case: Apply Variation Method
        //
	const double nomVal_Pl = getForwardBackwardRatioValue( nVar_FwPl.at("Val") , nVar_BwPl.at("Val") );
	const double nomVal_Mi = getForwardBackwardRatioValue( nVar_FwMi.at("Val") , nVar_BwMi.at("Val") );
	oVar_Pl.at("Val").at(b.first) = nomVal_Pl;
	oVar_Mi.at("Val").at(b.first) = nomVal_Mi;
        //
        // Case: Eta Correlated
        if (corrSyst==1 || corrSyst==3 || corrSyst==5) {
          const double varVal_Pl = getForwardBackwardRatioValue(iVar_FwPl.at("Val") , iVar_BwPl.at("Val"));
          const double varVal_Mi = getForwardBackwardRatioValue(iVar_FwMi.at("Val") , iVar_BwMi.at("Val"));
          for (const auto& t : systT) {
	    oVar_Pl.at(t).at(b.first) = std::abs(varVal_Pl - nomVal_Pl);
	    oVar_Mi.at(t).at(b.first) = std::abs(varVal_Mi - nomVal_Mi);
          }
          oVar_Pl.at("Val").at(b.first) = varVal_Pl;
          oVar_Mi.at("Val").at(b.first) = varVal_Mi;
        }
        // Case: Eta UnCorrelated
        else if (corrSyst==0 || corrSyst==2 || corrSyst==4) {
          const double vVar_FwPl = getForwardBackwardRatioValue(iVar_FwPl.at("Val"), nVar_BwPl.at("Val"));
          const double vVar_BwPl = getForwardBackwardRatioValue(nVar_FwPl.at("Val"), iVar_BwPl.at("Val"));
          const double vVar_FwMi = getForwardBackwardRatioValue(iVar_FwMi.at("Val"), nVar_BwMi.at("Val"));
          const double vVar_BwMi = getForwardBackwardRatioValue(nVar_FwMi.at("Val"), iVar_BwMi.at("Val"));
          for (const auto& t : systT) {
	    var.at(c.first).at("Pl")["ForwardBackward_Ratio_Fw"][t][b.first] = std::abs(vVar_FwPl - nomVal_Pl);
	    var.at(c.first).at("Pl")["ForwardBackward_Ratio_Bw"][t][b.first] = std::abs(vVar_BwPl - nomVal_Pl);
	    var.at(c.first).at("Mi")["ForwardBackward_Ratio_Fw"][t][b.first] = std::abs(vVar_FwMi - nomVal_Mi);
	    var.at(c.first).at("Mi")["ForwardBackward_Ratio_Bw"][t][b.first] = std::abs(vVar_BwMi - nomVal_Mi);
	  }
	}
	//
	// Inclusive Forward-Backward Ratio
	//
	// Case: Apply Variation Method
        //
	const double nomVal_Inc = getForwardBackwardRatioValue( ( nVar_FwPl.at("Val") + nVar_FwMi.at("Val") ) , ( nVar_BwPl.at("Val") + nVar_BwMi.at("Val") ) );
	oVar_Inc.at("Val").at(b.first) = nomVal_Inc;
        //
        // Case: Fully Correlated
        if (corrSyst==3 || corrSyst==5) {
          const double varVal_Inc = getForwardBackwardRatioValue((iVar_FwPl.at("Val")+iVar_FwMi.at("Val")) , (iVar_BwPl.at("Val")+iVar_BwMi.at("Val")));
          for (const auto& t : systT) { oVar_Inc.at(t).at(b.first) = std::abs(varVal_Inc - nomVal_Inc); }
          oVar_Inc.at("Val").at(b.first) = varVal_Inc;
        }
        // Case: Charge Correlated and Eta UnCorrelated
        if (corrSyst==2) {
          const double vVar_Fw = getForwardBackwardRatioValue((iVar_FwPl.at("Val")+iVar_FwMi.at("Val")), (nVar_BwPl.at("Val")+nVar_BwMi.at("Val")));
          const double vVar_Bw = getForwardBackwardRatioValue((nVar_FwPl.at("Val")+nVar_FwMi.at("Val")), (iVar_BwPl.at("Val")+iVar_BwMi.at("Val")));
          for (const auto& t : systT) {
            var.at(c.first).at("")["ForwardBackward_Ratio_Fw"][t][b.first] = std::abs(vVar_Fw - nomVal_Inc);
            var.at(c.first).at("")["ForwardBackward_Ratio_Bw"][t][b.first] = std::abs(vVar_Bw - nomVal_Inc);
	  }
        }
        // Case: Charge UnCorrelated and Eta Correlated
        if (corrSyst==1) {
          const double vVar_Pl = getForwardBackwardRatioValue((iVar_FwPl.at("Val")+nVar_FwMi.at("Val")), (iVar_BwPl.at("Val")+nVar_BwMi.at("Val")));
          const double vVar_Mi = getForwardBackwardRatioValue((nVar_FwPl.at("Val")+iVar_FwMi.at("Val")), (nVar_BwPl.at("Val")+iVar_BwMi.at("Val")));
          for (const auto& t : systT) {
	    var.at(c.first).at("")["ForwardBackward_Ratio_Pl"][t][b.first] = std::abs(vVar_Pl - nomVal_Inc);
	    var.at(c.first).at("")["ForwardBackward_Ratio_Mi"][t][b.first] = std::abs(vVar_Mi - nomVal_Inc);
	  }
        }
        // Case: Fully UnCorrelated
        if (corrSyst==0 || corrSyst==4) {
          const double vVar_Fw_Pl = getForwardBackwardRatioValue((iVar_FwPl.at("Val")+nVar_FwMi.at("Val")), (nVar_BwPl.at("Val")+nVar_BwMi.at("Val")));
          const double vVar_Bw_Pl = getForwardBackwardRatioValue((nVar_FwPl.at("Val")+nVar_FwMi.at("Val")), (iVar_BwPl.at("Val")+nVar_BwMi.at("Val")));
          const double vVar_Fw_Mi = getForwardBackwardRatioValue((nVar_FwPl.at("Val")+iVar_FwMi.at("Val")), (nVar_BwPl.at("Val")+nVar_BwMi.at("Val")));
          const double vVar_Bw_Mi = getForwardBackwardRatioValue((nVar_FwPl.at("Val")+nVar_FwMi.at("Val")), (nVar_BwPl.at("Val")+iVar_BwMi.at("Val")));
          for (const auto& t : systT) {
	    var.at(c.first).at("")["ForwardBackward_Ratio_Fw_Pl"][t][b.first] = std::abs(vVar_Fw_Pl - nomVal_Inc);
	    var.at(c.first).at("")["ForwardBackward_Ratio_Bw_Pl"][t][b.first] = std::abs(vVar_Bw_Pl - nomVal_Inc);
	    var.at(c.first).at("")["ForwardBackward_Ratio_Fw_Mi"][t][b.first] = std::abs(vVar_Fw_Mi - nomVal_Inc);
	    var.at(c.first).at("")["ForwardBackward_Ratio_Bw_Mi"][t][b.first] = std::abs(vVar_Bw_Mi - nomVal_Inc);
	  }
	}
      }
    }
  }
  return true;
};


double getCrossSectionValue( const double& N , const double& Luminosity , const double& BinWidth )
{
  return ( N / ( Luminosity * BinWidth ) );
};

  
double getCrossSectionError( const double& N , const double& Luminosity , const double& BinWidth , const double& Err_N , const double& Err_Luminosity )
{
  const double common  = getCrossSectionValue( N , Luminosity , BinWidth );
  const double error2N = std::pow( ( Err_N / N ) , 2.0 );
  const double error2L = std::pow( ( Err_Luminosity / Luminosity ) , 2.0 );
  return ( std::abs(common) * std::sqrt( error2N + error2L ) );
};


bool computeCrossSection( BinPentaMap& var , BinSextaMapVec& systVar , const VarBinMap& inputVar , const bool& doSyst , const VarBinMap& nomVar , const uint& corrSyst )
{
  //
  for (const auto& c : inputVar) {
    for (const auto& b : c.second) {
      // Check that everything is fine
      if (b.second.count("Pl")==0) { std::cout << "[ERROR] Plus charge in " << c.first << " is missing in bin [" << b.first.etabin().high() << " , " << b.first.etabin().low() << "]" << std::endl; return false; }
      if (b.second.count("Mi")==0) { std::cout << "[ERROR] Minus charge in " << c.first << " is missing in bin [" << b.first.etabin().high() << " , " << b.first.etabin().low() << "]" << std::endl; return false; }
      if (b.second.at("Mi").count("N_WToMu")==0 || b.second.at("Pl").count("N_WToMu")==0) {
        std::cout << "[ERROR] N_WToMu variable is missing in bin [" << b.first.etabin().high() << " , " << b.first.etabin().low() << "]" << std::endl; return false;
      }
      //
      auto&       oVar_Pl = var[c.first]["Pl"]["Cross_Section"];
      auto&       oVar_Mi = var[c.first]["Mi"]["Cross_Section"];
      const auto& iVar_Pl = b.second.at("Pl").at("N_WToMu");
      const auto& iVar_Mi = b.second.at("Mi").at("N_WToMu");
      //
      const double Luminosity = b.second.at("Pl").at("Luminosity").at("Val");
      std::map< std::string , double > Err_Luminosity;
      Err_Luminosity["Err_Syst_High"] = 0.0; // Not assigned
      Err_Luminosity["Err_Syst_Low" ] = 0.0; // Not assigned
      Err_Luminosity["Err_Stat_High"] = 0.0; // Not assigned
      Err_Luminosity["Err_Stat_Low" ] = 0.0; // Not assigned
      //
      const double BinWidth = ( b.first.etabin().high() - b.first.etabin().low() );
      //
      // Compute the cross-section value
      oVar_Pl["Val"][b.first] = getCrossSectionValue(iVar_Pl.at("Val"), Luminosity, BinWidth);
      oVar_Mi["Val"][b.first] = getCrossSectionValue(iVar_Mi.at("Val"), Luminosity, BinWidth);
      //
      // Compute the cross-section statistical error
      const std::vector< std::string > statT = { "Err_Stat_High" , "Err_Stat_Low" };
      for (const auto& t : statT) { oVar_Pl[t][b.first] = getCrossSectionError(iVar_Pl.at("Val"), Luminosity, BinWidth, iVar_Pl.at(t), Err_Luminosity.at(t)); }
      for (const auto& t : statT) { oVar_Mi[t][b.first] = getCrossSectionError(iVar_Mi.at("Val"), Luminosity, BinWidth, iVar_Mi.at(t), Err_Luminosity.at(t)); }
      //
      // Compute the cross-section systematic error (in this case from efficiency TnP)
      const std::vector< std::string > systT = { "Err_Syst_High" , "Err_Syst_Low" };
      // Initialize the Systematic Uncertainties
      for (const auto& t : systT) { oVar_Pl[t][b.first] = 0.0; }
      for (const auto& t : systT) { oVar_Mi[t][b.first] = 0.0; }
      //
      if (doSyst && nomVar.size()==0) {
        //
        std::string effType = ""; if (b.second.at("Pl").count("N_WToMu_Efficiency_TnP")>0) { effType = "TnP"; } else if (b.second.at("Pl").count("N_WToMu_Efficiency_MC")>0) { effType = "MC"; }
        //
        // If we are not using TnP efficieny, simply progate the uncertainties
        if (effType!="TnP") {
          for (const auto& t : systT) { oVar_Pl[t][b.first] = getCrossSectionError(iVar_Pl.at("Val"), Luminosity, BinWidth, iVar_Pl.at(t), Err_Luminosity.at(t)); }
          for (const auto& t : systT) { oVar_Mi[t][b.first] = getCrossSectionError(iVar_Mi.at("Val"), Luminosity, BinWidth, iVar_Mi.at(t), Err_Luminosity.at(t)); }
        }
        else {
          // Add the Statistical Error of the Efficiency
          const auto& iMCVar_Pl = b.second.at("Pl").at(Form("N_WToMu_Efficiency_%s", effType.c_str()));
          const auto& iMCVar_Mi = b.second.at("Mi").at(Form("N_WToMu_Efficiency_%s", effType.c_str()));
          //
          if (systVar["MC_Statistics"].size()<1) { systVar.at("MC_Statistics").resize(1); }
          // Plus
          auto& sMCVar_Pl = systVar.at("MC_Statistics")[0][c.first]["Pl"]["Cross_Section"];
          sMCVar_Pl["Val"][b.first] = oVar_Pl.at("Val").at(b.first);
          for (uint i=0; i<systT.size(); i++) {
            const double error_Pl = getCrossSectionError(iMCVar_Pl.at("Val"), Luminosity, BinWidth, iMCVar_Pl.at(statT[i]), Err_Luminosity.at(systT[i]));
            oVar_Pl.at(systT[i]).at(b.first) = sumErrors( oVar_Pl.at(systT[i]).at(b.first) , error_Pl );
            sMCVar_Pl[systT[i]][b.first] = error_Pl;
          }
          // Minus
          auto& sMCVar_Mi = systVar.at("MC_Statistics")[0][c.first]["Mi"]["Cross_Section"];
          sMCVar_Mi["Val"][b.first] = oVar_Mi.at("Val").at(b.first);
          for (uint i=0; i<systT.size(); i++) {
            const double error_Mi = getCrossSectionError(iMCVar_Mi.at("Val"), Luminosity, BinWidth, iMCVar_Mi.at(statT[i]), Err_Luminosity.at(systT[i]));
            oVar_Mi.at(systT[i]).at(b.first) = sumErrors( oVar_Mi.at(systT[i]).at(b.first) , error_Mi );
            sMCVar_Mi[systT[i]][b.first] = error_Mi;
          }
          // Add the other Systematic Errors of the Efficiency
          DoubleVecMap variation_Pl , variation_Mi;
          for (const auto& effT : effTnPType_) {
            //
            // Case: Apply Propagation Method
            if (effT.second.second<4) {  //false) { //
              const std::string vL = Form("N_WToMu_Efficiency_%s_PROP", effT.first.c_str());
              if (b.second.at("Mi").count(vL)==0 || b.second.at("Pl").count(vL)==0) {
                std::cout << "[ERROR] " << vL << " variable is missing in bin [" << b.first.etabin().high() << " , " << b.first.etabin().low() << "]" << std::endl; return false;
              }
              if (systVar[effT.first].size()<1) { systVar.at(effT.first).resize(1); }
              // Plus
              auto& sVar_Pl = systVar.at(effT.first)[0][c.first]["Pl"]["Cross_Section"];
              sVar_Pl["Val"][b.first] = oVar_Pl.at("Val").at(b.first);
              for (const auto& t : systT) {
                const double error_Pl = getCrossSectionError(b.second.at("Pl").at(vL).at("Val"), Luminosity, BinWidth, b.second.at("Pl").at(vL).at(t), 0.0);
                oVar_Pl.at(t).at(b.first) = sumErrors( oVar_Pl.at(t).at(b.first) , error_Pl );
                sVar_Pl[t][b.first] = error_Pl;
              }
              // Minus
              auto& sVar_Mi = systVar.at(effT.first)[0][c.first]["Mi"]["Cross_Section"];
              sVar_Mi["Val"][b.first] = oVar_Mi.at("Val").at(b.first);
              for (const auto& t : systT) {
                const double error_Mi = getCrossSectionError(b.second.at("Mi").at(vL).at("Val"), Luminosity, BinWidth, b.second.at("Mi").at(vL).at(t), 0.0);
                oVar_Mi.at(t).at(b.first) = sumErrors( oVar_Mi.at(t).at(b.first) , error_Mi );
                sVar_Mi[t][b.first] = error_Mi;
              }
            }
            // Case: Apply Variation Method
            if (effT.second.second>=4) { //true) { //
              if (systVar[effT.first].size()<effT.second.first) { systVar.at(effT.first).resize(effT.second.first); }
              for (uint i = 0; i < effT.second.first; i++) {
                const std::string vL = std::string("N_WToMu_Efficiency_") + effT.first + ( (effT.second.first>1) ? Form("_%d", i) : "" );
                if (b.second.at("Mi").count(vL)==0 || b.second.at("Pl").count(vL)==0) {
                  std::cout << "[ERROR] " << vL << " variable is missing in bin [" << b.first.etabin().high() << " , " << b.first.etabin().low() << "]" << std::endl; return false;
                }
                // Plus
                const double vVar_Pl = getCrossSectionValue(b.second.at("Pl").at(vL).at("Val"), Luminosity, BinWidth);
                variation_Pl[effT.first].push_back( vVar_Pl );
                auto& sVar_Pl = systVar.at(effT.first)[i][c.first]["Pl"]["Cross_Section"];
                for (const auto& t : systT) { sVar_Pl[t][b.first] = std::abs(vVar_Pl - oVar_Pl.at("Val").at(b.first)); }
                sVar_Pl["Val"][b.first] = vVar_Pl;
                // Minus
                const double vVar_Mi = getCrossSectionValue(b.second.at("Mi").at(vL).at("Val"), Luminosity, BinWidth);
                variation_Mi[effT.first].push_back( vVar_Mi );
                auto& sVar_Mi = systVar.at(effT.first)[i][c.first]["Mi"]["Cross_Section"];
                for (const auto& t : systT) { sVar_Mi[t][b.first] = std::abs(vVar_Mi - oVar_Mi.at("Val").at(b.first)); }
                sVar_Mi["Val"][b.first] = vVar_Mi;
              }
            }
          }
          for (const auto& vEr : variation_Pl) {
            computeTnPUncertainty(oVar_Pl.at("Err_Syst_High").at(b.first), oVar_Pl.at("Err_Syst_Low").at(b.first), variation_Pl.at(vEr.first), oVar_Pl.at("Val").at(b.first), vEr.first);
            computeTnPUncertainty(oVar_Mi.at("Err_Syst_High").at(b.first), oVar_Mi.at("Err_Syst_Low").at(b.first), variation_Mi.at(vEr.first), oVar_Mi.at("Val").at(b.first), vEr.first);
          }
        }
      }
      else if (doSyst && nomVar.size()>0) {
	//
	const auto& nVar_Pl = nomVar.at(c.first).at(b.first).at("Pl").at("N_WToMu");
	const auto& nVar_Mi = nomVar.at(c.first).at(b.first).at("Mi").at("N_WToMu");
	//
	// Add the Systematic Error
	//
	// Case: Apply Variation Method
        //
	const double nomVal_Pl = getCrossSectionValue(nVar_Pl.at("Val"), Luminosity, BinWidth);
	const double nomVal_Mi = getCrossSectionValue(nVar_Mi.at("Val"), Luminosity, BinWidth);
	const double varVal_Pl = getCrossSectionValue(iVar_Pl.at("Val"), Luminosity, BinWidth);
	const double varVal_Mi = getCrossSectionValue(iVar_Mi.at("Val"), Luminosity, BinWidth);
        //
	for (const auto& t : systT) {
	  oVar_Pl.at(t).at(b.first) = std::abs(varVal_Pl - nomVal_Pl);
	  oVar_Mi.at(t).at(b.first) = std::abs(varVal_Mi - nomVal_Mi);
	}
	oVar_Pl.at("Val").at(b.first) = varVal_Pl;
	oVar_Mi.at("Val").at(b.first) = varVal_Mi;
      }
    }
  }
  return true;
};


void procSyst(std::map<anabin<0>,double>& sys, const std::string type="Mean")
{
  if (type=="Mean") {
    double mean = 0.0;
    for (const auto& s : sys) { mean += std::abs(s.second); }
    mean /= sys.size();
    for (auto& s : sys) { s.second = mean; }
  }
  else if (type=="Max") {
    double max = -99999999.0;
    for (const auto& s : sys) { max = std::max(max, std::abs(s.second)); }
    for (auto& s : sys) { s.second = max; }
  }
  else if (type=="RMS") {
    double RMS = 0.0;
    for (const auto& s : sys) { RMS += (s.second*s.second); }
    RMS /= sys.size();
    RMS = std::sqrt(RMS);
    for (auto& s : sys) { s.second = RMS; }
  }
  else if (type=="Smooth") {
    const uint nBins = sys.size();
    double xx[nBins]; uint i = 0;
    i=0; for (const auto& s : sys) { xx[i] = s.second; i++; }
    TH1::SmoothArray(nBins, xx, 1);
    i=0; for (auto& s : sys) { s.second = xx[i]; i++; }
  }
};


void computeSystematic(BinSextaMap& varMap, BinSeptaMapVec& systVarMap)
{
  //
  auto& nomVar = varMap.at("Nominal");
  for (const auto& cat : systVarMap) {
    auto& var = varMap[cat.first];
    auto& systVarVec = systVarMap.at(cat.first);
    const auto origVarVec = systVarVec;
    var.clear(); systVarVec.clear();
    // Combine the sub-uncertainties of each systematic variation
    for (const auto& lbl : origVarVec) {
      systVarVec[lbl.first].clear();
      double maxErr = -9999999., minErr = 9999999.;
      for (uint iVr=0; iVr<lbl.second.size(); iVr++) {
	BinPentaMap systVar;
	for (const auto& c : lbl.second[iVr]) {
	  for (const auto& chg : c.second) {
	    // Combine the varied uncertainties
            auto& valMap = systVar[c.first][chg.first];
	    for (const auto& v : chg.second) {
	      std::string vLbl = v.first;
	      if (std::count(v.first.begin(), v.first.end(), '_')>1) {
		std::string tmp = vLbl.substr(vLbl.find("_")+1); vLbl = vLbl.substr(0, vLbl.find("_"))+"_"+tmp.substr(0, tmp.find("_"));
	      }
              valMap[vLbl]["Val"] = chg.second.at(vLbl).at("Val");
	      for (const auto& t : v.second) {
		if (t.first.find("Err_Syst_")!=std::string::npos) {
		  for (const auto& b : t.second) {
		    auto& val = valMap[vLbl][t.first][b.first];
		    val = sumErrors(val, b.second);
                    minErr = std::min( minErr , val ); maxErr = std::max( maxErr , val );
		  }
		}
              }
	      valMap[vLbl]["Nom"] = nomVar.at(c.first).at(chg.first).at(vLbl).at("Val");
              for (const auto& b : valMap.at(vLbl).at("Err_Syst_High")) {
                if (valMap.at(vLbl).at("Val").at(b.first)==valMap.at(vLbl).at("Nom").at(b.first)) {
                  valMap.at(vLbl).at("Val").at(b.first) += std::max(valMap.at(vLbl).at("Err_Syst_Low").at(b.first), valMap.at(vLbl).at("Err_Syst_High").at(b.first));
                }
              }
	    }
	  }
	}
        if (minErr==0.0 && maxErr==0.0) { std::cout << "[WARNING] Variation " << lbl.first << " index " << iVr << " is empty. Ignoring it!" << std::endl; break; }
	systVarVec.at(lbl.first).push_back(systVar);
      }
      if (minErr==0.0 && maxErr==0.0) { systVarVec.erase(lbl.first); continue; }
      // Compute the uncertainty of each systematic variation
      BinPentaMap systVar;
      for (const auto& c : lbl.second[0]) {
	for (const auto& chg : c.second) {
	  // Calculate the uncertainties
	  auto& valMap = systVar[c.first][chg.first];
	  for (const auto& v : chg.second) {
	    std::string vLbl = v.first;
	    if (std::count(v.first.begin(), v.first.end(), '_')>1) {
	      std::string tmp = vLbl.substr(vLbl.find("_")+1); vLbl = vLbl.substr(0, vLbl.find("_"))+"_"+tmp.substr(0, tmp.find("_"));
	    }
            valMap[vLbl]["Val"] = chg.second.at(vLbl).at("Val");
	    for (const auto& t : v.second) {
	      if (t.first.find("Err_Syst_")!=std::string::npos) {
                for (const auto& b : t.second) {
                  auto& val = valMap[vLbl][t.first][b.first];
                  uint nVariation = lbl.second.size();
                  if (nVariation>1 && ( (lbl.second[1].at(c.first).count(chg.first)==0) || (lbl.second[1].at(c.first).at(chg.first).count(v.first)==0) )) { nVariation = 1; }
                  double uncVal = 0.0;
                  //
                  if (nVariation == 1) {
                    const double diff = std::abs(lbl.second[0].at(c.first).at(chg.first).at(v.first).at(t.first).at(b.first));
                    uncVal = diff;
                  }
                  else if (lbl.first.find("MC_Syst_PDF")!=std::string::npos) {
                    // Variations (Use the offical EPPS16 approach)
                    double Err_High = 0.0 , Err_Low = 0.0;
                    for(uint i = 0; i < (nVariation/2); i++) {
                      const double diff_0 = lbl.second[(2*i)+0].at(c.first).at(chg.first).at(v.first).at(t.first).at(b.first);
                      const double diff_1 = lbl.second[(2*i)+1].at(c.first).at(chg.first).at(v.first).at(t.first).at(b.first);
                      Err_High += std::pow( std::max( std::max( diff_1 , diff_0 ) , 0.0 ) , 2.0 );
                      Err_Low  += std::pow( std::min( std::min( diff_1 , diff_0 ) , 0.0 ) , 2.0 );
                    }
                    // Convert from 90% CL to 68% CL
                    const double convFactor = TMath::ErfcInverse((1.-0.68))/TMath::ErfcInverse((1.-0.90));
                    uncVal = ( convFactor * std::max( std::sqrt(Err_High) , std::sqrt(Err_Low) ) );
                  }
                  else if (lbl.first.find("MC_Syst_Scale")!=std::string::npos) {
                    // Variations (Use the envelope approach)
                    double Err_High = 0.0 , Err_Low = 0.0;
                    for(uint i = 0; i < nVariation; i++) {
                      const double diff = lbl.second[i].at(c.first).at(chg.first).at(v.first).at(t.first).at(b.first);
                      Err_High = std::max( std::max( diff , Err_High ) , 0.0 );
                      Err_Low  = std::min( std::min( diff , Err_Low  ) , 0.0 );
                    }
                    uncVal = std::max( std::abs(Err_High) , std::abs(Err_Low) );
                  }
                  else if (lbl.first.find("MC_Syst_Alpha")!=std::string::npos) {
                    // Variations (Use the PDF4LHC15 approach)
                    const double diff_0 = lbl.second[0].at(c.first).at(chg.first).at(v.first).at(t.first).at(b.first);
                    const double diff_1 = lbl.second[1].at(c.first).at(chg.first).at(v.first).at(t.first).at(b.first);
                    const double err = ( ( std::abs(diff_1) + std::abs(diff_0) ) / 2.0 );
                    // Convert from deltaAlpha_s = 0.001 to 0.0015 (68% CL)
                    const double convFactor = (0.0015/0.0010);
                    uncVal = ( convFactor * std::abs( err ) );
                  }
                  else if (nVariation == 2) {
                    const double diff_0 = std::abs(lbl.second[0].at(c.first).at(chg.first).at(v.first).at(t.first).at(b.first));
                    const double diff_1 = std::abs(lbl.second[1].at(c.first).at(chg.first).at(v.first).at(t.first).at(b.first));
                    uncVal = std::max( diff_0 , diff_1 );
                  }
                  else if (nVariation > 2) {
                    double sum = 0.0;
                    for (uint i = 0; i < nVariation; i++) {
                      const double diff = std::abs(lbl.second[i].at(c.first).at(chg.first).at(v.first).at(t.first).at(b.first));
                      sum += ( diff * diff );
                    }
                    uncVal = std::sqrt( sum / nVariation );
                  }
                  if (vLbl==v.first) { val = uncVal; }
                  else { val = sumErrors(val, uncVal); }
                }
              }
	    }
            valMap[vLbl]["Nom"] = nomVar.at(c.first).at(chg.first).at(vLbl).at("Val");
            for (const auto& b : valMap.at(vLbl).at("Err_Syst_High")) {
              if (valMap.at(vLbl).at("Val").at(b.first)==valMap.at(vLbl).at("Nom").at(b.first)) {
                valMap.at(vLbl).at("Val").at(b.first) += std::max(valMap.at(vLbl).at("Err_Syst_High").at(b.first), valMap.at(vLbl).at("Err_Syst_High").at(b.first));
              }
            }
          }
          // Apply the mean to QCD variations
          if (lbl.first.find("_Avg")!=std::string::npos) {
            for (auto& sV : valMap) { for (auto& sT : sV.second) { if (sT.first.find("Err_Syst_")!=std::string::npos) { procSyst(valMap.at(sV.first).at(sT.first), "Mean"); } } }
          }
          // Apply to Recoil variations
          else if ((lbl.first.find("_Smearing")!=std::string::npos) || (lbl.first.find("_SystPtFunc")!=std::string::npos)
                   || (lbl.first.find("_SystJetEn")!=std::string::npos) || (lbl.first.find("_SystBWGauss")!=std::string::npos)) {
            for (auto& sV : valMap) { for (auto& sT : sV.second) { if (sT.first.find("Err_Syst_")!=std::string::npos) { procSyst(valMap.at(sV.first).at(sT.first), "Smooth"); } } }
          }
          // Apply to other variations
          else if ((lbl.first.find("_NTrackCorr")!=std::string::npos) || (lbl.first.find("BinWidth")!=std::string::npos)) {
            for (auto& sV : valMap) { for (auto& sT : sV.second) { if (sT.first.find("Err_Syst_")!=std::string::npos) { procSyst(valMap.at(sV.first).at(sT.first), "Smooth"); } } }
          }
          // Apply to No EWK Corr
          else if (lbl.first.find("_NoEWK")!=std::string::npos) {
            for (auto& sV : valMap) { for (auto& sT : sV.second) { if (sT.first.find("Err_Syst_")!=std::string::npos) { procSyst(valMap.at(sV.first).at(sT.first), "Smooth"); } } }
          }
	  for (const auto& sV : valMap) { for (const auto& sT : sV.second) { for (const auto& sB : sT.second) {
                auto& val = var[c.first][chg.first][sV.first][sT.first][sB.first];
                auto& nomVal = nomVar.at(c.first).at(chg.first).at(sV.first).at((sT.first=="Nom")?"Val":sT.first).at(sB.first);
		if (sT.first.find("Err_Syst_")!=std::string::npos) {
		  val = sumErrors(val, sB.second);
		  if (cat.first!="Efficiency") { nomVal = sumErrors(nomVal, sB.second); } // Needed, otherwise it double counts the TnP eff unc
		}
		else if (sT.first=="Nom") { val = nomVal; }
                else if (sT.first=="Val") { val = sB.second;
                  if (origVarVec.size()>1) { val = nomVal + std::max(var.at(c.first).at(chg.first).at(sV.first).at("Err_Syst_Low" ).at(sB.first),
                                                                     var.at(c.first).at(chg.first).at(sV.first).at("Err_Syst_High").at(sB.first)); }
                }
	      }
	    }
	  }
	}
      }
      systVarVec.at(lbl.first).push_back(systVar);
    }
  }
};


std::string formatResultVarName(const std::string varName, const bool useEtaCM, const bool isSyst = false, const bool useLATEX = false, const std::string chg="")
{
  std::string label = "";
  const std::string etaLbl = (useLATEX ? ( useEtaCM ? "\\etaMuCM" : "\\etaMuLAB" ) : ( useEtaCM ? "#eta^{#mu}_{CM}" : "#eta^{#mu}_{LAB}" ));
  if (isSyst) {
    if (varName == "Charge_Asymmetry"      ) {
      if (useLATEX) {
        label = "Abs. Unc. \\frac{( N_{\\mu}^{+} - N_{\\mu}^{-} )}{( N_{\\mu}^{+} + N_{\\mu}^{-} )}";
      }
      else {
        label = "Abs. Unc. ( N_{#mu}^{+} #font[122]{\55} N_{#mu}^{#font[122]{\55}} ) / ( N_{#mu}^{+} + N_{#mu}^{#font[122]{\55}} )";
      }
    }
    if (varName == "ForwardBackward_Ratio" ) {
      if (useLATEX) {
        if      (chg=="Pl" ) { label = Form("Abs. Unc.  \\frac{N_{\\mu}^{+}( +\\etaMuCM )}{N_{\\mu}^{+}( -\\etaMuCM )}"); }
        else if (chg=="Mi" ) { label = Form("Abs. Unc.  \\frac{N_{\\mu}^{-}( +\\etaMuCM )}{N_{\\mu}^{-}( -\\etaMuCM )}"); }
        else if (chg=="Inc") { label = Form("Abs. Unc.  \\frac{N_{\\mu}( +\\etaMuCM )}{N_{\\mu}( -\\etaMuCM )}");         }
        else { label = "Abs. Unc.  R_{FB}"; }
      }
      else {
        if      (chg=="Pl" ) { label = Form("Abs. Unc.  N_{#mu}^{+}(+#eta^{#mu}_{CM}) / N_{#mu}^{+}(#font[122]{\55}#eta^{#mu}_{CM})"); }
        else if (chg=="Mi" ) { label = Form("Abs. Unc.  N_{#mu}^{#font[122]{\55}}(+#eta^{#mu}_{CM}) / N_{#mu}^{#font[122]{\55}}(#font[122]{\55}#eta^{#mu}_{CM})"); }
        else if (chg=="Inc") { label = Form("Abs. Unc.  N_{#mu}(+#eta^{#mu}_{CM}) / N_{#mu}(#font[122]{\55}#eta^{#mu}_{CM})");         }
        else { label = "Abs. Unc.  R_{FB}"; }
      }
    }
    if (varName == "Cross_Section"         ) {
      if      (chg=="Pl") { label = Form("Rel. Unc. %s / d%s [nb]", (useLATEX ? "d\\sigma(\\WToMuNuPl)" : "d#sigma(W^{+}#kern[0.2]{#rightarrow}#kern[0.2]{#mu^{+}}#kern[0.2]{#nu_{#mu}})"), etaLbl.c_str()); }
      else if (chg=="Mi") { label = Form("Rel. Unc. %s / d%s [nb]", (useLATEX ? "d\\sigma(\\WToMuNuMi)" : "d#sigma(W^{#font[122]{\55}}#kern[0.2]{#rightarrow}#kern[0.2]{#mu^{#font[122]{\55}}}#kern[0.2]{#bar{#nu}_{#mu}})"), etaLbl.c_str()); }
      else { label = Form("Rel. Unc. B %s/d%s [nb]", (useLATEX ? "\\times d\\sigma" : "#times d#sigma"), etaLbl.c_str()); }
    }
  }
  else {
    if (varName == "Charge_Asymmetry"      ) {
      if (useLATEX) {
        label = "\\frac{( N_{\\mu}^{+} - N_{\\mu}^{-} )}{( N_{\\mu}^{+} + N_{\\mu}^{-} )}";
      }
      else {
        label = "( N_{#mu}^{+} #font[122]{\55} N_{#mu}^{#font[122]{\55}} ) / ( N_{#mu}^{+} + N_{#mu}^{#font[122]{\55}} )";
      }
    }
    if (varName == "ForwardBackward_Ratio" ) {
      if (useLATEX) {
        if      (chg=="Pl" ) { label = Form("\\frac{N_{\\mu}^{+}( +\\etaMuCM )}{ N_{\\mu}^{+}( -\\etaMuCM )}"); }
        else if (chg=="Mi" ) { label = Form("\\frac{N_{\\mu}^{-}( +\\etaMuCM )}{N_{\\mu}^{-}( -\\etaMuCM )}"); }
        else if (chg=="Inc") { label = Form("\\frac{N_{\\mu}( +\\etaMuCM )}{N_{\\mu}( -\\etaMuCM )}");         }
        else if (chg=="Com") { label = Form("\\frac{N_{\\mu}^{\\pm}( +\\etaMuCM )}{N_{\\mu}^{\\pm}( -\\etaMuCM )}"); }
        else { label = "R_{\\text{FB}}"; }
      }
      else {
        if      (chg=="Pl" ) { label = Form("N_{#mu}^{+}(+#eta^{#mu}_{CM}) / N_{#mu}^{+}(#font[122]{\55}#eta^{#mu}_{CM})"); }
        else if (chg=="Mi" ) { label = Form("N_{#mu}^{#font[122]{\55}}(+#eta^{#mu}_{CM}) / N_{#mu}^{#font[122]{\55}}(#font[122]{\55}#eta^{#mu}_{CM})"); }
        else if (chg=="Inc") { label = Form("N_{#mu}(+#eta^{#mu}_{CM}) / N_{#mu}(#font[122]{\55}#eta^{#mu}_{CM})");         }
        else if (chg=="Com") { label = Form("N_{#mu}^{#pm}(+#eta^{#mu}_{CM} ) / N_{#mu}^{#pm}(#font[122]{\55}#eta^{#mu}_{CM})"); }
        else { label = "R_{FB}"; }
      }
    }
    if (varName == "Cross_Section") {
      if      (chg=="Pl" ) { label = Form("%s / d%s [nb]", (useLATEX ? "d\\sigma(\\WToMuNuPl)" : "d#sigma(W^{+}#kern[0.2]{#rightarrow}#kern[0.2]{#mu^{+}}#kern[0.2]{#nu_{#mu}})"), etaLbl.c_str()); }
      else if (chg=="Mi" ) { label = Form("%s / d%s [nb]", (useLATEX ? "d\\sigma(\\WToMuNuMi)" : "d#sigma(W^{#font[122]{\55}}#kern[0.2]{#rightarrow}#kern[0.2]{#mu^{#font[122]{\55}}}#kern[0.2]{#bar{#nu}_{#mu}})"), etaLbl.c_str()); }
      else if (chg=="Inc") { label = Form("%s / d%s [nb]", (useLATEX ? "d\\sigma(\\WToMuNu)"   : "d#sigma(W#kern[0.2]{#rightarrow}#kern[0.2]{#mu}#kern[0.2]{#nu_{#mu}})"), etaLbl.c_str()); }
      else if (chg=="Com") { label = Form("%s / d%s", (useLATEX ? "d\\sigma(\\WToMuNu)"   : "d#sigma(W^{#pm}#kern[0.2]{#rightarrow}#kern[0.2]{#mu^{#pm}}#kern[0.2]{#nu_{#mu}})"), etaLbl.c_str()); }
      else { label = Form("B %s/d%s [nb]", (useLATEX ? "\\times d\\sigma" : "#times d#sigma"), etaLbl.c_str()); }
    }
    if (varName == "N_WToMu") { label = "Signal Yield"; }
  }
  return label;
};


bool iniResultsGraph(GraphPentaMap& graphMap, const BinSextaMapVec& var)
{
  //
  std::vector< std::string > nomGraphType = { "Err_Tot" , "Err_Stat" }; if (var.size()>1) { nomGraphType.push_back("Err_Syst"); }
  //
  for (const auto& c : var.begin()->second[0]) {
    for (const auto& ch : c.second) {
      for (const auto& v : ch.second) {
        const unsigned int nBins = v.second.begin()->second.size();
        for (const auto& lbl : var) {
          uint nLbl = lbl.second.size();
          if (nLbl>1 && ( (lbl.second[1].at(c.first).count(ch.first)==0) || (lbl.second[1].at(c.first).at(ch.first).count(v.first)==0) )) { nLbl = 2; }
          const uint nGraph = ( (lbl.first=="Nominal") ? nomGraphType.size() : nLbl );
          for (uint i = 0; i < nGraph; i++) {
            const std::string graphLbl = ( (lbl.first=="Nominal") ? nomGraphType[i] : ( (i<(nGraph-1)) ? Form("Variation_%d", i) : "Total" ) );
            auto& graph = graphMap[c.first][ch.first][v.first][lbl.first][graphLbl];
            graph.Set(nBins);
            // Set Graph Name
            const std::string name = Form("gr_WToMu%s_%s_%s_%s_%s", ch.first.c_str(), c.first.c_str(), v.first.c_str(), lbl.first.c_str(), graphLbl.c_str());
            graph.SetName(name.c_str());
          }
        }
      }
    }
  }
  // Return
  return true;
};


bool iniResultsGraph(GraphPentaMap& graphMap, const BinSextaMap& var)
{
  BinSextaMapVec varMap;
  for (const auto& cat : var) { varMap[cat.first].push_back(cat.second); varMap.at(cat.first).push_back(cat.second); }
  return iniResultsGraph(graphMap, varMap);
};


bool fillResultsGraph(GraphPentaMap& graphMap, const BinSextaMapVec& var)
{
  //
  std::cout << "[INFO] Filling the output graphs" << std::endl;
  std::vector< std::string > nomGraphType = { "Err_Tot" , "Err_Stat" }; if (var.size()>1) { nomGraphType.push_back("Err_Syst"); }
  //
  for (const auto& c : var.begin()->second[0]) {
    for (const auto& ch : c.second) {
      for (const auto& v : ch.second) {
	auto& graphMapMap = graphMap.at(c.first).at(ch.first).at(v.first);
	//
        // Determine index of bins
        auto& binMap = var.begin()->second[0].at(c.first).at(ch.first).at(v.first).begin()->second;
        std::map< anabin<0> , uint > binIdx;
        uint iBin = 0; for (const auto& b : binMap) { binIdx[b.first] = iBin; iBin++; }
        //
        // Compute systematic graphs
        for (const auto& lbl : var) {
          if (lbl.first=="Nominal") continue;
          uint nLbl = lbl.second.size();
          if (nLbl>1 && (var.at(lbl.first)[1].count(c.first)==0 || var.at(lbl.first)[1].at(c.first).count(ch.first)==0 || var.at(lbl.first)[1].at(c.first).at(ch.first).count(v.first)==0)) { nLbl = 2; }
          for (uint iGrP=0; iGrP<nLbl; iGrP++) {
            const uint iGr = ( (iGrP<(nLbl-1)) ? iGrP : (lbl.second.size()-1) ); // The last one in Var correspond to Total
            //
            // Compute the systematic variation graphs
            //
            if (var.count(lbl.first)==0) { std::cout << "[ERROR] Systematic result for " << lbl.first << " is missing" << std::endl; return false; }
            if (var.at(lbl.first).size()<lbl.second.size()) { std::cout << "[ERROR] Systematic result " << lbl.first << " is missing variations " << std::endl; return false; }
            if (var.at(lbl.first)[iGr].count(c.first)==0) { std::cout << "[ERROR] Systematic result " << lbl.first << " is missing " << c.first << std::endl; return false; }
            if (var.at(lbl.first)[iGr].at(c.first).count(ch.first)==0) { std::cout << "[ERROR] Systematic result " << lbl.first << " is missing " << ch.first << std::endl; return false; }
            if (var.at(lbl.first)[iGr].at(c.first).at(ch.first).count(v.first)==0) { std::cout << "[ERROR] Systematic result " << lbl.first << " is missing " << v.first << std::endl; return false; }
            if (var.at(lbl.first)[iGr].at(c.first).at(ch.first).at(v.first).count("Err_Syst_Low")==0) {
	      std::cout << "[ERROR] Systematic result " << lbl.first << "  " << v.first << " is missing the error" << std::endl; return false;
	    }
            if (var.at(lbl.first)[iGr].at(c.first).at(ch.first).at(v.first).count("Val")==0) {
	      std::cout << "[ERROR] Systematic result " << lbl.first << "  " << v.first << " is missing the value" << std::endl; return false;
	    }
            auto& vVar = var.at(lbl.first)[iGr].at(c.first).at(ch.first).at(v.first);
            auto& valBin = vVar.at("Val");
	    //
            const std::string graphLbl = ( (iGr<(lbl.second.size()-1)) ? Form("Variation_%d", iGr) : "Total" );
	    auto& graph = graphMapMap.at(lbl.first).at(graphLbl);
            //
            for (const auto& b : binMap) {
              //
              if (valBin.count(b.first)==0) {
                std::cout << "[ERROR] Systematic result " << lbl.first << " is missing bin [" << b.first.etabin().low() << " , " << b.first.etabin().high() << "]" << std::endl; return false;
              }
              //
              // Extract the parameters needed for each axis
              //
              // X Value
              const double X = ( (b.first.etabin().high() + b.first.etabin().low()) / 2.0 ); // Mean value of eta bin
              // X Error
              const double Err_X = ( (b.first.etabin().high() - b.first.etabin().low()) / 2.0 ); // Width of eta bin
              const double Err_X_High = Err_X;
              const double Err_X_Low  = Err_X;
              // Y Axis
              double Y = 0.0 , Err_Y_High = 0.0 , Err_Y_Low = 0.0;
	      if (graphLbl=="Total") {
                Y = vVar.at("Nom").at(b.first);
                Err_Y_High = vVar.at("Err_Syst_High").at(b.first);
		Err_Y_Low  = vVar.at("Err_Syst_Low").at(b.first);
	      }
	      else {
                Y = valBin.at(b.first);
		Err_Y_High = 0.0;
		Err_Y_Low  = 0.0;
	      }
              //
              // Fill the systematic variation graphs
              //
              const unsigned int iBin = binIdx.at(b.first);
              graph.SetPoint(iBin, X, Y);
              graph.SetPointError(iBin, Err_X_Low, Err_X_High, Err_Y_Low, Err_Y_High);
            }
          }
        }
        //
	// Compute nominal graph
	if (graphMapMap.count("Nominal")>0 && var.count("Nominal")>0) {
	  for (auto& gr : graphMapMap.at("Nominal")) {
	    auto& graph = gr.second;
	    auto& nomVar = var.at("Nominal")[0].at(c.first).at(ch.first).at(v.first);
	    for (const auto& b : binMap) {
	      const unsigned int iBin = binIdx.at(b.first);
	      //
	      // Extract the parameters needed for each axis
	      //
	      // X Value
	      const double X = ( (b.first.etabin().high() + b.first.etabin().low()) / 2.0 ); // Mean value of eta bin
	      // X Error
	      const double Err_X = ( (b.first.etabin().high() - b.first.etabin().low()) / 2.0 ); // Width of eta bin
	      double Err_X_High = Err_X;
	      double Err_X_Low  = Err_X;
	      // Y Value
	      const double Y = nomVar.at("Val").at(b.first);
	      //
	      // Compute total systematic error
	      const double Err_Y_Syst_High = nomVar.at("Err_Syst_High").at(b.first);
	      const double Err_Y_Syst_Low  = nomVar.at("Err_Syst_Low").at(b.first);
	      //
	      // Compute total statistic error
	      const double Err_Y_Stat_High = nomVar.at("Err_Stat_High").at(b.first);
	      const double Err_Y_Stat_Low  = nomVar.at("Err_Stat_Low" ).at(b.first);
	      //
	      // Y Error
	      double Err_Y_High = 0.0 , Err_Y_Low  = 0.0;
	      if (gr.first=="Err_Tot") {
		Err_Y_High = sumErrors( Err_Y_Stat_High , Err_Y_Syst_High );
		Err_Y_Low  = sumErrors( Err_Y_Stat_Low  , Err_Y_Syst_Low  );
	      }
	      if (gr.first=="Err_Stat") {
		Err_Y_High = Err_Y_Stat_High;
		Err_Y_Low  = Err_Y_Stat_Low;
	      }
	      if (gr.first=="Err_Syst") {
		Err_Y_High = Err_Y_Syst_High;
		Err_Y_Low  = Err_Y_Syst_Low;
	      }
	      //
	      // Fill the nominal graph
	      //
	      graph.SetPoint(iBin, X, Y);
	      graph.SetPointError(iBin, Err_X_Low, Err_X_High, Err_Y_Low, Err_Y_High);
	    }
	  }
	}
      }
    }
  }
  // Return
  return true;
};


bool fillResultsGraph(GraphPentaMap& graphMap, const BinSextaMap& var)
{
  BinSextaMapVec varMap;
  for (const auto& cat : var) { varMap[cat.first].push_back(cat.second); varMap.at(cat.first).push_back(cat.second); }
  return fillResultsGraph(graphMap, varMap);
};


void setStyle()
{
  // Set the CMS style
  setTDRStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TGaxis::SetMaxDigits(3); // to display powers of 10
  //
  // Set Palette
  gStyle->SetPalette(55);
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
  //
};


void formatLegendEntry(TLegendEntry& e, double size=0.040)
{
  e.SetTextSize(size);
};


void formatResultsGraph(TGraphAsymmErrors& graph, const std::string& col, const std::string& var, const std::string& chg, const bool& useEtaCM, const bool& incAcc, const bool isSyst = false, const int shiftEta=0)
{
  //
  // Set the Axis Titles
  std::string xLabel = "#eta^{#mu}";
  if (useEtaCM) { xLabel += "_{CM}"; }
  else { xLabel += "_{LAB}"; }
  if (var=="Cross_Section" && shiftEta>0) { xLabel = "#xi_{1} = (M^{W}/#sqrt{s_{NN}})#timese^{#eta^{#mu}_{CM}}";  }
  if (var=="Cross_Section" && shiftEta<0) { xLabel = "#xi_{2} = (M^{W}/#sqrt{s_{NN}})#timese^{-#eta^{#mu}_{CM}}"; }
  if (var=="Charge_Asymmetry" && shiftEta>0) { xLabel = "#eta^{#mu}_{ref}"; }
  std::string yLabel = formatResultVarName(var, useEtaCM, isSyst, false, chg);
  if (var=="Cross_Section" && shiftEta>0) { yLabel = "(#sqrt{s_{NN}} / GeV)^{-0.8} #times #xi_{1}d#sigma / d#xi_{1} [pb]"; }
  if (var=="Cross_Section" && shiftEta<0) { yLabel = "(#sqrt{s_{NN}} / GeV)^{-0.8} #times #xi_{2}d#sigma / d#xi_{2} [pb]"; }
  graph.SetTitle(Form(";%s;%s", xLabel.c_str(), yLabel.c_str()));
  //
  // General
  graph.SetMarkerColor(kBlack);
  graph.SetLineColor(kBlack); 
  graph.SetMarkerStyle(20);
  graph.SetMarkerSize(1.5);
  graph.SetLineWidth(3);
  graph.SetLineStyle(1);
  graph.SetFillStyle(0);
  // X-axis
  graph.GetXaxis()->CenterTitle(kTRUE);
  graph.GetXaxis()->SetTitleOffset(0.70);
  graph.GetXaxis()->SetTitleSize(0.065);
  if (var=="Cross_Section" && shiftEta!=0) {
    graph.GetXaxis()->SetTitleOffset(1.05);
    graph.GetXaxis()->SetTitleSize(0.050);
  }
  graph.GetXaxis()->SetLabelSize(0.035);
  double xMin=-2.5, xMax=2.5;
  if (useEtaCM) {
    if ( var == "Charge_Asymmetry"      ) { xMin = -3.0; xMax = 2.1; }
    if ( var == "Cross_Section"         ) { xMin = -3.0; xMax = 2.1; }
    if ( var == "ForwardBackward_Ratio" ) { xMin = -0.1; xMax = 2.0; }
  }
  if (shiftEta!=0) {
    if ( var == "Charge_Asymmetry"      ) { xMin = -3.5;  xMax = 2.5; }
    if ( var == "Cross_Section"         ) { xMin = 0.008; xMax = 0.4; }
  }
  graph.GetXaxis()->SetLimits(xMin , xMax);
  // Y-axis
  graph.GetYaxis()->CenterTitle(kTRUE);
  graph.GetYaxis()->SetTitleOffset(1.05);
  graph.GetYaxis()->SetTitleSize(0.065);
  graph.GetYaxis()->SetLabelSize(0.035);
  if ( var == "Charge_Asymmetry"      ) { graph.GetYaxis()->SetRangeUser(-0.12, 0.36); }
  if ( var == "ForwardBackward_Ratio" ) { graph.GetYaxis()->SetRangeUser( 0.50, 1.45); }
  if (incAcc){ if ( var == "Cross_Section" ) { graph.GetYaxis()->SetRangeUser(0.0, 300.0); } }
  else       { if ( var == "Cross_Section" ) { graph.GetYaxis()->SetRangeUser(65.0, 190.0); } }
  if ( var == "N_WToMu" ) { graph.GetYaxis()->SetRangeUser(0.0, 10000.); }
};


void setRangeYAxisGraph(TGraphAsymmErrors& graph, const double fracDo = 0.1, const double fracUp = 0.5)
{
  // Find maximum and minimum points of Plot to rescale Y axis
  double YMax = -1e99;
  for (int i=0; i<=graph.GetN(); i++) { double x, y; graph.GetPoint(i, x, y); y += graph.GetErrorY(i); YMax = std::max(YMax, y); }
  double YMin = 1e99;
  for (int i=0; i<=graph.GetN(); i++) { double x, y; graph.GetPoint(i, x, y); y -= graph.GetErrorY(i); YMin = std::min(YMin, y); }
  if (std::abs(YMax)>std::abs(YMin)) { YMin = -std::abs(YMax); } else { YMax = std::abs(YMin); }
  //
  const double YTot = (YMax - YMin)/(1.0 - fracUp - fracDo);
  const double Yup   = YMax + fracUp*YTot;
  const double Ydown = YMin - fracDo*YTot;
  //
  graph.GetYaxis()->SetRangeUser(Ydown, Yup);
};


void drawGraph( GraphPentaMap& graphMap , const std::string& outDir , const bool useEtaCM = true , const std::string accType = "MC" , const std::string effType = "TnP" )
{
  //
  // Set Style
  setStyle();
  //
  std::cout << "[INFO] Drawing the output graphs" << std::endl;
  //
  // Draw all graphs
  for (auto& c : graphMap) {
    for (auto& ch : c.second) {
      for (auto& v : ch.second) {
        for (auto& lbl : v.second) {
          //
          const std::string col  = c.first;
          const std::string chg  = ( (ch.first!="") ? ch.first : "Inc" );
          const std::string var  = v.first;
          const std::string type = lbl.first;
          auto& graph = lbl.second;
          if (std::count(var.begin(), var.end(), '_')>2) continue;
          //
          // Create Canvas
          TCanvas c("c", "c", 1000, 1000); c.cd();
          //
          // Create the Text Info
          TLatex tex; tex.SetNDC(); tex.SetTextSize(0.035); float dy = 0;
          std::vector< std::string > textToPrint;
          std::string sampleLabel = "W #rightarrow #mu + #nu_{#mu}";
          if (chg == "Pl") { sampleLabel = "W^{+} #rightarrow #mu^{+} + #nu_{#mu}"; }
          if (chg == "Mi") { sampleLabel = "W^{#font[122]{\55}} #rightarrow #mu^{#font[122]{\55}} + #bar{#nu}_{#mu}"; }
          textToPrint.push_back(sampleLabel);
          if (accType=="") { textToPrint.push_back("p^{#mu}_{T} > 25 GeV/c"); }
          //
          // Declare the graph vector (for drawing with markers)
          std::vector< TGraphAsymmErrors > grVec;
          // Initialize the Legend
          double legOff = 0.0 , legOff2 = 0.0; if (accType=="") { legOff = 0.05; }
          if ( type=="MC_Statistics" || type=="TnP_Syst_PU" || type=="TnP_Syst_STA" ) { legOff2 = 0.07; }
          TLegend leg(0.2, (0.68 - legOff + legOff2), 0.4, (0.81 - legOff));
          // Initialize the graph x range variables
          double xMin=0.0 , xMax=0.0 , yErrMin=9999999. , yErrMax=-1. , fbYLine = -1.;
          //
          std::string mainDir = "";
          // Draw graph
          if ( type=="Nominal" ) {
            mainDir = "Result";
            const bool incAcc = (accType!="");
            grVec.push_back(graph.at("Err_Stat"));
            if (graph.count("Err_Syst")>0) {
	      grVec.push_back(graph.at("Err_Syst"));
	      for (int j = 0; j < grVec.back().GetN(); j++) {
		grVec.back().SetPointError(j, grVec.back().GetErrorXlow(j)*0.4, grVec.back().GetErrorXhigh(j)*0.4, grVec.back().GetErrorYlow(j), grVec.back().GetErrorYhigh(j));
	      }
	      grVec.push_back(graph.at("Err_Tot"));
	      grVec.push_back(graph.at("Err_Tot"));
              for (int i=0; i<grVec[0].GetN(); i++) { double x, y; grVec[0].GetPoint(i, x, y); grVec[2].SetPoint(i, x, y+grVec[2].GetErrorYhigh(i)); grVec[3].SetPoint(i, x, y-grVec[3].GetErrorYlow(i)); }
	    }
            for (auto& gr : grVec) {
              formatResultsGraph(gr, col, var, chg, useEtaCM, incAcc);
              gr.SetFillStyle(1001);
            }
            if (graph.count("Err_Syst")>0) {
              for (uint j=2; j<=3; j++) {
                grVec[j].SetMarkerSize(0);
                for (int i=0; i<grVec[j].GetN(); i++) { grVec[j].SetPointEYhigh(i, 0.0); grVec[j].SetPointEYlow(i, 0.0); }
                for (int i=0; i<grVec[j].GetN(); i++) { grVec[j].SetPointEXhigh(i, 0.5*grVec[j].GetErrorXhigh(i)); grVec[j].SetPointEXlow(i, 0.5*grVec[j].GetErrorXlow(i)); }
              }
	    }
            for (int i=0; i<grVec[0].GetN(); i++) { grVec[0].SetPointEXhigh(i, 0.0); grVec[0].SetPointEXlow(i, 0.0); }
	    if (graph.count("Err_Syst")>0) { grVec[1].SetFillColor(kGreen+2); }
            // Create Legend
            //formatLegendEntry(*leg.AddEntry(&grVec[0], "Data", "pe"));
	    //if (graph.count("Err_Syst")>0) { formatLegendEntry(*leg.AddEntry(&grVec[2], "Systematic Uncertainty", "f"));  }
            // Draw the graphs
            grVec[0].Draw("ap");
            //if (graph.count("Err_Syst")>0) { grVec[1].Draw("same2"); }
            if (graph.count("Err_Syst")>0) { grVec[2].Draw("samep"); grVec[3].Draw("samep"); }
            grVec[0].Draw("samep");
            //
            xMin = grVec[0].GetXaxis()->GetXmin(); xMax = grVec[0].GetXaxis()->GetXmax(); fbYLine = 1.0;
          }
          else {
            mainDir = "Systematic";
            bool drawVar = true;
            if ( type=="MC_Statistics" || type=="TnP_Syst_PU" || type=="TnP_Syst_STA" ) { drawVar = false; }
            const auto ref = graph.at("Total");
            // Extract the Varied Efficiency graphs
            uint iMax = 0; double yMax = -9999999999.0;
            grVec.push_back(graph.at("Total"));
            for (auto gr : graph) { if (gr.first!="Total") { grVec.push_back(gr.second); } }
            for (uint i=0; i<grVec.size(); i++) {
              if (i==0) {
                for (int j = 0; j < grVec[i].GetN(); j++) {
                  double x, y; ref.GetPoint(j, x, y); grVec[i].SetPoint(j, x, 0.0);
                  double errLo = grVec[0].GetErrorYlow(j) , errHi = grVec[0].GetErrorYhigh(j);
                  if (var=="Cross_Section") { errLo /= y; errHi /= y; }
                  grVec[0].SetPointError(j, grVec[0].GetErrorXlow(j), grVec[0].GetErrorXhigh(j), errLo, errHi);
                  if (yMax < 0.90*std::max(errLo, errHi)) { yMax = std::max(errLo, errHi); }
                }
              }
              else if (drawVar) {
                for (int j = 0; j < grVec[i].GetN(); j++) {
                  double x=0., y1=0., y2=0.; ref.GetPoint(j, x, y1); grVec[i].GetPoint(j, x, y2);
                  double vVal = (y2 - y1);
                  if (var=="Cross_Section") { vVal /= y1; }
                  grVec[i].SetPoint(j, x, vVal);
                  grVec[i].SetPointError(j, grVec[i].GetErrorXlow(j)*0.5, grVec[i].GetErrorXhigh(j)*0.5, 0.0, 0.0);
                  if (yMax < 0.90*std::abs(vVal)) { yMax = std::abs(vVal); iMax = i; }
                }
              }
            }
            // Format the graphs
            const bool incAcc = (accType!="");
            for (auto& gr : grVec) { formatResultsGraph(gr, col, var, chg, useEtaCM, incAcc, true); }
            grVec[0].SetMarkerSize(0.0);
            setRangeYAxisGraph(grVec[iMax], 0.1, 0.4);
            std::vector<int> color;
            for (uint i=0; i<grVec.size(); i++) {  color.push_back(kRed); }
            if (grVec.size()==3) { color[1] = (kBlue+2); }
            for (uint i=0; i<grVec.size(); i++) { if (i!=0) { grVec[i].SetMarkerColor(kBlack); grVec[i].SetLineWidth(4); grVec[i].SetMarkerSize(0.0); grVec[i].SetLineColor(color[i]); } }
            // Create Legend
            formatLegendEntry(*leg.AddEntry(&grVec[0], Form("%s Uncertainty", type.c_str()), "f"), 0.028);
            if (drawVar) {
              if (grVec.size()==3) {
                formatLegendEntry(*leg.AddEntry(&grVec[2], Form("%s + Variation", type.c_str()), "l"), 0.028);
                formatLegendEntry(*leg.AddEntry(&grVec[1], Form("%s - Variation", type.c_str()), "l"), 0.028);
              }
              else { formatLegendEntry(*leg.AddEntry(&grVec[1], Form("%s Variation", type.c_str()), "l"), 0.028); }
            }
            // Draw the Graph
            grVec[iMax].Draw("apx");
            grVec[0].Draw("same2");
            if (drawVar) { for (uint i=0; i<grVec.size(); i++) { if (i!=0) { grVec[i].Draw("samep"); } } }
            if (iMax!=0) { grVec[0].Draw("same2"); }
            //
            xMin = grVec[0].GetXaxis()->GetXmin(); xMax = grVec[0].GetXaxis()->GetXmax(); fbYLine = 0.0;
            //
            for (int j = 0; j < ref.GetN(); j++) {
              double yErr = ref.GetErrorYhigh(j);
              if (var=="Cross_Section") { double x, yVal; ref.GetPoint(j, x, yVal); yErr /= std::abs(yVal); }
              if (yErrMax < yErr) { yErrMax = yErr; };  if (yErrMin > yErr) { yErrMin = yErr; }
            }
          }
          // Draw the Line
          TLine line_FB(xMin, fbYLine, xMax, fbYLine); line_FB.SetLineStyle(2);
          if (var=="ForwardBackward_Ratio") { line_FB.Draw("same"); }
          if (var=="Cross_Section" && type!="Nominal") { line_FB.Draw("same"); }
          TLine line_CA(xMin, 0.0, xMax, 0.0); line_CA.SetLineStyle(2);
          if (var=="Charge_Asymmetry") { line_CA.Draw("same"); }
          // Draw the Legend
          leg.Draw("same");
          // Update
          c.Modified(); c.Update();
          // Draw the text
          tex.SetTextSize(0.055); tex.DrawLatex(0.22, 0.84, textToPrint[0].c_str());
          tex.SetTextSize(0.058); tex.SetTextFont(61); tex.DrawLatex(0.78, 0.84, "CMS"); tex.SetTextFont(62);
          tex.SetTextSize(0.044); tex.SetTextFont(52); tex.DrawLatex(0.69, 0.79, "Preliminary"); tex.SetTextFont(62);
          if (textToPrint.size()>1) { tex.SetTextSize(0.040); tex.DrawLatex(0.22, 0.79, textToPrint[1].c_str()); }
          if ( type!="Nominal" )  {
            tex.SetTextSize(0.030);
            if (var=="Cross_Section") { tex.DrawLatex(0.20, 0.16, Form("Rel. Unc. [ %.2f%% , %.2f%% ]", (yErrMin*100.), (yErrMax*100.))); }
            else { tex.DrawLatex(0.20, 0.16, Form("Abs. Unc. [ %.2f%% , %.2f%% ]", (yErrMin*100.), (yErrMax*100.))); }
          }
          else {
            if (var=="Cross_Section") { tex.SetTextSize(0.030); tex.DrawLatex(0.25, 0.17, Form("Lumi. uncertainty (%.1f%%) not shown", LUMIUNC_*100.)); }
          }
          // Update
          c.Modified(); c.Update(); // Pure paranoia
          //
          // set the CMS style
          int option = 11830;
          if (col.find("pPb")!=std::string::npos) option = 115;
          if (col.find("Pbp")!=std::string::npos) option = 116;
          if (mainDir=="Systematic") { option = 117; }
          CMS_lumi(&c, option, 33, "", false, 0.6, false);
          // Update
          c.Modified(); c.Update(); // Pure paranoia
          //
          std::string label = "";
          if (accType==""   && effType==""   ) { label = "RAW";          }
          if (accType=="MC" && effType==""   ) { label = "AccMC";        }
          if (accType==""   && effType=="MC" ) { label = "EffMC";        }
          if (accType=="MC" && effType=="MC" ) { label = "AccMC_EffMC";  }
          if (accType==""   && effType=="TnP") { label = "EffTnP";       }
          if (accType=="MC" && effType=="TnP") { label = "AccMC_EffTnP"; }
          //
          // Create Output Directory
          std::string subDir = "";
          if (type=="Nominal") { subDir = "Nominal"; }
          else if (type.find("TnP_S")!=std::string::npos) { subDir = "TnP"; }
          else if (type.find("MC_Syst" )!=std::string::npos) { subDir = "PDF"; }
          else { subDir = "Variation"; }
          const std::string plotDir = outDir+"/Plots/" + mainDir+"/" + col+"/" + label+"/" + subDir+"/" + var;
          makeDir(plotDir + "/png/");
          makeDir(plotDir + "/pdf/");
          makeDir(plotDir + "/root/");
          //
          // Save Canvas
          const std::string name = Form("gr_WToMu%s_%s_%s_%s_%s", chg.c_str(), col.c_str(), var.c_str(), label.c_str(), type.c_str());
          c.SaveAs(( plotDir + "/png/"  + name + ".png"  ).c_str());
          c.SaveAs(( plotDir + "/pdf/"  + name + ".pdf"  ).c_str());
          c.SaveAs(( plotDir + "/root/" + name + ".root" ).c_str());
          //
          // Clean up memory
          c.Clear(); c.Close();
        }
      }
    }
  }
};


void drawSystematicGraph( const GraphPentaMap& graphMap , const std::string& outDir , const bool useEtaCM = true , const std::string accType = "MC" , const std::string effType = "TnP" )
{
  //
  // Set Style
  setStyle();
  //
  std::cout << "[INFO] Drawing the systematic graphs" << std::endl;
  //
  // Draw all graphs
  for (const auto& c : graphMap) {
    for (const auto& ch : c.second) {
      for (const auto& v : ch.second) {
        //
        // Declare the graph map
        std::map< std::string , TGraphAsymmErrors > grMap;
        for (auto& lbl : v.second) {
          if (lbl.first=="Nominal") continue;
          if (std::count(v.first.begin(), v.first.end(), '_')>2) continue;
          const auto& graph = lbl.second;
          grMap[lbl.first] = graph.at("Total");
        }
        if (v.second.count("Nominal")>0 && v.second.at("Nominal").count("Err_Syst")>0) { grMap["Systematic" ] = v.second.at("Nominal").at("Err_Syst"); }
        if (v.second.count("Nominal")>0 && v.second.at("Nominal").count("Err_Stat")>0) { grMap["Statistical"] = v.second.at("Nominal").at("Err_Stat"); }
        //
        for (auto& lbl : grMap) {
          if (lbl.first=="Nominal") continue;
          //
          const std::string col = c.first;
          const std::string chg = ( (ch.first!="") ? ch.first : "Inc" );
          const std::string var = v.first;
          const std::string type = lbl.first;
          auto& graph = lbl.second;
          //
          // Create Canvas
          TCanvas c("c", "c", 1000, 1000); c.cd();
          //
          // Create the Text Info
          TLatex tex; tex.SetNDC(); tex.SetTextSize(0.035); float dy = 0;
          std::vector< std::string > textToPrint;
          std::string sampleLabel = "W #rightarrow #mu + #nu_{#mu}";
          if (chg == "Pl") { sampleLabel = "W^{+} #rightarrow #mu^{+} + #nu_{#mu}"; }
          if (chg == "Mi") { sampleLabel = "W^{#font[122]{\55}} #rightarrow #mu^{#font[122]{\55}} + #bar{#nu}_{#mu}"; }
          textToPrint.push_back(sampleLabel);
          if (accType=="") { textToPrint.push_back("p^{#mu}_{T} > 25 GeV/c"); }
          //
          // Initialize the Legend
          double legOff = 0.0; if (accType=="") { legOff = 0.05; }
          TLegend leg(0.2, (0.76 - legOff), 0.4, (0.82 - legOff));
          double xMin = 0.0 , xMax = 0.0;
          //
          for (int j = 0; j < graph.GetN(); j++) {
            double X , Y; graph.GetPoint(j, X, Y); Y = std::abs(Y);
            graph.SetPoint(j, X, 0);
            //
            double Err_Y_Low  = graph.GetErrorYlow(j);
            double Err_Y_High = graph.GetErrorYhigh(j);
            if (var=="Cross_Section") { Err_Y_Low /= std::abs(Y); Err_Y_High /= std::abs(Y); }
            //
            graph.SetPointError(j,
                                graph.GetErrorXlow(j), graph.GetErrorXhigh(j),
                                Err_Y_Low, Err_Y_High
                                );
          }
          //
          // Proceed to draw the graphs
          //
          // Format the graphs
          const bool incAcc = (accType!="");
          formatResultsGraph(graph, col, var, chg, useEtaCM, incAcc, true);
          setRangeYAxisGraph(graph, 0.1, 0.4);
          // Create Legend
          formatLegendEntry(*leg.AddEntry(&graph, type.c_str(), "f"), 0.028);
          //
          graph.Draw("a2");
          //
          xMin = graph.GetXaxis()->GetXmin(); xMax = graph.GetXaxis()->GetXmax();
          //
          TLine line(xMin, 0.0, xMax, 0.0); line.SetLineStyle(2);
          line.Draw("same");
          //
          // Draw the Legend
          leg.Draw("same");
          // Update
          c.Modified(); c.Update();
          // Draw the text
          tex.SetTextSize(0.055); tex.DrawLatex(0.22, 0.84, textToPrint[0].c_str());
          tex.SetTextSize(0.058); tex.SetTextFont(61); tex.DrawLatex(0.78, 0.84, "CMS"); tex.SetTextFont(62);
          tex.SetTextSize(0.044); tex.SetTextFont(52); tex.DrawLatex(0.69, 0.80, "Preliminary"); tex.SetTextFont(62);
          if (textToPrint.size()>1) { tex.SetTextSize(0.035); tex.DrawLatex(0.22, 0.80, textToPrint[1].c_str()); }
          // Update
          c.Modified(); c.Update(); // Pure paranoia
          //
          // set the CMS style
          int option = 11830;
          if (col.find("pPb")!=std::string::npos) option = 115;
          if (col.find("Pbp")!=std::string::npos) option = 116;
          CMS_lumi(&c, option, 33, "", false, 0.6, false);
          // Update
          c.Modified(); c.Update(); // Pure paranoia
          //
          std::string label = "";
          if (accType==""   && effType==""   ) { label = "RAW";          }
          if (accType=="MC" && effType==""   ) { label = "AccMC";        }
          if (accType==""   && effType=="MC" ) { label = "EffMC";        }
          if (accType=="MC" && effType=="MC" ) { label = "AccMC_EffMC";  }
          if (accType==""   && effType=="TnP") { label = "EffTnP";       }
          if (accType=="MC" && effType=="TnP") { label = "AccMC_EffTnP"; }
          //
          // Create Output Directory
          std::string subDir = "";
          if (type=="Nominal") { subDir = "Nominal"; }
          else if (type.find("TnP_S")!=std::string::npos) { subDir = "TnP"; }
          else if (type.find("MC_Syst" )!=std::string::npos) { subDir = "PDF"; }
          else { subDir = "Variation"; }
          const std::string plotDir = outDir+"/Plots/Systematic/" + col+"/" + label+"/" + subDir+"/" + var;
          makeDir(plotDir + "/png/");
          makeDir(plotDir + "/pdf/");
          makeDir(plotDir + "/root/");
          //
          // Save Canvas
          const std::string name = Form("gr_WToMu%s_%s_%s_%s_%s", chg.c_str(), col.c_str(), var.c_str(), label.c_str(), type.c_str());
          c.SaveAs(( plotDir + "/png/"  + name + ".png"  ).c_str());
          c.SaveAs(( plotDir + "/pdf/"  + name + ".pdf"  ).c_str());
          c.SaveAs(( plotDir + "/root/" + name + ".root" ).c_str());
          //
          // Clean up memory
          c.Clear(); c.Close();
        }
      }
    }
  }
};

std::string frmLgn(std::string cl)
{
  if (cl=="Boson_pT") { return "Weak boson p_{T}"; }
  if (cl=="EWK_Background") { return "Electroweak and t#bar{t} backgrounds"; }
  if (cl=="EWK_Correction") { return "W-boson POWHEG BOX"; }
  if (cl=="Efficiency") { return "Signal efficiency"; }
  if (cl=="Event_Activity") { return "Event activity"; }
  if (cl=="QCD_Background") { return "QCD jet background"; }
  if (cl=="Recoil_Correction") { return "Recoil correction"; }
  if (cl=="Statistical") { return "Statistical uncertainty"; }
  if (cl=="Total_Systematic") { return "Total systematic uncertainty"; }
  return cl;
};

void drawCombineSystematicGraph( const GraphPentaMap& graphMap , const std::string& outDir , const bool useEtaCM = true ,
                                 const std::string accType = "MC" , const std::string effType = "TnP" , const std::string cmbType = "" )
{
  //
  // Set Style
  setStyle();
  //
  std::cout << "[INFO] Drawing the combine systematic graphs" << std::endl;
  //
  // Draw all graphs
  for (auto& c : graphMap) {
    for (auto& ch : c.second) {
      for (auto& v : ch.second) {
        //
        const std::string col = c.first;
        const std::string chg = ( (ch.first!="") ? ch.first : "Inc" );
        const std::string var = v.first;
        //
        // Create Canvas
        TCanvas c("c", "c", 1000, 1000); c.cd();
        //
        // Create the Text Info
        TLatex tex; tex.SetNDC(); tex.SetTextSize(0.035); float dy = 0;
        std::vector< std::string > textToPrint;
        std::string sampleLabel = "W #rightarrow #mu + #nu_{#mu}";
        if (chg == "Pl") { sampleLabel = "W^{+} #rightarrow #mu^{+} + #nu_{#mu}"; }
        if (chg == "Mi") { sampleLabel = "W^{#font[122]{\55}} #rightarrow #mu^{#font[122]{\55}} + #bar{#nu}_{#mu}"; }
        textToPrint.push_back(sampleLabel);
        if (accType=="") { textToPrint.push_back("p^{#mu}_{T} > 25 GeV/c"); }
        //
        // Declare the graph map (for drawing with markers)
        std::map< std::string , TGraphAsymmErrors > grMapTmP , grMap;
        //
        for (const auto& lbl : v.second) {
          if (lbl.first=="Nominal") continue;
          if (cmbType!="" && lbl.first.find(cmbType)==std::string::npos) continue;
          const std::string type = lbl.first;
          auto& graph = lbl.second;
          if (graph.count("Total")>0) { grMapTmP[lbl.first] = graph.at("Total"); }
        }
        if (v.second.count("Nominal")>0) {
          grMapTmP["Total_Systematic" ] = v.second.at("Nominal").at("Err_Syst");
          grMapTmP["Statistical"] = v.second.at("Nominal").at("Err_Stat");
        }
        //
        std::string lMax = ""; double vMax = -99999999.0;
        for (const auto& lbl : grMapTmP) {
          //
          std::string vLbl = lbl.first;
          if (std::count(lbl.first.begin(), lbl.first.end(), '_')>2) {
            std::string tmp = vLbl.substr(vLbl.find("_")+1) , tmp2 = tmp.substr(tmp.find("_")+1);
            vLbl = vLbl.substr(0, vLbl.find("_"))+"_"+tmp.substr(0, tmp.find("_"))+"_"+tmp2.substr(0, tmp2.find("_"));
          }
          //
          if (grMap.count(vLbl)==0) { grMap[vLbl].Set(lbl.second.GetN()); }
          for (int j = 0; j < lbl.second.GetN(); j++) {
            double X , Y; lbl.second.GetPoint(j, X, Y); Y = std::abs(Y);
            grMap.at(vLbl).SetPoint(j, X, 0);
            //
            double Err_Y_Low  = lbl.second.GetErrorYlow(j);
            double Err_Y_High = lbl.second.GetErrorYhigh(j);
            if (var=="Cross_Section") { Err_Y_Low /= std::abs(Y); Err_Y_High /= std::abs(Y); }
            Err_Y_Low  = std::sqrt( std::pow(Err_Y_Low , 2.0) + std::pow(grMap.at(vLbl).GetErrorYlow(j) , 2.0) );
            Err_Y_High = std::sqrt( std::pow(Err_Y_High, 2.0) + std::pow(grMap.at(vLbl).GetErrorYhigh(j), 2.0) );
            //
            grMap.at(vLbl).SetPointError(j, 
                                         lbl.second.GetErrorXlow(j), lbl.second.GetErrorXhigh(j), 
                                         Err_Y_Low, Err_Y_High
                                         );
            if (vMax < std::max(Err_Y_Low, Err_Y_High)) { vMax = std::max(Err_Y_Low, Err_Y_High); lMax = vLbl; }
          }
        }
        // Initialize the Legend
        double legOff = 0.0 , legOff2 = 0.0; if (accType=="") { legOff = 0.05; }
        if (grMap.size()==3) { legOff2 = 0.15; }
        TLegend leg(0.2, (0.54 - legOff + legOff2), 0.4, (0.82 - legOff));
        double xMin = 0.0 , xMax = 0.0;
        //
        // Draw the graphs in the same canvas
        bool firstPlot = true; uint iCnt = 1;
        std::vector<int> COLOR = { kBlack , kRed , kBlue , kGreen , kOrange , kViolet , kGreen+2 , kAzure-7 , kRed+3 , kBlue-2 };
        for (auto& gr : grMap) {
          const auto& grLbl = gr.first;
          auto& graph = gr.second;
          // Format the graphs
          const bool incAcc = (accType!="");
          formatResultsGraph(graph, col, var, chg, useEtaCM, incAcc, true);
          graph.SetMarkerSize(0.0);
          graph.SetLineColor((iCnt==uint(kWhite) ? uint(kOrange) : COLOR[iCnt-1]));
          if (var == "Charge_Asymmetry"     ) { graph.GetYaxis()->SetRangeUser(-0.03, 0.10); }
          if (var == "Cross_Section"        ) { graph.GetYaxis()->SetRangeUser(-0.1,  0.3); }
          if (var == "ForwardBackward_Ratio") { graph.GetYaxis()->SetRangeUser(-0.05,  0.1); }
          iCnt++;
          // Create Legend
          formatLegendEntry(*leg.AddEntry(&graph, frmLgn(grLbl).c_str(), "f"), 0.028);
          if (grLbl==lMax) { setRangeYAxisGraph(graph, 0.02, 0.55); }
          //
          if (firstPlot) { grMap.at(lMax).Draw("apx"); firstPlot = false; }
          graph.Draw("same2");
          //
          xMin = graph.GetXaxis()->GetXmin(); xMax = graph.GetXaxis()->GetXmax();
        }
        //
        TLine line(xMin, 0.0, xMax, 0.0); line.SetLineStyle(2);
        line.Draw("same");
        //
        // Draw the Legend
        leg.Draw("same");
        // Update
        c.Modified(); c.Update();
        // Draw the text
        tex.SetTextSize(0.055); tex.DrawLatex(0.22, 0.84, textToPrint[0].c_str());
        tex.SetTextSize(0.058); tex.SetTextFont(61); tex.DrawLatex(0.78, 0.84, "CMS"); tex.SetTextFont(62);
        tex.SetTextSize(0.044); tex.SetTextFont(52); tex.DrawLatex(0.69, 0.79, "Preliminary"); tex.SetTextFont(62);
        if (textToPrint.size()>1) { tex.SetTextSize(0.030); tex.DrawLatex(0.22, 0.79, textToPrint[1].c_str()); }
        // Update
        c.Modified(); c.Update(); // Pure paranoia
        //
        // set the CMS style
        int option = 11888;
        if (col.find("pPb")!=std::string::npos) option = 115;
        if (col.find("Pbp")!=std::string::npos) option = 116;
        CMS_lumi(&c, option, 33, "", false, 0.6, false);
        // Update
        c.Modified(); c.Update(); // Pure paranoia
        //
        std::string label = "";
        if (accType==""   && effType==""   ) { label = "RAW";          }
        if (accType=="MC" && effType==""   ) { label = "AccMC";        }
        if (accType==""   && effType=="MC" ) { label = "EffMC";        }
        if (accType=="MC" && effType=="MC" ) { label = "AccMC_EffMC";  }
        if (accType==""   && effType=="TnP") { label = "EffTnP";       }
        if (accType=="MC" && effType=="TnP") { label = "AccMC_EffTnP"; }
        //
        // Create Output Directory
        std::string subDir = "Combined";
        const std::string plotDir = outDir+"/Plots/Systematic/" + col+"/" + label+"/" + subDir+"/" + var;
        makeDir(plotDir + "/png/");
        makeDir(plotDir + "/pdf/");
        makeDir(plotDir + "/root/");
        //
        // Save Canvas
        if (cmbType!="") { label += ( "_" + cmbType ); }
        const std::string name = Form("gr_WToMu%s_%s_%s_%s", chg.c_str(), col.c_str(), var.c_str(), label.c_str());
        c.SaveAs(( plotDir + "/png/"  + name + ".png"  ).c_str());
        c.SaveAs(( plotDir + "/pdf/"  + name + ".pdf"  ).c_str());
        c.SaveAs(( plotDir + "/root/" + name + ".root" ).c_str());
        //
        // Clean up memory
        c.Clear(); c.Close();
      }
    }
  }
};


void drawGraphWithTheoryAndRatio( GraphPentaMap& graphMap , const std::string& outDir , const bool useEtaCM = true , const std::string accType = "MC" , const std::string effType = "TnP" , const bool forPaper = true , const bool addLumiBand = false )
{
  //
  // Set Style
  setStyle();
  //
  std::cout << "[INFO] Drawing the output graphs with Theory predictions" << std::endl;
  //
  // Add the Theory Predictions
  A1m(graphMap["PA"]["Mi"]["ForwardBackward_Ratio"]["Theory"]);
  A1p(graphMap["PA"]["Pl"]["ForwardBackward_Ratio"]["Theory"]);
  A3 (graphMap["PA"][  ""]["ForwardBackward_Ratio"]["Theory"]);
  Wp (graphMap["PA"]["Pl"]["Cross_Section"]["Theory"]);
  Wm (graphMap["PA"]["Mi"]["Cross_Section"]["Theory"]);
  chasym(graphMap["PA"][""]["Charge_Asymmetry"]["Theory"]);
  //
  // Draw all graphs
  for (const auto& c : graphMap) {
    for (const auto& ch : c.second) {
      for (const auto& v : ch.second) {
        if (v.first!="Cross_Section") continue;
	//
	const std::string col  = c.first;
	const std::string chg  = ( (ch.first!="") ? ch.first : "Inc" );
	const std::string var  = v.first;
	auto& nomGraph = graphMap.at("PA").at(ch.first).at(v.first).at("Nominal");
	auto& theGraph = graphMap.at("PA").at(ch.first).at(v.first).at("Theory");
	//
	// Create Canvas
	TCanvas c("c", "c", 1000, 1000); c.cd();
	//
	// Create the Text Info
	TLatex tex; tex.SetNDC(); tex.SetTextSize(0.045); float dy = 0;
	std::vector< std::string > textToPrint;
	std::string sampleLabel = "W#kern[0.2]{#rightarrow}#kern[0.2]{#mu}#kern[0.2]{#nu_{#mu}}";
	if (chg == "Pl") { sampleLabel = "W^{+}#kern[0.2]{#rightarrow}#kern[0.2]{#mu^{+}}#kern[0.2]{#nu_{#mu}}"; }
	if (chg == "Mi") { sampleLabel = "W^{#font[122]{\55}}#kern[0.2]{#rightarrow}#kern[0.2]{#mu^{#font[122]{\55}}}#kern[0.2]{#bar{#nu}_{#mu}}"; }
	textToPrint.push_back(sampleLabel);
	if (accType=="") { textToPrint.push_back("p^{#mu}_{T} > 25 GeV/c"); }
	//
	// Declare the graph vector (for drawing with markers)
	std::vector< TGraphAsymmErrors > grVec, grRatio;
	// Initialize the Legend
	double legOff = 0.0; if (accType=="") { legOff = 0.05; }
        double legXmin = 0.2 , legXmax = 0.4 , legYmin = (0.52 - legOff) , legYmax = (0.76 - legOff);
        if (var=="Cross_Section" && chg=="Pl") { legXmin = 0.35; legXmax = 0.55; legYmin = 0.12; legYmax = 0.36; }
        if (var=="ForwardBackward_Ratio") { legXmin = 0.2; legXmax = 0.4; legYmin = 0.15; legYmax = 0.39; }
	TLegend leg(legXmin, legYmin, legXmax, legYmax);
	// Initialize the graph x range variables
	double xMin=0.0 , xMax=0.0 , yErrMin=9999999. , yErrMax=-1.;
	// Format the graphs
	const bool incAcc = (accType!="");
	grVec.push_back(nomGraph.at("Err_Stat"));
	grVec.push_back(nomGraph.at("Err_Tot"));
	grVec.push_back(nomGraph.at("Err_Tot"));
	grVec.push_back(theGraph.at("CT14"));
	grVec.push_back(theGraph.at("EPPS16"));
	grVec.push_back(theGraph.at("nCTEQ15"));
        if (addLumiBand && var=="Cross_Section") { grVec.push_back(nomGraph.at("Err_Stat")); }
	for (auto& gr : grVec) { formatResultsGraph(gr, col, var, chg, useEtaCM, incAcc); }
        grVec[0].GetXaxis()->SetTitle("");
        grVec[0].GetXaxis()->SetLabelOffset(3);
        grVec[0].GetYaxis()->SetLabelSize(0.050);
        grVec[0].GetYaxis()->SetTitleSize(0.075); //0.07
        grVec[0].GetYaxis()->SetTitleOffset(0.9);
        //
        for (int i=0; i<grVec[0].GetN(); i++) { grVec[0].SetPointEXhigh(i, 0.0); grVec[0].SetPointEXlow(i, 0.0); }
        for (int i=0; i<grVec[0].GetN(); i++) { double x, y; grVec[0].GetPoint(i, x, y); grVec[1].SetPoint(i, x, y+grVec[1].GetErrorYhigh(i)); grVec[2].SetPoint(i, x, y-grVec[2].GetErrorYlow(i)); }
	grVec[1].SetMarkerSize(0);
        grVec[1].SetLineWidth(3);
        for (int i=0; i<grVec[1].GetN(); i++) { grVec[1].SetPointEYhigh(i, 0.0); grVec[1].SetPointEYlow(i, 0.0); }
        for (int i=0; i<grVec[1].GetN(); i++) { grVec[1].SetPointEXhigh(i, 0.5*grVec[1].GetErrorXhigh(i)); grVec[1].SetPointEXlow(i, 0.5*grVec[1].GetErrorXlow(i)); }
	grVec[2].SetMarkerSize(0);
        grVec[2].SetLineWidth(3);
        for (int i=0; i<grVec[2].GetN(); i++) { grVec[2].SetPointEYhigh(i, 0.0); grVec[2].SetPointEYlow(i, 0.0); }
        for (int i=0; i<grVec[2].GetN(); i++) { grVec[2].SetPointEXhigh(i, 0.5*grVec[2].GetErrorXhigh(i)); grVec[2].SetPointEXlow(i, 0.5*grVec[2].GetErrorXlow(i)); }
	grVec[3].SetFillColor(kYellow);
	grVec[3].SetLineColor(kRed);
	grVec[3].SetFillStyle(1001);
	grVec[3].SetLineWidth(4);
	grVec[3].SetMarkerSize(0);
	grVec[4].SetFillColor(kGreen+2);
	grVec[4].SetLineColor(kGreen+2);
	grVec[4].SetFillStyle(3275);
	grVec[4].SetLineStyle(7);
	grVec[4].SetLineWidth(4);
	grVec[4].SetMarkerSize(0);
	grVec[5].SetFillColor(kOrange+2);
	grVec[5].SetLineColor(kOrange+2);
	grVec[5].SetFillStyle(3257);
	grVec[5].SetLineStyle(7);
	grVec[5].SetLineWidth(4);
	grVec[5].SetMarkerSize(0);
        gStyle->SetHatchesSpacing(1.9);
        gStyle->SetHatchesLineWidth(2);
        if (addLumiBand && var=="Cross_Section") {
          grVec.back().SetLineColor(kAzure-7);
          grVec.back().SetLineWidth(4);
          for (int i=0; i<grVec.back().GetN(); i++) { double x,y; grVec.back().GetPoint(i, x, y); grVec.back().SetPointEYhigh(i, y*LUMIUNC_); grVec.back().SetPointEYlow(i, y*LUMIUNC_); }
        }
        //
        auto h1 = graphToHist(grVec[3]); for (int i=1; i<=grVec[3].GetN(); i++) { h1.SetBinError(i, 0.0002); }
        auto h2 = graphToHist(grVec[4]); for (int i=1; i<=grVec[4].GetN(); i++) { h2.SetBinError(i, 0.0002); }
        auto h3 = graphToHist(grVec[5]); for (int i=1; i<=grVec[5].GetN(); i++) { h3.SetBinError(i, 0.0002); }
        //
	// Create Legend
	formatLegendEntry(*leg.AddEntry(&grVec[0], "Data", "pe"), 0.050);
	formatLegendEntry(*leg.AddEntry(&grVec[3], "CT14 (68% CL)", "lf"), 0.050);
	formatLegendEntry(*leg.AddEntry(&grVec[4], "CT14+EPPS16 (68% CL)", "lf"), 0.050);
	formatLegendEntry(*leg.AddEntry(&grVec[5], "CT14+nCTEQ15 (68% CL)", "lf"), 0.050);
        //leg.SetHeader("NLO MCFM + NLO PDF, 68% CL");
        //TLegendEntry *header = (TLegendEntry*)leg.GetListOfPrimitives()->First();
        //header->SetTextSize(0.05);
        //
        // Define the plotting pads
        TPad *pad1  = new TPad("pad1", "", 0, 0.25, 1, 0.98);  // Unique Pointer does produce Segmentation Fault, so don't use it
        pad1->SetBottomMargin(0.00);
        //
        // Draw the Graphs
        //
        // Main Frame
        pad1->Draw();
        pad1->cd();
        //
	grVec[0].Draw("ap");
	grVec[3].Draw("same2"); h1.Draw("sameL");
	grVec[4].Draw("same2"); h2.Draw("sameL");
	grVec[5].Draw("same2"); h3.Draw("sameL");
        if (addLumiBand && var=="Cross_Section") { grVec.back().Draw("same2"); }
	grVec[1].Draw("samep"); grVec[2].Draw("samep");
	grVec[0].Draw("samep");
	//
	xMin = grVec[0].GetXaxis()->GetXmin(); xMax = grVec[0].GetXaxis()->GetXmax();
	// Draw the Line
	TLine line_FB(xMin, 1.0, xMax, 1.0); line_FB.SetLineStyle(2);
	if (var=="ForwardBackward_Ratio") { line_FB.Draw("same"); }
	TLine line_CA(xMin, 0.0, xMax, 0.0); line_CA.SetLineStyle(2);
	if (var=="Charge_Asymmetry") { line_CA.Draw("same"); }
	// Draw the Legend
	leg.Draw("same");
	// Update
	c.Modified(); c.Update();
	// Draw the text
	tex.SetTextSize(0.055*1.33); tex.DrawLatex(0.22, 0.82, textToPrint[0].c_str());
        tex.SetTextSize(0.058*1.33); tex.SetTextFont(61); tex.DrawLatex(0.78, 0.82, "CMS"); tex.SetTextFont(62);
        if (!forPaper) { tex.SetTextSize(0.044*1.33); tex.SetTextFont(52); tex.DrawLatex(0.69, 0.75, "Preliminary"); tex.SetTextFont(62); }
	if (textToPrint.size()>1) { tex.SetTextSize(0.035*1.33); tex.DrawLatex(0.22, 0.75, textToPrint[1].c_str()); }
        if (!addLumiBand && var=="Cross_Section") { tex.SetTextSize(0.030*1.33); tex.DrawLatex(0.30, 0.07, Form("Lumi. uncertainty (%.1f%%) not shown", LUMIUNC_*100.)); }
	// Update
	c.Modified(); c.Update(); // Pure paranoia
	//
	// set the CMS style
	int option = 11832;
	if (col.find("pPb")!=std::string::npos) option = 11830;
	if (col.find("Pbp")!=std::string::npos) option = 11831;
        CMS_lumi(pad1, option, 33, "", false, 0.8, false);
        pad1->SetFillStyle(4000); 
        pad1->SetFrameFillStyle(4000);
	// Update
	c.Modified(); c.Update(); // Pure paranoia
        //
        TPad *pad2  = new TPad("pad2", "", 0, 0, 1, 0.25); // Unique Pointer does produce Segmentation Fault, so don't use it
        pad2->SetTopMargin(0.00);
        pad2->SetBottomMargin(0.4);
        pad2->SetFillStyle(4000); 
        pad2->SetFrameFillStyle(4000);
        //
	grRatio.push_back(nomGraph.at("Err_Stat"));
	grRatio.push_back(nomGraph.at("Err_Tot"));
	grRatio.push_back(nomGraph.at("Err_Tot"));
	grRatio.push_back(theGraph.at("CT14"));
	grRatio.push_back(theGraph.at("EPPS16"));
	grRatio.push_back(theGraph.at("nCTEQ15"));
        for (auto& gr : grRatio) { gr.SetName(Form("%s_tmp", gr.GetName())); }
        for (int i=0; i<grRatio[3].GetN(); i++) {
          double x, y; grRatio[3].GetPoint(i, x, y);
          const double yH = grRatio[3].GetErrorYhigh(i);
          const double yL = grRatio[3].GetErrorYlow(i);
          for (uint iGr=0; iGr<grRatio.size(); iGr++) {
            auto& gr = grRatio[iGr];
            double x1, y1; gr.GetPoint(i, x1, y1); gr.SetPoint(i, x1, y1/y);
            const double yH1 = std::sqrt(std::pow(gr.GetErrorYhigh(i)/y1 , 2.0) + std::pow(yH/y , 2.0))*std::abs(y1/y);
            const double yL1 = std::sqrt(std::pow(gr.GetErrorYlow (i)/y1 , 2.0) + std::pow(yL/y , 2.0))*std::abs(y1/y);
            gr.SetPointEYhigh(i, yH1); gr.SetPointEYlow(i, yL1);
          }
        }
	for (auto& gr : grRatio) { formatResultsGraph(gr, col, var, chg, useEtaCM, incAcc); }
	grRatio[0].SetMarkerSize(1.0);
        grRatio[0].GetXaxis()->SetTitleOffset(0.6);
        grRatio[0].GetXaxis()->SetTitleSize(0.24);
        grRatio[0].GetXaxis()->SetLabelSize(0.14);
        grRatio[0].GetYaxis()->SetTitle("");
        c.cd(); tex.SetTextSize(0.032); tex.SetTextAngle(90); tex.DrawLatex(0.06, 0.07, "Ratio to CT14"); tex.SetTextAngle(0);
        grRatio[0].GetYaxis()->SetTitleOffset(0.55);
        grRatio[0].GetYaxis()->SetTitleSize(0.105);
        grRatio[0].GetYaxis()->SetLabelSize(0.14);
        grRatio[0].GetYaxis()->SetRangeUser(0.68, 1.26);
        grRatio[0].GetYaxis()->SetNdivisions(404);
        //
        for (int i=0; i<grRatio[0].GetN(); i++) { grRatio[0].SetPointEXhigh(i, 0.0); grRatio[0].SetPointEXlow(i, 0.0); }
        for (int i=0; i<grRatio[0].GetN(); i++) { double x, y; grRatio[0].GetPoint(i, x, y); grRatio[1].SetPoint(i, x, y+grRatio[1].GetErrorYhigh(i)); grRatio[2].SetPoint(i, x, y-grRatio[2].GetErrorYlow(i)); }
	grRatio[1].SetMarkerSize(0);
        for (int i=0; i<grRatio[1].GetN(); i++) { grRatio[1].SetPointEYhigh(i, 0.0); grRatio[1].SetPointEYlow(i, 0.0); }
        for (int i=0; i<grRatio[1].GetN(); i++) { grRatio[1].SetPointEXhigh(i, 0.5*grRatio[1].GetErrorXhigh(i)); grRatio[1].SetPointEXlow(i, 0.5*grRatio[1].GetErrorXlow(i)); }
	grRatio[2].SetMarkerSize(0);
        for (int i=0; i<grRatio[2].GetN(); i++) { grRatio[2].SetPointEYhigh(i, 0.0); grRatio[2].SetPointEYlow(i, 0.0); }
        for (int i=0; i<grRatio[2].GetN(); i++) { grRatio[2].SetPointEXhigh(i, 0.5*grRatio[2].GetErrorXhigh(i)); grRatio[2].SetPointEXlow(i, 0.5*grRatio[2].GetErrorXlow(i)); }
	grRatio[3].SetFillColor(kYellow);
	grRatio[3].SetLineColor(kRed);
	grRatio[3].SetFillStyle(1001);
	grRatio[3].SetLineWidth(4);
	grRatio[3].SetMarkerSize(0);
	grRatio[4].SetFillColor(kGreen+2);
	grRatio[4].SetLineColor(kGreen+2);
	grRatio[4].SetFillStyle(3275);
	grRatio[4].SetLineStyle(7);
	grRatio[4].SetLineWidth(4);
	grRatio[4].SetMarkerSize(0);
	grRatio[5].SetFillColor(kOrange+2);
	grRatio[5].SetLineColor(kOrange+2);
	grRatio[5].SetFillStyle(3257);
	grRatio[5].SetLineStyle(7);
	grRatio[5].SetLineWidth(4);
	grRatio[5].SetMarkerSize(0);
        gStyle->SetHatchesSpacing(1.9);
        gStyle->SetHatchesLineWidth(2);
        //
        auto h1Ratio = graphToHist(grRatio[3]); for (int i=1; i<=grRatio[3].GetN(); i++) { h1Ratio.SetBinError(i, 0.0002); }
        auto h2Ratio = graphToHist(grRatio[4]); for (int i=1; i<=grRatio[4].GetN(); i++) { h2Ratio.SetBinError(i, 0.0002); }
        auto h3Ratio = graphToHist(grRatio[5]); for (int i=1; i<=grRatio[5].GetN(); i++) { h3Ratio.SetBinError(i, 0.0002); }
        //
        // Ratio Frame
        c.cd();
        pad2->Draw();
        pad2->cd();
        //
	grRatio[0].Draw("ap");
	grRatio[3].Draw("same2"); h1Ratio.Draw("sameL");
	grRatio[4].Draw("same2"); h2Ratio.Draw("sameL");
	grRatio[5].Draw("same2"); h3Ratio.Draw("sameL");
	grRatio[1].Draw("samep"); grRatio[2].Draw("samep");
	grRatio[0].Draw("samep");
        //
        line_FB.Draw("same");
        pad2->Update();
	//
	std::string label = "";
	if (accType==""   && effType==""   ) { label = "RAW";          }
	if (accType=="MC" && effType==""   ) { label = "AccMC";        }
	if (accType==""   && effType=="MC" ) { label = "EffMC";        }
	if (accType=="MC" && effType=="MC" ) { label = "AccMC_EffMC";  }
	if (accType==""   && effType=="TnP") { label = "EffTnP";       }
	if (accType=="MC" && effType=="TnP") { label = "AccMC_EffTnP"; }
	//
	// Create Output Directory
	const std::string plotDir = outDir+"/Plots/Theory/" + (forPaper?"ForPaper/":"ForPAS/") + col+"/" + label+"/" + var;
	makeDir(plotDir + "/png/");
	makeDir(plotDir + "/pdf/");
	makeDir(plotDir + "/root/");
	//
	// Save Canvas
	const std::string name = Form("gr_WToMu%s_%s_%s_%s_%s", chg.c_str(), col.c_str(), var.c_str(), label.c_str(), "NominalWithTheoryAndRatio");
	c.SaveAs(( plotDir + "/png/"  + name + ".png"  ).c_str());
	c.SaveAs(( plotDir + "/pdf/"  + name + ".pdf"  ).c_str());
	c.SaveAs(( plotDir + "/root/" + name + ".root" ).c_str());
	//
	// Clean up memory
	c.Clear(); c.Close();
      }
    }
  }
};


void drawGraphWithTheory( GraphPentaMap& graphMap , const std::string& outDir , const std::string type="EPPS16", const bool useEtaCM = true , const std::string accType = "MC" , const std::string effType = "TnP" , const bool forPaper = true )
{
  //
  // Set Style
  setStyle();
  //
  std::cout << "[INFO] Drawing the output graphs with Theory predictions" << std::endl;
  //
  // Add the Theory Predictions
  A1m(graphMap["PA"]["Mi"]["ForwardBackward_Ratio"]["Theory"]);
  A1p(graphMap["PA"]["Pl"]["ForwardBackward_Ratio"]["Theory"]);
  A3 (graphMap["PA"][  ""]["ForwardBackward_Ratio"]["Theory"]);
  Wp (graphMap["PA"]["Pl"]["Cross_Section"]["Theory"]);
  Wm (graphMap["PA"]["Mi"]["Cross_Section"]["Theory"]);
  chasym(graphMap["PA"][""]["Charge_Asymmetry"]["Theory"]);
  //
  // Draw all graphs
  for (const auto& c : graphMap) {
    for (const auto& ch : c.second) {
      for (const auto& v : ch.second) {
	//
	const std::string col  = c.first;
	const std::string chg  = ( (ch.first!="") ? ch.first : "Inc" );
	const std::string var  = v.first;
	auto& nomGraph = graphMap.at("PA").at(ch.first).at(v.first).at("Nominal");
	auto& theGraph = graphMap.at("PA").at(ch.first).at(v.first).at("Theory");
	//
	// Create Canvas
	TCanvas c("c", "c", 1000, 1000); c.cd();
	//
	// Create the Text Info
	TLatex tex; tex.SetNDC(); tex.SetTextSize(0.035); float dy = 0;
	std::vector< std::string > textToPrint;
	std::string sampleLabel = "W#kern[0.2]{#rightarrow}#kern[0.2]{#mu}#kern[0.2]{#nu_{#mu}}";
	if (chg == "Pl") { sampleLabel = "W^{+}#kern[0.2]{#rightarrow}#kern[0.2]{#mu^{+}}#kern[0.2]{#nu_{#mu}}"; }
	if (chg == "Mi") { sampleLabel = "W^{#font[122]{\55}}#kern[0.2]{#rightarrow}#kern[0.2]{#mu^{#font[122]{\55}}}#kern[0.2]{#bar{#nu}_{#mu}}"; }
	textToPrint.push_back(sampleLabel);
	if (accType=="") { textToPrint.push_back("p^{#mu}_{T} > 25 GeV/c"); }
	//
	// Declare the graph vector (for drawing with markers)
	std::vector< TGraphAsymmErrors > grVec, grRatio;
	// Initialize the Legend
	double legOff = 0.0; if (accType=="") { legOff = 0.05; }
        double legXmin = 0.2 , legXmax = 0.4 , legYmin = (0.58 - legOff) , legYmax = (0.81 - legOff);
        if (var=="Cross_Section" && chg=="Pl") { legXmin = 0.35; legXmax = 0.55; legYmin = 0.22; legYmax = 0.45; }
        if (var=="ForwardBackward_Ratio") { legXmin = 0.2; legXmax = 0.4; legYmin = 0.15; legYmax = 0.38; }
	TLegend leg(legXmin, legYmin, legXmax, legYmax);
	// Initialize the graph x range variables
	double xMin=0.0 , xMax=0.0 , yErrMin=9999999. , yErrMax=-1.;
	// Format the graphs
	const bool incAcc = (accType!="");
	grVec.push_back(nomGraph.at("Err_Stat"));
	grVec.push_back(nomGraph.at("Err_Tot"));
	grVec.push_back(nomGraph.at("Err_Tot"));
        if (type=="EPPS16") {
          grVec.push_back(theGraph.at("CT14"));
          grVec.push_back(theGraph.at("EPPS16"));
          grVec.push_back(theGraph.at("nCTEQ15"));
        }
        else if (type=="EPS09") {
          grVec.push_back(theGraph.at("CT10"));
          grVec.push_back(theGraph.at("EPS09"));
        }
        else if (type=="EPS09_CT14") {
          grVec.push_back(theGraph.at("CT14"));
          grVec.push_back(theGraph.at("EPS09_CT14"));
        }
	for (auto& gr : grVec) { formatResultsGraph(gr, col, var, chg, useEtaCM, incAcc); }
        //
        for (int i=0; i<grVec[0].GetN(); i++) { grVec[0].SetPointEXhigh(i, 0.0); grVec[0].SetPointEXlow(i, 0.0); }
        for (int i=0; i<grVec[0].GetN(); i++) { double x, y; grVec[0].GetPoint(i, x, y); grVec[1].SetPoint(i, x, y+grVec[1].GetErrorYhigh(i)); grVec[2].SetPoint(i, x, y-grVec[2].GetErrorYlow(i)); }
	grVec[1].SetMarkerSize(0);
        for (int i=0; i<grVec[1].GetN(); i++) { grVec[1].SetPointEYhigh(i, 0.0); grVec[1].SetPointEYlow(i, 0.0); }
        for (int i=0; i<grVec[1].GetN(); i++) { grVec[1].SetPointEXhigh(i, 0.5*grVec[1].GetErrorXhigh(i)); grVec[1].SetPointEXlow(i, 0.5*grVec[1].GetErrorXlow(i)); }
	grVec[2].SetMarkerSize(0);
        for (int i=0; i<grVec[2].GetN(); i++) { grVec[2].SetPointEYhigh(i, 0.0); grVec[2].SetPointEYlow(i, 0.0); }
        for (int i=0; i<grVec[2].GetN(); i++) { grVec[2].SetPointEXhigh(i, 0.5*grVec[2].GetErrorXhigh(i)); grVec[2].SetPointEXlow(i, 0.5*grVec[2].GetErrorXlow(i)); }
        //
	grVec[3].SetFillColor(kYellow);
	grVec[3].SetLineColor(kRed);
	grVec[3].SetFillStyle(1001);
	grVec[3].SetLineWidth(4);
	grVec[3].SetMarkerSize(0);
        //
	grVec[4].SetFillColor(kGreen+2);
	grVec[4].SetLineColor(kGreen+2);
	grVec[4].SetFillStyle(3275);
	grVec[4].SetLineStyle(7);
	grVec[4].SetLineWidth(5);
	grVec[4].SetMarkerSize(0);
        //
        if (type=="EPPS16") {
          grVec[5].SetFillColor(kOrange+2);
          grVec[5].SetLineColor(kOrange+2);
          grVec[5].SetFillStyle(3257);
          grVec[5].SetLineStyle(7);
          grVec[5].SetLineWidth(5);
          grVec[5].SetMarkerSize(0);
        }
        //
        gStyle->SetHatchesSpacing(1.9);
        gStyle->SetHatchesLineWidth(2);
        //
        auto h1 = graphToHist(grVec[3]); for (int i=1; i<=grVec[3].GetN(); i++) { h1.SetBinError(i, 0.0002); }
        auto h2 = graphToHist(grVec[4]); for (int i=1; i<=grVec[4].GetN(); i++) { h2.SetBinError(i, 0.0002); }
        TH1D h3; if (type=="EPPS16") { h3 = graphToHist(grVec[5]); for (int i=1; i<=grVec[5].GetN(); i++) { h3.SetBinError(i, 0.0002); } }
        //
	// Create Legend
	formatLegendEntry(*leg.AddEntry(&grVec[0], "Data", "pe"));
        if (type=="EPPS16") {
          formatLegendEntry(*leg.AddEntry(&grVec[3], "CT14 (68% CL)", "lf"));
          formatLegendEntry(*leg.AddEntry(&grVec[4], "CT14+EPPS16 (68% CL)", "lf"));
          formatLegendEntry(*leg.AddEntry(&grVec[5], "CT14+nCTEQ15 (68% CL)", "lf"));
        }
        else if (type=="EPS09") {
          formatLegendEntry(*leg.AddEntry(&grVec[3], "CT10 (68% CL)", "lf"));
          formatLegendEntry(*leg.AddEntry(&grVec[4], "CT10+EPS09 (68% CL)", "lf"));
        }
        else if (type=="EPS09_CT14") {
          formatLegendEntry(*leg.AddEntry(&grVec[3], "CT14 (68% CL)", "lf"));
          formatLegendEntry(*leg.AddEntry(&grVec[4], "CT14+EPS09 (68% CL)", "lf"));
        }
        //leg.SetHeader("NLO MCFM + NLO PDF, 68% CL");
        //TLegendEntry *header = (TLegendEntry*)leg.GetListOfPrimitives()->First();
        //header->SetTextSize(0.04);
        //
        // Draw the Graphs
        //
	grVec[0].Draw("ap");
	grVec[3].Draw("same2"); h1.Draw("sameL");
	grVec[4].Draw("same2"); h2.Draw("sameL");
        if (type=="EPPS16") { grVec[5].Draw("same2"); h3.Draw("sameL"); }
	grVec[1].Draw("samep"); grVec[2].Draw("samep");
	grVec[0].Draw("samep");
	//
	xMin = grVec[0].GetXaxis()->GetXmin(); xMax = grVec[0].GetXaxis()->GetXmax();
	// Draw the Line
	TLine line_FB(xMin, 1.0, xMax, 1.0); line_FB.SetLineStyle(2);
	if (var=="ForwardBackward_Ratio") { line_FB.Draw("same"); }
	TLine line_CA(xMin, 0.0, xMax, 0.0); line_CA.SetLineStyle(2);
	if (var=="Charge_Asymmetry") { line_CA.Draw("same"); }
	// Draw the Legend
	leg.Draw("same");
	// Update
	c.Modified(); c.Update();
	// Draw the text
	tex.SetTextSize(0.055); tex.DrawLatex(0.22, 0.84, textToPrint[0].c_str());
        tex.SetTextSize(0.058); tex.SetTextFont(61); tex.DrawLatex(0.78, 0.84, "CMS"); tex.SetTextFont(62);
        if (!forPaper) { tex.SetTextSize(0.044); tex.SetTextFont(52); tex.DrawLatex(0.69, 0.79, "Preliminary"); tex.SetTextFont(62); }
	if (textToPrint.size()>1) { tex.SetTextSize(0.035); tex.DrawLatex(0.22, 0.78, textToPrint[1].c_str()); }
        if (var=="Cross_Section") { tex.SetTextSize(0.030); tex.DrawLatex(0.30, 0.17, Form("Lumi. uncertainty (%.1f%%) not shown", LUMIUNC_*100.)); }
	// Update
	c.Modified(); c.Update(); // Pure paranoia
	//
	// set the CMS style
	// Draw the text
	int option = 118;
	if (col.find("pPb")!=std::string::npos) option = 115;
	if (col.find("Pbp")!=std::string::npos) option = 116;
        CMS_lumi(&c, option, 33, "", false, 0.6, false);
	// Update
	c.Modified(); c.Update(); // Pure paranoia
	//
	std::string label = "";
	if (accType==""   && effType==""   ) { label = "RAW";          }
	if (accType=="MC" && effType==""   ) { label = "AccMC";        }
	if (accType==""   && effType=="MC" ) { label = "EffMC";        }
	if (accType=="MC" && effType=="MC" ) { label = "AccMC_EffMC";  }
	if (accType==""   && effType=="TnP") { label = "EffTnP";       }
	if (accType=="MC" && effType=="TnP") { label = "AccMC_EffTnP"; }
	//
	// Create Output Directory
	const std::string plotDir = outDir+"/Plots/Theory/" + (forPaper?"ForPaper/":"ForPAS/") + col+"/" + label+"/" + var;
	makeDir(plotDir + "/png/");
	makeDir(plotDir + "/pdf/");
	makeDir(plotDir + "/root/");
	//
	// Save Canvas
	const std::string name = Form("gr_WToMu%s_%s_%s_%s_%s_%s", chg.c_str(), col.c_str(), var.c_str(), label.c_str(), "NominalWithTheory", type.c_str());
	c.SaveAs(( plotDir + "/png/"  + name + ".png"  ).c_str());
	c.SaveAs(( plotDir + "/pdf/"  + name + ".pdf"  ).c_str());
	c.SaveAs(( plotDir + "/root/" + name + ".root" ).c_str());
	//
	// Clean up memory
	c.Clear(); c.Close();
      }
    }
  }
};


void drawGraphWithTheorySEP( GraphPentaMap& graphMap , const std::string& outDir , const std::string type="EPPS16", const bool useEtaCM = true , const std::string accType = "MC" , const std::string effType = "TnP" , const bool forPaper = true )
{
  //
  // Set Style
  setStyle();
  //
  std::cout << "[INFO] Drawing the output graphs with Theory predictions" << std::endl;
  //
  // Add the Theory Predictions
  std::map< std::string , std::map< std::string , std::map< std::string , std::map< std::string , std::vector< TGraphAsymmErrors > > > > > graphModelMap;
  //
  extractTheory(graphModelMap);
  //
  // Draw all graphs
  for (const auto& c : graphMap) {
    for (const auto& ch : c.second) {
      for (const auto& v : ch.second) {
	//
	const std::string col  = c.first;
	const std::string chg  = ( (ch.first!="") ? ch.first : "Inc" );
	const std::string var  = v.first;
	auto& nomGraph = graphMap.at("PA").at(ch.first).at(v.first).at("Nominal");
	auto& theGraph = graphModelMap.at(ch.first).at(v.first).at("Theory");
	//
	// Create Canvas
	TCanvas c("c", "c", 1000, 1000); c.cd();
	//
	// Create the Text Info
	TLatex tex; tex.SetNDC(); tex.SetTextSize(0.035); float dy = 0;
	std::vector< std::string > textToPrint;
	std::string sampleLabel = "W#kern[0.2]{#rightarrow}#kern[0.2]{#mu}#kern[0.2]{#nu_{#mu}}";
	if (chg == "Pl") { sampleLabel = "W^{+}#kern[0.2]{#rightarrow}#kern[0.2]{#mu^{+}}#kern[0.2]{#nu_{#mu}}"; }
	if (chg == "Mi") { sampleLabel = "W^{#font[122]{\55}}#kern[0.2]{#rightarrow}#kern[0.2]{#mu^{#font[122]{\55}}}#kern[0.2]{#bar{#nu}_{#mu}}"; }
	textToPrint.push_back(sampleLabel);
	if (accType=="") { textToPrint.push_back("p^{#mu}_{T} > 25 GeV/c"); }
	//
	// Declare the graph vector (for drawing with markers)
	std::vector< TGraphAsymmErrors > grVec, grRatio;
	// Initialize the Legend
	double legOff = 0.0; if (accType=="") { legOff = 0.05; }
        double legXmin = 0.2 , legXmax = 0.4 , legYmin = (0.68 - legOff) , legYmax = (0.81 - legOff);
        if (var=="Cross_Section" && chg=="Pl") { legXmin = 0.35; legXmax = 0.55; legYmin = 0.22; legYmax = 0.35; }
        if (var=="ForwardBackward_Ratio") { legXmin = 0.2; legXmax = 0.4; legYmin = 0.15; legYmax = 0.28; }
	TLegend leg(legXmin, legYmin, legXmax, legYmax);
	// Initialize the graph x range variables
	double xMin=0.0 , xMax=0.0 , yErrMin=9999999. , yErrMax=-1.;
	// Format the graphs
	const bool incAcc = (accType!="");
	grVec.push_back(nomGraph.at("Err_Stat"));
	grVec.push_back(nomGraph.at("Err_Tot"));
	grVec.push_back(nomGraph.at("Err_Tot"));
        if (type=="EPPS16") {
          for (const auto& gr : theGraph.at("CT14")   ) { grVec.push_back(gr); }
          for (const auto& gr : theGraph.at("EPPS16") ) { grVec.push_back(gr); }
          for (const auto& gr : theGraph.at("nCTEQ15")) { grVec.push_back(gr); }
        }
	for (auto& gr : grVec) { formatResultsGraph(gr, col, var, chg, useEtaCM, incAcc); }
        //
        for (int i=0; i<grVec[0].GetN(); i++) { grVec[0].SetPointEXhigh(i, 0.0); grVec[0].SetPointEXlow(i, 0.0); }
        for (int i=0; i<grVec[0].GetN(); i++) { double x, y; grVec[0].GetPoint(i, x, y); grVec[1].SetPoint(i, x, y+grVec[1].GetErrorYhigh(i)); grVec[2].SetPoint(i, x, y-grVec[2].GetErrorYlow(i)); }
	grVec[1].SetMarkerSize(0);
        for (int i=0; i<grVec[1].GetN(); i++) { grVec[1].SetPointEYhigh(i, 0.0); grVec[1].SetPointEYlow(i, 0.0); }
        for (int i=0; i<grVec[1].GetN(); i++) { grVec[1].SetPointEXhigh(i, 0.5*grVec[1].GetErrorXhigh(i)); grVec[1].SetPointEXlow(i, 0.5*grVec[1].GetErrorXlow(i)); }
	grVec[2].SetMarkerSize(0);
        for (int i=0; i<grVec[2].GetN(); i++) { grVec[2].SetPointEYhigh(i, 0.0); grVec[2].SetPointEYlow(i, 0.0); }
        for (int i=0; i<grVec[2].GetN(); i++) { grVec[2].SetPointEXhigh(i, 0.5*grVec[2].GetErrorXhigh(i)); grVec[2].SetPointEXlow(i, 0.5*grVec[2].GetErrorXlow(i)); }
        //
        int nIni = 0 , nEnd = (theGraph.at("CT14").size());
        for (int i=nIni; i<nEnd; i++) {
          grVec[3+i].SetLineColor(kRed);
          grVec[3+i].SetLineWidth(1);
          grVec[3+i].SetMarkerSize(0);
        }
        //
        nIni = theGraph.at("CT14").size() , nEnd = (theGraph.at("EPPS16").size()+theGraph.at("CT14").size());
        for (int i=nIni; i<nEnd; i++) {
          grVec[3+i].SetLineColor(kGreen+2);
          grVec[3+i].SetLineWidth(1);
          grVec[3+i].SetMarkerSize(0);
        }
        //
        nIni = (theGraph.at("EPPS16").size()+theGraph.at("CT14").size()) , nEnd = (theGraph.at("nCTEQ15").size()+theGraph.at("EPPS16").size()+theGraph.at("CT14").size());
        for (int i=nIni; i<nEnd; i++) {
          grVec[i+3].SetLineColor(kOrange+2);
          grVec[i+3].SetLineWidth(4);
          grVec[i+3].SetMarkerSize(0);
        }
        //
        std::vector< TH1D > hT;
        for (uint j=0; j<(grVec.size()-3); j++) { hT.push_back(graphToHist(grVec[j+3])); for (int i=1; i<=grVec[j+3].GetN(); i++) { hT[j].SetBinError(i, 0.0002); } }
        //
	// Create Legend
	formatLegendEntry(*leg.AddEntry(&grVec[0], "Data", "pe"));
        if (type=="EPPS16") {
          //formatLegendEntry(*leg.AddEntry(&grVec[3], "CT14", "l"));
          //formatLegendEntry(*leg.AddEntry(&grVec[3+theGraph.at("CT14").size()], "CT14+EPPS16", "l"));
          formatLegendEntry(*leg.AddEntry(&grVec[3+theGraph.at("CT14").size()+theGraph.at("EPPS16").size()], "CT14+nCTEQ15 error sets", "l"));
        }
        //
        // Draw the Graphs
        //
	grVec[0].Draw("ap");
        nIni = (theGraph.at("EPPS16").size()+theGraph.at("CT14").size());
        for (uint i=nIni; i<(grVec.size()-3); i++) { hT[i].Draw("sameL"); }
	grVec[1].Draw("samep"); grVec[2].Draw("samep");
	grVec[0].Draw("samep");
	//
	xMin = grVec[0].GetXaxis()->GetXmin(); xMax = grVec[0].GetXaxis()->GetXmax();
	// Draw the Line
	TLine line_FB(xMin, 1.0, xMax, 1.0); line_FB.SetLineStyle(2);
	if (var=="ForwardBackward_Ratio") { line_FB.Draw("same"); }
	TLine line_CA(xMin, 0.0, xMax, 0.0); line_CA.SetLineStyle(2);
	if (var=="Charge_Asymmetry") { line_CA.Draw("same"); }
	// Draw the Legend
	leg.Draw("same");
	// Update
	c.Modified(); c.Update();
	// Draw the text
	tex.SetTextSize(0.055); tex.DrawLatex(0.22, 0.84, textToPrint[0].c_str());
        tex.SetTextSize(0.058); tex.SetTextFont(61); tex.DrawLatex(0.78, 0.84, "CMS"); tex.SetTextFont(62);
        if (!forPaper) { tex.SetTextSize(0.044); tex.SetTextFont(52); tex.DrawLatex(0.69, 0.79, "Preliminary"); tex.SetTextFont(62); }
	if (textToPrint.size()>1) { tex.SetTextSize(0.035); tex.DrawLatex(0.22, 0.78, textToPrint[1].c_str()); }
        if (var=="Cross_Section") { tex.SetTextSize(0.030); tex.DrawLatex(0.30, 0.17, Form("Lumi. uncertainty (%.1f%%) not shown", LUMIUNC_*100.)); }
	// Update
	c.Modified(); c.Update(); // Pure paranoia
	//
	// set the CMS style
	// Draw the text
	int option = 118;
	if (col.find("pPb")!=std::string::npos) option = 115;
	if (col.find("Pbp")!=std::string::npos) option = 116;
        CMS_lumi(&c, option, 33, "", false, 0.6, false);
	// Update
	c.Modified(); c.Update(); // Pure paranoia
	//
	std::string label = "";
	if (accType==""   && effType==""   ) { label = "RAW";          }
	if (accType=="MC" && effType==""   ) { label = "AccMC";        }
	if (accType==""   && effType=="MC" ) { label = "EffMC";        }
	if (accType=="MC" && effType=="MC" ) { label = "AccMC_EffMC";  }
	if (accType==""   && effType=="TnP") { label = "EffTnP";       }
	if (accType=="MC" && effType=="TnP") { label = "AccMC_EffTnP"; }
	//
	// Create Output Directory
	const std::string plotDir = outDir+"/Plots/Theory/" + (forPaper?"ForPaper/":"ForPAS/") + col+"/" + label+"/" + var;
	makeDir(plotDir + "/png/");
	makeDir(plotDir + "/pdf/");
	makeDir(plotDir + "/root/");
	//
	// Save Canvas
	const std::string name = Form("gr_WToMu%s_%s_%s_%s_%s_%s", chg.c_str(), col.c_str(), var.c_str(), label.c_str(), "NominalWithTheorySEP", type.c_str());
	c.SaveAs(( plotDir + "/png/"  + name + ".png"  ).c_str());
	c.SaveAs(( plotDir + "/pdf/"  + name + ".pdf"  ).c_str());
	c.SaveAs(( plotDir + "/root/" + name + ".root" ).c_str());
	//
	// Clean up memory
	c.Clear(); c.Close();
      }
    }
  }
};


void drawGraphWithTheoryREW( GraphPentaMap& graphMap , const std::string& outDir , const bool useEtaCM = true , const std::string accType = "MC" , const std::string effType = "TnP" , const bool forPaper = true )
{
  //
  // Set Style
  setStyle();
  //
  std::cout << "[INFO] Drawing the output graphs with Theory predictions" << std::endl;
  //
  // Add the Theory Predictions
  A1m_REW(graphMap["PA"]["Mi"]["ForwardBackward_Ratio"]["Theory"]);
  A1p_REW(graphMap["PA"]["Pl"]["ForwardBackward_Ratio"]["Theory"]);
  A3_REW (graphMap["PA"][  ""]["ForwardBackward_Ratio"]["Theory"]);
  Wp_REW (graphMap["PA"]["Pl"]["Cross_Section"]["Theory"]);
  Wm_REW (graphMap["PA"]["Mi"]["Cross_Section"]["Theory"]);
  chasym_REW(graphMap["PA"][""]["Charge_Asymmetry"]["Theory"]);
  //
  // Draw all graphs
  for (const auto& c : graphMap) {
    for (const auto& ch : c.second) {
      for (const auto& v : ch.second) {
	//
	const std::string col  = c.first;
	const std::string chg  = ( (ch.first!="") ? ch.first : "Inc" );
	const std::string var  = v.first;
	auto& nomGraph = graphMap.at("PA").at(ch.first).at(v.first).at("Nominal");
	auto& theGraph = graphMap.at("PA").at(ch.first).at(v.first).at("Theory");
	//
	// Create Canvas
	TCanvas c("c", "c", 1000, 1000); c.cd();
	//
	// Create the Text Info
	TLatex tex; tex.SetNDC(); tex.SetTextSize(0.035); float dy = 0;
	std::vector< std::string > textToPrint;
	std::string sampleLabel = "W#kern[0.2]{#rightarrow}#kern[0.2]{#mu}#kern[0.2]{#nu_{#mu}}";
	if (chg == "Pl") { sampleLabel = "W^{+}#kern[0.2]{#rightarrow}#kern[0.2]{#mu^{+}}#kern[0.2]{#nu_{#mu}}"; }
	if (chg == "Mi") { sampleLabel = "W^{#font[122]{\55}}#kern[0.2]{#rightarrow}#kern[0.2]{#mu^{#font[122]{\55}}}#kern[0.2]{#bar{#nu}_{#mu}}"; }
	textToPrint.push_back(sampleLabel);
	if (accType=="") { textToPrint.push_back("p^{#mu}_{T} > 25 GeV/c"); }
	//
	// Declare the graph vector (for drawing with markers)
	std::vector< TGraphAsymmErrors > grVec, grRatio;
	// Initialize the Legend
	double legOff = 0.0; if (accType=="") { legOff = 0.05; }
        double legXmin = 0.2 , legXmax = 0.4 , legYmin = (0.58 - legOff) , legYmax = (0.81 - legOff);
        if (var=="Cross_Section" && chg=="Pl") { legXmin = 0.35; legXmax = 0.55; legYmin = 0.22; legYmax = 0.45; }
        if (var=="ForwardBackward_Ratio") { legXmin = 0.2; legXmax = 0.4; legYmin = 0.15; legYmax = 0.38; }
	TLegend leg(legXmin, legYmin, legXmax, legYmax);
	// Initialize the graph x range variables
	double xMin=0.0 , xMax=0.0 , yErrMin=9999999. , yErrMax=-1.;
	// Format the graphs
	const bool incAcc = (accType!="");
	grVec.push_back(nomGraph.at("Err_Stat"));
	grVec.push_back(nomGraph.at("Err_Tot"));
	grVec.push_back(nomGraph.at("Err_Tot"));
        grVec.push_back(theGraph.at("CT14"));
        grVec.push_back(theGraph.at("EPPS16"));
        grVec.push_back(theGraph.at("nCTEQ15"));
	for (auto& gr : grVec) { formatResultsGraph(gr, col, var, chg, useEtaCM, incAcc); }
        //
        for (int i=0; i<grVec[0].GetN(); i++) { grVec[0].SetPointEXhigh(i, 0.0); grVec[0].SetPointEXlow(i, 0.0); }
        for (int i=0; i<grVec[0].GetN(); i++) { double x, y; grVec[0].GetPoint(i, x, y); grVec[1].SetPoint(i, x, y+grVec[1].GetErrorYhigh(i)); grVec[2].SetPoint(i, x, y-grVec[2].GetErrorYlow(i)); }
	grVec[1].SetMarkerSize(0);
        for (int i=0; i<grVec[1].GetN(); i++) { grVec[1].SetPointEYhigh(i, 0.0); grVec[1].SetPointEYlow(i, 0.0); }
        for (int i=0; i<grVec[1].GetN(); i++) { grVec[1].SetPointEXhigh(i, 0.5*grVec[1].GetErrorXhigh(i)); grVec[1].SetPointEXlow(i, 0.5*grVec[1].GetErrorXlow(i)); }
	grVec[2].SetMarkerSize(0);
        for (int i=0; i<grVec[2].GetN(); i++) { grVec[2].SetPointEYhigh(i, 0.0); grVec[2].SetPointEYlow(i, 0.0); }
        for (int i=0; i<grVec[2].GetN(); i++) { grVec[2].SetPointEXhigh(i, 0.5*grVec[2].GetErrorXhigh(i)); grVec[2].SetPointEXlow(i, 0.5*grVec[2].GetErrorXlow(i)); }
        //
	grVec[3].SetFillColor(kYellow);
	grVec[3].SetLineColor(kRed);
	grVec[3].SetFillStyle(1001);
	grVec[3].SetLineWidth(4);
	grVec[3].SetMarkerSize(0);
        //
	grVec[4].SetFillColor(kGreen+2);
	grVec[4].SetLineColor(kGreen+2);
	grVec[4].SetFillStyle(3275);
	grVec[4].SetLineStyle(7);
	grVec[4].SetLineWidth(5);
	grVec[4].SetMarkerSize(0);
        //
        grVec[5].SetFillColor(kOrange+2);
        grVec[5].SetLineColor(kOrange+2);
        grVec[5].SetFillStyle(3257);
        grVec[5].SetLineStyle(7);
        grVec[5].SetLineWidth(5);
        grVec[5].SetMarkerSize(0);
        //
        gStyle->SetHatchesSpacing(1.9);
        gStyle->SetHatchesLineWidth(2);
        //
        auto h1 = graphToHist(grVec[3]); for (int i=1; i<=grVec[3].GetN(); i++) { h1.SetBinError(i, 0.0002); }
        auto h2 = graphToHist(grVec[4]); for (int i=1; i<=grVec[4].GetN(); i++) { h2.SetBinError(i, 0.0002); }
        TH1D h3; h3 = graphToHist(grVec[5]); for (int i=1; i<=grVec[5].GetN(); i++) { h3.SetBinError(i, 0.0002); }
        //
	// Create Legend
	formatLegendEntry(*leg.AddEntry(&grVec[0], "Data", "pe"));
        formatLegendEntry(*leg.AddEntry(&grVec[3], "CT14 (68% CL)", "lf"));
        formatLegendEntry(*leg.AddEntry(&grVec[4], "CT14+EPPS16 (68% CL)", "lf"));
        formatLegendEntry(*leg.AddEntry(&grVec[5], "CT14+nCTEQ15 (68% CL)", "lf"));
        //leg.SetHeader("NLO MCFM + NLO PDF, 68% CL");
        //TLegendEntry *header = (TLegendEntry*)leg.GetListOfPrimitives()->First();
        //header->SetTextSize(0.04);
        //
        // Draw the Graphs
        //
	grVec[0].Draw("ap");
	grVec[3].Draw("same2"); h1.Draw("sameL");
	grVec[4].Draw("same2"); h2.Draw("sameL");
        grVec[5].Draw("same2"); h3.Draw("sameL");
	grVec[1].Draw("samep"); grVec[2].Draw("samep");
	grVec[0].Draw("samep");
	//
	xMin = grVec[0].GetXaxis()->GetXmin(); xMax = grVec[0].GetXaxis()->GetXmax();
	// Draw the Line
	TLine line_FB(xMin, 1.0, xMax, 1.0); line_FB.SetLineStyle(2);
	if (var=="ForwardBackward_Ratio") { line_FB.Draw("same"); }
	TLine line_CA(xMin, 0.0, xMax, 0.0); line_CA.SetLineStyle(2);
	if (var=="Charge_Asymmetry") { line_CA.Draw("same"); }
	// Draw the Legend
	leg.Draw("same");
	// Update
	c.Modified(); c.Update();
	// Draw the text
	tex.SetTextSize(0.055); tex.DrawLatex(0.22, 0.84, textToPrint[0].c_str());
        tex.SetTextSize(0.058); tex.SetTextFont(61); tex.DrawLatex(0.78, 0.84, "CMS"); tex.SetTextFont(62);
        if (!forPaper) { tex.SetTextSize(0.044); tex.SetTextFont(52); tex.DrawLatex(0.69, 0.79, "Preliminary"); tex.SetTextFont(62); }
	if (textToPrint.size()>1) { tex.SetTextSize(0.035); tex.DrawLatex(0.22, 0.78, textToPrint[1].c_str()); }
        if (var=="Cross_Section") { tex.SetTextSize(0.030); tex.DrawLatex(0.30, 0.17, Form("Lumi. uncertainty (%.1f%%) not shown", LUMIUNC_*100.)); }
	// Update
	c.Modified(); c.Update(); // Pure paranoia
	//
	// set the CMS style
	// Draw the text
	int option = 118;
	if (col.find("pPb")!=std::string::npos) option = 115;
	if (col.find("Pbp")!=std::string::npos) option = 116;
        CMS_lumi(&c, option, 33, "", false, 0.6, false);
	// Update
	c.Modified(); c.Update(); // Pure paranoia
	//
	std::string label = "";
	if (accType==""   && effType==""   ) { label = "RAW";          }
	if (accType=="MC" && effType==""   ) { label = "AccMC";        }
	if (accType==""   && effType=="MC" ) { label = "EffMC";        }
	if (accType=="MC" && effType=="MC" ) { label = "AccMC_EffMC";  }
	if (accType==""   && effType=="TnP") { label = "EffTnP";       }
	if (accType=="MC" && effType=="TnP") { label = "AccMC_EffTnP"; }
	//
	// Create Output Directory
	const std::string plotDir = outDir+"/Plots/TheoryREW/" + (forPaper?"ForPaper/":"ForPAS/") + col+"/" + label+"/" + var;
	makeDir(plotDir + "/png/");
	makeDir(plotDir + "/pdf/");
	makeDir(plotDir + "/root/");
	//
	// Save Canvas
	const std::string name = Form("gr_WToMu%s_%s_%s_%s_%s_%s", chg.c_str(), col.c_str(), var.c_str(), label.c_str(), "NominalWithTheory", "REW");
	c.SaveAs(( plotDir + "/png/"  + name + ".png"  ).c_str());
	c.SaveAs(( plotDir + "/pdf/"  + name + ".pdf"  ).c_str());
	c.SaveAs(( plotDir + "/root/" + name + ".root" ).c_str());
	//
	// Clean up memory
	c.Clear(); c.Close();
      }
    }
  }
};


void graphEta5TeVto8TeV( TGraphAsymmErrors& gr )
{
  const double shift = TMath::Log(8.16/5.02);
  for (int i = 0; i < gr.GetN(); i++){
    double eta, y;
    gr.GetPoint(i, eta, y);
    double eta8TeV = 0.0;
    if (eta>=0.0) { eta += shift; }
    else          { eta -= shift; }
    gr.SetPoint(i, eta, y);
  }
};


void graphShiftChaAsym( TGraphAsymmErrors& gr , const double E=5020. )
{
  const double shift = TMath::Log(8160./E);
  for (int i = 0; i < gr.GetN(); i++){
    double eta, y;
    gr.GetPoint(i, eta, y);
    double eta8TeV = 0.0;
    if (eta>=0.0) { eta += shift; }
    else          { eta -= shift; }
    gr.SetPoint(i, eta, y);
  }
};


void graphShiftXSec( TGraphAsymmErrors& gr , const double E=8160. , const int addOnly=0 )
{
  auto tmp = gr;
  // Remove all unwanted points
  if (addOnly!=0) {
    for (int i = 0; i < tmp.GetN(); ){
      double eta, xSec;
      tmp.GetPoint(i, eta, xSec);
      if      (addOnly<0 && eta>=0) { tmp.RemovePoint(i); }
      else if (addOnly>0 && eta<=0) { tmp.RemovePoint(i); }
      else { i++; }
    }
  }
  // Make the change of variables
  gr = tmp;
  for (int i = 0; i < tmp.GetN(); i++){
    double eta, xSec;
    tmp.GetPoint(i, eta, xSec);
    //
    double etaErrLo = tmp.GetErrorXlow(i);
    double etaErrHi = tmp.GetErrorXhigh(i);;
    double x = (80.4/E)*std::exp(eta);
    double xErrLo = (x - (80.4/E)*std::exp(eta-etaErrLo));
    double xErrHi = ((80.4/E)*std::exp(eta+etaErrHi) - x);
    if (addOnly<0) {
      x = (80.4/E)*std::exp(-eta);
      xErrLo = (x - (80.4/E)*std::exp(-(eta+etaErrHi)));
      xErrHi = ((80.4/E)*std::exp(-(eta-etaErrLo)) - x);
    }
    //
    double xSecErrLo = tmp.GetErrorYlow(i);
    double xSecErrHi = tmp.GetErrorYhigh(i);
    double y = std::pow(E, -0.8)*1000.*xSec;
    double yErrLo = std::pow(E, -0.8)*1000.*xSecErrLo;
    double yErrHi = std::pow(E, -0.8)*1000.*xSecErrHi;
    //
    int iGr = i;
    if (addOnly<0) { iGr = (tmp.GetN()-i-1); }
    gr.SetPoint(iGr, x, y);
    gr.SetPointEXlow(iGr, xErrLo);
    gr.SetPointEXhigh(iGr, xErrHi);
    gr.SetPointEYlow(iGr, yErrLo);
    gr.SetPointEYhigh(iGr, yErrHi);
  }
};


void drawGraphWithpPb5TeV( GraphPentaMap& graphMap , const std::string& outDir , const bool useEtaCM = true , const std::string accType = "MC" , const std::string effType = "TnP" , bool addModel=false , int shiftEta=0 , const std::string model="EPS09", int plotStyle=1)
{
  //
  // Set Style
  setStyle();
  //
  std::cout << "[INFO] Drawing the output graphs with 5 TeV pPb results" << std::endl;
  //
  // Add the 5 TeV pPb results
  HIN_13007_Wp(graphMap["PA"]["Pl"]["Cross_Section"]["CMS"]);
  HIN_13007_Wm(graphMap["PA"]["Mi"]["Cross_Section"]["CMS"]);
  HIN_13007_chasym(graphMap["PA"][""]["Charge_Asymmetry"]["CMS"]);
  if (addModel) {
    Wp(graphMap["PA"]["Pl"]["Cross_Section"]["Theory_8TeV"]);
    Wm(graphMap["PA"]["Mi"]["Cross_Section"]["Theory_8TeV"]);
    chasym(graphMap["PA"][""]["Charge_Asymmetry"]["Theory_8TeV"]);
    HIN_13007_Theory_Wp(graphMap["PA"]["Pl"]["Cross_Section"]["Theory_5TeV"]);
    HIN_13007_Theory_Wm(graphMap["PA"]["Mi"]["Cross_Section"]["Theory_5TeV"]);
    HIN_13007_Theory_chasym(graphMap["PA"][""]["Charge_Asymmetry"]["Theory_5TeV"]);
  }
  //
  // Draw all graphs
  for (const auto& c : graphMap) {
    for (const auto& ch : c.second) {
      for (const auto& v : ch.second) {
        //
        if (v.first=="ForwardBackward_Ratio") continue;
        if (v.first=="Charge_Asymmetry" && shiftEta<0) continue;
	//
	const std::string col  = c.first;
	const std::string chg  = ( (ch.first!="") ? ch.first : "Inc" );
	const std::string var  = v.first;
	auto& nomGraph = graphMap.at("PA").at(ch.first).at(v.first).at("Nominal");
	auto& theGraph = graphMap.at("PA").at(ch.first).at(v.first).at("CMS");
	//
	// Create Canvas
	TCanvas c("c", "c", 1000, 1000); c.cd();
	//
	// Create the Text Info
	TLatex tex; tex.SetNDC(); tex.SetTextSize(0.035); float dy = 0;
	std::vector< std::string > textToPrint;
	std::string sampleLabel = "W#kern[0.2]{#rightarrow}#kern[0.2]{#mu}#kern[0.2]{#nu_{#mu}}";
	if (chg == "Pl") { sampleLabel = "W^{+}#kern[0.2]{#rightarrow}#kern[0.2]{#mu^{+}}#kern[0.2]{#nu_{#mu}}"; }
	if (chg == "Mi") { sampleLabel = "W^{#font[122]{\55}}#kern[0.2]{#rightarrow}#kern[0.2]{#mu^{#font[122]{\55}}}#kern[0.2]{#bar{#nu}_{#mu}}"; }
	textToPrint.push_back(sampleLabel);
	if (accType=="") { textToPrint.push_back("p^{#mu}_{T} > 25 GeV/c"); }
	//
	// Declare the graph vector (for drawing with markers)
	std::vector< TGraphAsymmErrors > grVec;
	// Initialize the Legend
	double legOff = 0.0; if (accType=="") { legOff = 0.05; }
        double legXmin = 0.2 , legXmax = 0.4 , legYmin = (0.58 - legOff) , legYmax = (0.81 - legOff);
        if (var=="Cross_Section") { legXmin = 0.45; legXmax = 0.65; legYmin = 0.15; legYmax = 0.38; }
        if (var=="Cross_Section" && shiftEta!=0 && chg=="Pl") { legXmin = 0.20; legXmax = 0.40; legYmin = 0.15; legYmax = 0.38; }
        if (var=="Cross_Section" && shiftEta!=0 && chg=="Mi") { legXmin = 0.20; legXmax = 0.40; legYmin = 0.45; legYmax = 0.68; }
        if (var=="Charge_Asymmetry") { legXmin = 0.45; legXmax = 0.65; legYmin = 0.17; legYmax = 0.39; }
	TLegend leg(legXmin, legYmin, legXmax, legYmax);
	// Initialize the graph x range variables
	double xMin=0.0 , xMax=0.0 , yErrMin=9999999. , yErrMax=-1.;
	// Format the graphs
	const bool incAcc = (accType!="");
	grVec.push_back(nomGraph.at("Err_Stat"));
	grVec.push_back(nomGraph.at("Err_Tot"));
	grVec.push_back(nomGraph.at("Err_Tot"));
	grVec.push_back(theGraph.at("HIN-13007_Stat"));
	grVec.push_back(theGraph.at("HIN-13007_Tot"));
	grVec.push_back(theGraph.at("HIN-13007_Tot"));
        if (var=="Cross_Section"    && shiftEta!=0) { for (uint i=0; i<3; i++) { graphShiftXSec(grVec[i], 8160., shiftEta); }; for (uint i=3; i<6; i++) { graphShiftXSec(grVec[i], 5020., shiftEta); } }
        if (var=="Charge_Asymmetry" && shiftEta >0) { for (uint i=3; i<6; i++) { graphShiftChaAsym(grVec[i], 5020.); } }
        if (addModel) {
          grVec.push_back(graphMap.at("PA").at(ch.first).at(v.first).at("Theory_8TeV").at(model));
          grVec.push_back(graphMap.at("PA").at(ch.first).at(v.first).at("Theory_5TeV").at(model));
          if (var=="Cross_Section" && shiftEta!=0) { graphShiftXSec(grVec[6], 8160., shiftEta); graphShiftXSec(grVec[7], 5020., shiftEta); }
          if (var=="Charge_Asymmetry" && shiftEta >0) { graphShiftChaAsym(grVec[7], 5020.); }
        }
	for (auto& gr : grVec) { formatResultsGraph(gr, col, var, chg, useEtaCM, incAcc, false, shiftEta); }
        if ( var == "Charge_Asymmetry" ) { grVec[0].GetYaxis()->SetRangeUser(-0.3, 0.45);   }
        if ( var == "Cross_Section" ) {
          if (shiftEta==0) { grVec[0].GetYaxis()->SetRangeUser(-10.0, 210.0); }
          else             { grVec[0].GetYaxis()->SetRangeUser( 40.0, 180.0); }
        }
        for (int i=0; i<grVec[0].GetN(); i++) { grVec[0].SetPointEXhigh(i, 0.0); grVec[0].SetPointEXlow(i, 0.0); }
        for (int i=0; i<grVec[0].GetN(); i++) { double x, y; grVec[0].GetPoint(i, x, y); grVec[1].SetPoint(i, x, y+grVec[1].GetErrorYhigh(i)); grVec[2].SetPoint(i, x, y-grVec[2].GetErrorYlow(i)); }
	grVec[1].SetMarkerSize(0);
        for (int i=0; i<grVec[1].GetN(); i++) { grVec[1].SetPointEYhigh(i, 0.0); grVec[1].SetPointEYlow(i, 0.0); }
        for (int i=0; i<grVec[1].GetN(); i++) { grVec[1].SetPointEXhigh(i, 0.5*grVec[1].GetErrorXhigh(i)); grVec[1].SetPointEXlow(i, 0.5*grVec[1].GetErrorXlow(i)); }
	grVec[2].SetMarkerSize(0);
        for (int i=0; i<grVec[2].GetN(); i++) { grVec[2].SetPointEYhigh(i, 0.0); grVec[2].SetPointEYlow(i, 0.0); }
        for (int i=0; i<grVec[2].GetN(); i++) { grVec[2].SetPointEXhigh(i, 0.5*grVec[2].GetErrorXhigh(i)); grVec[2].SetPointEXlow(i, 0.5*grVec[2].GetErrorXlow(i)); }
	grVec[3].SetMarkerStyle(21);
	grVec[3].SetMarkerColor(kBlue);
	grVec[3].SetLineColor(kBlue);
        for (int i=0; i<grVec[3].GetN(); i++) { grVec[3].SetPointEXhigh(i, 0.0); grVec[3].SetPointEXlow(i, 0.0); }
        for (int i=0; i<grVec[3].GetN(); i++) { double x, y; grVec[3].GetPoint(i, x, y); grVec[4].SetPoint(i, x, y+grVec[4].GetErrorYhigh(i)); grVec[5].SetPoint(i, x, y-grVec[5].GetErrorYlow(i)); }
	grVec[4].SetLineColor(kBlue);
	grVec[4].SetMarkerSize(0);
        for (int i=0; i<grVec[4].GetN(); i++) { grVec[4].SetPointEYhigh(i, 0.0); grVec[4].SetPointEYlow(i, 0.0); }
        for (int i=0; i<grVec[4].GetN(); i++) { grVec[4].SetPointEXhigh(i, 0.5*grVec[4].GetErrorXhigh(i)); grVec[4].SetPointEXlow(i, 0.5*grVec[4].GetErrorXlow(i)); }
	grVec[5].SetLineColor(kBlue);
	grVec[5].SetMarkerSize(0);
        for (int i=0; i<grVec[5].GetN(); i++) { grVec[5].SetPointEYhigh(i, 0.0); grVec[5].SetPointEYlow(i, 0.0); }
        for (int i=0; i<grVec[5].GetN(); i++) { grVec[5].SetPointEXhigh(i, 0.5*grVec[5].GetErrorXhigh(i)); grVec[5].SetPointEXlow(i, 0.5*grVec[5].GetErrorXlow(i)); }
        //
        TH1D h1, h2, h3;
        if (addModel) {
          grVec[6].SetFillColor(kGreen+2);
          grVec[6].SetLineColor(kGreen+2);
          grVec[6].SetFillStyle(3275);
          grVec[6].SetLineStyle(7);
          grVec[6].SetLineWidth(5);
          grVec[6].SetMarkerSize(0);
          //
          grVec[7].SetFillColor(kOrange+2);
          grVec[7].SetLineColor(kOrange+2);
          grVec[7].SetFillStyle(3257);
          grVec[7].SetLineStyle(7);
          grVec[7].SetLineWidth(5);
          grVec[7].SetMarkerSize(0);
          //
          gStyle->SetHatchesSpacing(1.9);
          gStyle->SetHatchesLineWidth(2);
          //
          h1 = graphToHist(grVec[6]); for (int i=1; i<=grVec[6].GetN(); i++) { h1.SetBinError(i, 0.0002); }
          if (var=="Charge_Asymmetry" && shiftEta>0) {
            h2 = graphToHist(grVec[7], -1); for (int i=1; i<=h2.GetNbinsX(); i++) { h2.SetBinError(i, 0.0002); }
            h3 = graphToHist(grVec[7], +1); for (int i=1; i<=h3.GetNbinsX(); i++) { h3.SetBinError(i, 0.0002); }
          }
          else { h2 = graphToHist(grVec[7]); for (int i=1; i<=grVec[7].GetN(); i++) { h2.SetBinError(i, 0.0002); } }
        }
	formatLegendEntry(*leg.AddEntry(&grVec[0], "Data 8.16 TeV", "pe"));
	formatLegendEntry(*leg.AddEntry(&grVec[3], "Data 5.02 TeV", "pe"));
        if (addModel) {
          if (model=="EPS09" ) { formatLegendEntry(*leg.AddEntry(&grVec[6], "CT10+EPS09 8.16 TeV" , "lf")); }
          if (model=="EPPS16") { formatLegendEntry(*leg.AddEntry(&grVec[6], "CT14+EPPS16 8.16 TeV", "lf")); }
        }
        if (addModel && grVec.size()==8) {
          if (model=="EPS09" ) { formatLegendEntry(*leg.AddEntry(&grVec[7], "CT10+EPS09 5.02 TeV" , "lf")); }
          if (model=="EPPS16") { formatLegendEntry(*leg.AddEntry(&grVec[7], "CT14+EPPS16 5.02 TeV", "lf")); }
        }
        //if (addModel) {
        //  leg.SetHeader("NLO MCFM + NLO PDF, 68% CL");
        //  TLegendEntry *header = (TLegendEntry*)leg.GetListOfPrimitives()->First();
        //  header->SetTextSize(0.04);
        //}
	// Draw the graphs
	grVec[0].Draw("ap");
        if (addModel) {
          grVec[6].Draw("same2"); h1.Draw("sameL");
          grVec[7].Draw("same2"); h2.Draw("sameL");
          if (var=="Charge_Asymmetry" && shiftEta>0) { h3.Draw("sameL"); }
        }
	grVec[4].Draw("samep"); grVec[5].Draw("samep");
	grVec[3].Draw("samep");
	grVec[1].Draw("samep"); grVec[2].Draw("samep");
	grVec[0].Draw("samep");
	//
	xMin = grVec[0].GetXaxis()->GetXmin(); xMax = grVec[0].GetXaxis()->GetXmax();
	// Draw the Line
	TLine line_FB(xMin, 1.0, xMax, 1.0); line_FB.SetLineStyle(2);
	if (var=="ForwardBackward_Ratio") { line_FB.Draw("same"); }
	TLine line_CA(xMin, 0.0, xMax, 0.0); line_CA.SetLineStyle(2);
	if (var=="Charge_Asymmetry") { line_CA.Draw("same"); }
	// Draw the Legend
	leg.Draw("same");
	// Update
	c.Modified(); c.Update();
	// Draw the texts
	tex.SetTextSize(0.055); tex.DrawLatex(0.22, 0.84, textToPrint[0].c_str());
        tex.SetTextSize(0.058); tex.SetTextFont(61); tex.DrawLatex(0.78, 0.84, "CMS"); tex.SetTextFont(62);
        if (plotStyle==1) { tex.SetTextSize(0.044); tex.SetTextFont(52); tex.DrawLatex(0.69, 0.79, "Preliminary"); tex.SetTextFont(62); }
        if (plotStyle==2) { tex.SetTextSize(0.044); tex.SetTextFont(52); tex.DrawLatex(0.67, 0.79, "Supplementary"); tex.SetTextFont(62); }
	if (textToPrint.size()>1) { tex.SetTextSize(0.035); tex.DrawLatex(0.22, 0.78, textToPrint[1].c_str()); }
        if (var=="Cross_Section") {
          if (shiftEta!=0) {
            if (shiftEta>0) { tex.SetTextSize(0.035); tex.DrawLatex(0.22, 0.74, "#eta^{#mu}_{CM} > 0"); }
            if (shiftEta<0) { tex.SetTextSize(0.035); tex.DrawLatex(0.22, 0.74, "#eta^{#mu}_{CM} < 0"); }
            tex.SetTextSize(0.030); tex.DrawLatex(0.22, 0.69, "Lumi. uncertainty not shown");
          }
          else { tex.SetTextSize(0.030); tex.DrawLatex(0.22, 0.74, "Lumi. uncertainty not shown"); }
        }
	// Update
	c.Modified(); c.Update(); // Pure paranoia
	//
	// set the CMS style
	int option = 118;
	if (col.find("pPb")!=std::string::npos) option = 115;
	if (col.find("Pbp")!=std::string::npos) option = 116;
        CMS_lumi(&c, option, 33, "", false, 0.6, false);
	// Update
	c.Modified(); c.Update(); // Pure paranoia
        if (var=="Cross_Section" && shiftEta!=0) { c.SetLogx(); c.Modified(); c.Update(); }
	//
	std::string label = "";
	if (accType==""   && effType==""   ) { label = "RAW";          }
	if (accType=="MC" && effType==""   ) { label = "AccMC";        }
	if (accType==""   && effType=="MC" ) { label = "EffMC";        }
	if (accType=="MC" && effType=="MC" ) { label = "AccMC_EffMC";  }
	if (accType==""   && effType=="TnP") { label = "EffTnP";       }
	if (accType=="MC" && effType=="TnP") { label = "AccMC_EffTnP"; }
	//
	// Create Output Directory
	const std::string plotDir = outDir+"/Plots/HIN_13007/" + (plotStyle==0?"ForPaper/":"ForPAS/") +  col+"/" + label+"/" + var;
	makeDir(plotDir + "/png/");
	makeDir(plotDir + "/pdf/");
	makeDir(plotDir + "/root/");
	//
	// Save Canvas
	std::string name = Form("gr_WToMu%s_%s_%s_%s_%s%s", chg.c_str(), col.c_str(), var.c_str(), label.c_str(), "HIN_13007", (addModel?("_With"+model).c_str():""));
        if (var=="Charge_Asymmetry" && shiftEta>0) { name += "_ShiftEta"; }
        if (var=="Cross_Section") { if (shiftEta>0) { name += "_ShiftEtaPos"; } else if (shiftEta<0) { name += "_ShiftEtaNeg"; } }
	c.SaveAs(( plotDir + "/png/"  + name + ".png"  ).c_str());
	c.SaveAs(( plotDir + "/pdf/"  + name + ".pdf"  ).c_str());
	c.SaveAs(( plotDir + "/root/" + name + ".root" ).c_str());
	//
	// Clean up memory
	c.Clear(); c.Close();
      }
    }
  }
};


const char*  sgn (const double n) {
  if (n >= 0.) { return "$+$"; } else { return "$-$"; }
};


void createYieldTable(std::vector< std::string >& texTable, const std::vector< std::string >& colVar, const std::vector< std::string >& colTitle1,
                      const VarBinMap& inputVar, const std::string& col, const std::string& chg)
{
  //
  const uint nCol = colVar.size();
  //
  texTable.push_back("  \\renewcommand{\\arraystretch}{1.5}");
  texTable.push_back(Form("  \\begin{tabular}{|c|*%dc|}", (nCol-1)));
  texTable.push_back("    \\hline");
  std::string tmp;
  tmp = ("    ");
  for (uint i = 0; i < nCol; i++) {
    tmp += colTitle1[i];
    if (i<(nCol-1)) { tmp += " & "; }
    else { tmp += "\\\\"; }
  }
  texTable.push_back(tmp);
  texTable.push_back("    \\hline\\hline");
  //
  for (const auto& b : inputVar.at(col)) {
    tmp = ("    ");
    for (uint i = 0; i < nCol; i++) {
      const auto& v = colVar[i];
      const auto& var = b.second.at(chg).at(v);
      std::string val;
      if (v=="Muon_Eta") {
        const double min = b.first.etabin().low();
        const double max = b.first.etabin().high();
        val = Form("%s%.2f , %s%.2f", sgn(min), abs(min), sgn(max) , abs(max));
      }
      else if (v=="N_DS_Entries") {
        val = Form("%g", var.at("Val"));
      }
      else if (v=="N_FIT_Entries") {
        val = Form("%.0f", var.at("Val"));
      }
      else if (v=="TEST_FIT") {
        const double pVal = var.at("Val");
        val = Form("%.2f", pVal);
      }
      else if (v=="TEST_FIT_BCChi2") {
        const double pVal = var.at("Val");
        //const double pVal = var.at("Chi2")/var.at("NDoF");
        val = Form("%.2f", pVal);
      }
      else if (v=="Acceptance_MC" || v=="Efficiency_MC" || v=="Efficiency_TnP") {
        const double errLow = var.at("Err_Stat_Low")*100.;
        const double errHi  = var.at("Err_Stat_High")*100.;
        if (std::abs(errLow-errHi)<0.005) {
          if (errLow<0.05) { val = Form("$%.2f \\pm %.2f$", var.at("Val")*100., errLow ); }
          else if (errLow<0.5) { val = Form("$%.1f \\pm %.1f$", var.at("Val")*100., errLow ); }
          else { val = Form("$%.0f \\pm %.0f$", var.at("Val")*100., errLow ); }
        }
        else { val = Form("$%.2f + %.2f - %.2f$", var.at("Val")*100., errHi, errLow ); }
      }
      else {
        const double errLow = var.at("Err_Stat_Low");
        const double errHi  = var.at("Err_Stat_High");
        if (errLow==errHi) {
          if (errLow < 0.5) { val = Form("$%.1f \\pm %.1f$", var.at("Val"), errLow ); }
          else { val = Form("$%.0f \\pm %.0f$", var.at("Val"), errLow ); }
        }
        else { val = Form("$%.2f + %.2f - %.2f$", var.at("Val"), errHi, errLow ); }
      }
      tmp += val;
      if (i<(nCol-1)) { tmp += " & "; } else { tmp += "\\\\"; }
    }
    texTable.push_back(tmp);
    texTable.push_back("    \\hline");
  }
  //
  texTable.push_back("  \\end{tabular}");
  //
};


void makeRawYieldsTable(std::ofstream& file, const VarBinMap& inputVar, const std::string& col, const std::string chg)
{
  //
  const auto& inVar = inputVar.at(col).begin()->second.at(chg);
  bool useEtaCM = false; if (inVar.count("useEtaCM")>0 && inVar.at("useEtaCM").at("Val")>0) { useEtaCM = true; }
  //
  // Extract kinematic info
  std::string pTCUT = "";
  if (inVar.count("Muon_Pt")>0) {
    if (inVar.at("Muon_Pt").at("Min") >0.0 && inVar.at("Muon_Pt").at("Max")>=1000.0) { pTCUT = Form("$p_{T} > %.0f$~GeV/c", inVar.at("Muon_Pt").at("Min")); }
    if (inVar.at("Muon_Pt").at("Min")<=0.0 && inVar.at("Muon_Pt").at("Max") <1000.0) { pTCUT = Form("$p_{T} < %.0f$~GeV/c", inVar.at("Muon_Pt").at("Max")); }
    if (inVar.at("Muon_Pt").at("Min") >0.0 && inVar.at("Muon_Pt").at("Max") <1000.0) { pTCUT = Form("$%.0f < p_{T} < %.0f$~GeV/c", inVar.at("Muon_Pt").at("Min"), inVar.at("Muon_Pt").at("Max")); }
  }
  //
  // Initialize the latex table
  std::vector< std::string > texTable;
  //
  // Determine number of columns
  const std::vector< std::string > allColVar    = { "Muon_Eta" , "N_DS_Entries" , "N_WToMu_RAW" , "N_DYToMu" , "N_WToTauToMu" , "N_DYToTauToMu" , "N_TTbarToMu" , "N_QCDToMu" }; // , "TEST_FIT_BCChi2"
  std::vector< std::string > allColTitle1 = { "$\\etaMuLAB$ Range" , "Total" , "Signal" , "\\DYToMuMu" , "\\WToTauNu" , "\\DYToTauTau" , "\\ttbar" , "QCD" }; // , "p-value ($\\chi^{2}$)"
  if (useEtaCM) { allColTitle1[0] = "$\\etaMuCM$ Range"; }
  //
  std::vector< std::string > colVar , colTitle1, colTitle2;
  for(uint i = 0; i < allColVar.size(); i++) {
    if (inVar.count(allColVar[i])>0) {
      if (allColVar[i]!="Muon_Eta" && inVar.at(allColVar[i]).at("Val")<0) continue;
      colVar.push_back(allColVar[i]);
      colTitle1.push_back(allColTitle1[i]);
    }
  }
  //
  texTable.push_back("\\begin{table}[htb!]");
  texTable.push_back("  \\centering");
  texTable.push_back("  \\resizebox{\\textwidth}{!}{");
  createYieldTable(texTable, colVar, colTitle1, inputVar, col, chg);
  texTable.push_back("  }");
  texTable.push_back(Form("  \\caption{%s}",
                          Form("Raw yields of %s and background processes, extracted from the nominal fits for each %s bin in the %s collision system. All analysis cuts are applied%s. All uncertainties shown are statistical only.",
                               (chg=="Pl" ? "\\WToMuNuPl" : "\\WToMuNuMi"),
                               (useEtaCM ? "muon $\\etaMuCM$" : "$\\etaMuLAB$"),
                               (col=="PA" ? "combined \\pPb and \\Pbp" : (col=="pPb" ? "\\pPb" : (col=="Pbp" ? "\\Pbp" : "????"))),
                               (pTCUT!="" ? Form(" including the muon %s cut", pTCUT.c_str()) : "")
                               )
                          )
                     );
  texTable.push_back(Form("  \\label{tab:RawYields_WToMu%s_%s}", chg.c_str(), col.c_str()));
  texTable.push_back("\\end{table}");
  //
  for (const auto& row : texTable) { file << row << std::endl; }
  file << std::endl; file << std::endl;
  //
};


void makeCorrYieldsTable(std::ofstream& file, const VarBinMap& inputVar, const std::string& col, const std::string chg)
{
  //
  const auto& inVar = inputVar.at(col).begin()->second.at(chg);
  bool useEtaCM = false; if (inVar.count("useEtaCM")>0 && inVar.at("useEtaCM").at("Val")>0) { useEtaCM = true; }
  //
  // Extract kinematic info
  std::string pTCUT = "";
  if (inVar.count("Muon_Pt")>0) {
    if (inVar.at("Muon_Pt").at("Min") >0.0 && inVar.at("Muon_Pt").at("Max")>=1000.0) { pTCUT = Form("$p_{T} > %.0f$~GeV/c", inVar.at("Muon_Pt").at("Min")); }
    if (inVar.at("Muon_Pt").at("Min")<=0.0 && inVar.at("Muon_Pt").at("Max") <1000.0) { pTCUT = Form("$p_{T} < %.0f$~GeV/c", inVar.at("Muon_Pt").at("Max")); }
    if (inVar.at("Muon_Pt").at("Min") >0.0 && inVar.at("Muon_Pt").at("Max") <1000.0) { pTCUT = Form("$%.0f < p_{T} < %.0f$~GeV/c", inVar.at("Muon_Pt").at("Min"), inVar.at("Muon_Pt").at("Max")); }
  }
  //
  // Initialize the latex table
  std::vector< std::string > texTable;
  //
  // Determine number of columns
  const std::vector< std::string > allColVar    = { "Muon_Eta" , "N_WToMu_RAW" , "Acceptance_MC" , "Efficiency_MC" , "Efficiency_TnP" , "N_WToMu" };
  std::vector< std::string > allColTitle1 = { "$\\etaMuLAB$ Range" , "Raw Yield"  , "MC Acceptance ($\\%$)" , "MC Efficiency ($\\%$)" , "Efficiency ($\\%$)" , "Corrected Yield" };
  if (useEtaCM) { allColTitle1[0] = "$\\etaMuCM$ Range"; }
  std::vector< std::string > colVar , colTitle1, colTitle2;
  for(uint i = 0; i < allColVar.size(); i++) {
    if (inVar.count(allColVar[i])>0) {
      colVar.push_back(allColVar[i]);
      colTitle1.push_back(allColTitle1[i]);
    }
  }
  //
  texTable.push_back("\\begin{table}[htb!]");
  texTable.push_back("  \\centering");
  //texTable.push_back("  \\resizebox{\\textwidth}{!}{");
  createYieldTable(texTable, colVar, colTitle1, inputVar, col, chg);
  //texTable.push_back("  }");
  texTable.push_back(Form("  \\caption{%s}",
                          Form("Corrected yields of %s, given for each %s bin in the %s collision system. All analysis cuts are applied%s.%s All uncertainties shown are statistical only.",
                               (chg=="Pl" ? "\\WToMuNuPl" : "\\WToMuNuMi"),
                               (useEtaCM ? "muon $\\etaMuCM$" : "$\\etaMuLAB$"),
                               (col=="PA" ? "combined \\pPb and \\Pbp" : (col=="pPb" ? "\\pPb" : (col=="Pbp" ? "\\Pbp" : "????"))),
                               (pTCUT!="" ? Form(" including the muon %s cut", pTCUT.c_str()) : ""),
                               (inVar.count("Efficiency_TnP")>0 ? " The muon efficiency has been corrected by applying the Tag and Probe scale factors, HF energy weights and vector boson \\pt weights, event by event." : "")
                               )
                          )
                     );
  texTable.push_back(Form("  \\label{tab:CorrYields_WToMu%s_%s}", chg.c_str(), col.c_str()));
  texTable.push_back("\\end{table}");
  //
  for (const auto& row : texTable) { file << row << std::endl; }
  file << std::endl; file << std::endl;
  //
};


bool printYieldsTables(const VarBinMap& inputVar, const std::string& outDir)
{
  //
  std::cout << "[INFO] Filling the yield tables" << std::endl;
  //
  for (const auto& c : inputVar) {
    // Create Output Directory
    const std::string tableDir = outDir + "/Tables/SignalExtraction/" + c.first;
    makeDir(tableDir);
    // Create Output Files for RAW Yields
    const std::string fileName_RAW = "rawYield_" + c.first;
    std::ofstream file_RAW((tableDir + "/" + fileName_RAW + ".tex").c_str());
    if (file_RAW.is_open()==false) { std::cout << "[ERROR] File " << fileName_RAW << " was not created!" << std::endl; return false; }
    //
    makeRawYieldsTable(file_RAW, inputVar, c.first, "Mi");
    makeRawYieldsTable(file_RAW, inputVar, c.first, "Pl");
    //
    // Create Output Files for CORRECTED Yields
    const std::string fileName_CORR = "corrYield_" + c.first;
    std::ofstream file_CORR((tableDir + "/" + fileName_CORR + ".tex").c_str());
    if (file_CORR.is_open()==false) { std::cout << "[ERROR] File " << fileName_CORR << " was not created!" << std::endl; return false; }
    //
    makeCorrYieldsTable(file_CORR, inputVar, c.first, "Mi");
    makeCorrYieldsTable(file_CORR, inputVar, c.first, "Pl");
  }
  //
  return true;
};


void createResultTable(std::vector< std::string >& texTable, const std::vector< std::string >& colVar, const std::vector< std::string >& colTitle1,
                       const std::vector< std::string >& colChg, const std::vector< std::string >& colEta,
                       const VarBinMap& inputVar, const BinPentaMap& resultVar, const std::string& col)
{
  //
  const uint nCol = colVar.size();
  //
  texTable.push_back("  \\renewcommand{\\arraystretch}{1.5}");
  texTable.push_back(Form("  \\begin{tabular}{|c|*%dc|}", (nCol-1)));
  texTable.push_back("    \\hline");
  std::string tmp;
  tmp = ("    ");
  for (uint i = 0; i < nCol; i++) {
    tmp += colTitle1[i];
    if (i<(nCol-1)) { tmp += " & "; }
    else { tmp += "\\\\"; }
  }
  texTable.push_back(tmp);
  texTable.push_back("    \\hline\\hline");
  //
  std::string varLbl = "Cross_Section"; for (const auto& v : colVar) { if (v=="ForwardBackward_Ratio") { varLbl = v; break; } }
  if (resultVar.at(col).at("Pl").count(varLbl)==0) { std::cout << "[INFO] Variable " << varLbl << " is missing, will skip the table!" << std::endl; return; }
  //
  for (const auto& b : resultVar.at(col).at("Pl").at(varLbl).at("Val")) {
    tmp = ("    ");
    for (uint i = 0; i < nCol; i++) {
      auto bin = b.first;
      if (colEta[i]=="Inv") { bin = anabin<0>(-1.0*bin.etabin().high() , ( (bin.etabin().low() == 0.0) ? 0.0 : -1.0*bin.etabin().low() )); }
      const auto&    v = colVar[i];
      const auto&  chg = colChg[i];
      if (resultVar.at(col).count(chg)>0) {
	std::string val;
	if (v=="Muon_Eta") {
	  const double min = b.first.etabin().low();
	  const double max = b.first.etabin().high();
	  val = Form("%s%.2f , %s%.2f", sgn(min), abs(min), sgn(max) , abs(max));
	}
	else if (v=="Cross_Section") {
	  const auto& rVar = resultVar.at(col).at(chg).at(v);
	  const double errLow = rVar.at("Err_Stat_Low").at(b.first);
	  const double errHi  = rVar.at("Err_Stat_High").at(b.first);
	  const double rVal   = rVar.at("Val").at(b.first);
	  if (errLow==errHi) { val = Form("$%.2f \\pm %.2f$", rVal, errLow ); }
	  else { val = Form("$%.2f + %.2f - %.2f$", rVal, errHi, errLow ); }
	}
	else if (v=="ForwardBackward_Ratio" || v=="Charge_Asymmetry") {
	  const auto& rVar = resultVar.at(col).at(chg).at(v);
	  const double errLow = rVar.at("Err_Stat_Low").at(b.first);
	  const double errHi  = rVar.at("Err_Stat_High").at(b.first);
	  const double rVal   = rVar.at("Val").at(b.first);
	  if (errLow==errHi) { val = Form("$%.3f \\pm %.3f$", rVal, errLow ); }
	  else { val = Form("$%.3f + %.3f - %.3f$", rVal, errHi, errLow ); }
	}
	else {
	  const auto& iVar = inputVar.at(col).at(bin).at(chg).at(v);
	  const double errLow = iVar.at("Err_Stat_Low");
	  const double errHi  = iVar.at("Err_Stat_High");
	  const double iVal   = iVar.at("Val");
	  if (errLow==errHi) { val = Form("$%.2f \\pm %.2f$", iVal, errLow ); }
	  else { val = Form("$%.2f + %.2f - %.2f$", iVal, errHi, errLow ); }
	}
	tmp += val;
      }
      if (i<(nCol-1)) { tmp += " & "; } else { tmp += "\\\\"; }
    }
    texTable.push_back(tmp);
    texTable.push_back("    \\hline");
  }
  //
  texTable.push_back("  \\end{tabular}");
};


void makeCrossSectionTable(std::ofstream& file, const BinPentaMap& resultVar, const VarBinMap& inputVar, const std::string& col)
{
  //
  const auto& inVar = inputVar.at(col).begin()->second.at("Pl");
  bool useEtaCM = false; if (inVar.count("useEtaCM")>0 && inVar.at("useEtaCM").at("Val")>0) { useEtaCM = true; }
  //
  // Extract kinematic info
  std::string pTCUT = "";
  if (inVar.count("Muon_Pt")>0) {
    if (inVar.at("Muon_Pt").at("Min") >0.0 && inVar.at("Muon_Pt").at("Max")>=1000.0) { pTCUT = Form("$p_{T} > %.0f$~GeV/c", inVar.at("Muon_Pt").at("Min")); }
    if (inVar.at("Muon_Pt").at("Min")<=0.0 && inVar.at("Muon_Pt").at("Max") <1000.0) { pTCUT = Form("$p_{T} < %.0f$~GeV/c", inVar.at("Muon_Pt").at("Max")); }
    if (inVar.at("Muon_Pt").at("Min") >0.0 && inVar.at("Muon_Pt").at("Max") <1000.0) { pTCUT = Form("$%.0f < p_{T} < %.0f$~GeV/c", inVar.at("Muon_Pt").at("Min"), inVar.at("Muon_Pt").at("Max")); }
  }
  const double lumiVal = ( (inVar.count("Luminosity")>0) ? inVar.at("Luminosity").at("Val")          : -1. );
  const double lumiErr = ( (inVar.count("Luminosity")>0) ? inVar.at("Luminosity").at("Val")*LUMIUNC_ : -1. ); // 5% error
  //
  // Initialize the latex table
  std::vector< std::string > texTable;
  //
  // Determine number of columns
  const std::vector< std::string > colChg = { "Mi" , "Mi" , "Mi" , "Pl" , "Pl" };
  const std::vector< std::string > colEta = { "" , "" , "" , "" , "" };
  const std::vector< std::string > colVar = { "Muon_Eta" , "N_WToMu" , "Cross_Section" , "N_WToMu" , "Cross_Section" };
  std::vector< std::string >    colTitle1 = { "$\\etaMuLAB$ Range" , "$\\W^{-} $ Yield"  , "" , "$\\W^{+} $ Yield" , "" };
  if (useEtaCM) { colTitle1[0] = "$\\etaMuCM$ Range"; }
  colTitle1[2] = Form("$%s$", formatResultVarName("Cross_Section", useEtaCM, false, true, "Mi").c_str());
  colTitle1[4] = Form("$%s$", formatResultVarName("Cross_Section", useEtaCM, false, true, "Pl").c_str());
  //
  texTable.push_back("\\begin{table}[htb!]");
  texTable.push_back("  \\centering");
  texTable.push_back("  \\resizebox{\\textwidth}{!}{");
  createResultTable(texTable, colVar, colTitle1, colChg, colEta, inputVar, resultVar, col);
  texTable.push_back("  }");
  texTable.push_back(Form("  \\caption{%s}",
                          Form("Differential cross section of \\WToMuNu, given for each %s bin in the %s collision system. The total integrated luminosity of the sample corresponds to $%.1f \\pm %.1f$~\\nbinv. All analysis cuts are applied%s. All uncertainties shown are statistical only.",
                               (useEtaCM ? "muon $\\etaMuCM$" : "$\\etaMuLAB$"),
                               (col=="PA" ? "combined \\pPb and \\Pbp" : (col=="pPb" ? "\\pPb" : (col=="Pbp" ? "\\Pbp" : "????"))),
                               lumiVal, lumiErr,
                               (pTCUT!="" ? Form(" including the muon %s cut", pTCUT.c_str()) : "")
                               )
                          )
                     );
  texTable.push_back(Form("  \\label{tab:CrossSection_WToMu_%s}", col.c_str()));
  texTable.push_back("\\end{table}");
  //
  for (const auto& row : texTable) { file << row << std::endl; }
  file << std::endl; file << std::endl;
  //
};


void makeForwardBackwardTable(std::ofstream& file, const BinPentaMap& resultVar, const VarBinMap& inputVar, const std::string& col)
{
  //
  const auto& inVar = inputVar.at(col).begin()->second.at("Pl");
  bool useEtaCM = false; if (inVar.count("useEtaCM")>0 && inVar.at("useEtaCM").at("Val")>0) { useEtaCM = true; }
  //
  // Extract kinematic info
  std::string pTCUT = "";
  if (inVar.count("Muon_Pt")>0) {
    if (inVar.at("Muon_Pt").at("Min") >0.0 && inVar.at("Muon_Pt").at("Max")>=1000.0) { pTCUT = Form("$p_{T} > %.0f$~GeV/c", inVar.at("Muon_Pt").at("Min")); }
    if (inVar.at("Muon_Pt").at("Min")<=0.0 && inVar.at("Muon_Pt").at("Max") <1000.0) { pTCUT = Form("$p_{T} < %.0f$~GeV/c", inVar.at("Muon_Pt").at("Max")); }
    if (inVar.at("Muon_Pt").at("Min") >0.0 && inVar.at("Muon_Pt").at("Max") <1000.0) { pTCUT = Form("$%.0f < p_{T} < %.0f$~GeV/c", inVar.at("Muon_Pt").at("Min"), inVar.at("Muon_Pt").at("Max")); }
  }
  //
  // Initialize the latex table
  std::vector< std::string > texTable;
  //
  // Determine number of columns
  const std::vector< std::string > colChg = { "Pl" , "Mi" , "Mi" , "Pl" , "Pl" , "Mi" , "Pl" , "" };
  const std::vector< std::string > colEta = { "" , "" , "Inv" , "" , "Inv" , "" , "" , "" };
  const std::vector< std::string > colVar = { "Muon_Eta" , "N_WToMu" , "N_WToMu" , "N_WToMu" , "N_WToMu" , "ForwardBackward_Ratio" , "ForwardBackward_Ratio" , "ForwardBackward_Ratio" };
  std::vector< std::string >    colTitle1 = { "$\\etaMuLAB$ Range" , "$\\W^{-} $ Fwd Yield" , "$\\W^{-} $ Bwd Yield" , "$\\W^{+} $ Fwd Yield" , "$\\W^{+} $ Bwd Yield"  , "" , "" , "" };
  if (useEtaCM) { colTitle1[0] = "$\\etaMuCM$ Range"; }
  colTitle1[5] = Form("$\\W^{-} %s$", formatResultVarName("ForwardBackward_Ratio", useEtaCM, false, true).c_str());
  colTitle1[6] = Form("$\\W^{+} %s$", formatResultVarName("ForwardBackward_Ratio", useEtaCM, false, true).c_str());
  colTitle1[7] = Form("$\\W %s$", formatResultVarName("ForwardBackward_Ratio", useEtaCM, false, true).c_str());
  //
  texTable.push_back("\\begin{table}[htb!]");
  texTable.push_back("  \\centering");
  texTable.push_back("  \\resizebox{\\textwidth}{!}{");
  createResultTable(texTable, colVar, colTitle1, colChg, colEta, inputVar, resultVar, col);
  texTable.push_back("  }");
  texTable.push_back(Form("  \\caption{%s}",
                          Form("Forward-Backward ratio of \\WToMuNu, given for each %s bin in the %s collision system. All analysis cuts are applied%s. All uncertainties shown are statistical only.",
                               (useEtaCM ? "muon $\\etaMuCM$" : "$\\etaMuLAB$"),
                               (col=="PA" ? "combined \\pPb and \\Pbp" : (col=="pPb" ? "\\pPb" : (col=="Pbp" ? "\\Pbp" : "????"))),
                               (pTCUT!="" ? Form(" including the muon %s cut", pTCUT.c_str()) : "")
                               )
                          )
                     );
  texTable.push_back(Form("  \\label{tab:ForwardBackwardRatio_WToMu_%s}", col.c_str()));
  texTable.push_back("\\end{table}");
  //
  for (const auto& row : texTable) { file << row << std::endl; }
  file << std::endl; file << std::endl;
  //
};


void makeChargeAsymmetryTable(std::ofstream& file, const BinPentaMap& resultVar, const VarBinMap& inputVar, const std::string& col)
{
  //
  const auto& inVar = inputVar.at(col).begin()->second.at("Pl");
  bool useEtaCM = false; if (inVar.count("useEtaCM")>0 && inVar.at("useEtaCM").at("Val")>0) { useEtaCM = true; }
  //
  // Extract kinematic info
  std::string pTCUT = "";
  if (inVar.count("Muon_Pt")>0) {
    if (inVar.at("Muon_Pt").at("Min") >0.0 && inVar.at("Muon_Pt").at("Max")>=1000.0) { pTCUT = Form("$p_{T} > %.0f$~GeV/c", inVar.at("Muon_Pt").at("Min")); }
    if (inVar.at("Muon_Pt").at("Min")<=0.0 && inVar.at("Muon_Pt").at("Max") <1000.0) { pTCUT = Form("$p_{T} < %.0f$~GeV/c", inVar.at("Muon_Pt").at("Max")); }
    if (inVar.at("Muon_Pt").at("Min") >0.0 && inVar.at("Muon_Pt").at("Max") <1000.0) { pTCUT = Form("$%.0f < p_{T} < %.0f$~GeV/c", inVar.at("Muon_Pt").at("Min"), inVar.at("Muon_Pt").at("Max")); }
  }
  //
  // Initialize the latex table
  std::vector< std::string > texTable;
  //
  // Determine number of columns
  const std::vector< std::string > colChg = { "Pl" , "Pl" , "Mi" , "" };
  const std::vector< std::string > colEta = { "" , "" , "" , "" };
  const std::vector< std::string > colVar = { "Muon_Eta" , "N_WToMu" , "N_WToMu" , "Charge_Asymmetry" };
  std::vector< std::string >    colTitle1 = { "$\\etaMuLAB$ Range" , "$\\W^{+} $ Yield" , "$\\W^{-} $ Yield" , "" };
  if (useEtaCM) { colTitle1[0] = "$\\etaMuCM$ Range"; }
  colTitle1[3] = Form("$%s$", formatResultVarName("Charge_Asymmetry", useEtaCM, false, true).c_str());
  //
  texTable.push_back("\\begin{table}[htb!]");
  texTable.push_back("  \\centering");
  //texTable.push_back("  \\resizebox{\\textwidth}{!}{");
  createResultTable(texTable, colVar, colTitle1, colChg, colEta, inputVar, resultVar, col);
  //texTable.push_back("  }");
  texTable.push_back(Form("  \\caption{%s}",
                          Form("\\W charge asymmetry, given for each %s bin in the %s collision system. All analysis cuts are applied%s. All uncertainties shown are statistical only.",
                               (useEtaCM ? "muon $\\etaMuCM$" : "$\\etaMuLAB$"),
                               (col=="PA" ? "combined \\pPb and \\Pbp" : (col=="pPb" ? "\\pPb" : (col=="Pbp" ? "\\Pbp" : "????"))),
                               (pTCUT!="" ? Form(" including the muon %s cut", pTCUT.c_str()) : "")
                               )
                          )
                     );
  texTable.push_back(Form("  \\label{tab:ChargeAsymmetry_WToMu_%s}", col.c_str()));
  texTable.push_back("\\end{table}");
  //
  for (const auto& row : texTable) { file << row << std::endl; }
  file << std::endl; file << std::endl;
  //
};


bool printResultsTables(const BinPentaMap& resultVar, const VarBinMap& inputVar, const std::string& outDir)
{
  //
  std::cout << "[INFO] Filling the result tables" << std::endl;
  //
  // Fill the tables
  //
  for (const auto& c : resultVar) {
    //
    // Create Output Directory
    const std::string tableDir = outDir + "/Tables/Results/" + c.first;
    makeDir(tableDir);
    // Create Output Files for Cross Sections
    const std::string fileName_XSEC = "crossSection_" + c.first;
    std::ofstream file_XSEC((tableDir + "/" + fileName_XSEC + ".tex").c_str());
    if (file_XSEC.is_open()==false) { std::cout << "[ERROR] File " << fileName_XSEC << " was not created!" << std::endl; return false; }
    makeCrossSectionTable(file_XSEC, resultVar, inputVar, c.first);
    //
    // Create Output Files for Forward-Backward Ratios
    const std::string fileName_FB = "forwardBackward_" + c.first;
    std::ofstream file_FB((tableDir + "/" + fileName_FB + ".tex").c_str());
    if (file_FB.is_open()==false) { std::cout << "[ERROR] File " << fileName_FB << " was not created!" << std::endl; return false; }
    makeForwardBackwardTable(file_FB, resultVar, inputVar, c.first);
    //
    // Create Output Files for Muon Charge Asymmetry Ratios
    const std::string fileName_CA = "chargeAsymmetry_" + c.first;
    std::ofstream file_CA((tableDir + "/" + fileName_CA + ".tex").c_str());
    if (file_CA.is_open()==false) { std::cout << "[ERROR] File " << fileName_CA << " was not created!" << std::endl; return false; }
    makeChargeAsymmetryTable(file_CA, resultVar, inputVar, c.first);
  }
  //
  return true;
};


void createSystematicTable(std::vector< std::string >& texTable, const std::vector< std::string >& colVar, const std::vector< std::string >& colTitle1,
                           const std::vector< std::string >& colChg, const BinPentaMap& varMap, const std::string& col)
{
  //
  const uint nCol = colVar.size();
  //
  texTable.push_back("  \\renewcommand{\\arraystretch}{1.5}");
  texTable.push_back(Form("  \\begin{tabular}{|c|*%dc|}", (nCol-1)));
  texTable.push_back("    \\hline");
  std::string tmp;
  tmp = ("    ");
  for (uint i = 0; i < nCol; i++) {
    tmp += colTitle1[i];
    if (i<(nCol-1)) { tmp += " & "; }
    else { tmp += "\\\\"; }
  }
  texTable.push_back(tmp);
  texTable.push_back("    \\hline\\hline");
  //
  std::string varLbl = "Cross_Section"; for (const auto& v : colVar) { if (v=="ForwardBackward_Ratio") { varLbl = v; break; } }
  //
  for (const auto& b : varMap.at(col).at("Pl").at(varLbl).begin()->second) {
    tmp = ("    ");
    for (uint i = 0; i < nCol; i++) {
      //
      const auto&    v = colVar[i];
      const auto&  chg = colChg[i];
      std::string val;
      //
      if (v=="Muon_Eta") {
        const double min = b.first.etabin().low();
        const double max = b.first.etabin().high();
        val = Form("%s%.2f , %s%.2f", sgn(min), abs(min), sgn(max) , abs(max));
      }
      else if (v=="Cross_Section" || v=="ForwardBackward_Ratio" || v=="Charge_Asymmetry") {
	const double statErrLow = varMap.at(col).at(chg).at(v).at("Err_Stat_Low").at(b.first);
	const double statErrHi  = varMap.at(col).at(chg).at(v).at("Err_Stat_High").at(b.first);
	const double systErrLow = varMap.at(col).at(chg).at(v).at("Err_Syst_Low").at(b.first);
	const double systErrHi  = varMap.at(col).at(chg).at(v).at("Err_Syst_High").at(b.first);
	const double rVal       = varMap.at(col).at(chg).at(v).at("Val").at(b.first);
	if (v=="Cross_Section") {
	  if (std::abs(statErrLow-statErrHi)<0.001) { val = Form("$%.2f \\pm %.2f$ (stat)", rVal, statErrLow ); }
	  else { val = Form("$%.2f + %.2f - %.2f$ (stat)", rVal, statErrHi, statErrLow ); }
	  if (std::abs(systErrLow-systErrHi)<0.001) { val += Form("$ \\pm %.2f$ (syst)", systErrLow ); }
	  else { val += Form("$ + %.2f - %.2f$ (syst)", systErrHi, systErrLow ); }
	}
	else if (v=="ForwardBackward_Ratio" || v=="Charge_Asymmetry") {
	  if (std::abs(statErrLow-statErrHi)<0.001) { val = Form("$%.3f \\pm %.3f$ (stat)", rVal, statErrLow ); }
	  else { val = Form("$%.2f + %.3f - %.3f$ (stat)", rVal, statErrHi, statErrLow ); }
	  if (std::abs(systErrLow-systErrHi)<0.001) { val += Form("$ \\pm %.3f$ (syst)", systErrLow ); }
	  else { val += Form("$ + %.3f - %.3f$ (syst)", systErrHi, systErrLow ); }
	}
      }
      tmp += val;
      if (i<(nCol-1)) { tmp += " & "; } else { tmp += "\\\\"; }
    }
    texTable.push_back(tmp);
    texTable.push_back("    \\hline");
  }
  //
  texTable.push_back("  \\end{tabular}");
  //
};


void makeCrossSectionSystTable(std::ofstream& file, const BinPentaMap& varMap, const std::string& col, const bool useEtaCM = true)
{
  //
  // Initialize the latex table
  std::vector< std::string > texTable;
  //
  const std::string pTCUT = Form("$p_{T} > %.0f$~GeV/c", 25.0);
  double lumiVal = -1.0;
  if (col=="PA" ) { lumiVal = PA::LUMI::Data_pPb + PA::LUMI::Data_Pbp; }
  if (col=="pPb") { lumiVal = PA::LUMI::Data_pPb; }
  if (col=="Pbp") { lumiVal = PA::LUMI::Data_Pbp; }
  const double lumiErr = LUMIUNC_*lumiVal;
  //
  // Determine number of columns
  const std::vector< std::string > colChg = { "Mi" , "Mi" , "Pl" };
  const std::vector< std::string > colVar = { "Muon_Eta" , "Cross_Section" , "Cross_Section" };
  std::vector< std::string >    colTitle1 = { "$\\etaMuLAB$ Range" , "" , "" };
  if (useEtaCM) { colTitle1[0] = "$\\etaMuCM$ Range"; }
  colTitle1[1] = Form("$%s$", formatResultVarName("Cross_Section", useEtaCM, false, true, "Mi").c_str());
  colTitle1[2] = Form("$%s$", formatResultVarName("Cross_Section", useEtaCM, false, true, "Pl").c_str());
  //
  texTable.push_back("\\begin{table}[htb!]");
  texTable.push_back("  \\centering");
  //texTable.push_back("  \\resizebox{\\textwidth}{!}{");
  createSystematicTable(texTable, colVar, colTitle1, colChg, varMap, col);
  //texTable.push_back("  }");
  texTable.push_back(Form("  \\caption{%s}",
                          Form("Differential cross section of \\WToMuNu, given for each %s bin in the %s collision system. The total integrated luminosity of the sample corresponds to $%.1f \\pm %.1f$~\nbinv. All analysis cuts are applied%s. The global luminosity uncertainty of $\\pm$%.1f$\\%%$ is not included.",
                               (useEtaCM ? "muon $\\etaMuCM$" : "$\\etaMuLAB$"),
                               (col=="PA" ? "combined \\pPb and \\Pbp" : (col=="pPb" ? "\\pPb" : (col=="Pbp" ? "\\Pbp" : "????"))),
                               lumiVal, lumiErr,
                               (pTCUT!="" ? Form(" including the muon %s cut", pTCUT.c_str()) : ""),
                               lumiErr
                               )
                          )
                     );
  texTable.push_back(Form("  \\label{tab:CrossSectionSyst_WToMu_%s}", col.c_str()));
  texTable.push_back("\\end{table}");
  //
  for (const auto& row : texTable) { file << row << std::endl; }
  file << std::endl; file << std::endl;
  //
};


void makeForwardBackwardSystTable(std::ofstream& file, const BinPentaMap& varMap, const std::string& col, const bool useEtaCM = true)
{
  //
  // Initialize the latex table
  std::vector< std::string > texTable;
  //
  const std::string pTCUT = Form("$p_{T} > %.0f$~GeV/c", 25.0);
  //
  // Determine number of columns
  const std::vector< std::string > colChg = { "Pl" , "Mi" , "Pl" , "" };
  const std::vector< std::string > colVar = { "Muon_Eta" , "ForwardBackward_Ratio" , "ForwardBackward_Ratio" , "ForwardBackward_Ratio" };
  std::vector< std::string >    colTitle1 = { "$\\etaMuLAB$ Range" , "" , "" , "" };
  if (useEtaCM) { colTitle1[0] = "$\\etaMuCM$ Range"; }
  colTitle1[1] = Form("$%s$", formatResultVarName("ForwardBackward_Ratio", useEtaCM, false, true, "Mi").c_str());
  colTitle1[2] = Form("$%s$", formatResultVarName("ForwardBackward_Ratio", useEtaCM, false, true, "Pl").c_str());
  colTitle1[3] = Form("$%s$", formatResultVarName("ForwardBackward_Ratio", useEtaCM, false, true, "Inc").c_str());
  //
  texTable.push_back("\\begin{table}[htb!]");
  texTable.push_back("  \\centering");
  texTable.push_back("  \\resizebox{\\textwidth}{!}{");
  createSystematicTable(texTable, colVar, colTitle1, colChg, varMap, col);
  texTable.push_back("  }");
  texTable.push_back(Form("  \\caption{%s}",
                          Form("Forward-Backward ratio of \\WToMuNu, given for each %s bin in the %s collision system. All analysis cuts are applied%s.",
                               (useEtaCM ? "muon $\\etaMuCM$" : "$\\etaMuLAB$"),
                               (col=="PA" ? "combined \\pPb and \\Pbp" : (col=="pPb" ? "\\pPb" : (col=="Pbp" ? "\\Pbp" : "????"))),
                               (pTCUT!="" ? Form(" including the muon %s cut", pTCUT.c_str()) : "")
                               )
                          )
                     );
  texTable.push_back(Form("  \\label{tab:ForwardBackwardRatioSyst_WToMu_%s}", col.c_str()));
  texTable.push_back("\\end{table}");
  //
  for (const auto& row : texTable) { file << row << std::endl; }
  file << std::endl; file << std::endl;
  //
};


void makeChargeAsymmetrySystTable(std::ofstream& file, const BinPentaMap& varMap, const std::string& col, const bool useEtaCM = true)
{
  //
  // Initialize the latex table
  std::vector< std::string > texTable;
  //
  const std::string pTCUT = Form("$p_{T} > %.0f$~GeV/c", 25.0);
  //
  // Determine number of columns
  const std::vector< std::string > colChg = { "Pl" , "" };
  const std::vector< std::string > colVar = { "Muon_Eta" , "Charge_Asymmetry" };
  std::vector< std::string >    colTitle1 = { "$\\etaMuLAB$ Range" , "" };
  if (useEtaCM) { colTitle1[0] = "$\\etaMuCM$ Range"; }
  colTitle1[1] = Form("$%s$", formatResultVarName("Charge_Asymmetry", useEtaCM, false, true).c_str());
  //
  texTable.push_back("\\begin{table}[htb!]");
  texTable.push_back("  \\centering");
  //texTable.push_back("  \\resizebox{\\textwidth}{!}{");
  createSystematicTable(texTable, colVar, colTitle1, colChg, varMap, col);
  //texTable.push_back("  }");
  texTable.push_back(Form("  \\caption{%s}",
                          Form("\\W charge asymmetry, given for each %s bin in the %s collision system. All analysis cuts are applied%s",
                               (useEtaCM ? "muon $\\etaMuCM$" : "$\\etaMuLAB$"),
                               (col=="PA" ? "combined \\pPb and \\Pbp" : (col=="pPb" ? "\\pPb" : (col=="Pbp" ? "\\Pbp" : "????"))),
                               (pTCUT!="" ? Form(" including the muon %s cut", pTCUT.c_str()) : "")
                               )
                          )
                     );
  texTable.push_back(Form("  \\label{tab:ChargeAsymmetrySyst_WToMu_%s}", col.c_str()));
  texTable.push_back("\\end{table}");
  //
  for (const auto& row : texTable) { file << row << std::endl; }
  file << std::endl; file << std::endl;
  //
};


bool printSystematicTables(const BinPentaMap& varMap , const std::string& outDir , const bool useEtaCM = true)
{
  //
  std::cout << "[INFO] Filling the systematic tables" << std::endl;
  //
  // Fill the tables
  //
  for (const auto& c : varMap) {
    //
    // Create Output Directory
    const std::string tableDir = outDir + "/Tables/Systematics/" + c.first;
    makeDir(tableDir);
    // Create Output Files for Cross Sections
    const std::string fileName_XSEC = "crossSection_" + c.first;
    std::ofstream file_XSEC((tableDir + "/" + fileName_XSEC + ".tex").c_str());
    if (file_XSEC.is_open()==false) { std::cout << "[ERROR] File " << fileName_XSEC << " was not created!" << std::endl; return false; }
    //
    makeCrossSectionSystTable(file_XSEC, varMap, c.first, useEtaCM);
    //
    // Create Output Files for Forward-Backward Ratios
    const std::string fileName_FB = "forwardBackward_" + c.first;
    std::ofstream file_FB((tableDir + "/" + fileName_FB + ".tex").c_str());
    if (file_FB.is_open()==false) { std::cout << "[ERROR] File " << fileName_FB << " was not created!" << std::endl; return false; }
    //
    makeForwardBackwardSystTable(file_FB, varMap, c.first, useEtaCM);
    //
    // Create Output Files for Muon Charge Asymmetry Ratios
    const std::string fileName_CA = "chargeAsymmetry_" + c.first;
    std::ofstream file_CA((tableDir + "/" + fileName_CA + ".tex").c_str());
    if (file_CA.is_open()==false) { std::cout << "[ERROR] File " << fileName_CA << " was not created!" << std::endl; return false; }
    //
    makeChargeAsymmetrySystTable(file_CA, varMap, c.first, useEtaCM);
  }
  //
  return true;
};


void createFullSystematicTable(std::vector< std::string >& texTable, const std::vector< std::string >& colVar, const std::vector< std::string >& colTitle1,
                               const std::vector< std::string >& colChg, const BinSextaMap& varMap, const std::string& col)
{
  //
  const uint nCol = colVar.size();
  //
  texTable.push_back("  \\renewcommand{\\arraystretch}{1.5}");
  texTable.push_back(Form("  \\begin{tabular}{|c|*%dc|}", (nCol-1)));
  texTable.push_back("    \\hline");
  std::string tmp;
  tmp = ("    ");
  for (uint i = 0; i < nCol; i++) {
    tmp += colTitle1[i];
    if (i<(nCol-1)) { tmp += " & "; }
    else { tmp += "\\\\"; }
  }
  texTable.push_back(tmp);
  texTable.push_back("    \\hline\\hline");
  //
  std::string varLbl = "Cross_Section"; for (const auto& v : colVar) { if (v=="ForwardBackward_Ratio") { varLbl = v; break; } }
  //
  for (const auto& sys : varMap) {
    //
    const auto& sysLbl = sys.first;
    if (sysLbl=="Nominal") continue;
    tmp = ("    ");
    for (uint i = 0; i < nCol; i++) {
      //
      const auto&   v = colVar[i];
      const auto& chg = colChg[i];
      std::string val;
      //
      if (v=="Title_Syst") {
        val = sysLbl; if (val.find("_")!=std::string::npos) { val.replace(val.find("_"), 1, " "); }
        if (val.find("_")!=std::string::npos) { val.replace(val.find("_"), 1, " "); }
      }
      else {
        const auto& value = sys.second.at(col).at(chg).at(v).at("Nom");
        const auto& errHi = sys.second.at(col).at(chg).at(v).at("Err_Syst_High");
        const auto& errLo = sys.second.at(col).at(chg).at(v).at("Err_Syst_Low");
        double eMax = -9999.0;
        for (const auto& b : errHi) {
          const double sVal = value.at(b.first);
          double sErrHi  = errHi.at(b.first);
          double sErrLow = errLo.at(b.first);
          if (v == "Cross_Section") { sErrHi  /= sVal;  sErrLow /= sVal; }
          if (sErrLow > eMax) { eMax = sErrLow; }
          if (sErrHi  > eMax) { eMax = sErrHi;  }
        }
        val = Form(" %.3f ", eMax );
      }
      tmp += val;
      if (i<(nCol-1)) { tmp += " & "; } else { tmp += "\\\\"; }
    }
    texTable.push_back(tmp);
    texTable.push_back("    \\hline");
  }
  if (true) { // ADD TOTAL SYSTEMATIC
    tmp = ("    ");
    for (uint i = 0; i < nCol; i++) {
      //
      const auto&   v = colVar[i];
      const auto& chg = colChg[i];
      std::string val;
      //
      if (v=="Title_Syst") { val = "Total Systematic Unc."; }
      else {
        const auto& value = varMap.at("Nominal").at(col).at(chg).at(v).at("Val");
        const auto& errHi = varMap.at("Nominal").at(col).at(chg).at(v).at("Err_Syst_High");
        const auto& errLo = varMap.at("Nominal").at(col).at(chg).at(v).at("Err_Syst_Low");
        double eMax = -9999.0;
        for (const auto& b : errHi) {
          const double sVal = value.at(b.first);
          double sErrHi  = errHi.at(b.first);
          double sErrLow = errLo.at(b.first);
          if (v == "Cross_Section") { sErrHi  /= sVal;  sErrLow /= sVal; }
          if (sErrLow > eMax) { eMax = sErrLow; }
          if (sErrHi  > eMax) { eMax = sErrHi;  }
        }
        val = Form(" %.3f ", eMax );
      }
      tmp += val;
      if (i<(nCol-1)) { tmp += " & "; } else { tmp += "\\\\"; }
    }
    texTable.push_back(tmp);
    texTable.push_back("    \\hline");
  }
  if (true) { // ADD STATISTICAL
    tmp = ("    ");
    for (uint i = 0; i < nCol; i++) {
      //
      const auto&   v = colVar[i];
      const auto& chg = colChg[i];
      std::string val;
      //
      if (v=="Title_Syst") { val = "Statistical  Unc."; }
      else {
        const auto& value = varMap.at("Nominal").at(col).at(chg).at(v).at("Val");
        const auto& errHi = varMap.at("Nominal").at(col).at(chg).at(v).at("Err_Stat_High");
        const auto& errLo = varMap.at("Nominal").at(col).at(chg).at(v).at("Err_Stat_Low");
        double eMax = -9999.0;
        for (const auto& b : errHi) {
          const double sVal = value.at(b.first);
          double sErrHi  = errHi.at(b.first);
          double sErrLow = errLo.at(b.first);
          if (v == "Cross_Section") { sErrHi  /= sVal;  sErrLow /= sVal; }
          if (sErrLow > eMax) { eMax = sErrLow; }
          if (sErrHi  > eMax) { eMax = sErrHi;  }
        }
        val = Form(" %.3f ", eMax );
      }
      tmp += val;
      if (i<(nCol-1)) { tmp += " & "; } else { tmp += "\\\\"; }
    }
    texTable.push_back(tmp);
    texTable.push_back("    \\hline");
  }
  //
  texTable.push_back("  \\end{tabular}");
  //
};


void makeFullSystematicTable(std::ofstream& file, const BinSextaMap& varMap, const std::string& col, const bool useEtaCM = true)
{
  //
  // Initialize the latex table
  std::vector< std::string > texTable;
  //
  // Determine number of columns
  const std::vector< std::string > colChg = { "" , "Mi" , "Pl" , "Mi" , "Pl" , "" , "" };
  const std::vector< std::string > colVar = { "Title_Syst" , "Cross_Section" , "Cross_Section" , "ForwardBackward_Ratio" , "ForwardBackward_Ratio" , "ForwardBackward_Ratio" , "Charge_Asymmetry" };
  std::vector< std::string >    colTitle1 = { "Systematic Variation" , "" , "" , "" , "" , "" , "" };
  colTitle1[1] = Form("$%s$", formatResultVarName("Cross_Section", useEtaCM, false, true, "Mi").c_str());
  colTitle1[2] = Form("$%s$", formatResultVarName("Cross_Section", useEtaCM, false, true, "Pl").c_str());
  //colTitle1[3] = Form("$\\W^{-} %s$", formatResultVarName("ForwardBackward_Ratio", useEtaCM, false, true).c_str());
  //colTitle1[4] = Form("$\\W^{+} %s$", formatResultVarName("ForwardBackward_Ratio", useEtaCM, false, true).c_str());
  //colTitle1[5] = Form("$\\W %s$", formatResultVarName("ForwardBackward_Ratio", useEtaCM, false, true).c_str());
  colTitle1[3] = Form("$%s$", formatResultVarName("ForwardBackward_Ratio", useEtaCM, false, true, "Mi").c_str());
  colTitle1[4] = Form("$%s$", formatResultVarName("ForwardBackward_Ratio", useEtaCM, false, true, "Pl").c_str());
  colTitle1[5] = Form("$%s$", formatResultVarName("ForwardBackward_Ratio", useEtaCM, false, true, "Inc").c_str());
  colTitle1[6] = Form("$%s$", formatResultVarName("Charge_Asymmetry", useEtaCM, false, true).c_str());
  //
  texTable.push_back("\\begin{table}[htb!]");
  texTable.push_back("  \\centering");
  texTable.push_back("  \\resizebox{\\textwidth}{!}{");
  createFullSystematicTable(texTable, colVar, colTitle1, colChg, varMap, col);
  texTable.push_back("  }");
  texTable.push_back(Form("  \\caption{%s}",
                          Form("Maximum error of the measured observables determined for each category in the %s collision system. The uncertainties of the differential cross sections are relative while for the asymmetries are absolute.",
                               (col=="PA" ? "combined \\pPb and \\Pbp" : (col=="pPb" ? "\\pPb" : (col=="Pbp" ? "\\Pbp" : "????")))
                               )
                          )
                     );
  texTable.push_back("  \\label{tab:Systematics}");
  texTable.push_back("\\end{table}");
  //
  for (const auto& row : texTable) { file << row << std::endl; }
  file << std::endl; file << std::endl;
  //
};


bool printFullSystematicTable(const BinSextaMap& varMap , const std::string& outDir , const bool useEtaCM = true)
{
  //
  std::cout << "[INFO] Filling the total systematic table" << std::endl;
  //
  // Fill the tables
  //
  for (const auto& c : varMap.at("Nominal")) {
    //
    // Create Output Directory
    const std::string tableDir = outDir + "/Tables/Systematics/" + c.first;
    makeDir(tableDir);
    // Create Output Files for Full Systematics
    const std::string fileName_SYST = "systematic_" + c.first;
    std::ofstream file_SYST((tableDir + "/" + fileName_SYST + ".tex").c_str());
    if (file_SYST.is_open()==false) { std::cout << "[ERROR] File " << fileName_SYST << " was not created!" << std::endl; return false; }
    //
    makeFullSystematicTable(file_SYST, varMap, c.first, useEtaCM);
  }
  //
  return true;
};


void createEffSystematicTable(std::vector< std::string >& texTable, const std::vector< std::string >& colVar, const std::vector< std::string >& colCor, const std::vector< std::string >& colTitle1,
                              const std::vector< std::string >& colChg, const BinSextaMapVec& systVar, const std::string& col)
{
  //
  const uint nCol = colVar.size();
  //
  texTable.push_back("  \\renewcommand{\\arraystretch}{1.5}");
  texTable.push_back(Form("  \\begin{tabular}{|c|*%dc|}", (nCol-1)));
  texTable.push_back("    \\hline");
  std::string tmp;
  tmp = ("    ");
  for (uint i = 0; i < nCol; i++) {
    tmp += colTitle1[i];
    if (i<(nCol-1)) { tmp += " & "; }
    else { tmp += "\\\\"; }
  }
  texTable.push_back(tmp);
  texTable.push_back("    \\hline\\hline");
  //
  std::string varLbl = "Cross_Section"; for (const auto& v : colVar) { if (v=="ForwardBackward_Ratio") { varLbl = v; break; } }
  //
  for (const auto& b : systVar.begin()->second[0].at(col).at("Pl").at(varLbl).begin()->second) {
    tmp = ("    ");
    for (uint i = 0; i < nCol; i++) {
      //
      const auto&  var = colVar[i];
      const auto&  cor = colCor[i];
      const auto&  chg = colChg[i];
      std::string val;
      //
      if (var=="Muon_Eta") {
        const double min = b.first.etabin().low();
        const double max = b.first.etabin().high();
        val = Form("%s%.2f , %s%.2f", sgn(min), abs(min), sgn(max) , abs(max));
      }
      else {
        double effErr = 0.0;
        const double& nom = systVar.begin()->second.back().at(col).at(chg).at(var).at("Nom").at(b.first);
        if (systVar.count(cor)>0) {
          const double& errLo = systVar.at(cor).back().at(col).at(chg).at(var).at("Err_Syst_Low" ).at(b.first);
          const double& errHi = systVar.at(cor).back().at(col).at(chg).at(var).at("Err_Syst_High").at(b.first);
          effErr = std::max( errLo , errHi );
        }
        else {
          double errLo=0.0 , errHi=0.0;
          for (const auto& lbl : systVar) {
            if (lbl.first.find(cor)!=std::string::npos) {
              errLo = sumErrors( errLo , lbl.second.back().at(col).at(chg).at(var).at("Err_Syst_Low" ).at(b.first) );
              errHi = sumErrors( errHi , lbl.second.back().at(col).at(chg).at(var).at("Err_Syst_High").at(b.first) );
            }
          }
          effErr = std::max( errLo , errHi );
        }
        if (var=="Cross_Section") { effErr /= nom; }
        val = Form("%.4f", effErr );
      }
      tmp += val;
      if (i<(nCol-1)) { tmp += " & "; } else { tmp += "\\\\"; }
    }
    texTable.push_back(tmp);
    texTable.push_back("    \\hline");
  }
  //
  texTable.push_back("  \\end{tabular}");
  //
};


void makeTnPSystematicTable(std::ofstream& file, const BinSextaMapVec& systVar, const std::string& col, const std::string& chg, const std::string& var, const bool useEtaCM = true)
{
  //
  // Initialize the latex table
  std::vector< std::string > texTable;
  //
  const std::string pTCUT = Form("$p_{T} > %.0f$~GeV/c", 25.0);
  //
  // Determine number of columns
  std::vector< std::string > colChg = { "" , chg , chg , chg , chg , chg , chg , chg , chg };
  std::vector< std::string > colCor = { "" , "TnP_Syst_MuID" , "TnP_Syst_Iso" , "TnP_Syst_Trig" , "TnP_Syst_BinIso" , "TnP_Syst_BinMuID" , "TnP_Syst_STA" , "TnP_Syst_PU" , "TnP_Syst" };
  std::vector< std::string > colVar = { "Muon_Eta" , var , var , var , var , var , var , var , var };
  std::vector< std::string > colTitle1 = { "$\\etaMuLAB$ Range" , "MuID" , "Iso" , "Trig" , "Iso Binned" , "MuID Binned" , "STA" , "PU" , "Total" };
  if (useEtaCM) { colTitle1[0] = "$\\etaMuCM$ Range"; }
  //
  texTable.push_back("\\begin{table}[htb!]");
  texTable.push_back("  \\centering");
  texTable.push_back("  \\resizebox{\\textwidth}{!}{");
  createEffSystematicTable(texTable, colVar, colCor, colTitle1, colChg, systVar, col);
  texTable.push_back("  }");
  texTable.push_back(Form("  \\caption{%s}",
                          Form("Systematic uncertainties of the %s corresponding to muon efficiency correction using the Tag and Probe method. The errors are shown as a function of the generated %s. %s",
                               (var=="Cross_Section" ? Form("%s differential cross section", (chg=="Pl" ? "$\\WToMuNuPl$" : "$\\WToMuNuMi$")) :
                                (var=="ForwardBackward_Ratio" ? Form("%s differential cross section", (chg=="" ? "$\\WToMuNu$" : (chg=="Pl" ? "$\\WToMuNuPl$" : "$\\WToMuNuMi$"))) :
                                 "$\\WToMuNu$ charge asymmetry")),
                               (useEtaCM ? "muon $\\etaMuCM$" : "$\\etaMuLAB$"),
                               ( (col=="PA") ? "The \\pPb and \\Pbp MC samples are combined as described in \\sect{sec:CombiningBeamDirection}. " : "")
                               )
                          )
                     );
  texTable.push_back(Form("  \\label{tab:tnpSystUncertainty_%s_%s_%s}", var.c_str(), chg.c_str(), col.c_str()));
  texTable.push_back("\\end{table}");
  //
  for (const auto& row : texTable) { file << row << std::endl; }
  file << std::endl; file << std::endl;
  //
  //
  // Initialize the latex table
  texTable.clear(); colChg.clear(); colCor.clear(); colVar.clear(); colTitle1.clear();
  //
  // Determine number of columns
  colChg = std::vector< std::string >({ "" , chg , chg , chg , chg });
  colCor = std::vector< std::string >({ "" , "TnP_Stat_MuID" , "TnP_Stat_Iso" , "TnP_Stat_Trig" , "TnP_Stat" });
  colVar = std::vector< std::string >({ "Muon_Eta" , var , var , var , var });
  colTitle1 = std::vector< std::string >({ "$\\etaMuLAB$ Range" , "MuID" , "Iso" , "Trig" , "Total" });
  if (useEtaCM) { colTitle1[0] = "$\\etaMuCM$ Range"; }
  //
  texTable.push_back("\\begin{table}[htb!]");
  texTable.push_back("  \\centering");
  //texTable.push_back("  \\resizebox{\\textwidth}{!}{");
  createEffSystematicTable(texTable, colVar, colCor, colTitle1, colChg, systVar, col);
  //texTable.push_back("  }");
  texTable.push_back(Form("  \\caption{%s}",
                          Form("Statistical uncertainties of the %s corresponding to muon efficiency correction using the Tag and Probe method. The errors are shown as a function of the generated %s. %s",
                               (var=="Cross_Section" ? Form("%s differential cross section", (chg=="Pl" ? "$\\WToMuNuPl$" : "$\\WToMuNuMi$")) :
                                (var=="ForwardBackward_Ratio" ? Form("%s differential cross section", (chg=="" ? "$\\WToMuNu$" : (chg=="Pl" ? "$\\WToMuNuPl$" : "$\\WToMuNuMi$"))) :
                                 "$\\WToMuNu$ charge asymmetry")),
                               (useEtaCM ? "muon $\\etaMuCM$" : "$\\etaMuLAB$"),
                               ( (col=="PA") ? "The \\pPb and \\Pbp MC samples are combined as described in \\sect{sec:CombiningBeamDirection}. " : "")
                               )
                          )
                     );
  texTable.push_back(Form("  \\label{tab:tnpStatUncertainty_%s_%s_%s}", var.c_str(), chg.c_str(), col.c_str()));
  texTable.push_back("\\end{table}");
  //
  for (const auto& row : texTable) { file << row << std::endl; }
  file << std::endl; file << std::endl;
  //
};


void makeMCSystematicTable(std::ofstream& file, const BinSextaMapVec& systVar, const std::string& col, const std::string& chg, const std::string& var, const bool useEtaCM = true)
{
  //
  // Initialize the latex table
  std::vector< std::string > texTable;
  //
  const std::string pTCUT = Form("$p_{T} > %.0f$~GeV/c", 25.0);
  //
  // Determine number of columns
  std::vector< std::string > colChg = { "" , chg , chg , chg , chg };
  std::vector< std::string > colCor = { "" , "MC_Statistics" , "MC_Syst_PDF" , "MC_Syst_Alpha" , "MC_Syst_Scale" };
  std::vector< std::string > colVar = { "Muon_Eta" , var , var , var , var };
  std::vector< std::string > colTitle1 = { "$\\etaMuLAB$ Range" , "MC Statistics" , "PDF" , "$\\Alpha_{s}$" , "($\\mu_{R},\\mu_{F}$) Scale" };
  if (useEtaCM) { colTitle1[0] = "$\\etaMuCM$ Range"; }
  //
  texTable.push_back("\\begin{table}[htb!]");
  texTable.push_back("  \\centering");
  texTable.push_back("  \\resizebox{\\textwidth}{!}{");
  createEffSystematicTable(texTable, colVar, colCor, colTitle1, colChg, systVar, col);
  texTable.push_back("  }");
  texTable.push_back(Form("  \\caption{%s}",
                          Form("Systematic uncertainties of the %s corresponding to the variation of the EPPS16+CT14 PDF, $\\alpha_{s}$ and ($\\mu_{R},\\mu_{F}$) scale variations. The errors are shown as a function of the generated %s. %s",
                               (var=="Cross_Section" ? Form("%s differential cross section", (chg=="Pl" ? "$\\WToMuNuPl$" : "$\\WToMuNuMi$")) :
                                (var=="ForwardBackward_Ratio" ? Form("%s differential cross section", (chg=="" ? "$\\WToMuNu$" : (chg=="Pl" ? "$\\WToMuNuPl$" : "$\\WToMuNuMi$"))) :
                                 "$\\WToMuNu$ charge asymmetry")),
                               (useEtaCM ? "muon $\\etaMuCM$" : "$\\etaMuLAB$"),
                               ( (col=="PA") ? "The \\pPb and \\Pbp MC samples are combined as described in \\sect{sec:CombiningBeamDirection}. " : "")
                               )
                          )
                     );
  texTable.push_back(Form("  \\label{tab:mcSystUncertainty_%s_%s_%s}", var.c_str(), chg.c_str(), col.c_str()));
  texTable.push_back("\\end{table}");
  //
  for (const auto& row : texTable) { file << row << std::endl; }
  file << std::endl; file << std::endl;
  //
};


bool printEffSystematicTables(const BinSextaMapVec& systVar , const std::string& outDir , const bool useEtaCM = true)
{
  //
  std::cout << "[INFO] Filling the TnP systematic tables" << std::endl;
  //
  for (const auto& c : systVar.begin()->second[0]) {
    // Create Output Directory
    const std::string tableDir = outDir + "/Tables/Systematics/" + c.first;
    makeDir(tableDir);
    // Create Output Files for TnP Uncertainties
    const std::string fileName_TnPUnc = "tnpUncertainty_" + c.first;
    std::ofstream file_TnPUnc((tableDir + "/" + fileName_TnPUnc + ".tex").c_str());
    if (file_TnPUnc.is_open()==false) { std::cout << "[ERROR] File " << fileName_TnPUnc << " was not created!" << std::endl; return false; }
    //
    makeTnPSystematicTable(file_TnPUnc, systVar, c.first, "Pl", "Cross_Section", useEtaCM);
    makeTnPSystematicTable(file_TnPUnc, systVar, c.first, "Mi", "Cross_Section", useEtaCM);
    //
    makeTnPSystematicTable(file_TnPUnc, systVar, c.first, "Pl", "ForwardBackward_Ratio", useEtaCM);
    makeTnPSystematicTable(file_TnPUnc, systVar, c.first, "Mi", "ForwardBackward_Ratio", useEtaCM);
    makeTnPSystematicTable(file_TnPUnc, systVar, c.first,   "", "ForwardBackward_Ratio", useEtaCM);
    //
    makeTnPSystematicTable(file_TnPUnc, systVar, c.first, "", "Charge_Asymmetry", useEtaCM);
    // Create Output Files for MC Uncertainties (PDF + Statistics)
    const std::string fileName_MCUnc = "mcUncertainty_" + c.first;
    std::ofstream file_MCUnc((tableDir + "/" + fileName_MCUnc + ".tex").c_str());
    if (file_MCUnc.is_open()==false) { std::cout << "[ERROR] File " << fileName_MCUnc << " was not created!" << std::endl; return false; }
    //
    makeMCSystematicTable(file_MCUnc, systVar, c.first, "Pl", "Cross_Section", useEtaCM);
    makeMCSystematicTable(file_MCUnc, systVar, c.first, "Mi", "Cross_Section", useEtaCM);
    //
    makeMCSystematicTable(file_MCUnc, systVar, c.first, "Pl", "ForwardBackward_Ratio", useEtaCM);
    makeMCSystematicTable(file_MCUnc, systVar, c.first, "Mi", "ForwardBackward_Ratio", useEtaCM);
    makeMCSystematicTable(file_MCUnc, systVar, c.first,   "", "ForwardBackward_Ratio", useEtaCM);
    //
    makeMCSystematicTable(file_MCUnc, systVar, c.first, "", "Charge_Asymmetry", useEtaCM);
    //
  }
  //
  return true;
};


void redrawBorder()
{
   gPad->Update();
   gPad->RedrawAxis();
   TLine l;
   l.DrawLine(gPad->GetUxmin(), gPad->GetUymax(), gPad->GetUxmax(), gPad->GetUymax());
   l.DrawLine(gPad->GetUxmax(), gPad->GetUymin(), gPad->GetUxmax(), gPad->GetUymax());
};


void makeCovarianceMatrix(std::ofstream& file , const WSDirMap& workDirNames, const BinSeptaMapVec& systVarMap, const BinPentaMap& varMap, const std::string& col, const std::string& var, const std::string& chg,
                          const std::string& outDir, const bool useEtaCM = true, const std::string type = "All")
{
  //
  // Initialize the latex table
  std::vector< std::string > texMatrix;
  //
  uint nCol = 0;
  const uint dCol = (var!="ForwardBackward_Ratio" ? 24 : 10);
  //
  // Fill the Nominal Value Vector
  //
  TVectorD nomVec;
  if (true) {
    std::vector<std::string> chgVec = { "Mi" , "Pl" };
    if (chg=="") { chgVec.clear(); chgVec.push_back(""); }
    const auto& sMap = systVarMap.begin()->second.begin()->second[0].at(col);
    const int nBin = sMap.at(chg).at(var).at("Nom").size();
    nomVec.ResizeTo(nBin*chgVec.size());
    ushort c1=0;
    for (const auto& chg : chgVec) {
      const auto& bMap = sMap.at(chg).at(var).at("Nom");
      ushort i=0;
      for (const auto& bin : bMap) {
        nomVec[i + nBin*c1] = bin.second;
        i++;
      }
      c1++;
    }
  }
  //
  // Fill the Systematic Covariance Matrices
  //
  CovMatrixDiMap covMatrixDiMap;
  //
  for (const auto& cat : systVarMap) {
    if (cat.first=="Nominal") continue;
    //
    auto& covMatrixMap = covMatrixDiMap[cat.first];
    //
    for (const auto& sys : cat.second) {
      //
      auto& covMatrix = covMatrixMap[sys.first];
      //
      const uint  nVar = systVarMap.at(cat.first).at(sys.first).size();
      const auto& vMap = systVarMap.at(cat.first).at(sys.first).back().at(col);
      const auto& sMap = systVarMap.at(cat.first).at(sys.first)[0].at(col);
      const auto& bMap = vMap.at(chg).at(var).at("Val");
      ushort corrType = 0;
      if (cat.first!="Efficiency") { corrType = workDirNames.at(cat.first).at(sys.first).second; }
      else if (sys.first!="MC_Statistics") { corrType = effTnPType_.at(sys.first).second; }
      nCol = bMap.size();
      nCol *= ( chg!="" ? 2.0 : 1.0 ); // If charged, increase the number of columns
      covMatrix.ResizeTo(nCol, nCol);
      std::vector<std::string> chgVec = { "Mi" , "Pl" };
      if (chg=="") { chgVec.clear(); chgVec.push_back(""); }
      //
      ushort c1=0;
      for (const auto& chg1 : chgVec) {
        ushort c2=0;
        for (const auto& chg2 : chgVec) {
          ushort i=0;
          for (const auto& bin1 : bMap) {
            ushort j=0;
            for (const auto& bin2 : bMap) {
              //
              // Determine the covariance
              const double errBin1  = vMap.at(chg1).at(var).at("Err_Syst_Low").at(bin1.first);
              const double errBin2  = vMap.at(chg2).at(var).at("Err_Syst_Low").at(bin2.first);
              const double covBin12 = ( errBin1 * errBin2 );
              //
              // Determine the sign
              const double nomBin1  = sMap.at(chg1).at(var).at("Nom").at(bin1.first);
              const double nomBin2  = sMap.at(chg2).at(var).at("Nom").at(bin2.first);
              const double varBin1  = sMap.at(chg1).at(var).at("Val").at(bin1.first);
              const double varBin2  = sMap.at(chg2).at(var).at("Val").at(bin2.first);
              double covSgn12 = 1.0;
              if (cat.first=="Efficiency" || cat.first=="EWK_Background") {
                if ( (((varBin1 - nomBin1)*(varBin2 - nomBin2)) < 0.0) && nVar <=3 ) { covSgn12 = -1.0; }
              }
              //
              bool fillCovariance = false;
              //
              // Case 1: Eta and Charge Uncorrelated
              if (corrType==0) {
                if (bin1.first==bin2.first && chg1==chg2) { fillCovariance = true; }
              }
              // Case 2: Eta Correlated and Charge Uncorrelated
              else if (corrType==1) {
                if (chg1==chg2) { fillCovariance = true; }
              }
              // Case 3: Eta UnCorrelated and Charge Correlated
              else if (corrType==2) {
                if (bin1.first==bin2.first) { fillCovariance = true; }
              }
              // Case 4: Eta Correlated and Charge Correlated
              else if (corrType==3 || corrType==5) { fillCovariance = true; }
              //
              if (fillCovariance) { covMatrix((i+bMap.size()*c1), (j+bMap.size()*c2)) = covSgn12*covBin12; }
              else { covMatrix((i+bMap.size()*c1), (j+bMap.size()*c2)) = 0.0; }
              //
              j++;
            }
            i++;
          }
          c2++;
        }
        c1++;
      }
    }
  }
  //
  // Fill the Statistical Covariance Matrices
  //
  if (true) {
    //
    auto& covMatrix = covMatrixDiMap["Statistical"]["Total"];
    //
    const auto& vMap = varMap.at(col);
    const auto& bMap = vMap.at(chg).at(var).at("Err_Stat_Low");
    nCol = bMap.size();
    nCol *= ( chg!="" ? 2.0 : 1.0 ); // If charged, increase the number of columns
    covMatrix.ResizeTo(nCol, nCol);
    std::vector<std::string> chgVec = { "Mi" , "Pl" };
    if (chg=="") { chgVec.clear(); chgVec.push_back(""); }
    //
    ushort c1=0;
    for (const auto& chg1 : chgVec) {
      ushort c2=0;
      for (const auto& chg2 : chgVec) {
        ushort i=0;
        for (const auto& bin1 : bMap) {
          ushort j=0;
          for (const auto& bin2 : bMap) {
            if (bin1.first==bin2.first && chg1==chg2) {
              // Determine the covariance
              const double errBin  = vMap.at(chg1).at(var).at("Err_Stat_Low").at(bin1.first);
              const double covBin = ( errBin * errBin );
              covMatrix((i+bMap.size()*c1), (j+bMap.size()*c2)) = covBin;
            }
            else { covMatrix((i+bMap.size()*c1), (j+bMap.size()*c2)) = 0.0; }
            //
            j++;
          }
          i++;
        }
        c2++;
      }
      c1++;
    }
  }
  //
  // Fill the Luminosity Covariance Matrices
  //
  if (var=="Cross_Section") {
    //
    auto& covMatrix = covMatrixDiMap["Luminosity"]["Total"];
    //
    const auto& vMap = varMap.at(col);
    const auto& bMap = vMap.at(chg).at(var).at("Err_Stat_Low");
    nCol = bMap.size();
    nCol *= ( chg!="" ? 2.0 : 1.0 ); // If charged, increase the number of columns
    covMatrix.ResizeTo(nCol, nCol);
    std::vector<std::string> chgVec = { "Mi" , "Pl" };
    if (chg=="") { chgVec.clear(); chgVec.push_back(""); }
    //
    ushort c1=0;
    for (const auto& chg1 : chgVec) {
      ushort c2=0;
      for (const auto& chg2 : chgVec) {
        ushort i=0;
        for (const auto& bin1 : bMap) {
          ushort j=0;
          for (const auto& bin2 : bMap) {
            // Determine the covariance
            const double errBin1  = LUMIUNC_*vMap.at(chg1).at(var).at("Val").at(bin1.first);
            const double errBin2  = LUMIUNC_*vMap.at(chg2).at(var).at("Val").at(bin2.first);
            const double covBin12 = ( errBin1 * errBin2 );
            covMatrix((i+bMap.size()*c1), (j+bMap.size()*c2)) = covBin12;
            //
            j++;
          }
          i++;
        }
        c2++;
      }
      c1++;
    }
  }
  //
  // Combine the Covariance Matrices
  //
  const auto covMatrixDiMapTmp = covMatrixDiMap;
  //
  auto& covMatrixTot = covMatrixDiMap["Total"]["Total"];
  covMatrixTot.ResizeTo(nCol, nCol);
  //
  for (const auto& cat : covMatrixDiMapTmp) {
    //
    auto& covMatrixCatTot = covMatrixDiMap.at(cat.first)["Total"];
    covMatrixCatTot.ResizeTo(nCol, nCol);
    for (const auto& sys : cat.second) { if (sys.first!="Total") { covMatrixCatTot += sys.second; } }
    covMatrixTot += covMatrixCatTot;
  }
  //
  // Transform from covanriance to correlation
  //
  auto corMatrixDiMap = covMatrixDiMap;
  for (const auto& cat : covMatrixDiMap) {
    for (const auto& sys : cat.second) {
            auto& corMatrix = corMatrixDiMap.at(cat.first).at(sys.first);
      const auto& covMatrix = covMatrixDiMap.at(cat.first).at(sys.first);
      for (ushort i=0; i<nCol; i++) {
        for (ushort j=0; j<nCol; j++) {
          corMatrix(i, j) = (covMatrix(i, j)/std::sqrt(covMatrix(i, i)*covMatrix(j, j)));
        }
      }
    }
  }
  //
  // Print the Covariance Matrices
  //
  texMatrix.push_back(Form("\\setcounter{MaxMatrixCols}{%d}", nCol));
  //
  for (const auto& cat : covMatrixDiMap) {
    if (type=="Total" && cat.first!="Total") continue;
    //
    texMatrix.push_back("\\verb|"+cat.first+"|");
    texMatrix.push_back("  "); texMatrix.push_back("  ");
    //
    for (const auto& sys : cat.second) {
      //
      const auto& covMatrix = sys.second;
      //
      texMatrix.push_back("\\verb|"+sys.first+"|");
      texMatrix.push_back("  "); texMatrix.push_back("  ");
      texMatrix.push_back("\\begin{displaymath}");
      texMatrix.push_back("\\resizebox{\\textwidth}{!}{$");
      texMatrix.push_back("  \\left[");
      texMatrix.push_back(Form("  \\begin{array}{@{}*{%d}c@{}}", nCol));
      //
      std::string tmp;
      for (ushort i=0; i<covMatrix.GetNrows(); i++) {
        tmp = ("    ");
        for (ushort j=0; j<covMatrix.GetNcols(); j++) {
          //
          std::string val;
          //
          if (j>=i) { val = Form("%.2f", covMatrix(i, j)); }
          else { val = " "; }
          //
          tmp += val;
          if (j<(covMatrix.GetNcols()-1)) { tmp += " & "; } else if (i<(covMatrix.GetNcols()-1)) { tmp += "\\\\"; }
        }
        texMatrix.push_back(tmp);
      }
      //
      texMatrix.push_back("  \\end{array}");
      texMatrix.push_back("  \\right]");
      texMatrix.push_back("$}");
      texMatrix.push_back("\\end{displaymath}");
      texMatrix.push_back("  "); texMatrix.push_back("  "); texMatrix.push_back("  "); texMatrix.push_back("  ");
    }
  }
  //
  for (const auto& row : texMatrix) { file << row << std::endl; }
  file << std::endl; file << std::endl;
  //
  // Draw the TH2 Plots
  //
  // Set Style
  setStyle();
  TGaxis::SetMaxDigits(2); // to display powers of 10
  //
  for (const auto& cat : covMatrixDiMap) {
    for (const auto& sys : cat.second) {
      //
      // Create Canvas
      TCanvas c("c", "c", 1000, 1000); c.cd();
      c.SetRightMargin(2.8);
      //
      auto& covMatrix = covMatrixDiMap.at(cat.first).at(sys.first);
      //
      TH2D graph("covMatrix", "", nCol, 0, nCol, std::ceil(nCol*1.4), 0, std::ceil(nCol*1.4));
      //
      for (ushort i=1; i<=nCol; i++) {
        for (ushort j=1; j<=nCol; j++) {
          graph.SetBinContent(j, (nCol+1-i), covMatrix(i-1,j-1) );
          graph.SetBinError(i, j, 0.0);
        }
      }
      //
      // Create the Text Info
      TLatex tex; tex.SetNDC(); tex.SetTextSize(0.035); float dy = 0;
      std::vector< std::string > textToPrint;
      std::string sampleLabel = "W #rightarrow #mu + #nu_{#mu}";
      if (chg == "Pl") { sampleLabel = "W^{#pm} #rightarrow #mu^{#pm} + #nu_{#mu}"; }
      textToPrint.push_back(sampleLabel);
      textToPrint.push_back("p^{#mu}_{T} > 25 GeV/c");
      //
      // Draw graph
      graph.Draw("COLZ");
      //
      // Format graph
      //
      // Set the Axis Titles
      std::string xLabel = "#eta^{#mu}";
      if (useEtaCM) { xLabel += "_{CM}"; }
      else { xLabel += "_{LAB}"; }
      std::string yLabel = formatResultVarName(var, useEtaCM, false, false, (chg=="Pl"?"Com":"Inc"));
      graph.SetTitle(Form("Covariance Matrix;%s;%s", "Analysis bin", "Analysis bin"));
      // X-axis
      graph.GetXaxis()->CenterTitle(kTRUE);
      graph.GetXaxis()->SetTitleOffset(0.6);
      graph.GetXaxis()->SetTitleSize(0.060);
      graph.GetXaxis()->SetLabelSize(0.0);
      // Y-axis
      graph.GetYaxis()->CenterTitle(kTRUE);
      graph.GetYaxis()->SetTitleOffset(0.7);
      graph.GetYaxis()->SetTitleSize(0.065);
      graph.GetYaxis()->SetLabelSize(0.0);
      // Z axis
      graph.GetZaxis()->SetLabelSize(0.030);
      //if (sys.first.find("Total")!=std::string::npos) { graph.GetZaxis()->SetRangeUser(-1.0 , 1.0); }
      c.Modified(); c.Update();
      // Draw the white box to hide under/overflow bins
      TBox box(c.GetFrame()->GetX1(), nCol, c.GetFrame()->GetX2(), c.GetFrame()->GetY2());
      box.SetFillColor(kWhite); box.SetLineColor(kWhite);
      box.Draw("same");
      redrawBorder();
      c.Modified(); c.Update();
      // Create line
      TLine line_1(nCol-dCol, 0      , nCol-dCol, nCol   ); line_1.SetLineWidth(3);
      TLine line_2(0      , nCol-dCol, nCol   , nCol-dCol); line_2.SetLineWidth(3);
      //
      if (nCol>dCol) {
        line_1.Draw("same");
        line_2.Draw("same");
      }
      //
      tex.SetTextSize(0.055); tex.DrawLatex(0.19, 0.85, textToPrint[0].c_str());
      tex.SetTextSize(0.058); tex.SetTextFont(61); tex.DrawLatex(0.75, 0.85, "CMS"); tex.SetTextFont(62);
      //tex.SetTextSize(0.044); tex.SetTextFont(52); tex.DrawLatex(0.66, 0.79, "Preliminary"); tex.SetTextFont(62);
      tex.SetTextSize(0.044); tex.SetTextFont(62); tex.DrawLatex(0.19, 0.79, yLabel.c_str()); tex.SetTextFont(62);
      if (textToPrint.size()>1) { tex.SetTextSize(0.035); tex.DrawLatex(0.19, 0.73, textToPrint[1].c_str()); }
      if (nCol>dCol) {
        tex.SetTextSize(0.035); tex.SetTextFont(62); tex.DrawLatex(0.42, 0.65, "Minus"); tex.SetTextFont(62);
        tex.SetTextSize(0.035); tex.SetTextFont(62); tex.DrawLatex(0.81, 0.65, "Plus" ); tex.SetTextFont(62);
        tex.SetTextSize(0.035); tex.SetTextFont(62); tex.SetTextAngle(90); tex.DrawLatex(0.21, 0.44, "Minus"); tex.SetTextFont(62); tex.SetTextAngle(0);
        tex.SetTextSize(0.035); tex.SetTextFont(62); tex.SetTextAngle(90); tex.DrawLatex(0.21, 0.15, "Plus" ); tex.SetTextFont(62); tex.SetTextAngle(0);
      }
      //
      // set the CMS style
      int option = 1118;
      CMS_lumi(&c, option, 33, "", false, 0.6, false);
      // Update
      c.Modified(); c.Update(); // Pure paranoia
      TPaletteAxis *palette = (TPaletteAxis*)graph.GetListOfFunctions()->FindObject("palette");
      palette->SetX2NDC(0.93);
      c.Modified();
      //
      // Save canvas
      //
      // Create Output Directory
      const std::string plotDir = outDir + "/Matrix/Covariance/" + col+"/Plot/"+var;
      makeDir(plotDir + "/png/");
      makeDir(plotDir + "/pdf/");
      makeDir(plotDir + "/root/");
      //
      // Save Canvas
      const std::string name = Form("covMatrix_WToMu%s_%s_%s_%s_%s", chg.c_str(), col.c_str(), var.c_str(), cat.first.c_str(), sys.first.c_str());
      c.SaveAs(( plotDir + "/png/"  + name + ".png"  ).c_str());
      c.SaveAs(( plotDir + "/pdf/"  + name + ".pdf"  ).c_str());
      c.SaveAs(( plotDir + "/root/" + name + ".root" ).c_str());
      //
      // Save Total
      if (cat.first=="Total" && sys.first=="Total") {
        const std::string dataDir = outDir + "/Matrix/Covariance/" + col+"/Data/"+var;
        makeDir(dataDir);
        const std::string fName = ( dataDir + "/"  + name + ".root"  );
        TFile f(fName.c_str(), "RECREATE"); f.cd();
        covMatrix.Write("covMatrix");
        nomVec.Write("nomVec");
        f.Write();
        f.Close();
      }
      //
      // Clean up memory
      c.Clear(); c.Close();
    }
  }
  //
  for (const auto& cat : corMatrixDiMap) {
    for (const auto& sys : cat.second) {
      //
      // Create Canvas
      TCanvas c("c", "c", 1000, 1000); c.cd();
      c.SetRightMargin(2.8);
      //
      auto& corMatrix = corMatrixDiMap.at(cat.first).at(sys.first);
      //
      TH2D graph("covMatrix", "", nCol, 0, nCol, std::ceil(nCol*1.4), 0, std::ceil(nCol*1.4));
      //
      for (ushort i=1; i<=nCol; i++) {
        for (ushort j=1; j<=nCol; j++) {
          graph.SetBinContent(j, (nCol+1-i), corMatrix(i-1,j-1) );
          graph.SetBinError(i, j, 0.0);
        }
      }
      //
      // Create the Text Info
      TLatex tex; tex.SetNDC(); tex.SetTextSize(0.035); float dy = 0;
      std::vector< std::string > textToPrint;
      std::string sampleLabel = "W #rightarrow #mu + #nu_{#mu}";
      if (chg == "Pl") { sampleLabel = "W^{#pm} #rightarrow #mu^{#pm} + #nu_{#mu}"; }
      textToPrint.push_back(sampleLabel);
      textToPrint.push_back("p^{#mu}_{T} > 25 GeV/c");
      //
      // Draw graph
      graph.Draw("COLZ");
      //
      // Format graph
      //
      // Set the Axis Titles
      std::string xLabel = "#eta^{#mu}";
      if (useEtaCM) { xLabel += "_{CM}"; }
      else { xLabel += "_{LAB}"; }
      std::string yLabel = formatResultVarName(var, useEtaCM, false, false, (chg=="Pl"?"Com":"Inc"));
      graph.SetTitle(Form("Correlation Matrix;%s;%s", "Analysis bin", "Analysis bin"));
      // X-axis
      graph.GetXaxis()->CenterTitle(kTRUE);
      graph.GetXaxis()->SetTitleOffset(0.6);
      graph.GetXaxis()->SetTitleSize(0.060);
      graph.GetXaxis()->SetLabelSize(0.0);
      // Y-axis
      graph.GetYaxis()->CenterTitle(kTRUE);
      graph.GetYaxis()->SetTitleOffset(0.7);
      graph.GetYaxis()->SetTitleSize(0.065);
      graph.GetYaxis()->SetLabelSize(0.0);
      // Z axis
      graph.GetZaxis()->SetLabelSize(0.030);
      graph.GetZaxis()->SetRangeUser(-1.0 , 1.0);
      c.Modified(); c.Update();
      // Draw the white box to hide under/overflow bins
      TBox box(c.GetFrame()->GetX1(), nCol, c.GetFrame()->GetX2(), c.GetFrame()->GetY2());
      box.SetFillColor(kWhite); box.SetLineColor(kWhite);
      box.Draw("same");
      redrawBorder();
      c.Modified(); c.Update();
      // Create line
      TLine line_1(nCol-dCol, 0      , nCol-dCol, nCol   ); line_1.SetLineWidth(3);
      TLine line_2(0      , nCol-dCol, nCol   , nCol-dCol); line_2.SetLineWidth(3);
      //
      if (nCol>dCol) {
        line_1.Draw("same");
        line_2.Draw("same");
      }
      //
      tex.SetTextSize(0.055); tex.DrawLatex(0.19, 0.85, textToPrint[0].c_str());
      tex.SetTextSize(0.058); tex.SetTextFont(61); tex.DrawLatex(0.75, 0.85, "CMS"); tex.SetTextFont(62);
      //tex.SetTextSize(0.044); tex.SetTextFont(52); tex.DrawLatex(0.66, 0.79, "Preliminary"); tex.SetTextFont(62);
      tex.SetTextSize(0.044); tex.SetTextFont(62); tex.DrawLatex(0.19, 0.79, yLabel.c_str()); tex.SetTextFont(62);
      if (textToPrint.size()>1) { tex.SetTextSize(0.035); tex.DrawLatex(0.19, 0.73, textToPrint[1].c_str()); }
      if (nCol>dCol) {
        tex.SetTextSize(0.035); tex.SetTextFont(62); tex.DrawLatex(0.42, 0.65, "Minus"); tex.SetTextFont(62);
        tex.SetTextSize(0.035); tex.SetTextFont(62); tex.DrawLatex(0.81, 0.65, "Plus" ); tex.SetTextFont(62);
        tex.SetTextSize(0.035); tex.SetTextFont(62); tex.SetTextAngle(90); tex.DrawLatex(0.21, 0.44, "Minus"); tex.SetTextFont(62); tex.SetTextAngle(0);
        tex.SetTextSize(0.035); tex.SetTextFont(62); tex.SetTextAngle(90); tex.DrawLatex(0.21, 0.15, "Plus" ); tex.SetTextFont(62); tex.SetTextAngle(0);
      }
      //
      // set the CMS style
      int option = 1118;
      CMS_lumi(&c, option, 33, "", false, 0.6, false);
      // Update
      c.Modified(); c.Update(); // Pure paranoia
      TPaletteAxis *palette = (TPaletteAxis*)graph.GetListOfFunctions()->FindObject("palette");
      palette->SetX2NDC(0.93);
      c.Modified();
      //
      // Save canvas
      //
      // Create Output Directory
      const std::string plotDir = outDir + "/Matrix/Correlation/" + col+"/Plot/"+var;
      makeDir(plotDir + "/png/");
      makeDir(plotDir + "/pdf/");
      makeDir(plotDir + "/root/");
      //
      // Save Canvas
      const std::string name = Form("corMatrix_WToMu%s_%s_%s_%s_%s", chg.c_str(), col.c_str(), var.c_str(), cat.first.c_str(), sys.first.c_str());
      c.SaveAs(( plotDir + "/png/"  + name + ".png"  ).c_str());
      c.SaveAs(( plotDir + "/pdf/"  + name + ".pdf"  ).c_str());
      c.SaveAs(( plotDir + "/root/" + name + ".root" ).c_str());
      //
      // Save Total
      if (cat.first=="Total" && sys.first=="Total") {
        const std::string dataDir = outDir + "/Matrix/Correlation/" + col+"/Data/"+var;
        makeDir(dataDir);
        const std::string fName = ( dataDir + "/"  + name + ".root"  );
        TFile f(fName.c_str(), "RECREATE"); f.cd();
        corMatrix.Write("corMatrix");
        nomVec.Write("nomVec");
        f.Write();
        f.Close();
      }
      //
      // Clean up memory
      c.Clear(); c.Close();
    }
  }
};


bool printCovarianceMatrix(const BinSeptaMapVec& systVarMap, const BinPentaMap& varMap, const WSDirMap& workDirNames, const std::string& outDir, const bool useEtaCM = true)
{
  //
  std::cout << "[INFO] Creating the covariance matrix for systematic uncertainties" << std::endl;
  //
  // Fill the tables
  //
  for (const auto& c : systVarMap.begin()->second.begin()->second[0]) {
    //
    // For Cross Section
    //
    // Create Output Directory
    const std::string tableDir_XSec = outDir + "/Matrix/Covariance/Table/" + c.first + "/Cross_Section";
    makeDir(tableDir_XSec);
    // Create Output Files for all the convariance matrices
    const std::string fileName_CovMatrix_XSec = Form("covMatrix_WToMu%s_%s_%s", "Pl", c.first.c_str(), "Cross_Section");
    std::ofstream file_CovMatrix_XSec((tableDir_XSec + "/" + fileName_CovMatrix_XSec + ".tex").c_str());
    if (file_CovMatrix_XSec.is_open()==false) { std::cout << "[ERROR] File " << fileName_CovMatrix_XSec << " was not created!" << std::endl; return false; }
    makeCovarianceMatrix(file_CovMatrix_XSec, workDirNames, systVarMap, varMap, c.first, "Cross_Section", "Pl", outDir, useEtaCM, "Total");
    //
    // For Forward Backward Ratio
    //
    // Create Output Directory
    const std::string tableDir_RFB = outDir + "/Matrix/Covariance/Table/" + c.first + "/ForwardBackward_Ratio";
    makeDir(tableDir_RFB);
    // Create Output Files for all the convariance matrices
    const std::string fileName_CovMatrix_RFB = Form("covMatrix_WToMu%s_%s_%s", "Pl", c.first.c_str(), "ForwardBackward_Ratio");
    std::ofstream file_CovMatrix_RFB((tableDir_RFB + "/" + fileName_CovMatrix_RFB + ".tex").c_str());
    if (file_CovMatrix_RFB.is_open()==false) { std::cout << "[ERROR] File " << fileName_CovMatrix_RFB << " was not created!" << std::endl; return false; }
    makeCovarianceMatrix(file_CovMatrix_RFB, workDirNames, systVarMap, varMap, c.first, "ForwardBackward_Ratio", "Pl", outDir, useEtaCM, "Total");
    //
    // Create Output Files for all the convariance matrices
    const std::string fileName_CovMatrix_RFB2 = Form("covMatrix_WToMu%s_%s_%s", "", c.first.c_str(), "ForwardBackward_Ratio");
    std::ofstream file_CovMatrix_RFB2((tableDir_RFB + "/" + fileName_CovMatrix_RFB2 + ".tex").c_str());
    if (file_CovMatrix_RFB2.is_open()==false) { std::cout << "[ERROR] File " << fileName_CovMatrix_RFB2 << " was not created!" << std::endl; return false; }
    makeCovarianceMatrix(file_CovMatrix_RFB2, workDirNames, systVarMap, varMap, c.first, "ForwardBackward_Ratio", "", outDir, useEtaCM, "Total");
    //
    // For Charge Asymmetry
    //
    // Create Output Directory
    const std::string tableDir_C = outDir + "/Matrix/Covariance/Table/" + c.first + "/Charge_Asymmetry";
    makeDir(tableDir_C);
    // Create Output Files for all the convariance matrices
    const std::string fileName_CovMatrix_C = Form("covMatrix_WToMu%s_%s_%s", "", c.first.c_str(), "Charge_Asymmetry");
    std::ofstream file_CovMatrix_C((tableDir_C + "/" + fileName_CovMatrix_C + ".tex").c_str());
    if (file_CovMatrix_C.is_open()==false) { std::cout << "[ERROR] File " << fileName_CovMatrix_C << " was not created!" << std::endl; return false; }
    makeCovarianceMatrix(file_CovMatrix_C, workDirNames, systVarMap, varMap, c.first, "Charge_Asymmetry", "", outDir, useEtaCM, "Total");
    //
  }
  //
  return true;
};


bool printCrossSectionTableForPaper(const BinPentaMap& varMap , const std::string& outDir)
{
  //
  // Create Output Directory
  const std::string tableDir = outDir + "/Tables/Paper/PA";
  makeDir(tableDir);
  // Create Output Files for Cross Sections
  const std::string fileName_XSEC = "crossSection_PA";
  std::ofstream file((tableDir + "/" + fileName_XSEC + ".tex").c_str());
  if (file.is_open()==false) { std::cout << "[ERROR] File " << fileName_XSEC << " was not created!" << std::endl; return false; }
  //
  // Initialize the latex table
  std::vector< std::string > texTable;
  //
  texTable.push_back("\\begin{table*}[htb!]");
  texTable.push_back(Form("  \\topcaption{\\label{tab:Cross-sectionMuPlus} Production cross section for $\\pPb \\to \\WToMuNu + X$ for positively (top) and negatively (bottom) charged muons of \\pt larger than 25\\GeVc, in nanobarns, as a function of the muon pseudorapidity in the center-of-mass frame. Quoted uncertainties are first statistical, then systematic. The global normalization uncertainty of %.1f\\%% is not included in the listed uncertainties.}", LUMIUNC_*100.));
  texTable.push_back("  \\centering");
  texTable.push_back("  \\resizebox{\\textwidth}{!}{");
  texTable.push_back("  \\begin{tabular}{c|cccccc}");
  texTable.push_back("    \\hline");
  //
  std::string tmp1 = "    $B\\times\\frac{\\rd\\sigma}{\\rd\\etaCM}$ (nb) [\\etaCM bin]";
  std::string tmp2 = "    $\\Pgmm$";
  std::string tmp3 = "    $\\Pgmp$";
  //
  int iCnt = 0;
  for (const auto& b : varMap.at("PA").at("Pl").at("Cross_Section").begin()->second) {
    // Add
    if (iCnt >= 6) {
      tmp1 += "\\\\"; tmp2 += "\\\\"; tmp3 += "\\\\";
      texTable.push_back(tmp1);
      texTable.push_back("    \\hline");
      texTable.push_back(tmp2);
      texTable.push_back(tmp3);
      texTable.push_back("    \\hline \\hline");
      //
      tmp1 = "    $B\\times\\frac{\\rd\\sigma}{\\rd\\etaCM}$ (nb) [\\etaCM bin]";
      tmp2 = "    $\\Pgmm$";
      tmp3 = "    $\\Pgmp$";
      //
      iCnt = 0;
    }
    // First Line
    const double min = b.first.etabin().low();
    const double max = b.first.etabin().high();
    tmp1 += Form(" & $[%g,%g]$", min, max);
    // Second Line
    const double statErrLow_Mi = varMap.at("PA").at("Mi").at("Cross_Section").at("Err_Stat_Low").at(b.first);
    const double statErrHi_Mi  = varMap.at("PA").at("Mi").at("Cross_Section").at("Err_Stat_High").at(b.first);
    const double systErrLow_Mi = varMap.at("PA").at("Mi").at("Cross_Section").at("Err_Syst_Low").at(b.first);
    const double systErrHi_Mi  = varMap.at("PA").at("Mi").at("Cross_Section").at("Err_Syst_High").at(b.first);
    const double rVal_Mi       = varMap.at("PA").at("Mi").at("Cross_Section").at("Val").at(b.first);
    tmp2 += Form(" & $%.1f\\pm%.1f\\pm%.1f$", rVal_Mi, statErrLow_Mi, systErrLow_Mi);
    // Third Line
    const double statErrLow_Pl = varMap.at("PA").at("Pl").at("Cross_Section").at("Err_Stat_Low").at(b.first);
    const double statErrHi_Pl  = varMap.at("PA").at("Pl").at("Cross_Section").at("Err_Stat_High").at(b.first);
    const double systErrLow_Pl = varMap.at("PA").at("Pl").at("Cross_Section").at("Err_Syst_Low").at(b.first);
    const double systErrHi_Pl  = varMap.at("PA").at("Pl").at("Cross_Section").at("Err_Syst_High").at(b.first);
    const double rVal_Pl       = varMap.at("PA").at("Pl").at("Cross_Section").at("Val").at(b.first);
    tmp3 += Form(" & $%.1f\\pm%.1f\\pm%.1f$", rVal_Pl, statErrLow_Pl, systErrLow_Pl);
    //
    iCnt++;
  }
  //
  tmp1 += "\\\\"; tmp2 += "\\\\"; tmp3 += "\\\\";
  texTable.push_back(tmp1);
  texTable.push_back("    \\hline");
  texTable.push_back(tmp2);
  texTable.push_back(tmp3);
  texTable.push_back("    \\hline");
  //
  texTable.push_back("  \\end{tabular}");
  texTable.push_back("  }");
  texTable.push_back("\\end{table*}");
  //
  for (const auto& row : texTable) { file << row << std::endl; }
  file << std::endl; file << std::endl;
  //
  return true;
};


#endif // ifndef resultUtils_h

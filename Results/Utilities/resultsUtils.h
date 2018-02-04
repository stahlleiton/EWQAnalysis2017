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
#include "TCanvas.h"
#include "TPad.h"
#include "TFrame.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TPaletteAxis.h"
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
#include "HIN13007.h"


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
//typedef std::vector< std::pair< std::string , uint > > CorrMap;
typedef std::map< std::string , std::pair< uint , uint > > CorrMap;


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
// Variation   Method: 4: Fully Uncorrelated , 1: Eta Correlated , 2: Charge Correlated , 3: Fully Correlated
CorrMap effTnPType_ = {
  { "TnP_Stat_Iso"    , { 100 , 1 } } , { "TnP_Stat_MuID"    , { 100 , 1 } } , { "TnP_Stat_Trig" , { 2 , 1 } } ,
  { "TnP_Syst_BinIso" , {   1 , 2 } } , { "TnP_Syst_BinMuID" , {   1 , 2 } } , { "TnP_Syst_Iso"  , { 2 , 2 } } , { "TnP_Syst_MuID" , { 2 , 2 } } , { "TnP_Syst_Trig" , { 2 , 2 } } ,
  { "TnP_Syst_PU"     , {   2 , 2 } } , { "TnP_Syst_STA"     , {   2 , 2 } } ,
  { "MC_Syst_PDF"     , {  96 , 3 } } , { "MC_Syst_Scale"    , {   6 , 3 } } , { "MC_Syst_Alpha" , { 2 , 3 } }
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
  if (type.find("TnP_")!=std::string::npos) {
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
    else if (variation.size() == 1) {
      Err_High = std::abs(variation[0] - nominal);
      Err_Low  = Err_High;
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
    const double convFactor = TMath::ErfcInverse((1.-0.68))/TMath::ErfcInverse((1.-0.90));
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
    // Convert from deltaAlpha_s = 0.002 to 0.0015 (68% CL)
    const double convFactor = (0.0015/0.0020);
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


bool computeChargeAsymmetry( BinPentaMap& var , const VarBinMap& inputVar , const bool& doSyst , const VarBinMap& nomVar , const uint& corrSyst )
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
          for (uint i=0; i<systT.size(); i++) {
            oVar.at(systT[i]).at(b.first) = sumErrors( oVar.at(systT[i]).at(b.first) , getChargeAsymmetryError(iMCVar_Pl.at("Val"), iMCVar_Mi.at("Val"), iMCVar_Pl.at(statT[i]), iMCVar_Mi.at(statT[i])) );
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
              for (const auto& t : systT) {
                // Case: Charge UnCorrelated
                if (effT.second.second==0 || effT.second.second==1) {
                  oVar.at(t).at(b.first) = sumErrors( oVar.at(t).at(b.first) , getChargeAsymmetryError(b.second.at("Pl").at(vL).at("Val"), b.second.at("Mi").at(vL).at("Val"),
                                                                                                       b.second.at("Pl").at(vL).at(t)    , b.second.at("Mi").at(vL).at(t)) );
                }
              }
            }
            // Case: Apply Variation Method
            if (effT.second.second>=2) {
              for (uint i = 0; i < effT.second.first; i++) {
                const std::string vL = std::string("N_WToMu_Efficiency_") + effT.first + ( (effT.second.first>1) ? Form("_%d", i) : "" );
                if (b.second.at("Mi").count(vL)==0 || b.second.at("Pl").count(vL)==0) {
                  std::cout << "[ERROR] " << vL << " variable is missing in bin [" << b.first.etabin().high() << " , " << b.first.etabin().low() << "]" << std::endl; return false;
                }
                // Case: Charge Correlated
                if (effT.second.second==2 || effT.second.second==3) {
                  variation[effT.first].push_back( getChargeAsymmetryValue(b.second.at("Pl").at(vL).at("Val"), b.second.at("Mi").at(vL).at("Val")) );
                }
                // Case: Charge Uncorrelated
                if (effT.second.second==4) {
                  variation[(effT.first+"_Pl")].push_back( getChargeAsymmetryValue(b.second.at("Pl").at(vL).at("Val") ,                  iVar_Mi.at("Val")) );
                  variation[(effT.first+"_Mi")].push_back( getChargeAsymmetryValue(                 iVar_Pl.at("Val") , b.second.at("Mi").at(vL).at("Val")) );
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
	for (const auto& t : systT) {
	  // Case: Charge Correlated
	  if (corrSyst==2 || corrSyst==3) {
	    for (const auto& t : systT) {
	      const double varVal = getChargeAsymmetryValue(iVar_Pl.at("Val"), iVar_Mi.at("Val"));
	      oVar.at(t).at(b.first) = (varVal - nomVal);
              oVar.at("Val").at(b.first) = varVal;
	    }
	  }
	  // Case: Charge Uncorrelated
	  if (corrSyst==0 || corrSyst==1 ||  corrSyst==4) {
	    var.at(c.first).at("")["Charge_Asymmetry_Pl"][t][b.first] = std::abs(getChargeAsymmetryValue(iVar_Pl.at("Val"), nVar_Mi.at("Val")) - nomVal);
	    var.at(c.first).at("")["Charge_Asymmetry_Mi"][t][b.first] = std::abs(getChargeAsymmetryValue(nVar_Pl.at("Val"), iVar_Mi.at("Val")) - nomVal);
	  }
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


bool computeForwardBackwardRatio( BinPentaMap& var , const VarBinMap& inputVar , const bool& doSyst , const VarBinMap& nomVar , const uint& corrSyst )
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
          for (uint i=0; i<systT.size(); i++) {
            oVar_Inc.at(systT[i]).at(b.first) = sumErrors( oVar_Inc.at(systT[i]).at(b.first) ,
                                                           getForwardBackwardRatioError(
                                                                                        iMCVar_FwPl.at("Val")    , iMCVar_BwPl.at("Val")    , iMCVar_FwMi.at("Val")    , iMCVar_BwMi.at("Val"), 
                                                                                        iMCVar_FwPl.at(statT[i]) , iMCVar_BwPl.at(statT[i]) , iMCVar_FwMi.at(statT[i]) , iMCVar_BwMi.at(statT[i])
                                                                                        ) );
            oVar_Pl.at(systT[i]).at(b.first) =  sumErrors( oVar_Pl.at(systT[i]).at(b.first) ,
                                                           getForwardBackwardRatioError(iMCVar_FwPl.at("Val"), iMCVar_BwPl.at("Val"), iMCVar_FwPl.at(statT[i]), iMCVar_BwPl.at(statT[i])) );
            oVar_Mi.at(systT[i]).at(b.first) =  sumErrors( oVar_Mi.at(systT[i]).at(b.first) ,
                                                           getForwardBackwardRatioError(iMCVar_FwMi.at("Val"), iMCVar_BwMi.at("Val"), iMCVar_FwMi.at(statT[i]), iMCVar_BwMi.at(statT[i])) );
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
              //
              for (const auto& t : systT) {
                // Case: Eta UnCorrelated
                if (effT.second.second==0 || effT.second.second==2) {
                  oVar_Pl.at(t).at(b.first) = sumErrors( oVar_Pl.at(t).at(b.first) ,
                                                         getForwardBackwardRatioError(iEVar_FwPl.at("Val"), iEVar_BwPl.at("Val"), iEVar_FwPl.at(t), iEVar_BwPl.at(t)) );
                  oVar_Mi.at(t).at(b.first) = sumErrors( oVar_Mi.at(t).at(b.first) ,
                                                         getForwardBackwardRatioError(iEVar_FwMi.at("Val"), iEVar_BwMi.at("Val"), iEVar_FwMi.at(t), iEVar_BwMi.at(t)) );
                }
              }
            }
            // Case: Apply Variation Method
            if (effT.second.second==1 || effT.second.second==3 || effT.second.second==4) {
              for (uint i = 0; i < effT.second.first; i++) {
                const std::string vL = std::string("N_WToMu_Efficiency_") + effT.first + ( (effT.second.first>1) ? Form("_%d", i) : "" );
                if (c.second.at(binFw).at("Pl").count(vL)==0 || c.second.at(binBw).at("Pl").count(vL)==0 || c.second.at(binFw).at("Mi").count(vL)==0 || c.second.at(binBw).at("Mi").count(vL)==0) {
                  std::cout << "[ERROR] " << vL << " variable is missing in bin [" << b.first.etabin().high() << " , " << b.first.etabin().low() << "]" << std::endl; return false;
                }
                const auto& iEVar_FwPl = c.second.at(binFw).at("Pl").at(vL);
                const auto& iEVar_FwMi = c.second.at(binFw).at("Mi").at(vL);
                const auto& iEVar_BwPl = c.second.at(binBw).at("Pl").at(vL);
                const auto& iEVar_BwMi = c.second.at(binBw).at("Mi").at(vL);
                //
                // Case: Eta Correlated
                if (effT.second.second==1 || effT.second.second==3) {
                  variation_Pl [effT.first].push_back( getForwardBackwardRatioValue( iEVar_FwPl.at("Val") , iEVar_BwPl.at("Val") ) );
                  variation_Mi [effT.first].push_back( getForwardBackwardRatioValue( iEVar_FwMi.at("Val") , iEVar_BwMi.at("Val") ) );
                }
                // Case: Eta UnCorrelated
                else if (effT.second.second==4) {
                  variation_Pl [(effT.first+"_Fw")].push_back( getForwardBackwardRatioValue( iEVar_FwPl.at("Val") , iVar_BwPl.at("Val") ) );
                  variation_Mi [(effT.first+"_Fw")].push_back( getForwardBackwardRatioValue( iEVar_FwMi.at("Val") , iVar_BwMi.at("Val") ) );
                  variation_Pl [(effT.first+"_Bw")].push_back( getForwardBackwardRatioValue( iVar_FwPl.at("Val") , iEVar_BwPl.at("Val") ) );
                  variation_Mi [(effT.first+"_Bw")].push_back( getForwardBackwardRatioValue( iVar_FwMi.at("Val") , iEVar_BwMi.at("Val") ) );
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
              //
              for (const auto& t : systT) {
                // Case: Eta and Charge UnCorrelated
                if (effT.second.second==0) {
                  oVar_Inc.at(t).at(b.first) = sumErrors( oVar_Inc.at(t).at(b.first) ,
                                                          getForwardBackwardRatioError(
                                                                                       iEVar_FwPl.at("Val") , iEVar_BwPl.at("Val") , iEVar_FwMi.at("Val") , iEVar_BwMi.at("Val"), 
                                                                                       iEVar_FwPl.at(t)     , iEVar_BwPl.at(t)     , iEVar_FwMi.at(t)     , iEVar_BwMi.at(t)
                                                                                       ) );
                }
              }
            }
            // Case: Apply Variation Method
            if (effT.second.second>0) {
              for (uint i = 0; i < effT.second.first; i++) {
                const std::string vL = std::string("N_WToMu_Efficiency_") + effT.first + ( (effT.second.first>1) ? Form("_%d", i) : "" );
                if (c.second.at(binFw).at("Pl").count(vL)==0 || c.second.at(binBw).at("Pl").count(vL)==0 || c.second.at(binFw).at("Mi").count(vL)==0 || c.second.at(binBw).at("Mi").count(vL)==0) {
                  std::cout << "[ERROR] " << vL << " variable is missing in bin [" << b.first.etabin().high() << " , " << b.first.etabin().low() << "]" << std::endl; return false;
                }
                const auto& iEVar_FwPl = c.second.at(binFw).at("Pl").at(vL);
                const auto& iEVar_FwMi = c.second.at(binFw).at("Mi").at(vL);
                const auto& iEVar_BwPl = c.second.at(binBw).at("Pl").at(vL);
                const auto& iEVar_BwMi = c.second.at(binBw).at("Mi").at(vL);
                //
                // Case: Fully Correlated
                if (effT.second.second==3) {
                  variation_Inc[effT.first].push_back( getForwardBackwardRatioValue( ( iEVar_FwPl.at("Val") + iEVar_FwMi.at("Val") ) , ( iEVar_BwPl.at("Val") + iEVar_BwMi.at("Val") ) ) );
                }
                // Case: Charge Correlated and Eta UnCorrelated
                if (effT.second.second==2) {
                  variation_Inc[(effT.first+"_Fw")].push_back( getForwardBackwardRatioValue( ( iEVar_FwPl.at("Val") + iEVar_FwMi.at("Val") ) , (  iVar_BwPl.at("Val") +  iVar_BwMi.at("Val") ) ) );
                  variation_Inc[(effT.first+"_Bw")].push_back( getForwardBackwardRatioValue( (  iVar_FwPl.at("Val") +  iVar_FwMi.at("Val") ) , ( iEVar_BwPl.at("Val") + iEVar_BwMi.at("Val") ) ) );
                }
                // Case: Charge UnCorrelated and Eta Correlated
                if (effT.second.second==1) {
                  variation_Inc[(effT.first+"_Pl")].push_back( getForwardBackwardRatioValue( ( iEVar_FwPl.at("Val") +  iVar_FwMi.at("Val") ) , ( iEVar_BwPl.at("Val") +  iVar_BwMi.at("Val") ) ) );
                  variation_Inc[(effT.first+"_Mi")].push_back( getForwardBackwardRatioValue( (  iVar_FwPl.at("Val") + iEVar_FwMi.at("Val") ) , (  iVar_BwPl.at("Val") + iEVar_BwMi.at("Val") ) ) );
                }
                // Case: Fully UnCorrelated
                if (effT.second.second==4) {
                  variation_Inc[(effT.first+"_Fw_Pl")].push_back( getForwardBackwardRatioValue( ( iEVar_FwPl.at("Val") +  iVar_FwMi.at("Val") ) , (  iVar_BwPl.at("Val") +  iVar_BwMi.at("Val") ) ) );
                  variation_Inc[(effT.first+"_Bw_Pl")].push_back( getForwardBackwardRatioValue( (  iVar_FwPl.at("Val") +  iVar_FwMi.at("Val") ) , ( iEVar_BwPl.at("Val") +  iVar_BwMi.at("Val") ) ) );
                  variation_Inc[(effT.first+"_Fw_Mi")].push_back( getForwardBackwardRatioValue( (  iVar_FwPl.at("Val") + iEVar_FwMi.at("Val") ) , (  iVar_BwPl.at("Val") +  iVar_BwMi.at("Val") ) ) );
                  variation_Inc[(effT.first+"_Bw_Mi")].push_back( getForwardBackwardRatioValue( (  iVar_FwPl.at("Val") +  iVar_FwMi.at("Val") ) , (  iVar_BwPl.at("Val") + iEVar_BwMi.at("Val") ) ) );
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
	for (const auto& t : systT) {
	  // Case: Eta Correlated
	  if (corrSyst==1 || corrSyst==3) {
	    const double varVal_Pl = getForwardBackwardRatioValue(iVar_FwPl.at("Val") , iVar_BwPl.at("Val"));
	    const double varVal_Mi = getForwardBackwardRatioValue(iVar_FwMi.at("Val") , iVar_BwMi.at("Val"));
	    oVar_Pl.at(t).at(b.first) = (varVal_Pl - nomVal_Pl);
	    oVar_Mi.at(t).at(b.first) = (varVal_Mi - nomVal_Mi);
            oVar_Pl.at("Val").at(b.first) = varVal_Pl;
            oVar_Mi.at("Val").at(b.first) = varVal_Mi;
	  }
	  // Case: Eta UnCorrelated
	  else if (corrSyst==0 || corrSyst==2 || corrSyst==4) {
	    var.at(c.first).at("Pl")["ForwardBackward_Ratio_Fw"][t][b.first] = std::abs(getForwardBackwardRatioValue(iVar_FwPl.at("Val"), nVar_BwPl.at("Val")) - nomVal_Pl);
	    var.at(c.first).at("Pl")["ForwardBackward_Ratio_Bw"][t][b.first] = std::abs(getForwardBackwardRatioValue(nVar_FwPl.at("Val"), iVar_BwPl.at("Val")) - nomVal_Pl);
	    var.at(c.first).at("Mi")["ForwardBackward_Ratio_Fw"][t][b.first] = std::abs(getForwardBackwardRatioValue(iVar_FwMi.at("Val"), nVar_BwMi.at("Val")) - nomVal_Mi);
	    var.at(c.first).at("Mi")["ForwardBackward_Ratio_Bw"][t][b.first] = std::abs(getForwardBackwardRatioValue(nVar_FwMi.at("Val"), iVar_BwMi.at("Val")) - nomVal_Mi);
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
	for (const auto& t : systT) {
	  // Case: Fully Correlated
	  if (corrSyst==3) {
            const double varVal_Inc = getForwardBackwardRatioValue((iVar_FwPl.at("Val")+iVar_FwMi.at("Val")) , (iVar_BwPl.at("Val")+iVar_BwMi.at("Val")));
	    oVar_Inc.at(t).at(b.first) = (varVal_Inc - nomVal_Inc);
            oVar_Inc.at("Val").at(b.first) = varVal_Inc;
	  }
	  // Case: Charge Correlated and Eta UnCorrelated
	  if (corrSyst==2) {
	    var.at(c.first).at("")["ForwardBackward_Ratio_Fw"][t][b.first] =
	      std::abs(getForwardBackwardRatioValue((iVar_FwPl.at("Val")+iVar_FwMi.at("Val")), (nVar_BwPl.at("Val")+nVar_BwMi.at("Val"))) - nomVal_Inc);
	    var.at(c.first).at("")["ForwardBackward_Ratio_Bw"][t][b.first] =
	      std::abs(getForwardBackwardRatioValue((nVar_FwPl.at("Val")+nVar_FwMi.at("Val")), (iVar_BwPl.at("Val")+iVar_BwMi.at("Val"))) - nomVal_Inc);
	  }
	  // Case: Charge UnCorrelated and Eta Correlated
	  if (corrSyst==1) {
	    var.at(c.first).at("")["ForwardBackward_Ratio_Pl"][t][b.first] =
	      std::abs(getForwardBackwardRatioValue((iVar_FwPl.at("Val")+nVar_FwMi.at("Val")), (iVar_BwPl.at("Val")+nVar_BwMi.at("Val"))) - nomVal_Inc);
	    var.at(c.first).at("")["ForwardBackward_Ratio_Mi"][t][b.first] =
	      std::abs(getForwardBackwardRatioValue((nVar_FwPl.at("Val")+iVar_FwMi.at("Val")), (nVar_BwPl.at("Val")+iVar_BwMi.at("Val"))) - nomVal_Inc);
	  }
	  // Case: Fully UnCorrelated
	  if (corrSyst==0 || corrSyst==4) {
	    var.at(c.first).at("")["ForwardBackward_Ratio_Fw_Pl"][t][b.first] =
	      std::abs(getForwardBackwardRatioValue((iVar_FwPl.at("Val")+nVar_FwMi.at("Val")), (nVar_BwPl.at("Val")+nVar_BwMi.at("Val"))) - nomVal_Inc);
	    var.at(c.first).at("")["ForwardBackward_Ratio_Bw_Pl"][t][b.first] =
		std::abs(getForwardBackwardRatioValue((nVar_FwPl.at("Val")+nVar_FwMi.at("Val")), (iVar_BwPl.at("Val")+nVar_BwMi.at("Val"))) - nomVal_Inc);
	    var.at(c.first).at("")["ForwardBackward_Ratio_Fw_Mi"][t][b.first] =
	      std::abs(getForwardBackwardRatioValue((nVar_FwPl.at("Val")+iVar_FwMi.at("Val")), (nVar_BwPl.at("Val")+nVar_BwMi.at("Val"))) - nomVal_Inc);
	    var.at(c.first).at("")["ForwardBackward_Ratio_Bw_Mi"][t][b.first] =
	      std::abs(getForwardBackwardRatioValue((nVar_FwPl.at("Val")+nVar_FwMi.at("Val")), (nVar_BwPl.at("Val")+iVar_BwMi.at("Val"))) - nomVal_Inc);
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


bool computeCrossSection( BinPentaMap& var , const VarBinMap& inputVar , const bool& doSyst , const VarBinMap& nomVar , const uint& corrSyst )
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
          for (uint i=0; i<systT.size(); i++) {
            oVar_Pl.at(systT[i]).at(b.first) = sumErrors( oVar_Pl.at(systT[i]).at(b.first) , getCrossSectionError(iMCVar_Pl.at("Val"), Luminosity, BinWidth, iMCVar_Pl.at(statT[i]), Err_Luminosity.at(systT[i])) );
          }
          for (uint i=0; i<systT.size(); i++) {
            oVar_Mi.at(systT[i]).at(b.first) = sumErrors( oVar_Mi.at(systT[i]).at(b.first) , getCrossSectionError(iMCVar_Mi.at("Val"), Luminosity, BinWidth, iMCVar_Mi.at(statT[i]), Err_Luminosity.at(systT[i])) );
          }
          // Add the other Systematic Errors of the Efficiency
          DoubleVecMap variation_Pl , variation_Mi;
          for (const auto& effT : effTnPType_) {
            //
            // Case: Apply Propagation Method
            if (effT.second.second<4) {
              const std::string vL = Form("N_WToMu_Efficiency_%s_PROP", effT.first.c_str());
              if (b.second.at("Mi").count(vL)==0 || b.second.at("Pl").count(vL)==0) {
                std::cout << "[ERROR] " << vL << " variable is missing in bin [" << b.first.etabin().high() << " , " << b.first.etabin().low() << "]" << std::endl; return false;
                }
              for (const auto& t : systT) {
                oVar_Pl.at(t).at(b.first) = sumErrors( oVar_Pl.at(t).at(b.first) , getCrossSectionError(b.second.at("Pl").at(vL).at("Val"), Luminosity, BinWidth, b.second.at("Pl").at(vL).at(t), 0.0) );
              }
              for (const auto& t : systT) {
                oVar_Mi.at(t).at(b.first) = sumErrors( oVar_Mi.at(t).at(b.first) , getCrossSectionError(b.second.at("Mi").at(vL).at("Val"), Luminosity, BinWidth, b.second.at("Mi").at(vL).at(t), 0.0) );
              }
            }
            // Case: Apply Variation Method
            if (effT.second.second>=4) {
              for (uint i = 0; i < effT.second.first; i++) {
                const std::string vL = std::string("N_WToMu_Efficiency_") + effT.first + ( (effT.second.first>1) ? Form("_%d", i) : "" );
                if (b.second.at("Mi").count(vL)==0 || b.second.at("Pl").count(vL)==0) {
                  std::cout << "[ERROR] " << vL << " variable is missing in bin [" << b.first.etabin().high() << " , " << b.first.etabin().low() << "]" << std::endl; return false;
                }
                variation_Pl[effT.first].push_back( getCrossSectionValue(b.second.at("Pl").at(vL).at("Val"), Luminosity, BinWidth) );
                variation_Mi[effT.first].push_back( getCrossSectionValue(b.second.at("Mi").at(vL).at("Val"), Luminosity, BinWidth) );
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
	  oVar_Pl.at(t).at(b.first) = (varVal_Pl - nomVal_Pl);
	  oVar_Mi.at(t).at(b.first) = (varVal_Mi - nomVal_Mi);
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
    for (const auto& s : sys) { mean += s.second; }
    mean /= sys.size();
    for (auto& s : sys) { s.second = mean; }
  }
  if (type=="Max") {
    double max = -99999999.0;
    for (const auto& s : sys) { if (max < s.second) { max = s.second; } }
    for (auto& s : sys) { s.second = max; }
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
		  }
		}
              }
	      valMap[vLbl]["Nom"] = nomVar.at(c.first).at(chg.first).at(vLbl).at("Val");
              for (const auto& b : valMap.at(vLbl).at("Err_Syst_High")) {
                if (valMap.at(vLbl).at("Val").at(b.first)==valMap.at(vLbl).at("Nom").at(b.first)) {
                  valMap.at(vLbl).at("Val").at(b.first) += valMap.at(vLbl).at("Err_Syst_High").at(b.first);
                }
              }
	    }
	  }
	}
	systVarVec.at(lbl.first).push_back(systVar);
      }
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
		  const uint nVariation = lbl.second.size();
		  double uncVal = 0.0;
		  if (nVariation > 2) {
		    double sum = 0.0;
		    for (uint i = 0; i < nVariation; i++) {
		      const double diff = std::abs(lbl.second[i].at(c.first).at(chg.first).at(v.first).at(t.first).at(b.first));
		      sum += ( diff * diff );
		    }
		    uncVal = std::sqrt( sum / nVariation );
		  }
		  else if (nVariation == 2) {
		    const double diff_0 = std::abs(lbl.second[0].at(c.first).at(chg.first).at(v.first).at(t.first).at(b.first));
		    const double diff_1 = std::abs(lbl.second[1].at(c.first).at(chg.first).at(v.first).at(t.first).at(b.first));
		    uncVal = std::max( diff_0 , diff_1 );
		  }
		  else if (nVariation == 1) {
		    const double diff = std::abs(lbl.second[0].at(c.first).at(chg.first).at(v.first).at(t.first).at(b.first));
		    uncVal = diff;
		  }
		  if (vLbl==v.first) { val = uncVal; }
		  else { val = sumErrors(val, uncVal); }
		}
              }
	    }
            valMap[vLbl]["Nom"] = nomVar.at(c.first).at(chg.first).at(vLbl).at("Val");
            for (const auto& b : valMap.at(vLbl).at("Err_Syst_High")) {
              if (valMap.at(vLbl).at("Val").at(b.first)==valMap.at(vLbl).at("Nom").at(b.first)) {
                valMap.at(vLbl).at("Val").at(b.first) += valMap.at(vLbl).at("Err_Syst_High").at(b.first);
              }
            }
          }
          // Apply the mean to QCD Mean
          if (lbl.first.find("_Avg")!=std::string::npos) {
            for (auto& sV : valMap) { for (auto& sT : sV.second) { if (sT.first.find("Err_Syst_")!=std::string::npos) { procSyst(valMap.at(sV.first).at(sT.first), "Mean"); } } }
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
                  if (origVarVec.size()>1) { val = nomVal + var.at(c.first).at(chg.first).at(sV.first).at("Err_Syst_High").at(sB.first); }
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


std::string formatResultVarName(const std::string varName, const bool useEtaCM, const bool isSyst = false, const bool useLATEX = false)
{
  std::string label = "";
  const std::string etaLbl = std::string(useLATEX ? "\\" : "#") + ( useEtaCM ? "eta_{CM}" : "eta_{LAB}" );
  if (isSyst) {
    if (varName == "Charge_Asymmetry"      ) { label = "Abs. Unc. ( N^{+} - N^{-} ) / ( N^{+} + N^{-} )"; }
    if (varName == "ForwardBackward_Ratio" ) { label = "Abs. Unc. R_{FB}"; }
    if (varName == "Cross_Section"         ) { label = Form("Rel. Unc. B %s/d%s", (useLATEX ? "\\times d\\sigma" : "#times d#sigma"), etaLbl.c_str()); }
  }
  else {
    if (varName == "Charge_Asymmetry"      ) { label = "( N^{+} - N^{-} ) / ( N^{+} + N^{-} )"; }
    if (varName == "ForwardBackward_Ratio" ) { label = "R_{FB}"; }
    if (varName == "Cross_Section"         ) { label = Form("B %s/d%s (nb)", (useLATEX ? "\\times d\\sigma" : "#times d#sigma"), etaLbl.c_str()); }
    if (varName == "N_WToMu"               ) { label = "Signal Yield"; }
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
          const uint nGraph = ( (lbl.first=="Nominal") ? nomGraphType.size() : lbl.second.size() );
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
          for (uint iGr=0; iGr<lbl.second.size(); iGr++) {
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
  //
};


void formatLegendEntry(TLegendEntry& e)
{
  e.SetTextSize(0.028);
};


void formatResultsGraph(TGraphAsymmErrors& graph, const std::string& col, const std::string& var, const std::string& chg, const bool& useEtaCM, const bool& incAcc, const bool isSyst = false)
{
  //
  // Set the Axis Titles
  std::string xLabel = "#mu"; if (chg == "Pl") { xLabel += "^{+}"; }; if (chg == "Mi") { xLabel += "^{-}"; }; xLabel += " #eta";
  if (useEtaCM) { xLabel += "_{CM}"; }
  else { xLabel += "_{LAB}"; }
  std::string yLabel = formatResultVarName(var, useEtaCM, isSyst);
  graph.SetTitle(Form(";%s;%s", xLabel.c_str(), yLabel.c_str()));
  //
  // General
  graph.SetMarkerColor(kBlue);
  graph.SetMarkerStyle(20);
  graph.SetMarkerSize(1.0);
  graph.SetFillStyle(1001);
  // X-axis
  graph.GetXaxis()->CenterTitle(kFALSE);
  graph.GetXaxis()->SetTitleOffset(0.9);
  graph.GetXaxis()->SetTitleSize(0.050);
  graph.GetXaxis()->SetLabelSize(0.035);
  double xMin=-2.5, xMax=2.5;
  if (useEtaCM) {
    if ( var == "Charge_Asymmetry"      ) { xMin = -3.0; xMax = 2.1; }
    if ( var == "Cross_Section"         ) { xMin = -3.0; xMax = 2.1; }
    if ( var == "ForwardBackward_Ratio" ) { xMin = -0.1; xMax = 2.1; }
  }
  graph.GetXaxis()->SetLimits(xMin , xMax);
  // Y-axis
  graph.GetYaxis()->CenterTitle(kFALSE);
  graph.GetYaxis()->SetTitleOffset(1.38);
  graph.GetYaxis()->SetTitleSize(0.050);
  if ( var == "Charge_Asymmetry" ) { graph.GetYaxis()->SetTitleSize(0.040); }
  graph.GetYaxis()->SetLabelSize(0.035);
  if ( var == "Charge_Asymmetry"      ) { graph.GetYaxis()->SetRangeUser(-0.2, 0.4); }
  if ( var == "ForwardBackward_Ratio" ) { graph.GetYaxis()->SetRangeUser( 0.6, 1.5); }
  if (incAcc){ if ( var == "Cross_Section" ) { graph.GetYaxis()->SetRangeUser(0.0, 300.0); } }
  else       { if ( var == "Cross_Section" ) { graph.GetYaxis()->SetRangeUser(0.0, 250.0); } }
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
          if (chg == "Mi") { sampleLabel = "W^{-} #rightarrow #mu^{-} + #nu_{#mu}"; }
          textToPrint.push_back(sampleLabel);
          if (accType=="") { textToPrint.push_back("p^{#mu}_{T} > 25 GeV/c"); }
          //
          // Declare the graph vector (for drawing with markers)
          std::vector< TGraphAsymmErrors > grVec;
          // Initialize the Legend
          double legOff = 0.0; if (accType=="") { legOff = 0.05; }
          TLegend leg(0.2, (0.71 - legOff), 0.4, (0.84 - legOff));
          // Initialize the graph x range variables
          double xMin=0.0 , xMax=0.0 , yErrMin=9999999. , yErrMax=-1. , fbYLine = -1.;
          // Draw graph
          if ( type=="Nominal" ) {
            // Format the graphs
            const bool incAcc = (accType!="");
            grVec.push_back(graph.at("Err_Tot"));
            grVec.push_back(graph.at("Err_Stat"));
            for (int j = 0; j < grVec.back().GetN(); j++) {
              grVec.back().SetPointError(j, grVec.back().GetErrorXlow(j)*0.8, grVec.back().GetErrorXhigh(j)*0.8, grVec.back().GetErrorYlow(j), grVec.back().GetErrorYhigh(j));
            }
            if (graph.count("Err_Syst")>0) {
	      grVec.push_back(graph.at("Err_Syst"));
	      for (int j = 0; j < grVec.back().GetN(); j++) {
		grVec.back().SetPointError(j, grVec.back().GetErrorXlow(j)*0.4, grVec.back().GetErrorXhigh(j)*0.4, grVec.back().GetErrorYlow(j), grVec.back().GetErrorYhigh(j));
	      }
	    }
            for (auto& gr : grVec) { formatResultsGraph(gr, col, var, chg, useEtaCM, incAcc); }
            grVec[0].SetMarkerColor(kBlack);
            grVec[1].SetFillColor(kGreen+3);
	    if (graph.count("Err_Syst")>0) { grVec[2].SetFillColor(kOrange); }
            // Create Legend
            formatLegendEntry(*leg.AddEntry(&grVec[0], "Data", "pe"));
            formatLegendEntry(*leg.AddEntry(&grVec[1], "Statistical Uncertainty", "f"));
	    if (graph.count("Err_Syst")>0) { formatLegendEntry(*leg.AddEntry(&grVec[2], "Systematic Uncertainty", "f")); }
            // Draw the graphs
            grVec[0].Draw("ap");
            grVec[1].Draw("same2");
            if (graph.count("Err_Syst")>0) { grVec[2].Draw("same2"); }
            grVec[0].Draw("samep");
            //
            xMin = grVec[0].GetXaxis()->GetXmin(); xMax = grVec[0].GetXaxis()->GetXmax(); fbYLine = 1.0;
          }
          else {
            const auto ref = graph.at("Total");
            // Extract the Varied Efficiency graphs
            uint iMax = 0; double yMax = -9999999999.0;
            grVec.push_back(graph.at("Total"));
            for (auto gr : graph) { if (gr.first!="Total") { grVec.push_back(gr.second); } }
            for (uint i=0; i<grVec.size(); i++) {
              if (i==0) {
                for (int j = 0; j < grVec[i].GetN(); j++) {
                  double x, y; ref.GetPoint(j, x, y); grVec[i].SetPoint(j, x, 0.0);
                  double errLo = grVec[i].GetErrorYlow(j) , errHi = grVec[i].GetErrorYhigh(j);
                  if (var=="Cross_Section") { errLo /= y; errHi /= y; }
                  grVec[i].SetPointError(j, grVec[i].GetErrorXlow(j), grVec[i].GetErrorXhigh(j), errLo, errHi);
                  if (yMax < 0.90*std::max(errLo, errHi)) { yMax = std::max(errLo, errHi); }
                }
              }
              else {
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
            grVec[0].SetMarkerColor(kBlack);
            grVec[0].SetMarkerSize(1.0);
            grVec[0].SetLineColor(kBlack);
            grVec[0].SetLineWidth(3);
            grVec[0].SetFillStyle(0);
            grVec[0].SetMarkerSize(0.0);
            setRangeYAxisGraph(grVec[iMax], 0.1, 0.4);
            for (uint i=0; i<grVec.size(); i++) { if (i!=0) { grVec[i].SetMarkerColor(kBlack); grVec[i].SetLineWidth(3); grVec[i].SetMarkerSize(0.0); grVec[i].SetLineColor(kRed); } }
            TGaxis::SetMaxDigits(3); // to display powers of 10
            // Create Legend
            formatLegendEntry(*leg.AddEntry(&grVec[0], Form("%s Uncertainty", type.c_str()), "f"));
            formatLegendEntry(*leg.AddEntry(&grVec[1], Form("%s Variation", type.c_str()), "l"));
            // Draw the Graph
            grVec[iMax].Draw("apx");
            grVec[0].Draw("same2");
            for (uint i=0; i<grVec.size(); i++) { if (i!=0) { grVec[i].Draw("samep"); } }
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
          for (const auto& s: textToPrint) { tex.DrawLatex(0.22, 0.86-dy, s.c_str()); dy+=0.045; }
          if ( type!="Nominal" )  {
            tex.SetTextSize(0.030);
            if (var=="Cross_Section") { tex.DrawLatex(0.20, 0.16, Form("Rel. Unc. [ %.2f%% , %.2f%% ]", (yErrMin*100.), (yErrMax*100.))); }
            else { tex.DrawLatex(0.20, 0.16, Form("Abs. Unc. [ %.2f%% , %.2f%% ]", (yErrMin*100.), (yErrMax*100.))); }
          }
          else {
            if (var=="Cross_Section") { tex.SetTextSize(0.025); tex.DrawLatex(0.58, 0.72, "Luminosity uncertainty: 5.0%"); }
          }
          // Update
          c.Modified(); c.Update(); // Pure paranoia
          //
          // set the CMS style
          int option = 118;
          if (col.find("pPb")!=std::string::npos) option = 115;
          if (col.find("Pbp")!=std::string::npos) option = 116;
          CMS_lumi(&c, option, 33, "");
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
          const std::string plotDir = outDir+"/Plots/Result/" + col+"/" + label+"/" + var;
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
          if (chg == "Mi") { sampleLabel = "W^{-} #rightarrow #mu^{-} + #nu_{#mu}"; }
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
          graph.SetMarkerColor(kBlack);
          graph.SetMarkerSize(1.0);
          graph.SetLineColor(kBlack);
          graph.SetLineWidth(3);
          graph.SetFillStyle(0);
          setRangeYAxisGraph(graph, 0.1, 0.4);
          TGaxis::SetMaxDigits(3); // to display powers of 10
          // Create Legend
          formatLegendEntry(*leg.AddEntry(&graph, type.c_str(), "f"));
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
          for (const auto& s: textToPrint) { tex.DrawLatex(0.22, 0.86-dy, s.c_str()); dy+=0.045; }
          // Update
          c.Modified(); c.Update(); // Pure paranoia
          //
          // set the CMS style
          int option = 118;
          if (col.find("pPb")!=std::string::npos) option = 115;
          if (col.find("Pbp")!=std::string::npos) option = 116;
          CMS_lumi(&c, option, 33, "");
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
          const std::string plotDir = outDir+"/Plots/Systematic/" + col+"/" + label+"/" + var;
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


void drawCombineSystematicGraph( const GraphPentaMap& graphMap , const std::string& outDir , const bool useEtaCM = true , const std::string accType = "MC" , const std::string effType = "TnP" )
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
        if (chg == "Mi") { sampleLabel = "W^{-} #rightarrow #mu^{-} + #nu_{#mu}"; }
        textToPrint.push_back(sampleLabel);
        if (accType=="") { textToPrint.push_back("p^{#mu}_{T} > 25 GeV/c"); }
        //
        // Declare the graph map (for drawing with markers)
        std::map< std::string , TGraphAsymmErrors > grMap;
        // Initialize the Legend
        double legOff = 0.0; if (accType=="") { legOff = 0.05; }
        TLegend leg(0.2, (0.54 - legOff), 0.4, (0.82 - legOff));
        double xMin = 0.0 , xMax = 0.0;
        //
        for (const auto& lbl : v.second) {
          if (lbl.first=="Nominal") continue;
          const std::string type = lbl.first;
          auto& graph = lbl.second;
          grMap[lbl.first] = graph.at("Total");
        }
        grMap["Systematic" ] = v.second.at("Nominal").at("Err_Syst");
        grMap["Statistical"] = v.second.at("Nominal").at("Err_Stat");
        //
        std::string lMax = ""; double vMax = -99999999.0;
        for (auto& lbl : grMap) {
          for (int j = 0; j < grMap.at(lbl.first).GetN(); j++) {
            double X , Y; grMap.at(lbl.first).GetPoint(j, X, Y); Y = std::abs(Y);
            grMap.at(lbl.first).SetPoint(j, X, 0);
            //
            double Err_Y_Low  = grMap.at(lbl.first).GetErrorYlow(j);
            double Err_Y_High = grMap.at(lbl.first).GetErrorYhigh(j);
            if (var=="Cross_Section") { Err_Y_Low /= std::abs(Y); Err_Y_High /= std::abs(Y); }
            //
            grMap.at(lbl.first).SetPointError(j, 
                                              grMap.at(lbl.first).GetErrorXlow(j), grMap.at(lbl.first).GetErrorXhigh(j), 
                                              Err_Y_Low, Err_Y_High
                                              );
            if (vMax < std::max(Err_Y_Low, Err_Y_High)) { vMax = std::max(Err_Y_Low, Err_Y_High); lMax = lbl.first; }
          }
        }
        //
        // Draw the graphs in the same canvas
        bool firstPlot = true; uint iCnt = 1;
        for (auto& gr : grMap) {
          const auto& grLbl = gr.first;
          auto& graph = gr.second;
          // Format the graphs
          const bool incAcc = (accType!="");
          formatResultsGraph(graph, col, var, chg, useEtaCM, incAcc, true);
          graph.SetMarkerColor(kBlack);
          graph.SetMarkerSize(0.0);
          graph.SetLineColor((iCnt==uint(kWhite) ? uint(kOrange) : iCnt));
          graph.SetLineWidth(3);
          graph.SetFillStyle(0);
          if (var == "Charge_Asymmetry"     ) { graph.GetYaxis()->SetRangeUser(-0.03, 0.10); }
          if (var == "Cross_Section"        ) { graph.GetYaxis()->SetRangeUser(-0.1,  0.3); }
          if (var == "ForwardBackward_Ratio") { graph.GetYaxis()->SetRangeUser(-0.05,  0.1); }
          iCnt++;
          // Create Legend
          formatLegendEntry(*leg.AddEntry(&graph, grLbl.c_str(), "f"));
          if (grLbl==lMax) { setRangeYAxisGraph(graph, 0.02, 0.55); }
          TGaxis::SetMaxDigits(3); // to display powers of 10
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
        for (const auto& s: textToPrint) { tex.DrawLatex(0.22, 0.86-dy, s.c_str()); dy+=0.045; }
        if (var=="Cross_Section") { tex.SetTextSize(0.025); tex.DrawLatex(0.58, 0.72, "Luminosity uncertainty: 5.0%"); }
        // Update
        c.Modified(); c.Update(); // Pure paranoia
        //
        // set the CMS style
        int option = 118;
        if (col.find("pPb")!=std::string::npos) option = 115;
        if (col.find("Pbp")!=std::string::npos) option = 116;
        CMS_lumi(&c, option, 33, "");
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
        const std::string plotDir = outDir+"/Plots/Systematic/" + col+"/" + label+"/" + var;
        makeDir(plotDir + "/png/");
        makeDir(plotDir + "/pdf/");
        makeDir(plotDir + "/root/");
        //
        // Save Canvas
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


void drawGraphWithTheory( GraphPentaMap& graphMap , const std::string& outDir , const bool useEtaCM = true , const std::string accType = "MC" , const std::string effType = "TnP" )
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
	std::string sampleLabel = "W #rightarrow #mu + #nu_{#mu}";
	if (chg == "Pl") { sampleLabel = "W^{+} #rightarrow #mu^{+} + #nu_{#mu}"; }
	if (chg == "Mi") { sampleLabel = "W^{-} #rightarrow #mu^{-} + #nu_{#mu}"; }
	textToPrint.push_back(sampleLabel);
	if (accType=="") { textToPrint.push_back("p^{#mu}_{T} > 25 GeV/c"); }
	//
	// Declare the graph vector (for drawing with markers)
	std::vector< TGraphAsymmErrors > grVec;
	// Initialize the Legend
	double legOff = 0.0; if (accType=="") { legOff = 0.05; }
	TLegend leg(0.2, (0.71 - legOff), 0.4, (0.84 - legOff));
	// Initialize the graph x range variables
	double xMin=0.0 , xMax=0.0 , yErrMin=9999999. , yErrMax=-1.;
	// Format the graphs
	const bool incAcc = (accType!="");
	grVec.push_back(nomGraph.at("Err_Tot"));
	grVec.push_back(theGraph.at("CT14"));
	grVec.push_back(theGraph.at("EPPS16"));
	for (auto& gr : grVec) { formatResultsGraph(gr, col, var, chg, useEtaCM, incAcc); }
	grVec[0].SetMarkerColor(kBlack);
	grVec[1].SetFillColor(kYellow);
	grVec[1].SetLineColor(kRed);
	grVec[1].SetFillStyle(1001);
	grVec[1].SetLineStyle(1);
	grVec[1].SetLineWidth(3);
	grVec[1].SetMarkerSize(0);
	grVec[2].SetFillColor(kGreen+2);
	grVec[2].SetLineColor(kGreen+2);
	grVec[2].SetFillStyle(3005);
	grVec[2].SetLineStyle(7);
	grVec[2].SetLineWidth(3);
	grVec[2].SetMarkerSize(0);
        //
        auto h1 = graphToHist(grVec[1]); for (int i=1; i<=grVec[1].GetN(); i++) { h1.SetBinError(i, 0.0001); }
        auto h2 = graphToHist(grVec[2]); for (int i=1; i<=grVec[2].GetN(); i++) { h2.SetBinError(i, 0.0001); }
        //
	// Create Legend
	formatLegendEntry(*leg.AddEntry(&grVec[0], "Data", "pe"));
	formatLegendEntry(*leg.AddEntry(&grVec[1], "CT14", "lf"));
	formatLegendEntry(*leg.AddEntry(&grVec[2], "CT14+EPPS16", "lf"));
	// Draw the graphs
	grVec[0].Draw("ap");
	grVec[1].Draw("same2"); h1.Draw("sameL");
	grVec[2].Draw("same2"); h2.Draw("sameL");
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
	for (const auto& s: textToPrint) { tex.DrawLatex(0.22, 0.86-dy, s.c_str()); dy+=0.045; }
        if (var=="Cross_Section") { tex.SetTextSize(0.025); tex.DrawLatex(0.58, 0.72, "Luminosity uncertainty: 5.0%"); }
	// Update
	c.Modified(); c.Update(); // Pure paranoia
	//
	// set the CMS style
	int option = 118;
	if (col.find("pPb")!=std::string::npos) option = 115;
	if (col.find("Pbp")!=std::string::npos) option = 116;
	CMS_lumi(&c, option, 33, "");
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
	const std::string plotDir = outDir+"/Plots/Theory/" + col+"/" + label+"/" + var;
	makeDir(plotDir + "/png/");
	makeDir(plotDir + "/pdf/");
	makeDir(plotDir + "/root/");
	//
	// Save Canvas
	const std::string name = Form("gr_WToMu%s_%s_%s_%s_%s", chg.c_str(), col.c_str(), var.c_str(), label.c_str(), "NominalWithTheory");
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


void drawGraphWithpPb5TeV( GraphPentaMap& graphMap , const std::string& outDir , const bool useEtaCM = true , const std::string accType = "MC" , const std::string effType = "TnP" )
{
  //
  // Set Style
  setStyle();
  //
  std::cout << "[INFO] Drawing the output graphs with 5 TeV pPb results" << std::endl;
  //
  // Add the 5 TeV pPb results
  HIN_13007_Wp (graphMap["PA"]["Pl"]["Cross_Section"]["CMS"]);
  HIN_13007_Wm (graphMap["PA"]["Mi"]["Cross_Section"]["CMS"]);
  HIN_13007_chasym(graphMap["PA"][""]["Charge_Asymmetry"]["CMS"]);
  //
  // Draw all graphs
  for (const auto& c : graphMap) {
    for (const auto& ch : c.second) {
      for (const auto& v : ch.second) {
        //
        if (v.first=="ForwardBackward_Ratio") continue;
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
	std::string sampleLabel = "W #rightarrow #mu + #nu_{#mu}";
	if (chg == "Pl") { sampleLabel = "W^{+} #rightarrow #mu^{+} + #nu_{#mu}"; }
	if (chg == "Mi") { sampleLabel = "W^{-} #rightarrow #mu^{-} + #nu_{#mu}"; }
	textToPrint.push_back(sampleLabel);
	if (accType=="") { textToPrint.push_back("p^{#mu}_{T} > 25 GeV/c"); }
	//
	// Declare the graph vector (for drawing with markers)
	std::vector< TGraphAsymmErrors > grVec;
	// Initialize the Legend
	double legOff = 0.0; if (accType=="") { legOff = 0.05; }
	TLegend leg(0.2, (0.71 - legOff), 0.4, (0.84 - legOff));
	// Initialize the graph x range variables
	double xMin=0.0 , xMax=0.0 , yErrMin=9999999. , yErrMax=-1.;
	// Format the graphs
	const bool incAcc = (accType!="");
	grVec.push_back(nomGraph.at("Err_Tot"));
	grVec.push_back(theGraph.at("HIN-13007"));
	for (auto& gr : grVec) { formatResultsGraph(gr, col, var, chg, useEtaCM, incAcc); }
	grVec[0].SetMarkerColor(kBlack);
	grVec[1].SetMarkerColor(kRed);
	grVec[1].SetMarkerStyle(21);
	grVec[1].SetMarkerSize(2.);
	// Create Legend
	formatLegendEntry(*leg.AddEntry(&grVec[0], "Data", "pe"));
	formatLegendEntry(*leg.AddEntry(&grVec[1], "HIN-13007", "pe"));
	// Draw the graphs
	grVec[0].Draw("ap");
	grVec[1].Draw("samep");
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
	for (const auto& s: textToPrint) { tex.DrawLatex(0.22, 0.86-dy, s.c_str()); dy+=0.045; }
        if (var=="Cross_Section") { tex.SetTextSize(0.025); tex.DrawLatex(0.58, 0.72, "Luminosity uncertainty: 5.0%"); }
	// Update
	c.Modified(); c.Update(); // Pure paranoia
	//
	// set the CMS style
	int option = 118;
	if (col.find("pPb")!=std::string::npos) option = 115;
	if (col.find("Pbp")!=std::string::npos) option = 116;
	CMS_lumi(&c, option, 33, "");
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
	const std::string plotDir = outDir+"/Plots/HIN_13007/" + col+"/" + label+"/" + var;
	makeDir(plotDir + "/png/");
	makeDir(plotDir + "/pdf/");
	makeDir(plotDir + "/root/");
	//
	// Save Canvas
	const std::string name = Form("gr_WToMu%s_%s_%s_%s_%s", chg.c_str(), col.c_str(), var.c_str(), label.c_str(), "HIN_13007");
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
  if (n >= 0.) { return "+"; } else { return ""; }
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
        val = Form("%s%.2f , %s%.2f", sgn(min), min, sgn(max) , max);
      }
      else if (v=="N_DS_Entries") {
        val = Form("%g", var.at("Val"));
      }
      else if (v=="N_FIT_Entries") {
        val = Form("%.0f", var.at("Val"));
      }
      else if (v=="TEST_FIT") {
        //const double Chi2 = var.at("Chi2");
        //const double NDoF = var.at("NDoF");
        //const double pVal = TMath::Prob(Chi2, NDoF);
        const double pVal = var.at("Val");
        val = Form("%.2f", pVal);
      }
      else if (v=="Acceptance_MC" || v=="Efficiency_MC" || v=="Efficiency_TnP") {
        const double errLow = var.at("Err_Stat_Low")*100.;
        const double errHi  = var.at("Err_Stat_High")*100.;
        if (errLow==errHi) { val = Form("$%.2f \\pm %.2f$", var.at("Val")*100., errLow ); }
        else { val = Form("$%.2f + %.2f - %.2f$", var.at("Val")*100., errHi, errLow ); }
      }
      else {
        const double errLow = var.at("Err_Stat_Low");
        const double errHi  = var.at("Err_Stat_High");
        if (errLow==errHi) { val = Form("$%.0f \\pm %.0f$", var.at("Val"), errLow ); }
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
  const std::vector< std::string > allColVar    = { "Muon_Eta" , "N_DS_Entries" , "N_FIT_Entries" , "N_WToMu_RAW" , "N_DYToMu" , "N_WToTauToMu" , "N_DYToTauToMu" , "N_TTbarToMu" , "N_QCDToMu" , "TEST_FIT" };
  std::vector< std::string > allColTitle1 = { "$\\eta_{LAB}$ Range" , "Total"  , "Fitted" , "Signal" , "$Z/\\gamma*\\to\\mu^{+}+\\mu^{-}$" , "$W\\to\\tau+\\nu_{\\tau}$" ,
                                              "$Z/\\gamma*\\to\\tau^{+}+\\tau^{-}$" , "$t\\bar{t}$" , "QCD" , "$\\chi^{2}/NDoF$" };
  if (useEtaCM) { allColTitle1[0] = "$\\eta_{CM}$ Range"; }
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
  texTable.push_back("\\begin{table}[h!]");
  texTable.push_back("  \\centering");
  texTable.push_back("  \\resizebox{\\textwidth}{!}{");
  createYieldTable(texTable, colVar, colTitle1, inputVar, col, chg);
  texTable.push_back("  }");
  texTable.push_back(Form("  \\caption{%s}",
                          Form("Raw yields of %s and background processes, extracted from the nominal fits for each %s bin in the %s collision system. All analysis cuts are applied%s. All uncertainties shown are statistical only.",
                               (chg=="Pl" ? "\\WToMuNuPl" : "\\WToMuNuMi"),
                               (useEtaCM ? "muon $\\eta_{CM}$" : "$\\eta_{LAB}$"),
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
  std::vector< std::string > allColTitle1 = { "$\\eta_{LAB}$ Range" , "Raw Yield"  , "MC Acceptance ($\\%$)" , "MC Efficiency ($\\%$)" , "Efficiency ($\\%$)" , "Corrected Yield" };
  if (useEtaCM) { allColTitle1[0] = "$\\eta_{CM}$ Range"; }
  std::vector< std::string > colVar , colTitle1, colTitle2;
  for(uint i = 0; i < allColVar.size(); i++) {
    if (inVar.count(allColVar[i])>0) {
      colVar.push_back(allColVar[i]);
      colTitle1.push_back(allColTitle1[i]);
    }
  }
  //
  texTable.push_back("\\begin{table}[h!]");
  texTable.push_back("  \\centering");
  //texTable.push_back("  \\resizebox{\\textwidth}{!}{");
  createYieldTable(texTable, colVar, colTitle1, inputVar, col, chg);
  //texTable.push_back("  }");
  texTable.push_back(Form("  \\caption{%s}",
                          Form("Corrected yields of %s, given for each %s bin in the %s collision system. All analysis cuts are applied%s.%s All uncertainties shown are statistical only.",
                               (chg=="Pl" ? "\\WToMuNuPl" : "\\WToMuNuMi"),
                               (useEtaCM ? "muon $\\eta_{CM}$" : "$\\eta_{LAB}$"),
                               (col=="PA" ? "combined \\pPb and \\Pbp" : (col=="pPb" ? "\\pPb" : (col=="Pbp" ? "\\Pbp" : "????"))),
                               (pTCUT!="" ? Form(" including the muon %s cut", pTCUT.c_str()) : ""),
                               (inVar.count("Efficiency_TnP")>0 ? " The muon efficiency has been corrected by applying the Tag and Probe scale factors event by event." : "")
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
	  val = Form("%s%.2f , %s%.2f", sgn(min), min, sgn(max) , max);
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
  const double lumiVal = ( (inVar.count("Luminosity")>0) ? inVar.at("Luminosity").at("Val")      : -1. );
  const double lumiErr = ( (inVar.count("Luminosity")>0) ? inVar.at("Luminosity").at("Val")*0.05 : -1. ); // 5% error
  //
  // Initialize the latex table
  std::vector< std::string > texTable;
  //
  // Determine number of columns
  const std::vector< std::string > colChg = { "Mi" , "Mi" , "Mi" , "Pl" , "Pl" };
  const std::vector< std::string > colEta = { "" , "" , "" , "" , "" };
  const std::vector< std::string > colVar = { "Muon_Eta" , "N_WToMu" , "Cross_Section" , "N_WToMu" , "Cross_Section" };
  std::vector< std::string >    colTitle1 = { "$\\eta_{LAB}$ Range" , "$\\W^{-} $ Yield"  , "" , "$\\W^{+} $ Yield" , "" };
  if (useEtaCM) { colTitle1[0] = "$\\eta_{CM}$ Range"; }
  colTitle1[2] = Form("$\\W^{-} %s$", formatResultVarName("Cross_Section", useEtaCM, false, true).c_str());
  colTitle1[4] = Form("$\\W^{+} %s$", formatResultVarName("Cross_Section", useEtaCM, false, true).c_str());
  //
  texTable.push_back("\\begin{table}[h!]");
  texTable.push_back("  \\centering");
  texTable.push_back("  \\resizebox{\\textwidth}{!}{");
  createResultTable(texTable, colVar, colTitle1, colChg, colEta, inputVar, resultVar, col);
  texTable.push_back("  }");
  texTable.push_back(Form("  \\caption{%s}",
                          Form("Cross-Section of \\WToMuNu, given for each %s bin in the %s collision system. The total integrated luminosity of the sample corresponds to $%.1f \\pm %.1f$~\nbinv. All analysis cuts are applied%s. All uncertainties shown are statistical only.",
                               (useEtaCM ? "muon $\\eta_{CM}$" : "$\\eta_{LAB}$"),
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
  std::vector< std::string >    colTitle1 = { "$\\eta_{LAB}$ Range" , "$\\W^{-} $ Fwd Yield" , "$\\W^{-} $ Bwd Yield" , "$\\W^{+} $ Fwd Yield" , "$\\W^{+} $ Bwd Yield"  , "" , "" , "" };
  if (useEtaCM) { colTitle1[0] = "$\\eta_{CM}$ Range"; }
  colTitle1[5] = Form("$\\W^{-} %s$", formatResultVarName("ForwardBackward_Ratio", useEtaCM, false, true).c_str());
  colTitle1[6] = Form("$\\W^{+} %s$", formatResultVarName("ForwardBackward_Ratio", useEtaCM, false, true).c_str());
  colTitle1[7] = Form("$\\W %s$", formatResultVarName("ForwardBackward_Ratio", useEtaCM, false, true).c_str());
  //
  texTable.push_back("\\begin{table}[h!]");
  texTable.push_back("  \\centering");
  texTable.push_back("  \\resizebox{\\textwidth}{!}{");
  createResultTable(texTable, colVar, colTitle1, colChg, colEta, inputVar, resultVar, col);
  texTable.push_back("  }");
  texTable.push_back(Form("  \\caption{%s}",
                          Form("Forward-Backward ratio of \\WToMuNu, given for each %s bin in the %s collision system. All analysis cuts are applied%s. All uncertainties shown are statistical only.",
                               (useEtaCM ? "muon $\\eta_{CM}$" : "$\\eta_{LAB}$"),
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
  std::vector< std::string >    colTitle1 = { "$\\eta_{LAB}$ Range" , "$\\W^{+} $ Yield" , "$\\W^{-} $ Yield" , "" };
  if (useEtaCM) { colTitle1[0] = "$\\eta_{CM}$ Range"; }
  colTitle1[3] = Form("$%s$", formatResultVarName("Charge_Asymmetry", useEtaCM, false, true).c_str());
  //
  texTable.push_back("\\begin{table}[h!]");
  texTable.push_back("  \\centering");
  //texTable.push_back("  \\resizebox{\\textwidth}{!}{");
  createResultTable(texTable, colVar, colTitle1, colChg, colEta, inputVar, resultVar, col);
  //texTable.push_back("  }");
  texTable.push_back(Form("  \\caption{%s}",
                          Form("\\W charge asymmetry, given for each %s bin in the %s collision system. All analysis cuts are applied%s. All uncertainties shown are statistical only.",
                               (useEtaCM ? "muon $\\eta_{CM}$" : "$\\eta_{LAB}$"),
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
        val = Form("%s%.2f , %s%.2f", sgn(min), min, sgn(max) , max);
      }
      else if (v=="Cross_Section" || v=="ForwardBackward_Ratio" || v=="Charge_Asymmetry") {
	const double statErrLow = varMap.at(col).at(chg).at(v).at("Err_Stat_Low").at(b.first);
	const double statErrHi  = varMap.at(col).at(chg).at(v).at("Err_Stat_Low").at(b.first);
	const double systErrLow = varMap.at(col).at(chg).at(v).at("Err_Syst_Low").at(b.first);
	const double systErrHi  = varMap.at(col).at(chg).at(v).at("Err_Syst_High").at(b.first);
	const double rVal       = varMap.at(col).at(chg).at(v).at("Val").at(b.first);
	if (v=="Cross_Section") {
	  if (statErrLow==statErrHi) { val = Form("$%.2f \\pm %.2f$ (stat)", rVal, statErrLow ); }
	  else { val = Form("$%.2f + %.2f - %.2f$ (stat)", rVal, statErrHi, statErrLow ); }
	  if (systErrLow==systErrHi) { val += Form("$ \\pm %.2f$ (syst)", systErrLow ); }
	  else { val += Form("$ + %.2f - %.2f$ (syst)", systErrHi, systErrLow ); }
	}
	else if (v=="ForwardBackward_Ratio" || v=="Charge_Asymmetry") {
	  if (statErrLow==statErrHi) { val = Form("$%.3f \\pm %.3f$ (stat)", rVal, statErrLow ); }
	  else { val = Form("$%.2f + %.3f - %.3f$ (stat)", rVal, statErrHi, statErrLow ); }
	  if (systErrLow==systErrHi) { val += Form("$ \\pm %.3f$ (syst)", systErrLow ); }
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
  const double lumiErr = 0.05*lumiVal;
  //
  // Determine number of columns
  const std::vector< std::string > colChg = { "Mi" , "Mi" , "Pl" };
  const std::vector< std::string > colVar = { "Muon_Eta" , "Cross_Section" , "Cross_Section" };
  std::vector< std::string >    colTitle1 = { "$\\eta_{LAB}$ Range" , "" , "" };
  if (useEtaCM) { colTitle1[0] = "$\\eta_{CM}$ Range"; }
  colTitle1[1] = Form("$\\W^{-} %s$", formatResultVarName("Cross_Section", useEtaCM, false, true).c_str());
  colTitle1[2] = Form("$\\W^{+} %s$", formatResultVarName("Cross_Section", useEtaCM, false, true).c_str());
  //
  texTable.push_back("\\begin{table}[h!]");
  texTable.push_back("  \\centering");
  //texTable.push_back("  \\resizebox{\\textwidth}{!}{");
  createSystematicTable(texTable, colVar, colTitle1, colChg, varMap, col);
  //texTable.push_back("  }");
  texTable.push_back(Form("  \\caption{%s}",
                          Form("Cross-Section of \\WToMuNu, given for each %s bin in the %s collision system. The total integrated luminosity of the sample corresponds to $%.1f \\pm %.1f$~\nbinv. All analysis cuts are applied%s. The global luminosity uncertainty of $\\pm$5.0$\\%%$ is not included.",
                               (useEtaCM ? "muon $\\eta_{CM}$" : "$\\eta_{LAB}$"),
                               (col=="PA" ? "combined \\pPb and \\Pbp" : (col=="pPb" ? "\\pPb" : (col=="Pbp" ? "\\Pbp" : "????"))),
                               lumiVal, lumiErr,
                               (pTCUT!="" ? Form(" including the muon %s cut", pTCUT.c_str()) : "")
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
  std::vector< std::string >    colTitle1 = { "$\\eta_{LAB}$ Range" , "" , "" , "" };
  if (useEtaCM) { colTitle1[0] = "$\\eta_{CM}$ Range"; }
  colTitle1[1] = Form("$\\W^{-} %s$", formatResultVarName("ForwardBackward_Ratio", useEtaCM, false, true).c_str());
  colTitle1[2] = Form("$\\W^{+} %s$", formatResultVarName("ForwardBackward_Ratio", useEtaCM, false, true).c_str());
  colTitle1[3] = Form("$\\W %s$", formatResultVarName("ForwardBackward_Ratio", useEtaCM, false, true).c_str());
  //
  texTable.push_back("\\begin{table}[h!]");
  texTable.push_back("  \\centering");
  texTable.push_back("  \\resizebox{\\textwidth}{!}{");
  createSystematicTable(texTable, colVar, colTitle1, colChg, varMap, col);
  texTable.push_back("  }");
  texTable.push_back(Form("  \\caption{%s}",
                          Form("Forward-Backward ratio of \\WToMuNu, given for each %s bin in the %s collision system. All analysis cuts are applied%s.",
                               (useEtaCM ? "muon $\\eta_{CM}$" : "$\\eta_{LAB}$"),
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
  std::vector< std::string >    colTitle1 = { "$\\eta_{LAB}$ Range" , "" };
  if (useEtaCM) { colTitle1[0] = "$\\eta_{CM}$ Range"; }
  colTitle1[1] = Form("$%s$", formatResultVarName("Charge_Asymmetry", useEtaCM, false, true).c_str());
  //
  texTable.push_back("\\begin{table}[h!]");
  texTable.push_back("  \\centering");
  //texTable.push_back("  \\resizebox{\\textwidth}{!}{");
  createSystematicTable(texTable, colVar, colTitle1, colChg, varMap, col);
  //texTable.push_back("  }");
  texTable.push_back(Form("  \\caption{%s}",
                          Form("\\W charge asymmetry, given for each %s bin in the %s collision system. All analysis cuts are applied%s",
                               (useEtaCM ? "muon $\\eta_{CM}$" : "$\\eta_{LAB}$"),
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
  //if (true) { // ADD LUMINOSITY
  //  texTable.push_back("    Luminosity & 0.050 & 0.050  & 0.000 & 0.000 & 0.000 & 0.000 \\\\");
  //  texTable.push_back("    \\hline");
  //}
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
  const std::string pTCUT = Form("$p_{T} > %.0f$~GeV/c", 25.0);
  //
  // Determine number of columns
  const std::vector< std::string > colChg = { "" , "Mi" , "Pl" , "Mi" , "Pl" , "" , "" };
  const std::vector< std::string > colVar = { "Title_Syst" , "Cross_Section" , "Cross_Section" , "ForwardBackward_Ratio" , "ForwardBackward_Ratio" , "ForwardBackward_Ratio" , "Charge_Asymmetry" };
  std::vector< std::string >    colTitle1 = { "Systematic Variation" , "" , "" , "" , "" , "" , "" };
  colTitle1[1] = Form("$\\W^{-} %s$", formatResultVarName("Cross_Section", useEtaCM, false, true).c_str());
  colTitle1[2] = Form("$\\W^{+} %s$", formatResultVarName("Cross_Section", useEtaCM, false, true).c_str());
  colTitle1[3] = Form("$\\W^{-} %s$", formatResultVarName("ForwardBackward_Ratio", useEtaCM, false, true).c_str());
  colTitle1[4] = Form("$\\W^{+} %s$", formatResultVarName("ForwardBackward_Ratio", useEtaCM, false, true).c_str());
  colTitle1[5] = Form("$\\W %s$", formatResultVarName("ForwardBackward_Ratio", useEtaCM, false, true).c_str());
  colTitle1[6] = Form("$%s$", formatResultVarName("Charge_Asymmetry", useEtaCM, false, true).c_str());
  //
  texTable.push_back("\\begin{table}[h!]");
  texTable.push_back("  \\centering");
  texTable.push_back("  \\resizebox{\\textwidth}{!}{");
  createFullSystematicTable(texTable, colVar, colTitle1, colChg, varMap, col);
  texTable.push_back("  }");
  texTable.push_back(Form("  \\caption{%s}",
                          Form("Maximum error of the measured observables determined for each category in the %s collision system. The uncertainties of the cross-sections are relative while for the asymmetries are absolute. All analysis cuts are applied%s.",
                               (col=="PA" ? "combined \\pPb and \\Pbp" : (col=="pPb" ? "\\pPb" : (col=="Pbp" ? "\\Pbp" : "????"))),
                               (pTCUT!="" ? Form(" including the muon %s cut", pTCUT.c_str()) : "")
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


#endif // ifndef resultUtils_h
  

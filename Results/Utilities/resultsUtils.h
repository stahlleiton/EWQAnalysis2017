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
#include "TVector.h"
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
#include <map>
#include <array>
// CMS headers
#include "../../Utilities/CMS/tdrstyle.C"
#include "../../Utilities/CMS/CMS_lumi.C"



// ------------------ TYPE -------------------------------
typedef std::map< std::string , std::string       > StringMap;
typedef std::map< std::string , std::string*      > StringPMap;
typedef std::map< std::string , bool              > BoolMap;
typedef std::map< std::string , TGraphAsymmErrors > GraphMap;
typedef std::map< std::string , GraphMap          > GraphMapMap;
typedef std::map< std::string , GraphMapMap       > GraphTriMap;
typedef std::map< std::string , GraphTriMap       > GraphQuadMap;
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
        std::cout << "[INFO] Adding file: " << epdf->d_name << std::endl;
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
  if (thePoiNamesStr == "all") { thePoiNamesStr = "N_WToMu,N_WToTauToMu,N_DYToMu,N_QCDToMu,N_TTbarToMu,Alpha_QCDToMu,Beta_QCDToMu,x0_QCDToMu"; }
  thePoiNamesStr = thePoiNamesStr + ","; // Just to include the last name XD
  std::vector< std::string > thePoiNamesVec;
  while (thePoiNamesStr.find(",")!=std::string::npos) {
    thePoiNamesVec.push_back( thePoiNamesStr.substr(0, thePoiNamesStr.find(",")) );
    thePoiNamesStr = thePoiNamesStr.substr(thePoiNamesStr.find(",")+1);
  }
  const std::vector< std::string > poiVarType = { "Min" , "Max" , "Val" , "Err" , "parIni_Val" , "parIni_Err" };
  //
  for (const auto& p : thePoiNamesVec) {
    for (const auto& t : poiVarType) {
      info.Var["POI_"+p][t] = -99.0;
    }
  }
  //
  // Initialize the model name containers
  const std::vector< std::string > objType = { "W" , "WToTau" , "DY" , "TTbar" , "QCD" };
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
  for (const auto& c : inputVar) {
    for (const auto& b : c.second) {
      for (const auto& ch : b.second) {
        // For MC Acceptance
        eff[c.first][ch.first]["Acceptance_MC"]["Val"][b.first] = -1.0;
        eff[c.first][ch.first]["Acceptance_MC"]["Err_Stat_High"][b.first] = 0.0;
        eff[c.first][ch.first]["Acceptance_MC"]["Err_Stat_Low" ][b.first] = 0.0;
        eff[c.first][ch.first]["Acceptance_MC"]["Err_Syst_High"][b.first] = 0.0;
        eff[c.first][ch.first]["Acceptance_MC"]["Err_Syst_Low" ][b.first] = 0.0;
        // For MC Efficiency
        eff[c.first][ch.first]["Efficiency_MC"]["Val"][b.first] = -1.0;
        eff[c.first][ch.first]["Efficiency_MC"]["Err_Stat_High"][b.first] = 0.0;
        eff[c.first][ch.first]["Efficiency_MC"]["Err_Stat_Low" ][b.first] = 0.0;
        eff[c.first][ch.first]["Efficiency_MC"]["Err_Syst_High"][b.first] = 0.0;
        eff[c.first][ch.first]["Efficiency_MC"]["Err_Syst_Low" ][b.first] = 0.0;
        // For Corrected Efficiency
        eff[c.first][ch.first]["Efficiency_TnP"]["Val"][b.first] = -1.0;
        eff[c.first][ch.first]["Efficiency_TnP"]["Err_Stat_High"][b.first] = 0.0;
        eff[c.first][ch.first]["Efficiency_TnP"]["Err_Stat_Low" ][b.first] = 0.0;
        eff[c.first][ch.first]["Efficiency_TnP"]["Err_Syst_High"][b.first] = 0.0;
        eff[c.first][ch.first]["Efficiency_TnP"]["Err_Syst_Low" ][b.first] = 0.0;
      }
    }
  }
};


void getEffContent( double& val , double& err_High , double& err_Low , const TEfficiency& eff , const double& binVal )
{
  auto hist = eff.GetTotalHistogram();
  if (hist==NULL) { return; }
  const int iBin = hist->GetXaxis()->FindBin(binVal);
  val      = eff.GetEfficiency(iBin);
  err_High = eff.GetEfficiencyErrorUp(iBin);
  err_Low  = eff.GetEfficiencyErrorLow(iBin);
};


void getUncContent( double& err_High , double& err_Low , const TVector& unc , const TEfficiency& eff , const double& binVal )
{
  auto hist = eff.GetTotalHistogram();
  if (hist==NULL) { return; }
  const int idx = ( hist->GetXaxis()->FindBin(binVal) - 1 ); // index is bin number - 1
  const double unc_High = unc[idx];
  const double unc_Low  = unc[idx];
  err_High = std::sqrt( std::pow( err_High , 2.0 ) + std::pow( unc_High , 2.0 ) );
  err_Low  = std::sqrt( std::pow( err_Low  , 2.0 ) + std::pow( unc_Low  , 2.0 ) );
};


bool getAcceptanceAndEfficiency( BinPentaMap& effMap , const std::string& inputFilePath , const bool useEtaCM = true )
{
  //
  // Open the input file
  TFile inputFile(inputFilePath.c_str(), "READ");
  if (inputFile.IsOpen()==false || inputFile.IsZombie()==true) { std::cout << "[ERROR] The input efficiency file " << inputFilePath << " was not found!" << std::endl; return false; }
  inputFile.cd();
  //
  std::string var = "Eta"; if (useEtaCM) { var = "EtaCM"; }
  const std::string sample = "MC_WToMuNu";
  //
  for (const auto& c : effMap) {
    for (const auto& ch : c.second) {
      if (ch.second.count("Acceptance_MC")==0 || ch.second.at("Acceptance_MC").count("Val")==0) { std::cout << "[ERROR] The efficiency container is not valid" << std::endl; return false; }
      auto& eff = effMap.at(c.first).at(ch.first);
      //
      const std::string col = c.first;
      std::string charge = ""; if (ch.first=="Pl") { charge = "Plus"; } if (ch.first=="Mi") { charge = "Minus"; }
      //
      const std::string accDir    = Form("TnPEfficiency1D/%s/%s/%s/Acceptance/NoCorr" , var.c_str(), sample.c_str(), col.c_str());
      const std::string effDir    = Form("TnPEfficiency1D/%s/%s/%s/Total/NoCorr"      , var.c_str(), sample.c_str(), col.c_str());
      const std::string effTnPDir = Form("TnPEfficiency1D/%s/%s/%s/Total/TnP_Nominal" , var.c_str(), sample.c_str(), col.c_str());
      const std::string uncTnPDir = Form("TnPEfficiency1D/%s/%s/%s/Total"             , var.c_str(), sample.c_str(), col.c_str());
      //
      const std::string accName        = Form("eff1D_%s_%s_%s_%s_Acceptance_NoCorr" , var.c_str(), sample.c_str(), col.c_str(), charge.c_str());
      const std::string effName        = Form("eff1D_%s_%s_%s_%s_Total_NoCorr"      , var.c_str(), sample.c_str(), col.c_str(), charge.c_str());
      const std::string effTnPName     = Form("eff1D_%s_%s_%s_%s_Total_TnP_Nominal" , var.c_str(), sample.c_str(), col.c_str(), charge.c_str());
      const std::string uncStatTnPName = Form("unc1D_%s_%s_%s_%s_Total_TnP_Stat"    , var.c_str(), sample.c_str(), col.c_str(), charge.c_str());
      const std::string uncSystTnPName = Form("unc1D_%s_%s_%s_%s_Total_TnP_Syst"    , var.c_str(), sample.c_str(), col.c_str(), charge.c_str());
      //
      auto accObj        = (TEfficiency*) inputFile.Get(Form("%s/%s", accDir.c_str(), accName.c_str()));
      auto effObj        = (TEfficiency*) inputFile.Get(Form("%s/%s", effDir.c_str(), effName.c_str()));
      auto effTnPObj     = (TEfficiency*) inputFile.Get(Form("%s/%s", effTnPDir.c_str(), effTnPName.c_str()));
      auto uncStatTnPObj = (TVector*    ) inputFile.Get(Form("%s/%s", uncTnPDir.c_str(), uncStatTnPName.c_str()));
      auto uncSystTnPObj = (TVector*    ) inputFile.Get(Form("%s/%s", uncTnPDir.c_str(), uncSystTnPName.c_str()));
      //
      if (accObj==NULL       ) { std::cout << "[ERROR] Acceptance object in "                  << accDir    << " was not found in file " << inputFilePath << std::endl; return false; }
      if (effObj==NULL       ) { std::cout << "[ERROR] Efficiency object in "                  << effDir    << " was not found in file " << inputFilePath << std::endl; return false; }
      if (effTnPObj==NULL    ) { std::cout << "[ERROR] TnP Efficiency object in "              << effTnPDir << " was not found in file " << inputFilePath << std::endl; return false; }
      if (uncStatTnPObj==NULL) { std::cout << "[ERROR] TnP Statistical Uncertainty object in " << uncTnPDir << " was not found in file " << inputFilePath << std::endl; return false; }
      if (uncSystTnPObj==NULL) { std::cout << "[ERROR] TnP Systematic Uncertainty object in "  << uncTnPDir << " was not found in file " << inputFilePath << std::endl; return false; }
      //
      for (const auto& b : ch.second.at("Acceptance_MC").at("Val")) {
        // Compute the center value of the bin
        const double etaVal = ( ( b.first.etabin().high() + b.first.etabin().low() ) / 2.0 );
        // Fill Acceptance
        getEffContent(eff.at("Acceptance_MC").at("Val").at(b.first) , eff.at("Acceptance_MC").at("Err_Stat_High").at(b.first) , eff.at("Acceptance_MC").at("Err_Stat_Low").at(b.first) , *accObj , etaVal);
        // Fill Efficiency
        getEffContent(eff.at("Efficiency_MC").at("Val").at(b.first) , eff.at("Efficiency_MC").at("Err_Stat_High").at(b.first) , eff.at("Efficiency_MC").at("Err_Stat_Low").at(b.first) , *effObj , etaVal);
        // Fill TnP Efficiency
        getEffContent(eff.at("Efficiency_TnP").at("Val").at(b.first) , eff.at("Efficiency_TnP").at("Err_Stat_High").at(b.first) , eff.at("Efficiency_TnP").at("Err_Stat_Low").at(b.first) , *effTnPObj , etaVal);
        // Fill TnP Statistical Uncertainty
        getUncContent(eff.at("Efficiency_TnP").at("Err_Stat_High").at(b.first) , eff.at("Efficiency_TnP").at("Err_Stat_Low").at(b.first) , *uncStatTnPObj , *effTnPObj , etaVal);
        // Fill TnP Systematic Uncertainty
        getUncContent(eff.at("Efficiency_TnP").at("Err_Syst_High").at(b.first) , eff.at("Efficiency_TnP").at("Err_Syst_Low").at(b.first) , *uncSystTnPObj , *effTnPObj , etaVal);
      }
    }
  }
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


bool correctRawYields( VarBinMap& inputVar , const BinPentaMap& effMap , const std::string accType = "MC" , const std::string effType = "TnP" )
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
          Acceptance["Val"] = 1.0;
          for (const auto& t : ch.second.at("N_WToMu")) { if (t.first=="Val") continue; Acceptance[t.first] = 0.0; }
        }
        if (effName!="") { for (const auto& t : eff.at(effName)) { Efficiency[t.first] = t.second.at(b.first); } }
        else {
          Efficiency["Val"] = 1.0;
          for (const auto& t : ch.second.at("N_WToMu")) { if (t.first=="Val") continue; Efficiency[t.first] = 0.0; }
        }
        //
        const auto& N_Raw  = ch.second.at("N_WToMu");
        auto&       N_Corr = inputVar.at(c.first).at(b.first).at(ch.first).at("N_WToMu");
        // Create a backup
        for (const auto& t : N_Raw) { inputVar.at(c.first).at(b.first).at(ch.first)["N_WToMu_RAW"][t.first] = t.second; }
        //
        // Fill with the corrected values
        N_Corr.at("Val") = getCorrectedYieldValue(N_Raw.at("Val") , Acceptance.at("Val") , Efficiency.at("Val"));
        for (const auto& t : N_Raw) {
          if (t.first=="Val") continue;
          N_Corr.at(t.first) = getCorrectedYieldError(N_Raw.at("Val") , Acceptance.at("Val") , Efficiency.at("Val") , N_Raw.at(t.first) , Acceptance.at(t.first) , Efficiency.at(t.first));
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
  const double common   = ( 2.0 / (sum * sum) );
  const double error2Pl = std::pow( ( Err_Plus * N_Minus ) , 2.0 );
  const double error2Mi = std::pow( ( Err_Minus * N_Plus ) , 2.0 );
  return ( std::abs(common) * std::sqrt( error2Pl + error2Mi ) );
};


bool computeChargeAsymmetry( BinPentaMap& var , const VarBinMap& inputVar )
{
  //
  for (const auto& c : inputVar) {
    for (const auto& b : c.second) {
      // Check that everything is fine
      if (b.second.count("Pl")==0) { std::cout << "[ERROR] Plus charge is missing in bin [" << b.first.etabin().high() << " , " << b.first.etabin().low() << "]" << std::endl; return false; }
      if (b.second.count("Mi")==0) { std::cout << "[ERROR] Minus charge is missing in bin [" << b.first.etabin().high() << " , " << b.first.etabin().low() << "]" << std::endl; return false; }
      if (b.second.at("Mi").count("N_WToMu")==0 || b.second.at("Pl").count("N_WToMu")==0) {
        std::cout << "[ERROR] N_WToMu variable is missing in bin [" << b.first.etabin().high() << " , " << b.first.etabin().low() << "]" << std::endl; return false;
      }
      //
      auto&       oVar    = var[c.first][""]["Charge_Asymmetry"];
      const auto& iVar_Pl = b.second.at("Pl").at("N_WToMu");
      const auto& iVar_Mi = b.second.at("Mi").at("N_WToMu");
      // Compute charge asymmetry value
      oVar["Val"][b.first] = getChargeAsymmetryValue(iVar_Pl.at("Val"), iVar_Mi.at("Val"));
      // Compute charge asymmetry error
      const std::vector< std::string > errType = { "Err_Stat_High" , "Err_Stat_Low" , "Err_Syst_High" , "Err_Syst_Low" };
      for (const auto& t : errType) { oVar[t][b.first] = getChargeAsymmetryError(iVar_Pl.at("Val"), iVar_Mi.at("Val"), iVar_Pl.at(t), iVar_Mi.at(t)); }
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


bool computeForwardBackwardRatio( BinPentaMap& var , const VarBinMap& inputVar )
{
  //
  for (const auto& c : inputVar) {
    for (const auto& b : c.second) {
      if (b.first.etabin().low() < 0.0 ) continue; // Ignore the backward bins when looping
      // Build the backward bin
      const auto& binFw = b.first;
      const auto  binBw = anabin<0>(-1.0*binFw.etabin().high() , ( (binFw.etabin().low() == 0.0) ? 0.0 : -1.0*binFw.etabin().low() ));
      // Check that everything is fine
      if (c.second.count(binBw)==0) { std::cout << "[ERROR] Backward bin [" << binBw.etabin().low() << " , " << binBw.etabin().high() << "] is not defined" << std::endl; return false; }
      if (c.second.count(binFw)==0) { std::cout << "[ERROR] Forward bin [" << binBw.etabin().low() << " , " << binBw.etabin().high() << "] is not defined" << std::endl; return false; }
      if (c.second.at(binFw).count("Pl")==0) { std::cout << "[ERROR] Plus charge is missing in bin [" << binFw.etabin().low() << " , " << binFw.etabin().high() << "]" << std::endl; return false; }
      if (c.second.at(binFw).count("Mi")==0) { std::cout << "[ERROR] Minus charge is missing in bin [" << binFw.etabin().low() << " , " << binFw.etabin().high() << "]" << std::endl; return false; }
      if (c.second.at(binBw).count("Pl")==0) { std::cout << "[ERROR] Plus charge is missing in bin [" << binBw.etabin().low() << " , " << binBw.etabin().high() << "]" << std::endl; return false; }
      if (c.second.at(binBw).count("Mi")==0) { std::cout << "[ERROR] Minus charge is missing in bin [" << binBw.etabin().low() << " , " << binBw.etabin().high() << "]" << std::endl; return false; }
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
      // Compute forward-backward ratio value
      oVar_Inc["Val"][b.first] = getForwardBackwardRatioValue( ( iVar_FwPl.at("Val") + iVar_FwMi.at("Val") ) , ( iVar_BwPl.at("Val") + iVar_BwMi.at("Val") ) );
      oVar_Pl["Val"][b.first]  = getForwardBackwardRatioValue( iVar_FwPl.at("Val") , iVar_BwPl.at("Val") );
      oVar_Mi["Val"][b.first]  = getForwardBackwardRatioValue( iVar_FwMi.at("Val") , iVar_BwMi.at("Val") );
      // Compute forward-backward ratio error
      const std::vector< std::string > errType = { "Err_Stat_High" , "Err_Stat_Low" , "Err_Syst_High" , "Err_Syst_Low" };
      for (const auto& t : errType) {
        oVar_Inc[t][b.first] = getForwardBackwardRatioError(
                                                            iVar_FwPl.at("Val"), iVar_BwPl.at("Val"), iVar_FwMi.at("Val"), iVar_BwMi.at("Val"), 
                                                            iVar_FwPl.at(t)    , iVar_BwPl.at(t)    , iVar_FwMi.at(t)    , iVar_BwMi.at(t)
                                                            );
      }
      for (const auto& t : errType) { oVar_Pl[t][b.first] = getForwardBackwardRatioError(iVar_FwPl.at("Val"), iVar_BwPl.at("Val"), iVar_FwPl.at(t), iVar_BwPl.at(t)); }
      for (const auto& t : errType) { oVar_Mi[t][b.first] = getForwardBackwardRatioError(iVar_FwMi.at("Val"), iVar_BwMi.at("Val"), iVar_FwMi.at(t), iVar_BwMi.at(t)); }
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


bool computeCrossSection( BinPentaMap& var , const VarBinMap& inputVar )
{
  //
  for (const auto& c : inputVar) {
    for (const auto& b : c.second) {
      // Check that everything is fine
      if (b.second.count("Pl")==0) { std::cout << "[ERROR] Plus charge is missing in bin [" << b.first.etabin().high() << " , " << b.first.etabin().low() << "]" << std::endl; return false; }
      if (b.second.count("Mi")==0) { std::cout << "[ERROR] Minus charge is missing in bin [" << b.first.etabin().high() << " , " << b.first.etabin().low() << "]" << std::endl; return false; }
      if (b.second.at("Mi").count("N_WToMu")==0 || b.second.at("Pl").count("N_WToMu")==0) {
        std::cout << "[ERROR] N_WToMu variable is missing in bin [" << b.first.etabin().high() << " , " << b.first.etabin().low() << "]" << std::endl; return false;
      }
      //
      auto&       oVar_Pl = var[c.first]["Pl"]["Cross_Section"];
      auto&       oVar_Mi = var[c.first]["Mi"]["Cross_Section"];
      const auto& iVar_Pl = b.second.at("Pl").at("N_WToMu");
      const auto& iVar_Mi = b.second.at("Mi").at("N_WToMu");
      //
      double Luminosity = -1.0;
      if (c.first=="PA" ) { Luminosity = ( PA::LUMI::Data_pPb + PA::LUMI::Data_Pbp ); }
      if (c.first=="pPb") { Luminosity = PA::LUMI::Data_pPb; }
      if (c.first=="Pbp") { Luminosity = PA::LUMI::Data_Pbp; }
      std::map< std::string , double > Err_Luminosity;
      Err_Luminosity["Err_Stat_High"] = ( 0.05 * Luminosity ); // 5% relative error
      Err_Luminosity["Err_Stat_Low" ] = ( 0.05 * Luminosity ); // 5% relative error
      Err_Luminosity["Err_Syst_High"] = 0.0; // Not assigned
      Err_Luminosity["Err_Syst_Low" ] = 0.0; // Not assigned
      //
      const double BinWidth = ( b.first.etabin().high() - b.first.etabin().low() );
      //
      // Compute the cross-section value
      oVar_Pl["Val"][b.first] = getCrossSectionValue(iVar_Pl.at("Val"), Luminosity, BinWidth);
      oVar_Mi["Val"][b.first] = getCrossSectionValue(iVar_Mi.at("Val"), Luminosity, BinWidth);
      // Compute the cross-section error
      const std::vector< std::string > errType = { "Err_Stat_High" , "Err_Stat_Low" , "Err_Syst_High" , "Err_Syst_Low" };
      for (const auto& t : errType) { oVar_Pl[t][b.first] = getCrossSectionError(iVar_Pl.at("Val"), Luminosity, BinWidth, iVar_Pl.at(t), Err_Luminosity.at(t)); }
      for (const auto& t : errType) { oVar_Mi[t][b.first] = getCrossSectionError(iVar_Mi.at("Val"), Luminosity, BinWidth, iVar_Mi.at(t), Err_Luminosity.at(t)); }
    }
  }
  return true;
};


std::string formatResultVarName(const std::string varName)
{
  std::string label = "";
  if (varName == "Charge_Asymmetry"      ) { label = "( N^{+} - N^{-} ) / ( N^{+} + N^{-} )"; }
  if (varName == "ForwardBackward_Ratio" ) { label = "R_{FB}"; }
  if (varName == "Cross_Section"         ) { label = "B #times d#sigma/d#eta (nb)"; }
  return label;
};


void iniResultsGraph(GraphQuadMap& graphMap, const BinPentaMap& var)
{
  //
  std::cout << "[INFO] Initializing the output graphs" << std::endl;
  const std::vector< std::string > graphType = { "Err_Tot" , "Err_Stat" , "Err_Syst" };
  //
  // Initialize the graphs
  for (const auto& c : var) {
    for (const auto& ch : c.second) {
      for (const auto& v : ch.second) {
        const unsigned int nBins = v.second.at("Val").size();
        for (const auto& t : graphType) {
          auto& graph = graphMap[c.first][ch.first][v.first][t];
          graph.Set(nBins);
          //
          // Set Graph Name
          const std::string name = Form("gr_WToMu%s_%s_%s_%s", ch.first.c_str(), c.first.c_str(), v.first.c_str(), t.c_str());
          graph.SetName(name.c_str());
        }
      }
    }
  }
};


bool fillResultsGraph(GraphQuadMap& graphMap, const BinPentaMap& var)
{
  //
  std::cout << "[INFO] Filling the output graphs" << std::endl;
  const std::vector< std::string > graphType = { "Err_Tot" , "Err_Stat" , "Err_Syst" };
  //
  for (const auto& c : var) {
    for (const auto& ch : c.second) {
      for (const auto& v : ch.second) {
        if (v.second.count("Val")==0) { std::cout << "[ERROR] Value is missing for " << v.first << std::endl; return false; }
        if ( (v.second.count("Err_Stat_High")==0) || (v.second.count("Err_Stat_Low")==0) ) { std::cout << "[ERROR] Statisticial errors are missing for " << v.first << std::endl; return false; }
        if ( (v.second.count("Err_Syst_High")==0) || (v.second.count("Err_Syst_Low")==0) ) { std::cout << "[ERROR] Systematic errors are missing for "   << v.first << std::endl; return false; }
        for (const auto& t : graphType) {
          // Get the graph
          auto& graph = graphMap.at(c.first).at(ch.first).at(v.first).at(t);
          //
          // Determine index of bins
          std::map< anabin<0> , unsigned int > binIdx;
          unsigned int iBin = 0; for (const auto& b : v.second.at("Val")) { binIdx[b.first] = iBin; iBin++; }
          //
          for (const auto& b : v.second.at("Val")) {
            //
            // Extract the parameters needed for each axis
            //
            // X Value
            const double X = ( (b.first.etabin().high() + b.first.etabin().low()) / 2.0 ); // Mean value of eta bin
            // X Error
            const double Err_X      = ( (b.first.etabin().high() - b.first.etabin().low()) / 2.0 ); // Width of eta bin
            double Err_X_High = Err_X;
            double Err_X_Low  = Err_X;
            if (t=="Err_Stat") { Err_X_High *= 0.4; Err_X_Low *= 0.4; }
            if (t=="Err_Syst") { Err_X_High *= 0.6; Err_X_Low *= 0.6; }
            // Y Value
            const double Y = v.second.at("Val").at(b.first);
            // Y Error
            double Err_Y_High = -1.0;
            double Err_Y_Low  = -1.0;
            if (t=="Err_Tot") {
              Err_Y_High = std::sqrt( std::pow( v.second.at("Err_Stat_High").at(b.first) , 2.0 ) + std::pow( v.second.at("Err_Syst_High").at(b.first) , 2.0 ) );
              Err_Y_Low  = std::sqrt( std::pow( v.second.at("Err_Stat_Low" ).at(b.first) , 2.0 ) + std::pow( v.second.at("Err_Syst_Low" ).at(b.first) , 2.0 ) );
            }
            else {
              Err_Y_High = v.second.at(Form("%s_High", t.c_str())).at(b.first);
              Err_Y_Low  = v.second.at(Form("%s_Low" , t.c_str())).at(b.first);
            }
            //
            // Fill the graphs
            //
            const unsigned int iBin = binIdx.at(b.first);
            graph.SetPoint(iBin, X, Y);
            graph.SetPointError(iBin, Err_X_Low, Err_X_High, Err_Y_Low, Err_Y_High);
          }
        }
      }
    }
  }
  //
  return true;
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


void formatResultsGraph(TGraphAsymmErrors& graph, const std::string& var, const std::string& chg, const bool& useEtaCM, const bool& isEffCorr)
{
  //
  // Set the Axis Titles
  std::string xLabel = "#mu"; if (chg == "Pl") { xLabel += "^{+}"; }; if (chg == "Mi") { xLabel += "^{-}"; }; xLabel += " #eta";
  if (useEtaCM) { xLabel += "_{CM}"; }
  else { xLabel += "_{LAB}"; }
  std::string yLabel = formatResultVarName(var);
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
  if (useEtaCM) {
    if ( var == "Charge_Asymmetry"      ) { graph.GetXaxis()->SetLimits(-2.0, 2.0); }
    if ( var == "ForwardBackward_Ratio" ) { graph.GetXaxis()->SetLimits( 0.0, 2.0); }
    if ( var == "Cross_Section"         ) { graph.GetXaxis()->SetLimits(-2.0, 2.0); }
  }
  else {
    if ( var == "Charge_Asymmetry"      ) { graph.GetXaxis()->SetLimits(-2.5, 2.5); }
    if ( var == "ForwardBackward_Ratio" ) { graph.GetXaxis()->SetLimits( 0.0, 2.5); }
    if ( var == "Cross_Section"         ) { graph.GetXaxis()->SetLimits(-2.5, 2.5); }
  }
  // Y-axis
  graph.GetYaxis()->CenterTitle(kFALSE);
  graph.GetYaxis()->SetTitleOffset(1.45);
  graph.GetYaxis()->SetTitleSize(0.050);
  if ( var == "Charge_Asymmetry" ) { graph.GetYaxis()->SetTitleSize(0.040); }
  graph.GetYaxis()->SetLabelSize(0.035);
  if ( var == "Charge_Asymmetry"      ) { graph.GetYaxis()->SetRangeUser(   0.0,   0.5); }
  if ( var == "ForwardBackward_Ratio" ) { graph.GetYaxis()->SetRangeUser(   0.6,   1.5); }
  if ( var == "Cross_Section"         ) { graph.GetYaxis()->SetRangeUser(  50.0, 350.0); }
};


void drawGraph( GraphQuadMap& graphMap , const std::string& outDir , const bool useEtaCM = true , const std::string accType = "MC" , const std::string effType = "TnP" )
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
        //
        const std::string col = c.first;
        const std::string chg = ( (ch.first!="") ? ch.first : "Inc" );
        const std::string var = v.first;
        auto& graph = v.second;
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
        //
        // Create Legend
        TLegend leg(0.2, 0.71, 0.4, 0.84);
        formatLegendEntry(*leg.AddEntry(&graph.at("Err_Tot") , "Data", "pe"));
        formatLegendEntry(*leg.AddEntry(&graph.at("Err_Stat"), "Statistical Uncertainty", "f"));
        formatLegendEntry(*leg.AddEntry(&graph.at("Err_Syst"), "Systematic Uncertainty", "f"));
        //
        // Format the graphs
        const bool isEffCorr = (accType!="" || effType!="");
        for (auto& gr : graph) { formatResultsGraph(gr.second, var, chg, useEtaCM, isEffCorr); }
        graph.at("Err_Tot").SetMarkerColor(kBlack);
        graph.at("Err_Stat").SetFillColor(kOrange);
        graph.at("Err_Syst").SetFillColor(kGreen+3);
        //
        // Draw the graphs
        graph.at("Err_Syst").Draw("a2");
        graph.at("Err_Stat").Draw("same2");
        graph.at("Err_Tot").Draw("samep");
        // Draw the Line
        double etaMin = -2.5; if (useEtaCM) { etaMin = -2.0; }
        double etaMax = 2.5; if (useEtaCM) { etaMax = 2.0; }
        TLine line_FB(0.0, 1.0, etaMax, 1.0); line_FB.SetLineStyle(2);
        if (var=="ForwardBackward_Ratio") { line_FB.Draw("same"); }
        TLine line_CA(etaMin, 0.0, etaMax, 0.0); line_CA.SetLineStyle(2);
        if (var=="Charge_Asymmetry") { line_CA.Draw("same"); }
        // Draw the Legend
        leg.Draw("same");
        // Update
        c.Modified(); c.Update();
        // Draw the text
        for (const auto& s: textToPrint) { tex.DrawLatex(0.22, 0.86-dy, s.c_str()); dy+=0.04; }
        // Update
        c.Modified(); c.Update(); // Pure paranoia
        //
        // set the CMS style
        int option = 117;
        if (col.find("pPb")!=std::string::npos) option = 115;
        if (col.find("Pbp")!=std::string::npos) option = 116;
        CMS_lumi(&c, option, 33, "");
        // Update
        c.Modified(); c.Update(); // Pure paranoia
        //
        // Create Output Directory
        const std::string plotDir = outDir + "/Output/" + col+"/" + var;
        makeDir(plotDir + "/png/");
        makeDir(plotDir + "/pdf/");
        makeDir(plotDir + "/root/");
        //
        std::string label = "";
        if (accType==""   && effType==""   ) { label = "RAW";          }
        if (accType=="MC" && effType==""   ) { label = "AccMC";        }
        if (accType==""   && effType=="MC" ) { label = "EffMC";        }
        if (accType=="MC" && effType=="MC" ) { label = "AccMC_EffMC";  }
        if (accType==""   && effType=="TnP") { label = "EffTnP";       }
        if (accType=="MC" && effType=="TnP") { label = "AccMC_EffTnP"; }
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


#endif // ifndef resultUtils_h
  

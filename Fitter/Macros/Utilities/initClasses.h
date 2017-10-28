#ifndef initClasses_h
#define initClasses_h

#include "TSystem.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TSystemFile.h"
#include "TSystemDirectory.h"

#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TH1.h"
#include "TF1.h"
#include "TLine.h"
#include "TPad.h"
#include "TObjString.h"

#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooCategory.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooHistPdf.h"
#include "RooClassFactory.h"
#include "RooFitResult.h"
#include "RooArgList.h"
#include "RooPlot.h"
#include "RooHist.h"

#include "../../../Utilities/CMS/tdrstyle.C"
#include "../../../Utilities/CMS/CMS_lumi.C"
#include "../../../Utilities/EVENTUTILS.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <bitset>
#include <algorithm>


typedef std::vector< std::string > StringVector;
typedef std::map< std::string , RooWorkspace                    > RooWorkspaceMap_t;
typedef std::map< std::string , std::vector< std::string >      > StringVectorMap_t;
typedef std::map< std::string , StringVectorMap_t               > StringVectorDiMap_t;
typedef std::map< std::string , StringVectorDiMap_t             > StringVectorTriMap_t;
typedef std::map< std::string , std::map< std::string , double >> DoubleMapMap_t;
typedef std::map< std::string , std::string                     > StringMap_t;
typedef std::map< std::string , int                             > IntMap_t;
typedef std::map< std::string , IntMap_t                        > IntMapMap_t;
typedef std::map< std::string , bool                            > BoolMap_t;
typedef std::map< std::string , StringMap_t                     > StrMapMap_t;
typedef std::map< std::string , std::map< std::string , std::map< std::string , int > > > ModelMap;
typedef std::map< std::string , std::map< std::string , std::vector< int > > >  IntVecMapMap_t;

// MET Model
enum class METModel 
{
    InvalidModel,
    CutAndCount,
    Template,
    MultiJetBkg,
    ModifiedRayleigh,
    Size
};
const IntMap_t METModelDictionary = {
  {"InvalidModel",     int(METModel::InvalidModel)},
  {"CutAndCount",      int(METModel::CutAndCount)},
  {"Template",         int(METModel::Template)},
  {"MultiJetBkg",      int(METModel::MultiJetBkg)},
  {"ModifiedRayleigh", int(METModel::ModifiedRayleigh)}
};
// Add the Models to the Main Dictionary
const IntMapMap_t ModelDictionary = {
  { "MET" , METModelDictionary }
};

const StringMap_t varEWQLabel = {
  { "MET"        , "|#slash{E}_{T}|"   },
  { "Muon_Iso"   , "Iso_{PFR03}^{#mu}" },
  { "Muon_Pt"    , "p_{T}^{#mu}"       },
  { "Muon_Eta"   , "#eta_{LAB}^{#mu}"  },
  { "Muon_EtaCM" , "#eta_{CM}^{#mu}"   },
  { "Muon_MT"    , "M_{T}^{#mu}"       },
  { "Centrality" , "Cent."             }
};

// Global Info Structure (wrapper to carry information around)
typedef struct GlobalInfo {
  DoubleMapMap_t    Var;
  StringMap_t       Par;
  IntMap_t          Int;
  StringVectorMap_t StrV;
  BoolMap_t         Flag;
  void              Clear() { this->Var.clear(); this->Par.clear(); this->Int.clear(); this->StrV.clear(); this->Flag.clear(); }
  GlobalInfo() {}
  GlobalInfo(const GlobalInfo &ref, bool keep = false) {
    this->Copy(ref, keep);
  }
  ~GlobalInfo() {
    this->Clear();
  }
  void Copy(const DoubleMapMap_t &ref, bool keep = true) {
    if (!keep) this->Var.clear();
    for (const auto& var : ref) {
      for (const auto& ele : var.second) {
        this->Var[var.first][ele.first] = ele.second;
      }
    }
  }
  void Copy(const StringMap_t &ref, bool keep = true) {
    if (!keep) this->Par.clear();
    for (const auto& par : ref) {
      this->Par[par.first] = par.second;
    }
  }
  void Copy(const IntMap_t &ref, bool keep = true) {
    if (!keep) this->Int.clear();
    for (const auto& i : ref) {
      this->Int[i.first] = i.second;
    }
  }
  void Copy(const StringVectorMap_t &ref, bool keep = true) {
    if (!keep) this->StrV.clear();
    for (const auto& i : ref) {
      this->StrV[i.first] = i.second;
    }
  }
  void Copy(const BoolMap_t &ref, bool keep = true) {
    if (!keep) this->Flag.clear();
    for (const auto& flag : ref) {
      this->Flag[flag.first] = flag.second;
    }
  }
  void Copy(const GlobalInfo &ref, bool keep = true) {
    this->Copy(ref.Var, keep);
    this->Copy(ref.Par, keep);
    this->Copy(ref.Int, keep);
    this->Copy(ref.StrV, keep);
    this->Copy(ref.Flag, keep);
  }
  bool operator == (const DoubleMapMap_t &ref) const 
  {
    if (ref.size() != this->Var.size()) return false;
    for (const auto& var : this->Var) {
      if (ref.count(var.first)==0 || ref.at(var.first).count("Min")==0 || ref.at(var.first).count("Max")==0) return false;
      if (var.second.at("Min") != ref.at(var.first).at("Min")) return false;
      if (var.second.at("Max") != ref.at(var.first).at("Max")) return false;
    }
    return true;
  }
  bool operator == (const StringMap_t &ref) const 
  {
    if (ref.size() != this->Par.size()) return false;
    for (const auto& par : this->Par) {
      if (ref.count(par.first)==0) return false;
      if (par.second != ref.at(par.first)) return false;
    }
    return true;
  }
  bool operator == (const IntMap_t &ref) const 
  {
    if (ref.size() != this->Int.size()) return false;
    for (const auto& i : this->Int) {
      if (ref.count(i.first)==0) return false;
      if (i.second != ref.at(i.first)) return false;
    }
    return true;
  }
  bool operator == (const StringVectorMap_t &ref) const 
  {
    if (ref.size() != this->StrV.size()) return false;
    for (const auto& i : this->StrV) {
      if (ref.count(i.first)==0) return false;
      if (i.second != ref.at(i.first)) return false;
    }
    return true;
  }
  bool operator == (const BoolMap_t &ref) const 
  {
    if (ref.size() != this->Flag.size()) return false;
    for (const auto& flag : this->Flag) {
      if (ref.count(flag.first)==0) return false;
      if (flag.second != ref.at(flag.first)) return false;
    }
    return true;
  }
  bool operator != ( const DoubleMapMap_t     &ref ) const { if (*this == ref) { return false; } return true; }
  bool operator != ( const StringMap_t        &ref ) const { if (*this == ref) { return false; } return true; }
  bool operator != ( const IntMap_t           &ref ) const { if (*this == ref) { return false; } return true; }
  bool operator != ( const StringVectorMap_t  &ref ) const { if (*this == ref) { return false; } return true; }
  bool operator != ( const BoolMap_t          &ref ) const { if (*this == ref) { return false; } return true; }
  bool operator == ( const GlobalInfo         &ref) const 
  {
    if ( *this != ref.Var  ) return false;
    if ( *this != ref.Par  ) return false;
    if ( *this != ref.Int  ) return false;
    if ( *this != ref.StrV ) return false;
    if ( *this != ref.Flag ) return false;
    return true;
  }
} GlobalInfo;


bool splitString(std::string input, std::string delimiter, std::vector< std::string >& output)
{
  // remove spaces from input string 
  input.erase(std::remove(input.begin(), input.end(), ' '), input.end());
  // proceed to parse input string
  while(input!="") {
    std::string d = input.substr(0, input.find(delimiter));
    output.push_back(d);
    if(input.find(delimiter.c_str())!= std::string::npos){ input.erase(0, input.find(delimiter) + delimiter.length()); }
    else { input = ""; }
  }
  return true;
};

void setFixedVarsToContantVars(RooWorkspace& ws)
{
  RooArgSet listVar = ws.allVars();
  std::unique_ptr<TIterator> parIt = std::unique_ptr<TIterator>(listVar.createIterator());
  for (RooRealVar* it = (RooRealVar*)parIt->Next(); it!=NULL; it = (RooRealVar*)parIt->Next() ) {
    if ( it->getMin()==it->getMax() && !it->isConstant() ) {
      std::cout << "[INFO] Setting " << it->GetName() << " constant!" << std::endl;
      it->setConstant(kTRUE);
    }
  }
};

bool saveWorkSpace(const RooWorkspace& ws, const string& outputDir, const string& fileName)
{
  // Save the workspace
  gSystem->mkdir(outputDir.c_str(), kTRUE);
  std::cout << (outputDir+fileName) << std::endl;
  std::unique_ptr<TFile> file = std::unique_ptr<TFile>(new TFile((outputDir+fileName).c_str(), "RECREATE"));
  if (!file) {
    file->Close();
    std::cout << "[ERROR] Output root file with fit results could not be created!" << std::endl; return false; 
  } else {
    file->cd();    
    ws.Write("workspace"); 
    file->Write(); file->Close();
  }
  return true;
};

bool compareSnapshots(const RooArgSet& pars1, const RooArgSet& pars2) {
  std::unique_ptr<TIterator> parIt = std::unique_ptr<TIterator>(pars1.createIterator());
  for (RooRealVar* it = (RooRealVar*)parIt->Next(); it!=NULL; it = (RooRealVar*)parIt->Next() ) {
    double val = pars2.getRealValue(it->GetName(),-1e99);
    const std::string name = it->GetName();
    if ( (name=="MET") || (name=="Muon_Pt") || (name=="Muon_Eta") || (name=="Muon_Iso") || (name=="Muon_MT") || (name=="Centrality") || (name=="Event_Type") ) continue;
    if (val==-1e99) return false;          // the parameter was not found!
    if (val != it->getVal()) return false;  // the parameter was found, but with a different value!
    if ( ((RooRealVar&)pars2[it->GetName()]).getMin() != it->getMin() ) return false;  // the parameter has different lower limit
    if ( ((RooRealVar&)pars2[it->GetName()]).getMax() != it->getMax() ) return false;  // the parameter has different upper limit
  }
  return true;
};

bool isFitAlreadyFound(const RooArgSet& newpars, const string& fileName, const string& pdfName) 
{
  std::cout << "[INFO] Checking if fit was already done!" << std::endl;
  if (gSystem->AccessPathName(fileName.c_str())) {
    std::cout << "[INFO] FileName: " << fileName << " was not found" << std::endl;
    return false; // File was not found
  }
  std::unique_ptr<TFile> file = std::unique_ptr<TFile>(new TFile(fileName.c_str()));
  if (!file) return false;
  RooWorkspace* ws = (RooWorkspace*) file->Get("workspace");
  if (!ws) {
    std::cout << "[INFO] Workspace was not found" << std::endl;
    file->Close();
    return false;
  }
  const RooArgSet* params = ( (ws->getSnapshot("initialParameters")) ? ws->getSnapshot("initialParameters") :  ws->getSnapshot(Form("%s_parIni", pdfName.c_str())) );
  if (!params) {
    std::cout << "[INFO] Snapshot of initial parameters was not found!" << std::endl;
    file->Close();
    return false;
  }
  bool result = compareSnapshots(newpars, *params);
  file->Close();
  return result;
};

void stringReplace(std::string& txt, const std::string& from, const std::string& to)
{
  std::string::size_type n = 0;
  while ( ( n = txt.find( from, n ) ) != std::string::npos )
    {
      txt.replace( n, from.size(), to );
      n += to.size();
    }
};

std::string formatCut(const std::string& cut, const StringMap_t& map = StringMap_t())
{
  std::string str = cut;
  str.erase(std::remove(str.begin(), str.end(), ' '), str.end());
  stringReplace( str, "<=MET", " GeV/c<=MET" ); stringReplace( str, "<MET", " GeV/c<MET" );
  stringReplace( str, "<=Muon_MT", " GeV/c^{2}<=Muon_MT" ); stringReplace( str, "<Muon_MT", " GeV/c^{2}<Muon_MT" );
  for (auto& elem : map) { stringReplace( str, elem.first, elem.second ); }
  stringReplace( str, "Pl", "^{+}+x" ); stringReplace( str, "Mi", "^{-}+x" );
  stringReplace( str, "Mu", "#mu" ); stringReplace( str, "Tau", "#tau" ); stringReplace( str, "DY", "Z/#gamma*" ); stringReplace( str, "TTbar", "t#bar{t}" );
  stringReplace( str, "To", "#rightarrow" ); stringReplace( str, "&&", " & " ); stringReplace( str, "(", "" ); stringReplace( str, ")", "" );
  stringReplace( str, "<=", " #leq " ); stringReplace( str, "<", " < " ); stringReplace( str, ">=", " #geq " ); stringReplace( str, ">", " > " );
  return str;
};


RooRealVar getVar(const RooArgSet& set, const std::string varName)
{
  if (set.find(varName.c_str()) != NULL) {
    return *(RooRealVar*)set.find(varName.c_str());
  }
  std::cout << "[WARNING] Variable " << varName << " was not found in the set, return empty RooRealVar" << std::endl;
  return RooRealVar();
};


void getPoint(const RooHist& rH, const uint& i, double& x, double& y, double& exl, double& exh, double& eyl, double& eyh)
{
  //
  rH.GetPoint(i, x, y);
  //
  eyl = rH.GetErrorYlow(i);
  eyh = rH.GetErrorYhigh(i);
  //
  exl = rH.GetErrorXlow(i);
  exh = rH.GetErrorXhigh(i);
  if (exl<=0.0 ) { exl = rH.GetErrorX(i); }
  if (exh<=0.0 ) { exh = rH.GetErrorX(i); }
  if (exl<=0.0 ) { exl = 0.5*rH.getNominalBinWidth(); }
  if (exh<=0.0 ) { exh = 0.5*rH.getNominalBinWidth(); }
};


bool rooPlotToTH1(TH1D& hData, TH1D& hFit, const RooPlot& frame, const bool useAverage = true)
{
  // Find curve object
  auto rFit = (RooCurve*) frame.findObject(0, RooCurve::Class());
  if (!rFit) { std::cout << "[ERROR] The latest RooCurve was not found" << std::endl; return false; }
  // Find histogram object
  auto rData = (RooHist*) frame.findObject(0, RooHist::Class());
  if (!rData) { std::cout << "[ERROR] The latest RooHist was not found" << std::endl; return false; }
  // Determine range of curve
  double xstart, xstop, yDummy;
  rFit->GetPoint(0, xstart, yDummy);
  rFit->GetPoint((rFit->GetN()-1), xstop, yDummy);
  // Get Binning
  std::vector<double> binV;
  for (int i = 0; i < rData->GetN(); i++) {
    double x, y; rData->GetPoint(i, x, y);
    // Only consider bins inside curve range
    if (x<xstart || x>xstop) continue;
    const double binW = rData->getNominalBinWidth();
    binV.push_back(x - (binW*0.5));
    if (i==(rData->GetN()-1)) { binV.push_back(x + (binW*0.5)); }
  }
  const uint nBin = (binV.size()-1);
  Double_t bin[nBin+1];
  for (uint i = 0; i < binV.size(); i++) { bin[i] = binV[i]; }
  //
  hData.Reset(); hData = TH1D(Form("hData_%s", rData->GetName()), rData->GetTitle(), nBin, bin); hData.Sumw2();
  hFit.Reset();  hFit  = TH1D(Form("hFit_%s" , rFit->GetName()) , rFit->GetTitle() , nBin, bin); hFit.Sumw2();
  // Set Histogram entries
  for (uint i = 0; i < nBin; i++) {
    double x, dataVal, exl, exh, eyl, eyh;
    getPoint(*rData, i, x, dataVal, exl, exh, eyl, eyh);
    double fitVal = 0.0;
    if (useAverage) { fitVal = rFit->average(x-exl, x+exh); }
    else            { fitVal = rFit->interpolate(x);        }
    hData.SetBinContent((i+1), dataVal);
    hData.SetBinError((i+1), std::sqrt((eyl*eyl + eyh*eyh)/2.0));
    hFit.SetBinContent((i+1), fitVal);
    hFit.SetBinError((i+1), std::sqrt(fitVal));
  }
  return true;
};


bool makePullHist(RooHist& pHist, const RooPlot& frame, const std::string histname, const std::string curvename, const bool useAverage)
{
  // Find curve object
  auto rFit = (RooCurve*) frame.findObject(((curvename=="") ? 0 : curvename.c_str()), RooCurve::Class());
  if (!rFit) { std::cout << "[ERROR] makePullHist(" << curvename << ") cannot find curve" << std::endl; return false; }
  // Find histogram object
  auto rData = (RooHist*) frame.findObject(((histname=="") ? 0 : histname.c_str()), RooHist::Class());
  if (!rData) { std::cout << "[ERROR] makePullHist(" << histname  << ") cannot find histogram" << std::endl; return false; }
  // Determine range of curve
  double xstart, xstop, yDummy;
  rFit->GetPoint(0, xstart, yDummy);
  rFit->GetPoint((rFit->GetN()-1), xstop, yDummy);
  // Add histograms, calculate Poisson confidence interval on sum value
  for (int i = 0; i < rData->GetN(); i++) {
    // Get Data Value and Error
    double x, dataVal, exl, exh, eyl, eyh;
    getPoint(*rData, i, x, dataVal, exl, exh, eyl, eyh);
    // Only calculate pull for bins inside curve range
    if (x<xstart || x>xstop) continue ;
    // Get Fit Value
    double fitVal = 0.0;
    if (useAverage) { fitVal = rFit->average(x-exl, x+exh); }
    else            { fitVal = rFit->interpolate(x);        }
    // Compute the Residual
    double y = (dataVal - fitVal);
    // Get the Norm factor
    const double norm = ( (y>0) ? eyl : eyh );
    // Set Values
    if (norm>0.0) {
      y   /= norm;
      eyl /= norm;
      eyh /= norm;
      pHist.addBinWithXYError(x, y, exl, exh, eyl, eyh);
    }
    else {
      pHist.addBinWithXYError(x, 0, exl, exh, 0, 0);
    }
  }
  return true;
};


double getErrorHi(const RooRealVar& var)
{
  if (var.getErrorLo()==0.0 && var.getErrorHi()==0.0) { return var.getError(); }
  return std::abs(var.getErrorHi());
};


double getErrorLo(const RooRealVar& var)
{
  if (var.getErrorLo()==0.0 && var.getErrorHi()==0.0) { return var.getError(); }
  return std::abs(var.getErrorLo());
};  


bool isParAtLimit(const RooRealVar& var)
{
  if (
      ( (std::abs(var.getValV() - var.getMin())/getErrorLo(var)) <= 3.0 ) ||
      ( (std::abs(var.getValV() - var.getMax())/getErrorHi(var)) <= 3.0 )
      )
    { return true; }
  return false;
};


#endif // #ifndef initClasses_h

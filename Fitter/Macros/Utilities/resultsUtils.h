#ifndef Utilities_resultUtils_h
#define Utilities_resultUtils_h

#include "bin.h"

#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooAbsData.h"
#include "RooWorkspace.h"
#include "TString.h"
#include "TFile.h"
#include "TSystemFile.h"
#include "TSystemDirectory.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "RooFitResult.h"
#include "TFitResult.h"
#include "TTree.h"

#include <string>
#include <vector>
#include <iostream>

///////////////
// CONSTANTS //
///////////////

bool binok(std::vector< anabin<0> > thecats, const std::string& xaxis, anabin<0> &tocheck, bool override=true);
bool binok(anabin<0> thecat, const std::string& xaxis, anabin<0> &tocheck, bool override=true);
TString treeFileName(const char* workDirName, const char* DSTag, const char* prependPath, const char* fitType);
RooRealVar* poiFromWS(const RooWorkspace& ws, const char* token, const char* thepoiname);
std::vector<TString> fileList(const char* input, const char* token, const char* DSTag, const char* prependPath, const char* fitType);
void resultsEWQ2tree(const char* workDirName = "Test", const char* DSTag = "DATA",  const char* prependPath = "",  const char* fitType = "MET", const char* thePoiNames = "all");

#include "../../resultsEWQ2tree.C"

TString treeFileName(const char* workDirName, const char* DSTag, const char* prependPath, const char* fitType)
{
  TString outputFileName("");
  if (!strcmp(fitType,"")) outputFileName = Form("Output/%s/result/%s/tree_allvars.root",workDirName,DSTag);
  else outputFileName = Form("Output/%s/%s/%s/result/tree_allvars.root",workDirName,fitType,DSTag);
  if ( strcmp(prependPath,"") ) outputFileName.Prepend(Form("%s/",prependPath));
  return outputFileName;
};

RooRealVar* poiFromWS(const RooWorkspace& ws, const char* token, const char* thepoiname)
{
  TString poiname_and_token = TString(thepoiname) + TString(token);
  RooRealVar* ans = ((RooRealVar*) ws.var(poiname_and_token));
  if (!ans) { return NULL; }
  return ans;
};

RooAbsPdf* pdfFromWS(const RooWorkspace& ws, const char* token, const char* thepdfname)
{
  TString pdfname_and_token = TString(thepdfname) + TString(token);
  RooAbsPdf* ans = ((RooAbsPdf*) ws.pdf(pdfname_and_token));
  if (!ans) { return NULL; }
  return ans;
};

RooAbsData* dataFromWS(const RooWorkspace& ws, const char* token, const char* thedataname)
{
  TString dataname_and_token = TString(thedataname) + TString(token);
  RooAbsData* ans = ((RooAbsData*) ws.data(dataname_and_token));
  if (!ans) { return NULL; }
  return ans;
};

RooRealVar poiFromFile(const char* thepoiname, const char* filename, const char* token)
{
  TFile *f = TFile::Open(filename);
  if (!f) { std::cout << "[Error] File: " << filename << " does not exist!" << std::endl; return RooRealVar(); }
  if (!f->IsOpen()) { std::cout << "[Error] File: " << filename << " fail to open!" << std::endl; delete f; return RooRealVar(); }
  RooWorkspace* ws = (RooWorkspace*) f->Get("workspace");
  if (!ws) { std::cout << "[Error] File: " << filename << " does not have the workspace!" << std::endl; f->Close(); delete f; return RooRealVar(); }
  RooRealVar* ans = poiFromWS(*ws, token, thepoiname); 
  if (!ans) { std::cout << "[Error] File: " << filename << " does not have the variable: " << thepoiname << " !" << std::endl; f->Close(); delete f; return RooRealVar(); }
  TString poiname_and_token = TString(thepoiname) + TString(token);
  RooRealVar ansc = RooRealVar(*ans, poiname_and_token + Form("_from_%s",filename));
  delete ws; f->Close(); delete f;
  return ansc;
};

std::vector<TString> fileList(const char* input, const char* token, const char* DSTag, const char* prependPath, const char* fitType)
{
  std::vector<TString> ans;
  TString basedir("");
  if (!strcmp(fitType,"")) basedir = Form("./Output/%s/%s/result/",input,DSTag);
  else basedir = Form("./Output/%s/%s/%s/result/",input,fitType,DSTag);
  if ( strcmp(prependPath,"") ) basedir.Prepend(Form("%s/",prependPath));
  TSystemDirectory dir(input,basedir);
  TList *files = dir.GetListOfFiles();
  if (files) {
    TIter next(files);
    TSystemFile *file;
    TString fname;
    while ((file=(TSystemFile*)next())) {
      fname = file->GetName();
      if (fname.EndsWith(".root") && fname.Index("FIT") != kNPOS
          && (TString(token) == "" || fname.Index(token) != kNPOS)) {
        ans.push_back(basedir+fname);
      }
    }
  }
  delete files;
  return ans;
};

RooRealVar ratioVar( const RooRealVar& num, const RooRealVar& den, const bool& usedenerror )
{
  double n = num.getVal();
  double d = den.getVal();
  double dn = num.getError();
  double dd = den.getError();
  double r = d!=0 ? n/d : 0;
  double dr = n!=0 && d!=0 ? fabs(r * sqrt(pow(dn/n,2) + pow(dd/d,2))) : 0;
  if (!usedenerror && n!=0) dr = fabs((dn/n)*r);
  RooRealVar ans = RooRealVar(Form("%s_over_%s", num.GetName(), den.GetName()), Form("%s / %s", num.GetTitle(), den.GetTitle()), r);
  ans.setError(dr);
  return ans;
};

bool isSameCollisionSystem( TString name1, TString name2 )
{
  std::vector<std::string> col = { "_PP_", "_PbPb_", "_Pbp_", "_pPb_", "_PA_" };
  for (auto& c : col) { if (name1.Contains(c.c_str()) && name2.Contains(c.c_str())) return true; }
  return false;
};
bool isSameCharge( TString name1, TString name2 )
{
  std::vector<std::string> chg = { "Pl_", "Mi_" };
  for (auto& c : chg) { if (name1.Contains(c.c_str()) && name2.Contains(c.c_str())) return true; }
  return false;
};

bool binFromFile(const char* filename, anabin<0>& ans)
{
  TFile *f = TFile::Open(filename);
  if (!f) { std::cout << "[Error] File: " << filename << " does not exist!" << std::endl; return false; }
  if (!f->IsOpen()) { std::cout << "[Error] File: " << filename << " fail to open!" << std::endl; delete f; return false; }
  RooWorkspace* ws = (RooWorkspace*) f->Get("workspace");
  if (!ws) { std::cout << "[Error] File: " << filename << " does not have the workspace!" << std::endl; f->Close(); delete f; return false; }
  RooRealVar *mueta = (RooRealVar*) ws->var("Muon_Eta");
  RooRealVar *muiso = (RooRealVar*) ws->var("Muon_Iso");
  RooRealVar *cent  = (RooRealVar*) ws->var("Centrality");
  if (!mueta || !muiso || !cent) { std::cout << "[Error] File: " << filename << " does not have the EWQ variables." << std::endl; f->Close(); delete f;return false; }
  ans.setmuetabin(binF(mueta->getMin(),mueta->getMax())); ans.setmuisobin(binF(muiso->getMin(),muiso->getMax())); ans.setcentbin(binI(cent->getMin(),cent->getMax()));
  delete ws; f->Close(); delete f;
  return true;
};
bool binFromFile(const char* filename, anabin<1>& ans)
{
  TFile *f = TFile::Open(filename);
  if (!f) { std::cout << "[Error] File: " << filename << " does not exist!" << std::endl; return false; }
  if (!f->IsOpen()) { std::cout << "[Error] File: " << filename << " fail to open!" << std::endl; delete f; return false; }
  RooWorkspace* ws = (RooWorkspace*) f->Get("workspace");
  if (!ws) { std::cout << "[Error] File: " << filename << " does not have the workspace!" << std::endl; f->Close(); delete f; return false; }
  RooRealVar *pt   = (RooRealVar*) ws->var("rap");
  RooRealVar *rap  = (RooRealVar*) ws->var("pt");
  RooRealVar *cent = (RooRealVar*) ws->var("cent");
  if (!pt || !rap || !cent) { std::cout << "[Error] File: " << filename << " does not have the Onia variables." << std::endl; f->Close(); delete f; return false; }
  ans.setrapbin(binF(rap->getMin(),rap->getMax())); ans.setptbin(binF(pt->getMin(),pt->getMax())); ans.setcentbin(binI(cent->getMin(),cent->getMax()));
  delete ws; f->Close(); delete f;
  return true;
};

bool binok(std::vector< anabin<0> > thecats, const std::string& xaxis, anabin<0> &tocheck, bool override)
{
   bool ok=false;
   for (auto it : thecats) {
      if (xaxis=="eta" && it.muisobin()==tocheck.muisobin() && it.centbin()==tocheck.centbin()
            && ! (it.muetabin()==tocheck.muetabin())) {
         ok=true;
         if (override) tocheck.setmuetabin(it.muetabin());
         break;
      } else if (xaxis=="cent" && it.muisobin()==tocheck.muisobin() && it.muetabin()==tocheck.muetabin()
            && ! (it.centbin()==tocheck.centbin())) {
         ok=true;
         if (override) tocheck.setcentbin(it.centbin());
         break;
      } else if (xaxis=="iso" && it.centbin()==tocheck.centbin() && it.muetabin()==tocheck.muetabin()
            && ! (it.muisobin()==tocheck.muisobin())) {
         ok=true;
         if (override) tocheck.setmuisobin(it.muisobin());
         break;
      } else if ((it.centbin().low()<=0 && it.centbin().high()<=0)
            && it.muisobin()==tocheck.muisobin() && it.muetabin()==tocheck.muetabin()
            &&  (abs(it.centbin().low())==abs(tocheck.centbin().low()) && abs(it.centbin().high())==abs(tocheck.centbin().high()))) {
         ok=true;
         break;
      }
   }
   return ok;
};

bool binok(anabin<0> thecat, const std::string& xaxis, anabin<0> &tocheck, bool override)
{
  std::vector< anabin<0> > thecats; thecats.push_back(thecat);
   return binok(thecats, xaxis, tocheck, override);
};

bool isSameBinEWQ(const char* filename1, const char* filename2)
{
  TString fN1(filename1); TString fN2(filename2);
  if ( isSameCollisionSystem(filename1, filename2) && isSameCharge(filename1, filename2) ) { std::cout << "[ERROR] Comparing bins of same collision system and charge" << std::endl; return false; }
  anabin<0> thebin1(0,0,0,0,0,0), thebin2(0,0,0,0,0,0);
  binFromFile(fN1.Data(), thebin1);
  binFromFile(fN2.Data(), thebin2);
  if ( (thebin1.muetabin().low() == thebin2.muetabin().low()) && (thebin1.muetabin().high() == thebin2.muetabin().high())  && (thebin1.muisobin().low() == thebin2.muisobin().low()) && (thebin1.muisobin().high() == thebin2.muisobin().high()) ) return true;
  else return false;
};

std::pair<float, float> 
poiFromBin( const char* thepoiname, anabin<0> thebin, const char* theCollSystem, const char* theCharge, const char* workDirName, const char* DSTag, const char* prependPath )
{
  // Get the Content of the Result Tree
  TString tfname = treeFileName(workDirName, DSTag, prependPath, "MET");
  TFile *f = TFile::Open(tfname);
  if (!f || !f->IsOpen()) {
    resultsEWQ2tree(workDirName, DSTag, prependPath);
    f = new TFile(tfname);
    if (!f) return std::make_pair(-999., -999.);
    if (!f->IsOpen()) { delete f; return std::make_pair(-999., -999.); }
  }
  TTree *tr = (TTree*) f->Get("fitresults");
  if (!tr) { f->Close(); delete f; return std::make_pair(-999., -999.); }

  // fix centrality for pp and PA
  if (TString(theCollSystem) != "PbPb") thebin.setcentbin(binI(0,200));
  float muetamin, muetamax, muisomin, muisomax, centmin, centmax, val=-999., err=-999.;
  int valI=-999;
  char collSystem[5], charge[5];
  tr->SetBranchAddress("Muon_Eta_Min",&muetamin);
  tr->SetBranchAddress("Muon_Eta_Max",&muetamax);
  tr->SetBranchAddress("Muon_Iso_Min",&muisomin);
  tr->SetBranchAddress("Muon_Iso_Max",&muisomax);
  tr->SetBranchAddress("Centrality_Min",&centmin);
  tr->SetBranchAddress("Centrality_Max",&centmax);
  if (!TString(thepoiname).Contains("Chi2") && !TString(thepoiname).Contains("NDoF")) { tr->SetBranchAddress(Form("%s_Val",thepoiname),&val); tr->SetBranchAddress(Form("%s_Err",thepoiname),&err); }
  else if (!TString(thepoiname).Contains("NDoF")) { tr->SetBranchAddress(thepoiname,&val); err = 0.; }
  else { tr->SetBranchAddress(thepoiname,&valI); err = 0.; }
  tr->SetBranchAddress("collSystem",collSystem);
  tr->SetBranchAddress("charge",charge);
  int ntr = tr->GetEntries();
  bool found=false;
  for (int i=0; i<ntr; i++) {
    tr->GetEntry(i);
    if ((anabin<0>(muetamin, muetamax, muisomin, muisomax, centmin, centmax) == thebin) && (TString(collSystem) == TString(theCollSystem)) && (TString(charge) == TString(theCharge))) {
      found=true;
      break;
    }
  }
  if (!found) { val = -999.; err = -999.; }
  if (TString(thepoiname).Contains("NDoF")) val = float(valI);
  f->Close(); delete f;
  return std::make_pair(val, err);
};
std::pair<float, float> 
poiFromBin( const char* thepoiname, anabin<1> thebin, const char* theCollSystem, const char* workDirName, const char* DSTag, const char* prependPath )
{
  // Get the Content of the Result Tree
  TString tfname = treeFileName(workDirName, DSTag, prependPath, "MET");
  TFile *f = TFile::Open(tfname);
  if (!f || !f->IsOpen()) {
    //results2tree(workDirName, DSTag, prependPath); // NEEDS TO BE IMPLEMENTED
    f = new TFile(tfname);
    if (!f) return std::make_pair(-999., -999.);
    if (!f->IsOpen()) { delete f; return std::make_pair(-999., -999.); }
  }
  TTree *tr = (TTree*) f->Get("fitresults");
  if (!tr) { f->Close(); delete f; return std::make_pair(-999., -999.); }

  // fix centrality for pp and PA
  if (TString(theCollSystem) != "PbPb") thebin.setcentbin(binI(0,200));
  float ptmin, ptmax, ymin, ymax, centmin, centmax, val=-999., err=-999.;
  int valI=-999;
  char collSystem[5], charge[5];
  tr->SetBranchAddress("ptmin",&ptmin);
  tr->SetBranchAddress("ptmax",&ptmax);
  tr->SetBranchAddress("ymin",&ymin);
  tr->SetBranchAddress("ymax",&ymax);
  tr->SetBranchAddress("centmin",&centmin);
  tr->SetBranchAddress("centmax",&centmax);
  if (!TString(thepoiname).Contains("Chi2") && !TString(thepoiname).Contains("NDoF")) { tr->SetBranchAddress(Form("%s_Val",thepoiname),&val); tr->SetBranchAddress(Form("%s_Err",thepoiname),&err); }
  else if (!TString(thepoiname).Contains("NDoF")) { tr->SetBranchAddress(thepoiname,&val); err = 0.; }
  else { tr->SetBranchAddress(thepoiname,&valI); err = 0.; }
  tr->SetBranchAddress("collSystem",collSystem);
  int ntr = tr->GetEntries();
  bool found=false;
  for (int i=0; i<ntr; i++) {
    tr->GetEntry(i);
    if ((anabin<1>(ymin, ymax, ptmin, ptmax, centmin, centmax) == thebin) && (TString(collSystem) == TString(theCollSystem))) {
      found=true;
      break;
    }
  }
  if (!found) { val = -999.; err = -999.; }
  if (TString(thepoiname).Contains("NDoF")) val = float(valI);
  f->Close(); delete f; 
  return std::make_pair(val, err);
};

// For analysis bins
void prune( std::vector<anabin<0>>& v, const bool& keepshort )
{
  std::vector<anabin<0>> ans;
  for (auto it1 : v) {
    bool binok=true;
    for (auto it2 : v) {
      if (it1==it2) continue;
      if (it1.muetabin()==it2.muetabin() && it1.muisobin()==it2.muisobin()) {
        binI cb1 = it1.centbin();
        binI cb2 = it2.centbin();
        if (!(cb1==binI(0,200)) && cb1.low()==cb2.low()) { // the bin is not MB and there is another bin with the same lower edge
          if (keepshort && cb1.high()>cb2.high()) binok=false;
          if (!keepshort && cb1.high()<cb2.high()) binok=false;
        }
      } // same pt and rap bins
    } // for it2
    if (binok) ans.push_back(it1);
  } // for it1
  v = ans;
};
void prune( std::vector<anabin<1>>& v, const bool& keepshort )
{
   vector<anabin<1>> ans;
   for (auto it1 : v) {
      bool binok=true;
      for (auto it2 : v) {
         if (it1==it2) continue;
         if (it1.rapbin()==it2.rapbin() && it1.ptbin()==it2.ptbin()) {
            binI cb1 = it1.centbin();
            binI cb2 = it2.centbin();
            if (!(cb1==binI(0,200)) && cb1.low()==cb2.low()) { // the bin is not MB and there is another bin with the same lower edge
               if (keepshort && cb1.high()>cb2.high()) binok=false;
               if (!keepshort && cb1.high()<cb2.high()) binok=false;
            }
         } // same pt and rap bins
      } // for it2
      if (binok) ans.push_back(it1);
   } // for it1
   v = ans;
};

// For TGraph
void prune( TGraphAsymmErrors *g, TGraphAsymmErrors *gsyst, const bool& keepshort )
{
  int n = g->GetN();
  for (int i1=0; i1<n; i1++) {
    double xl1 = g->GetX()[i1]-g->GetErrorXlow(i1);
    double xh1 = g->GetX()[i1]+g->GetErrorXhigh(i1);
    bool binok=true;
    for (int i2=0; i2<n; i2++) {
      if (i2==i1) continue;
      double xl2 = g->GetX()[i2]-g->GetErrorXlow(i2);
      double xh2 = g->GetX()[i2]+g->GetErrorXhigh(i2);
      if (fabs(xl1-xl2)<1e-3 && keepshort && xh1>xh2) binok=false;
      if (fabs(xl1-xl2)<1e-3 && !keepshort && xh1<xh2) binok=false;
    } // for i2
    if (!binok) {
      g->SetPoint(i1,g->GetX()[i1],-g->GetY()[i1]);
      if (gsyst) gsyst->SetPoint(i1,gsyst->GetX()[i1],-gsyst->GetY()[i1]);
    }
  } // for i1
};
void prune( TGraphErrors *g, const bool& keepshort )
{
  int n = g->GetN();
  for (int i1=0; i1<n; i1++) {
    double xl1 = g->GetX()[i1]-g->GetErrorX(i1);
    double xh1 = g->GetX()[i1]+g->GetErrorX(i1);
    bool binok=true;
    for (int i2=0; i2<n; i2++) {
      if (i2==i1) continue;
      double xl2 = g->GetX()[i2]-g->GetErrorX(i2);
      double xh2 = g->GetX()[i2]+g->GetErrorX(i2);
      if (fabs(xl1-xl2)<1e-3 && keepshort && xh1>xh2) binok=false;
      if (fabs(xl1-xl2)<1e-3 && !keepshort && xh1<xh2) binok=false;
    } // for i2
    if (!binok) g->SetPoint(i1,g->GetX()[i1],-g->GetY()[i1]);
  } // for i1
};


#endif // ifndef resultUtils_h
  

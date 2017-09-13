#ifndef resultsEWQ2tree_C
#define resultsEWQ2tree_C

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH1.h"
#include "TMath.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooAbsData.h"
#include "RooStringVar.h"
#include "RooWorkspace.h"
#include "RooPlot.h"
#include "RooHist.h"

#include <vector>
#include <cstring>

#include "Macros/Utilities/resultsUtils.h"


struct poi {
  Char_t Name[64];
  float Val;
  float Err;
  float Min;
  float Max;
  float parIni_Val;
  float parIni_Err;
};

const int nBins = 46;

void resultsEWQ2tree(
                     const char* workDirName,
                     const char* DSTag,
                     const char* prependPath,
                     const char* fitType,
                     const char* thePoiNames
                     ) {
  // workDirName: usual tag where to look for files in Output
  // thePoiNames: comma-separated list of parameters to store ("par1,par2,par3"). Default: all

  TFile *f  = new TFile(treeFileName(workDirName,DSTag,prependPath,fitType), "RECREATE");
  TTree *tr = new TTree("fitresults", "fit results");

  // Event Variables
  std::map< std::string , std::map < std::string , float > > evtVar = 
    {
      { "MET",        {{"Min" , -99.}, {"Max" , -99.}, {"Val" , -99.}, {"Err" , -99.}} },
      { "Muon_Pt",    {{"Min" , -99.}, {"Max" , -99.}, {"Val" , -99.}, {"Err" , -99.}} },
      { "Muon_Eta",   {{"Min" , -99.}, {"Max" , -99.}, {"Val" , -99.}, {"Err" , -99.}} },
      { "Muon_Iso",   {{"Min" , -99.}, {"Max" , -99.}, {"Val" , -99.}, {"Err" , -99.}} },
      { "Muon_MT",    {{"Min" , -99.}, {"Max" , -99.}, {"Val" , -99.}, {"Err" , -99.}} },
      { "Centrality", {{"Min" , -99.}, {"Max" , -99.}, {"Val" , -99.}, {"Err" , -99.}} }
    };
  // model names
  Char_t W_Model[128]="" , WToTau_Model[128]="", DYZ_Model[128]="", QCD_Model[128]="";
  // Cut Value
  Char_t CutAndCount_WToMu[200]="";
  // system information
  Char_t collSystem[8]="", charge[5]="";
  // goodness of fit
  std::map< std::string , std::map < std::string , float > > fitPar = 
    { 
      { "MET", {{"Chi2" , -99.}, {"NDoF" , -99.}, {"NLL" , -99.}, {"NPar" , -99.}, {"NormChi2" , -99.}, {"Chi2Prob" , -99.}} }
    };
  // parame
  TString thePoiNamesStr(thePoiNames);
  if (thePoiNamesStr.EqualTo("all")) {
    thePoiNamesStr = TString("N_WToMu,N_WToTauToMu,N_DYZToMu,N_QCDToMu,Alpha_QCDToMu,Beta_QCDToMu,x0_QCDToMu");
  }
  std::vector<poi> thePois;
  TString t; Int_t from = 0;
  while (thePoiNamesStr.Tokenize(t, from , ",")) {
    poi p; strcpy(p.Name, t.Data());
    cout << p.Name << endl;
    thePois.push_back(p);
  }

  // create tree branches
  // Event Variables
  for (auto& v : evtVar) {
    for (auto& p : v.second) {
      tr->Branch(Form("%s_%s", v.first.c_str(), p.first.c_str()) ,&p.second, Form("%s_%s/F", v.first.c_str(), p.first.c_str()));
    }
  }
  // Goodness of fit
  for (auto& v : fitPar) {
    for (auto& p : v.second) {
      tr->Branch(Form("%s_%s", v.first.c_str(), p.first.c_str()) ,&p.second, Form("%s_%s/F", v.first.c_str(), p.first.c_str()));
    }
  }
  // model names
  tr->Branch("Model_W",W_Model,"Model_W/C");
  tr->Branch("Model_WToTau",WToTau_Model,"Model_WToTau/C");
  tr->Branch("Model_DYZ",DYZ_Model,"Model_DYZ/C");
  tr->Branch("Model_QCD",QCD_Model,"Model_QCD/C");
  // system information
  tr->Branch("collSystem",collSystem,"collSystem/C");
  tr->Branch("charge",charge,"collSystem/C");
  // parameters to store
  tr->Branch("CutAndCount_WToMu",CutAndCount_WToMu,"CutAndCount_WToMu/C");
  for (vector<poi>::iterator it=thePois.begin(); it!=thePois.end(); it++) {
    tr->Branch(Form("%s_Val",it->Name),&(it->Val),Form("%s_Val/F",it->Name));
    tr->Branch(Form("%s_Err",it->Name),&(it->Err),Form("%s_Err/F",it->Name));
    tr->Branch(Form("%s_Min",it->Name),&(it->Min),Form("%s_Min/F",it->Name));
    tr->Branch(Form("%s_Max",it->Name),&(it->Max),Form("%s_Max/F",it->Name));
    tr->Branch(Form("%s_parIni_Val",it->Name),&(it->parIni_Val),Form("%s_parIni_Val/F",it->Name));
    tr->Branch(Form("%s_parIni_Err",it->Name),&(it->parIni_Err),Form("%s_parIni_Err/F",it->Name));
  }

  // list of files
  std::vector<TString> theFiles = fileList(workDirName,"",DSTag,"",fitType);
  std::cout << theFiles.size() << std::endl;

  int cnt=0;
  for (auto& fileN : theFiles) {
    cout << "Parsing file " << cnt << " / " << theFiles.size() << ": " << fileN << endl;
    // parse the file name to get info
    if (fileN.Contains("Pbp_"))  { strcpy(collSystem, "Pbp");  }
    if (fileN.Contains("pPb_"))  { strcpy(collSystem, "pPb");  }
    if (fileN.Contains("PA_"))   { strcpy(collSystem, "PA");   }
    if (fileN.Contains("PbPb_")) { strcpy(collSystem, "PbPb"); }
    if (fileN.Contains("PP_"))   { strcpy(collSystem, "PP");   }
    if (fileN.Contains("Pl_"))   { strcpy(charge,     "Pl");   }
    if (fileN.Contains("Mi_"))   { strcpy(charge,     "Mi");   }
    const char* dsName = Form("d%s_%s_MUON_%s", charge, DSTag, collSystem);
    const char* token = Form("%s_%s", charge, collSystem);
    
    TFile *f = TFile::Open(fileN.Data());
    if (!f) { std::cout << "[Error] File: " << fileN << " does not exist!" << std::endl; return; }
    if (!f->IsOpen()) { std::cout << "[Error] File: " << fileN << " fail to open!" << std::endl; delete f; return; }
    RooWorkspace* ws = (RooWorkspace*) f->Get("workspace");
    if (!ws) { std::cout << "[Error] File: " << fileN << " does not have the workspace!" << std::endl; f->Close(); delete f; return; }

    // get the snapshots
    const RooArgSet *parIni = ws->getSnapshot("initialParameters");

    // get the POIs
    for (auto& itpoi : thePois) {
      RooRealVar *thevar = poiFromWS(*ws, token, itpoi.Name);
      RooRealVar *thevar_parIni = parIni ? (RooRealVar*) parIni->find(Form("%s%s", itpoi.Name, token)) : NULL;
      itpoi.Val = thevar ? thevar->getVal()   : 0;
      itpoi.Err = thevar ? thevar->getError() : 0;
      itpoi.Min = thevar ? thevar->getMin()   : 0;
      itpoi.Max = thevar ? thevar->getMax()   : 0;
      itpoi.parIni_Val = thevar_parIni ? thevar_parIni->getVal()   : 0;
      itpoi.parIni_Err = thevar_parIni ? thevar_parIni->getError() : 0;
    }

    RooAbsData * ds = (RooAbsData*) ws->data(Form("CutAndCount_%s", dsName));
    if (!ds) { ds = (RooAbsData*) ws->data(Form("%s", dsName)); }
    for (auto& v : evtVar) {
      RooRealVar *thevar  = parIni ? (RooRealVar*) parIni->find(v.first.c_str()) : NULL;
      RooRealVar *themean = (thevar && ds) ? (RooRealVar*) ds->meanVar(*thevar) : NULL;
      if (v.second.count("Min")) { v.second.at("Min") = thevar ? thevar->getMin() : 0; }
      if (v.second.count("Max")) { v.second.at("Max") = thevar ? thevar->getMax() : 0; }
      if (v.second.count("Val")) { v.second.at("Val") = themean ? themean->getVal() : 0; }
      if (v.second.count("Err")) { v.second.at("Err") = themean ? themean->getError() : 0; }
    }
    RooStringVar *thecut  = (RooStringVar*) ws->obj(Form("CutAndCount_WToMu%s", token));
    strcpy(CutAndCount_WToMu, thecut ? thecut->getVal() : "");

    delete ws; f->Close(); delete f;

    // fill the tree
    tr->Fill();
    cnt++;
  } // loop on the files

  f->Write();
  f->Close();
}


#endif // #ifndef resultsEWQ2tree_C

#ifndef Histogram2_h
#define Histogram2_h

// Header file for ROOT classes
#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <THStack.h>
#include <TPad.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TCanvas.h>

// Header file for c++ classes
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <iterator>

// CMS STYLE
#include "CMS/tdrstyle.C"
#include "CMS/CMS_lumi.C"


typedef struct TypeInfo {
  std::vector<std::string>   sample;
  std::vector<std::string>   cutSelection;
} TypeInfo;

struct DrawInfo {
  bool          yLogScale;
  DrawInfo& operator=(const DrawInfo& a)
  {
    this->yLogScale = a.yLogScale;
    return *this;
  }
};

typedef struct VarInfo {
  std::string   label;
  unsigned int  nBin;
  float         min;
  float         max;
  DrawInfo drawInfo;
} VarInfo;

template<typename Out>
void split(const std::string &s, char delim, Out result) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
      *(result++) = item;
    }
}

std::vector<int> COLOR  = { kRed, kGreen+2, kBlue+2, kOrange+2, kViolet+2, kMagenta+2, kAzure-7, kBlack, kViolet };
std::vector<int> MARKER = { 20, 24, 21, 34, 22 };
std::vector<int> LINE   = { 2, 1, 3, 4, 5, 6 };

class Histogram2 {

 public :

  Histogram2();
  virtual ~Histogram2();
  
  virtual void         Book    ( const std::string&, const std::string&, const std::map< std::string , struct VarInfo >& );
  virtual void         Fill    ( const std::string&, const std::string&, const std::map< std::string , float >& );
  virtual void         Fill    ( const std::string&, const std::string&, const std::map< std::string , std::pair< float , float > >& );
  virtual void         Book    ( const std::string& sample, const std::map< std::string , struct VarInfo >& varMap               ) { Book(sample, "", varMap);   }
  virtual void         Fill    ( const std::string& sample, const std::map< std::string , float >& valueMap                      ) { Fill(sample, "", valueMap); }
  virtual void         Fill    ( const std::string& sample, const std::map< std::string , std::pair< float , float > >& valueMap ) { Fill(sample, "", valueMap); }
  virtual void         Draw    ( const std::string& tag="", const std::map< std::string , struct TypeInfo >& typeInfo={});
  virtual void         Delete  ( void );
  virtual void         Norm    ( void );
  virtual void         Fit     ( const std::string& fName="gaus" );
  virtual void         Save    ( const std::string& );
  virtual void         ComputeRatio ( const std::string& rName="DATAvsMC" );
  virtual void         PrintMean ( void );
  virtual void         PrintFinalMean ( void );
  virtual void         PrintFitResults ( void );
  
 private :
  virtual void         MakeStackHistogram ( const std::string& tag="" );
  virtual void         setYRange ( TH1D*    , TH1D* ref=NULL, bool logScale=false);
  virtual void         setYRange ( THStack* , TH1D* ref=NULL, bool logScale=false);
  virtual std::string  getBeam   ( const std::string& );
  virtual std::string  getDecayChannel ( const std::string& );
  virtual std::string  GetString ( const TF1&  );
  virtual std::string  GetString ( const TH1D& , bool Gauss=false );

  std::map< std::string ,std::map< std::string ,std::map< std::string , DrawInfo > > > drawInfo_;
  std::map< std::string ,std::map< std::string ,std::map< std::string , TH1D*    > > > TH1D_;
  std::map< std::string ,std::map< std::string ,std::map< std::string , TH1D*    > > > TH1D_RATIO_;
  std::map< std::string ,std::map< std::string ,std::map< std::string , TH1D*    > > > TH1D_RATIO_NORM_;
  std::map< std::string ,std::map< std::string ,std::map< std::string , TF1*     > > > TF1_;
  std::map< std::string ,std::map< std::string ,std::map< std::string , TF1*     > > > TF1_RATIO_;
  std::map< std::string ,std::map< std::string ,std::map< std::string , THStack* > > > THStack_;

};

Histogram2::Histogram2()
{
  // set the CMS style
  setTDRStyle();
}
 
Histogram2::~Histogram2()
{
   Delete();
}

void 
Histogram2::Book(const std::string& sample, const std::string& type, const std::map< std::string , struct VarInfo >& varMap)
{
  for (auto& var : varMap) {
    std::string    varName = var.first;
    struct VarInfo varInfo = var.second;
    std::string histName = (std::string("h_") + (type!="" ? (sample + "_" + type) : sample) + "_" + varName);
    if (TH1D_.count(sample)==0 || TH1D_[sample].count(type)==0 || TH1D_[sample][type].count(varName)==0) {
      // Create the Histogram
      if (TH1D_[sample][type][varName]) delete TH1D_[sample][type][varName];
      TH1D_[sample][type][varName] = new TH1D(histName.c_str(), histName.c_str(), varInfo.nBin, varInfo.min, varInfo.max);
      drawInfo_[sample][type][varName] = varInfo.drawInfo;
      // Initialize the Histogram
      TH1D_[sample][type][varName]->GetYaxis()->SetTitle("Number of Entries");
      TH1D_[sample][type][varName]->GetXaxis()->SetTitle(varInfo.label.c_str());
      TH1D_[sample][type][varName]->Sumw2(kTRUE);
      std::cout << "[INFO] Added histogram: " << histName <<  std::endl;
    }
  }
}

void 
Histogram2::Fill(const std::string& sample, const std::string& type, const std::map< std::string , float >& valueMap)
{
  for (auto& value : valueMap) {
    std::string varName  = value.first;
    if (TH1D_.count(sample)>0 && TH1D_.at(sample).count(type)>0 && TH1D_.at(sample).at(type).count(varName)>0 && TH1D_.at(sample).at(type).at(varName)) {
      TH1D_.at(sample).at(type).at(varName)->Fill(value.second);
    }
  }
}

void 
Histogram2::Fill(const std::string& sample, const std::string& type, const std::map< std::string , std::pair< float , float > >& valueMap)
{
  for (auto& value : valueMap) {
    std::string varName  = value.first;
    if (TH1D_.count(sample)>0 && TH1D_.at(sample).count(type)>0 && TH1D_.at(sample).at(type).count(varName)>0 && TH1D_.at(sample).at(type).at(varName)) {
      TH1D_.at(sample).at(type).at(varName)->Fill(value.second.first, value.second.second);
    }
  }
}

void 
Histogram2::Norm(void)
{
  for (auto& s : TH1D_) {
    for (auto& t : s.second) {
      for (auto& v : t.second) {
        if (v.second) {
          if (v.second->GetEntries()>0.) {
            v.second->ClearUnderflowAndOverflow();
            v.second->Scale(1.0/v.second->GetSumOfWeights());
          }
        }
      }
    }
  }
}

void 
Histogram2::Fit(const std::string& fName)
{
  for (auto& s : TH1D_) {
    std::string sample = s.first;
    for (auto& t : TH1D_[sample]) {
      std::string type = t.first;
      for (auto& elem : TH1D_[sample][type]) {
        std::string varName  = elem.first;
        std::string funcName = (std::string("f_") + sample + "_" + type + "_" + varName);
        if (type=="") { funcName = (std::string("f_") + sample + "_" + varName); }
        std::string funcRATIOName = (std::string("f_RATIO_") + sample + "_" + type + "_" + varName);
        if (type=="") { funcRATIOName = (std::string("f_RATIO_") + sample + "_" + varName); }
        double minX = TH1D_[sample][type][varName]->GetXaxis()->GetXmin();
        double maxX = TH1D_[sample][type][varName]->GetXaxis()->GetXmax();
        if (fName=="gaus") {
          if (TF1_[sample][type][varName]) delete TF1_[sample][type][varName];
          TF1_[sample][type][varName] = new TF1(funcName.c_str(), "gaus", minX, maxX);
          double maxY = TH1D_[sample][type][varName]->GetBinContent(TH1D_[sample][type][varName]->GetMaximumBin());
          TF1_[sample][type][varName]->SetParameter(0,maxY);
          //TF1_[sample][type][varName]->FixParameter(0,maxY);
          TF1_[sample][type][varName]->SetParameter(1,TH1D_[sample][type][varName]->GetMean());
          TF1_[sample][type][varName]->SetParameter(2,TH1D_[sample][type][varName]->GetStdDev());
          if (TF1_RATIO_[sample][type][varName]) delete TF1_RATIO_[sample][type][varName];
          TF1_RATIO_[sample][type][varName] = new TF1(funcRATIOName.c_str(), "gaus(0)/(gaus(3))", minX, maxX);
        }
        if (TF1_[sample][type][varName]) {
          TH1D_[sample][type][varName]->Fit(funcName.c_str(),"0");
        }
      }
    }
  }
  // Apply General Settings to all TF1
  for (auto& s : TF1_) { for (auto& t : s.second) { for (auto& v : t.second) { if (v.second) {
          v.second->SetLineColor(kBlack);
          v.second->SetLineWidth(3);
        }
      }
    }
  }
}

void 
Histogram2::ComputeRatio(const std::string& rName)
{
  // For TH1D
  if (rName=="DATAvsFIT" || rName=="DATAvsMC") {
    for (auto& s : TH1D_) {
      std::string sample = s.first;
      for (auto& t : s.second) {
        std::string type = t.first;
        for (auto& v : t.second) {
          std::string varName  = v.first;
          std::string histName = (std::string("hRATIO_") + (type!="" ? (sample + "_" + type) : sample) + "_" + varName);
          std::string histNormName = (std::string("hRATIO_NORM_") + (type!="" ? (sample + "_" + type) : sample) + "_" + varName);
          double minX = v.second->GetXaxis()->GetXmin();
          double maxX = v.second->GetXaxis()->GetXmax();
          if (rName=="DATAvsMC" && sample.find("MC")!=std::string::npos) {
            // Create Data Histogram
            if (TH1D_RATIO_[sample][type][varName]!=NULL) delete TH1D_RATIO_[sample][type][varName];
            TH1D_RATIO_[sample][type][varName] = new TH1D(histName.c_str(), histName.c_str(), v.second->GetNbinsX(), minX, maxX);
            TH1D_RATIO_.at(sample).at(type).at(varName)->SetDirectory(0);
            if (TH1D_RATIO_NORM_[sample][type][varName]!=NULL) delete TH1D_RATIO_NORM_[sample][type][varName];
            TH1D_RATIO_NORM_[sample][type][varName] = new TH1D(histNormName.c_str(), histNormName.c_str(), v.second->GetNbinsX(), minX, maxX);
            TH1D_RATIO_NORM_.at(sample).at(type).at(varName)->SetDirectory(0);
            // Initialize the Histogram
            TH1D_RATIO_.at(sample).at(type).at(varName)->GetYaxis()->SetTitle("DATA/MC");
            TH1D_RATIO_NORM_.at(sample).at(type).at(varName)->GetYaxis()->SetTitle("DATA/MC");
            // Make the ratio
            std::string sampleData = "DATA_" + getBeam(sample);
            if (TH1D_.count(sampleData)>0 && TH1D_.at(sampleData).count(type)>0 && TH1D_.at(sampleData).at(type).count(varName)>0 && 
                TH1D_.at(sampleData).at(type).at(varName) && TH1D_.at(sampleData).at(type).at(varName)->GetEntries()>0.) {
              TH1D_.at(sampleData).at(type).at(varName)->ClearUnderflowAndOverflow();
              bool isNorm = false; //if (v.second->GetSumOfWeights()<=1.5) { isNorm = true; }
              if (isNorm) { v.second->Scale(v.second->GetSumOfWeights()); }
              if (isNorm) { TH1D_.at(sampleData).at(type).at(varName)->Scale(TH1D_.at(sampleData).at(type).at(varName)->GetSumOfWeights()); }
              if (isNorm) { 
                TH1D_RATIO_.at(sample).at(type).at(varName)->Divide(TH1D_.at(sampleData).at(type).at(varName), v.second, 
                                                                    (1./TH1D_.at(sampleData).at(type).at(varName)->GetSumOfWeights()), 
                                                                    (1./v.second->GetSumOfWeights()), 
                                                                    "b");
              }
              else {
                TH1D_RATIO_.at(sample).at(type).at(varName)->Divide(TH1D_.at(sampleData).at(type).at(varName), v.second, 1.0, 1.0, "b");
              }
              TH1D_RATIO_NORM_.at(sample).at(type).at(varName)->Divide(TH1D_.at(sampleData).at(type).at(varName), v.second, 
                                                                       (1./TH1D_.at(sampleData).at(type).at(varName)->GetSumOfWeights()), 
                                                                       (1./v.second->GetSumOfWeights()), 
                                                                       "b");
              if (isNorm) { v.second->Scale(1./v.second->GetSumOfWeights()); }
              if (isNorm) { TH1D_.at(sampleData).at(type).at(varName)->Scale(1./TH1D_.at(sampleData).at(type).at(varName)->GetSumOfWeights()); }
            }
            else {
              TH1D* tmp = new TH1D("TMP", "TMP", v.second->GetNbinsX(), minX, maxX);
              TH1D_RATIO_.at(sample).at(type).at(varName)->Divide(tmp, v.second, 1.0, 1.0, "b");
              delete tmp;
            }
            if (TF1_.count(sampleData)>0 && TF1_.at(sampleData).count(type)>0 && TF1_.at(sampleData).at(type).count(varName)>0 && TF1_.at(sampleData).at(type).at(varName)) {        
              TF1_RATIO_.at(sample).at(type).at(varName)->SetParameter(0, TF1_.at(sampleData).at(type).at(varName)->GetParameter(0));
              TF1_RATIO_.at(sample).at(type).at(varName)->SetParameter(1, TF1_.at(sampleData).at(type).at(varName)->GetParameter(1));
              TF1_RATIO_.at(sample).at(type).at(varName)->SetParameter(2, TF1_.at(sampleData).at(type).at(varName)->GetParameter(2));
              TF1_RATIO_.at(sample).at(type).at(varName)->SetParameter(3, TF1_.at(sample).at(type).at(varName)->GetParameter(0));
              TF1_RATIO_.at(sample).at(type).at(varName)->SetParameter(4, TF1_.at(sample).at(type).at(varName)->GetParameter(1));
              TF1_RATIO_.at(sample).at(type).at(varName)->SetParameter(5, TF1_.at(sample).at(type).at(varName)->GetParameter(2));
            }
          }
          if (rName=="DATAvsFIT" && TF1_RATIO_.count(sample)>0 && TF1_RATIO_.at(sample).count(type)>0 && TF1_RATIO_.at(sample).at(type).count(varName)>0 && TF1_RATIO_.at(sample).at(type).at(varName)) {
            // Create Data Histogram
            if (TH1D_RATIO_.at(sample).at(type).at(varName)) delete TH1D_RATIO_.at(sample).at(type).at(varName);
            TH1D_RATIO_.at(sample).at(type).at(varName) = new TH1D(histName.c_str(), histName.c_str(), v.second->GetNbinsX(), minX, maxX);
            TH1D_RATIO_.at(sample).at(type).at(varName)->SetDirectory(0);
            // Create Fit Histogram
            TH1D* Ref = (TH1D *) TH1D_RATIO_.at(sample).at(type).at(varName)->Clone("REF");
            Ref->FillRandom(TF1_.at(sample).at(type).at(varName)->GetName(), v.second->GetEffectiveEntries());
            if (v.second->GetSumOfWeights()<=1.5) { Ref->Scale(1.0/Ref->Integral()); }
            // Initialize the Histogram
            TH1D_RATIO_.at(sample).at(type).at(varName)->GetYaxis()->SetTitle("DATA/FIT");
            TH1D_RATIO_.at(sample).at(type).at(varName)->Divide(v.second, Ref, 1., 1., "b");
            delete Ref;
          }
        }
      }
    }
  }
  // For THStack
  if (rName=="DATAvsMCStack") {
    for (auto& s : THStack_) {
      std::string sample = s.first;
      for (auto& t : s.second) {
        std::string type = t.first;
        for (auto& v : t.second) {
          if (v.second) {
            std::string varName  = v.first;
            std::string histName = (std::string("hRATIO_") + (type!="" ? (sample + "_" + type) : sample) + "_" + varName);
            std::string histNormName = (std::string("hRATIO_NORM_") + (type!="" ? (sample + "_" + type) : sample) + "_" + varName);
            TH1D* hStack = ((TH1D*)v.second->GetStack()->At(v.second->GetNhists()-1));
            double minX = hStack->GetXaxis()->GetXmin();
            double maxX = hStack->GetXaxis()->GetXmax();
            if (rName=="DATAvsMCStack" && sample.find("MC")!=std::string::npos) {
              // Create Data Histogram
              if (TH1D_RATIO_[sample][type][varName]!=NULL) delete TH1D_RATIO_[sample][type][varName];
              TH1D_RATIO_[sample][type][varName] = new TH1D(histName.c_str(), histName.c_str(), hStack->GetNbinsX(), minX, maxX);
              TH1D_RATIO_.at(sample).at(type).at(varName)->SetDirectory(0);
              if (TH1D_RATIO_NORM_[sample][type][varName]!=NULL) delete TH1D_RATIO_NORM_[sample][type][varName];
              TH1D_RATIO_NORM_[sample][type][varName] = new TH1D(histNormName.c_str(), histNormName.c_str(), hStack->GetNbinsX(), minX, maxX);
              TH1D_RATIO_NORM_.at(sample).at(type).at(varName)->SetDirectory(0);
              // Initialize the Histogram
              TH1D_RATIO_.at(sample).at(type).at(varName)->GetYaxis()->SetTitle("DATA/MC");
              TH1D_RATIO_NORM_.at(sample).at(type).at(varName)->GetYaxis()->SetTitle("DATA/MC");
              // Make the ratio
              std::string sampleData = "DATA_" + getBeam(sample);
              if (TH1D_.count(sampleData)>0 && TH1D_.at(sampleData).count(type)>0 && TH1D_.at(sampleData).at(type).count(varName)>0 && 
                  TH1D_.at(sampleData).at(type).at(varName) && TH1D_.at(sampleData).at(type).at(varName)->GetEntries()>0.) {
                TH1D_.at(sampleData).at(type).at(varName)->ClearUnderflowAndOverflow();
                TH1D_RATIO_.at(sample).at(type).at(varName)->Divide(TH1D_.at(sampleData).at(type).at(varName), hStack, 1.0, 1.0, "b");
                TH1D_RATIO_NORM_.at(sample).at(type).at(varName)->Divide(TH1D_.at(sampleData).at(type).at(varName), hStack, 
                                                                         1.0/TH1D_.at(sampleData).at(type).at(varName)->GetSumOfWeights(), 
                                                                         1.0/hStack->GetSumOfWeights(), 
                                                                         "b");
              }
              else {
                TH1D* tmp = new TH1D("TMP", "TMP", hStack->GetNbinsX(), minX, maxX);
                TH1D_RATIO_.at(sample).at(type).at(varName)->Divide(tmp, hStack, 1.0, 1.0, "b");
                delete tmp;
              }
            }
          }
        }
      }
    }
  }
  // For TH1D
  if (rName=="MCNLOvsMCLO") {
    for (auto& s : TH1D_) {
      std::string sample = s.first;
      for (auto& t : s.second) {
        std::string type = t.first;
        for (auto& v : t.second) {
          std::string varName  = v.first;
          std::string histName = (std::string("hRATIO_") + (type!="" ? (sample + "_" + type) : sample) + "_" + varName);
          std::string histNormName = (std::string("hRATIO_NORM_") + (type!="" ? (sample + "_" + type) : sample) + "_" + varName);
          double minX = v.second->GetXaxis()->GetXmin();
          double maxX = v.second->GetXaxis()->GetXmax();
          if (rName=="MCNLOvsMCLO" && sample.find("MCNLO")!=std::string::npos) {
            // Create Data Histogram
            if (TH1D_RATIO_[sample][type][varName]!=NULL) delete TH1D_RATIO_[sample][type][varName];
            TH1D_RATIO_[sample][type][varName] = new TH1D(histName.c_str(), histName.c_str(), v.second->GetNbinsX(), minX, maxX);
            TH1D_RATIO_.at(sample).at(type).at(varName)->SetDirectory(0);
            if (TH1D_RATIO_NORM_[sample][type][varName]!=NULL) delete TH1D_RATIO_NORM_[sample][type][varName];
            TH1D_RATIO_NORM_[sample][type][varName] = new TH1D(histNormName.c_str(), histNormName.c_str(), v.second->GetNbinsX(), minX, maxX);
            TH1D_RATIO_NORM_.at(sample).at(type).at(varName)->SetDirectory(0);
            // Initialize the Histogram
            TH1D_RATIO_.at(sample).at(type).at(varName)->GetYaxis()->SetTitle("MC_NLO/MC_LO");
            TH1D_RATIO_NORM_.at(sample).at(type).at(varName)->GetYaxis()->SetTitle("MC_NLO/MC_LO");
            // Make the ratio
            std::string sampleMCLO = "MCLO_" + getDecayChannel(sample) + getBeam(sample);
            if (TH1D_.count(sampleMCLO)>0 && TH1D_.at(sampleMCLO).count(type)>0 && TH1D_.at(sampleMCLO).at(type).count(varName)>0 && 
                TH1D_.at(sampleMCLO).at(type).at(varName) && TH1D_.at(sampleMCLO).at(type).at(varName)->GetEntries()>0.) {
              TH1D_.at(sampleMCLO).at(type).at(varName)->ClearUnderflowAndOverflow();
              bool isNorm = false; //if (v.second->GetSumOfWeights()<=1.5) { isNorm = true; }
              if (isNorm) { v.second->Scale(v.second->GetSumOfWeights()); }
              if (isNorm) { TH1D_.at(sampleMCLO).at(type).at(varName)->Scale(TH1D_.at(sampleMCLO).at(type).at(varName)->GetSumOfWeights()); }
              if (isNorm) { 
                TH1D_RATIO_.at(sample).at(type).at(varName)->Divide(v.second, TH1D_.at(sampleMCLO).at(type).at(varName), 
                                                                    (1./v.second->GetSumOfWeights()), 
                                                                    (1./TH1D_.at(sampleMCLO).at(type).at(varName)->GetSumOfWeights()), 
                                                                    "b");
              }
              else {
                TH1D_RATIO_.at(sample).at(type).at(varName)->Divide(v.second, TH1D_.at(sampleMCLO).at(type).at(varName), 1.0, 1.0, "b");
              }
              TH1D_RATIO_NORM_.at(sample).at(type).at(varName)->Divide(v.second, TH1D_.at(sampleMCLO).at(type).at(varName), 
                                                                       (1./v.second->GetSumOfWeights()), 
                                                                       (1./TH1D_.at(sampleMCLO).at(type).at(varName)->GetSumOfWeights()), 
                                                                       "b");
              if (isNorm) { v.second->Scale(1./v.second->GetSumOfWeights()); }
              if (isNorm) { TH1D_.at(sampleMCLO).at(type).at(varName)->Scale(1./TH1D_.at(sampleMCLO).at(type).at(varName)->GetSumOfWeights()); }
            }
            else {
              TH1D* tmp = new TH1D("TMP", "TMP", v.second->GetNbinsX(), minX, maxX);
              TH1D_RATIO_.at(sample).at(type).at(varName)->Divide(tmp, v.second, 1.0, 1.0, "b");
              delete tmp;
            }
          }
        }
      }
    }
  }
  // Apply General Settings to all TH1D_RATIO_
  for (auto& s : TH1D_RATIO_) { for (auto& t : s.second) { for (auto& v : t.second) { if (v.second) {
          v.second->GetXaxis()->SetTitle("");
          v.second->SetTitle("");
          v.second->GetYaxis()->CenterTitle(kTRUE);
          v.second->GetYaxis()->SetTitleOffset(0.4);
          v.second->GetYaxis()->SetTitleSize(0.15);
          v.second->GetYaxis()->SetLabelSize(0.13);
          v.second->GetYaxis()->SetNdivisions(204);
          v.second->GetYaxis()->SetTitle("");
          v.second->GetXaxis()->SetTitleOffset(1);
          v.second->GetXaxis()->SetTitleSize(0.15);
          v.second->GetXaxis()->SetLabelSize(0.15);
          v.second->GetYaxis()->SetRangeUser(0.0, 2.0);
          v.second->SetMarkerColor(kBlue);
        }
      }
    }
  }
  // Apply General Settings to all TF1_RATIO
  for (auto& s : TF1_RATIO_) { for (auto& t : s.second) { for (auto& v : t.second) { if (v.second) {
          v.second->SetLineColor(kBlack);
          v.second->SetLineWidth(3);
        }
      }
    }
  }
}

void 
Histogram2::MakeStackHistogram(const std::string& tag)
{
  std::map< std::string ,std::map< std::string ,std::map< std::string , int> > > iSH;
  for (auto& s : TH1D_) {
    if (s.first.find("MC_")==std::string::npos) continue;
    for (auto& t : s.second) {
      for (auto& v : t.second) {
        if (v.second && v.second->GetEntries()>0) {
          const std::string sample  = "MC_" + getBeam(s.first);
          const std::string type    = t.first;
          const std::string varName = v.first;
          const std::string name    = "sh_" + (type!="" ? (sample + "_" + type) : sample) + "_" + varName;
          if (THStack_[sample][type][varName]==NULL) {
            THStack_[sample][type][varName] = new THStack(name.c_str(), name.c_str());
            iSH[sample][type][varName] = 0;
          }
          v.second->SetMarkerColor(COLOR[iSH.at(sample).at(type).at(varName)]);
          v.second->SetFillColor(COLOR[iSH.at(sample).at(type).at(varName)]);
          v.second->SetMarkerStyle(MARKER[0]);
          v.second->SetFillStyle(3001);
          THStack_.at(sample).at(type).at(varName)->Add(v.second);
          iSH.at(sample).at(type).at(varName)++;
        }
      }
    }
  }
  // Apply General Settings to all THStack
  for (auto& s : THStack_) { for (auto& t : s.second) { for (auto& v : t.second) { if (v.second) {
          v.second->Draw();
          v.second->GetYaxis()->SetTitleOffset(1.7);
          v.second->GetYaxis()->SetTitleSize(0.04);
          v.second->GetYaxis()->SetLabelSize(0.04);
          v.second->GetXaxis()->SetTitle(((TH1D*)v.second->GetHists()->At(0))->GetXaxis()->GetTitle());
          v.second->GetYaxis()->SetTitle(((TH1D*)v.second->GetHists()->At(0))->GetYaxis()->GetTitle());
          if (v.second->GetMaximum()<=1.5) { v.second->GetYaxis()->SetTitle("Normalized Entries"); }
          v.second->GetXaxis()->SetTitleOffset(1.0);
          v.second->GetXaxis()->SetTitleSize(0.048);
          v.second->GetXaxis()->SetLabelSize(0.035);
          setYRange(v.second, NULL, drawInfo_.begin()->second.at(t.first).at(v.first).yLogScale);
          if (tag=="DATAvsMCStack" || TH1D_RATIO_.size()>0) {
            v.second->GetYaxis()->SetTitleSize(0.04*(1./0.8));
            v.second->GetYaxis()->SetLabelSize(0.04*(1./0.8));
            v.second->GetXaxis()->SetTitleOffset(3);
            v.second->GetXaxis()->SetLabelOffset(3);
          }
        }
      }
    }
  }
}

void 
Histogram2::PrintMean(void)
{
  std::vector< std::string > sampleV;
  std::vector< std::string > typeV;
  std::vector< std::string > varNameV;
  for (auto& s : TH1D_) { sampleV.push_back(s.first); }
  for (auto& t : TH1D_[sampleV[0]]) { typeV.push_back(t.first); }
  for (auto& e : TH1D_[sampleV[0]][typeV[0]]) { varNameV.push_back(e.first); }
  for (auto& varName : varNameV) {
    std::cout << "  " << std::endl;
    std::cout << " Histrogram Statistics " << std::endl;
    std::cout << "  " << std::endl;
    for (auto& type : typeV) {
      for (auto& sample : sampleV) {
        if (TH1D_[sample][type][varName]) {
          double maxY = TH1D_[sample][type][varName]->GetBinContent(TH1D_[sample][type][varName]->GetMaximumBin());
          std::cout << (((type!="")?(sample+"_"+type):sample)+"_"+varName) 
                    << " has Norm:  " << maxY << "+-" << 0.0 
                    << " has Mean:  " << TH1D_[sample][type][varName]->GetMean()    << "+-" << TH1D_[sample][type][varName]->GetMeanError()
                    << " has Sigma: " << TH1D_[sample][type][varName]->GetStdDev()  << "+-" << TH1D_[sample][type][varName]->GetStdDevError()
                    << std::endl;
        }
      }
    }
    std::cout << "  " << std::endl;
    std::cout << "  " << std::endl;
  }
}

std::string 
Histogram2::GetString(const TH1D& h, bool Gauss)
{
  std::string out = "";
  if (Gauss) {
    double maxY = h.GetBinContent(h.GetMaximumBin());
    double mean = h.GetMean();
    double sigm = h.GetStdDev();;
    double maxYErr = (sqrt(maxY*h.Integral())/h.Integral());
    double meanErr = h.GetMeanError();
    double sigmErr = h.GetStdDevError();;
    out = Form("N: %.4f#pm%.4f , #mu: %.4f#pm%.4f , #sigma: %.4f#pm%.4f", 
               maxY , maxYErr ,
               mean , meanErr ,
               sigm , sigmErr);
  }
  else {
    double numEntries = h.GetSumOfWeights();
    double numEntriesErr = sqrt(numEntries);
    double mean = h.GetMean();
    double meanErr = h.GetMeanError();
    out = Form("Entries: %.3f#pm%.3f , Mean: %.3f#pm%.3f", 
               numEntries , numEntriesErr ,
               mean , meanErr);
  }
  return out;
}

std::string 
Histogram2::GetString(const TF1& f)
{
  std::string out = "";
  double maxY = f.GetParameter(0);
  double mean = f.GetParameter(1);
  double sigm = f.GetParameter(2);
  double maxYErr = f.GetParError(0);
  double meanErr = f.GetParError(1);
  double sigmErr = f.GetParError(2);
  out = Form("N: %.4f#pm%.4f , #mu: %.4f#pm%.4f , #sigma: %.4f#pm%.4f", 
             maxY , maxYErr ,
             mean , meanErr ,
             sigm , sigmErr);
  return out;
}

void 
Histogram2::PrintFinalMean(void)
{
  std::vector< std::string > sampleV;
  std::vector< std::string > typeV;
  std::vector< std::string > varNameV;
  for (auto& s : TH1D_) { sampleV.push_back(s.first); }
  for (auto& t : TH1D_[sampleV[0]]) { typeV.push_back(t.first); }
  for (auto& e : TH1D_[sampleV[0]][typeV[0]]) { varNameV.push_back(e.first); }
  for (auto& varName : varNameV) {
    double n_norm_Data=0, n_mean_Data=0, n_sigma_Data=0, n_norm_MC=0, n_mean_MC=0, n_sigma_MC=0;
    double norm_Data=0, mean_Data=0, sigma_Data=0, norm_MC=0, mean_MC=0, sigma_MC=0;
    std::cout << "  " << std::endl;
    std::cout << " Final Results " << std::endl;
    std::cout << "  " << std::endl;
    for (auto& type : typeV) {
      for (auto& sample : sampleV) {
        if (TH1D_[sample][type][varName]) {
          if (sample.find("DATA")!=std::string::npos) {
            norm_Data  += TF1_[sample][type][varName]->GetParameter(0)*((TF1_[sample][type][varName]->GetParError(0)>0.) ? TF1_[sample][type][varName]->GetParError(0) : 1.0);
            mean_Data  += TF1_[sample][type][varName]->GetParameter(1)*((TF1_[sample][type][varName]->GetParError(1)>0.) ? TF1_[sample][type][varName]->GetParError(1) : 1.0);
            sigma_Data += TF1_[sample][type][varName]->GetParameter(2)*((TF1_[sample][type][varName]->GetParError(2)>0.) ? TF1_[sample][type][varName]->GetParError(2) : 1.0);
            n_norm_Data  += ((TF1_[sample][type][varName]->GetParError(0)>0.) ? TF1_[sample][type][varName]->GetParError(0) : 1.0);;
            n_mean_Data  += ((TF1_[sample][type][varName]->GetParError(1)>0.) ? TF1_[sample][type][varName]->GetParError(1) : 1.0);;
            n_sigma_Data += ((TF1_[sample][type][varName]->GetParError(2)>0.) ? TF1_[sample][type][varName]->GetParError(2) : 1.0);;
          }
          if (sample.find("MC")!=std::string::npos) {
            norm_MC  += TF1_[sample][type][varName]->GetParameter(0)*((TF1_[sample][type][varName]->GetParError(0)>0.) ? TF1_[sample][type][varName]->GetParError(0) : 1.0);
            mean_MC  += TF1_[sample][type][varName]->GetParameter(1)*((TF1_[sample][type][varName]->GetParError(1)>0.) ? TF1_[sample][type][varName]->GetParError(1) : 1.0);
            sigma_MC += TF1_[sample][type][varName]->GetParameter(2)*((TF1_[sample][type][varName]->GetParError(2)>0.) ? TF1_[sample][type][varName]->GetParError(2) : 1.0);
            n_norm_MC  += ((TF1_[sample][type][varName]->GetParError(0)>0.) ? TF1_[sample][type][varName]->GetParError(0) : 1.0);
            n_mean_MC  += ((TF1_[sample][type][varName]->GetParError(1)>0.) ? TF1_[sample][type][varName]->GetParError(1) : 1.0);
            n_sigma_MC += ((TF1_[sample][type][varName]->GetParError(2)>0.) ? TF1_[sample][type][varName]->GetParError(2) : 1.0);
          }
        }
      }
      std::cout << ("DATA_"+varName) 
                << " has Norm:  " << (norm_Data/n_norm_Data) << "+-" << 0.0 
                << " has Mean:  " << (mean_Data/n_mean_Data)    << "+-" << 0.0
                << " has Sigma: " << (sigma_Data/n_sigma_Data)  << "+-" << 0.0
                << std::endl;
      std::cout << "  " << std::endl;
      std::cout << ("MC_"+varName) 
                << " has Norm:  " << (norm_MC/n_norm_MC) << "+-" << 0.0 
                << " has Mean:  " << (mean_MC/n_mean_MC)    << "+-" << 0.0
                << " has Sigma: " << (sigma_MC/n_sigma_MC)  << "+-" << 0.0
                << std::endl;
      std::cout << "  " << std::endl;
      std::cout << "  " << std::endl;
    }
  }
}

void 
Histogram2::PrintFitResults(void)
{
  std::vector< std::string > sampleV;
  std::vector< std::string > typeV;
  std::vector< std::string > varNameV;
  for (auto& s : TH1D_) { sampleV.push_back(s.first); }
  for (auto& t : TH1D_[sampleV[0]]) { typeV.push_back(t.first); }
  for (auto& e : TH1D_[sampleV[0]][typeV[0]]) { varNameV.push_back(e.first); }
  for (auto& varName : varNameV) {
    std::cout << " Fitting Results " << std::endl;
    std::cout << "  " << std::endl;
    for (auto& type : typeV) {
      for (auto& sample : sampleV) {
        if (TF1_[sample][type][varName]) {
          std::cout << (((type!="")?(sample+"_"+type):sample)+"_"+varName) 
                    << " has Norm:  " << TF1_[sample][type][varName]->GetParameter(0) << "+-" << TF1_[sample][type][varName]->GetParError(0) 
                    << " has Mean:  " << TF1_[sample][type][varName]->GetParameter(1) << "+-" << TF1_[sample][type][varName]->GetParError(1) 
                    << " has Sigma: " << TF1_[sample][type][varName]->GetParameter(2) << "+-" << TF1_[sample][type][varName]->GetParError(2)
                    << std::endl;
        }
      }
    }
    std::cout << "  " << std::endl;
    std::cout << "  " << std::endl;
  }
}

void 
Histogram2::Draw(const std::string& tag, const std::map< std::string , struct TypeInfo >& typeInfo)
{
  if (tag=="") return;

  // set the CMS style
  setTDRStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gSystem->mkdir("Plots", kTRUE);

  // Apply General Settings to all TH1D
  for (auto& s : TH1D_) { for (auto& t : s.second) { for (auto& v : t.second) { if (v.second) {
          v.second->GetYaxis()->SetTitleOffset(1.7);
          v.second->GetYaxis()->SetTitleSize(0.04);
          v.second->GetYaxis()->SetLabelSize(0.04);
          if (v.second->GetSumOfWeights()<=1.5) { v.second->GetYaxis()->SetTitle("Normalized Entries"); }
          v.second->GetXaxis()->SetTitleOffset(1.0);
          v.second->GetXaxis()->SetTitleSize(0.048);
          v.second->GetXaxis()->SetLabelSize(0.035);
          v.second->SetMarkerColor(kBlack);
          v.second->SetFillStyle(3001);
          v.second->SetFillColor(kBlue);
          setYRange(v.second, NULL, drawInfo_.at(s.first).at(t.first).at(v.first).yLogScale);
          if (tag=="DATAvsMC" || tag=="DATAvsFIT" || tag!="DATAvsMCStack" || TH1D_RATIO_.size()>0) {
            v.second->GetYaxis()->SetTitleSize(0.04*(1./0.8));
            v.second->GetYaxis()->SetLabelSize(0.04*(1./0.8));
            v.second->GetXaxis()->SetTitleOffset(3);
            v.second->GetXaxis()->SetLabelOffset(3);
          }
        }
      }
    }
  }

  // Case: Separate With Data -> Use one canvas for each histogram
  if (tag=="separate") {
    Double_t xl1=.20, yl1=0.62, xl2=xl1+.3, yl2=yl1+.250;
    for (const auto& s : TH1D_) {
      const std::string sample = s.first;
      for (const auto& t : s.second) {
        const std::string type = t.first;
        gSystem->mkdir(Form("Plots/%s/%s/%s/png/", tag.c_str(), type.c_str(), sample.c_str()), kTRUE);
        gSystem->mkdir(Form("Plots/%s/%s/%s/pdf/", tag.c_str(), type.c_str(), sample.c_str()), kTRUE);
        for (const auto& v : t.second) {
          const std::string varName = v.first;
          if (v.second) {
            std::string cName = (std::string("c_") + (type!="" ? (sample + "_" + type) : sample) + "_" + varName);
            TCanvas* c = new TCanvas(cName.c_str(), cName.c_str(), 1000, 1000);
            TLegend *leg = new TLegend(xl1,yl1,xl2,yl2);
            c->cd();
            if (v.second->GetEntries()>0) {
              c->SetLogy(drawInfo_.at(sample).at(type).at(varName).yLogScale);
              v.second->Draw("p");
              leg->AddEntry(v.second, (type!="" ? (sample+"_"+type) : sample).c_str(), "p");
              if (TF1_.count(sample)>0 && TF1_.at(sample).count(type)>0 && TF1_.at(sample).at(type).count(varName)>0 && TF1_.at(sample).at(type).at(varName)) {
                leg->AddEntry(TF1_.at(sample).at(type).at(varName), "Fit", "l");
              }
              leg->Draw("SAME");
              c->Update();
              int option = 111;
              if (sample.find("pPb")!=std::string::npos) option = 109;
              if (sample.find("Pbp")!=std::string::npos) option = 110;
              CMS_lumi(c, option, 33, "");
              c->Update();
              c->SaveAs(Form("Plots/%s/%s/%s/png/%s.png", tag.c_str(), type.c_str(), sample.c_str(), cName.c_str()));
              c->SaveAs(Form("Plots/%s/%s/%s/pdf/%s.pdf", tag.c_str(), type.c_str(), sample.c_str(), cName.c_str()));
            }
            c->Clear();
            c->Close();
            delete c;
            delete leg;
          }
        }
      }
    }
    return;
  }
  // Case: Data vs StackMC -> Use one canvas for each stack histogram comparing MC vs Data
  if (tag=="DATAvsMCStack") {
    MakeStackHistogram("DATAvsMCStack");
    ComputeRatio("DATAvsMCStack");
    Double_t xl1=.20, yl1=0.75, xl2=xl1+.3, yl2=yl1+.135;
    for (const auto& s : THStack_) {
      const std::string sample = s.first;
      for (auto& t : s.second) {
        const std::string type = t.first;
        gSystem->mkdir(Form("Plots/%s/%s/%s/png/", tag.c_str(), type.c_str(), sample.c_str()), kTRUE);
        gSystem->mkdir(Form("Plots/%s/%s/%s/pdf/", tag.c_str(), type.c_str(), sample.c_str()), kTRUE);
        for (auto& v : t.second) {
          const std::string varName = v.first;
          if (v.second) {
            std::string cName = (std::string("c_") + "DATAvsMCStack" + "_" + (type!="" ? (sample + "_" + type) : sample) + "_" + varName);
            TCanvas* c = new TCanvas(cName.c_str(), cName.c_str(), 1000, 1000);
            TLegend *leg = new TLegend(xl1,yl1,xl2,yl2);
            TLatex *tex = new TLatex(); tex->SetNDC(); tex->SetTextSize(0.025); float dy = 0; 
            TPad* padRATIO = NULL; TPad* padMAIN = NULL;
            bool draw = false;
            TH1D* hStack = ((TH1D*)v.second->GetStack()->At(v.second->GetNhists()-1));
            c->cd();
            if (v.second->GetMaximum()>0.) {
              padMAIN = NULL;
              if (TH1D_RATIO_.count(sample)>0 && TH1D_RATIO_.at(sample).count(type)>0 && TH1D_RATIO_.at(sample).at(type).count(varName)>0 && TH1D_RATIO_.at(sample).at(type).at(varName)) {
                padMAIN = new TPad( "padMAIN", "", 0, 0.2, 1, 1 );
                padMAIN->SetFixedAspectRatio(kTRUE);
                padMAIN->SetBottomMargin(0.015);
              }
              else { padMAIN = new TPad( "padMAIN", "", 0, 0, 1, 1 ); }
              padMAIN->Draw();
              padMAIN->cd();
              padMAIN->SetLogy(drawInfo_.begin()->second.at(type).at(varName).yLogScale);
              v.second->Draw("HISTF");
              std::map<std::string , double> nEntries;
              for (int i=0; i<v.second->GetNhists(); i++) {
                std::string name = ((TH1D*)v.second->GetHists()->At(i))->GetTitle();
                std::string tmp = name.substr(name.find("MC_")+3, name.length()); tmp = tmp.substr(0, tmp.find("_"));
                name = "MC_" + tmp + "_" + getBeam(name);
                leg->AddEntry(v.second->GetHists()->At(i), name.c_str(), "f");
                nEntries[tmp] = ((TH1D*)v.second->GetHists()->At(i))->GetSumOfWeights();
              }
              if (typeInfo.count(type)>0) {
                double dy = 0.08; for (const auto& s: typeInfo.at(type).cutSelection) { tex->DrawLatex(0.66, 0.72-dy, s.c_str()); dy+=0.04; }
              }
              // Plot Data
              std::string sampleData = "DATA_" + getBeam(sample);
              tex->DrawLatex(0.20, 0.72-dy, ("MC    : " + GetString(*hStack)).c_str()); dy+=0.040;
              if (TH1D_.count(sampleData)>0 && TH1D_.at(sampleData).count(type)>0 && TH1D_.at(sampleData).at(type).count(varName)>0 && 
                  TH1D_.at(sampleData).at(type).at(varName) && TH1D_.at(sampleData).at(type).at(varName)->GetEntries()>0.) {
                TH1D_.at(sampleData).at(type).at(varName)->SetMarkerColor(kRed); TH1D_.at(sampleData).at(type).at(varName)->Draw("samep"); 
                leg->AddEntry(TH1D_.at(sampleData).at(type).at(varName), sampleData.c_str(), "p");
                tex->DrawLatex(0.20, 0.72-dy, ("DATA: " + GetString(*TH1D_.at(sampleData).at(type).at(varName))).c_str()); dy+=0.040;
                setYRange(v.second, TH1D_.at(sampleData).at(type).at(varName), drawInfo_.at(sampleData).at(type).at(varName).yLogScale);
                double ratio = TH1D_.at(sampleData).at(type).at(varName)->GetSumOfWeights()/hStack->GetSumOfWeights();
                tex->DrawLatex(0.20, 0.72-dy, Form("DATA/MC: %.4f", ratio));
                if (typeInfo.count(type)>0) {
                  std::vector<std::string> typeSample = typeInfo.at(type).sample;
                  double bkgN = 0.; for (const auto& n : nEntries) { bool isSig=false; for (const auto& s : typeSample) { if (n.first==s) {isSig=true; break;} } if (!isSig) bkgN += n.second; }
                  double sigN = 0.; for (const auto& s : typeSample) { if (nEntries.count(s)>0) sigN += nEntries.at(s); }
                  if (sigN>0.) {
                    double ratio = ( (TH1D_.at(sampleData).at(type).at(varName)->GetSumOfWeights() - bkgN) / sigN );
                    tex->DrawLatex(0.37, 0.72-dy, Form("%s: %.4f", typeSample.at(0).c_str(), ratio));
                  }
                }                  
              }
              leg->Draw("SAME");
              padMAIN->Modified();
              padMAIN->Update();
              int option = 111;
              if (sample.find("pPb")!=std::string::npos) option = 109;
              if (sample.find("Pbp")!=std::string::npos) option = 110;
              CMS_lumi(padMAIN, option, 33, "");
              padMAIN->Modified();
              padMAIN->Update();
              draw = true;
            }
            c->cd();
            if (draw) {
              if (sample.find("MC")!=std::string::npos) {
                padRATIO = new TPad( "padRATIO", "", 0, 0, 1, 0.20 );      
                padRATIO->SetFixedAspectRatio(kTRUE);
                padRATIO->SetTopMargin(0.02);
                padRATIO->SetBottomMargin(0.4);
                padRATIO->SetFillStyle(4000);
                padRATIO->SetFrameFillStyle(4000);
                padRATIO->SetGridx(kTRUE);
                padRATIO->SetGridy(kTRUE);
                padRATIO->Draw();
                padRATIO->cd();
                if (TH1D_RATIO_.count(sample)>0 && TH1D_RATIO_.at(sample).count(type)>0 && TH1D_RATIO_.at(sample).at(type).count(varName)>0 && TH1D_RATIO_.at(sample).at(type).at(varName)) {
                  TH1D_RATIO_.at(sample).at(type).at(varName)->GetYaxis()->SetTitle("#frac{DATA}{MC}");
                  TH1D_RATIO_.at(sample).at(type).at(varName)->Draw("p");
                  if (TF1_RATIO_[sample][type][varName]!=NULL) TF1_RATIO_[sample][type][varName]->Draw("samel");
                }
                padRATIO->Modified();
                padRATIO->Update();
              }
              c->Modified();
              c->Update();
              c->SaveAs(Form("Plots/%s/%s/%s/png/%s.png", tag.c_str(), type.c_str(), sample.c_str(), cName.c_str()));
              c->SaveAs(Form("Plots/%s/%s/%s/pdf/%s.pdf", tag.c_str(), type.c_str(), sample.c_str(), cName.c_str()));
            }
            c->Clear();
            c->Close();
            delete c;
            delete leg;
            delete tex;
          }
        }
      }
    }
    return;
  }
  // Case: Data vs Fit -> Use one canvas for each histogram comparing Fit vs Data
  if (tag=="DATAvsFIT") {
    ComputeRatio("DATAvsFIT");
    Double_t xl1=.20, yl1=0.75, xl2=xl1+.3, yl2=yl1+.125;
    for (const auto& s : TH1D_) {
      const std::string sample = s.first;
      for (auto& t : s.second) {
        const std::string type = t.first;
        gSystem->mkdir(Form("Plots/%s/%s/%s/png/", tag.c_str(), type.c_str(), sample.c_str()), kTRUE);
        gSystem->mkdir(Form("Plots/%s/%s/%s/pdf/", tag.c_str(), type.c_str(), sample.c_str()), kTRUE);
        for (auto& v : t.second) {
          const std::string varName = v.first;
          if (v.second) {
            std::string cName = (std::string("c_") + "DATAvsFIT" + "_" + (type!="" ? (sample + "_" + type) : sample) + "_" + varName);
            TCanvas* c = new TCanvas(cName.c_str(), cName.c_str(), 1000, 1000);
            TLegend *leg = new TLegend(xl1,yl1,xl2,yl2);
            TLatex *tex = new TLatex(); tex->SetNDC(); tex->SetTextSize(0.025); float dy = 0; 
            TPad* padRATIO = NULL; TPad* padMAIN = NULL;
            bool draw = false;
            c->cd();
            if (v.second->GetEntries()>0.) {
              padMAIN = NULL;
              if (TH1D_RATIO_.count(sample)>0 && TH1D_RATIO_.at(sample).count(type)>0 && TH1D_RATIO_.at(sample).at(type).count(varName)>0 && TH1D_RATIO_.at(sample).at(type).at(varName)) {
                padMAIN = new TPad( "padMAIN", "", 0, 0.2, 1, 1 );
                padMAIN->SetFixedAspectRatio(kTRUE);
                padMAIN->SetBottomMargin(0.015);
              }
              else { padMAIN = new TPad( "padMAIN", "", 0, 0, 1, 1 ); }
              padMAIN->Draw();
              padMAIN->cd();
              padMAIN->SetLogy(drawInfo_.at(sample).at(type).at(varName).yLogScale);
              v.second->Draw("p");
              leg->AddEntry(v.second, (type!="" ? (sample+"_"+type) : sample).c_str(), "p");
              tex->DrawLatex(0.20, 0.72-dy, ("HIST: " + GetString(*v.second)).c_str()); dy+=0.040;
              if (TF1_.count(sample)>0 && TF1_.at(sample).count(type)>0 && TF1_.at(sample).at(type).count(varName)>0 && TF1_.at(sample).at(type).at(varName)) {
                TF1_.at(sample).at(type).at(varName)->Draw("samel");
                if (sample.find("DATA")!=std::string::npos) { leg->AddEntry(TF1_.at(sample).at(type).at(varName), "Fit DATA", "l"); }
                if (sample.find("MC")!=std::string::npos  ) { leg->AddEntry(TF1_.at(sample).at(type).at(varName), "Fit MC", "l");   }
                tex->DrawLatex(0.20, 0.72-dy, ("FIT : " + GetString(*TF1_.at(sample).at(type).at(varName))).c_str()); dy+=0.040;
              }
              leg->Draw("SAME");
              padMAIN->Update();
              int option = 111;
              if (sample.find("pPb")!=std::string::npos) option = 109;
              if (sample.find("Pbp")!=std::string::npos) option = 110;
              CMS_lumi(padMAIN, option, 33, "");
              padMAIN->Update();
              draw = true;
            }
            c->cd();
            if (draw) {
              padRATIO = new TPad( "padRATIO", "", 0, 0, 1, 0.20 );      
              padRATIO->SetFixedAspectRatio(kTRUE);
              padRATIO->SetTopMargin(0.02);
              padRATIO->SetBottomMargin(0.4);
              padRATIO->SetFillStyle(4000);
              padRATIO->SetFrameFillStyle(4000);
              padRATIO->SetGridx(kTRUE);
              padRATIO->SetGridy(kTRUE);
              padRATIO->Draw();
              padRATIO->cd();
              if (TH1D_RATIO_.count(sample)>0 && TH1D_RATIO_.at(sample).count(type)>0 && TH1D_RATIO_.at(sample).at(type).count(varName)>0 && TH1D_RATIO_.at(sample).at(type).at(varName)) {
                TH1D_RATIO_[sample][type][varName]->GetYaxis()->SetTitle("#frac{DATA}{FIT}");
                TH1D_RATIO_[sample][type][varName]->Draw("p");
              }
              padRATIO->Update();
              c->Update();
              c->SaveAs(Form("Plots/%s/%s/%s/png/%s.png", tag.c_str(), type.c_str(), sample.c_str(), cName.c_str()));
              c->SaveAs(Form("Plots/%s/%s/%s/pdf/%s.pdf", tag.c_str(), type.c_str(), sample.c_str(), cName.c_str()));
            }
            c->Clear();
            c->Close();
            delete c;
            delete leg;
            delete tex;
          }
        }
      }
    }
    return;
  }
  // Case: Data vs MC -> Use one canvas for each histogram comparing MC vs Data
  if (tag=="DATAvsMC") {
    ComputeRatio("DATAvsMC");
    Double_t xl1=.20, yl1=0.75, xl2=xl1+.3, yl2=yl1+.125;
    for (const auto& s : TH1D_) {
      const std::string sample = s.first;
      if ( tag=="DATAvsMC" && sample.find("MC")==std::string::npos) continue;
      for (auto& t : s.second) {
        const std::string type = t.first;
        gSystem->mkdir(Form("Plots/%s/%s/%s/png/", tag.c_str(), type.c_str(), sample.c_str()), kTRUE);
        gSystem->mkdir(Form("Plots/%s/%s/%s/pdf/", tag.c_str(), type.c_str(), sample.c_str()), kTRUE);
        for (auto& v : t.second) {
          const std::string varName = v.first;
          if (v.second) {
            std::string cName = (std::string("c_") + "DATAvsMC" + "_" + (type!="" ? (sample + "_" + type) : sample) + "_" + varName);
            TCanvas* c = new TCanvas(cName.c_str(), cName.c_str(), 1000, 1000);
            TLegend *leg = new TLegend(xl1,yl1,xl2,yl2);
            TLatex *tex = new TLatex(); tex->SetNDC(); tex->SetTextSize(0.025); float dy = 0; 
            TPad* padRATIO = NULL; TPad* padMAIN = NULL;
            bool draw = false;
            c->cd();
            if (v.second->GetEntries()>0.) {
              padMAIN = NULL;
              if (TH1D_RATIO_.count(sample)>0 && TH1D_RATIO_.at(sample).count(type)>0 && TH1D_RATIO_.at(sample).at(type).count(varName)>0 && TH1D_RATIO_.at(sample).at(type).at(varName)) {
                padMAIN = new TPad( "padMAIN", "", 0, 0.2, 1, 1 );
                padMAIN->SetFixedAspectRatio(kTRUE);
                padMAIN->SetBottomMargin(0.015);
              }
              else { padMAIN = new TPad( "padMAIN", "", 0, 0, 1, 1 ); }
              padMAIN->Draw();
              padMAIN->cd();
              padMAIN->SetLogy(drawInfo_.at(sample).at(type).at(varName).yLogScale);
              v.second->Draw("histf");
              leg->AddEntry(v.second, (type!="" ? (sample+"_"+type) : sample).c_str(), "f");
              if (typeInfo.count(type)>0) {
                double dy = 0.08; for (const auto& s: typeInfo.at(type).cutSelection) { tex->DrawLatex(0.66, 0.72-dy, s.c_str()); dy+=0.04; }
              }
              if (TF1_.count(sample)>0 && TF1_.at(sample).count(type)>0 && TF1_.at(sample).at(type).count(varName)>0 && TF1_.at(sample).at(type).at(varName)) {
                TF1_.at(sample).at(type).at(varName)->Draw("samel");
                if (sample.find("DATA")!=std::string::npos) { leg->AddEntry(TF1_.at(sample).at(type).at(varName), "Fit DATA", "l"); }
                if (sample.find("MC")!=std::string::npos  ) { leg->AddEntry(TF1_.at(sample).at(type).at(varName), "Fit MC", "l");   }
                tex->DrawLatex(0.20, 0.72-dy, ("FIT : " + GetString(*TF1_.at(sample).at(type).at(varName))).c_str()); dy+=0.040;
              }
              // Plot Data
              std::string sampleData = "DATA_" + getBeam(sample);
              if (TF1_[sample][type][varName]==NULL) { tex->DrawLatex(0.20, 0.72-dy, ("MC    : " + GetString(*v.second)).c_str()); dy+=0.040; }
              if (TH1D_.count(sampleData)>0 && TH1D_.at(sampleData).count(type)>0 && TH1D_.at(sampleData).at(type).count(varName)>0 && 
                  TH1D_.at(sampleData).at(type).at(varName) && TH1D_.at(sampleData).at(type).at(varName)->GetEntries()>0.) {
                TH1D_.at(sampleData).at(type).at(varName)->SetMarkerColor(kRed); TH1D_.at(sampleData).at(type).at(varName)->Draw("samep"); 
                leg->AddEntry(TH1D_.at(sampleData).at(type).at(varName), sampleData.c_str(), "p");
                if (TF1_[sampleData][type][varName]==NULL) { tex->DrawLatex(0.20, 0.72-dy, ("DATA: " + GetString(*TH1D_[sampleData][type][varName])).c_str()); dy+=0.040; }
                setYRange(v.second, TH1D_.at(sampleData).at(type).at(varName), drawInfo_.at(sample).at(type).at(varName).yLogScale);
                double ratio = TH1D_.at(sampleData).at(type).at(varName)->GetSumOfWeights()/v.second->GetSumOfWeights();
                tex->DrawLatex(0.20, 0.72-dy, Form("DATA/MC: %.3f", ratio));
              }
              if (TF1_[sampleData][type][varName]!=NULL) { 
                TF1_[sampleData][type][varName]->SetLineColor(kRed); TF1_[sampleData][type][varName]->Draw("samel"); leg->AddEntry(TF1_[sampleData][type][varName], "Fit Data", "l");
                tex->DrawLatex(0.20, 0.72-dy, ("DATA: " + GetString(*TF1_[sampleData][type][varName])).c_str()); dy+=0.040;
              }
              leg->Draw("SAME");
              padMAIN->Update();
              int option = 111;
              if (sample.find("pPb")!=std::string::npos) option = 109;
              if (sample.find("Pbp")!=std::string::npos) option = 110;
              CMS_lumi(padMAIN, option, 33, "");
              padMAIN->Update();
              draw = true;
            }
            c->cd();
            if (draw) {
              if (sample.find("MC")!=std::string::npos) {
                padRATIO = new TPad( "padRATIO", "", 0, 0, 1, 0.20 );      
                padRATIO->SetFixedAspectRatio(kTRUE);
                padRATIO->SetTopMargin(0.02);
                padRATIO->SetBottomMargin(0.4);
                padRATIO->SetFillStyle(4000);
                padRATIO->SetFrameFillStyle(4000);
                padRATIO->SetGridx(kTRUE);
                padRATIO->SetGridy(kTRUE);
                padRATIO->Draw();
                padRATIO->cd();
                if (TH1D_RATIO_.count(sample)>0 && TH1D_RATIO_.at(sample).count(type)>0 && TH1D_RATIO_.at(sample).at(type).count(varName)>0 && TH1D_RATIO_.at(sample).at(type).at(varName)) {
                  TH1D_RATIO_.at(sample).at(type).at(varName)->GetYaxis()->SetTitle("#frac{DATA}{MC}");
                  TH1D_RATIO_.at(sample).at(type).at(varName)->Draw("p");
                  if (TF1_RATIO_[sample][type][varName]!=NULL) TF1_RATIO_[sample][type][varName]->Draw("samel");
                }
                padRATIO->Update();
              }
              c->Update();
              c->SaveAs(Form("Plots/%s/%s/%s/png/%s.png", tag.c_str(), type.c_str(), sample.c_str(), cName.c_str()));
              c->SaveAs(Form("Plots/%s/%s/%s/pdf/%s.pdf", tag.c_str(), type.c_str(), sample.c_str(), cName.c_str()));
            }
            c->Clear();
            c->Close();
            delete c;
            delete leg;
            delete tex;
          }
        }
      }
    }
    return;
  }
  // Case: MC_NLO vs MC_LO -> Use one canvas for each histogram comparing MC_NLO vs MC_LO
  if (tag=="MCNLOvsMCLO") {
    ComputeRatio("MCNLOvsMCLO");
    Double_t xl1=.20, yl1=0.75, xl2=xl1+.3, yl2=yl1+.125;
    for (const auto& s : TH1D_) {
      const std::string sample = s.first;
      if ( tag=="MCNLOvsMCLO" && sample.find("MCNLO")==std::string::npos) continue;
      for (auto& t : s.second) {
        const std::string type = t.first;
        gSystem->mkdir(Form("Plots/%s/%s/%s/png/", tag.c_str(), type.c_str(), sample.c_str()), kTRUE);
        gSystem->mkdir(Form("Plots/%s/%s/%s/pdf/", tag.c_str(), type.c_str(), sample.c_str()), kTRUE);
        for (auto& v : t.second) {
          const std::string varName = v.first;
          if (v.second) {
            std::string cName = (std::string("c_") + "MCNLOvsMCLO" + "_" + (type!="" ? (sample + "_" + type) : sample) + "_" + varName);
            TCanvas* c = new TCanvas(cName.c_str(), cName.c_str(), 1000, 1000);
            TLegend *leg = new TLegend(xl1,yl1,xl2,yl2);
            TLatex *tex = new TLatex(); tex->SetNDC(); tex->SetTextSize(0.025); float dy = 0; 
            TPad* padRATIO = NULL; TPad* padMAIN = NULL;
            bool draw = false;
            c->cd();
            if (v.second->GetEntries()>0.) {
              padMAIN = NULL;
              if (TH1D_RATIO_.count(sample)>0 && TH1D_RATIO_.at(sample).count(type)>0 && TH1D_RATIO_.at(sample).at(type).count(varName)>0 && TH1D_RATIO_.at(sample).at(type).at(varName)) {
                padMAIN = new TPad( "padMAIN", "", 0, 0.2, 1, 1 );
                padMAIN->SetFixedAspectRatio(kTRUE);
                padMAIN->SetBottomMargin(0.015);
              }
              else { padMAIN = new TPad( "padMAIN", "", 0, 0, 1, 1 ); }
              padMAIN->Draw();
              padMAIN->cd();
              padMAIN->SetLogy(drawInfo_.at(sample).at(type).at(varName).yLogScale);
              v.second->Draw("histf");
              leg->AddEntry(v.second, sample.c_str(), "f");
              if (typeInfo.count(type)>0) {
                double dy = 0.08; for (const auto& s: typeInfo.at(type).cutSelection) { tex->DrawLatex(0.66, 0.72-dy, s.c_str()); dy+=0.04; }
              }
              // Plot Leading Order MC
              std::string sampleMCLO = "MCLO_" + getDecayChannel(sample) + getBeam(sample);
              tex->DrawLatex(0.20, 0.72-dy, ("MCNLO    : " + GetString(*v.second)).c_str()); dy+=0.040;
              if (TH1D_.count(sampleMCLO)>0 && TH1D_.at(sampleMCLO).count(type)>0 && TH1D_.at(sampleMCLO).at(type).count(varName)>0 && 
                  TH1D_.at(sampleMCLO).at(type).at(varName) && TH1D_.at(sampleMCLO).at(type).at(varName)->GetEntries()>0.) {
                TH1D_.at(sampleMCLO).at(type).at(varName)->SetMarkerColor(kRed); TH1D_.at(sampleMCLO).at(type).at(varName)->Draw("samep"); 
                leg->AddEntry(TH1D_.at(sampleMCLO).at(type).at(varName), sampleMCLO.c_str(), "p");
                tex->DrawLatex(0.20, 0.72-dy, ("MCLO: " + GetString(*TH1D_[sampleMCLO][type][varName])).c_str()); dy+=0.040;
                setYRange(v.second, TH1D_.at(sampleMCLO).at(type).at(varName), drawInfo_.at(sample).at(type).at(varName).yLogScale);
                double ratio = v.second->GetSumOfWeights()/TH1D_.at(sampleMCLO).at(type).at(varName)->GetSumOfWeights();
                tex->DrawLatex(0.20, 0.72-dy, Form("MCNLO/MCLO: %.3f", ratio));
              }
              leg->Draw("SAME");
              padMAIN->Update();
              int option = 111;
              if (sample.find("pPb")!=std::string::npos) option = 109;
              if (sample.find("Pbp")!=std::string::npos) option = 110;
              CMS_lumi(padMAIN, option, 33, "");
              padMAIN->Update();
              draw = true;
            }
            c->cd();
            if (draw) {
              if (sample.find("MCNLO")!=std::string::npos) {
                padRATIO = new TPad( "padRATIO", "", 0, 0, 1, 0.20 );      
                padRATIO->SetFixedAspectRatio(kTRUE);
                padRATIO->SetTopMargin(0.02);
                padRATIO->SetBottomMargin(0.4);
                padRATIO->SetFillStyle(4000);
                padRATIO->SetFrameFillStyle(4000);
                padRATIO->SetGridx(kTRUE);
                padRATIO->SetGridy(kTRUE);
                padRATIO->Draw();
                padRATIO->cd();
                if (TH1D_RATIO_.count(sample)>0 && TH1D_RATIO_.at(sample).count(type)>0 && TH1D_RATIO_.at(sample).at(type).count(varName)>0 && TH1D_RATIO_.at(sample).at(type).at(varName)) {
                  TH1D_RATIO_.at(sample).at(type).at(varName)->GetYaxis()->SetTitle("#frac{MC[NLO]}{MC[LO]}");
                  TH1D_RATIO_.at(sample).at(type).at(varName)->Draw("p");
                  if (TF1_RATIO_[sample][type][varName]!=NULL) TF1_RATIO_[sample][type][varName]->Draw("samel");
                }
                padRATIO->Update();
              }
              c->Update();
              c->SaveAs(Form("Plots/%s/%s/%s/png/%s.png", tag.c_str(), type.c_str(), sample.c_str(), cName.c_str()));
              c->SaveAs(Form("Plots/%s/%s/%s/pdf/%s.pdf", tag.c_str(), type.c_str(), sample.c_str(), cName.c_str()));
            }
            c->Clear();
            c->Close();
            delete c;
            delete leg;
            delete tex;
          }
        }
      }
    }
    return;
  }
  // Initialize the labels
  std::vector< std::string > sampleV;
  std::vector< std::string > typeV;
  std::vector< std::string > varNameV;
  for (auto& s : TH1D_) { sampleV.push_back(s.first); }
  for (auto& t : TH1D_[sampleV[0]]) { typeV.push_back(t.first); }
  for (auto& e : TH1D_[sampleV[0]][typeV[0]]) { varNameV.push_back(e.first); }
  // Case: join -> Draw all histograms in the same canvas for each variable
  if (tag=="join") {
    Double_t xl1=.20, yl1=0.75, xl2=xl1+.5, yl2=yl1+.065;
    std::vector<std::string> beamDir = {"pPb", "Pbp", "PA"};
    for (const auto& beam : beamDir) {
      for (const auto& varName : varNameV) {
        for (const auto& type : typeV) {      
          gSystem->mkdir(Form("Plots/%s/%s/%s/png/", tag.c_str(), type.c_str(), beam.c_str()), kTRUE);
          gSystem->mkdir(Form("Plots/%s/%s/%s/pdf/", tag.c_str(), type.c_str(), beam.c_str()), kTRUE);
          std::string cName = (std::string("c_") + beam + "_" + (type!="" ? (type + "_" + varName) : varName));
          TCanvas* c = new TCanvas(cName.c_str(), cName.c_str(), 1000, 1000);
          TLegend *leg = new TLegend(xl1,yl1,xl2,yl2);
          TLatex *tex = new TLatex(); tex->SetNDC(); tex->SetTextSize(0.025); float dy = 0; 
          TPad* padRATIO = NULL; TPad* padMAIN = NULL;
          bool draw = false;
          c->cd();
          bool found = false;
          if (TH1D_.size()>0) {
            padMAIN = NULL;
            if (TH1D_RATIO_.size()>0) { 
              padMAIN = new TPad( "padMAIN", "", 0, 0.2, 1, 1 );
              padMAIN->SetFixedAspectRatio(kTRUE);
              padMAIN->SetBottomMargin(0.015);
            }
            else { padMAIN = new TPad( "padMAIN", "", 0, 0, 1, 1 ); }
            padMAIN->Draw();
            padMAIN->cd();
            bool firstDraw = true;
            uint i = 0;
            std::string prevDS = "";
            for (const auto& sample : sampleV) {
              if (beam!="PA" && sample.find(beam)==std::string::npos) continue;
              if ( TH1D_.count(sample)>0 && TH1D_.at(sample).count(type)>0 && TH1D_.at(sample).at(type).count(varName)>0 && TH1D_.at(sample).at(type).at(varName) && 
                   TH1D_.at(sample).at(type).at(varName)->GetEntries()>0. ) {
                if (sample.find("DATA")!=std::string::npos) { TH1D_.at(sample).at(type).at(varName)->SetMarkerStyle(MARKER[0]); if (prevDS!="DATA") { i=0; prevDS="DATA"; } }
                if (sample.find("MC")!=std::string::npos)   { TH1D_.at(sample).at(type).at(varName)->SetMarkerStyle(MARKER[1]); if (prevDS!="MC"  ) { i=0; prevDS="MC";   } }
                TH1D_.at(sample).at(type).at(varName)->SetMarkerColor(COLOR[i]);
                if (firstDraw) { TH1D_.at(sample).at(type).at(varName)->Draw("P"); firstDraw = false; }
                else  { TH1D_.at(sample).at(type).at(varName)->Draw("SAMEP"); }
                if (TF1_[sample][type][varName]) {
                  TF1_.at(sample).at(type).at(varName)->SetLineColor(COLOR[i]);
                  if (sample.find("DATA")!=std::string::npos) { TF1_.at(sample).at(type).at(varName)->SetLineStyle(LINE[0]); }
                  if (sample.find("MC")!=std::string::npos)   { TF1_.at(sample).at(type).at(varName)->SetLineStyle(LINE[1]); }
                  TF1_.at(sample).at(type).at(varName)->Draw("SAMEL");
                }
                leg->AddEntry(TH1D_.at(sample).at(type).at(varName), sample.c_str(), "p");
                i++;
                found = true;
              }
            }
            if (typeInfo.count(type)>0) {
              double dy = 0.08; for (const auto& s: typeInfo.at(type).cutSelection) { tex->DrawLatex(0.66, 0.72-dy, s.c_str()); dy+=0.04; }
            }
          }
          if (found) {
            leg->Draw("SAME");
            padMAIN->Update();
            int option = 111;
            if (beam=="pPb") option = 109;
            if (beam=="Pbp") option = 110;
            CMS_lumi(padMAIN, option, 33, "");
            padMAIN->Update();
            draw = true;
          }
          c->cd();
          if (TH1D_RATIO_.size()>0) {
            padRATIO = new TPad( "padRATIO", "", 0, 0, 1, 0.20 );      
            padRATIO->SetFixedAspectRatio(kTRUE);
            padRATIO->SetTopMargin(0.02);
            padRATIO->SetBottomMargin(0.4);
            padRATIO->SetFillStyle(4000);
            padRATIO->SetFrameFillStyle(4000);
            padRATIO->SetGridx(kTRUE);
            padRATIO->SetGridy(kTRUE);
            padRATIO->Draw();
            padRATIO->cd();
            bool firstDraw = true;
            uint i = 0;
            std::string prevDS = "";
            for (const auto& sample : sampleV) {
              if (beam!="PA" && sample.find(beam)==std::string::npos) continue;
              if (sample.find("MC")==std::string::npos) continue;
              if ( TH1D_RATIO_.count(sample)>0 && TH1D_RATIO_.at(sample).count(type)>0 && TH1D_RATIO_.at(sample).at(type).count(varName)>0 && TH1D_RATIO_.at(sample).at(type).at(varName) && 
                   TH1D_RATIO_.at(sample).at(type).at(varName)->GetEntries()>0. ) {
                if (tag=="join") TH1D_RATIO_.at(sample).at(type).at(varName)->GetYaxis()->SetTitle("#frac{DATA}{MC}");
                if (sample.find("DATA")!=std::string::npos) { TH1D_RATIO_.at(sample).at(type).at(varName)->SetMarkerStyle(MARKER[0]); if (prevDS!="DATA") { i=0; prevDS="DATA"; }  }
                if (sample.find("MC")!=std::string::npos)   { TH1D_RATIO_.at(sample).at(type).at(varName)->SetMarkerStyle(MARKER[1]); if (prevDS!="MC"  ) { i=0; prevDS="MC";   }  }
                TH1D_RATIO_.at(sample).at(type).at(varName)->SetMarkerColor(COLOR[i]);
                if (firstDraw) { TH1D_RATIO_.at(sample).at(type).at(varName)->Draw("P"); firstDraw = false; }
                else  {
                  TH1D_RATIO_.at(sample).at(type).at(varName)->Draw("SAMEP");
                }
                if ( TF1_RATIO_.count(sample)>0 && TF1_RATIO_.at(sample).count(type)>0 && TF1_RATIO_.at(sample).at(type).count(varName)>0 && TF1_RATIO_.at(sample).at(type).at(varName) ) {
                  TF1_RATIO_.at(sample).at(type).at(varName)->SetLineColor(COLOR[i]);
                  if (sample.find("DATA")!=std::string::npos) { TF1_RATIO_.at(sample).at(type).at(varName)->SetLineStyle(LINE[0]); }
                  if (sample.find("MC")!=std::string::npos)   { TF1_RATIO_.at(sample).at(type).at(varName)->SetLineStyle(LINE[1]); }
                  TF1_RATIO_.at(sample).at(type).at(varName)->Draw("SAMEL");
                }
                i++;
              }
            }
            padRATIO->Update();
            draw = true;
          }
          if (draw) {
            c->Update();
            c->SaveAs(Form("Plots/%s/%s/%s/png/%s.png", tag.c_str(), type.c_str(), beam.c_str(), cName.c_str()));
            c->SaveAs(Form("Plots/%s/%s/%s/pdf/%s.pdf", tag.c_str(), type.c_str(), beam.c_str(), cName.c_str()));
          }
          c->Clear();
          c->Close();
          delete c;
          delete leg;
          delete tex;
        }
      }
      return;
    }
  }
}

void 
Histogram2::setYRange(TH1D* h, TH1D* ref, bool logScale)
{ 
  // Find maximum and minimum points of Plot to rescale Y axis
  Double_t YMax = h->GetBinContent(h->GetMaximumBin());
  Double_t YMin = 1e99;
  for (int i=1; i<=h->GetNbinsX(); i++) if (h->GetBinContent(i)>0) YMin = min(YMin, h->GetBinContent(i));
  if (ref!=NULL && ref->Integral() > h->Integral()) {
    YMax = ref->GetBinContent(ref->GetMaximumBin());
    YMin = 1e99;
    for (int i=1; i<=ref->GetNbinsX(); i++) if (ref->GetBinContent(i)>0) YMin = min(YMin, ref->GetBinContent(i));
  }
  // Set Range
  Double_t Yup(0.),Ydown(0.);
  Double_t fup(0.4),fdown(0.0);
  if (logScale) {
    YMin = max(YMin, 0.0001);
    YMax = max(YMax, 0.0001);
    Ydown = YMin/(std::pow((YMax/YMin), (fdown/(1.0-fdown-fup))));
    Yup = YMax*std::pow((YMax/YMin), (fup/(1.0-fdown-fup)));
    Ydown = max(Ydown, 0.001);
    Yup = max(Yup, 0.001);
  }
  else {
    Ydown = 0.0;
    Yup   = YMax/(1.0-fup);
  }
  h->GetYaxis()->SetRangeUser(Ydown,Yup);
}

void 
Histogram2::setYRange(THStack* h, TH1D* ref, bool logScale)
{ 
  // Find maximum and minimum points of Plot to rescale Y axis
  Double_t YMax = h->GetMaximum();
  Double_t YMin = h->GetMinimum();
  if (ref!=NULL && ref->GetBinContent(ref->GetMaximumBin()) > YMax) {
    YMax = ref->GetBinContent(ref->GetMaximumBin());
    YMin = 1e99;
    for (int i=1; i<=ref->GetNbinsX(); i++) if (ref->GetBinContent(i)>0) YMin = min(YMin, ref->GetBinContent(i));
  }
  // Set Range
  Double_t Yup(0.),Ydown(0.);
  Double_t fup(0.4),fdown(0.0);
  if (logScale) {
    YMin = max(YMin, 0.0001);
    YMax = max(YMax, 0.0001);
    Ydown = YMin/(std::pow((YMax/YMin), (fdown/(1.0-fdown-fup))));
    Yup = YMax*std::pow((YMax/YMin), (fup/(1.0-fdown-fup)));
    Ydown = max(Ydown, 0.001);
    Yup = max(Yup, 0.001);
  }
  else {
    Ydown = 0.0;
    Yup   = YMax/(1.0-fup);
  }
  h->SetMinimum(Ydown);
  h->SetMaximum(Yup);
  h->GetYaxis()->SetRangeUser(Ydown,Yup);
}

std::string  
Histogram2::getBeam(const std::string& sample)
{ 
  std::string coll = "";
  if (sample.find("pPb")!=std::string::npos)  coll = "pPb";
  if (sample.find("Pbp")!=std::string::npos)  coll = "Pbp";
  if (sample.find("PA")!=std::string::npos)   coll = "PA";
  if (sample.find("PP")!=std::string::npos)   coll = "PP";
  if (sample.find("PbPb")!=std::string::npos) coll = "PbPb";
  return coll;
}

std::string  
Histogram2::getDecayChannel(const std::string& sample)
{ 
  std::string dec = "";
  if (sample.find("WToMuNu")!=std::string::npos)  dec = "WToMuNu_";
  if (sample.find("WToTauNu")!=std::string::npos) dec = "WToTauNu_";
  if (sample.find("QCDToMu")!=std::string::npos)  dec = "QCDToMu_";
  if (sample.find("ZToMuMu")!=std::string::npos)  dec = "ZToMuMu_";
  return dec;
}

void 
Histogram2::Save(const std::string& label)
{ 
  gSystem->mkdir("Histos", kTRUE);
  if (label=="NORM_RATIO") {
    gSystem->mkdir("Histos/NormRatios", kTRUE);
    for (auto& s : TH1D_RATIO_NORM_) { for (auto& t : s.second) { for (auto& v : t.second) { if (v.second && v.second->GetSumOfWeights()>0.) {
            TFile *file = new TFile(Form("Histos/NormRatios/%s.root", v.second->GetName()),"RECREATE");
            v.second->Write();
            file->Write();
            file->Close();
            delete file;
          }
        }
      }
    }
  }
}

void 
Histogram2::Delete(void)
{
  for (auto& s : TH1D_)       { for (auto& elem : s.second) { for (auto& hist : elem.second) { if (hist.second) delete hist.second; } } }
  for (auto& s : TF1_)        { for (auto& elem : s.second) { for (auto& func : elem.second) { if (func.second) delete func.second; } } }
  for (auto& s : THStack_)    { for (auto& elem : s.second) { for (auto& stac : elem.second) { if (stac.second) delete stac.second; } } }
  for (auto& s : TH1D_RATIO_) { for (auto& elem : s.second) { for (auto& hist : elem.second) { if (hist.second) delete hist.second; } } }
  for (auto& s : TH1D_RATIO_NORM_) { for (auto& elem : s.second) { for (auto& hist : elem.second) { if (hist.second) delete hist.second; } } }
  for (auto& s : TF1_RATIO_)  { for (auto& elem : s.second) { for (auto& func : elem.second) { if (func.second) delete func.second; } } }
}

#endif

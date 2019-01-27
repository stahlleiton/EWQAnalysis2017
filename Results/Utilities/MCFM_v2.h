#include "TH1F.h"
#include "TGraphAsymmErrors.h"
#include "TROOT.h"

#include <vector>
#include <string>
#include <map>
#include <iostream>

#define NCT10 53
#define NEPS09 30
#define NCT14 57
#define NEPPS16 97
#define NnCTEQ15 33

typedef std::map< std::string , std::map< std::string , std::map< std::string , std::map< std::string , std::vector< TGraphAsymmErrors > > > > > GraphTheory;


void extractTheory(GraphTheory& graphObs)
{
  //
  const std::string DIR = "/home/llr/cms/stahl/ElectroWeakAnalysis/NOMINAL/BOSONPTCORR/EWQAnalysis2017/Theory/npdf_theory";
  //
  // Extract the observables
  //
  std::map< std::string , int > modelVec = { { "CT14" , NCT14 } , { "nCTEQ15" , (NnCTEQ15 + NCT14 - 1) } , { "EPPS16" , NEPPS16 } };
  std::map< std::string , std::vector< std::string > > varMap = { { "Cross_Section" , { "Pl" , "Mi" } } , { "ForwardBackward_Ratio" , { "Pl" , "Mi" , "Inc" } } , { "Charge_Asymmetry" , { "Inc" } } };
  //
  std::map< std::string , std::map< std::string , std::map< std::string , std::vector< TH1F > > > > PDFSetObs;
  //
  for (const auto& mod : modelVec) {
    for (const auto& var : varMap) {
      for (const auto& chg : var.second) {
        const std::string fileName = Form("%s/BACKUP/%s_%s_%s.root", DIR.c_str(), mod.first.c_str(), var.first.c_str(), chg.c_str());
        // Open the input file
        auto f = std::unique_ptr<TFile>( TFile::Open(fileName.c_str()) );
        if (f!=NULL && f->IsOpen() && !f->IsZombie()) {
          for (int iSet = 0; iSet < mod.second; iSet++) {
            const std::string hName = Form("%s_%s_%s_%d", mod.first.c_str(), var.first.c_str(), chg.c_str(), iSet);
            auto hSet = (TH1F*)f->Get(hName.c_str());
            if (hSet==NULL) { std::cout << "[ERROR] " << hName << " histogram was not found in " << fileName << std::endl; }
            else {
              PDFSetObs[mod.first][var.first][chg].push_back( *hSet );
            }
          }
        }
        else { std::cout << "[ERROR] File " << fileName << " was not found!" << std::endl; }
        if (f!=NULL) { f->Close(); }
      }
    }
  }
  //
  // Make the graphs
  //
  for (const auto& mod : PDFSetObs) {
    for (const auto& var : mod.second) {
      for (const auto& chg : var.second) {
        const std::string ch = (chg.first!="Inc" ? chg.first : "");
        auto& gr = graphObs[ch][var.first]["Theory"][mod.first];
        for (uint iSet = 0; iSet < chg.second.size(); iSet++) {
          gr.push_back(TGraphAsymmErrors());
          // Initialize the TGraph
          const uint nBin = chg.second[iSet].GetNbinsX();
          gr[iSet].Set(nBin);
          gr[iSet].SetName(Form("%s_pdf", chg.second[iSet].GetName()));
          //
          if (mod.first=="nCTEQ15" && (iSet!=24 && iSet!=13)) continue;
          // Loop on the bins
          for (uint i = 0; i < nBin; i++) {
            // Get X value and error
            const double xVal   = chg.second[iSet].GetBinCenter(i+1);
            const double xErrLo = (chg.second[iSet].GetBinWidth(i+1)/2.);
            const double xErrHi = xErrLo;
            // Extract the PDF set value and error
            const double yVal = chg.second[iSet].GetBinContent(i+1);
            double yErrLo = 0., yErrHi = 0.;
            // Set the point
            gr[iSet].SetPoint(i, xVal, yVal);
            // Set the error
            gr[iSet].SetPointError(i, xErrLo, xErrHi, yErrLo, yErrHi);
          }
        }
      }
    }
  }
};


#include "../../Results/Utilities/resultsUtils.h"
#include "modelUtils_5TeV.h"
#include "TH1F.h"
#include "TMatrixDSym.h"
#include "TVectorD.h"
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


typedef std::map< std::string , std::map< std::string , std::map< std::string , TMatrixDSym               > > > MatrixMap;
typedef std::map< std::string , std::map< std::string , std::map< std::string , std::pair< double , int > > > > Chi2Map;
typedef std::map< std::string , std::map< std::string , std::map< std::string , std::vector< std::pair< double , int > > > > > Chi2VecMap;


void makeCovarianceMatrix_5TeV()
{
  //
  // Initialize the W Boson MCFM PDF Set names
  //
  std::map< std::string , std::map< std::string , std::vector< std::string > > > PDFSetName;
  for (uint i = 0; i < NEPPS16; i++) {
    const int i1 = (i<=40) ? 0 : i-40;
    PDFSetName["EPPS16"]["Pl"].push_back(Form("MCFM_5TeV/EPPS16nlo_5TeV/W_only_nlo_CT14nlo_80___80___W1_nlo_EPPS16nlo_%d_%d.root",i1,i));
    PDFSetName["EPPS16"]["Mi"].push_back(Form("MCFM_5TeV/EPPS16nlo_5TeV/W_only_nlo_CT14nlo_80___80___W6_nlo_EPPS16nlo_%d_%d.root",i1,i));
  }
  //
  // Extract the PDF Yields
  //
  std::map< std::string , std::map< std::string , std::vector< TH1F > > > PDFSetYield;
  for (const auto& mod : PDFSetName) {
    for (const auto& chg : mod.second) {
      PDFSetYield[mod.first][chg.first].resize(chg.second.size());
      for (uint iSet = 0; iSet < chg.second.size(); iSet++) {
        auto& h = PDFSetYield.at(mod.first).at(chg.first)[iSet];
        if (!extractYield(h, chg.second[iSet])) { return; }
      }
    }
  }
  //
  // Compute the observables
  //
  std::map< std::string , std::map< std::string , std::map< std::string , std::vector< TH1F > > > > PDFSetObs;
  for (const auto& mod : PDFSetYield) {
    const uint nSet = mod.second.at("Pl").size();
    const bool isEPS09 = (mod.first.find("EPS09")!=std::string::npos);
    // Initialize the Observables
    PDFSetObs[mod.first]["Cross_Section"]["Pl"].resize(nSet);
    PDFSetObs[mod.first]["Cross_Section"]["Mi"].resize(nSet);
    PDFSetObs[mod.first]["ForwardBackward_Ratio"]["Pl" ].resize(nSet);
    PDFSetObs[mod.first]["ForwardBackward_Ratio"]["Mi" ].resize(nSet);
    PDFSetObs[mod.first]["ForwardBackward_Ratio"]["Inc"].resize(nSet);
    PDFSetObs[mod.first]["Charge_Asymmetry"]["Inc"].resize(nSet);
    // Compute the Observables
    for (uint iSet = 0; iSet < nSet; iSet++) {
      // Extract the yields
      auto& yield_Pl = PDFSetYield.at(mod.first).at("Pl")[iSet];
      auto& yield_Mi = PDFSetYield.at(mod.first).at("Mi")[iSet];
      // Get the Cross Section
      getCrossSection( PDFSetObs.at(mod.first).at("Cross_Section").at("Pl")[iSet] , yield_Pl , isEPS09);
      getCrossSection( PDFSetObs.at(mod.first).at("Cross_Section").at("Mi")[iSet] , yield_Mi , isEPS09 );
      // Get the Forward Backward Ratios
      getForwardBackwardRatio( PDFSetObs.at(mod.first).at("ForwardBackward_Ratio").at("Pl" )[iSet] , yield_Pl , isEPS09 );
      getForwardBackwardRatio( PDFSetObs.at(mod.first).at("ForwardBackward_Ratio").at("Mi" )[iSet] , yield_Mi , isEPS09 );
      getForwardBackwardRatio( PDFSetObs.at(mod.first).at("ForwardBackward_Ratio").at("Inc")[iSet] , yield_Pl , yield_Mi , isEPS09 );
      // Get the Charge Asymmetry
      getChargeAsymmetry( PDFSetObs.at(mod.first).at("Charge_Asymmetry").at("Inc")[iSet] , yield_Pl , yield_Mi , isEPS09 );
    }
  }
  //
  const std::string CWD = getcwd(NULL, 0);
  makeDir(CWD+"/BACKUP_5TeV");
  //
  // Save the observables
  //
  for (const auto& mod : PDFSetObs) {
    for (const auto& var : mod.second) {
      for (const auto& chg : var.second) {
        const std::string fileName = Form("%s/BACKUP_5TeV/%s_%s_%s.root", CWD.c_str(), mod.first.c_str(), var.first.c_str(), chg.first.c_str());
        TFile file(fileName.c_str(), "RECREATE"); file.cd();
        for (uint iSet=0; iSet<chg.second.size(); iSet++) {
          const auto& h = chg.second[iSet];
          const std::string hName = Form("%s_%s_%s_%d", mod.first.c_str(), var.first.c_str(), chg.first.c_str(), iSet);
          h.Write(hName.c_str());
        }
        file.Write();
        file.Close();
      }
    }
  }  
  //
  // Compute the PDF uncertainties
  //
  std::map< std::string , std::map< std::string , std::map< std::string , TGraphAsymmErrors > > > GraphObs;
  for (const auto& mod : PDFSetObs) {
    for (const auto& var : mod.second) {
      for (const auto& chg : var.second) {
        auto& gr = GraphObs[mod.first][var.first][chg.first];
        if (!getPDFHessianError(gr, chg.second)) { return; }
      }
    }
  }
  //
  // Extract Graphs (draw)
  //  
  for (auto& mod : GraphObs) {
    for (auto& var : mod.second) {
      for (auto& chg : var.second) {
        TCanvas c("c", "", 600, 600); c.cd();
        chg.second.Draw();
        const std::string cname = Form("%s_%s_%s", mod.first.c_str(), var.first.c_str(), chg.first.c_str());
        c.SaveAs(Form("Plot_5TeV/%s.C", cname.c_str()));
      }
    }
  }
};

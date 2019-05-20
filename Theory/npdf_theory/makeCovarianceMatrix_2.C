#include "../../Results/Utilities/resultsUtils.h"
#include "modelUtils.h"
#include "TH1F.h"
#include "TMatrixDSym.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TROOT.h"

#include "../../Results/Utilities/MCFM.h"

#include <vector>
#include <string>
#include <map>
#include <iostream>

#define NCT10 53
#define NEPS09 30
#define NCT14 57
#define NEPPS16 97
#define NnCTEQ15 33
const bool addCT14VAR = true;
const bool modelUncorr = false;
const bool dataUncorr = false;
const bool excLumiUnc = false;


typedef std::map< std::string , std::map< std::string , std::map< std::string , TMatrixDSym               > > > MatrixMap;
typedef std::map< std::string , std::map< std::string , std::map< std::string , std::pair< double , int > > > > Chi2Map;
typedef std::map< std::string , std::map< std::string , std::map< std::string , std::vector< std::pair< double , int > > > > > Chi2VecMap;
typedef std::map< std::string , std::map< std::string , std::map< std::string , TGraph > > > Chi2GraphMap;


void drawMatrix ( const MatrixMap& , const std::string& , const std::string& );
void printChi2  ( const Chi2Map&   , const std::vector<std::string> , const std::string& );
void drawGraph  ( const Chi2GraphMap& , const std::string& );
void drawChi2Hist ( const Chi2VecMap& , const std::string& );


void makeCovarianceMatrix()
{
  //
  // Initialize the W Boson MCFM PDF Set names
  //
  std::map< std::string , std::map< std::string , std::vector< std::string > > > PDFSetName;
  for (uint i = 0; i < NCT14; i++) {
    PDFSetName["CT14"]["Pl"].push_back(Form("MCFM/CT14nlo/W_only_nlo_CT14nlo_80___80___W1_nlo_CT14nlo_%d.root", i));
    PDFSetName["CT14"]["Mi"].push_back(Form("MCFM/CT14nlo/W_only_nlo_CT14nlo_80___80___W6_nlo_CT14nlo_%d.root", i));
  }
  for (uint i = 0; i < NCT10; i++) {
    PDFSetName["CT10"]["Pl"].push_back(Form("MCFM/CT10nlo/W_only_nlo_CT10nlo_80___80___W1_nlo_CT10nlo_%d_1.root",i));
    PDFSetName["CT10"]["Mi"].push_back(Form("MCFM/CT10nlo/W_only_nlo_CT10nlo_80___80___W6_nlo_CT10nlo_%d_1.root",i));
  }
  PDFSetName["nCTEQ15"]["Pl"].push_back(Form("MCFM/nCTEQ15nlo/W_only_nlo_CT14nlo_80___80___W1_nlo_nCTEQ15nlo_%d_0.root",0));
  PDFSetName["nCTEQ15"]["Mi"].push_back(Form("MCFM/nCTEQ15nlo/W_only_nlo_CT14nlo_80___80___W6_nlo_nCTEQ15nlo_%d_0.root",0));
  for (uint i = 1; i < NnCTEQ15; i++) {
    PDFSetName["nCTEQ15"]["Pl"].push_back(Form("MCFM/nCTEQ15nlo/W_only_nlo_CT14nlo_80___80___W1_nlo_nCTEQ15nlo_0_%d.root",i));
    PDFSetName["nCTEQ15"]["Mi"].push_back(Form("MCFM/nCTEQ15nlo/W_only_nlo_CT14nlo_80___80___W6_nlo_nCTEQ15nlo_0_%d.root",i));
  }
  if (addCT14VAR) {
    for (uint i = 1; i < NCT14; i++) {
      PDFSetName["nCTEQ15"]["Pl"].push_back(Form("MCFM/nCTEQ15nlo/W_only_nlo_CT14nlo_80___80___W1_nlo_nCTEQ15nlo_%d_0.root",i));
      PDFSetName["nCTEQ15"]["Mi"].push_back(Form("MCFM/nCTEQ15nlo/W_only_nlo_CT14nlo_80___80___W6_nlo_nCTEQ15nlo_%d_0.root",i));
    }
  }
  PDFSetName["EPPS16"]["Pl"].push_back(Form("MCFM/EPPS16nlo/W_only_nlo_CT14nlo_80___80___W1_nlo_EPPS16nlo_%d_%d.root",0,0));
  PDFSetName["EPPS16"]["Mi"].push_back(Form("MCFM/EPPS16nlo/W_only_nlo_CT14nlo_80___80___W6_nlo_EPPS16nlo_%d_%d.root",0,0));
  for (uint i = 1; i < (NEPPS16-NCT14+1); i++) {
    const int i1 = (i<=40) ? 0 : i-40;
    PDFSetName["EPPS16"]["Pl"].push_back(Form("MCFM/EPPS16nlo/W_only_nlo_CT14nlo_80___80___W1_nlo_EPPS16nlo_%d_%d.root",i1,i));
    PDFSetName["EPPS16"]["Mi"].push_back(Form("MCFM/EPPS16nlo/W_only_nlo_CT14nlo_80___80___W6_nlo_EPPS16nlo_%d_%d.root",i1,i));
  }
  if (addCT14VAR) {
    for (uint i = (NEPPS16-NCT14+1); i < NEPPS16; i++) {
      const int i1 = (i<=40) ? 0 : i-40;
      PDFSetName["EPPS16"]["Pl"].push_back(Form("MCFM/EPPS16nlo/W_only_nlo_CT14nlo_80___80___W1_nlo_EPPS16nlo_%d_%d.root",i1,i));
      PDFSetName["EPPS16"]["Mi"].push_back(Form("MCFM/EPPS16nlo/W_only_nlo_CT14nlo_80___80___W6_nlo_EPPS16nlo_%d_%d.root",i1,i));
    }
  }
  for (int i = 0; i < NCT10; i++) {
    PDFSetName["EPS09"]["Pl"].push_back(Form("MCFM/EPS09nlo/W_only_nlo_CT10nlo_80___80___W1_nlo_EPS09_%d_%d.root",i,1));
    PDFSetName["EPS09"]["Mi"].push_back(Form("MCFM/EPS09nlo/W_only_nlo_CT10nlo_80___80___W6_nlo_EPS09_%d_%d.root",i,1));
  }
  for (int i = 2; i < NEPS09+2; i++) {
    PDFSetName["EPS09"]["Pl"].push_back(Form("MCFM/EPS09nlo/W_only_nlo_CT10nlo_80___80___W1_nlo_EPS09_%d_%d.root",0,i));
    PDFSetName["EPS09"]["Mi"].push_back(Form("MCFM/EPS09nlo/W_only_nlo_CT10nlo_80___80___W6_nlo_EPS09_%d_%d.root",0,i));
  }
  for (int i = 0; i < NCT14; i++) {
    PDFSetName["EPS09_CT14"]["Pl"].push_back(Form("MCFM/EPS09nlo_CT14nlo/W_only_nlo_CT14nlo_80___80___W1_nlo_EPS09_%d_%d.root",i,1));
    PDFSetName["EPS09_CT14"]["Mi"].push_back(Form("MCFM/EPS09nlo_CT14nlo/W_only_nlo_CT14nlo_80___80___W6_nlo_EPS09_%d_%d.root",i,1));
  }
  for (int i = 2; i < NEPS09+2; i++) {
    PDFSetName["EPS09_CT14"]["Pl"].push_back(Form("MCFM/EPS09nlo_CT14nlo/W_only_nlo_CT14nlo_80___80___W1_nlo_EPS09_%d_%d.root",0,i));
    PDFSetName["EPS09_CT14"]["Mi"].push_back(Form("MCFM/EPS09nlo_CT14nlo/W_only_nlo_CT14nlo_80___80___W6_nlo_EPS09_%d_%d.root",0,i));
  }
  //
  std::map< std::string , std::string > LHAPDFModel = { {"CT14", "CT14nlo"}, {"CT10", "CT10nlo"}, {"nCTEQ15", "nCTEQ15_208_82"}, {"EPPS16", "EPPS16nlo_CT14nlo_Pb208"}, {"EPS09", "NONLHA"}, {"EPS09_CT14", "NONLHA"} };
  if (addCT14VAR==true ) { LHAPDFModel.at("nCTEQ15") = "NONLHA"; }
  if (addCT14VAR==false) { LHAPDFModel.at("EPPS16")  = "NONLHA"; }
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
      getCrossSection( PDFSetObs.at(mod.first).at("Cross_Section").at("Pl")[iSet] , yield_Pl );
      getCrossSection( PDFSetObs.at(mod.first).at("Cross_Section").at("Mi")[iSet] , yield_Mi );
      // Get the Forward Backward Ratios
      getForwardBackwardRatio( PDFSetObs.at(mod.first).at("ForwardBackward_Ratio").at("Pl" )[iSet] , yield_Pl );
      getForwardBackwardRatio( PDFSetObs.at(mod.first).at("ForwardBackward_Ratio").at("Mi" )[iSet] , yield_Mi );
      getForwardBackwardRatio( PDFSetObs.at(mod.first).at("ForwardBackward_Ratio").at("Inc")[iSet] , yield_Pl , yield_Mi );
      // Get the Charge Asymmetry
      getChargeAsymmetry( PDFSetObs.at(mod.first).at("Charge_Asymmetry").at("Inc")[iSet] , yield_Pl , yield_Mi );
    }
  }
  //
  const std::string CWD = getcwd(NULL, 0);
  makeDir(CWD+"/BACKUP");
  //
  // Save the observables
  //
  for (const auto& mod : PDFSetObs) {
    for (const auto& var : mod.second) {
      for (const auto& chg : var.second) {
        const std::string fileName = Form("%s/BACKUP/%s_%s_%s.root", CWD.c_str(), mod.first.c_str(), var.first.c_str(), chg.first.c_str());
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
        if (!getPDFHessianError(gr, chg.second, LHAPDFModel.at(mod.first))) { return; }
      }
    }
  }
  //
  // Check the PDF uncertainties
  //
  if (addCT14VAR) {
    std::map< std::string , std::map< std::string , std::map< std::string , TGraphAsymmErrors > > > graphMap;
    A1m(graphMap["Mi"]["ForwardBackward_Ratio"]);
    A1p(graphMap["Pl"]["ForwardBackward_Ratio"]);
    A3 (graphMap["Inc"]["ForwardBackward_Ratio"]);
    Wp (graphMap["Pl"]["Cross_Section"]);
    Wm (graphMap["Mi"]["Cross_Section"]);
    chasym(graphMap["Inc"]["Charge_Asymmetry"]);
    for (const auto& mod : GraphObs) {
      for (const auto& var : mod.second) {
        for (const auto& chg : var.second) {
	  const auto& graphI = chg.second;
	  const auto& graphT = graphMap.at(chg.first).at(var.first).at(mod.first);
	  if (graphI.GetN()!=graphT.GetN()) { std::cout << "[ERROR] The number of entries on graphI and graphT are different" << std::endl; return; }
	  for (int i=0; i<graphI.GetN(); i++) {
	    if (std::abs(graphI.GetErrorXlow(i)-graphT.GetErrorXlow(i))>0.00001) { std::cout << "[ERROR] Different low X error: " << graphI.GetErrorXlow(i) << " , " << graphT.GetErrorXlow(i) << std::endl; return; }
            if (std::abs(graphI.GetErrorXhigh(i)-graphT.GetErrorXhigh(i))>0.00001) { std::cout << "[ERROR] Different high X error: " << graphI.GetErrorXhigh(i) << " , " << graphT.GetErrorXhigh(i) << std::endl; return; }
            if (std::abs(graphI.GetErrorYlow(i)-graphT.GetErrorYlow(i))>0.00001) { std::cout << "[ERROR] Different low Y error: " << graphI.GetErrorYlow(i) << " , " << graphT.GetErrorYlow(i) << std::endl; return; }
            if (std::abs(graphI.GetErrorYhigh(i)-graphT.GetErrorYhigh(i))>0.00001) { std::cout << "[ERROR] Different high Y error: " << graphI.GetErrorYhigh(i) << " , " << graphT.GetErrorYhigh(i) << std::endl; return; }
	    double X1, Y1, X2, Y2;
	    graphI.GetPoint(i, X1, Y1);
	    graphT.GetPoint(i, X2, Y2);
	    if (std::abs(X1-X2)>0.0001) { std::cout << "[ERROR] Different X value: " << X1 << " , " << X2 << std::endl; return; }
	    if (std::abs(Y1-Y2)>0.0001) { std::cout << "[ERROR] Different Y value: " << Y1 << " , " << Y2 << std::endl; return; }
  	  }
        }
      }
    }
  }
  //
  // Compute Covariance Matrix for Theory
  //
  MatrixMap CovMatrixObs;
  std::map< std::string , std::map< std::string , std::map< std::string , TVectorD > > > NomVecObs;
  std::map< std::string , std::map< std::string , std::map< std::string , std::vector< TVectorD > > > > PDFNomVecObs;
  for (const auto& mod : GraphObs) {
    for (const auto& var : mod.second) {
      for (const auto& chg : var.second) {
        if (chg.first=="Mi") continue;
        //
        std::vector< std::string > chgVec = { "Mi" , "Pl" };
        if (chg.first=="Inc") { chgVec.clear(); chgVec.push_back("Inc"); }
        const int nBin = GraphObs.at(mod.first).at(var.first).at(chg.first).GetN();
        const int nChg = chgVec.size();
        //
        auto& PDFSet    = PDFSetObs.at(mod.first).at(var.first).at(chg.first);
        auto& PDFNomVec = PDFNomVecObs[mod.first][var.first][chg.first];
        auto& nomVec    = NomVecObs[mod.first][var.first][chg.first];
        auto& covMatrix = CovMatrixObs[mod.first][var.first][chg.first];
        covMatrix.ResizeTo( nChg*nBin , nChg*nBin );
        nomVec.ResizeTo( nChg*nBin );
        PDFNomVec.resize(PDFSet.size());
        //
        for (int i1 = 0; i1 < nChg; i1++) {
          const auto& gra1 = GraphObs.at(mod.first).at(var.first).at(chgVec[i1]);
          const auto& pdf1 = PDFSetObs.at(mod.first).at(var.first).at(chgVec[i1]);
          for (int i2 = 0; i2 < nChg; i2++) {
            const auto& gra2 = GraphObs.at(mod.first).at(var.first).at(chgVec[i2]);
            const auto& pdf2 = PDFSetObs.at(mod.first).at(var.first).at(chgVec[i2]);
            for (int j1 = 0; j1 < nBin; j1++) {
              for (int j2 = 0; j2 < nBin; j2++) {
                const double errBin1 = 0.5*( gra1.GetErrorYlow(j1) + gra1.GetErrorYhigh(j1) );
                const double errBin2 = 0.5*( gra2.GetErrorYlow(j2) + gra2.GetErrorYhigh(j2) );
                double corr12  = 1.0;
                if (j1!=j2 || i1!=i2) { if (!getPDFHessianCorrelation(corr12, j1, j2, pdf1, pdf2, LHAPDFModel.at(mod.first))) { return; } }
		if (modelUncorr && (j1!=j2 || i1!=i2)) { corr12 = 0.0; }
                const double covBin12 ( corr12 * errBin1 * errBin2 );
                covMatrix(j1 + nBin*i1 , j2 + nBin*i2) = covBin12;
              }
            }
          }
        }
        //
        for (int i = 0; i < nChg; i++) {
          const auto& gra = GraphObs.at(mod.first).at(var.first).at(chgVec[i]);
          for (int j = 0; j < nBin; j++) {
            double x, y; gra.GetPoint(j, x, y);
            nomVec[j + nBin*i] = y;
          }
        }
        //
        for (uint s = 0; s < PDFSet.size(); s++) {
          auto& pdfNom = PDFNomVec[s];
          pdfNom.ResizeTo( nChg*nBin );
          for (int i = 0; i < nChg; i++) {
            const auto& pdf = PDFSetObs.at(mod.first).at(var.first).at(chgVec[i])[s];
            for (int j = 0; j < nBin; j++) {
              pdfNom[j + nBin*i] = pdf.GetBinContent(j+1);
            }
          }
        }
      }
    }
  }
  //
  // Extract Covariance Matrix from Data
  //
  const std::string path = "/home/llr/cms/stahl/ElectroWeakAnalysis/NOMINAL/BOSONPTCORR/EWQAnalysis2017/Results/Output/NominalCM/METPF_RAW/DATA/Matrix/Covariance/PA/Data";
  for (const auto& var : CovMatrixObs.at("EPPS16")) {
    for (const auto& chg : var.second) {
      if (chg.first=="Mi") continue;
      auto& nomVec    = NomVecObs["DATA"][var.first][chg.first];
      auto& covMatrix = CovMatrixObs["DATA"][var.first][chg.first];
      const std::string ch = (chg.first!="Inc" ? chg.first : "");
      const std::string fName = Form("%s/%s/covMatrix_WToMu%s_PA_%s_Total_Total.root", path.c_str(), var.first.c_str(), ch.c_str(), var.first.c_str());
      TFile f(fName.c_str(), "READ");
      if (!f.IsOpen() && f.IsZombie()) { std::cout << "[ERROR] File " << fName << " was not found!" << std::endl; return; }
      // Extract the covariance matrix
      auto matrix = (TMatrixDSym*) f.Get("covMatrix");
      if (matrix==NULL) { std::cout << "[ERROR] Covariance matrix was not found in " << fName << std::endl; return; }
      covMatrix.ResizeTo(*matrix);
      covMatrix = *matrix;
      // Extract the nominal value vector
      auto vec = (TVectorD*) f.Get("nomVec");
      if (vec==NULL) { std::cout << "[ERROR] Nominal value vector was not found in " << fName << std::endl; return; }
      nomVec.ResizeTo(*vec);
      nomVec = *vec;
      // Close file
      f.Close();
      // Remove luminosity uncertainty if not wanted
      if (excLumiUnc && var.first=="Cross_Section") {
      	TMatrixDSym lumiMatrix; lumiMatrix.ResizeTo(covMatrix);
	for (int i1=0; i1<covMatrix.GetNcols(); i1++) { for (int i2=0; i2<covMatrix.GetNcols(); i2++) { lumiMatrix[i1][i2] = (nomVec[i1]*0.035)*(nomVec[i2]*0.035); } }
	covMatrix -= lumiMatrix;
      }
      if (dataUncorr) { for (int i1=0; i1<covMatrix.GetNcols(); i1++) { for (int i2=0; i2<covMatrix.GetNcols(); i2++) { if (i1!=i2) covMatrix[i1][i2] = 0.0; } } }
    }
  }
  //
  for (const auto& mod : CovMatrixObs) {
    for (const auto& var : mod.second) {
      for (const auto& chg : var.second) {
        const auto& covMatrix  = CovMatrixObs.at(mod.first).at(var.first).at(chg.first);
        const auto& nomVec     = NomVecObs.at(mod.first).at(var.first).at(chg.first);
        const auto& covMatrixX = CovMatrixObs.at(mod.first).at("Cross_Section").at("Pl");
        const auto& nomVecX    = NomVecObs.at(mod.first).at("Cross_Section").at("Pl");
        TMatrixDSym covMatrixT; covMatrixT.ResizeTo(covMatrix);
        TVectorD   nomVecT; nomVecT.ResizeTo(nomVec);
        const uint nColT = covMatrixT.GetNcols();
        const uint nColX = covMatrixX.GetNcols();
        TMatrixD   jacMatrixT; jacMatrixT.ResizeTo(nColT, nColX);
        std::vector<TMatrixD>  hesMatrixT;
        if (var.first=="ForwardBackward_Ratio") {
          if (chg.first!="Inc") {
            for (uint iT=0; iT<nColT; iT++) {
              const uint iXFW = ((iT<nColT/2) ? 14 : 28) + iT;
              const uint iXBW = ((iT<nColT/2) ? 13 : 47) - iT;
              nomVecT[iT] = nomVecX[iXFW]/nomVecX[iXBW];
              TMatrixD hesM; hesM.ResizeTo(nColX, nColX);
              for (uint iX=0; iX<nColX; iX++) {
                double val = 0.0;
                if (iX==iXFW) { val = +nomVecT[iT]/nomVecX[iXFW]; }
                if (iX==iXBW) { val = -nomVecT[iT]/nomVecX[iXBW]; }
                jacMatrixT[iT][iX] = val;
                for (uint iX2=0; iX2<nColX; iX2++) {
                  double val = 0.0;
                  if (iX==iXFW && iX2==iXFW) { val = 0.0; }
                  if (iX==iXFW && iX2==iXBW) { val = -nomVecT[iT]/(nomVecX[iXFW]*nomVecX[iXBW]); }
                  if (iX==iXBW && iX2==iXFW) { val = -nomVecT[iT]/(nomVecX[iXFW]*nomVecX[iXBW]); }
                  if (iX==iXBW && iX2==iXBW) { val = 2.*nomVecT[iT]/(nomVecX[iXBW]*nomVecX[iXBW]); }
                  hesM[iX][iX2] = val;
                } 
              }
              hesMatrixT.push_back(hesM);
            }
          }
          else {
            for (uint iT=0; iT<nColT; iT++) {
              const uint iXFWMi = 14 + iT;
              const uint iXBWMi = 13 - iT;
              const uint iXFWPl = 38 + iT;
              const uint iXBWPl = 37 - iT;
              nomVecT[iT] = (nomVecX[iXFWPl]+nomVecX[iXFWMi])/(nomVecX[iXBWPl]+nomVecX[iXBWMi]);
	      TMatrixD hesM; hesM.ResizeTo(nColX, nColX);
              for (uint iX=0; iX<nColX; iX++) {
                double val = 0.0;
                if (iX==iXFWMi) { val = +nomVecT[iT]/(nomVecX[iXFWPl]+nomVecX[iXFWMi]); }
                if (iX==iXFWPl) { val = +nomVecT[iT]/(nomVecX[iXFWPl]+nomVecX[iXFWMi]); }
                if (iX==iXBWMi) { val = -nomVecT[iT]/(nomVecX[iXBWPl]+nomVecX[iXBWMi]); }
                if (iX==iXBWPl) { val = -nomVecT[iT]/(nomVecX[iXBWPl]+nomVecX[iXBWMi]); }
                jacMatrixT[iT][iX] = val;
		for (uint iX2=0; iX2<nColX; iX2++) {
                  double val = 0.0;
                  if (iX==iXFWMi && iX2==iXFWMi) { val = 0.0; }
		  if (iX==iXFWMi && iX2==iXFWPl) { val = 0.0; }
		  if (iX==iXFWPl && iX2==iXFWMi) { val = 0.0; }
                  if (iX==iXFWPl && iX2==iXFWPl) { val = 0.0; }
		  if (iX==iXFWMi && iX2==iXBWMi) { val = -nomVecT[iT]/((nomVecX[iXBWPl]+nomVecX[iXBWMi])*(nomVecX[iXFWPl]+nomVecX[iXFWMi])); }
                  if (iX==iXFWMi && iX2==iXBWPl) { val = -nomVecT[iT]/((nomVecX[iXBWPl]+nomVecX[iXBWMi])*(nomVecX[iXFWPl]+nomVecX[iXFWMi])); }
		  if (iX==iXFWPl && iX2==iXBWMi) { val = -nomVecT[iT]/((nomVecX[iXBWPl]+nomVecX[iXBWMi])*(nomVecX[iXFWPl]+nomVecX[iXFWMi])); }
                  if (iX==iXFWPl && iX2==iXBWPl) { val = -nomVecT[iT]/((nomVecX[iXBWPl]+nomVecX[iXBWMi])*(nomVecX[iXFWPl]+nomVecX[iXFWMi])); }
		  if (iX==iXBWMi && iX2==iXFWMi) { val = -nomVecT[iT]/((nomVecX[iXBWPl]+nomVecX[iXBWMi])*(nomVecX[iXFWPl]+nomVecX[iXFWMi])); }
		  if (iX==iXBWMi && iX2==iXFWPl) { val = -nomVecT[iT]/((nomVecX[iXBWPl]+nomVecX[iXBWMi])*(nomVecX[iXFWPl]+nomVecX[iXFWMi])); }
		  if (iX==iXBWPl && iX2==iXFWMi) { val = -nomVecT[iT]/((nomVecX[iXBWPl]+nomVecX[iXBWMi])*(nomVecX[iXFWPl]+nomVecX[iXFWMi])); }
                  if (iX==iXBWPl && iX2==iXFWPl) { val = -nomVecT[iT]/((nomVecX[iXBWPl]+nomVecX[iXBWMi])*(nomVecX[iXFWPl]+nomVecX[iXFWMi])); }
		  if (iX==iXBWMi && iX2==iXBWMi) { val = 2.*nomVecT[iT]/((nomVecX[iXBWPl]+nomVecX[iXBWMi])*(nomVecX[iXBWPl]+nomVecX[iXBWMi])); }
		  if (iX==iXBWMi && iX2==iXBWPl) { val = 2.*nomVecT[iT]/((nomVecX[iXBWPl]+nomVecX[iXBWMi])*(nomVecX[iXBWPl]+nomVecX[iXBWMi])); }
		  if (iX==iXBWPl && iX2==iXBWMi) { val = 2.*nomVecT[iT]/((nomVecX[iXBWPl]+nomVecX[iXBWMi])*(nomVecX[iXBWPl]+nomVecX[iXBWMi])); }
		  if (iX==iXBWPl && iX2==iXBWPl) { val = 2.*nomVecT[iT]/((nomVecX[iXBWPl]+nomVecX[iXBWMi])*(nomVecX[iXBWPl]+nomVecX[iXBWMi])); }
                  hesM[iX][iX2] = val;
                }
              }
	      hesMatrixT.push_back(hesM);
            }
          }
        }
        else if (var.first=="Charge_Asymmetry") {
          for (uint iT=0; iT<nColT; iT++) {
            const uint iXMi = 0  + iT;
            const uint iXPl = 24 + iT;
            nomVecT[iT] = (nomVecX[iXPl]-nomVecX[iXMi])/(nomVecX[iXPl]+nomVecX[iXMi]);
            for (uint iX=0; iX<nColX; iX++) {
              double val = 0.0;
              if (iX==iXPl) { val = (+nomVecT[iT]/(nomVecX[iXPl]-nomVecX[iXMi]) - nomVecT[iT]/(nomVecX[iXPl]+nomVecX[iXMi])); }
              if (iX==iXMi) { val = (-nomVecT[iT]/(nomVecX[iXPl]-nomVecX[iXMi]) - nomVecT[iT]/(nomVecX[iXPl]+nomVecX[iXMi])); }
              jacMatrixT[iT][iX] = val;
            }
          }
        }
        else continue;
        //
        // Check the nominal value
        for (uint iT=0; iT<nColT; iT++) {
          if (std::abs(nomVecT[iT]-nomVec[iT])>0.0000001) {
            std::cout << "[ERROR] Test value: " << nomVecT[iT] << " , " << nomVec[iT] << " are different for " << mod.first << "  " << var.first << "  " << chg.first << std::endl; return;
          }
        }
        // Compute the covariance matrix using error propagation
        auto tmpMatrix = covMatrixX;
        covMatrixT = tmpMatrix.Similarity(jacMatrixT);
        for (uint i1=0; i1<hesMatrixT.size(); i1++) {
          for (uint i2=0; i2<hesMatrixT.size(); i2++) {
            const auto tmpM (covMatrixX*hesMatrixT[i1]*covMatrixX*hesMatrixT[i2]);
            for (uint i3=0; i3<nColX; i3++) { covMatrixT[i1][i2] += 0.5*tmpM[i3][i3]; }
          }
        }
        // Check the covariance
        for (uint i1=0; i1<nColT; i1++) {
          for (uint i2=0; i2<nColT; i2++) {
            if (std::sqrt(std::abs(covMatrixT[i1][i2]-covMatrix[i1][i2])/(var.first=="Charge_Asymmetry"?1.0:nomVecT[i1]*nomVecT[i2]))>0.035) {
              std::cout << "[ERROR] Test value: " << (covMatrixT[i1][i2]/(var.first=="Charge_Asymmetry"?1.0:nomVecT[i1]*nomVecT[i2])) << " , " << (covMatrix[i1][i2]/(var.first=="Charge_Asymmetry"?1.0:nomVecT[i1]*nomVecT[i2])) << " are different for " << mod.first << "  " << var.first << "  " << chg.first << std::endl;
            }
            if (var.first=="ForwardBackward_Ratio") {
              //std::cout << "[INFO] Test value: " << covMatrixT[i1][i2] << " -> " << covMatrix[i1][i2] << " for " << mod.first << "  " << var.first << "  " << chg.first << std::endl;
            }
          }
        }
      }
    }
  }
  return;
  //
  // Compute the Correlation Matrix
  //
  MatrixMap CorMatrixObs;
  for (const auto& mod : CovMatrixObs) {
    for (const auto& var : mod.second) {
      for (const auto& chg : var.second) {
        auto& corMatrix = CorMatrixObs[mod.first][var.first][chg.first];
        const auto& covMatrix = CovMatrixObs.at(mod.first).at(var.first).at(chg.first);
        corMatrix.ResizeTo(covMatrix);
        uint nCol = covMatrix.GetNcols();
        for (uint i1=0; i1 < nCol; i1++) {
          for (uint i2=0; i2 < nCol; i2++) {
            corMatrix(i1 , i2) = covMatrix(i1 , i2)/(std::sqrt(covMatrix(i1 , i1)*covMatrix(i2 , i2)));
          }
        }
      }
    }
  }
  //
  // Draw the Matrices
  //
  drawMatrix(CovMatrixObs, "Covariance" , CWD);
  drawMatrix(CorMatrixObs, "Correlation", CWD);
  //
  // Compute the Global Chi2
  //
  Chi2Map Chi2Obs;
  for (const auto& mod : CovMatrixObs) {
    if (mod.first=="DATA") continue;
    for (const auto& var : mod.second) {
      for (const auto& chg : var.second) {
        auto& chi2Result = Chi2Obs[mod.first][var.first][chg.first];
        const auto& nomVec_Theo    = NomVecObs.at(mod.first).at(var.first).at(chg.first);
        const auto& nomVec_Data    = NomVecObs.at("DATA").at(var.first).at(chg.first);
        const auto& covMatrix_Theo = CovMatrixObs.at(mod.first).at(var.first).at(chg.first);
        const auto& covMatrix_Data = CovMatrixObs.at("DATA").at(var.first).at(chg.first);
        const uint nCol = covMatrix_Data.GetNcols();
        // Compute Total Covariance Matrix
        TMatrixDSym covMatrix(nCol);
        covMatrix += covMatrix_Theo;
        covMatrix += covMatrix_Data;
        // Compute the inverse
        TMatrixDSym invCovMatrix(nCol); invCovMatrix = covMatrix; invCovMatrix.Invert();
        // Compute the chi2
        double chi2 = 0.0;
        for (uint i1=0; i1 < nCol; i1++) {
          for (uint i2=0; i2 < nCol; i2++) {
            // Extract Theory points
            const double t1 = nomVec_Theo[i1];
            const double t2 = nomVec_Theo[i2];
            // Extract Data points
            const double d1 = nomVec_Data[i1];
            const double d2 = nomVec_Data[i2];
            // Sum the chi2
            chi2 += (d1-t1)*invCovMatrix[i1][i2]*(d2-t2);
	    if (modelUncorr && dataUncorr) {
	      if (covMatrix[i1][i2]!=0.0 && std::abs(invCovMatrix[i1][i2]*covMatrix[i1][i2]-1.0)>0.000001) { std::cout << mod.first << " " << var.first << " has inconsistent covariance matrix " << invCovMatrix[i1][i2] << " -> " << 1./covMatrix[i1][i2] << std::endl; return; }
	    }
          }
        }
	TMatrixDSym invCovMatrixMi(nCol/2), invCovMatrixPl(nCol/2);
	for (uint i1=0; i1 < nCol/2; i1++) { for (uint i2=0; i2 < nCol/2; i2++) { invCovMatrixMi[i1][i2] = covMatrix[i1][i2]; } }
	for (uint i1=0; i1 < nCol/2; i1++) { for (uint i2=0; i2 < nCol/2; i2++) { invCovMatrixPl[i1][i2] = covMatrix[i1+nCol/2][i2+nCol/2]; } }
	invCovMatrixMi.Invert(); invCovMatrixPl.Invert();
        double chiP2 = 0.0;
        for (uint i1=0; i1 < nCol/2; i1++) {
          for (uint i2=0; i2 < nCol/2; i2++) {
            // Extract Theory points
            const double t1 = nomVec_Theo[i1];
            const double t2 = nomVec_Theo[i2];
            // Extract Data points
            const double d1 = nomVec_Data[i1];
            const double d2 = nomVec_Data[i2];
            // Sum the chi2
            chiP2 += (d1-t1)*invCovMatrixMi[i1][i2]*(d2-t2);
          }
        }
        double chiM2 = 0.0;
        for (uint i1=0; i1 < nCol/2; i1++) {
          for (uint i2=0; i2 < nCol/2; i2++) {
            // Extract Theory points
            const double t1 = nomVec_Theo[i1+nCol/2];
            const double t2 = nomVec_Theo[i2+nCol/2];
            // Extract Data points
            const double d1 = nomVec_Data[i1+nCol/2];
            const double d2 = nomVec_Data[i2+nCol/2];
            // Sum the chi2
            chiM2 += (d1-t1)*invCovMatrixPl[i1][i2]*(d2-t2);
          }
        }
	if (var.first=="Cross_Section") std::cout << mod.first << "  " << var.first << "  " << chg.first << "  " << chi2 << " , " << chiP2 << " , " << chiM2 << std::endl;
        // Save the results
        chi2Result = std::pair<double,int>({chi2, nCol});
      }
    }
  }
return;
  //
  // Print the Global Chi2
  //
  printChi2(Chi2Obs, {"CT10", "EPS09"}, CWD);
  printChi2(Chi2Obs, {"CT14", "EPPS16", "nCTEQ15"}, CWD);
  //
  // Extract Graphs (draw)
  //  
  for (auto& mod : GraphObs) {
    for (auto& var : mod.second) {
      for (auto& chg : var.second) {
        TCanvas c("c", "", 600, 600); c.cd();
        chg.second.Draw();
        const std::string cname = Form("%s_%s_%s", mod.first.c_str(), var.first.c_str(), chg.first.c_str());
        c.SaveAs(Form("Plot/%s.C", cname.c_str()));
      }
    }
  }
  //
  // Compute the Local Chi2
  //
  Chi2VecMap Chi2ObsVec;
  for (const auto& mod : PDFNomVecObs) {
    for (const auto& var : mod.second) {
      for (const auto& chg : var.second) {
        auto& chi2Result = Chi2ObsVec[mod.first][var.first][chg.first];
        const auto& nomVec_Data    = NomVecObs.at("DATA").at(var.first).at(chg.first);
        const auto& covMatrix_Data = CovMatrixObs.at("DATA").at(var.first).at(chg.first);
        const uint nCol = covMatrix_Data.GetNcols();
        // Compute Total Covariance Matrix
        TMatrixDSym covMatrix(nCol);
        covMatrix += covMatrix_Data;
        // Compute the inverse
        auto invCovMatrix = covMatrix.Invert();
        // Compute the chi2
        for (const auto& set : chg.second) {
          const auto& nomVec_Theo = set;
          double chi2 = 0.0;
          for (uint i1=0; i1 < nCol; i1++) {
            for (uint i2=0; i2 < nCol; i2++) {
              // Extract Theory points
              const double t1 = nomVec_Theo[i1];
              const double t2 = nomVec_Theo[i2];
              // Extract Data points
              const double d1 = nomVec_Data[i1];
              const double d2 = nomVec_Data[i2];
              // Sum the chi2
              chi2 += (d1-t1)*invCovMatrix[i1][i2]*(d2-t2);
            }
          }
          // Save the results
          chi2Result.push_back( std::pair<double,int>({chi2, nCol}) );
        }          
      }
    }
  }
  //
  // Print the Local Chi2 Rejected Sets
  //
  for (const auto& mod : Chi2ObsVec) {
    for (const auto& var : mod.second) {
      for (const auto& chg : var.second) {
        std::vector< int > badPDFs;
        for (uint s = 0; s < chg.second.size(); s++) {
          const double chi2 = chg.second[s].first;
          const double ndof = chg.second[s].second;
          const double prob = TMath::Prob(chi2, ndof);
          if (TMath::ErfcInverse(prob)*std::sqrt(2.0) > 5.0) { badPDFs.push_back(s); }
        }
        std::cout << "MODEL: " << mod.first << "  " << var.first << "  " << chg.first << "  has " << badPDFs.size() << " rejected by the data of " << chg.second.size() << std::endl;
      }
    }
  }
  // DO FOR nCTEQ15 : 24 (WORST) and 13 (BEST)
  //
  // Draw the Local Chi2
  //
  drawChi2Hist(Chi2ObsVec, CWD);
  //
  Chi2GraphMap Chi2Graph;
  for (const auto& var : Chi2ObsVec.at("EPPS16")) {
    for (const auto& chg : var.second) {
      for (const auto& mod : Chi2ObsVec) {
        auto& graph = Chi2Graph[var.first][chg.first][mod.first];
        const auto& chi2 = Chi2ObsVec.at(mod.first).at(var.first).at(chg.first);
        graph.Set(4*chg.second.size());
        for (uint s = 0; s < chg.second.size(); s++) {
          const double chi2N  = ( (s < chi2.size())  ? (chi2[s].first / chi2[s].second) : -1.0 );
          double x = double(s);
          graph.SetPoint(4*s+0, x, 0.0); graph.SetPoint(4*s+1, x+0.25, 0.0); graph.SetPoint(4*s+2, x+0.50, 0.0); graph.SetPoint(4*s+3, x+0.75, 0.0);
          int iM = 0; if (mod.first=="CT14") { iM = 0; } else if (mod.first=="nCTEQ15") { iM = 1; } else if (mod.first=="EPPS16") { iM = 2; }
          graph.SetPoint(4*s+iM, x+0.25*iM, chi2N);
        }
      }
    }
  }
  drawGraph(Chi2Graph, CWD);
  //
};


void printChi2(const Chi2Map& Chi2Obs, const std::vector<std::string> modelList, const std::string& outDir)
{
  //
  // Create Output Directory
  const std::string tableDir = outDir + "/Chi2/Table";
  makeDir(tableDir);
  // Create Output Files for all the convariance matrices
  const std::string fileName = Form("chi2_%s", modelList[0].c_str());
  std::ofstream file((tableDir + "/" + fileName + ".tex").c_str());
  if (file.is_open()==false) { std::cout << "[ERROR] File " << fileName << " was not created!" << std::endl; return; }
  //
  // Initialize the latex table
  std::vector< std::string > texTable;
  const uint nCol = (modelList.size() + 1);
  //
  // Determine number of columns
  std::vector< std::string > rowChg =  { "Pl" , "Inc" , "Pl" , "Inc" };
  std::vector< std::string > rowVar =  { "Cross_Section" , "Charge_Asymmetry" , "ForwardBackward_Ratio" , "ForwardBackward_Ratio" };
  //
  texTable.push_back("\\begin{table}[h!]");
  texTable.push_back("  \\centering");
  texTable.push_back("  \\resizebox{\\textwidth}{!}{");
  //
  texTable.push_back("  \\renewcommand{\\arraystretch}{1.5}");
  std::string tmp = Form("  \\begin{tabular}{c *%dc", 3);
  for (uint i = 0; i < (nCol-2); i++) { tmp += Form(" | *%dc", 3); }
  tmp += "}";
  texTable.push_back(tmp);
  texTable.push_back("    \\hline");
  tmp = ("    ");
  tmp += "\\multirow{2}{*}{Observable} ";
  for (uint i = 0; i < nCol-1; i++) {
    if      (modelList[i]=="EPPS16" ) { tmp += " & \\multicolumn{3}{c}{CT14+EPPS16}";  }
    else if (modelList[i]=="nCTEQ15") { tmp += " & \\multicolumn{3}{c}{CT14+nCTEQ15}"; }
    else if (modelList[i]=="EPS09"  ) { tmp += " & \\multicolumn{3}{c}{CT10+EPS09}";   }
    else { tmp += Form(" & \\multicolumn{3}{c}{%s}", modelList[i].c_str()); }
    if (i==(nCol-2)) { tmp += Form(" \\\\ \\cline{2-%d}", 3*i+4); }
  }
  texTable.push_back(tmp);
  tmp = ("    ");
  for (uint i = 0; i < nCol-1; i++) {
    tmp += " & $\\chi^{2}$ & NDoF & Prob.($\\%$)";
    if (i==(nCol-2)) { tmp += " \\\\"; }
  }
  texTable.push_back(tmp);
  texTable.push_back("    \\hline");
  //
  for (uint i = 0; i < rowVar.size(); i++) {
    const auto& var = rowVar[i];
    const auto& chg = rowChg[i];
    tmp = ("    ");
    tmp += "$"+formatResultVarName(var, true, false, true, (chg=="Pl"?"Com":chg))+"$";
    for (uint j = 0; j < modelList.size(); j++) {
      const auto& mod = modelList[j];
      const auto& chi2Result = Chi2Obs.at(mod).at(var).at(chg);
      const double chi2 = chi2Result.first;
      const double ndof = chi2Result.second;
      const double prob = TMath::Prob(chi2, ndof)*100.;
      tmp += Form(" & %.0f & %.0f & %.0f", chi2, ndof, prob);
      if (j==(modelList.size()-1)) { tmp += "  \\\\"; }
    }
    texTable.push_back(tmp);
  }
  texTable.push_back("    \\hline");
  //
  texTable.push_back("  \\end{tabular}");
  //
  texTable.push_back("  }");
  texTable.push_back(Form("  \\caption{%s}",
                          Form("Bla Bla")
                          )
                     );
  texTable.push_back(Form("  \\label{tab:modelComparison_%s}", modelList[0].c_str()));
  texTable.push_back("\\end{table}");
  //
  for (const auto& row : texTable) { file << row << std::endl; }
  file << std::endl; file << std::endl;
  //
};


void drawMatrix(const MatrixMap& MatrixObs, const std::string& type, const std::string& outDir)
{
  // Set Style
  setStyle();
  //
  for (const auto& mod : MatrixObs) {
    for (const auto& var : mod.second) {
      for (const auto& chg : var.second) {
        //
        // Create Canvas
        TCanvas c("c", "c", 1000, 1000); c.cd();
        c.SetRightMargin(2.8);
        //
        auto& matrix = chg.second;
        const uint nCol = matrix.GetNcols();
        const uint dCol = (var.first!="ForwardBackward_Ratio" ? 24 : 10);
        //
        TH2D graph("matrix", "", nCol, 0, nCol, std::ceil(nCol*1.3), 0, std::ceil(nCol*1.3));
        //
        for (ushort i=1; i<=nCol; i++) {
          for (ushort j=1; j<=nCol; j++) {
            graph.SetBinContent(j, (nCol+1-i), matrix(i-1,j-1) );
            graph.SetBinError(i, j, 0.0);
          }
        }
        //
        // Draw graph
        graph.Draw("COLZ");
        //
        // Format graph
        //
        //
        // Set the Axis Titles
        std::string yLabel = formatResultVarName(var.first, true, false, false, (chg.first=="Pl"?"Com":chg.first));
        if      (type=="Correlation") { graph.SetTitle(Form("Correlation Matrix;%s;%s", "", yLabel.c_str())); }
        else if (type=="Covariance" ) { graph.SetTitle(Form("Covariance Matrix;%s;%s" , "", yLabel.c_str())); }
        // X-axis
        graph.GetXaxis()->SetLabelSize(0.0);
        // Y-axis
        graph.GetYaxis()->SetTitleOffset(0.7);
        graph.GetYaxis()->SetLabelSize(0.0);
        // Z axis
        graph.GetZaxis()->SetLabelSize(0.030);
        if (type=="Correlation") { graph.GetZaxis()->SetRangeUser(-1.0 , 1.0); }
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
        // set the CMS style
        int option = 120;
        CMS_lumi(&c, option, 33, "");
        // Update
        c.Modified(); c.Update(); // Pure paranoia
        TPaletteAxis *palette = (TPaletteAxis*)graph.GetListOfFunctions()->FindObject("palette");
        palette->SetX2NDC(0.93);
        c.Modified(); c.Update(); // Pure paranoia
        //
        // Save canvas
        //
        // Create Output Directory
        const std::string plotDir = outDir + "/Matrix/" + type + "/PA/Plot/" + var.first;
        makeDir(plotDir + "/png/");
        makeDir(plotDir + "/pdf/");
        makeDir(plotDir + "/root/");
        makeDir(plotDir + "/C/");
        //
        // Save Canvas
        const std::string name = Form("%sMatrix_WToMu%s_PA_%s_%s", (type=="Correlation"?"cor":"cov"), chg.first.c_str(), var.first.c_str(), mod.first.c_str());
        c.SaveAs(( plotDir + "/png/"  + name + ".png"  ).c_str());
        c.SaveAs(( plotDir + "/pdf/"  + name + ".pdf"  ).c_str());
        c.SaveAs(( plotDir + "/root/" + name + ".root" ).c_str());
        c.SaveAs(( plotDir + "/C/"    + name + ".C"    ).c_str());
        //
        // Clean up memory
        c.Clear(); c.Close();
      }
    }
  }
};


void drawGraph(const Chi2GraphMap& Chi2Graph, const std::string& outDir)
{
  // Set Style
  setStyle();
  //
  //
  for (const auto& var : Chi2Graph) {
    for (const auto& chg : var.second) {
      //
      // Create Canvas
      TCanvas c("c", "c", 4000, 1000); c.cd();
      //
      std::map< std::string , TGraph > graph;
      for (const auto& mod : chg.second) { graph[mod.first] = mod.second; }
      //
      // Legend
      TLegend leg(0.22, 0.64, 0.42, 0.84);
      //
      // Format graph
      //
      // Set the Axis Titles
      std::string yLabel = formatResultVarName(var.first, true, false, false, (chg.first=="Pl"?"Com":chg.first));
      graph.at("CT14").SetTitle(Form("Correlation Matrix;%s;%s", "PDF set number", "#Chi^{2}/ndf"));
      graph.at("CT14").SetFillColor(kRed);
      graph.at("nCTEQ15").SetFillColor(kBlue);
      graph.at("EPPS16").SetFillColor(kGreen+2);
      // X-axis
      graph.at("CT14").GetXaxis()->SetRangeUser(0.0, 100.);
      graph.at("CT14").GetXaxis()->CenterTitle(kTRUE);
      // Y-axis
      graph.at("CT14").GetYaxis()->SetTitleOffset(0.7);
      graph.at("CT14").GetYaxis()->SetRangeUser(0.0, 10.0);
      graph.at("CT14").GetYaxis()->CenterTitle(kTRUE);
      c.Modified(); c.Update();
      //
      // Draw graph
      graph.at("CT14").Draw("AB");
      graph.at("nCTEQ15").Draw("sameB");
      graph.at("EPPS16").Draw("sameB");
      //
      leg.AddEntry(&graph.at("CT14"), "CT14");
      leg.AddEntry(&graph.at("EPPS16"), "CT14+EPPS16");
      leg.AddEntry(&graph.at("nCTEQ15"), "CT14+nCTEQ15");
      leg.Draw("same");
      //
      // Create line
      TLine line(0, 1.0, 100, 1.0);
      line.SetLineWidth(3);
      line.Draw("same");
      //
      // set the CMS style
      int option = 118;
      CMS_lumi(&c, option, 33, "");
      c.Modified(); c.Update(); // Pure paranoia
      //
      // Save canvas
      //
      // Create Output Directory
      const std::string plotDir = outDir + "/Chi2/PA/Plot/" + var.first;
      makeDir(plotDir + "/png/");
      makeDir(plotDir + "/pdf/");
      makeDir(plotDir + "/root/");
      makeDir(plotDir + "/C/");
      //
      // Save Canvas
      const std::string name = Form("Chi2_WToMu%s_PA_%s_ALL", chg.first.c_str(), var.first.c_str());
      c.SaveAs(( plotDir + "/png/"  + name + ".png"  ).c_str());
      c.SaveAs(( plotDir + "/pdf/"  + name + ".pdf"  ).c_str());
      c.SaveAs(( plotDir + "/root/" + name + ".root" ).c_str());
      c.SaveAs(( plotDir + "/C/"    + name + ".C"    ).c_str());
      //
      // Clean up memory
      c.Clear(); c.Close();
    }
  }
};


void drawChi2Hist(const Chi2VecMap& Chi2ObsVec, const std::string& outDir)
{
  // Set Style
  setStyle();
  //
  //
  for (const auto& var : Chi2ObsVec.at("EPPS16")) {
    for (const auto& chg : var.second) {
      //
      // Create Canvas
      TCanvas c("c", "c", 1000, 1000); c.cd();
      //
      // Create the Text Info
      TLatex tex; tex.SetNDC(); tex.SetTextSize(0.035); float dy = 0;
      std::vector< std::string > textToPrint;
      std::string sampleLabel = "W#kern[0.2]{#rightarrow}#kern[0.2]{#mu}#kern[0.2]{#nu_{#mu}}";
      if (chg.first == "Pl") { sampleLabel = "W^{#pm}#kern[0.2]{#rightarrow}#kern[0.2]{#mu^{#pm}}#kern[0.2]{#nu_{#mu}}"; }
      textToPrint.push_back(sampleLabel);
      textToPrint.push_back("p^{#mu}_{T} > 25 GeV/c");
      textToPrint.push_back(formatResultVarName(var.first, true, false, false, (chg.first=="Pl"?"Com":chg.first)));
      //
      // Make histogram
      //
      std::map< std::string , TH1D > graph;
      for (const auto& mod : Chi2ObsVec) {
        graph[mod.first] = TH1D((mod.first+"_TEMP").c_str(), "", 60, 0, 12);
        const auto& chi2 = mod.second.at(var.first).at(chg.first);
        for (uint s = 0; s < chi2.size(); s++) {
          graph.at(mod.first).Fill(chi2[s].first/chi2[s].second);
        }
      }
      //
      // Legend
      double xMin = 0.58, xMax = 0.78, yMin = 0.54, yMax = 0.74;
      if (var.first=="ForwardBackward_Ratio" && chg.first=="Inc") { xMin = 0.18; xMax = 0.38; yMin = 0.49; yMax = 0.69; }
      TLegend leg(xMin, yMin, xMax, yMax);
      //
      // Format graph
      //
      for (auto& gr : graph) {
        // General
        gr.second.SetTitle(";#chi^{2}/ndf;Number of error sets");
        gr.second.SetLineColor(kBlack);
        gr.second.SetMarkerSize(0);
        gr.second.SetLineStyle(1);
        // X-axis
        gr.second.GetXaxis()->CenterTitle(kTRUE);
        gr.second.GetXaxis()->SetTitleOffset(0.70);
        gr.second.GetXaxis()->SetTitleSize(0.065);
        gr.second.GetXaxis()->SetLabelSize(0.035);
        gr.second.GetXaxis()->SetLimits(0.0 , 12.0);
        // Y-axis
        gr.second.GetYaxis()->CenterTitle(kTRUE);
        gr.second.GetYaxis()->SetTitleOffset(1.05);
        gr.second.GetYaxis()->SetTitleSize(0.065);
        gr.second.GetYaxis()->SetLabelSize(0.035);
        if (addCT14VAR) { gr.second.GetYaxis()->SetRangeUser(0.0, 110.); }
        else {
          if (var.first=="Cross_Section") { gr.second.GetYaxis()->SetRangeUser(0.0, 45.); }
          if (var.first=="ForwardBackward_Ratio" && chg.first=="Pl") { gr.second.GetYaxis()->SetRangeUser(0.0, 35.); }
          if (var.first=="ForwardBackward_Ratio" && chg.first!="Pl") { gr.second.GetYaxis()->SetRangeUser(0.0, 40.); }
          if (var.first=="Charge_Asymmetry") { gr.second.GetYaxis()->SetRangeUser(0.0, 70.); }
        }
      }
      //
      graph.at("CT14").SetFillColorAlpha(kRed, 0.2);
      graph.at("nCTEQ15").SetFillColorAlpha(kBlue, 0.2);
      graph.at("EPPS16").SetFillColorAlpha(kGreen+2, 0.2);
      //
      c.Modified(); c.Update();
      //
      // Draw graph
      graph.at("EPPS16").Draw("HIST");
      graph.at("nCTEQ15").Draw("HIST same");
      graph.at("CT14").Draw("HIST same");
      //
      formatLegendEntry(*leg.AddEntry(&graph.at("CT14"), "CT14", "f"));
      formatLegendEntry(*leg.AddEntry(&graph.at("EPPS16"), "(CT14+)EPPS16", "f"));
      formatLegendEntry(*leg.AddEntry(&graph.at("nCTEQ15"), "(CT14+)nCTEQ15", "f"));
      leg.Draw("same");
      //
      // Create line
      const auto& chi0_CT14    = Chi2ObsVec.at("CT14").at(var.first).at(chg.first)[0].first/Chi2ObsVec.at("CT14").at(var.first).at(chg.first)[0].second;
      const auto& chi0_EPPS16  = Chi2ObsVec.at("EPPS16").at(var.first).at(chg.first)[0].first/Chi2ObsVec.at("EPPS16").at(var.first).at(chg.first)[0].second;
      const auto& chi0_nCTEQ15 = Chi2ObsVec.at("nCTEQ15").at(var.first).at(chg.first)[0].first/Chi2ObsVec.at("nCTEQ15").at(var.first).at(chg.first)[0].second;
      double lMax = 0.;
      if (var.first=="Cross_Section") { lMax = 30.; }
      if (var.first=="Charge_Asymmetry") { lMax = 50.; }
      if (var.first=="ForwardBackward_Ratio" && chg.first=="Pl") { lMax = 25.; }
      if (var.first=="ForwardBackward_Ratio" && chg.first!="Pl") { lMax = 15.; }
      TLine line_CT14(chi0_CT14, 0., chi0_CT14, lMax);
      TLine line_EPPS16(chi0_EPPS16, 0., chi0_EPPS16, lMax);
      TLine line_nCTEQ15(chi0_nCTEQ15, 0., chi0_nCTEQ15, lMax);
      line_CT14.SetLineColor(kRed);
      line_EPPS16.SetLineColor(kGreen+2);
      line_nCTEQ15.SetLineColor(kBlue);
      line_CT14.SetLineWidth(3);
      line_EPPS16.SetLineWidth(3);
      line_nCTEQ15.SetLineWidth(3);
      line_CT14.SetLineStyle(7);
      line_EPPS16.SetLineStyle(7);
      line_nCTEQ15.SetLineStyle(7);
      line_CT14.Draw("same");
      line_EPPS16.Draw("same");
      line_nCTEQ15.Draw("same");
      //
      // Draw the text
      tex.SetTextSize(0.055); tex.DrawLatex(0.22, 0.84, textToPrint[0].c_str());
      tex.SetTextSize(0.058); tex.SetTextFont(61); tex.DrawLatex(0.78, 0.84, "CMS"); tex.SetTextFont(62);
      bool forPaper = true;
      //if (var.first=="Cross_Section") { forPaper = true; }
      if (!forPaper) { tex.SetTextSize(0.044); tex.SetTextFont(52); tex.DrawLatex(0.69, 0.79, "Preliminary"); tex.SetTextFont(62); }
      if (textToPrint.size()>1) { tex.SetTextSize(0.035); tex.DrawLatex(0.22, 0.725, textToPrint[1].c_str()); }
      if (textToPrint.size()>2) { tex.SetTextSize(0.040); tex.DrawLatex(0.22, 0.775, textToPrint[2].c_str()); }
      //
      // set the CMS style
      // Draw the text
      int option = 118;
      CMS_lumi(&c, option, 33, "", false, 0.6, false);
      c.Modified(); c.Update(); // Pure paranoia
      //
      // Save canvas
      //
      // Create Output Directory
      const std::string plotDir = outDir + "/Chi2/PA/Plot/" + var.first;
      makeDir(plotDir + "/png/");
      makeDir(plotDir + "/pdf/");
      makeDir(plotDir + "/root/");
      makeDir(plotDir + "/C/");
      //
      // Save Canvas
      const std::string name = Form("Chi2Hist_WToMu%s_PA_%s_ALL", chg.first.c_str(), var.first.c_str());
      c.SaveAs(( plotDir + "/png/"  + name + ".png"  ).c_str());
      c.SaveAs(( plotDir + "/pdf/"  + name + ".pdf"  ).c_str());
      c.SaveAs(( plotDir + "/root/" + name + ".root" ).c_str());
      c.SaveAs(( plotDir + "/C/"    + name + ".C"    ).c_str());
      //
      // Clean up memory
      c.Clear(); c.Close();
    }
  }
};

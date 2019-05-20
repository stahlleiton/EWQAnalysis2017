#include "TFile.h"
#include "TH1.h"
#include "TRandom.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"

#include "LHAPDF/PDFSet.h"

#include <string>
#include <iostream>


const double CL68      = (100.*std::erf(1./std::sqrt(2.)));
const int    nbins     = 24;
const int    nbins2    = 10;
const double binsEP[25]  = {-2.8 , -2.6 , -2.4 , -2.2 , -2.0 , -1.8 , -1.6 , -1.4 , -1.2 , -1.0 , -0.8 , -0.6 , -0.4 , -0.2 , 0.0 , 0.2 , 0.4 , 0.6 , 0.8 , 1.0 , 1.2 , 1.4 , 1.6 , 1.8 , 2.0};
const double binsEP2[11] = {0.0 , 0.2 , 0.4 , 0.6 , 0.8 , 1.0 , 1.2 , 1.4 , 1.6 , 1.8 , 2.0};
const double bins[25]  = {-2.86 , -2.6 , -2.4 , -2.2 , -1.93 , -1.8 , -1.6 , -1.4 , -1.2 , -1.0 , -0.8 , -0.6 , -0.4 , -0.2 , 0.0 , 0.2 , 0.4 , 0.6 , 0.8 , 1.0 , 1.2 , 1.4 , 1.6 , 1.8 , 1.93};
const double bins2[11] = {0.0 , 0.2 , 0.4 , 0.6 , 0.8 , 1.0 , 1.2 , 1.4 , 1.6 , 1.8 , 1.93};


const char* RandomString(int length)
{
  char* ct = new char[length+1];
  for (int i=0; i<length; i++) { ct[i] = (char) (97+26.*(gRandom->Rndm())); }
  ct[length] = '\0';
  string result(ct);
  delete[] ct;
  return result.c_str();
};


bool extractYield(TH1F& h, const std::string& fileName, const bool doReBin=true)
{
  // Open the input file
  auto f = std::unique_ptr<TFile>( TFile::Open(fileName.c_str()) );
  if (f==NULL || !f->IsOpen() || f->IsZombie()) { std::cout << "[ERROR] File " << fileName << " was not found!" << std::endl; return false; }
  //
  bool isEPS09 = false; if (fileName.find("EPS09_")!=std::string::npos) { isEPS09 = true; }
  //
  // Initialize the histogram name
  const std::string hName = Form("hReBin_%s", RandomString(8));
  // Rebin the histogram
  if (doReBin) {
    // Initialize the input histogram
    h = TH1F( hName.c_str() , "", nbins, bins );
    // Loop over the bins
    for (uint i = 1; i <= nbins; i++) {
      // Get the bin conent
      double bincontent=0;
      double binerr=0;
      double cnt=0;
      //
      std::vector<int> iHV = { 4, 5, 6, 7, 8};
      if (isEPS09 && f->Get("id8")==NULL) { iHV = std::vector<int>({3}); }  //|| std::abs(bins[i]-bins[i-1]-0.2)<0.01
      // Extract the histogram
      for (const auto& iH : iHV) {
        if (iHV.size()>1 && bins[i-1]>=double(iH-6)) continue;
        if (iHV.size()>1 && bins[i]  <=double(iH-7)) break;
        auto hSet = (TH1F*)f->Get(Form("id%d", iH));
        if (hSet==NULL) { std::cout << "[ERROR] id" << iH << " histogram was not found in " << fileName << std::endl; f->Close(); return false; }
        for (int j = 1; j <= hSet->GetNbinsX(); j++) {
          if ( (hSet->GetBinCenter(j) > bins[i-1]) && (hSet->GetBinCenter(j) < bins[i]) ) {
            bincontent += hSet->GetBinContent(j);
            binerr += std::pow( hSet->GetBinError(j) , 2 );
            cnt++;
          }
        }
      }
      if ((iHV.size()>1 && cnt<12) || (iHV.size()==1 && cnt>3)) { std::cout << "[ERROR] cnt is " << cnt << " in " << fileName << " for bin " << i << "  " << bins[i-1] << " " << bins[i] << std::endl; return false; }
      //
      binerr = std::sqrt(binerr);
      // Set the bin content
      h.SetBinContent(i , bincontent/(double(cnt)) );
      h.SetBinError(i, binerr/(double(cnt)) );
    }
    // Scale histogram
    h.Scale(208.*1e-6);
  }
  else {
    auto hSet = (TH1F*)f->Get("id3");
    if (hSet==NULL) { std::cout << "[ERROR] id3 histogram was not found in " << fileName << std::endl; f->Close(); return false; }
    h = *((TH1F*)hSet->Clone(hName.c_str()));
  }
  // Close the input file
  f->Close();
  // Return
  return true;
};


void getCrossSection(TH1F& h, const TH1F& hYl)
{
  const std::string hName = string("hCrossSection_") + string(RandomString(8));
  h = TH1F( hName.c_str() , "" , nbins , bins );
  for (int j = 1; j <= nbins; j++) {
    const double y  = hYl.GetBinContent(j);
    const double ye = hYl.GetBinError(j);
    h.SetBinContent( j , y );
    h.SetBinError( j , ye );
  }
};


void getChargeAsymmetry(TH1F& h, const TH1F& hPl, const TH1F& hMi)
{
  const std::string hName = string("hChargeAsymmetry_") + string(RandomString(8));
  h = TH1F( hName.c_str() , "" , nbins , bins );
  for (int j = 1; j <= nbins; j++) {
    const double yp  = hPl.GetBinContent(j);
    const double ype = hPl.GetBinError(j);
    const double ym  = hMi.GetBinContent(j);
    const double yme = hMi.GetBinError(j);
    double val = 0.0;
    double err = 0.0;
    if ((yp+ym) != 0.) {
      val = ((yp - ym)/(yp + ym));
      err = std::sqrt((4*ym*ym/std::pow(yp+ym,4))*ype*ype + (4*yp*yp/std::pow(yp+ym,4)*yme*yme));
    }
    h.SetBinContent( j , val );
    h.SetBinError( j , err );
  }
};


void getForwardBackwardRatio(TH1F& h, TH1F& hYl)
{
  const std::string hName = string("hForwardBackwardRatio_") + string(RandomString(8));
  h = TH1F( hName.c_str() , "" , nbins2 , bins2 );
  for (int j = 1; j <= nbins2; j++) {
    const double x = h.GetBinCenter(j);
    const int iFw = hYl.FindBin(x);
    const int iBw = hYl.FindBin(-x);
    const double yFwV = hYl.GetBinContent( iFw );
    const double yFwE = hYl.GetBinError( iFw );
    const double yBwV = hYl.GetBinContent( iBw );
    const double yBwE = hYl.GetBinError( iBw );
    double val = 0.0;
    double err = 0.0;
    if (yBwV != 0.0) {
      val = (yFwV / yBwV);
      err = std::abs(val)*std::sqrt( std::pow( (yFwE/yFwV) , 2 ) + std::pow( (yBwE/yBwV) , 2 ) );
    }
    h.SetBinContent( j , val );
    h.SetBinError( j , err );
  }
};


void getForwardBackwardRatio(TH1F& h, TH1F& hPl, TH1F& hMi)
{
  const std::string hName = string("hForwardBackwardRatio_Inc_") + string(RandomString(8));
  h = TH1F( hName.c_str() , "" , nbins2 , bins2 );
  for (int j = 1; j <= nbins2; j++) {
    const double x = h.GetBinCenter(j);
    const int iPlFw = hPl.FindBin(x);
    const int iPlBw = hPl.FindBin(-x);
    const double yPlFwV = hPl.GetBinContent( iPlFw );
    const double yPlFwE = hPl.GetBinError( iPlFw );
    const double yPlBwV = hPl.GetBinContent( iPlBw );
    const double yPlBwE = hPl.GetBinError( iPlBw );
    const int iMiFw = hMi.FindBin(x);
    const int iMiBw = hMi.FindBin(-x);
    const double yMiFwV = hMi.GetBinContent( iMiFw );
    const double yMiFwE = hMi.GetBinError( iMiFw );
    const double yMiBwV = hMi.GetBinContent( iMiBw );
    const double yMiBwE = hMi.GetBinError( iMiBw );
    //
    const double yFwV = ( yPlFwV + yMiFwV );
    const double yFwE = std::sqrt( std::pow( yPlFwE , 2 ) + std::pow( yMiFwE , 2 ) );
    const double yBwV = ( yPlBwV + yMiBwV );
    const double yBwE = std::sqrt( std::pow( yPlBwE , 2 ) + std::pow( yMiBwE , 2 ) );
    //
    double val = 0.0;
    double err = 0.0;
    if (yBwV != 0.0) {
      val = (yFwV / yBwV);
      err = std::abs(val)*std::sqrt( std::pow( (yFwE/yFwV) , 2 ) + std::pow( (yBwE/yBwV) , 2 ) );
    }
    h.SetBinContent( j , val );
    h.SetBinError( j , err );
  }
};


void getHessianError(double& errLo, double& errHi, double& errSy, const std::vector< double >& y, const double inCL=90., const double ouCL=CL68)
{
  // Extract the central value
  const auto cntVal = y[0];
  // Extract the uncertainty
  double eLo = 0.0 , eHi = 0.0 , eSy = 0.0;
  for (uint j=1; j<=(y.size()-1)/2; j++) {
    const int j1 = 2*j-1;
    const int j2 = 2*j;
    eSy += std::pow( ( y[j1] - y[j2] ) , 2 );
    eHi += std::pow( std::max( std::max( y[j1] - cntVal , y[j2] - cntVal ) , 0. ) , 2 );
    eLo += std::pow( std::max( std::max( cntVal - y[j1] , cntVal - y[j2] ) , 0. ) , 2 );
  }
  eSy = 0.5*std::sqrt( eSy );
  eHi = std::sqrt( eHi );
  eLo = std::sqrt( eLo );
  // Convert from inCL to ouCL
  const double convFactor = TMath::ErfcInverse((1.-(ouCL/100.)))/TMath::ErfcInverse((1.-(inCL/100.)));
  eSy *= convFactor;
  eHi *= convFactor;
  eLo *= convFactor;
  // Return the result
  errSy = eSy;
  errHi = eHi;
  errLo = eLo;
};


bool getPDFHessianError(TGraphAsymmErrors& gr, const std::vector< TH1F >& hVec, const std::string& pdfName="NONLHA", const double CL=90.)
{
  if (hVec.size()==0) { std::cout << "[ERROR] getPDFHessianError: The vector of histograms is empty" << std::endl; return false; }
  // Initialize the TGraph
  const uint nBin = hVec[0].GetNbinsX();
  gr.Set(nBin);
  gr.SetName(Form("%s_pdf", hVec[0].GetName()));
  // Loop on the bins
  for (uint i = 0; i < nBin; i++) {
    // Get X value and error
    const double xVal   = hVec[0].GetBinCenter(i+1);
    const double xErrLo = (hVec[0].GetBinWidth(i+1)/2.);
    const double xErrHi = xErrLo;
    // Extract the PDF set values
    std::vector< double > y;
    for (uint j=0; j < hVec.size(); j++) { y.push_back(hVec[j].GetBinContent(i+1)); }
    // Get Y value and error
    const double yVal = y[0];
    double yErrLo = -1., yErrHi = -1., yErrSy = -1.;
    if (pdfName!="NONLHA") {
      LHAPDF::PDFSet pdfSet(pdfName.c_str());
      const auto nm = pdfSet.size();
      if (hVec.size()!=nm) { std::cout << "[ERROR] Expected " << nm << " histos for " << pdfName << ", but got << " << hVec.size() << std::endl; return false; }
      LHAPDF::PDFUncertainty u = pdfSet.uncertainty(y, CL68);
      yErrLo = u.errminus;
      yErrHi = u.errplus;
      yErrSy = u.errsymm;
      double yErrLoOld = -1., yErrHiOld = -1., yErrSyOld = -1.;
      getHessianError(yErrLoOld, yErrHiOld, yErrSyOld, y, CL, CL68);
      if (pdfSet.errorConfLevel()!=CL) { std::cout << "[ERROR] Assumed CL: " << CL << " is different from LHAPDF " << pdfName << " CL: " << pdfSet.errorConfLevel() << std::endl; return false; }
      if (std::abs(yErrLo-yErrLoOld)*100./yVal>0.0001) { std::cout << "[ERROR] Old ErrMinus: " << yErrLoOld << " is different from LHAPDF " << pdfName << " ErrMinus: " << yErrLo << std::endl; return false; }
      if (std::abs(yErrHi-yErrHiOld)*100./yVal>0.0001) { std::cout << "[ERROR] Old ErrPlus: "  << yErrHiOld << " is different from LHAPDF " << pdfName << " ErrPlus: "  << yErrHi << std::endl; return false; }
      if (std::abs(yErrSy-yErrSyOld)*100./yVal>0.0001) { std::cout << "[ERROR] Old ErrSymm: "  << yErrSyOld << " is different from LHAPDF " << pdfName << " ErrSymm: "  << yErrSy << std::endl; return false; }
    }
    else { getHessianError(yErrLo, yErrHi, yErrSy, y, CL, CL68); }
    // Set the point
    gr.SetPoint(i, xVal, yVal);
    // Set the error
    gr.SetPointError(i, xErrLo, xErrHi, yErrLo, yErrHi);
  }
  return true;
};


double getHessianCorrelation(const std::vector<double>& valuesA, const std::vector<double>& valuesB)
{
  // Get the uncertainty for A
  double errLoA, errHiA, errSyA;
  getHessianError(errLoA, errHiA, errSyA, valuesA, 90., 90.);
  // Get the uncertainty for B
  double errLoB, errHiB, errSyB;
  getHessianError(errLoB, errHiB, errSyB, valuesB, 90., 90.);
  //
  const int nmem = valuesA.size()-1;
  double cor = 0.0;
  // Calculate the correlation using Eq. (2.5) of arXiv:1106.5788v2.
  for (int ieigen = 1; ieigen <= nmem/2; ieigen++) {
    const int j1 = 2*ieigen-1;
    const int j2 = 2*ieigen;
    cor += (valuesA[j1]-valuesA[j2]) * (valuesB[j1]-valuesB[j2]);
  }
  if (errSyA==0.) { std::cout << "[ERROR] The error A is " << errSyA << " Low " << errLoA << " High " << errHiA << std::endl; return 0.0; }
  if (errSyB==0.) { std::cout << "[ERROR] The error B is " << errSyB << " Low " << errLoB << " High " << errHiB << std::endl; return 0.0; }
  // Normalize the correlation
  cor /= 4.0*errSyA*errSyB;
  if (cor<-1. || cor>1.) { std::cout << "[ERROR] The correlation: " << cor << " is wrong!" << std::endl; return 0.0; }
  // Return
  return cor;
};


bool getPDFHessianCorrelation(double& corr, const int binA, const int binB, const std::vector< TH1F >& hVecA, const std::vector< TH1F >& hVecB, const std::string& pdfName="NONLHA")
{
  if (hVecA.size()==0 || hVecB.size()==0) { std::cout << "[ERROR] getPDFHessianCorrelation: The vector of histograms is empty" << std::endl; return false; }
  if (hVecA.size()!=hVecB.size()) { std::cout << "[ERROR] getPDFHessianCorrelation: The vectors of histograms are different" << std::endl; return false; }
  // Extract the PDF set values
  std::vector< double > valuesA;
  for (uint j=0; j < hVecA.size(); j++) { valuesA.push_back(hVecA[j].GetBinContent(binA+1)); }
  std::vector< double > valuesB;
  for (uint j=0; j < hVecB.size(); j++) { valuesB.push_back(hVecB[j].GetBinContent(binB+1)); }
  // Get correlation
  if (pdfName!="NONLHA") {
    LHAPDF::PDFSet pdfSet(pdfName.c_str());
    const auto nm = pdfSet.size();
    if (hVecA.size()!=nm || hVecB.size()!=nm) { std::cout << "[ERROR] Expected " << nm << " histos for " << pdfName << ", but got << " << hVecA.size() << "," << hVecB.size() << std::endl; return false; }
    corr = pdfSet.correlation(valuesA, valuesB);
    const double corrOld = getHessianCorrelation(valuesA, valuesB);
    if (corr!=corrOld) { std::cout << "[ERROR] Old correlation: " << corrOld << " is different from LHAPDF " << pdfName << " correction: " << corr << std::endl; return false; }
  }
  else { corr = getHessianCorrelation(valuesA, valuesB); }
  // Return
  return true;
};

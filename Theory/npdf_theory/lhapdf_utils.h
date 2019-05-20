#ifndef lhapdf_utils_h
#define lhapdf_utils_h

#include "TH1.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"
#include <iostream>

const double CL68 = (100*std::erf(1./std::sqrt(2.)));

void getHessianError(double& errLo, double& errHi, const std::vector< double >& y, const double inCL=90., const double ouCL=CL68)
{
  // Extract the central value
  const auto cntVal = y[0];
  // Extract the uncertainty
  double eLo = 0.0 , eHi = 0.0;
  for (uint j=1; j<=(y.size()-1)/2; j++) {
    const int j1 = 2*j-1;
    const int j2 = 2*j;
    eLo += std::pow( std::min( std::min( y[j1] - cntVal , y[j2] - cntVal ) , 0. ) , 2 );
    eHi += std::pow( std::max( std::max( y[j1] - cntVal , y[j2] - cntVal ) , 0. ) , 2 );
  }
  eLo = std::sqrt( eLo );
  eHi = std::sqrt( eHi );
  // Convert from inCL to ouCL
  const double convFactor = TMath::ErfcInverse((1.-(ouCL/100.)))/TMath::ErfcInverse((1.-(inCL/100.)));
  eLo *= convFactor;
  eHi *= convFactor;
  // Return the result
  errLo = eLo;
  errHi = eHi;
};


TGraphAsymmErrors* getPDFHessianError(const std::vector< TH1F* >& h, const double CL=90.)
{
  if (h.size()==0) { std::cout << "[ERROR] getPDFHessianError: The vector of histograms is empty" << std::endl; return NULL; }
  for (uint j=0; j<h.size(); j++) { if (h[j]==NULL ) { std::cout << "[ERROR] getPDFHessianError: The " << j << " histogram is NULL"     << std::endl; return NULL; } }
  // Initialize the TGraph
  const uint nBin = h[0]->GetNbinsX();
  auto gr = new TGraphAsymmErrors(nBin);
  gr->SetName(Form("%s_pdf", h[0]->GetName()));
  // Loop on the bins
  for (uint i = 0; i < nBin; i++) {
    // Get X value and error
    const double xVal   = h[0]->GetBinCenter(i+1);
    const double xErrLo = (h[0]->GetBinWidth(i+1)/2.);
    const double xErrHi = xErrLo;
    // Extract the PDF set values
    std::vector< double > y;
    for (uint j=0; j<h.size(); j++) { y.push_back(h[j]->GetBinContent(i+1)); }
    // Get Y value and error
    const double yVal = y[0];
    double yErrLo = -1., yErrHi = -1.;
    getHessianError(yErrLo, yErrHi, y, CL);
    // Set the point
    gr->SetPoint(i, xVal, yVal);
    // Set the error
    gr->SetPointError(i, xErrLo, xErrHi, yErrLo, yErrHi);
  }
  return gr;
};


TGraphAsymmErrors* hist2graph(const TH1* hist, const double syst=0.)
{
   const int nBin = hist->GetNbinsX();
   auto gr = new TGraphAsymmErrors(nBin);
   gr->Set(nBin);
   gr->SetName(Form("%s_graph", hist->GetName()));
   for (int i=0; i<nBin; i++) {
     const double x = hist->GetBinCenter(i+1);
     const double ex = hist->GetBinWidth(i+1)/2.;
     const double y  = hist->GetBinContent(i+1);
     const double ey = std::sqrt( std::pow( hist->GetBinError(i+1) , 2 ) + (syst*syst) );
     gr->SetPoint(i, x, y);
     gr->SetPointError(i, ex, ex, ey, ey);
   }
   return gr;
};


#endif // #ifndef lhapdf_utils_h

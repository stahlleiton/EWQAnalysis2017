#include <iostream>
#include <vector>
#include <map>
#include <iterator>
#include <algorithm>

#include <TSystem.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1D.h>
#include <TGraphAsymmErrors.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLine.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>


using namespace std;

// Binning
const double bins_eta[] = {-2.4,-2.2,-2.0,-1.8,-1.6,-1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4};
const int nbins_eta = sizeof(bins_eta)/sizeof(double) -1;
const double bins_pt[] = {0, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90};
const int nbins_pt = sizeof(bins_pt)/sizeof(double) -1;



void SetHistStyle(TH1 *h, const int color, const int marker, const int line,
    const double xmin=0, const double xmax=0, const double ymin=0, const double ymax=1) {
  h->SetLineColor(color);
  h->SetLineWidth(1);
  h->SetLineStyle(line);
  h->SetMarkerStyle(marker);
  h->SetMarkerColor(color);

  if (xmin!=0 && xmax!=0) {
    h->GetXaxis()->SetRangeUser(xmin,xmax);
  }
  if (ymin!=0 && ymax!=0) {
    h->GetYaxis()->SetRangeUser(ymin,ymax);
  }
}

void SetGraphStyle(TGraph *h, const int color, const int marker, const int line,
    const double xmin=0, const double xmax=0, const double ymin=0, const double ymax=1) {
  h->SetLineColor(color);
  h->SetLineWidth(1);
  h->SetLineStyle(line);
  h->SetMarkerStyle(marker);
  h->SetMarkerColor(color);

  if (!(xmin==0 && xmax==0)) {
    h->GetXaxis()->SetLimits(xmin,xmax);
  }
  if (!(ymin==0 && ymax==0)) {
    h->SetMinimum(ymin);
    h->SetMaximum(ymax);
  }
}


void CheckUnderFlow(TH1 &hnum, TH1 &hden) {
  for (int j=0; j<=hnum.GetNbinsX()+1; j++) {
    double num0 = hnum.GetBinContent(j);
    double den0 = hden.GetBinContent(j);
    if ((j==0 && num0>den0) || (j==hnum.GetNbinsX()+1 && num0>den0)) {
      hnum.SetBinContent(j,0);
      hden.SetBinContent(j,0);
      hnum.SetBinError(j,0);
      hden.SetBinError(j,0);
      cout << "CheckUnderFlow: " << hnum.GetName() << endl;
    }   
  }
}


unsigned short trackDownMothers(
    const std::vector <int> &v_pdgid,
    const std::vector <unsigned char> &v_status,
    const std::vector <std::vector<unsigned short>> &vv_mother,
    const int ipar, const int pdgid, const unsigned char status) {
  
  unsigned short mother = 9999; // obtain mother's idx
  if (vv_mother.size()==0) return mother;
  vector <unsigned short> v_mother = vv_mother.at(ipar);
  if (v_mother.size()==0) return mother;

  for (auto v_int=v_mother.begin(); v_int!=v_mother.end(); ++v_int) {
    mother = *v_int; // obtain mother's idx
    std::cout << "trackDownMothers() " << std::distance(v_mother.begin(),v_int) << ": " << pdgid << "/" << static_cast<unsigned short>(status) << " "
              << mother << " " << v_pdgid.at(mother) << "/" << static_cast<unsigned short>(v_status.at(mother)) << "?" << (status == v_status.at(mother)) << std::endl;
    if (pdgid == v_pdgid.at(mother) && status == v_status.at(mother)) {
      std::cout << "trackDownMothers() found!! " << *v_int << " " << pdgid << "/" << static_cast<unsigned short>(status) << std::endl;
      return mother;
    } else {
      mother = trackDownMothers(v_pdgid, v_status, vv_mother, mother, pdgid, status);
      if (mother == *v_int) { // mother's idx is idx of itself -> infinite loop, skip this particle
        return 8888;
      }
      // matched case found, return the number (idx)
      if (mother!=9999) return mother;
    }
  }

  return mother;
}



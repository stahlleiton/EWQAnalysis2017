#include "TFile.h"
#include "TH1.h"
#include "TString.h"


bool kint = true;
const int nvars = 45;
TString varnames[nvars] = {
   "PF_Muon_IsoPFR03","PF_Muon_IsoPFR03NoPUCorr","PF_Muon_IsoPFR04","PF_Muon_IsoPFR04NoPUCorr","PF_Muon_EM_Chg_sumR03Pt",
   "PF_Muon_EM_Chg_sumR04Pt","PF_Muon_EM_Neu_sumR03Pt","PF_Muon_EM_Neu_sumR04Pt","PF_Muon_Had_Chg_sumR03Pt","PF_Muon_Had_Chg_sumR04Pt",
   "PF_Muon_Had_Neu_sumR03Pt","PF_Muon_Had_Neu_sumR04Pt","PF_Muon_Had_PU_sumR03Pt","PF_Muon_Had_PU_sumR04Pt","PF_DiMuon_VtxProb",
   "PF_Muon_IsoPFR03_reliso","PF_Muon_IsoPFR03NoPUCorr_reliso","PF_Muon_IsoPFR04_reliso","PF_Muon_IsoPFR04NoPUCorr_reliso","PF_Muon_EM_Chg_sumR03Pt_reliso",
   "PF_Muon_EM_Chg_sumR04Pt_reliso","PF_Muon_EM_Neu_sumR03Pt_reliso","PF_Muon_EM_Neu_sumR04Pt_reliso","PF_Muon_Had_Chg_sumR03Pt_reliso","PF_Muon_Had_Chg_sumR04Pt_reliso",
   "PF_Muon_Had_Neu_sumR03Pt_reliso","PF_Muon_Had_Neu_sumR04Pt_reliso","PF_Muon_Had_PU_sumR03Pt_reliso","PF_Muon_Had_PU_sumR04Pt_reliso", "Reco_Muon_IsoR03_reliso",
   "Reco_Muon_IsoR05", "Reco_Muon_Trk_sumR03Pt", "Reco_Muon_Trk_sumR05Pt", "Reco_Muon_IsoR03", "Reco_Muon_IsoR05_reliso",
   "Reco_Muon_Trk_sumR03Pt_reliso", "Reco_Muon_Trk_sumR05Pt_reliso", "PF_Muon_myIsoPFR010", "PF_Muon_myIsoPFR015", "PF_Muon_myIsoPFR020",
   "PF_Muon_myIsoPFR025", "PF_Muon_myIsoPFR030", "PF_Muon_myIsoPFR035", "PF_Muon_myIsoPFR040", "PF_Muon_myIsoPFR045"
};

void cumulative(TH1F* hist);
TGraph* roc(TH1F* OS, TH1F* SS);
double area(TGraph *g);
double integralUFOF(TH1 *h);

void plot(const char* dataname="data.root", const char* mcname="mc.root") {
   TCanvas *c1 = new TCanvas();
   TFile *fdata = TFile::Open(dataname);
   TFile *fmc = TFile::Open(mcname);

   TH2F *haxes_roc = new TH2F("haxes_roc",";Sig. eff;1-Bkg. eff.",1,0,1,1,0,1);

   float maxroc_SS=-1;
   float maxroc_bkg=-1;
   TString maxroc_SS_name;
   TString maxroc_bkg_name;

   for (int i=0; i<nvars; i++) {
      TH1F *hdata_OS = (TH1F*) fdata->Get("hist_" + varnames[i] + "_OS");
      TH1F *hdata_SS = (TH1F*) fdata->Get("hist_" + varnames[i] + "_SS");
      TH1F *hdata_bkg = (TH1F*) fdata->Get("hist_" + varnames[i] + "_bkg");
      TH1F *hmc_OS = (TH1F*) fmc->Get("hist_" + varnames[i] + "_OS");
      TH1F *hmc_SS = (TH1F*) fmc->Get("hist_" + varnames[i] + "_SS");

      hdata_OS->Scale(1./integralUFOF(hdata_OS));
      hdata_SS->Scale(1./integralUFOF(hdata_SS));
      hdata_bkg->Scale(1./integralUFOF(hdata_bkg));
      hmc_OS->Scale(1./integralUFOF(hmc_OS));
      hmc_SS->Scale(1./integralUFOF(hmc_SS));

      if (kint) {
         cumulative(hdata_OS);
         cumulative(hdata_SS);
         cumulative(hdata_bkg);
         cumulative(hmc_OS);
         cumulative(hmc_SS);
      }

      hdata_OS->SetLineColor(kBlack);
      hdata_OS->SetMarkerColor(kBlack);
      hdata_OS->GetXaxis()->SetTitle(varnames[i]);
      hdata_SS->SetLineColor(kBlue);
      hdata_SS->SetMarkerColor(kBlue);
      hdata_bkg->SetLineColor(kGreen+2);
      hdata_bkg->SetMarkerColor(kGreen+2);
      hmc_OS->SetLineColor(kRed);
      hmc_OS->SetMarkerColor(kRed);

      hdata_OS->Draw();
      if (kint) hdata_OS->GetYaxis()->SetRangeUser(0,1);
      hdata_SS->Draw("same");
      hdata_bkg->Draw("same");
      hmc_OS->Draw("same");

      TLegend *tleg = new TLegend(0.5,0.2,0.9,0.5);
      tleg->SetBorderSize(0);
      tleg->AddEntry(hdata_OS,"Data OS", "p");
      tleg->AddEntry(hmc_OS,"MC OS", "p");
      tleg->AddEntry(hdata_SS,"Data SS", "p");
      tleg->AddEntry(hdata_bkg,"Data low MET", "p");
      tleg->Draw();

      c1->SaveAs(varnames[i] + ".pdf");

      if (kint) {
         TText tt;
         cout << varnames[i] << ": ";
         TGraph *thegraph = roc(hdata_OS, hdata_SS);
         haxes_roc->Draw();
         thegraph->Draw("P");
         float a = area(thegraph);
         cout << a << " ";
         tt.DrawText(0.1,0.6,Form("Area: %.3f",a));
         tt.DrawText(0.1,0.7,varnames[i]);
         c1->SaveAs(varnames[i] + "_rocSS.pdf");
         if (a>maxroc_SS) {
            maxroc_SS = a;
            maxroc_SS_name = varnames[i];
         }

         thegraph = roc(hdata_OS, hdata_bkg);
         haxes_roc->Draw();
         thegraph->Draw("P");
         a = area(thegraph);
         cout << area(thegraph) << endl;
         tt.DrawText(0.1,0.6,Form("Area: %.3f",a));
         tt.DrawText(0.1,0.7,varnames[i]);
         c1->SaveAs(varnames[i] + "_roc.pdf");
         if (a>maxroc_bkg) {
            maxroc_bkg = a;
            maxroc_bkg_name = varnames[i];
         }
      } // if kint
   } // varname loop

   cout << "Max area for SS ROC: " << maxroc_SS << " (" << maxroc_SS_name << ")" << endl;
   cout << "Max area for bkg ROC: " << maxroc_bkg << " (" << maxroc_bkg_name << ")" << endl;
}

void cumulative(TH1F* hist) {
   int nentries = hist->GetNbinsX();
   TH1F *historig = (TH1F*) hist->Clone("historig");
   for (int i=1; i<nentries+2; i++) {
      double val, err;
      val = historig->IntegralAndError(0,i,err);
      hist->SetBinContent(i,val);
      hist->SetBinError(i,err);
   }
   hist->GetYaxis()->SetTitle("Efficiency");
   delete historig;
}

TGraph* roc(TH1F* sig, TH1F* bkg) {
   int n = sig->GetNbinsX();
   double *x = new double[n+2];
   double *y = new double[n+2];

   // what is the working point that makes us closest to (1,1)?
   double fom = 0;
   double maxfom = 0;
   double cut_maxfom = 0;

   x[0]=0; y[0]=1;
   x[n+1]=1; y[n+1]=0;
   for (int i=1; i<n+1; i++) {
      x[i] = sig->GetBinContent(i);
      y[i] = 1.-bkg->GetBinContent(i);

      fom = x[i]*y[i];
      if (fom>maxfom) {
         maxfom = fom;
         cut_maxfom = sig->GetBinLowEdge(i);
      }
   }

   TGraph *ans = new TGraph(n+2,x,y);
   delete[] x; delete[] y;

   cout << "Optimal cut for " << sig->GetTitle() << " obtained for iso = " << cut_maxfom << endl;
   return ans;
}

double area(TGraph *g) {
   double area=0;
   int n = g->GetN();
   double *x = g->GetX();
   double *y = g->GetY();

   for (int i=0; i<g->GetN()-1; i++) {
      area += (y[i+1]+y[i])*(x[i+1]-x[i])/2.;
   }

   return area;
}

double integralUFOF(TH1 *h) {
   return h->Integral() + h->GetBinContent(0) + h->GetBinContent(h->GetNbinsX()+1);
}

#include <vector>
#include <map>
#include <string>
#include <iostream>

#include "TFile.h"
#include "TH1.h"
#include "TRandom.h"

using namespace std;

const int nbins=24;
const int nbins2=10;
const double bins[25] = {-2.8,-2.6,-2.4 , -2.2 , -2.0 , -1.8 , -1.6 , -1.4 , -1.2 , -1.0 , -0.8 , -0.6 , -0.4 , -0.2 , 0.0 , 
   0.2 , 0.4 , 0.6 , 0.8 , 1.0 , 1.2 , 1.4 , 1.6 , 1.8 , 2.0};
const double bins2[11] = {0.0, 0.2 , 0.4 , 0.6 , 0.8 , 1.0 , 1.2 , 1.4 , 1.6 , 1.8 , 2.0};

const char* RandomString(int length)
{
	char* ct = new char[length+1];
	for (int i=0; i<length; i++)
		ct[i] = (char) (97+26.*(gRandom->Rndm()));
	ct[length] = '\0';
	string result(ct);
   delete[] ct;
	return result.c_str();
}

TH1F* rebin(const char* filename, bool dorebin=true)
{
   string hname = string("hrebin_") + string(RandomString(8));
   TH1F *hrebin = new TH1F(hname.c_str(),hname.c_str(),nbins,bins);

   TFile *f = TFile::Open(filename);
   if (!(f->IsOpen())) {cout << "Error, couldn't open " << filename << endl; return NULL;}
   TH1F *id3 = (TH1F*) f->Get("id3");
   if (!id3) {cout << "Error, couldn't find id3 in " << filename << endl; return NULL;}

   if (!dorebin) {
      TH1F *ans = (TH1F*) id3->Clone(TString("hrebin_") + TString(RandomString(8)));
      f->Close();
      return ans;
   }

   for (int i=1; i<nbins+1; i++)
   {
      double bincontent=0;
      double binerr=0;
      double cnt=0;

      for (int j=1; j<id3->GetNbinsX()+1; j++)
      {
         // cout << i << " " << j << " " << id3->GetBinCenter(j) << " " << bins[i-1] << " " << bins[i] << endl;
         if (id3->GetBinCenter(j)>bins[i-1]&&id3->GetBinCenter(j)<bins[i])
         {
            bincontent+=id3->GetBinContent(j);
            binerr+=pow(id3->GetBinError(j),2);
            cnt++;
         }
      }
      binerr = sqrt(binerr);

      hrebin->SetBinContent(i,bincontent/((double) cnt));
      hrebin->SetBinError(i,binerr/((double) cnt));
      // cout << hrebin->GetBinCenter(i) << " " <<  hrebin->GetBinContent(i) << endl;
   }

   hrebin->Scale(208.*1e-6);

   f->Close();

   return hrebin;
}

TH1F* chasym(TH1F *hplus, TH1F *hminus)
{
   string hname = string("hasym_") + string(RandomString(8));

   TH1F *hasym = new TH1F(hname.c_str(),hname.c_str(),nbins,bins);

   for (int j=1; j<nbins+1; j++)
   {
      double yp = hplus->GetBinContent(j);
      double ype = hplus->GetBinError(j);
      double ym = hminus->GetBinContent(j);
      double yme = hminus->GetBinError(j);
      hasym->SetBinContent(j,(yp+ym!=0) ? (yp-ym)/(yp+ym) : 0);
      hasym->SetBinError(j,(yp+ym!=0) ? sqrt((4*ym*ym/pow(yp+ym,4))*ype*ype + (4*yp*yp/pow(yp+ym,4)*yme*yme)) : 0);
   }

   // hasym->Draw();

   return hasym;
}

TH1F* chasym(const char* nameplus, const char* nameminus)
{
   TH1F *hplus = rebin(nameplus);
   TH1F *hminus = rebin(nameminus);

   return chasym(hplus,hminus);
}

TH1F* a1plus(TH1F *hplus)
{
   string hname = string("ha1p_") + string(RandomString(8));

   TH1F *hasym = new TH1F(hname.c_str(),hname.c_str(),nbins2,bins2);

   for (int j=1; j<nbins2+1; j++)
   {
      // int nh = nbins/2;
      // int jm = nh+1-j;
      // int jp = nh+j;
      double x = hasym->GetBinCenter(j);
      int jm = hplus->FindBin(-x);
      int jp = hplus->FindBin(x);
      double ym = hplus->GetBinContent(jm);
      double yme = hplus->GetBinError(jm);
      double yp = hplus->GetBinContent(jp);
      double ype = hplus->GetBinError(jp);
      hasym->SetBinContent(j,(ym!=0) ? yp/ym : 0);
      hasym->SetBinError(j,(ym!=0) ? fabs(yp/ym)*sqrt(pow(ype/yp,2) + pow(yme/ym,2)) : 0);
   }

   // hasym->Draw();

   return hasym;
}

TH1F* a1plus(const char* nameplus)
{
   TH1F *hplus = rebin(nameplus);
   return a1plus(hplus);
}

TH1F* a3(TH1F *hplus, TH1F *hminus)
{
   string hname = string("ha1p_") + string(RandomString(8));

   TH1F *hasym = new TH1F(hname.c_str(),hname.c_str(),nbins2,bins2);

   for (int j=1; j<nbins2+1; j++)
   {
      // int nh = nbins/2;
      // int jm = nh+1-j;
      // int jp = nh+j;
      double x = hasym->GetBinCenter(j);
      int jm = hplus->FindBin(-x);
      int jp = hplus->FindBin(x);
      double ypm = hplus->GetBinContent(jm);
      double ypme = hplus->GetBinError(jm);
      double ypp = hplus->GetBinContent(jp);
      double yppe = hplus->GetBinError(jp);
      double ymm = hminus->GetBinContent(jm);
      double ymme = hminus->GetBinError(jm);
      double ymp = hminus->GetBinContent(jp);
      double ympe = hminus->GetBinError(jp);
      hasym->SetBinContent(j,(ypm+ymm!=0) ? (ypp+ymp)/(ypm+ymm) : 0);
      hasym->SetBinError(j,(ypm+ymm!=0) ? fabs((ypp+ymp)/(ypm+ymm))*sqrt((yppe*yppe+ympe*ympe)/pow(ypp+ymp,2) + (ypme*ypme+ymme*ymme)/pow(ypm+ymm,2)) : 0);
   }

   // hasym->Draw();

   return hasym;
}

TH1F* a3(const char* nameplus, const char* nameminus)
{
   TH1F *hplus = rebin(nameplus);
   TH1F *hminus = rebin(nameminus);

   return a3(hplus, hminus);
}

void reverse(TH1F *hist) {
   double n = hist->GetNbinsX();
   for (int i=1; i<=n/2; i++) {
      double c1 = hist->GetBinContent(i);
      double e1 = hist->GetBinError(i);
      double c2 = hist->GetBinContent(n+1-i);
      double e2 = hist->GetBinError(n+1-i);
      hist->SetBinContent(i,c2);
      hist->SetBinError(i,e2);
      hist->SetBinContent(n+1-i,c1);
      hist->SetBinError(n+1-i,e1);
   }
}

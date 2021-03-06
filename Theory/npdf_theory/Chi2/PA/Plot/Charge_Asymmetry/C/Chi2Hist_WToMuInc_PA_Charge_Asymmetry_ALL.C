void Chi2Hist_WToMuInc_PA_Charge_Asymmetry_ALL()
{
//=========Macro generated from canvas: c/c
//=========  (Thu Mar 28 20:02:04 2019) by ROOT version 6.12/07
   TCanvas *c = new TCanvas("c", "c",0,0,1000,1000);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c->Range(-2.4,-16.5,12.6,121);
   c->SetFillColor(0);
   c->SetBorderMode(0);
   c->SetBorderSize(2);
   c->SetTickx(1);
   c->SetTicky(1);
   c->SetLeftMargin(0.16);
   c->SetRightMargin(0.04);
   c->SetTopMargin(0.08);
   c->SetBottomMargin(0.12);
   c->SetFrameFillStyle(0);
   c->SetFrameBorderMode(0);
   c->SetFrameFillStyle(0);
   c->SetFrameBorderMode(0);
   
   TH1D *EPPS16_TEMP__57 = new TH1D("EPPS16_TEMP__57","",60,0,12);
   EPPS16_TEMP__57->SetBinContent(4,3);
   EPPS16_TEMP__57->SetBinContent(5,71);
   EPPS16_TEMP__57->SetBinContent(6,17);
   EPPS16_TEMP__57->SetBinContent(7,5);
   EPPS16_TEMP__57->SetBinContent(8,1);
   EPPS16_TEMP__57->SetBinError(4,1.732051);
   EPPS16_TEMP__57->SetBinError(5,8.42615);
   EPPS16_TEMP__57->SetBinError(6,4.123106);
   EPPS16_TEMP__57->SetBinError(7,2.236068);
   EPPS16_TEMP__57->SetBinError(8,1);
   EPPS16_TEMP__57->SetMinimum(0);
   EPPS16_TEMP__57->SetMaximum(110);
   EPPS16_TEMP__57->SetEntries(97);
   EPPS16_TEMP__57->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = 2401;
   color = new TColor(ci, 0, 0.6, 0, " ", 0.2);
   EPPS16_TEMP__57->SetFillColor(ci);
   EPPS16_TEMP__57->SetMarkerStyle(20);
   EPPS16_TEMP__57->SetMarkerSize(0);
   EPPS16_TEMP__57->GetXaxis()->SetTitle("#chi^{2}/ndf");
   EPPS16_TEMP__57->GetXaxis()->CenterTitle(true);
   EPPS16_TEMP__57->GetXaxis()->SetLabelFont(42);
   EPPS16_TEMP__57->GetXaxis()->SetLabelOffset(0.007);
   EPPS16_TEMP__57->GetXaxis()->SetLabelSize(0.035);
   EPPS16_TEMP__57->GetXaxis()->SetTitleSize(0.065);
   EPPS16_TEMP__57->GetXaxis()->SetTitleOffset(0.7);
   EPPS16_TEMP__57->GetXaxis()->SetTitleFont(42);
   EPPS16_TEMP__57->GetYaxis()->SetTitle("Number of error sets");
   EPPS16_TEMP__57->GetYaxis()->CenterTitle(true);
   EPPS16_TEMP__57->GetYaxis()->SetLabelFont(42);
   EPPS16_TEMP__57->GetYaxis()->SetLabelOffset(0.007);
   EPPS16_TEMP__57->GetYaxis()->SetLabelSize(0.035);
   EPPS16_TEMP__57->GetYaxis()->SetTitleSize(0.065);
   EPPS16_TEMP__57->GetYaxis()->SetTitleOffset(1.05);
   EPPS16_TEMP__57->GetYaxis()->SetTitleFont(42);
   EPPS16_TEMP__57->GetZaxis()->SetLabelFont(42);
   EPPS16_TEMP__57->GetZaxis()->SetLabelOffset(0.007);
   EPPS16_TEMP__57->GetZaxis()->SetLabelSize(0.05);
   EPPS16_TEMP__57->GetZaxis()->SetTitleSize(0.06);
   EPPS16_TEMP__57->GetZaxis()->SetTitleFont(42);
   EPPS16_TEMP__57->Draw("HIST");
   
   TH1D *nCTEQ15_TEMP__58 = new TH1D("nCTEQ15_TEMP__58","",60,0,12);
   nCTEQ15_TEMP__58->SetBinContent(7,1);
   nCTEQ15_TEMP__58->SetBinContent(8,7);
   nCTEQ15_TEMP__58->SetBinContent(9,52);
   nCTEQ15_TEMP__58->SetBinContent(10,24);
   nCTEQ15_TEMP__58->SetBinContent(11,4);
   nCTEQ15_TEMP__58->SetBinContent(12,1);
   nCTEQ15_TEMP__58->SetBinError(7,1);
   nCTEQ15_TEMP__58->SetBinError(8,2.645751);
   nCTEQ15_TEMP__58->SetBinError(9,7.211103);
   nCTEQ15_TEMP__58->SetBinError(10,4.898979);
   nCTEQ15_TEMP__58->SetBinError(11,2);
   nCTEQ15_TEMP__58->SetBinError(12,1);
   nCTEQ15_TEMP__58->SetMinimum(0);
   nCTEQ15_TEMP__58->SetMaximum(110);
   nCTEQ15_TEMP__58->SetEntries(89);
   nCTEQ15_TEMP__58->SetStats(0);

   ci = 2400;
   color = new TColor(ci, 0, 0, 1, " ", 0.2);
   nCTEQ15_TEMP__58->SetFillColor(ci);
   nCTEQ15_TEMP__58->SetMarkerStyle(20);
   nCTEQ15_TEMP__58->SetMarkerSize(0);
   nCTEQ15_TEMP__58->GetXaxis()->SetTitle("#chi^{2}/ndf");
   nCTEQ15_TEMP__58->GetXaxis()->CenterTitle(true);
   nCTEQ15_TEMP__58->GetXaxis()->SetLabelFont(42);
   nCTEQ15_TEMP__58->GetXaxis()->SetLabelOffset(0.007);
   nCTEQ15_TEMP__58->GetXaxis()->SetLabelSize(0.035);
   nCTEQ15_TEMP__58->GetXaxis()->SetTitleSize(0.065);
   nCTEQ15_TEMP__58->GetXaxis()->SetTitleOffset(0.7);
   nCTEQ15_TEMP__58->GetXaxis()->SetTitleFont(42);
   nCTEQ15_TEMP__58->GetYaxis()->SetTitle("Number of error sets");
   nCTEQ15_TEMP__58->GetYaxis()->CenterTitle(true);
   nCTEQ15_TEMP__58->GetYaxis()->SetLabelFont(42);
   nCTEQ15_TEMP__58->GetYaxis()->SetLabelOffset(0.007);
   nCTEQ15_TEMP__58->GetYaxis()->SetLabelSize(0.035);
   nCTEQ15_TEMP__58->GetYaxis()->SetTitleSize(0.065);
   nCTEQ15_TEMP__58->GetYaxis()->SetTitleOffset(1.05);
   nCTEQ15_TEMP__58->GetYaxis()->SetTitleFont(42);
   nCTEQ15_TEMP__58->GetZaxis()->SetLabelFont(42);
   nCTEQ15_TEMP__58->GetZaxis()->SetLabelOffset(0.007);
   nCTEQ15_TEMP__58->GetZaxis()->SetLabelSize(0.05);
   nCTEQ15_TEMP__58->GetZaxis()->SetTitleSize(0.06);
   nCTEQ15_TEMP__58->GetZaxis()->SetTitleFont(42);
   nCTEQ15_TEMP__58->Draw("HIST same");
   
   TH1D *CT14_TEMP__59 = new TH1D("CT14_TEMP__59","",60,0,12);
   CT14_TEMP__59->SetBinContent(5,13);
   CT14_TEMP__59->SetBinContent(6,41);
   CT14_TEMP__59->SetBinContent(7,3);
   CT14_TEMP__59->SetBinError(5,3.605551);
   CT14_TEMP__59->SetBinError(6,6.403124);
   CT14_TEMP__59->SetBinError(7,1.732051);
   CT14_TEMP__59->SetMinimum(0);
   CT14_TEMP__59->SetMaximum(110);
   CT14_TEMP__59->SetEntries(57);
   CT14_TEMP__59->SetStats(0);

   ci = 2399;
   color = new TColor(ci, 1, 0, 0, " ", 0.2);
   CT14_TEMP__59->SetFillColor(ci);
   CT14_TEMP__59->SetMarkerStyle(20);
   CT14_TEMP__59->SetMarkerSize(0);
   CT14_TEMP__59->GetXaxis()->SetTitle("#chi^{2}/ndf");
   CT14_TEMP__59->GetXaxis()->CenterTitle(true);
   CT14_TEMP__59->GetXaxis()->SetLabelFont(42);
   CT14_TEMP__59->GetXaxis()->SetLabelOffset(0.007);
   CT14_TEMP__59->GetXaxis()->SetLabelSize(0.035);
   CT14_TEMP__59->GetXaxis()->SetTitleSize(0.065);
   CT14_TEMP__59->GetXaxis()->SetTitleOffset(0.7);
   CT14_TEMP__59->GetXaxis()->SetTitleFont(42);
   CT14_TEMP__59->GetYaxis()->SetTitle("Number of error sets");
   CT14_TEMP__59->GetYaxis()->CenterTitle(true);
   CT14_TEMP__59->GetYaxis()->SetLabelFont(42);
   CT14_TEMP__59->GetYaxis()->SetLabelOffset(0.007);
   CT14_TEMP__59->GetYaxis()->SetLabelSize(0.035);
   CT14_TEMP__59->GetYaxis()->SetTitleSize(0.065);
   CT14_TEMP__59->GetYaxis()->SetTitleOffset(1.05);
   CT14_TEMP__59->GetYaxis()->SetTitleFont(42);
   CT14_TEMP__59->GetZaxis()->SetLabelFont(42);
   CT14_TEMP__59->GetZaxis()->SetLabelOffset(0.007);
   CT14_TEMP__59->GetZaxis()->SetLabelSize(0.05);
   CT14_TEMP__59->GetZaxis()->SetTitleSize(0.06);
   CT14_TEMP__59->GetZaxis()->SetTitleFont(42);
   CT14_TEMP__59->Draw("HIST same");
   
   TLegend *leg = new TLegend(0.58,0.54,0.78,0.74,NULL,"brNDC");
   leg->SetTextFont(62);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("CT14_TEMP","CT14","f");

   ci = 2399;
   color = new TColor(ci, 1, 0, 0, " ", 0.2);
   entry->SetFillColor(ci);
   entry->SetFillStyle(1001);
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(62);
   entry->SetTextSize(0.04);
   entry=leg->AddEntry("EPPS16_TEMP","(CT14+)EPPS16","f");

   ci = 2401;
   color = new TColor(ci, 0, 0.6, 0, " ", 0.2);
   entry->SetFillColor(ci);
   entry->SetFillStyle(1001);
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(62);
   entry->SetTextSize(0.04);
   entry=leg->AddEntry("nCTEQ15_TEMP","(CT14+)nCTEQ15","f");

   ci = 2400;
   color = new TColor(ci, 0, 0, 1, " ", 0.2);
   entry->SetFillColor(ci);
   entry->SetFillStyle(1001);
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(62);
   entry->SetTextSize(0.04);
   leg->Draw();
   TLine *line = new TLine(1.034392,0,1.034392,50);

   ci = TColor::GetColor("#ff0000");
   line->SetLineColor(ci);
   line->SetLineStyle(7);
   line->SetLineWidth(3);
   line->Draw();
   line = new TLine(0.8984559,0,0.8984559,50);

   ci = TColor::GetColor("#009900");
   line->SetLineColor(ci);
   line->SetLineStyle(7);
   line->SetLineWidth(3);
   line->Draw();
   line = new TLine(1.793194,0,1.793194,50);

   ci = TColor::GetColor("#0000ff");
   line->SetLineColor(ci);
   line->SetLineStyle(7);
   line->SetLineWidth(3);
   line->Draw();
   TLatex *   tex = new TLatex(0.22,0.84,"W#kern[0.2]{#rightarrow}#kern[0.2]{#mu}#kern[0.2]{#nu_{#mu}}");
tex->SetNDC();
   tex->SetTextSize(0.055);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.78,0.84,"CMS");
tex->SetNDC();
   tex->SetTextFont(61);
   tex->SetTextSize(0.058);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.22,0.725,"p^{#mu}_{T} > 25 GeV/c");
tex->SetNDC();
   tex->SetTextSize(0.035);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.22,0.775,"( N_{#mu}^{+} #font[122]{-} N_{#mu}^{#font[122]{-}} ) / ( N_{#mu}^{+} + N_{#mu}^{#font[122]{-}} )");
tex->SetNDC();
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.96,0.9424,"pPb 173.4 nb^{-1}             #sqrt{s_{NN}} = 8.16 TeV");
tex->SetNDC();
   tex->SetTextAlign(31);
   tex->SetTextFont(42);
   tex->SetTextSize(0.048);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.16,0.9424,"");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.048);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.914,0.874,"");
tex->SetNDC();
   tex->SetTextAlign(33);
   tex->SetTextFont(61);
   tex->SetTextSize(0.06);
   tex->SetLineWidth(2);
   tex->Draw();
   c->Modified();
   c->cd();
   c->SetSelected(c);
}

void Wm()
{
//=========Macro generated from canvas: Wm/
//=========  (Fri Apr  6 14:52:00 2018) by ROOT version6.06/00
   TCanvas *Wm = new TCanvas("Wm", "",0,0,600,600);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   Wm->Range(-3.78688,-30,2.26112,220);
   Wm->SetFillColor(0);
   Wm->SetBorderMode(0);
   Wm->SetBorderSize(2);
   Wm->SetTickx(1);
   Wm->SetTicky(1);
   Wm->SetLeftMargin(0.16);
   Wm->SetRightMargin(0.04);
   Wm->SetTopMargin(0.08);
   Wm->SetBottomMargin(0.12);
   Wm->SetFrameFillStyle(0);
   Wm->SetFrameBorderMode(0);
   Wm->SetFrameFillStyle(0);
   Wm->SetFrameBorderMode(0);
   
   Double_t hrebin_ottriune_pdf_fx3002[24] = {
   -2.7,
   -2.5,
   -2.3,
   -2.1,
   -1.9,
   -1.7,
   -1.5,
   -1.3,
   -1.1,
   -0.9,
   -0.7,
   -0.5,
   -0.3,
   -0.1,
   0.1,
   0.3,
   0.5,
   0.7,
   0.9,
   1.1,
   1.3,
   1.5,
   1.7,
   1.9};
   Double_t hrebin_ottriune_pdf_fy3002[24] = {
   109.5147,
   113.2587,
   115.9851,
   117.7109,
   119.5455,
   120.9316,
   121.5356,
   121.5387,
   121.4287,
   120.4911,
   119.2498,
   117.8252,
   115.3809,
   113.1418,
   110.5728,
   107.0026,
   103.644,
   100.3832,
   96.94868,
   92.74773,
   88.98379,
   85.26247,
   81.17089,
   77.03313};
   Double_t hrebin_ottriune_pdf_felx3002[24] = {
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1};
   Double_t hrebin_ottriune_pdf_fely3002[24] = {
   3.842444,
   4.500608,
   3.666764,
   3.327953,
   3.209269,
   3.794825,
   3.695264,
   3.49067,
   4.962109,
   5.032679,
   6.153435,
   6.870237,
   7.110426,
   7.823986,
   7.971544,
   8.425935,
   8.372711,
   8.275661,
   8.46904,
   7.96289,
   8.303458,
   7.411268,
   7.364643,
   6.84928};
   Double_t hrebin_ottriune_pdf_fehx3002[24] = {
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1};
   Double_t hrebin_ottriune_pdf_fehy3002[24] = {
   4.41122,
   3.911194,
   4.599979,
   3.671759,
   3.179768,
   1.959384,
   2.114329,
   2.439649,
   2.12181,
   3.165043,
   3.068359,
   3.7715,
   4.327594,
   4.377491,
   5.181433,
   5.690239,
   5.263571,
   5.793194,
   5.635074,
   6.076509,
   5.48455,
   5.749824,
   5.469996,
   5.785023};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(24,hrebin_ottriune_pdf_fx3002,hrebin_ottriune_pdf_fy3002,hrebin_ottriune_pdf_felx3002,hrebin_ottriune_pdf_fehx3002,hrebin_ottriune_pdf_fely3002,hrebin_ottriune_pdf_fehy3002);
   grae->SetName("hrebin_ottriune_pdf");
   grae->SetTitle("Graph");

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#ff00ff");
   grae->SetFillColor(ci);
   grae->SetFillStyle(3007);

   ci = TColor::GetColor("#ff00ff");
   grae->SetLineColor(ci);

   ci = TColor::GetColor("#ff00ff");
   grae->SetMarkerColor(ci);
   
   TH1F *Graph_hrebin_ottriune_pdf3002 = new TH1F("Graph_hrebin_ottriune_pdf3002","Graph",100,-3.28,2.48);
   Graph_hrebin_ottriune_pdf3002->SetMinimum(0);
   Graph_hrebin_ottriune_pdf3002->SetMaximum(200);
   Graph_hrebin_ottriune_pdf3002->SetDirectory(0);
   Graph_hrebin_ottriune_pdf3002->SetStats(0);
   Graph_hrebin_ottriune_pdf3002->SetLineStyle(0);
   Graph_hrebin_ottriune_pdf3002->SetMarkerStyle(20);
   Graph_hrebin_ottriune_pdf3002->GetXaxis()->SetTitle("#eta_{cm}");
   Graph_hrebin_ottriune_pdf3002->GetXaxis()->SetRange(9,92);
   Graph_hrebin_ottriune_pdf3002->GetXaxis()->SetLabelFont(42);
   Graph_hrebin_ottriune_pdf3002->GetXaxis()->SetLabelOffset(0.007);
   Graph_hrebin_ottriune_pdf3002->GetXaxis()->SetLabelSize(0.05);
   Graph_hrebin_ottriune_pdf3002->GetXaxis()->SetTitleSize(0.06);
   Graph_hrebin_ottriune_pdf3002->GetXaxis()->SetTitleOffset(0.9);
   Graph_hrebin_ottriune_pdf3002->GetXaxis()->SetTitleFont(42);
   Graph_hrebin_ottriune_pdf3002->GetYaxis()->SetTitle("d#sigma(W^{-} #rightarrow l^{-}#nu / d#eta_{cm} [nb]");
   Graph_hrebin_ottriune_pdf3002->GetYaxis()->SetLabelFont(42);
   Graph_hrebin_ottriune_pdf3002->GetYaxis()->SetLabelOffset(0.007);
   Graph_hrebin_ottriune_pdf3002->GetYaxis()->SetLabelSize(0.05);
   Graph_hrebin_ottriune_pdf3002->GetYaxis()->SetTitleSize(0.06);
   Graph_hrebin_ottriune_pdf3002->GetYaxis()->SetTitleOffset(1.25);
   Graph_hrebin_ottriune_pdf3002->GetYaxis()->SetTitleFont(42);
   Graph_hrebin_ottriune_pdf3002->GetZaxis()->SetLabelFont(42);
   Graph_hrebin_ottriune_pdf3002->GetZaxis()->SetLabelOffset(0.007);
   Graph_hrebin_ottriune_pdf3002->GetZaxis()->SetLabelSize(0.05);
   Graph_hrebin_ottriune_pdf3002->GetZaxis()->SetTitleSize(0.06);
   Graph_hrebin_ottriune_pdf3002->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_hrebin_ottriune_pdf3002);
   
   grae->Draw("a5");
   
   TLegend *leg = new TLegend(0.55,0.17,0.88,0.4,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(62);
   leg->SetTextSize(0.04);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("NULL","MCFM nlo","h");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(62);
   entry=leg->AddEntry("hrebin_ottriune_pdf","nCTEQ15_208_82","lpf");

   ci = TColor::GetColor("#ff00ff");
   entry->SetFillColor(ci);
   entry->SetFillStyle(3007);

   ci = TColor::GetColor("#ff00ff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#ff00ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(1);
   entry->SetMarkerSize(1);
   entry->SetTextFont(62);
   leg->Draw();
   TLatex *   tex = new TLatex(0.96,0.9424,"pA [285410-286504] (8.16 TeV)");
tex->SetNDC();
   tex->SetTextAlign(31);
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.16,0.9424,"");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.186,0.874,"CMS");
tex->SetNDC();
   tex->SetTextAlign(13);
   tex->SetTextFont(61);
   tex->SetTextSize(0.06);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.186,0.802,"Projection");
tex->SetNDC();
   tex->SetTextAlign(13);
   tex->SetTextFont(52);
   tex->SetTextSize(0.0456);
   tex->SetLineWidth(2);
   tex->Draw();
   Wm->Modified();
   Wm->cd();
   Wm->SetSelected(Wm);
}

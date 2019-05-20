void Wp()
{
//=========Macro generated from canvas: Wp/
//=========  (Fri Apr  6 14:52:00 2018) by ROOT version6.06/00
   TCanvas *Wp = new TCanvas("Wp", "",0,0,600,600);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   Wp->Range(-3.78688,-30,2.26112,220);
   Wp->SetFillColor(0);
   Wp->SetBorderMode(0);
   Wp->SetBorderSize(2);
   Wp->SetTickx(1);
   Wp->SetTicky(1);
   Wp->SetLeftMargin(0.16);
   Wp->SetRightMargin(0.04);
   Wp->SetTopMargin(0.08);
   Wp->SetBottomMargin(0.12);
   Wp->SetFrameFillStyle(0);
   Wp->SetFrameBorderMode(0);
   Wp->SetFrameFillStyle(0);
   Wp->SetFrameBorderMode(0);
   
   Double_t hrebin_zehygmyt_pdf_fx3001[24] = {
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
   Double_t hrebin_zehygmyt_pdf_fy3001[24] = {
   101.2194,
   114.303,
   123.5385,
   131.5611,
   135.5926,
   138.9257,
   140.2984,
   140.8792,
   140.4074,
   139.8834,
   138.9021,
   136.6738,
   135.6048,
   132.8683,
   131.641,
   129.8863,
   128.5495,
   127.3229,
   125.8323,
   125.491,
   124.3651,
   123.5164,
   122.7035,
   120.5999};
   Double_t hrebin_zehygmyt_pdf_felx3001[24] = {
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
   Double_t hrebin_zehygmyt_pdf_fely3001[24] = {
   3.670377,
   4.927807,
   4.050295,
   5.679569,
   2.769592,
   4.22753,
   3.839958,
   5.123603,
   5.145812,
   5.583397,
   7.356064,
   7.011968,
   8.8553,
   8.509213,
   9.802973,
   10.01182,
   9.943668,
   10.76481,
   10.04339,
   10.83006,
   10.89878,
   10.89474,
   11.2579,
   10.52425};
   Double_t hrebin_zehygmyt_pdf_fehx3001[24] = {
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
   Double_t hrebin_zehygmyt_pdf_fehy3001[24] = {
   4.518036,
   3.820741,
   4.804236,
   2.939747,
   5.024149,
   3.196404,
   3.75059,
   2.509411,
   3.211668,
   3.593604,
   3.652106,
   4.826737,
   4.684279,
   6.110307,
   5.792998,
   6.344376,
   6.580945,
   7.314928,
   7.807191,
   7.971424,
   8.151022,
   8.304945,
   8.408161,
   8.946451};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(24,hrebin_zehygmyt_pdf_fx3001,hrebin_zehygmyt_pdf_fy3001,hrebin_zehygmyt_pdf_felx3001,hrebin_zehygmyt_pdf_fehx3001,hrebin_zehygmyt_pdf_fely3001,hrebin_zehygmyt_pdf_fehy3001);
   grae->SetName("hrebin_zehygmyt_pdf");
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
   
   TH1F *Graph_hrebin_zehygmyt_pdf3001 = new TH1F("Graph_hrebin_zehygmyt_pdf3001","Graph",100,-3.28,2.48);
   Graph_hrebin_zehygmyt_pdf3001->SetMinimum(0);
   Graph_hrebin_zehygmyt_pdf3001->SetMaximum(200);
   Graph_hrebin_zehygmyt_pdf3001->SetDirectory(0);
   Graph_hrebin_zehygmyt_pdf3001->SetStats(0);
   Graph_hrebin_zehygmyt_pdf3001->SetLineStyle(0);
   Graph_hrebin_zehygmyt_pdf3001->SetMarkerStyle(20);
   Graph_hrebin_zehygmyt_pdf3001->GetXaxis()->SetTitle("#eta_{cm}");
   Graph_hrebin_zehygmyt_pdf3001->GetXaxis()->SetRange(9,92);
   Graph_hrebin_zehygmyt_pdf3001->GetXaxis()->SetLabelFont(42);
   Graph_hrebin_zehygmyt_pdf3001->GetXaxis()->SetLabelOffset(0.007);
   Graph_hrebin_zehygmyt_pdf3001->GetXaxis()->SetLabelSize(0.05);
   Graph_hrebin_zehygmyt_pdf3001->GetXaxis()->SetTitleSize(0.06);
   Graph_hrebin_zehygmyt_pdf3001->GetXaxis()->SetTitleOffset(0.9);
   Graph_hrebin_zehygmyt_pdf3001->GetXaxis()->SetTitleFont(42);
   Graph_hrebin_zehygmyt_pdf3001->GetYaxis()->SetTitle("d#sigma(W^{+} #rightarrow l^{+}#nu / d#eta_{cm} [nb]");
   Graph_hrebin_zehygmyt_pdf3001->GetYaxis()->SetLabelFont(42);
   Graph_hrebin_zehygmyt_pdf3001->GetYaxis()->SetLabelOffset(0.007);
   Graph_hrebin_zehygmyt_pdf3001->GetYaxis()->SetLabelSize(0.05);
   Graph_hrebin_zehygmyt_pdf3001->GetYaxis()->SetTitleSize(0.06);
   Graph_hrebin_zehygmyt_pdf3001->GetYaxis()->SetTitleOffset(1.25);
   Graph_hrebin_zehygmyt_pdf3001->GetYaxis()->SetTitleFont(42);
   Graph_hrebin_zehygmyt_pdf3001->GetZaxis()->SetLabelFont(42);
   Graph_hrebin_zehygmyt_pdf3001->GetZaxis()->SetLabelOffset(0.007);
   Graph_hrebin_zehygmyt_pdf3001->GetZaxis()->SetLabelSize(0.05);
   Graph_hrebin_zehygmyt_pdf3001->GetZaxis()->SetTitleSize(0.06);
   Graph_hrebin_zehygmyt_pdf3001->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_hrebin_zehygmyt_pdf3001);
   
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
   entry=leg->AddEntry("hrebin_zehygmyt_pdf","nCTEQ15_208_82","lpf");

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
   Wp->Modified();
   Wp->cd();
   Wp->SetSelected(Wp);
}

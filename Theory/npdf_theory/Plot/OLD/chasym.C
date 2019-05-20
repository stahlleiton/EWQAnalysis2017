void chasym()
{
//=========Macro generated from canvas: chasym/
//=========  (Fri Apr  6 14:52:00 2018) by ROOT version6.06/00
   TCanvas *chasym = new TCanvas("chasym", "",0,0,600,600);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   chasym->Range(-3.78688,-0.275,2.26112,0.35);
   chasym->SetFillColor(0);
   chasym->SetBorderMode(0);
   chasym->SetBorderSize(2);
   chasym->SetTickx(1);
   chasym->SetTicky(1);
   chasym->SetLeftMargin(0.16);
   chasym->SetRightMargin(0.04);
   chasym->SetTopMargin(0.08);
   chasym->SetBottomMargin(0.12);
   chasym->SetFrameFillStyle(0);
   chasym->SetFrameBorderMode(0);
   chasym->SetFrameFillStyle(0);
   chasym->SetFrameBorderMode(0);
   
   Double_t hasym_bqsqscsc_pdf_fx3003[24] = {
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
   Double_t hasym_bqsqscsc_pdf_fy3003[24] = {
   -0.03936382,
   0.004588801,
   0.03153526,
   0.05556284,
   0.06289559,
   0.06924605,
   0.07165899,
   0.07370086,
   0.07248313,
   0.07447848,
   0.07612707,
   0.07406124,
   0.08057776,
   0.08018577,
   0.08698174,
   0.09660129,
   0.1072621,
   0.1183089,
   0.1296505,
   0.150034,
   0.1658379,
   0.183227,
   0.2037166,
   0.2204427};
   Double_t hasym_bqsqscsc_pdf_felx3003[24] = {
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
   Double_t hasym_bqsqscsc_pdf_fely3003[24] = {
   0.01493304,
   0.01919298,
   0.01747401,
   0.0183904,
   0.006019536,
   0.008177934,
   0.00582284,
   0.009691705,
   0.004749256,
   0.006587645,
   0.006337597,
   0.002532932,
   0.005545626,
   0.001224936,
   0.005054263,
   0.005242697,
   0.002143763,
   0.006265048,
   0.002824897,
   0.00710641,
   0.003588946,
   0.006863122,
   0.01197798,
   0.007960514};
   Double_t hasym_bqsqscsc_pdf_fehx3003[24] = {
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
   Double_t hasym_bqsqscsc_pdf_fehy3003[24] = {
   0.01684826,
   0.01578628,
   0.01434179,
   0.006474419,
   0.01325006,
   0.009910128,
   0.01029839,
   0.003635874,
   0.00777095,
   0.005923162,
   0.004594933,
   0.006419378,
   0.003067948,
   0.007800263,
   0.002395113,
   0.003062959,
   0.004959777,
   0.003705798,
   0.008247906,
   0.004945033,
   0.01029107,
   0.007429249,
   0.008985005,
   0.007577535};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(24,hasym_bqsqscsc_pdf_fx3003,hasym_bqsqscsc_pdf_fy3003,hasym_bqsqscsc_pdf_felx3003,hasym_bqsqscsc_pdf_fehx3003,hasym_bqsqscsc_pdf_fely3003,hasym_bqsqscsc_pdf_fehy3003);
   grae->SetName("hasym_bqsqscsc_pdf");
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
   
   TH1F *Graph_hasym_bqsqscsc_pdf3003 = new TH1F("Graph_hasym_bqsqscsc_pdf3003","Graph",100,-3.28,2.48);
   Graph_hasym_bqsqscsc_pdf3003->SetMinimum(-0.2);
   Graph_hasym_bqsqscsc_pdf3003->SetMaximum(0.3);
   Graph_hasym_bqsqscsc_pdf3003->SetDirectory(0);
   Graph_hasym_bqsqscsc_pdf3003->SetStats(0);
   Graph_hasym_bqsqscsc_pdf3003->SetLineStyle(0);
   Graph_hasym_bqsqscsc_pdf3003->SetMarkerStyle(20);
   Graph_hasym_bqsqscsc_pdf3003->GetXaxis()->SetTitle("#eta_{cm}");
   Graph_hasym_bqsqscsc_pdf3003->GetXaxis()->SetRange(9,92);
   Graph_hasym_bqsqscsc_pdf3003->GetXaxis()->SetLabelFont(42);
   Graph_hasym_bqsqscsc_pdf3003->GetXaxis()->SetLabelOffset(0.007);
   Graph_hasym_bqsqscsc_pdf3003->GetXaxis()->SetLabelSize(0.05);
   Graph_hasym_bqsqscsc_pdf3003->GetXaxis()->SetTitleSize(0.06);
   Graph_hasym_bqsqscsc_pdf3003->GetXaxis()->SetTitleOffset(0.9);
   Graph_hasym_bqsqscsc_pdf3003->GetXaxis()->SetTitleFont(42);
   Graph_hasym_bqsqscsc_pdf3003->GetYaxis()->SetTitle("(N^{+} - N^{-}) / (N^{+} + N^{-})");
   Graph_hasym_bqsqscsc_pdf3003->GetYaxis()->SetLabelFont(42);
   Graph_hasym_bqsqscsc_pdf3003->GetYaxis()->SetLabelOffset(0.007);
   Graph_hasym_bqsqscsc_pdf3003->GetYaxis()->SetLabelSize(0.05);
   Graph_hasym_bqsqscsc_pdf3003->GetYaxis()->SetTitleSize(0.06);
   Graph_hasym_bqsqscsc_pdf3003->GetYaxis()->SetTitleOffset(1.25);
   Graph_hasym_bqsqscsc_pdf3003->GetYaxis()->SetTitleFont(42);
   Graph_hasym_bqsqscsc_pdf3003->GetZaxis()->SetLabelFont(42);
   Graph_hasym_bqsqscsc_pdf3003->GetZaxis()->SetLabelOffset(0.007);
   Graph_hasym_bqsqscsc_pdf3003->GetZaxis()->SetLabelSize(0.05);
   Graph_hasym_bqsqscsc_pdf3003->GetZaxis()->SetTitleSize(0.06);
   Graph_hasym_bqsqscsc_pdf3003->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_hasym_bqsqscsc_pdf3003);
   
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
   entry=leg->AddEntry("hasym_bqsqscsc_pdf","nCTEQ15_208_82","lpf");

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
   chasym->Modified();
   chasym->cd();
   chasym->SetSelected(chasym);
}

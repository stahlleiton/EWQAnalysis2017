void A1p()
{
//=========Macro generated from canvas: A1p/
//=========  (Fri Apr  6 14:52:00 2018) by ROOT version6.06/00
   TCanvas *A1p = new TCanvas("A1p", "",0,0,600,600);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   A1p->Range(-0.44,0.68,2.31,1.68);
   A1p->SetFillColor(0);
   A1p->SetBorderMode(0);
   A1p->SetBorderSize(2);
   A1p->SetTickx(1);
   A1p->SetTicky(1);
   A1p->SetLeftMargin(0.16);
   A1p->SetRightMargin(0.04);
   A1p->SetTopMargin(0.08);
   A1p->SetBottomMargin(0.12);
   A1p->SetFrameFillStyle(0);
   A1p->SetFrameBorderMode(0);
   A1p->SetFrameFillStyle(0);
   A1p->SetFrameBorderMode(0);
   
   Double_t ha1p_hxfldkre_pdf_fx3004[10] = {
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
   Double_t ha1p_hxfldkre_pdf_fy3004[10] = {
   0.9907628,
   0.9578303,
   0.9405574,
   0.9166374,
   0.8995518,
   0.8937631,
   0.8827789,
   0.8803835,
   0.8832308,
   0.8894285};
   Double_t ha1p_hxfldkre_pdf_felx3004[10] = {
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
   Double_t ha1p_hxfldkre_pdf_fely3004[10] = {
   0.01525227,
   0.01670076,
   0.0311632,
   0.04049268,
   0.04699116,
   0.062318,
   0.06395599,
   0.08116489,
   0.08450638,
   0.09725654};
   Double_t ha1p_hxfldkre_pdf_fehx3004[10] = {
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
   Double_t ha1p_hxfldkre_pdf_fehy3004[10] = {
   0.003344126,
   0.02265263,
   0.02763441,
   0.04459985,
   0.05268218,
   0.06131672,
   0.06684606,
   0.06982223,
   0.07513433,
   0.07587482};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(10,ha1p_hxfldkre_pdf_fx3004,ha1p_hxfldkre_pdf_fy3004,ha1p_hxfldkre_pdf_felx3004,ha1p_hxfldkre_pdf_fehx3004,ha1p_hxfldkre_pdf_fely3004,ha1p_hxfldkre_pdf_fehy3004);
   grae->SetName("ha1p_hxfldkre_pdf");
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
   
   TH1F *Graph_ha1p_hxfldkre_pdf3004 = new TH1F("Graph_ha1p_hxfldkre_pdf3004","Graph",100,0,2.2);
   Graph_ha1p_hxfldkre_pdf3004->SetMinimum(0.8);
   Graph_ha1p_hxfldkre_pdf3004->SetMaximum(1.6);
   Graph_ha1p_hxfldkre_pdf3004->SetDirectory(0);
   Graph_ha1p_hxfldkre_pdf3004->SetStats(0);
   Graph_ha1p_hxfldkre_pdf3004->SetLineStyle(0);
   Graph_ha1p_hxfldkre_pdf3004->SetMarkerStyle(20);
   Graph_ha1p_hxfldkre_pdf3004->GetXaxis()->SetTitle("#eta_{cm}");
   Graph_ha1p_hxfldkre_pdf3004->GetXaxis()->SetRange(1,100);
   Graph_ha1p_hxfldkre_pdf3004->GetXaxis()->SetLabelFont(42);
   Graph_ha1p_hxfldkre_pdf3004->GetXaxis()->SetLabelOffset(0.007);
   Graph_ha1p_hxfldkre_pdf3004->GetXaxis()->SetLabelSize(0.05);
   Graph_ha1p_hxfldkre_pdf3004->GetXaxis()->SetTitleSize(0.06);
   Graph_ha1p_hxfldkre_pdf3004->GetXaxis()->SetTitleOffset(0.9);
   Graph_ha1p_hxfldkre_pdf3004->GetXaxis()->SetTitleFont(42);
   Graph_ha1p_hxfldkre_pdf3004->GetYaxis()->SetTitle("N^{+} (+#eta_{cm}) / N^{+} (-#eta_{cm})");
   Graph_ha1p_hxfldkre_pdf3004->GetYaxis()->SetLabelFont(42);
   Graph_ha1p_hxfldkre_pdf3004->GetYaxis()->SetLabelOffset(0.007);
   Graph_ha1p_hxfldkre_pdf3004->GetYaxis()->SetLabelSize(0.05);
   Graph_ha1p_hxfldkre_pdf3004->GetYaxis()->SetTitleSize(0.06);
   Graph_ha1p_hxfldkre_pdf3004->GetYaxis()->SetTitleOffset(1.25);
   Graph_ha1p_hxfldkre_pdf3004->GetYaxis()->SetTitleFont(42);
   Graph_ha1p_hxfldkre_pdf3004->GetZaxis()->SetLabelFont(42);
   Graph_ha1p_hxfldkre_pdf3004->GetZaxis()->SetLabelOffset(0.007);
   Graph_ha1p_hxfldkre_pdf3004->GetZaxis()->SetLabelSize(0.05);
   Graph_ha1p_hxfldkre_pdf3004->GetZaxis()->SetTitleSize(0.06);
   Graph_ha1p_hxfldkre_pdf3004->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_ha1p_hxfldkre_pdf3004);
   
   grae->Draw("a5");
   
   TLegend *leg = new TLegend(0.2,0.63,0.55,0.86,NULL,"brNDC");
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
   entry=leg->AddEntry("ha1p_hxfldkre_pdf","nCTEQ15_208_82","lpf");

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
   A1p->Modified();
   A1p->cd();
   A1p->SetSelected(A1p);
}

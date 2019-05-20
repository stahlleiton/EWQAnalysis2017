void A1m()
{
//=========Macro generated from canvas: A1m/
//=========  (Fri Apr  6 14:52:00 2018) by ROOT version6.06/00
   TCanvas *A1m = new TCanvas("A1m", "",0,0,600,600);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   A1m->Range(-0.44,0.525,2.31,1.15);
   A1m->SetFillColor(0);
   A1m->SetBorderMode(0);
   A1m->SetBorderSize(2);
   A1m->SetTickx(1);
   A1m->SetTicky(1);
   A1m->SetLeftMargin(0.16);
   A1m->SetRightMargin(0.04);
   A1m->SetTopMargin(0.08);
   A1m->SetBottomMargin(0.12);
   A1m->SetFrameFillStyle(0);
   A1m->SetFrameBorderMode(0);
   A1m->SetFrameFillStyle(0);
   A1m->SetFrameBorderMode(0);
   
   Double_t ha1p_elqrsmci_pdf_fx3005[10] = {
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
   Double_t ha1p_elqrsmci_pdf_fy3005[10] = {
   0.977294,
   0.9273855,
   0.8796414,
   0.8417896,
   0.8046129,
   0.7638041,
   0.7321435,
   0.7015429,
   0.6712131,
   0.6443833};
   Double_t ha1p_elqrsmci_pdf_felx3005[10] = {
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
   Double_t ha1p_elqrsmci_pdf_fely3005[10] = {
   0.005928844,
   0.02116911,
   0.02623285,
   0.03602371,
   0.05091286,
   0.04893767,
   0.06275518,
   0.05945962,
   0.06342108,
   0.06795588};
   Double_t ha1p_elqrsmci_pdf_fehx3005[10] = {
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
   Double_t ha1p_elqrsmci_pdf_fehy3005[10] = {
   0.01215519,
   0.02393009,
   0.02869467,
   0.0393064,
   0.04673432,
   0.05599812,
   0.05091258,
   0.05897219,
   0.05917078,
   0.06084942};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(10,ha1p_elqrsmci_pdf_fx3005,ha1p_elqrsmci_pdf_fy3005,ha1p_elqrsmci_pdf_felx3005,ha1p_elqrsmci_pdf_fehx3005,ha1p_elqrsmci_pdf_fely3005,ha1p_elqrsmci_pdf_fehy3005);
   grae->SetName("ha1p_elqrsmci_pdf");
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
   
   TH1F *Graph_ha1p_elqrsmci_pdf3005 = new TH1F("Graph_ha1p_elqrsmci_pdf3005","Graph",100,0,2.2);
   Graph_ha1p_elqrsmci_pdf3005->SetMinimum(0.6);
   Graph_ha1p_elqrsmci_pdf3005->SetMaximum(1.1);
   Graph_ha1p_elqrsmci_pdf3005->SetDirectory(0);
   Graph_ha1p_elqrsmci_pdf3005->SetStats(0);
   Graph_ha1p_elqrsmci_pdf3005->SetLineStyle(0);
   Graph_ha1p_elqrsmci_pdf3005->SetMarkerStyle(20);
   Graph_ha1p_elqrsmci_pdf3005->GetXaxis()->SetTitle("#eta_{cm}");
   Graph_ha1p_elqrsmci_pdf3005->GetXaxis()->SetRange(1,100);
   Graph_ha1p_elqrsmci_pdf3005->GetXaxis()->SetLabelFont(42);
   Graph_ha1p_elqrsmci_pdf3005->GetXaxis()->SetLabelOffset(0.007);
   Graph_ha1p_elqrsmci_pdf3005->GetXaxis()->SetLabelSize(0.05);
   Graph_ha1p_elqrsmci_pdf3005->GetXaxis()->SetTitleSize(0.06);
   Graph_ha1p_elqrsmci_pdf3005->GetXaxis()->SetTitleOffset(0.9);
   Graph_ha1p_elqrsmci_pdf3005->GetXaxis()->SetTitleFont(42);
   Graph_ha1p_elqrsmci_pdf3005->GetYaxis()->SetTitle("N^{-} (+#eta_{cm}) / N^{-} (-#eta_{cm})");
   Graph_ha1p_elqrsmci_pdf3005->GetYaxis()->SetLabelFont(42);
   Graph_ha1p_elqrsmci_pdf3005->GetYaxis()->SetLabelOffset(0.007);
   Graph_ha1p_elqrsmci_pdf3005->GetYaxis()->SetLabelSize(0.05);
   Graph_ha1p_elqrsmci_pdf3005->GetYaxis()->SetTitleSize(0.06);
   Graph_ha1p_elqrsmci_pdf3005->GetYaxis()->SetTitleOffset(1.25);
   Graph_ha1p_elqrsmci_pdf3005->GetYaxis()->SetTitleFont(42);
   Graph_ha1p_elqrsmci_pdf3005->GetZaxis()->SetLabelFont(42);
   Graph_ha1p_elqrsmci_pdf3005->GetZaxis()->SetLabelOffset(0.007);
   Graph_ha1p_elqrsmci_pdf3005->GetZaxis()->SetLabelSize(0.05);
   Graph_ha1p_elqrsmci_pdf3005->GetZaxis()->SetTitleSize(0.06);
   Graph_ha1p_elqrsmci_pdf3005->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_ha1p_elqrsmci_pdf3005);
   
   grae->Draw("a5");
   
   TLegend *leg = new TLegend(0.2,0.17,0.55,0.4,NULL,"brNDC");
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
   entry=leg->AddEntry("ha1p_elqrsmci_pdf","nCTEQ15_208_82","lpf");

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
   A1m->Modified();
   A1m->cd();
   A1m->SetSelected(A1m);
}

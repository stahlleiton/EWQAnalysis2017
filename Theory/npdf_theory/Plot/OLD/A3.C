void A3()
{
//=========Macro generated from canvas: A3/
//=========  (Fri Apr  6 14:52:00 2018) by ROOT version6.06/00
   TCanvas *A3 = new TCanvas("A3", "",0,0,600,600);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   A3->Range(-0.44,0.6175,2.31,1.305);
   A3->SetFillColor(0);
   A3->SetBorderMode(0);
   A3->SetBorderSize(2);
   A3->SetTickx(1);
   A3->SetTicky(1);
   A3->SetLeftMargin(0.16);
   A3->SetRightMargin(0.04);
   A3->SetTopMargin(0.08);
   A3->SetBottomMargin(0.12);
   A3->SetFrameFillStyle(0);
   A3->SetFrameBorderMode(0);
   A3->SetFrameFillStyle(0);
   A3->SetFrameBorderMode(0);
   
   Double_t ha1p_pomsdbid_pdf_fx3006[10] = {
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
   Double_t ha1p_pomsdbid_pdf_fy3006[10] = {
   0.9845684,
   0.9438345,
   0.9123552,
   0.8820625,
   0.8556178,
   0.8334935,
   0.8130121,
   0.797371,
   0.7845626,
   0.774612};
   Double_t ha1p_pomsdbid_pdf_felx3006[10] = {
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
   Double_t ha1p_pomsdbid_pdf_fely3006[10] = {
   0.009692578,
   0.01817533,
   0.02837802,
   0.03780693,
   0.04815222,
   0.05518327,
   0.06280921,
   0.06984143,
   0.07294086,
   0.08214058};
   Double_t ha1p_pomsdbid_pdf_fehx3006[10] = {
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
   Double_t ha1p_pomsdbid_pdf_fehy3006[10] = {
   0.006132722,
   0.02284595,
   0.02798617,
   0.04154682,
   0.04945012,
   0.05825721,
   0.05850316,
   0.06408042,
   0.06674231,
   0.06793937};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(10,ha1p_pomsdbid_pdf_fx3006,ha1p_pomsdbid_pdf_fy3006,ha1p_pomsdbid_pdf_felx3006,ha1p_pomsdbid_pdf_fehx3006,ha1p_pomsdbid_pdf_fely3006,ha1p_pomsdbid_pdf_fehy3006);
   grae->SetName("ha1p_pomsdbid_pdf");
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
   
   TH1F *Graph_ha1p_pomsdbid_pdf3006 = new TH1F("Graph_ha1p_pomsdbid_pdf3006","Graph",100,0,2.2);
   Graph_ha1p_pomsdbid_pdf3006->SetMinimum(0.7);
   Graph_ha1p_pomsdbid_pdf3006->SetMaximum(1.25);
   Graph_ha1p_pomsdbid_pdf3006->SetDirectory(0);
   Graph_ha1p_pomsdbid_pdf3006->SetStats(0);
   Graph_ha1p_pomsdbid_pdf3006->SetLineStyle(0);
   Graph_ha1p_pomsdbid_pdf3006->SetMarkerStyle(20);
   Graph_ha1p_pomsdbid_pdf3006->GetXaxis()->SetTitle("#eta_{cm}");
   Graph_ha1p_pomsdbid_pdf3006->GetXaxis()->SetRange(1,100);
   Graph_ha1p_pomsdbid_pdf3006->GetXaxis()->SetLabelFont(42);
   Graph_ha1p_pomsdbid_pdf3006->GetXaxis()->SetLabelOffset(0.007);
   Graph_ha1p_pomsdbid_pdf3006->GetXaxis()->SetLabelSize(0.05);
   Graph_ha1p_pomsdbid_pdf3006->GetXaxis()->SetTitleSize(0.06);
   Graph_ha1p_pomsdbid_pdf3006->GetXaxis()->SetTitleOffset(0.9);
   Graph_ha1p_pomsdbid_pdf3006->GetXaxis()->SetTitleFont(42);
   Graph_ha1p_pomsdbid_pdf3006->GetYaxis()->SetTitle("N (+#eta_{cm}) / N (-#eta_{cm})");
   Graph_ha1p_pomsdbid_pdf3006->GetYaxis()->SetLabelFont(42);
   Graph_ha1p_pomsdbid_pdf3006->GetYaxis()->SetLabelOffset(0.007);
   Graph_ha1p_pomsdbid_pdf3006->GetYaxis()->SetLabelSize(0.05);
   Graph_ha1p_pomsdbid_pdf3006->GetYaxis()->SetTitleSize(0.06);
   Graph_ha1p_pomsdbid_pdf3006->GetYaxis()->SetTitleOffset(1.25);
   Graph_ha1p_pomsdbid_pdf3006->GetYaxis()->SetTitleFont(42);
   Graph_ha1p_pomsdbid_pdf3006->GetZaxis()->SetLabelFont(42);
   Graph_ha1p_pomsdbid_pdf3006->GetZaxis()->SetLabelOffset(0.007);
   Graph_ha1p_pomsdbid_pdf3006->GetZaxis()->SetLabelSize(0.05);
   Graph_ha1p_pomsdbid_pdf3006->GetZaxis()->SetTitleSize(0.06);
   Graph_ha1p_pomsdbid_pdf3006->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_ha1p_pomsdbid_pdf3006);
   
   grae->Draw("a5");
   
   TLegend *leg = new TLegend(0.55,0.68,0.88,0.91,NULL,"brNDC");
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
   entry=leg->AddEntry("ha1p_pomsdbid_pdf","nCTEQ15_208_82","lpf");

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
   A3->Modified();
   A3->cd();
   A3->SetSelected(A3);
}

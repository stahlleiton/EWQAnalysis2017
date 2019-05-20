void CT14_ForwardBackward_Ratio_Inc()
{
//=========Macro generated from canvas: c/
//=========  (Tue May 29 12:18:46 2018) by ROOT version6.06/00
   TCanvas *c = new TCanvas("c", "",0,0,600,600);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c->Range(0,0,1,1);
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
   
   Double_t hForwardBackwardRatio_Inc_wiecfvkv_pdf_fx3004[10] = {
   0.1,
   0.3,
   0.5,
   0.7,
   0.9,
   1.1,
   1.3,
   1.5,
   1.7,
   1.865};
   Double_t hForwardBackwardRatio_Inc_wiecfvkv_pdf_fy3004[10] = {
   1.003327,
   1.001404,
   0.9987258,
   1.000469,
   1.008213,
   1.008673,
   1.008204,
   1.020619,
   1.02406,
   1.024763};
   Double_t hForwardBackwardRatio_Inc_wiecfvkv_pdf_felx3004[10] = {
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.065};
   Double_t hForwardBackwardRatio_Inc_wiecfvkv_pdf_fely3004[10] = {
   0.006430418,
   0.01008219,
   0.00132236,
   0.004568841,
   0.00838858,
   0.002609202,
   0.001423,
   0.01027559,
   0.004574841,
   0};
   Double_t hForwardBackwardRatio_Inc_wiecfvkv_pdf_fehx3004[10] = {
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.065};
   Double_t hForwardBackwardRatio_Inc_wiecfvkv_pdf_fehy3004[10] = {
   0.002708477,
   0.003371278,
   0.007693347,
   0.004865333,
   0.001966853,
   0.006169183,
   0.006884357,
   0.002159401,
   0.003340935,
   0.02617781};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(10,hForwardBackwardRatio_Inc_wiecfvkv_pdf_fx3004,hForwardBackwardRatio_Inc_wiecfvkv_pdf_fy3004,hForwardBackwardRatio_Inc_wiecfvkv_pdf_felx3004,hForwardBackwardRatio_Inc_wiecfvkv_pdf_fehx3004,hForwardBackwardRatio_Inc_wiecfvkv_pdf_fely3004,hForwardBackwardRatio_Inc_wiecfvkv_pdf_fehy3004);
   grae->SetName("hForwardBackwardRatio_Inc_wiecfvkv_pdf");
   grae->SetTitle("");
   grae->SetFillColor(1);
   grae->Draw("alp");
   c->Modified();
   c->cd();
   c->SetSelected(c);
}

void nCTEQ15_ForwardBackward_Ratio_Inc()
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
   
   Double_t hForwardBackwardRatio_Inc_afpituuh_pdf_fx3016[10] = {
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
   Double_t hForwardBackwardRatio_Inc_afpituuh_pdf_fy3016[10] = {
   0.9901973,
   0.9826617,
   0.9651726,
   0.954473,
   0.9463504,
   0.9387919,
   0.9264547,
   0.9192091,
   0.9125227,
   0.9088137};
   Double_t hForwardBackwardRatio_Inc_afpituuh_pdf_felx3016[10] = {
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
   Double_t hForwardBackwardRatio_Inc_afpituuh_pdf_fely3016[10] = {
   0.003673593,
   0.0193474,
   0.02223016,
   0.03124732,
   0.04199121,
   0.04916281,
   0.0484406,
   0.05527663,
   0.05516189,
   0.0560003};
   Double_t hForwardBackwardRatio_Inc_afpituuh_pdf_fehx3016[10] = {
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
   Double_t hForwardBackwardRatio_Inc_afpituuh_pdf_fehy3016[10] = {
   0.01349868,
   0.02191956,
   0.03997965,
   0.04866539,
   0.06016177,
   0.06878158,
   0.07577826,
   0.07711235,
   0.08342245,
   0.07958377};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(10,hForwardBackwardRatio_Inc_afpituuh_pdf_fx3016,hForwardBackwardRatio_Inc_afpituuh_pdf_fy3016,hForwardBackwardRatio_Inc_afpituuh_pdf_felx3016,hForwardBackwardRatio_Inc_afpituuh_pdf_fehx3016,hForwardBackwardRatio_Inc_afpituuh_pdf_fely3016,hForwardBackwardRatio_Inc_afpituuh_pdf_fehy3016);
   grae->SetName("hForwardBackwardRatio_Inc_afpituuh_pdf");
   grae->SetTitle("");
   grae->SetFillColor(1);
   grae->Draw("alp");
   c->Modified();
   c->cd();
   c->SetSelected(c);
}

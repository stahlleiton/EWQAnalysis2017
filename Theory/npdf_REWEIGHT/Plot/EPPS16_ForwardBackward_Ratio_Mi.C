void EPPS16_ForwardBackward_Ratio_Mi()
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
   
   Double_t hForwardBackwardRatio_yrhgsqxi_pdf_fx3011[10] = {
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
   Double_t hForwardBackwardRatio_yrhgsqxi_pdf_fy3011[10] = {
   0.9864246,
   0.959353,
   0.924716,
   0.8991607,
   0.8709249,
   0.8447637,
   0.8179572,
   0.7896617,
   0.7649657,
   0.738382};
   Double_t hForwardBackwardRatio_yrhgsqxi_pdf_felx3011[10] = {
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
   Double_t hForwardBackwardRatio_yrhgsqxi_pdf_fely3011[10] = {
   0.004472798,
   0.009542538,
   0.008752259,
   0.01462903,
   0.01108498,
   0.01742957,
   0.02438653,
   0.02309526,
   0.02371212,
   0.02106005};
   Double_t hForwardBackwardRatio_yrhgsqxi_pdf_fehx3011[10] = {
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
   Double_t hForwardBackwardRatio_yrhgsqxi_pdf_fehy3011[10] = {
   0.006883433,
   0.004467838,
   0.01127895,
   0.008995091,
   0.01707501,
   0.01907378,
   0.0101898,
   0.01789254,
   0.02635515,
   0.02223025};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(10,hForwardBackwardRatio_yrhgsqxi_pdf_fx3011,hForwardBackwardRatio_yrhgsqxi_pdf_fy3011,hForwardBackwardRatio_yrhgsqxi_pdf_felx3011,hForwardBackwardRatio_yrhgsqxi_pdf_fehx3011,hForwardBackwardRatio_yrhgsqxi_pdf_fely3011,hForwardBackwardRatio_yrhgsqxi_pdf_fehy3011);
   grae->SetName("hForwardBackwardRatio_yrhgsqxi_pdf");
   grae->SetTitle("");
   grae->SetFillColor(1);
   grae->Draw("alp");
   c->Modified();
   c->cd();
   c->SetSelected(c);
}

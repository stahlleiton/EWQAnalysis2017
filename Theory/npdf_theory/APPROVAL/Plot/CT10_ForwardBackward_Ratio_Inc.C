void CT10_ForwardBackward_Ratio_Inc()
{
//=========Macro generated from canvas: c/
//=========  (Fri May 11 16:18:31 2018) by ROOT version6.06/00
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
   
   Double_t hForwardBackwardRatio_Inc_qdinrwlw_pdf_fx3004[10] = {
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
   Double_t hForwardBackwardRatio_Inc_qdinrwlw_pdf_fy3004[10] = {
   1.003443,
   1.003671,
   1.006122,
   1.009587,
   1.017421,
   1.020384,
   1.020993,
   1.02925,
   1.03341,
   1.03847};
   Double_t hForwardBackwardRatio_Inc_qdinrwlw_pdf_felx3004[10] = {
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
   Double_t hForwardBackwardRatio_Inc_qdinrwlw_pdf_fely3004[10] = {
   0.003523965,
   0.003247276,
   0.002005988,
   0.002894526,
   0.003813249,
   0.005022102,
   0.001663841,
   0.006181681,
   0.007510675,
   0.001006151};
   Double_t hForwardBackwardRatio_Inc_qdinrwlw_pdf_fehx3004[10] = {
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
   Double_t hForwardBackwardRatio_Inc_qdinrwlw_pdf_fehy3004[10] = {
   0.004297403,
   0.001987226,
   0.004371281,
   0.003515773,
   0.002660478,
   0.002486827,
   0.007700758,
   0.00259338,
   0.003134803,
   0.009879413};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(10,hForwardBackwardRatio_Inc_qdinrwlw_pdf_fx3004,hForwardBackwardRatio_Inc_qdinrwlw_pdf_fy3004,hForwardBackwardRatio_Inc_qdinrwlw_pdf_felx3004,hForwardBackwardRatio_Inc_qdinrwlw_pdf_fehx3004,hForwardBackwardRatio_Inc_qdinrwlw_pdf_fely3004,hForwardBackwardRatio_Inc_qdinrwlw_pdf_fehy3004);
   grae->SetName("hForwardBackwardRatio_Inc_qdinrwlw_pdf");
   grae->SetTitle("");
   grae->SetFillColor(1);
   grae->Draw("alp");
   c->Modified();
   c->cd();
   c->SetSelected(c);
}

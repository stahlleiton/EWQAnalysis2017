void CT14_ForwardBackward_Ratio_Inc()
{
//=========Macro generated from canvas: c/
//=========  (Thu Mar 28 20:02:04 2019) by ROOT version 6.12/07
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
   
   Double_t hForwardBackwardRatio_Inc_dynqfrwm_pdf_fx3010[10] = {
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
   Double_t hForwardBackwardRatio_Inc_dynqfrwm_pdf_fy3010[10] = {
   1.007122,
   1.004488,
   1.00971,
   1.011673,
   1.020019,
   1.023663,
   1.026357,
   1.033182,
   1.035446,
   1.044762};
   Double_t hForwardBackwardRatio_Inc_dynqfrwm_pdf_felx3010[10] = {
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
   Double_t hForwardBackwardRatio_Inc_dynqfrwm_pdf_fely3010[10] = {
   0.01254379,
   0.001737757,
   0.008338488,
   0.002563678,
   0.004222031,
   0.006077461,
   0.007092925,
   0.006009441,
   0.002233938,
   0.008751382};
   Double_t hForwardBackwardRatio_Inc_dynqfrwm_pdf_fehx3010[10] = {
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
   Double_t hForwardBackwardRatio_Inc_dynqfrwm_pdf_fehy3010[10] = {
   0,
   0.006041244,
   0.001557884,
   0.008675031,
   0.007765526,
   0.004396419,
   0.004549892,
   0.004607572,
   0.008266665,
   0.003308993};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(10,hForwardBackwardRatio_Inc_dynqfrwm_pdf_fx3010,hForwardBackwardRatio_Inc_dynqfrwm_pdf_fy3010,hForwardBackwardRatio_Inc_dynqfrwm_pdf_felx3010,hForwardBackwardRatio_Inc_dynqfrwm_pdf_fehx3010,hForwardBackwardRatio_Inc_dynqfrwm_pdf_fely3010,hForwardBackwardRatio_Inc_dynqfrwm_pdf_fehy3010);
   grae->SetName("hForwardBackwardRatio_Inc_dynqfrwm_pdf");
   grae->SetTitle("");
   grae->SetFillStyle(1000);
   grae->Draw("alp");
   c->Modified();
   c->cd();
   c->SetSelected(c);
}

void EPS09_ForwardBackward_Ratio_Pl()
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
   
   Double_t hForwardBackwardRatio_ucrrxogi_pdf_fx3024[10] = {
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
   Double_t hForwardBackwardRatio_ucrrxogi_pdf_fy3024[10] = {
   0.992147,
   0.9869962,
   0.9806539,
   0.975882,
   0.976548,
   0.9831443,
   0.9895453,
   0.9963073,
   1.012973,
   1.036189};
   Double_t hForwardBackwardRatio_ucrrxogi_pdf_felx3024[10] = {
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
   Double_t hForwardBackwardRatio_ucrrxogi_pdf_fely3024[10] = {
   0.009440419,
   0.01202409,
   0.0199038,
   0.02520132,
   0.0342673,
   0.03755428,
   0.0467946,
   0.04534121,
   0.0506278,
   0.05311413};
   Double_t hForwardBackwardRatio_ucrrxogi_pdf_fehx3024[10] = {
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
   Double_t hForwardBackwardRatio_ucrrxogi_pdf_fehy3024[10] = {
   0.005126722,
   0.0170737,
   0.01849931,
   0.02753288,
   0.03461488,
   0.03761363,
   0.03786639,
   0.04557599,
   0.0450975,
   0.0472572};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(10,hForwardBackwardRatio_ucrrxogi_pdf_fx3024,hForwardBackwardRatio_ucrrxogi_pdf_fy3024,hForwardBackwardRatio_ucrrxogi_pdf_felx3024,hForwardBackwardRatio_ucrrxogi_pdf_fehx3024,hForwardBackwardRatio_ucrrxogi_pdf_fely3024,hForwardBackwardRatio_ucrrxogi_pdf_fehy3024);
   grae->SetName("hForwardBackwardRatio_ucrrxogi_pdf");
   grae->SetTitle("");
   grae->SetFillColor(1);
   grae->Draw("alp");
   c->Modified();
   c->cd();
   c->SetSelected(c);
}

void nCTEQ15_ForwardBackward_Ratio_Inc()
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
   
   Double_t hForwardBackwardRatio_Inc_kmhujhto_pdf_fx3034[10] = {
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
   Double_t hForwardBackwardRatio_Inc_kmhujhto_pdf_fy3034[10] = {
   0.9845684,
   0.9438345,
   0.9123552,
   0.8820625,
   0.8556178,
   0.8334935,
   0.8130121,
   0.797371,
   0.7845626,
   0.7768694};
   Double_t hForwardBackwardRatio_Inc_kmhujhto_pdf_felx3034[10] = {
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
   Double_t hForwardBackwardRatio_Inc_kmhujhto_pdf_fely3034[10] = {
   0.00796442,
   0.01749709,
   0.02659184,
   0.03612525,
   0.04562042,
   0.05333753,
   0.06095085,
   0.06674404,
   0.06986845,
   0.07865652};
   Double_t hForwardBackwardRatio_Inc_kmhujhto_pdf_fehx3034[10] = {
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
   Double_t hForwardBackwardRatio_Inc_kmhujhto_pdf_fehy3034[10] = {
   0.005997128,
   0.02148806,
   0.02737058,
   0.03994131,
   0.04777681,
   0.05580508,
   0.05680441,
   0.06158351,
   0.06410674,
   0.06789915};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(10,hForwardBackwardRatio_Inc_kmhujhto_pdf_fx3034,hForwardBackwardRatio_Inc_kmhujhto_pdf_fy3034,hForwardBackwardRatio_Inc_kmhujhto_pdf_felx3034,hForwardBackwardRatio_Inc_kmhujhto_pdf_fehx3034,hForwardBackwardRatio_Inc_kmhujhto_pdf_fely3034,hForwardBackwardRatio_Inc_kmhujhto_pdf_fehy3034);
   grae->SetName("hForwardBackwardRatio_Inc_kmhujhto_pdf");
   grae->SetTitle("");
   grae->SetFillColor(1);
   grae->Draw("alp");
   c->Modified();
   c->cd();
   c->SetSelected(c);
}

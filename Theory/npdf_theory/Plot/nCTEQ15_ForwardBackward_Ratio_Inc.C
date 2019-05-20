void nCTEQ15_ForwardBackward_Ratio_Inc()
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
   
   Double_t hForwardBackwardRatio_Inc_sdpkzqhw_pdf_fx3034[10] = {
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
   Double_t hForwardBackwardRatio_Inc_sdpkzqhw_pdf_fy3034[10] = {
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
   Double_t hForwardBackwardRatio_Inc_sdpkzqhw_pdf_felx3034[10] = {
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
   Double_t hForwardBackwardRatio_Inc_sdpkzqhw_pdf_fely3034[10] = {
   0.009844003,
   0.01828313,
   0.02832079,
   0.03808075,
   0.04779101,
   0.05521187,
   0.06331587,
   0.06924665,
   0.07303691,
   0.08167949};
   Double_t hForwardBackwardRatio_Inc_sdpkzqhw_pdf_fehx3034[10] = {
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
   Double_t hForwardBackwardRatio_Inc_sdpkzqhw_pdf_fehy3034[10] = {
   0.006368638,
   0.02274531,
   0.0280755,
   0.04156621,
   0.04933252,
   0.05831828,
   0.05870321,
   0.06395784,
   0.06685638,
   0.07051238};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(10,hForwardBackwardRatio_Inc_sdpkzqhw_pdf_fx3034,hForwardBackwardRatio_Inc_sdpkzqhw_pdf_fy3034,hForwardBackwardRatio_Inc_sdpkzqhw_pdf_felx3034,hForwardBackwardRatio_Inc_sdpkzqhw_pdf_fehx3034,hForwardBackwardRatio_Inc_sdpkzqhw_pdf_fely3034,hForwardBackwardRatio_Inc_sdpkzqhw_pdf_fehy3034);
   grae->SetName("hForwardBackwardRatio_Inc_sdpkzqhw_pdf");
   grae->SetTitle("");
   grae->SetFillStyle(1000);
   grae->Draw("alp");
   c->Modified();
   c->cd();
   c->SetSelected(c);
}

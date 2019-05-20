void EPPS16_ForwardBackward_Ratio_Mi()
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
   
   Double_t hForwardBackwardRatio_dmewoimi_pdf_fx3017[10] = {
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
   Double_t hForwardBackwardRatio_dmewoimi_pdf_fy3017[10] = {
   0.9872179,
   0.9569585,
   0.9199166,
   0.890173,
   0.86306,
   0.8320614,
   0.7990534,
   0.7710296,
   0.7455351,
   0.7184044};
   Double_t hForwardBackwardRatio_dmewoimi_pdf_felx3017[10] = {
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
   Double_t hForwardBackwardRatio_dmewoimi_pdf_fely3017[10] = {
   0.01116718,
   0.02895549,
   0.04028672,
   0.04482742,
   0.06167399,
   0.06864932,
   0.07192289,
   0.07426386,
   0.08457957,
   0.0748669};
   Double_t hForwardBackwardRatio_dmewoimi_pdf_fehx3017[10] = {
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
   Double_t hForwardBackwardRatio_dmewoimi_pdf_fehy3017[10] = {
   0.00628199,
   0.01880561,
   0.0367885,
   0.04727695,
   0.05074465,
   0.06486515,
   0.06916108,
   0.06889333,
   0.07174585,
   0.07237423};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(10,hForwardBackwardRatio_dmewoimi_pdf_fx3017,hForwardBackwardRatio_dmewoimi_pdf_fy3017,hForwardBackwardRatio_dmewoimi_pdf_felx3017,hForwardBackwardRatio_dmewoimi_pdf_fehx3017,hForwardBackwardRatio_dmewoimi_pdf_fely3017,hForwardBackwardRatio_dmewoimi_pdf_fehy3017);
   grae->SetName("hForwardBackwardRatio_dmewoimi_pdf");
   grae->SetTitle("");
   grae->SetFillColor(1);
   grae->Draw("alp");
   c->Modified();
   c->cd();
   c->SetSelected(c);
}

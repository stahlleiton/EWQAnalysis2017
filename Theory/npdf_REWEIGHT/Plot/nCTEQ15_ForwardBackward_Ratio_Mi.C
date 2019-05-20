void nCTEQ15_ForwardBackward_Ratio_Mi()
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
   
   Double_t hForwardBackwardRatio_dxgxgqxz_pdf_fx3017[10] = {
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
   Double_t hForwardBackwardRatio_dxgxgqxz_pdf_fy3017[10] = {
   0.9870969,
   0.9639618,
   0.9340343,
   0.9086681,
   0.8846139,
   0.8607687,
   0.8298706,
   0.8048885,
   0.7788555,
   0.7541152};
   Double_t hForwardBackwardRatio_dxgxgqxz_pdf_felx3017[10] = {
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
   Double_t hForwardBackwardRatio_dxgxgqxz_pdf_fely3017[10] = {
   0.003795548,
   0.01602055,
   0.0255984,
   0.03121961,
   0.03734758,
   0.04980297,
   0.04435054,
   0.05059588,
   0.05077331,
   0.05355548};
   Double_t hForwardBackwardRatio_dxgxgqxz_pdf_fehx3017[10] = {
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
   Double_t hForwardBackwardRatio_dxgxgqxz_pdf_fehy3017[10] = {
   0.01635067,
   0.02277128,
   0.0382728,
   0.04724899,
   0.05640717,
   0.0650639,
   0.0673348,
   0.06807763,
   0.07115218,
   0.06468352};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(10,hForwardBackwardRatio_dxgxgqxz_pdf_fx3017,hForwardBackwardRatio_dxgxgqxz_pdf_fy3017,hForwardBackwardRatio_dxgxgqxz_pdf_felx3017,hForwardBackwardRatio_dxgxgqxz_pdf_fehx3017,hForwardBackwardRatio_dxgxgqxz_pdf_fely3017,hForwardBackwardRatio_dxgxgqxz_pdf_fehy3017);
   grae->SetName("hForwardBackwardRatio_dxgxgqxz_pdf");
   grae->SetTitle("");
   grae->SetFillColor(1);
   grae->Draw("alp");
   c->Modified();
   c->cd();
   c->SetSelected(c);
}

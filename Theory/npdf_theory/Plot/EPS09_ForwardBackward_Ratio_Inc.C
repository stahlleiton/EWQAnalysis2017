void EPS09_ForwardBackward_Ratio_Inc()
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
   
   Double_t hForwardBackwardRatio_Inc_jgberysm_pdf_fx3022[10] = {
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
   Double_t hForwardBackwardRatio_Inc_jgberysm_pdf_fy3022[10] = {
   0.9883355,
   0.9713978,
   0.9526381,
   0.9365803,
   0.9239611,
   0.9145011,
   0.9047774,
   0.8960387,
   0.894389,
   0.8945876};
   Double_t hForwardBackwardRatio_Inc_jgberysm_pdf_felx3022[10] = {
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
   Double_t hForwardBackwardRatio_Inc_jgberysm_pdf_fely3022[10] = {
   0.004637817,
   0.009679135,
   0.01892279,
   0.02363861,
   0.02971647,
   0.03357101,
   0.03958737,
   0.03857385,
   0.04348525,
   0.04396508};
   Double_t hForwardBackwardRatio_Inc_jgberysm_pdf_fehx3022[10] = {
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
   Double_t hForwardBackwardRatio_Inc_jgberysm_pdf_fehy3022[10] = {
   0.006427302,
   0.0162868,
   0.01800076,
   0.02401029,
   0.03022497,
   0.03462988,
   0.03358292,
   0.03848627,
   0.03854744,
   0.03906115};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(10,hForwardBackwardRatio_Inc_jgberysm_pdf_fx3022,hForwardBackwardRatio_Inc_jgberysm_pdf_fy3022,hForwardBackwardRatio_Inc_jgberysm_pdf_felx3022,hForwardBackwardRatio_Inc_jgberysm_pdf_fehx3022,hForwardBackwardRatio_Inc_jgberysm_pdf_fely3022,hForwardBackwardRatio_Inc_jgberysm_pdf_fehy3022);
   grae->SetName("hForwardBackwardRatio_Inc_jgberysm_pdf");
   grae->SetTitle("");
   grae->SetFillStyle(1000);
   grae->Draw("alp");
   c->Modified();
   c->cd();
   c->SetSelected(c);
}

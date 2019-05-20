void EPPS16_ForwardBackward_Ratio_Inc()
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
   
   Double_t hForwardBackwardRatio_Inc_wheryblw_pdf_fx3010[10] = {
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
   Double_t hForwardBackwardRatio_Inc_wheryblw_pdf_fy3010[10] = {
   0.9889374,
   0.9754338,
   0.955188,
   0.9422884,
   0.9301074,
   0.9192874,
   0.9122632,
   0.9023745,
   0.9000344,
   0.9003616};
   Double_t hForwardBackwardRatio_Inc_wheryblw_pdf_felx3010[10] = {
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
   Double_t hForwardBackwardRatio_Inc_wheryblw_pdf_fely3010[10] = {
   0.003124578,
   0.007628103,
   0.009448523,
   0.01301366,
   0.01535325,
   0.01499199,
   0.02256035,
   0.01953106,
   0.02168845,
   0.0244586};
   Double_t hForwardBackwardRatio_Inc_wheryblw_pdf_fehx3010[10] = {
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
   Double_t hForwardBackwardRatio_Inc_wheryblw_pdf_fehy3010[10] = {
   0.005964894,
   0.004664016,
   0.009519514,
   0.008700403,
   0.01321926,
   0.01805065,
   0.01113375,
   0.01773265,
   0.01940334,
   0.01788538};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(10,hForwardBackwardRatio_Inc_wheryblw_pdf_fx3010,hForwardBackwardRatio_Inc_wheryblw_pdf_fy3010,hForwardBackwardRatio_Inc_wheryblw_pdf_felx3010,hForwardBackwardRatio_Inc_wheryblw_pdf_fehx3010,hForwardBackwardRatio_Inc_wheryblw_pdf_fely3010,hForwardBackwardRatio_Inc_wheryblw_pdf_fehy3010);
   grae->SetName("hForwardBackwardRatio_Inc_wheryblw_pdf");
   grae->SetTitle("");
   grae->SetFillColor(1);
   grae->Draw("alp");
   c->Modified();
   c->cd();
   c->SetSelected(c);
}

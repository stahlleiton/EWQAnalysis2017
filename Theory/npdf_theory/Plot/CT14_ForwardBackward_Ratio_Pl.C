void CT14_ForwardBackward_Ratio_Pl()
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
   
   Double_t hForwardBackwardRatio_oizhqykx_pdf_fx3012[10] = {
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
   Double_t hForwardBackwardRatio_oizhqykx_pdf_fy3012[10] = {
   1.015028,
   1.022036,
   1.041326,
   1.062626,
   1.087,
   1.107904,
   1.131902,
   1.164073,
   1.185143,
   1.222182};
   Double_t hForwardBackwardRatio_oizhqykx_pdf_felx3012[10] = {
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
   Double_t hForwardBackwardRatio_oizhqykx_pdf_fely3012[10] = {
   0.01782196,
   0.004422112,
   0.008649755,
   0.009327171,
   0.006877089,
   0.01492304,
   0.01640024,
   0.02117836,
   0.01144065,
   0.03011314};
   Double_t hForwardBackwardRatio_oizhqykx_pdf_fehx3012[10] = {
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
   Double_t hForwardBackwardRatio_oizhqykx_pdf_fehy3012[10] = {
   0,
   0.007801223,
   0.003944037,
   0.01127936,
   0.0173735,
   0.01251444,
   0.01611152,
   0.01667588,
   0.02880735,
   0.01728268};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(10,hForwardBackwardRatio_oizhqykx_pdf_fx3012,hForwardBackwardRatio_oizhqykx_pdf_fy3012,hForwardBackwardRatio_oizhqykx_pdf_felx3012,hForwardBackwardRatio_oizhqykx_pdf_fehx3012,hForwardBackwardRatio_oizhqykx_pdf_fely3012,hForwardBackwardRatio_oizhqykx_pdf_fehy3012);
   grae->SetName("hForwardBackwardRatio_oizhqykx_pdf");
   grae->SetTitle("");
   grae->SetFillStyle(1000);
   grae->Draw("alp");
   c->Modified();
   c->cd();
   c->SetSelected(c);
}

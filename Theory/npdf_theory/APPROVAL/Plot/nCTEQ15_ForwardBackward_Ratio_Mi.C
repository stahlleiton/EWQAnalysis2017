void nCTEQ15_ForwardBackward_Ratio_Mi()
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
   
   Double_t hForwardBackwardRatio_nkgzeiec_pdf_fx3035[10] = {
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
   Double_t hForwardBackwardRatio_nkgzeiec_pdf_fy3035[10] = {
   0.977294,
   0.9273855,
   0.8796414,
   0.8417896,
   0.8046129,
   0.7638041,
   0.7321435,
   0.701543,
   0.6712132,
   0.6514753};
   Double_t hForwardBackwardRatio_nkgzeiec_pdf_felx3035[10] = {
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
   Double_t hForwardBackwardRatio_nkgzeiec_pdf_fely3035[10] = {
   0.005112137,
   0.02007458,
   0.02464102,
   0.03454406,
   0.04712218,
   0.04755358,
   0.05904302,
   0.05681418,
   0.06107747,
   0.06791261};
   Double_t hForwardBackwardRatio_nkgzeiec_pdf_fehx3035[10] = {
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
   Double_t hForwardBackwardRatio_nkgzeiec_pdf_fehy3035[10] = {
   0.01024153,
   0.02214663,
   0.02802245,
   0.03728843,
   0.04576887,
   0.05272866,
   0.0497084,
   0.05617803,
   0.05599345,
   0.05581494};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(10,hForwardBackwardRatio_nkgzeiec_pdf_fx3035,hForwardBackwardRatio_nkgzeiec_pdf_fy3035,hForwardBackwardRatio_nkgzeiec_pdf_felx3035,hForwardBackwardRatio_nkgzeiec_pdf_fehx3035,hForwardBackwardRatio_nkgzeiec_pdf_fely3035,hForwardBackwardRatio_nkgzeiec_pdf_fehy3035);
   grae->SetName("hForwardBackwardRatio_nkgzeiec_pdf");
   grae->SetTitle("");
   grae->SetFillColor(1);
   grae->Draw("alp");
   c->Modified();
   c->cd();
   c->SetSelected(c);
}

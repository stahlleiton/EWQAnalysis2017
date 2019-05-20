void EPPS16_ForwardBackward_Ratio_Pl()
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
   
   Double_t hForwardBackwardRatio_nycjvurh_pdf_fx3018[10] = {
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
   Double_t hForwardBackwardRatio_nycjvurh_pdf_fy3018[10] = {
   0.9914179,
   0.9868066,
   0.9787798,
   0.9770872,
   0.9762217,
   0.9775951,
   0.9839736,
   0.9915134,
   1.001143,
   1.02921};
   Double_t hForwardBackwardRatio_nycjvurh_pdf_felx3018[10] = {
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
   Double_t hForwardBackwardRatio_nycjvurh_pdf_fely3018[10] = {
   0.006350575,
   0.0271097,
   0.03719965,
   0.06193689,
   0.07689448,
   0.07519107,
   0.09971078,
   0.1081256,
   0.1040074,
   0.1297722};
   Double_t hForwardBackwardRatio_nycjvurh_pdf_fehx3018[10] = {
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
   Double_t hForwardBackwardRatio_nycjvurh_pdf_fehy3018[10] = {
   0.009271281,
   0.02938228,
   0.04716056,
   0.05039133,
   0.07290863,
   0.08669505,
   0.09262917,
   0.09841089,
   0.1188775,
   0.1081478};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(10,hForwardBackwardRatio_nycjvurh_pdf_fx3018,hForwardBackwardRatio_nycjvurh_pdf_fy3018,hForwardBackwardRatio_nycjvurh_pdf_felx3018,hForwardBackwardRatio_nycjvurh_pdf_fehx3018,hForwardBackwardRatio_nycjvurh_pdf_fely3018,hForwardBackwardRatio_nycjvurh_pdf_fehy3018);
   grae->SetName("hForwardBackwardRatio_nycjvurh_pdf");
   grae->SetTitle("");
   grae->SetFillColor(1);
   grae->Draw("alp");
   c->Modified();
   c->cd();
   c->SetSelected(c);
}

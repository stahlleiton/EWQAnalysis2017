void EPS09_CT14_ForwardBackward_Ratio_Pl()
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
   
   Double_t hForwardBackwardRatio_aszfemqo_pdf_fx3030[10] = {
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
   Double_t hForwardBackwardRatio_aszfemqo_pdf_fy3030[10] = {
   0.9912542,
   0.9919463,
   0.9823563,
   0.9804887,
   0.9798509,
   0.9889419,
   0.9939959,
   1.00211,
   1.022002,
   1.04671};
   Double_t hForwardBackwardRatio_aszfemqo_pdf_felx3030[10] = {
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
   Double_t hForwardBackwardRatio_aszfemqo_pdf_fely3030[10] = {
   0.006079832,
   0.02172382,
   0.02158856,
   0.02914376,
   0.03281024,
   0.04234206,
   0.04584438,
   0.04461716,
   0.05290224,
   0.05643688};
   Double_t hForwardBackwardRatio_aszfemqo_pdf_fehx3030[10] = {
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
   Double_t hForwardBackwardRatio_aszfemqo_pdf_fehy3030[10] = {
   0.008107374,
   0.007804418,
   0.02193774,
   0.02301271,
   0.0387465,
   0.03363148,
   0.04237019,
   0.04942561,
   0.05007132,
   0.04981931};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(10,hForwardBackwardRatio_aszfemqo_pdf_fx3030,hForwardBackwardRatio_aszfemqo_pdf_fy3030,hForwardBackwardRatio_aszfemqo_pdf_felx3030,hForwardBackwardRatio_aszfemqo_pdf_fehx3030,hForwardBackwardRatio_aszfemqo_pdf_fely3030,hForwardBackwardRatio_aszfemqo_pdf_fehy3030);
   grae->SetName("hForwardBackwardRatio_aszfemqo_pdf");
   grae->SetTitle("");
   grae->SetFillStyle(1000);
   grae->Draw("alp");
   c->Modified();
   c->cd();
   c->SetSelected(c);
}

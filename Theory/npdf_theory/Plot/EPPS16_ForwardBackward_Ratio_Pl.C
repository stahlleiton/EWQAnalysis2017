void EPPS16_ForwardBackward_Ratio_Pl()
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
   
   Double_t hForwardBackwardRatio_unhauupb_pdf_fx3018[10] = {
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
   Double_t hForwardBackwardRatio_unhauupb_pdf_fy3018[10] = {
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
   Double_t hForwardBackwardRatio_unhauupb_pdf_felx3018[10] = {
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
   Double_t hForwardBackwardRatio_unhauupb_pdf_fely3018[10] = {
   0.007407233,
   0.02732907,
   0.03743752,
   0.06283646,
   0.0779527,
   0.07546897,
   0.100857,
   0.1095102,
   0.1043695,
   0.1323223};
   Double_t hForwardBackwardRatio_unhauupb_pdf_fehx3018[10] = {
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
   Double_t hForwardBackwardRatio_unhauupb_pdf_fehy3018[10] = {
   0.01066103,
   0.03111952,
   0.04812549,
   0.05079643,
   0.07427136,
   0.08846538,
   0.09356613,
   0.09923948,
   0.1213726,
   0.1090706};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(10,hForwardBackwardRatio_unhauupb_pdf_fx3018,hForwardBackwardRatio_unhauupb_pdf_fy3018,hForwardBackwardRatio_unhauupb_pdf_felx3018,hForwardBackwardRatio_unhauupb_pdf_fehx3018,hForwardBackwardRatio_unhauupb_pdf_fely3018,hForwardBackwardRatio_unhauupb_pdf_fehy3018);
   grae->SetName("hForwardBackwardRatio_unhauupb_pdf");
   grae->SetTitle("");
   grae->SetFillStyle(1000);
   grae->Draw("alp");
   c->Modified();
   c->cd();
   c->SetSelected(c);
}

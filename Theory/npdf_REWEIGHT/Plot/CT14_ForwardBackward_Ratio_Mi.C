void CT14_ForwardBackward_Ratio_Mi()
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
   
   Double_t hForwardBackwardRatio_yukguuqo_pdf_fx3005[10] = {
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
   Double_t hForwardBackwardRatio_yukguuqo_pdf_fy3005[10] = {
   0.994778,
   0.9799527,
   0.9627461,
   0.9441959,
   0.92962,
   0.9140208,
   0.891353,
   0.877198,
   0.8579851,
   0.8360456};
   Double_t hForwardBackwardRatio_yukguuqo_pdf_felx3005[10] = {
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
   Double_t hForwardBackwardRatio_yukguuqo_pdf_fely3005[10] = {
   0.004088883,
   0.007833604,
   0.004816524,
   0.005768666,
   0.008335281,
   0.006111959,
   0.00492779,
   0.01079028,
   0.005479252,
   0};
   Double_t hForwardBackwardRatio_yukguuqo_pdf_fehx3005[10] = {
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
   Double_t hForwardBackwardRatio_yukguuqo_pdf_fehy3005[10] = {
   0.00464276,
   0.001026529,
   0.005628894,
   0.006276887,
   0.004455541,
   0.006227309,
   0.008059448,
   0.005485804,
   0.005904425,
   0.03371735};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(10,hForwardBackwardRatio_yukguuqo_pdf_fx3005,hForwardBackwardRatio_yukguuqo_pdf_fy3005,hForwardBackwardRatio_yukguuqo_pdf_felx3005,hForwardBackwardRatio_yukguuqo_pdf_fehx3005,hForwardBackwardRatio_yukguuqo_pdf_fely3005,hForwardBackwardRatio_yukguuqo_pdf_fehy3005);
   grae->SetName("hForwardBackwardRatio_yukguuqo_pdf");
   grae->SetTitle("");
   grae->SetFillColor(1);
   grae->Draw("alp");
   c->Modified();
   c->cd();
   c->SetSelected(c);
}

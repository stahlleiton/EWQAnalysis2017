void nCTEQ15_ForwardBackward_Ratio_Mi()
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
   
   Double_t hForwardBackwardRatio_svljlewy_pdf_fx3035[10] = {
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
   Double_t hForwardBackwardRatio_svljlewy_pdf_fy3035[10] = {
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
   Double_t hForwardBackwardRatio_svljlewy_pdf_felx3035[10] = {
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
   Double_t hForwardBackwardRatio_svljlewy_pdf_fely3035[10] = {
   0.006133503,
   0.02114631,
   0.02627601,
   0.03617737,
   0.05052946,
   0.04904786,
   0.0628552,
   0.05925915,
   0.06365228,
   0.07249567};
   Double_t hForwardBackwardRatio_svljlewy_pdf_fehx3035[10] = {
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
   Double_t hForwardBackwardRatio_svljlewy_pdf_fehy3035[10] = {
   0.01186996,
   0.02395573,
   0.02869467,
   0.03931907,
   0.04673181,
   0.0560989,
   0.0509045,
   0.05898911,
   0.05926756,
   0.05796841};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(10,hForwardBackwardRatio_svljlewy_pdf_fx3035,hForwardBackwardRatio_svljlewy_pdf_fy3035,hForwardBackwardRatio_svljlewy_pdf_felx3035,hForwardBackwardRatio_svljlewy_pdf_fehx3035,hForwardBackwardRatio_svljlewy_pdf_fely3035,hForwardBackwardRatio_svljlewy_pdf_fehy3035);
   grae->SetName("hForwardBackwardRatio_svljlewy_pdf");
   grae->SetTitle("");
   grae->SetFillStyle(1000);
   grae->Draw("alp");
   c->Modified();
   c->cd();
   c->SetSelected(c);
}

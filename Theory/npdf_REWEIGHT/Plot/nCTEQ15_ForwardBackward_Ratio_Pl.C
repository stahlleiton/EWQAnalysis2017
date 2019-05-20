void nCTEQ15_ForwardBackward_Ratio_Pl()
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
   
   Double_t hForwardBackwardRatio_bazgdkdq_pdf_fx3018[10] = {
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
   Double_t hForwardBackwardRatio_bazgdkdq_pdf_fy3018[10] = {
   0.9928014,
   0.9985462,
   0.9916614,
   0.9936767,
   0.9993475,
   1.005908,
   1.009299,
   1.017655,
   1.028742,
   1.044168};
   Double_t hForwardBackwardRatio_bazgdkdq_pdf_felx3018[10] = {
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
   Double_t hForwardBackwardRatio_bazgdkdq_pdf_fely3018[10] = {
   0.005154578,
   0.02370682,
   0.02082634,
   0.03253144,
   0.04833412,
   0.05078744,
   0.05393891,
   0.06219027,
   0.06266794,
   0.06153664};
   Double_t hForwardBackwardRatio_bazgdkdq_pdf_fehx3018[10] = {
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
   Double_t hForwardBackwardRatio_bazgdkdq_pdf_fehy3018[10] = {
   0.01340751,
   0.02185509,
   0.04209625,
   0.05048933,
   0.06385982,
   0.07340176,
   0.08462635,
   0.08609589,
   0.09557784,
   0.09613815};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(10,hForwardBackwardRatio_bazgdkdq_pdf_fx3018,hForwardBackwardRatio_bazgdkdq_pdf_fy3018,hForwardBackwardRatio_bazgdkdq_pdf_felx3018,hForwardBackwardRatio_bazgdkdq_pdf_fehx3018,hForwardBackwardRatio_bazgdkdq_pdf_fely3018,hForwardBackwardRatio_bazgdkdq_pdf_fehy3018);
   grae->SetName("hForwardBackwardRatio_bazgdkdq_pdf");
   grae->SetTitle("");
   grae->SetFillColor(1);
   grae->Draw("alp");
   c->Modified();
   c->cd();
   c->SetSelected(c);
}

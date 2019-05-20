void CT10_ForwardBackward_Ratio_Mi()
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
   
   Double_t hForwardBackwardRatio_yotzplcw_pdf_fx3005[10] = {
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
   Double_t hForwardBackwardRatio_yotzplcw_pdf_fy3005[10] = {
   0.9949959,
   0.9820575,
   0.9678833,
   0.9509708,
   0.9389985,
   0.9239438,
   0.9007124,
   0.8838656,
   0.8663865,
   0.8493971};
   Double_t hForwardBackwardRatio_yotzplcw_pdf_felx3005[10] = {
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
   Double_t hForwardBackwardRatio_yotzplcw_pdf_fely3005[10] = {
   0.002297881,
   0.001698065,
   0.004176787,
   0.004779116,
   0.007624872,
   0.006964364,
   0.006002745,
   0.008498983,
   0.01105284,
   0.007743882};
   Double_t hForwardBackwardRatio_yotzplcw_pdf_fehx3005[10] = {
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
   Double_t hForwardBackwardRatio_yotzplcw_pdf_fehy3005[10] = {
   0.004371081,
   0.004508948,
   0.006150061,
   0.006676645,
   0.002914558,
   0.006453269,
   0.01202591,
   0.008587717,
   0.008470457,
   0.01621984};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(10,hForwardBackwardRatio_yotzplcw_pdf_fx3005,hForwardBackwardRatio_yotzplcw_pdf_fy3005,hForwardBackwardRatio_yotzplcw_pdf_felx3005,hForwardBackwardRatio_yotzplcw_pdf_fehx3005,hForwardBackwardRatio_yotzplcw_pdf_fely3005,hForwardBackwardRatio_yotzplcw_pdf_fehy3005);
   grae->SetName("hForwardBackwardRatio_yotzplcw_pdf");
   grae->SetTitle("");
   grae->SetFillStyle(1000);
   grae->Draw("alp");
   c->Modified();
   c->cd();
   c->SetSelected(c);
}

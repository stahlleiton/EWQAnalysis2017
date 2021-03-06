void nCTEQ15_Cross_Section_Mi()
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
   
   Double_t hCrossSection_avngcllc_pdf_fx3032[24] = {
   -2.73,
   -2.5,
   -2.3,
   -2.065,
   -1.865,
   -1.7,
   -1.5,
   -1.3,
   -1.1,
   -0.9,
   -0.7,
   -0.5,
   -0.3,
   -0.1,
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
   Double_t hCrossSection_avngcllc_pdf_fy3032[24] = {
   108.8755,
   113.2587,
   115.9851,
   118.1488,
   119.6237,
   120.9316,
   121.5356,
   121.5387,
   121.4287,
   120.4911,
   119.2498,
   117.8252,
   115.3809,
   113.1418,
   110.5728,
   107.0026,
   103.644,
   100.3832,
   96.94868,
   92.74773,
   88.98379,
   85.26247,
   81.17089,
   77.93192};
   Double_t hCrossSection_avngcllc_pdf_felx3032[24] = {
   0.13,
   0.1,
   0.1,
   0.135,
   0.065,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
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
   Double_t hCrossSection_avngcllc_pdf_fely3032[24] = {
   2.341979,
   2.52162,
   2.376871,
   2.348348,
   1.89202,
   2.464501,
   2.526851,
   2.723527,
   4.110825,
   4.626629,
   5.628146,
   6.550445,
   6.74709,
   7.521067,
   7.761679,
   8.227557,
   8.134235,
   8.106723,
   8.254486,
   7.820482,
   8.026167,
   7.211555,
   7.148624,
   6.93642};
   Double_t hCrossSection_avngcllc_pdf_fehx3032[24] = {
   0.13,
   0.1,
   0.1,
   0.135,
   0.065,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
   0.1,
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
   Double_t hCrossSection_avngcllc_pdf_fehy3032[24] = {
   3.286419,
   2.925383,
   2.579169,
   2.345777,
   2.175018,
   0.9545738,
   0.9595166,
   1.351334,
   1.348877,
   2.329006,
   2.770139,
   3.353949,
   4.041926,
   4.221523,
   4.851196,
   5.593352,
   5.125644,
   5.649625,
   5.512534,
   5.888555,
   5.392422,
   5.60182,
   5.327667,
   5.535366};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(24,hCrossSection_avngcllc_pdf_fx3032,hCrossSection_avngcllc_pdf_fy3032,hCrossSection_avngcllc_pdf_felx3032,hCrossSection_avngcllc_pdf_fehx3032,hCrossSection_avngcllc_pdf_fely3032,hCrossSection_avngcllc_pdf_fehy3032);
   grae->SetName("hCrossSection_avngcllc_pdf");
   grae->SetTitle("");
   grae->SetFillColor(1);
   grae->Draw("alp");
   c->Modified();
   c->cd();
   c->SetSelected(c);
}

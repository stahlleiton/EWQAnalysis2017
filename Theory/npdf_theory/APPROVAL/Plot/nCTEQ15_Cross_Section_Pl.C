void nCTEQ15_Cross_Section_Pl()
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
   
   Double_t hCrossSection_ojpkuoit_pdf_fx3033[24] = {
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
   Double_t hCrossSection_ojpkuoit_pdf_fy3033[24] = {
   99.00308,
   114.303,
   123.5385,
   132.1217,
   136.599,
   138.9257,
   140.2984,
   140.8792,
   140.4074,
   139.8834,
   138.9021,
   136.6738,
   135.6048,
   132.8683,
   131.641,
   129.8864,
   128.5495,
   127.3229,
   125.8323,
   125.491,
   124.3651,
   123.5164,
   122.7035,
   121.1197};
   Double_t hCrossSection_ojpkuoit_pdf_felx3033[24] = {
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
   Double_t hCrossSection_ojpkuoit_pdf_fely3033[24] = {
   2.114805,
   2.952918,
   2.426383,
   2.519658,
   1.978664,
   2.201836,
   1.993019,
   3.701495,
   4.010041,
   4.851867,
   6.578301,
   6.600268,
   8.38121,
   8.28836,
   9.418864,
   9.789748,
   9.741861,
   10.49932,
   9.858199,
   10.68366,
   10.79498,
   10.68496,
   10.74592,
   10.29808};
   Double_t hCrossSection_ojpkuoit_pdf_fehx3033[24] = {
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
   Double_t hCrossSection_ojpkuoit_pdf_fehy3033[24] = {
   3.214863,
   2.778856,
   3.119621,
   1.987158,
   2.447481,
   0.993077,
   1.371119,
   1.083989,
   1.432366,
   2.423457,
   2.973131,
   4.073876,
   4.303735,
   5.680022,
   5.645771,
   6.195652,
   6.339476,
   7.208841,
   7.638982,
   7.859199,
   8.048733,
   8.186892,
   8.31633,
   9.208969};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(24,hCrossSection_ojpkuoit_pdf_fx3033,hCrossSection_ojpkuoit_pdf_fy3033,hCrossSection_ojpkuoit_pdf_felx3033,hCrossSection_ojpkuoit_pdf_fehx3033,hCrossSection_ojpkuoit_pdf_fely3033,hCrossSection_ojpkuoit_pdf_fehy3033);
   grae->SetName("hCrossSection_ojpkuoit_pdf");
   grae->SetTitle("");
   grae->SetFillColor(1);
   grae->Draw("alp");
   c->Modified();
   c->cd();
   c->SetSelected(c);
}

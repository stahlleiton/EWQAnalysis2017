void EPPS16_Cross_Section_Mi()
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
   
   Double_t hCrossSection_hjkrougt_pdf_fx3014[24] = {
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
   Double_t hCrossSection_hjkrougt_pdf_fy3014[24] = {
   110.4337,
   114.162,
   116.3075,
   117.8456,
   119.1601,
   119.4978,
   120.193,
   120.2718,
   120.2248,
   119.5666,
   119.0128,
   118.061,
   116.6198,
   115.1478,
   113.6759,
   111.6003,
   108.6063,
   105.9419,
   103.1932,
   100.0344,
   96.10363,
   92.67233,
   89.08978,
   85.60516};
   Double_t hCrossSection_hjkrougt_pdf_felx3014[24] = {
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
   Double_t hCrossSection_hjkrougt_pdf_fely3014[24] = {
   3.947335,
   4.013794,
   3.811035,
   3.883252,
   3.794598,
   3.287636,
   4.297704,
   3.95138,
   4.435338,
   3.710859,
   4.226774,
   4.434425,
   4.311683,
   4.350871,
   5.10673,
   6.053109,
   5.960967,
   5.797005,
   6.933669,
   7.973524,
   7.585964,
   7.931172,
   8.594581,
   8.058616};
   Double_t hCrossSection_hjkrougt_pdf_fehx3014[24] = {
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
   Double_t hCrossSection_hjkrougt_pdf_fehy3014[24] = {
   4.34022,
   4.006952,
   4.167271,
   3.91298,
   3.861195,
   4.487789,
   3.325002,
   3.60876,
   3.14887,
   3.944057,
   3.815982,
   3.840733,
   4.436207,
   4.810897,
   4.884736,
   5.30371,
   6.098854,
   6.88146,
   6.655506,
   7.19386,
   7.793773,
   7.275551,
   8.08046,
   7.88116};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(24,hCrossSection_hjkrougt_pdf_fx3014,hCrossSection_hjkrougt_pdf_fy3014,hCrossSection_hjkrougt_pdf_felx3014,hCrossSection_hjkrougt_pdf_fehx3014,hCrossSection_hjkrougt_pdf_fely3014,hCrossSection_hjkrougt_pdf_fehy3014);
   grae->SetName("hCrossSection_hjkrougt_pdf");
   grae->SetTitle("");
   grae->SetFillStyle(1000);
   grae->Draw("alp");
   c->Modified();
   c->cd();
   c->SetSelected(c);
}

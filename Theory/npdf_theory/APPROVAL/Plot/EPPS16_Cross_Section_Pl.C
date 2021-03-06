void EPPS16_Cross_Section_Pl()
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
   
   Double_t hCrossSection_ethisyzm_pdf_fx3015[24] = {
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
   Double_t hCrossSection_ethisyzm_pdf_fy3015[24] = {
   92.81647,
   108.0727,
   118.5696,
   127.0627,
   132.9537,
   136.045,
   138.1191,
   139.5036,
   139.4451,
   139.8309,
   139.8183,
   139.3929,
   138.3454,
   137.9691,
   136.785,
   136.5202,
   136.435,
   136.6146,
   136.5059,
   136.3208,
   137.2679,
   136.9469,
   136.2006,
   136.8373};
   Double_t hCrossSection_ethisyzm_pdf_felx3015[24] = {
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
   Double_t hCrossSection_ethisyzm_pdf_fely3015[24] = {
   1.831753,
   2.093265,
   2.155807,
   2.250467,
   2.458154,
   2.914499,
   2.415197,
   2.148763,
   2.102714,
   1.792382,
   1.308308,
   2.466321,
   2.242537,
   3.305921,
   3.582239,
   4.843116,
   6.303669,
   7.74123,
   9.310299,
   9.045404,
   12.06497,
   12.84881,
   12.37461,
   15.27241};
   Double_t hCrossSection_ethisyzm_pdf_fehx3015[24] = {
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
   Double_t hCrossSection_ethisyzm_pdf_fehy3015[24] = {
   2.179169,
   2.530953,
   2.573235,
   3.368663,
   2.802371,
   2.634636,
   2.943098,
   2.871488,
   2.76848,
   2.941822,
   3.411826,
   2.719567,
   3.71463,
   4.827818,
   5.44204,
   6.542002,
   7.577282,
   7.844176,
   9.404215,
   10.85935,
   11.30666,
   11.61708,
   13.5261,
   12.20154};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(24,hCrossSection_ethisyzm_pdf_fx3015,hCrossSection_ethisyzm_pdf_fy3015,hCrossSection_ethisyzm_pdf_felx3015,hCrossSection_ethisyzm_pdf_fehx3015,hCrossSection_ethisyzm_pdf_fely3015,hCrossSection_ethisyzm_pdf_fehy3015);
   grae->SetName("hCrossSection_ethisyzm_pdf");
   grae->SetTitle("");
   grae->SetFillColor(1);
   grae->Draw("alp");
   c->Modified();
   c->cd();
   c->SetSelected(c);
}

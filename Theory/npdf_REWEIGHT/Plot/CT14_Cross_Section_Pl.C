void CT14_Cross_Section_Pl()
{
//=========Macro generated from canvas: c/
//=========  (Tue May 29 12:18:45 2018) by ROOT version6.06/00
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
   
   Double_t hCrossSection_imgocddn_pdf_fx3003[24] = {
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
   Double_t hCrossSection_imgocddn_pdf_fy3003[24] = {
   90.99908,
   104.7654,
   113.7499,
   121.0072,
   126.3165,
   128.5147,
   130.6319,
   132.761,
   133.0481,
   134.2398,
   135.2877,
   136.542,
   136.6821,
   137.6319,
   139.0698,
   139.3454,
   140.5239,
   141.8487,
   144.3964,
   145.0434,
   147.1871,
   149.5057,
   150.1938,
   150.6102};
   Double_t hCrossSection_imgocddn_pdf_felx3003[24] = {
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
   Double_t hCrossSection_imgocddn_pdf_fely3003[24] = {
   1.931395,
   2.364285,
   2.724539,
   2.640958,
   3.495732,
   2.683105,
   3.063488,
   3.425991,
   2.882663,
   3.78596,
   3.150544,
   4.108088,
   3.468202,
   3.536312,
   3.838897,
   3.98106,
   3.266144,
   3.170014,
   4.744156,
   2.676656,
   3.281492,
   3.845813,
   3.107038,
   3.188731};
   Double_t hCrossSection_imgocddn_pdf_fehx3003[24] = {
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
   Double_t hCrossSection_imgocddn_pdf_fehy3003[24] = {
   2.124265,
   2.032804,
   2.214006,
   2.384451,
   1.989127,
   2.738026,
   2.444016,
   2.363435,
   2.932806,
   2.365515,
   2.869634,
   2.025828,
   3.240723,
   2.779882,
   2.497327,
   2.410547,
   2.964431,
   2.963954,
   1.930999,
   3.577031,
   2.783249,
   2.368338,
   3.052334,
   3.242757};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(24,hCrossSection_imgocddn_pdf_fx3003,hCrossSection_imgocddn_pdf_fy3003,hCrossSection_imgocddn_pdf_felx3003,hCrossSection_imgocddn_pdf_fehx3003,hCrossSection_imgocddn_pdf_fely3003,hCrossSection_imgocddn_pdf_fehy3003);
   grae->SetName("hCrossSection_imgocddn_pdf");
   grae->SetTitle("");
   grae->SetFillColor(1);
   grae->Draw("alp");
   c->Modified();
   c->cd();
   c->SetSelected(c);
}

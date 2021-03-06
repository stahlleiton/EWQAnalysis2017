void EPS09_CT14_Cross_Section_Pl()
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
   
   Double_t hCrossSection_omzbvnfn_pdf_fx3027[24] = {
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
   Double_t hCrossSection_omzbvnfn_pdf_fy3027[24] = {
   92.00321,
   107.1899,
   117.8809,
   126.6708,
   132.93,
   135.8058,
   138.707,
   140.222,
   140.1942,
   141.2365,
   141.0639,
   140.5965,
   139.5442,
   139.746,
   138.5238,
   138.4204,
   138.1159,
   138.3116,
   138.3907,
   138.6439,
   139.3801,
   138.9997,
   138.7938,
   139.1391};
   Double_t hCrossSection_omzbvnfn_pdf_felx3027[24] = {
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
   Double_t hCrossSection_omzbvnfn_pdf_fely3027[24] = {
   3.210821,
   3.670204,
   3.765902,
   4.154048,
   4.681032,
   4.209417,
   5.111104,
   4.48955,
   4.437352,
   5.147828,
   4.331086,
   4.64153,
   4.431787,
   5.320336,
   5.531833,
   5.720541,
   5.170595,
   5.575625,
   5.682618,
   6.047601,
   6.430415,
   6.251223,
   6.437518,
   7.19755};
   Double_t hCrossSection_omzbvnfn_pdf_fehx3027[24] = {
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
   Double_t hCrossSection_omzbvnfn_pdf_fehy3027[24] = {
   2.780295,
   3.485311,
   3.579281,
   3.706779,
   3.802949,
   4.042195,
   3.331601,
   3.900102,
   4.252507,
   3.434723,
   4.280733,
   4.017307,
   4.740935,
   3.453858,
   3.933702,
   4.134253,
   4.538164,
   4.714659,
   5.08449,
   5.091054,
   5.600517,
   5.866365,
   6.103996,
   5.595522};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(24,hCrossSection_omzbvnfn_pdf_fx3027,hCrossSection_omzbvnfn_pdf_fy3027,hCrossSection_omzbvnfn_pdf_felx3027,hCrossSection_omzbvnfn_pdf_fehx3027,hCrossSection_omzbvnfn_pdf_fely3027,hCrossSection_omzbvnfn_pdf_fehy3027);
   grae->SetName("hCrossSection_omzbvnfn_pdf");
   grae->SetTitle("");
   grae->SetFillStyle(1000);
   grae->Draw("alp");
   c->Modified();
   c->cd();
   c->SetSelected(c);
}

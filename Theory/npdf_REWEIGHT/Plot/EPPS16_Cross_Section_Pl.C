void EPPS16_Cross_Section_Pl()
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
   
   Double_t hCrossSection_uqzrwqfp_pdf_fx3009[24] = {
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
   Double_t hCrossSection_uqzrwqfp_pdf_fy3009[24] = {
   91.20517,
   106.6266,
   117.6953,
   126.7655,
   133.2896,
   136.7156,
   139.6231,
   141.6359,
   142.3342,
   143.0906,
   143.6685,
   143.3433,
   142.6732,
   142.5168,
   141.2386,
   141.1002,
   140.6044,
   140.6221,
   140.2828,
   139.8824,
   140.5672,
   139.4816,
   139.0391,
   139.0837};
   Double_t hCrossSection_uqzrwqfp_pdf_felx3009[24] = {
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
   Double_t hCrossSection_uqzrwqfp_pdf_fely3009[24] = {
   3.262782,
   3.64101,
   4.110421,
   4.07721,
   4.209622,
   4.515383,
   4.361284,
   4.643046,
   5.075213,
   4.433076,
   4.86279,
   4.791157,
   4.540607,
   4.871522,
   4.632049,
   4.830961,
   4.804475,
   4.785526,
   5.108533,
   4.244768,
   5.304133,
   4.026819,
   4.909916,
   5.373549};
   Double_t hCrossSection_uqzrwqfp_pdf_fehx3009[24] = {
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
   Double_t hCrossSection_uqzrwqfp_pdf_fehy3009[24] = {
   3.154872,
   3.471046,
   3.383264,
   3.680681,
   4.164664,
   3.507601,
   3.888197,
   3.82759,
   3.195037,
   4.086492,
   3.822557,
   3.860963,
   4.185489,
   3.814691,
   4.274632,
   3.964384,
   3.92797,
   3.722661,
   3.551473,
   4.250724,
   3.477753,
   4.674976,
   4.09863,
   3.82896};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(24,hCrossSection_uqzrwqfp_pdf_fx3009,hCrossSection_uqzrwqfp_pdf_fy3009,hCrossSection_uqzrwqfp_pdf_felx3009,hCrossSection_uqzrwqfp_pdf_fehx3009,hCrossSection_uqzrwqfp_pdf_fely3009,hCrossSection_uqzrwqfp_pdf_fehy3009);
   grae->SetName("hCrossSection_uqzrwqfp_pdf");
   grae->SetTitle("");
   grae->SetFillColor(1);
   grae->Draw("alp");
   c->Modified();
   c->cd();
   c->SetSelected(c);
}

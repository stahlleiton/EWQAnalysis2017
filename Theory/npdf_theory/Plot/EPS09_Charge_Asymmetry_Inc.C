void EPS09_Charge_Asymmetry_Inc()
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
   
   Double_t hChargeAsymmetry_tnhochwg_pdf_fx3019[24] = {
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
   Double_t hChargeAsymmetry_tnhochwg_pdf_fy3019[24] = {
   -0.07267374,
   -0.01545634,
   0.02137616,
   0.04696335,
   0.06298075,
   0.07005699,
   0.07526799,
   0.07891734,
   0.07852183,
   0.0814039,
   0.08395101,
   0.08464652,
   0.08682069,
   0.09162742,
   0.09583712,
   0.1042726,
   0.1165445,
   0.1294368,
   0.1429516,
   0.1594765,
   0.1800003,
   0.1955927,
   0.2119318,
   0.2312363};
   Double_t hChargeAsymmetry_tnhochwg_pdf_felx3019[24] = {
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
   Double_t hChargeAsymmetry_tnhochwg_pdf_fely3019[24] = {
   0.0177547,
   0.005418745,
   0.01092488,
   0.004993188,
   0.00680966,
   0.004758178,
   0.004372587,
   0.003853377,
   0.005821295,
   0.003726355,
   0.007356726,
   0.004565965,
   0.005321917,
   0.00405479,
   0.006957087,
   0.006138517,
   0.005460634,
   0.006996384,
   0.007116556,
   0.009666016,
   0.008557266,
   0.006548771,
   0.009024151,
   0.008777784};
   Double_t hChargeAsymmetry_tnhochwg_pdf_fehx3019[24] = {
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
   Double_t hChargeAsymmetry_tnhochwg_pdf_fehy3019[24] = {
   0.006347949,
   0.009050014,
   0.004410186,
   0.009665799,
   0.00928624,
   0.005786504,
   0.005672104,
   0.006243966,
   0.003147399,
   0.0056536,
   0.003768523,
   0.005176983,
   0.006058768,
   0.008255532,
   0.004628761,
   0.006885335,
   0.006131568,
   0.006494814,
   0.007258941,
   0.004690806,
   0.007082197,
   0.009200171,
   0.008488748,
   0.01069308};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(24,hChargeAsymmetry_tnhochwg_pdf_fx3019,hChargeAsymmetry_tnhochwg_pdf_fy3019,hChargeAsymmetry_tnhochwg_pdf_felx3019,hChargeAsymmetry_tnhochwg_pdf_fehx3019,hChargeAsymmetry_tnhochwg_pdf_fely3019,hChargeAsymmetry_tnhochwg_pdf_fehy3019);
   grae->SetName("hChargeAsymmetry_tnhochwg_pdf");
   grae->SetTitle("");
   grae->SetFillStyle(1000);
   grae->Draw("alp");
   c->Modified();
   c->cd();
   c->SetSelected(c);
}

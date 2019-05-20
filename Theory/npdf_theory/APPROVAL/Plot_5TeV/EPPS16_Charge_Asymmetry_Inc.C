void EPPS16_Charge_Asymmetry_Inc()
{
//=========Macro generated from canvas: c/
//=========  (Thu May  3 16:25:11 2018) by ROOT version6.06/00
   TCanvas *c = new TCanvas("c", "",0,0,600,600);
   c->Range(0,0,1,1);
   c->SetFillColor(0);
   c->SetBorderMode(0);
   c->SetBorderSize(2);
   c->SetFrameBorderMode(0);
   
   Double_t hChargeAsymmetry_usfbjrib_pdf_fx3001[10] = {
   -2.2,
   -1.75,
   -1.25,
   -0.75,
   -0.25,
   0.25,
   0.75,
   1.25,
   1.75,
   2.2};
   Double_t hChargeAsymmetry_usfbjrib_pdf_fy3001[10] = {
   -0.2336219,
   -0.06014597,
   0.04517266,
   0.09224602,
   0.1132514,
   0.1284726,
   0.1503874,
   0.1867246,
   0.2307891,
   0.2603281};
   Double_t hChargeAsymmetry_usfbjrib_pdf_felx3001[10] = {
   0.2,
   0.25,
   0.25,
   0.25,
   0.25,
   0.25,
   0.25,
   0.25,
   0.25,
   0.2};
   Double_t hChargeAsymmetry_usfbjrib_pdf_fely3001[10] = {
   0.0244727,
   0.02410125,
   0.01962496,
   0.01396378,
   0.008932523,
   0.005300171,
   0.006142093,
   0.009793986,
   0.01204817,
   0.01241369};
   Double_t hChargeAsymmetry_usfbjrib_pdf_fehx3001[10] = {
   0.2,
   0.25,
   0.25,
   0.25,
   0.25,
   0.25,
   0.25,
   0.25,
   0.25,
   0.2};
   Double_t hChargeAsymmetry_usfbjrib_pdf_fehy3001[10] = {
   0.02524164,
   0.02342904,
   0.02015705,
   0.01373849,
   0.008144228,
   0.007470627,
   0.009097115,
   0.01005905,
   0.01298757,
   0.01487545};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(10,hChargeAsymmetry_usfbjrib_pdf_fx3001,hChargeAsymmetry_usfbjrib_pdf_fy3001,hChargeAsymmetry_usfbjrib_pdf_felx3001,hChargeAsymmetry_usfbjrib_pdf_fehx3001,hChargeAsymmetry_usfbjrib_pdf_fely3001,hChargeAsymmetry_usfbjrib_pdf_fehy3001);
   grae->SetName("hChargeAsymmetry_usfbjrib_pdf");
   grae->SetTitle("");
   grae->SetFillColor(1);
   grae->Draw("alp");
   c->Modified();
   c->cd();
   c->SetSelected(c);
}

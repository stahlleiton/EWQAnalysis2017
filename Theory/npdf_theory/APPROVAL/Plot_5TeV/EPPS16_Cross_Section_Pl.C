void EPPS16_Cross_Section_Pl()
{
//=========Macro generated from canvas: c/
//=========  (Thu May  3 16:25:11 2018) by ROOT version6.06/00
   TCanvas *c = new TCanvas("c", "",0,0,600,600);
   c->Range(0,0,1,1);
   c->SetFillColor(0);
   c->SetBorderMode(0);
   c->SetBorderSize(2);
   c->SetFrameBorderMode(0);
   
   Double_t hCrossSection_ppeguzkf_pdf_fx3003[10] = {
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
   Double_t hCrossSection_ppeguzkf_pdf_fy3003[10] = {
   42.18692,
   66.15761,
   85.01131,
   94.07329,
   96.98965,
   97.20647,
   96.67413,
   96.76384,
   96.13419,
   92.13782};
   Double_t hCrossSection_ppeguzkf_pdf_felx3003[10] = {
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
   Double_t hCrossSection_ppeguzkf_pdf_fely3003[10] = {
   1.653754,
   2.493004,
   3.125351,
   3.250007,
   3.354755,
   3.279073,
   3.52704,
   4.943047,
   6.65084,
   7.836764};
   Double_t hCrossSection_ppeguzkf_pdf_fehx3003[10] = {
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
   Double_t hCrossSection_ppeguzkf_pdf_fehy3003[10] = {
   1.768937,
   2.495407,
   2.909781,
   3.005679,
   3.077289,
   3.04223,
   3.751912,
   5.106617,
   7.007804,
   8.226671};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(10,hCrossSection_ppeguzkf_pdf_fx3003,hCrossSection_ppeguzkf_pdf_fy3003,hCrossSection_ppeguzkf_pdf_felx3003,hCrossSection_ppeguzkf_pdf_fehx3003,hCrossSection_ppeguzkf_pdf_fely3003,hCrossSection_ppeguzkf_pdf_fehy3003);
   grae->SetName("hCrossSection_ppeguzkf_pdf");
   grae->SetTitle("");
   grae->SetFillColor(1);
   grae->Draw("alp");
   c->Modified();
   c->cd();
   c->SetSelected(c);
}

void EPPS16_ForwardBackward_Ratio_Pl()
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
   
   Double_t hForwardBackwardRatio_gidigbmu_pdf_fx3012[10] = {
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
   Double_t hForwardBackwardRatio_gidigbmu_pdf_fy3012[10] = {
   0.9910318,
   0.9889746,
   0.9808924,
   0.9787956,
   0.9803774,
   0.9827746,
   0.9924549,
   0.9989862,
   1.016995,
   1.04347};
   Double_t hForwardBackwardRatio_gidigbmu_pdf_felx3012[10] = {
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
   Double_t hForwardBackwardRatio_gidigbmu_pdf_fely3012[10] = {
   0.004142725,
   0.008436369,
   0.0126311,
   0.01436482,
   0.02372564,
   0.01661723,
   0.02921637,
   0.02497345,
   0.0288189,
   0.04034488};
   Double_t hForwardBackwardRatio_gidigbmu_pdf_fehx3012[10] = {
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
   Double_t hForwardBackwardRatio_gidigbmu_pdf_fehy3012[10] = {
   0.008231878,
   0.006858576,
   0.01212631,
   0.01247248,
   0.0160127,
   0.02705874,
   0.01904992,
   0.02975524,
   0.02922635,
   0.02677191};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(10,hForwardBackwardRatio_gidigbmu_pdf_fx3012,hForwardBackwardRatio_gidigbmu_pdf_fy3012,hForwardBackwardRatio_gidigbmu_pdf_felx3012,hForwardBackwardRatio_gidigbmu_pdf_fehx3012,hForwardBackwardRatio_gidigbmu_pdf_fely3012,hForwardBackwardRatio_gidigbmu_pdf_fehy3012);
   grae->SetName("hForwardBackwardRatio_gidigbmu_pdf");
   grae->SetTitle("");
   grae->SetFillColor(1);
   grae->Draw("alp");
   c->Modified();
   c->cd();
   c->SetSelected(c);
}

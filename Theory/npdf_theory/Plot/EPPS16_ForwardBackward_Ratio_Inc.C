void EPPS16_ForwardBackward_Ratio_Inc()
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
   
   Double_t hForwardBackwardRatio_Inc_ekmoecna_pdf_fx3016[10] = {
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
   Double_t hForwardBackwardRatio_Inc_ekmoecna_pdf_fy3016[10] = {
   0.9895072,
   0.9731542,
   0.9517868,
   0.9371233,
   0.924061,
   0.9102143,
   0.8983585,
   0.888922,
   0.881615,
   0.8823097};
   Double_t hForwardBackwardRatio_Inc_ekmoecna_pdf_felx3016[10] = {
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
   Double_t hForwardBackwardRatio_Inc_ekmoecna_pdf_fely3016[10] = {
   0.0085249,
   0.02804175,
   0.03859414,
   0.05390298,
   0.07009402,
   0.07190552,
   0.08657513,
   0.0919108,
   0.09479006,
   0.1033093};
   Double_t hForwardBackwardRatio_Inc_ekmoecna_pdf_fehx3016[10] = {
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
   Double_t hForwardBackwardRatio_Inc_ekmoecna_pdf_fehy3016[10] = {
   0.007731235,
   0.02467296,
   0.04263093,
   0.04903309,
   0.0624812,
   0.07655125,
   0.08128097,
   0.08408839,
   0.0955562,
   0.09006462};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(10,hForwardBackwardRatio_Inc_ekmoecna_pdf_fx3016,hForwardBackwardRatio_Inc_ekmoecna_pdf_fy3016,hForwardBackwardRatio_Inc_ekmoecna_pdf_felx3016,hForwardBackwardRatio_Inc_ekmoecna_pdf_fehx3016,hForwardBackwardRatio_Inc_ekmoecna_pdf_fely3016,hForwardBackwardRatio_Inc_ekmoecna_pdf_fehy3016);
   grae->SetName("hForwardBackwardRatio_Inc_ekmoecna_pdf");
   grae->SetTitle("");
   grae->SetFillStyle(1000);
   grae->Draw("alp");
   c->Modified();
   c->cd();
   c->SetSelected(c);
}

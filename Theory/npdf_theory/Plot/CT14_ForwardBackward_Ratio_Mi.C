void CT14_ForwardBackward_Ratio_Mi()
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
   
   Double_t hForwardBackwardRatio_tthzdffw_pdf_fx3011[10] = {
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
   Double_t hForwardBackwardRatio_tthzdffw_pdf_fy3011[10] = {
   0.9976748,
   0.9837714,
   0.9724824,
   0.9523835,
   0.9424112,
   0.9265168,
   0.9049227,
   0.883846,
   0.8652766,
   0.8466122};
   Double_t hForwardBackwardRatio_tthzdffw_pdf_felx3011[10] = {
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
   Double_t hForwardBackwardRatio_tthzdffw_pdf_fely3011[10] = {
   0.009088434,
   0.001989348,
   0.01223937,
   0.002768657,
   0.009338371,
   0.01159086,
   0.0143496,
   0.0120747,
   0.01507562,
   0.01717654};
   Double_t hForwardBackwardRatio_tthzdffw_pdf_fehx3011[10] = {
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
   Double_t hForwardBackwardRatio_tthzdffw_pdf_fehy3011[10] = {
   0.001173679,
   0.006370311,
   0.001980949,
   0.01358325,
   0.006850326,
   0.009747164,
   0.01004261,
   0.01370684,
   0.01322472,
   0.01986967};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(10,hForwardBackwardRatio_tthzdffw_pdf_fx3011,hForwardBackwardRatio_tthzdffw_pdf_fy3011,hForwardBackwardRatio_tthzdffw_pdf_felx3011,hForwardBackwardRatio_tthzdffw_pdf_fehx3011,hForwardBackwardRatio_tthzdffw_pdf_fely3011,hForwardBackwardRatio_tthzdffw_pdf_fehy3011);
   grae->SetName("hForwardBackwardRatio_tthzdffw_pdf");
   grae->SetTitle("");
   grae->SetFillStyle(1000);
   grae->Draw("alp");
   c->Modified();
   c->cd();
   c->SetSelected(c);
}

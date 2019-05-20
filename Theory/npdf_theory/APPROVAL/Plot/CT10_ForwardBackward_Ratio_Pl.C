void CT10_ForwardBackward_Ratio_Pl()
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
   
   Double_t hForwardBackwardRatio_jqpoimmj_pdf_fx3006[10] = {
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
   Double_t hForwardBackwardRatio_jqpoimmj_pdf_fy3006[10] = {
   1.010469,
   1.021855,
   1.038422,
   1.059533,
   1.084677,
   1.103164,
   1.12444,
   1.155112,
   1.178784,
   1.204964};
   Double_t hForwardBackwardRatio_jqpoimmj_pdf_felx3006[10] = {
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
   Double_t hForwardBackwardRatio_jqpoimmj_pdf_fely3006[10] = {
   0.006252135,
   0.007565259,
   0.003886278,
   0.00594442,
   0.006688132,
   0.01162589,
   0.008920767,
   0.01530979,
   0.01715712,
   0.01064654};
   Double_t hForwardBackwardRatio_jqpoimmj_pdf_fehx3006[10] = {
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
   Double_t hForwardBackwardRatio_jqpoimmj_pdf_fehy3006[10] = {
   0.005355547,
   0.002420648,
   0.007231001,
   0.006522823,
   0.00743943,
   0.008206959,
   0.01510528,
   0.009511108,
   0.01145872,
   0.01921176};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(10,hForwardBackwardRatio_jqpoimmj_pdf_fx3006,hForwardBackwardRatio_jqpoimmj_pdf_fy3006,hForwardBackwardRatio_jqpoimmj_pdf_felx3006,hForwardBackwardRatio_jqpoimmj_pdf_fehx3006,hForwardBackwardRatio_jqpoimmj_pdf_fely3006,hForwardBackwardRatio_jqpoimmj_pdf_fehy3006);
   grae->SetName("hForwardBackwardRatio_jqpoimmj_pdf");
   grae->SetTitle("");
   grae->SetFillColor(1);
   grae->Draw("alp");
   c->Modified();
   c->cd();
   c->SetSelected(c);
}

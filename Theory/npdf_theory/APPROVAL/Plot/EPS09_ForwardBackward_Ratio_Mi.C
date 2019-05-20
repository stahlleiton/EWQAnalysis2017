void EPS09_ForwardBackward_Ratio_Mi()
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
   
   Double_t hForwardBackwardRatio_ilaybwhf_pdf_fx3023[10] = {
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
   Double_t hForwardBackwardRatio_ilaybwhf_pdf_fy3023[10] = {
   0.9837553,
   0.9528332,
   0.9194409,
   0.890075,
   0.862054,
   0.8341593,
   0.8054838,
   0.7794476,
   0.7579384,
   0.7339513};
   Double_t hForwardBackwardRatio_ilaybwhf_pdf_felx3023[10] = {
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
   Double_t hForwardBackwardRatio_ilaybwhf_pdf_fely3023[10] = {
   0.002386336,
   0.008295354,
   0.01923491,
   0.02273036,
   0.02546945,
   0.03049337,
   0.03339399,
   0.03272972,
   0.03814735,
   0.03633109};
   Double_t hForwardBackwardRatio_ilaybwhf_pdf_fehx3023[10] = {
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
   Double_t hForwardBackwardRatio_ilaybwhf_pdf_fehy3023[10] = {
   0.01122874,
   0.01696362,
   0.01781185,
   0.0217923,
   0.02666901,
   0.03282793,
   0.03025761,
   0.03217068,
   0.03309561,
   0.03361428};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(10,hForwardBackwardRatio_ilaybwhf_pdf_fx3023,hForwardBackwardRatio_ilaybwhf_pdf_fy3023,hForwardBackwardRatio_ilaybwhf_pdf_felx3023,hForwardBackwardRatio_ilaybwhf_pdf_fehx3023,hForwardBackwardRatio_ilaybwhf_pdf_fely3023,hForwardBackwardRatio_ilaybwhf_pdf_fehy3023);
   grae->SetName("hForwardBackwardRatio_ilaybwhf_pdf");
   grae->SetTitle("");
   grae->SetFillColor(1);
   grae->Draw("alp");
   c->Modified();
   c->cd();
   c->SetSelected(c);
}

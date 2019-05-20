void CT14_ForwardBackward_Ratio_Pl()
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
   
   Double_t hForwardBackwardRatio_beumhmzo_pdf_fx3006[10] = {
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
   Double_t hForwardBackwardRatio_beumhmzo_pdf_fy3006[10] = {
   1.010447,
   1.019486,
   1.029163,
   1.048496,
   1.07566,
   1.090158,
   1.108663,
   1.144481,
   1.168689,
   1.192324};
   Double_t hForwardBackwardRatio_beumhmzo_pdf_felx3006[10] = {
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
   Double_t hForwardBackwardRatio_beumhmzo_pdf_fely3006[10] = {
   0.0112952,
   0.01382538,
   0.000593273,
   0.005932192,
   0.01143617,
   0.003735916,
   0.004130421,
   0.01438591,
   0.008118435,
   0.002863365};
   Double_t hForwardBackwardRatio_beumhmzo_pdf_fehx3006[10] = {
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
   Double_t hForwardBackwardRatio_beumhmzo_pdf_fehy3006[10] = {
   0.003055784,
   0.007372066,
   0.01230471,
   0.005338841,
   0.001480129,
   0.009886814,
   0.009495266,
   0.005059516,
   0.005325548,
   0.01952472};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(10,hForwardBackwardRatio_beumhmzo_pdf_fx3006,hForwardBackwardRatio_beumhmzo_pdf_fy3006,hForwardBackwardRatio_beumhmzo_pdf_felx3006,hForwardBackwardRatio_beumhmzo_pdf_fehx3006,hForwardBackwardRatio_beumhmzo_pdf_fely3006,hForwardBackwardRatio_beumhmzo_pdf_fehy3006);
   grae->SetName("hForwardBackwardRatio_beumhmzo_pdf");
   grae->SetTitle("");
   grae->SetFillColor(1);
   grae->Draw("alp");
   c->Modified();
   c->cd();
   c->SetSelected(c);
}

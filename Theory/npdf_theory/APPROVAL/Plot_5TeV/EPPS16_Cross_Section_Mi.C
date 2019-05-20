void EPPS16_Cross_Section_Mi()
{
//=========Macro generated from canvas: c/
//=========  (Thu May  3 16:25:11 2018) by ROOT version6.06/00
   TCanvas *c = new TCanvas("c", "",0,0,600,600);
   c->Range(0,0,1,1);
   c->SetFillColor(0);
   c->SetBorderMode(0);
   c->SetBorderSize(2);
   c->SetFrameBorderMode(0);
   
   Double_t hCrossSection_okmjehgq_pdf_fx3002[10] = {
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
   Double_t hCrossSection_okmjehgq_pdf_fy3002[10] = {
   67.90736,
   74.62512,
   77.66289,
   78.1833,
   77.25607,
   75.07324,
   71.39817,
   66.31332,
   60.08134,
   54.07461};
   Double_t hCrossSection_okmjehgq_pdf_felx3002[10] = {
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
   Double_t hCrossSection_okmjehgq_pdf_fely3002[10] = {
   2.456047,
   2.547761,
   2.634702,
   2.487857,
   2.573399,
   2.80536,
   2.868756,
   3.004656,
   3.58028,
   4.08787};
   Double_t hCrossSection_okmjehgq_pdf_fehx3002[10] = {
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
   Double_t hCrossSection_okmjehgq_pdf_fehy3002[10] = {
   2.673406,
   2.758138,
   2.462283,
   2.346367,
   2.261493,
   2.168716,
   2.523086,
   3.031841,
   3.67035,
   4.067623};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(10,hCrossSection_okmjehgq_pdf_fx3002,hCrossSection_okmjehgq_pdf_fy3002,hCrossSection_okmjehgq_pdf_felx3002,hCrossSection_okmjehgq_pdf_fehx3002,hCrossSection_okmjehgq_pdf_fely3002,hCrossSection_okmjehgq_pdf_fehy3002);
   grae->SetName("hCrossSection_okmjehgq_pdf");
   grae->SetTitle("");
   grae->SetFillColor(1);
   grae->Draw("alp");
   c->Modified();
   c->cd();
   c->SetSelected(c);
}

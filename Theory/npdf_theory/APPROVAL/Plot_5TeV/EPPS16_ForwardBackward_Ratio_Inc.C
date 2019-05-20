void EPPS16_ForwardBackward_Ratio_Inc()
{
//=========Macro generated from canvas: c/
//=========  (Thu May  3 16:25:11 2018) by ROOT version6.06/00
   TCanvas *c = new TCanvas("c", "",0,0,600,600);
   c->Range(0,0,1,1);
   c->SetFillColor(0);
   c->SetBorderMode(0);
   c->SetBorderSize(2);
   c->SetFrameBorderMode(0);
   
   Double_t hForwardBackwardRatio_Inc_ziobvoay_pdf_fx3004[5] = {
   0.25,
   0.75,
   1.25,
   1.75,
   2.2};
   Double_t hForwardBackwardRatio_Inc_ziobvoay_pdf_fy3004[5] = {
   0.988717,
   0.975709,
   1.002477,
   1.109621,
   1.328066};
   Double_t hForwardBackwardRatio_Inc_ziobvoay_pdf_felx3004[5] = {
   0.25,
   0.25,
   0.25,
   0.25,
   0.2};
   Double_t hForwardBackwardRatio_Inc_ziobvoay_pdf_fely3004[5] = {
   0.01335378,
   0.03228396,
   0.05401767,
   0.08152985,
   0.1181682};
   Double_t hForwardBackwardRatio_Inc_ziobvoay_pdf_fehx3004[5] = {
   0.25,
   0.25,
   0.25,
   0.25,
   0.2};
   Double_t hForwardBackwardRatio_Inc_ziobvoay_pdf_fehy3004[5] = {
   0.01433975,
   0.03588067,
   0.0582416,
   0.0814554,
   0.1145925};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(5,hForwardBackwardRatio_Inc_ziobvoay_pdf_fx3004,hForwardBackwardRatio_Inc_ziobvoay_pdf_fy3004,hForwardBackwardRatio_Inc_ziobvoay_pdf_felx3004,hForwardBackwardRatio_Inc_ziobvoay_pdf_fehx3004,hForwardBackwardRatio_Inc_ziobvoay_pdf_fely3004,hForwardBackwardRatio_Inc_ziobvoay_pdf_fehy3004);
   grae->SetName("hForwardBackwardRatio_Inc_ziobvoay_pdf");
   grae->SetTitle("");
   grae->SetFillColor(1);
   grae->Draw("alp");
   c->Modified();
   c->cd();
   c->SetSelected(c);
}

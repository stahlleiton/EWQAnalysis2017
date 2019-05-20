void EPPS16_ForwardBackward_Ratio_Mi()
{
//=========Macro generated from canvas: c/
//=========  (Thu May  3 16:25:11 2018) by ROOT version6.06/00
   TCanvas *c = new TCanvas("c", "",0,0,600,600);
   c->Range(0,0,1,1);
   c->SetFillColor(0);
   c->SetBorderMode(0);
   c->SetBorderSize(2);
   c->SetFrameBorderMode(0);
   
   Double_t hForwardBackwardRatio_wxclstcy_pdf_fx3005[5] = {
   0.25,
   0.75,
   1.25,
   1.75,
   2.2};
   Double_t hForwardBackwardRatio_wxclstcy_pdf_fy3005[5] = {
   0.9717455,
   0.913215,
   0.8538611,
   0.8051088,
   0.7962998};
   Double_t hForwardBackwardRatio_wxclstcy_pdf_felx3005[5] = {
   0.25,
   0.25,
   0.25,
   0.25,
   0.2};
   Double_t hForwardBackwardRatio_wxclstcy_pdf_fely3005[5] = {
   0.01404678,
   0.0344003,
   0.04630017,
   0.05957893,
   0.07066451};
   Double_t hForwardBackwardRatio_wxclstcy_pdf_fehx3005[5] = {
   0.25,
   0.25,
   0.25,
   0.25,
   0.2};
   Double_t hForwardBackwardRatio_wxclstcy_pdf_fehy3005[5] = {
   0.0113514,
   0.0348643,
   0.05025444,
   0.05781033,
   0.06726458};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(5,hForwardBackwardRatio_wxclstcy_pdf_fx3005,hForwardBackwardRatio_wxclstcy_pdf_fy3005,hForwardBackwardRatio_wxclstcy_pdf_felx3005,hForwardBackwardRatio_wxclstcy_pdf_fehx3005,hForwardBackwardRatio_wxclstcy_pdf_fely3005,hForwardBackwardRatio_wxclstcy_pdf_fehy3005);
   grae->SetName("hForwardBackwardRatio_wxclstcy_pdf");
   grae->SetTitle("");
   grae->SetFillColor(1);
   grae->Draw("alp");
   c->Modified();
   c->cd();
   c->SetSelected(c);
}

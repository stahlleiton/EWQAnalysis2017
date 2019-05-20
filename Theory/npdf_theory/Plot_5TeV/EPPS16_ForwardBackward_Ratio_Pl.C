void EPPS16_ForwardBackward_Ratio_Pl()
{
//=========Macro generated from canvas: c/
//=========  (Thu May  3 16:25:11 2018) by ROOT version6.06/00
   TCanvas *c = new TCanvas("c", "",0,0,600,600);
   c->Range(0,0,1,1);
   c->SetFillColor(0);
   c->SetBorderMode(0);
   c->SetBorderSize(2);
   c->SetFrameBorderMode(0);
   
   Double_t hForwardBackwardRatio_havggtzh_pdf_fx3006[5] = {
   0.25,
   0.75,
   1.25,
   1.75,
   2.2};
   Double_t hForwardBackwardRatio_havggtzh_pdf_fy3006[5] = {
   1.002235,
   1.027647,
   1.138247,
   1.453109,
   2.184037};
   Double_t hForwardBackwardRatio_havggtzh_pdf_felx3006[5] = {
   0.25,
   0.25,
   0.25,
   0.25,
   0.2};
   Double_t hForwardBackwardRatio_havggtzh_pdf_fely3006[5] = {
   0.01543372,
   0.03820225,
   0.072812,
   0.1221599,
   0.2166946};
   Double_t hForwardBackwardRatio_havggtzh_pdf_fehx3006[5] = {
   0.25,
   0.25,
   0.25,
   0.25,
   0.2};
   Double_t hForwardBackwardRatio_havggtzh_pdf_fehy3006[5] = {
   0.01983011,
   0.04553022,
   0.07795836,
   0.1259766,
   0.2151431};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(5,hForwardBackwardRatio_havggtzh_pdf_fx3006,hForwardBackwardRatio_havggtzh_pdf_fy3006,hForwardBackwardRatio_havggtzh_pdf_felx3006,hForwardBackwardRatio_havggtzh_pdf_fehx3006,hForwardBackwardRatio_havggtzh_pdf_fely3006,hForwardBackwardRatio_havggtzh_pdf_fehy3006);
   grae->SetName("hForwardBackwardRatio_havggtzh_pdf");
   grae->SetTitle("");
   grae->SetFillColor(1);
   grae->Draw("alp");
   c->Modified();
   c->cd();
   c->SetSelected(c);
}

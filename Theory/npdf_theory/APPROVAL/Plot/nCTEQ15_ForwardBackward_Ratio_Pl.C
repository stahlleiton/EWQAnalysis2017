void nCTEQ15_ForwardBackward_Ratio_Pl()
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
   
   Double_t hForwardBackwardRatio_ekjfctkc_pdf_fx3036[10] = {
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
   Double_t hForwardBackwardRatio_ekjfctkc_pdf_fy3036[10] = {
   0.9907628,
   0.9578304,
   0.9405574,
   0.9166374,
   0.8995518,
   0.8937631,
   0.8827789,
   0.8803835,
   0.8832308,
   0.8866808};
   Double_t hForwardBackwardRatio_ekjfctkc_pdf_felx3036[10] = {
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
   Double_t hForwardBackwardRatio_ekjfctkc_pdf_fely3036[10] = {
   0.01165769,
   0.01553326,
   0.02876073,
   0.03742515,
   0.04430688,
   0.05862906,
   0.06216219,
   0.07583989,
   0.07779715,
   0.08817004};
   Double_t hForwardBackwardRatio_ekjfctkc_pdf_fehx3036[10] = {
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
   Double_t hForwardBackwardRatio_ekjfctkc_pdf_fehy3036[10] = {
   0.003285852,
   0.02098782,
   0.02680482,
   0.04256571,
   0.04950818,
   0.05836574,
   0.0629399,
   0.06619166,
   0.07114692,
   0.07903761};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(10,hForwardBackwardRatio_ekjfctkc_pdf_fx3036,hForwardBackwardRatio_ekjfctkc_pdf_fy3036,hForwardBackwardRatio_ekjfctkc_pdf_felx3036,hForwardBackwardRatio_ekjfctkc_pdf_fehx3036,hForwardBackwardRatio_ekjfctkc_pdf_fely3036,hForwardBackwardRatio_ekjfctkc_pdf_fehy3036);
   grae->SetName("hForwardBackwardRatio_ekjfctkc_pdf");
   grae->SetTitle("");
   grae->SetFillColor(1);
   grae->Draw("alp");
   c->Modified();
   c->cd();
   c->SetSelected(c);
}

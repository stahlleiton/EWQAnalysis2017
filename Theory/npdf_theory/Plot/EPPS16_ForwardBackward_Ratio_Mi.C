void EPPS16_ForwardBackward_Ratio_Mi()
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
   
   Double_t hForwardBackwardRatio_ixxeubkg_pdf_fx3017[10] = {
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
   Double_t hForwardBackwardRatio_ixxeubkg_pdf_fy3017[10] = {
   0.9872179,
   0.9569585,
   0.9199166,
   0.890173,
   0.86306,
   0.8320614,
   0.7990534,
   0.7710296,
   0.7455351,
   0.7184044};
   Double_t hForwardBackwardRatio_ixxeubkg_pdf_felx3017[10] = {
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
   Double_t hForwardBackwardRatio_ixxeubkg_pdf_fely3017[10] = {
   0.01300169,
   0.03030172,
   0.04069332,
   0.04497633,
   0.06258313,
   0.06930585,
   0.07271978,
   0.07499482,
   0.08638873,
   0.0761247};
   Double_t hForwardBackwardRatio_ixxeubkg_pdf_fehx3017[10] = {
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
   Double_t hForwardBackwardRatio_ixxeubkg_pdf_fehy3017[10] = {
   0.006432106,
   0.0189076,
   0.03730151,
   0.04823601,
   0.05098133,
   0.06535122,
   0.07016492,
   0.06988732,
   0.07254623,
   0.07385724};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(10,hForwardBackwardRatio_ixxeubkg_pdf_fx3017,hForwardBackwardRatio_ixxeubkg_pdf_fy3017,hForwardBackwardRatio_ixxeubkg_pdf_felx3017,hForwardBackwardRatio_ixxeubkg_pdf_fehx3017,hForwardBackwardRatio_ixxeubkg_pdf_fely3017,hForwardBackwardRatio_ixxeubkg_pdf_fehy3017);
   grae->SetName("hForwardBackwardRatio_ixxeubkg_pdf");
   grae->SetTitle("");
   grae->SetFillStyle(1000);
   grae->Draw("alp");
   c->Modified();
   c->cd();
   c->SetSelected(c);
}

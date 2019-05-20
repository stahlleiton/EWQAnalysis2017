void EPS09_CT14_ForwardBackward_Ratio_Mi()
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
   
   Double_t hForwardBackwardRatio_pzstpxwe_pdf_fx3029[10] = {
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
   Double_t hForwardBackwardRatio_pzstpxwe_pdf_fy3029[10] = {
   0.9918966,
   0.9540738,
   0.9171107,
   0.8916354,
   0.8665023,
   0.8356045,
   0.8056865,
   0.7801462,
   0.7561284,
   0.7291713};
   Double_t hForwardBackwardRatio_pzstpxwe_pdf_felx3029[10] = {
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
   Double_t hForwardBackwardRatio_pzstpxwe_pdf_fely3029[10] = {
   0.02619931,
   0.008564445,
   0.01331265,
   0.01942191,
   0.03155084,
   0.02928062,
   0.03188445,
   0.03622788,
   0.0415365,
   0.0332725};
   Double_t hForwardBackwardRatio_pzstpxwe_pdf_fehx3029[10] = {
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
   Double_t hForwardBackwardRatio_pzstpxwe_pdf_fehy3029[10] = {
   0.00038125,
   0.01747923,
   0.02870738,
   0.02558463,
   0.02301984,
   0.0343626,
   0.03563939,
   0.04467325,
   0.03341935,
   0.03927142};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(10,hForwardBackwardRatio_pzstpxwe_pdf_fx3029,hForwardBackwardRatio_pzstpxwe_pdf_fy3029,hForwardBackwardRatio_pzstpxwe_pdf_felx3029,hForwardBackwardRatio_pzstpxwe_pdf_fehx3029,hForwardBackwardRatio_pzstpxwe_pdf_fely3029,hForwardBackwardRatio_pzstpxwe_pdf_fehy3029);
   grae->SetName("hForwardBackwardRatio_pzstpxwe_pdf");
   grae->SetTitle("");
   grae->SetFillStyle(1000);
   grae->Draw("alp");
   c->Modified();
   c->cd();
   c->SetSelected(c);
}

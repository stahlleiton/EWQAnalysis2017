void nCTEQ15_ForwardBackward_Ratio_Pl()
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
   
   Double_t hForwardBackwardRatio_cwrmdfds_pdf_fx3036[10] = {
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
   Double_t hForwardBackwardRatio_cwrmdfds_pdf_fy3036[10] = {
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
   Double_t hForwardBackwardRatio_cwrmdfds_pdf_felx3036[10] = {
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
   Double_t hForwardBackwardRatio_cwrmdfds_pdf_fely3036[10] = {
   0.01549142,
   0.01690008,
   0.03103619,
   0.04083446,
   0.04661636,
   0.06227018,
   0.06487022,
   0.08019096,
   0.08445629,
   0.09188672};
   Double_t hForwardBackwardRatio_cwrmdfds_pdf_fehx3036[10] = {
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
   Double_t hForwardBackwardRatio_cwrmdfds_pdf_fehy3036[10] = {
   0.003567235,
   0.02234386,
   0.02789695,
   0.04446537,
   0.05250926,
   0.06142352,
   0.06723356,
   0.06954773,
   0.07526133,
   0.08390271};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(10,hForwardBackwardRatio_cwrmdfds_pdf_fx3036,hForwardBackwardRatio_cwrmdfds_pdf_fy3036,hForwardBackwardRatio_cwrmdfds_pdf_felx3036,hForwardBackwardRatio_cwrmdfds_pdf_fehx3036,hForwardBackwardRatio_cwrmdfds_pdf_fely3036,hForwardBackwardRatio_cwrmdfds_pdf_fehy3036);
   grae->SetName("hForwardBackwardRatio_cwrmdfds_pdf");
   grae->SetTitle("");
   grae->SetFillStyle(1000);
   grae->Draw("alp");
   c->Modified();
   c->cd();
   c->SetSelected(c);
}

void EPS09_CT14_ForwardBackward_Ratio_Inc()
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
   
   Double_t hForwardBackwardRatio_Inc_uvutiivw_pdf_fx3028[10] = {
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
   Double_t hForwardBackwardRatio_Inc_uvutiivw_pdf_fy3028[10] = {
   0.9915463,
   0.9745867,
   0.9523693,
   0.9395897,
   0.9275965,
   0.917851,
   0.906818,
   0.8989404,
   0.8975152,
   0.8968738};
   Double_t hForwardBackwardRatio_Inc_uvutiivw_pdf_felx3028[10] = {
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
   Double_t hForwardBackwardRatio_Inc_uvutiivw_pdf_fely3028[10] = {
   0.01331945,
   0.0137121,
   0.01609379,
   0.02376349,
   0.0309862,
   0.03433627,
   0.03792077,
   0.03877762,
   0.04339108,
   0.04312554};
   Double_t hForwardBackwardRatio_Inc_uvutiivw_pdf_fehx3028[10] = {
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
   Double_t hForwardBackwardRatio_Inc_uvutiivw_pdf_fehy3028[10] = {
   0.00259253,
   0.01099324,
   0.02360109,
   0.02353227,
   0.02964918,
   0.03265796,
   0.037134,
   0.04258141,
   0.03904847,
   0.04026014};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(10,hForwardBackwardRatio_Inc_uvutiivw_pdf_fx3028,hForwardBackwardRatio_Inc_uvutiivw_pdf_fy3028,hForwardBackwardRatio_Inc_uvutiivw_pdf_felx3028,hForwardBackwardRatio_Inc_uvutiivw_pdf_fehx3028,hForwardBackwardRatio_Inc_uvutiivw_pdf_fely3028,hForwardBackwardRatio_Inc_uvutiivw_pdf_fehy3028);
   grae->SetName("hForwardBackwardRatio_Inc_uvutiivw_pdf");
   grae->SetTitle("");
   grae->SetFillStyle(1000);
   grae->Draw("alp");
   c->Modified();
   c->cd();
   c->SetSelected(c);
}

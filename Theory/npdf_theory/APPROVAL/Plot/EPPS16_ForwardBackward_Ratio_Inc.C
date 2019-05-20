void EPPS16_ForwardBackward_Ratio_Inc()
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
   
   Double_t hForwardBackwardRatio_Inc_ltiseevo_pdf_fx3016[10] = {
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
   Double_t hForwardBackwardRatio_Inc_ltiseevo_pdf_fy3016[10] = {
   0.9895072,
   0.9731542,
   0.9517868,
   0.9371233,
   0.924061,
   0.9102143,
   0.8983585,
   0.888922,
   0.881615,
   0.8823097};
   Double_t hForwardBackwardRatio_Inc_ltiseevo_pdf_felx3016[10] = {
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
   Double_t hForwardBackwardRatio_Inc_ltiseevo_pdf_fely3016[10] = {
   0.007806329,
   0.02780043,
   0.03853858,
   0.05374167,
   0.06966656,
   0.07186573,
   0.08643008,
   0.09176217,
   0.09476115,
   0.1029079};
   Double_t hForwardBackwardRatio_Inc_ltiseevo_pdf_fehx3016[10] = {
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
   Double_t hForwardBackwardRatio_Inc_ltiseevo_pdf_fehy3016[10] = {
   0.007417356,
   0.02432385,
   0.04225177,
   0.04875558,
   0.06224935,
   0.07621699,
   0.08115126,
   0.08401089,
   0.09536141,
   0.09001049};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(10,hForwardBackwardRatio_Inc_ltiseevo_pdf_fx3016,hForwardBackwardRatio_Inc_ltiseevo_pdf_fy3016,hForwardBackwardRatio_Inc_ltiseevo_pdf_felx3016,hForwardBackwardRatio_Inc_ltiseevo_pdf_fehx3016,hForwardBackwardRatio_Inc_ltiseevo_pdf_fely3016,hForwardBackwardRatio_Inc_ltiseevo_pdf_fehy3016);
   grae->SetName("hForwardBackwardRatio_Inc_ltiseevo_pdf");
   grae->SetTitle("");
   grae->SetFillColor(1);
   grae->Draw("alp");
   c->Modified();
   c->cd();
   c->SetSelected(c);
}

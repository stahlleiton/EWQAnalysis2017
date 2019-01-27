// ROOT headers
#include "TGraphAsymmErrors.h"
// c++ headers
#include <iostream>
#include <map>

typedef std::map< std::string , TGraphAsymmErrors > GraphMap;

void etaLABtoCM(double& eta, const bool ispPb=true)
{
  const double etaLAB = eta;
  const double shift = ( ispPb ? 0.465 : -0.465 );
  double etaCM = etaLAB - shift;
  roundValue(etaCM, 4);
  eta = etaCM;
};

void unc90to68CL(double& unc)
{
  const double convFactor = (1./std::sqrt(2.))/TMath::ErfcInverse(1.-0.90);
  unc *= convFactor;
};

void graphInvEta( TGraphAsymmErrors& gr )
{
  const auto tmp = gr;
  for (int i = 0; i < tmp.GetN(); i++){
    double eta, y;
    tmp.GetPoint(i, eta, y);
    double etaErrLo = tmp.GetErrorXlow(i);
    double etaErrHi = tmp.GetErrorXhigh(i);
    double yErrLo = tmp.GetErrorYlow(i);
    double yErrHi = tmp.GetErrorYhigh(i);
    gr.SetPoint((tmp.GetN()-i-1), eta, y);
    gr.SetPointEXlow((tmp.GetN()-i-1), etaErrLo);
    gr.SetPointEXhigh((tmp.GetN()-i-1), etaErrHi);
    gr.SetPointEYlow((tmp.GetN()-i-1), yErrLo);
    gr.SetPointEYhigh((tmp.GetN()-i-1), yErrHi);
  }
};

void HIN_13007_chasym(GraphMap& graph)
{
  //
  // HIN-13007
  Double_t Wasym_fx[10] = {
    -2.2,
    -1.75,
    -1.25,
    -0.75,
    -0.25,
    0.25,
    0.75,
    1.250,
    1.75,
    2.2
  };
  for (uint i=0; i<10; i++) { etaLABtoCM(Wasym_fx[i]); }
  Double_t Wasym_fex[10] = {
    0.2,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.2
  };
  Double_t Wasym_fy[10] = {
    -0.23,
    -0.11,
    0.01,
    0.1,
    0.11,
    0.16,
    0.14,
    0.17,
    0.27,
    0.26
  };
  Double_t Wasym_fey_stat[10] = {
    0.02,
    0.02,
    0.02,
    0.01,
    0.02,
    0.01,
    0.01,
    0.02,
    0.02,
    0.02
  };
  Double_t Wasym_fey_syst[10] = {
    0.03,
    0.02,
    0.02,
    0.02,
    0.02,
    0.02,
    0.02,
    0.02,
    0.02,
    0.03
  };
  Double_t Wasym_fey[10];
  for (uint i=0; i<10; i++) {
    Wasym_fey[i] = std::sqrt( std::pow(Wasym_fey_stat[i], 2.0) + std::pow(Wasym_fey_syst[i], 2.0) );
  };
  //
  graph["HIN-13007_Stat"] = TGraphAsymmErrors(10,Wasym_fx,Wasym_fy,Wasym_fex,Wasym_fex,Wasym_fey_stat,Wasym_fey_stat);
  graph.at("HIN-13007_Stat").SetName("HIN-13007_Stat"); graph.at("HIN-13007_Stat").SetTitle("HIN-13007_Stat");
  graph["HIN-13007_Tot"] = TGraphAsymmErrors(10,Wasym_fx,Wasym_fy,Wasym_fex,Wasym_fex,Wasym_fey,Wasym_fey);
  graph.at("HIN-13007_Tot").SetName("HIN-13007_Tot"); graph.at("HIN-13007_Tot").SetTitle("HIN-13007_Tot");
  //
};


void HIN_13007_Wm(GraphMap& graph)
{
  //
  // HIN-13007
  Double_t Wm_fx[10] = {
    -2.2,
    -1.75,
    -1.25,
    -0.75,
    -0.25,
    0.25,
    0.75,
    1.250,
    1.75,
    2.2
  };
  for (uint i=0; i<10; i++) { etaLABtoCM(Wm_fx[i]); }
  Double_t Wm_fex[10] = {
    0.2,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.2
  };
  Double_t Wm_fy[10] = {
    72.1,
    79.9,
    85.4,
    81.1,
    80.8,
    78.0,
    76.5,
    70.1,
    59.8,
    63.7
  };
  Double_t Wm_fey_stat[10] = {
    2.2,
    2.1,
    2.0,
    1.8,
    1.9,
    1.8,
    1.8,
    1.8,
    1.7,
    2.1
  };
  Double_t Wm_fey_syst[10] = {
    3.7,
    3.3,
    2.7,
    2.1,
    2.2,
    2.3,
    2.4,
    2.5,
    2.6,
    3.9
  };
  Double_t Wm_fey[10];
  for (uint i=0; i<10; i++) {
    Wm_fey[i] = std::sqrt( std::pow(Wm_fey_stat[i], 2.0) + std::pow(Wm_fey_syst[i], 2.0) );
  };
  //
  graph["HIN-13007_Stat"] = TGraphAsymmErrors(10,Wm_fx,Wm_fy,Wm_fex,Wm_fex,Wm_fey_stat,Wm_fey_stat);
  graph.at("HIN-13007_Stat").SetName("HIN-13007_Stat"); graph.at("HIN-13007_Stat").SetTitle("HIN-13007_Stat");
  graph["HIN-13007_Tot"] = TGraphAsymmErrors(10,Wm_fx,Wm_fy,Wm_fex,Wm_fex,Wm_fey,Wm_fey);
  graph.at("HIN-13007_Tot").SetName("HIN-13007_Tot"); graph.at("HIN-13007_Tot").SetTitle("HIN-13007_Tot");
  //
};


void HIN_13007_Wp(GraphMap& graph)
{
  //
  // HIN-13007
  Double_t Wp_fx[10] = {
    -2.2,
    -1.75,
    -1.25,
    -0.75,
    -0.25,
    0.25,
    0.75,
    1.250,
    1.75,
    2.2
  };
  for (uint i=0; i<10; i++) { etaLABtoCM(Wp_fx[i]); }
  Double_t Wp_fex[10] = {
    0.2,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.2
  };
  Double_t Wp_fy[10] = {
    44.5,
    62.9,
    85.9,
    98.6,
    99.7,
    105.3,
    101.8,
    100.2,
    102.3,
    108.1
  };
  Double_t Wp_fey_stat[10] = {
    1.7,
    1.8,
    2.0,
    2.1,
    2.1,
    2.1,
    2.1,
    2.2,
    2.3,
    2.8
  };
  Double_t Wp_fey_syst[10] = {
    2.3,
    2.2,
    2.7,
    2.3,
    2.7,
    2.8,
    2.5,
    3.1,
    4.2,
    6.0
  };
  Double_t Wp_fey[10];
  for (uint i=0; i<10; i++) {
    Wp_fey[i] = std::sqrt( std::pow(Wp_fey_stat[i], 2.0) + std::pow(Wp_fey_syst[i], 2.0) );
  }
  //
  graph["HIN-13007_Stat"] = TGraphAsymmErrors(10,Wp_fx,Wp_fy,Wp_fex,Wp_fex,Wp_fey_stat,Wp_fey_stat);
  graph.at("HIN-13007_Stat").SetName("HIN-13007_Stat"); graph.at("HIN-13007_Stat").SetTitle("HIN-13007_Stat");
  graph["HIN-13007_Tot"] = TGraphAsymmErrors(10,Wp_fx,Wp_fy,Wp_fex,Wp_fex,Wp_fey,Wp_fey);
  graph.at("HIN-13007_Tot").SetName("HIN-13007_Tot"); graph.at("HIN-13007_Tot").SetTitle("HIN-13007_Tot");
  //
};


void HIN_13007_Theory_chasym(GraphMap& graph)
{
  //
  // HIN-13007 CT10
  Double_t Graph0_fx3001[10] = {
    2.2,
    1.75,
    1.25,
    0.75,
    0.25,
    -0.25,
    -0.75,
    -1.25,
    -1.75,
    -2.2
  };
  for (uint i=0; i<10; i++) { etaLABtoCM(Graph0_fx3001[i]); }
  Double_t Graph0_fy3001[10] = {
    0.263377,
    0.234428,
    0.19231,
    0.156091,
    0.13209,
    0.115851,
    0.0964275,
    0.0539073,
    -0.0465022,
    -0.215376
  };
  Double_t Graph0_felx3001[10] = {
    0.2,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.2
  };
  Double_t Graph0_fely3001[10] = {
    0.01276083,
    0.01766791,
    0.02757422,
    0.01753697,
    0.01045832,
    0.01296183,
    0.0161304,
    0.01688416,
    0.009157931,
    0.01377135
  };
  for (uint i=0; i<10; i++) { unc90to68CL(Graph0_fely3001[i]); }
  Double_t Graph0_fehx3001[10] = {
    0.2,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.2
  };
  Double_t Graph0_fehy3001[10] = {
    0.02858905,
    0.01664689,
    0.004525079,
    0.01285556,
    0.01411928,
    0.01166375,
    0.006458504,
    0.007979965,
    0.02122952,
    0.0218889
  };
  for (uint i=0; i<10; i++) { unc90to68CL(Graph0_fehy3001[i]); }
  //
  graph["CT10"] = TGraphAsymmErrors(10,Graph0_fx3001,Graph0_fy3001,Graph0_felx3001,Graph0_fehx3001,Graph0_fely3001,Graph0_fehy3001);
  graph.at("CT10").SetName("HIN-13007_CT10"); graph.at("CT10").SetTitle("HIN-13007_CT10");
  graphInvEta(graph.at("CT10"));
  //
  // HIN-13007 EPS09
  Double_t Graph0_fx3002[10] = {
    2.2,
    1.75,
    1.25,
    0.75,
    0.25,
    -0.25,
    -0.75,
    -1.25,
    -1.75,
    -2.2
  };
  for (uint i=0; i<10; i++) { etaLABtoCM(Graph0_fx3002[i]); }
  Double_t Graph0_fy3002[10] = {
    0.263377,
    0.234428,
    0.19231,
    0.156091,
    0.13209,
    0.115851,
    0.0964275,
    0.0539073,
    -0.0465022,
    -0.215376
  };
  Double_t Graph0_felx3002[10] = {
    0.2,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.2
  };
  Double_t Graph0_fely3002[10] = {
    0.01276083,
    0.01766791,
    0.02757422,
    0.01753697,
    0.01045832,
    0.01296183,
    0.0161304,
    0.01688416,
    0.009157931,
    0.01377135
  };
  for (uint i=0; i<10; i++) { unc90to68CL(Graph0_fely3002[i]); }
  Double_t Graph0_fehx3002[10] = {
    0.2,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.2
  };
  Double_t Graph0_fehy3002[10] = {
    0.02858905,
    0.01664689,
    0.004525079,
    0.01285556,
    0.01411928,
    0.01166375,
    0.006458504,
    0.007979965,
    0.02122952,
    0.0218889
  };
  for (uint i=0; i<10; i++) { unc90to68CL(Graph0_fehy3002[i]); }
  //
  graph["EPS09"] = TGraphAsymmErrors(10,Graph0_fx3002,Graph0_fy3002,Graph0_felx3002,Graph0_fehx3002,Graph0_fely3002,Graph0_fehy3002);
  graph.at("EPS09").SetName("HIN-13007_EPS09"); graph.at("EPS09").SetTitle("HIN-13007_EPS09");
  graphInvEta(graph.at("EPS09"));
  //
  // HIN-13007 EPPS16
  Double_t hChargeAsymmetry_usfbjrib_pdf_fx3001[10] = {
    -2.2,
    -1.75,
    -1.25,
    -0.75,
    -0.25,
    0.25,
    0.75,
    1.25,
    1.75,
    2.2
  };
  for (uint i=0; i<10; i++) { etaLABtoCM(hChargeAsymmetry_usfbjrib_pdf_fx3001[i]); }
  Double_t hChargeAsymmetry_usfbjrib_pdf_fy3001[10] = {
    -0.2336219,
    -0.06014597,
    0.04517266,
    0.09224602,
    0.1132514,
    0.1284726,
    0.1503874,
    0.1867246,
    0.2307891,
    0.2603281
  };
  Double_t hChargeAsymmetry_usfbjrib_pdf_felx3001[10] = {
    0.2,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.2
  };
  Double_t hChargeAsymmetry_usfbjrib_pdf_fely3001[10] = {
    0.0244727,
    0.02410125,
    0.01962496,
    0.01396378,
    0.008932523,
    0.005300171,
    0.006142093,
    0.009793986,
    0.01204817,
    0.01241369
  };
  Double_t hChargeAsymmetry_usfbjrib_pdf_fehx3001[10] = {
    0.2,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.2
  };
  Double_t hChargeAsymmetry_usfbjrib_pdf_fehy3001[10] = {
    0.02524164,
    0.02342904,
    0.02015705,
    0.01373849,
    0.008144228,
    0.007470627,
    0.009097115,
    0.01005905,
    0.01298757,
    0.01487545
  };
  //
  graph["EPPS16"] = TGraphAsymmErrors(10,hChargeAsymmetry_usfbjrib_pdf_fx3001,hChargeAsymmetry_usfbjrib_pdf_fy3001,hChargeAsymmetry_usfbjrib_pdf_felx3001,hChargeAsymmetry_usfbjrib_pdf_fehx3001,hChargeAsymmetry_usfbjrib_pdf_fely3001,hChargeAsymmetry_usfbjrib_pdf_fehy3001);
  graph.at("EPPS16").SetName("HIN-13007_EPPS16_ChgAsy"); graph.at("EPPS16").SetTitle("HIN-13007_EPPS16_ChgAsy");
  //
};


void HIN_13007_Theory_Wm(GraphMap& graph)
{
  //
  // HIN-13007 CT10
  Double_t Graph0_fx3003[10] = {
    2.2,
    1.75,
    1.25,
    0.75,
    0.25,
    -0.25,
    -0.75,
    -1.25,
    -1.75,
    -2.2
  };
  for (uint i=0; i<10; i++) { etaLABtoCM(Graph0_fx3003[i]); }
  Double_t Graph0_fy3003[10] = {
    59.32539,
    64.47146,
    69.07048,
    72.39233,
    74.34135,
    75.07561,
    74.77608,
    73.3034,
    70.13963,
    64.15945
  };
  Double_t Graph0_felx3003[10] = {
    0.2,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.2
  };
  Double_t Graph0_fely3003[10] = {
    4.049307,
    4.291061,
    4.910586,
    5.11398,
    5.163876,
    5.085138,
    4.844703,
    4.758015,
    4.724943,
    4.441151
  };
  for (uint i=0; i<10; i++) { unc90to68CL(Graph0_fely3003[i]); }
  Double_t Graph0_fehx3003[10] = {
    0.2,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.2
  };
  Double_t Graph0_fehy3003[10] = {
    3.565033,
    3.783372,
    4.186292,
    4.760494,
    4.941513,
    4.68498,
    4.43467,
    4.34438,
    4.178857,
    3.863407
  };
  for (uint i=0; i<10; i++) { unc90to68CL(Graph0_fehy3003[i]); }
  //
  graph["CT10"] = TGraphAsymmErrors(10,Graph0_fx3003,Graph0_fy3003,Graph0_felx3003,Graph0_fehx3003,Graph0_fely3003,Graph0_fehy3003);
  graph.at("CT10").SetName("HIN-13007_CT10"); graph.at("CT10").SetTitle("HIN-13007_CT10");
  graphInvEta(graph.at("CT10"));
  //
  // HIN-13007 EPS09
  Double_t Graph0_fx3004[10] = {
    2.2,
    1.75,
    1.25,
    0.75,
    0.25,
    -0.25,
    -0.75,
    -1.25,
    -1.75,
    -2.2
  };
  for (uint i=0; i<10; i++) { etaLABtoCM(Graph0_fx3004[i]); }
  Double_t Graph0_fy3004[10] = {
    53.40554,
    59.08411,
    64.79595,
    69.69657,
    73.55716,
    76.08652,
    76.94142,
    75.91387,
    72.39233,
    65.239
  };
  Double_t Graph0_felx3004[10] = {
    0.2,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.2
  };
  Double_t Graph0_fely3004[10] = {
    4.840949,
    4.675164,
    4.500348,
    4.806465,
    5.240032,
    5.356636,
    4.781826,
    4.688571,
    6.164546,
    4.613319
  };
  for (uint i=0; i<10; i++) { unc90to68CL(Graph0_fely3004[i]); }
  Double_t Graph0_fehx3004[10] = {
    0.2,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.2
  };
  Double_t Graph0_fehy3004[10] = {
    3.981867,
    4.47893,
    5.20512,
    4.926794,
    5.167483,
    4.712881,
    4.983093,
    4.985132,
    3.578494,
    4.19526
  };
  for (uint i=0; i<10; i++) { unc90to68CL(Graph0_fehy3004[i]); }
  //
  graph["EPS09"] = TGraphAsymmErrors(10,Graph0_fx3004,Graph0_fy3004,Graph0_felx3004,Graph0_fehx3004,Graph0_fely3004,Graph0_fehy3004);
  graph.at("EPS09").SetName("HIN-13007_EPS09"); graph.at("EPS09").SetTitle("HIN-13007_EPS09");
  graphInvEta(graph.at("EPS09"));
  //
  // HIN-13007 EPPS16
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
    2.2
  };
  for (uint i=0; i<10; i++) { etaLABtoCM(hCrossSection_okmjehgq_pdf_fx3002[i]); }
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
    54.07461
  };
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
    0.2
  };
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
    4.08787
  };
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
    0.2
  };
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
    4.067623
  };
  //
  graph["EPPS16"] = TGraphAsymmErrors(10,hCrossSection_okmjehgq_pdf_fx3002,hCrossSection_okmjehgq_pdf_fy3002,hCrossSection_okmjehgq_pdf_felx3002,hCrossSection_okmjehgq_pdf_fehx3002,hCrossSection_okmjehgq_pdf_fely3002,hCrossSection_okmjehgq_pdf_fehy3002);
  graph.at("EPPS16").SetName("HIN-13007_EPPS16_XSec"); graph.at("EPPS16").SetTitle("HIN-13007_EPPS16_XSec");
  //
};


void HIN_13007_Theory_Wp(GraphMap& graph)
{
  //
  // HIN-13007 CT10
  Double_t Graph0_fx3005[10] = {
    2.2,
    1.75,
    1.25,
    0.75,
    0.25,
    -0.25,
    -0.75,
    -1.25,
    -1.75,
    -2.2
  };
  for (uint i=0; i<10; i++) { etaLABtoCM(Graph0_fx3005[i]); }
  Double_t Graph0_fy3005[10] = {
    103.56,
    105.744,
    103.4081,
    100.0093,
    97.03688,
    94.36609,
    90.5471,
    82.43278,
    65.82142,
    43.43165
  };
  Double_t Graph0_felx3005[10] = {
    0.2,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.2
  };
  Double_t Graph0_fely3005[10] = {
    7.299733,
    7.206486,
    6.824632,
    6.594841,
    6.406375,
    6.079494,
    5.746822,
    5.225105,
    4.165823,
    2.787035
  };
  for (uint i=0; i<10; i++) { unc90to68CL(Graph0_fely3005[i]); }
  Double_t Graph0_fehx3005[10] = {
    0.2,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.2
  };
  Double_t Graph0_fehy3005[10] = {
    5.988963,
    6.218694,
    6.307709,
    6.368316,
    5.847823,
    5.655323,
    5.083818,
    4.561256,
    3.722445,
    2.542957
  };
  for (uint i=0; i<10; i++) { unc90to68CL(Graph0_fehy3005[i]); }
  //
  graph["CT10"] = TGraphAsymmErrors(10,Graph0_fx3005,Graph0_fy3005,Graph0_felx3005,Graph0_fehx3005,Graph0_fely3005,Graph0_fehy3005);
  graph.at("CT10").SetName("HIN-13007_CT10"); graph.at("CT10").SetTitle("HIN-13007_CT10");
  graphInvEta(graph.at("CT10"));
  //
  // HIN-13007 EPS09
  Double_t Graph0_fx3006[10] = {
    2.2,
    1.75,
    1.25,
    0.75,
    0.25,
    -0.25,
    -0.75,
    -1.25,
    -1.75,
    -2.2
  };
  for (uint i=0; i<10; i++) { etaLABtoCM(Graph0_fx3006[i]); }
  Double_t Graph0_fy3006[10] = {
    91.59545,
    95.26883,
    95.65156,
    95.47892,
    95.94693,
    96.02597,
    93.3635,
    84.56484,
    65.95871,
    42.11706
  };
  Double_t Graph0_felx3006[10] = {
    0.2,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.2
  };
  Double_t Graph0_fely3006[10] = {
    8.093772,
    8.102706,
    8.186826,
    6.574564,
    6.279891,
    6.620883,
    6.246574,
    5.584516,
    4.554484,
    2.654461
  };
  for (uint i=0; i<10; i++) { unc90to68CL(Graph0_fely3006[i]); }
  Double_t Graph0_fehx3006[10] = {
    0.2,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.2
  };
  Double_t Graph0_fehy3006[10] = {
    9.277223,
    7.867945,
    5.692544,
    6.530371,
    6.577697,
    5.769801,
    5.218286,
    4.471689,
    3.667556,
    2.965252
  };
  for (uint i=0; i<10; i++) { unc90to68CL(Graph0_fehy3006[i]); }
  //
  graph["EPS09"] = TGraphAsymmErrors(10,Graph0_fx3006,Graph0_fy3006,Graph0_felx3006,Graph0_fehx3006,Graph0_fely3006,Graph0_fehy3006);
  graph.at("EPS09").SetName("HIN-13007_EPS09"); graph.at("EPS09").SetTitle("HIN-13007_EPS09");
  graphInvEta(graph.at("EPS09"));
  //
  // HIN-13007 EPPS16
  Double_t hCrossSection_ppeguzkf_pdf_fx3003[10] = {
    -2.2,
    -1.75,
    -1.25,
    -0.75,
    -0.25,
    0.25,
    0.75,
    1.25,
    1.75,
    2.2
  };
  for (uint i=0; i<10; i++) { etaLABtoCM(hCrossSection_ppeguzkf_pdf_fx3003[i]); }
  Double_t hCrossSection_ppeguzkf_pdf_fy3003[10] = {
    42.18692,
    66.15761,
    85.01131,
    94.07329,
    96.98965,
    97.20647,
    96.67413,
    96.76384,
    96.13419,
    92.13782
  };
  Double_t hCrossSection_ppeguzkf_pdf_felx3003[10] = {
    0.2,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.2
  };
  Double_t hCrossSection_ppeguzkf_pdf_fely3003[10] = {
    1.653754,
    2.493004,
    3.125351,
    3.250007,
    3.354755,
    3.279073,
    3.52704,
    4.943047,
    6.65084,
    7.836764
  };
  Double_t hCrossSection_ppeguzkf_pdf_fehx3003[10] = {
    0.2,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.25,
    0.2
  };
  Double_t hCrossSection_ppeguzkf_pdf_fehy3003[10] = {
    1.768937,
    2.495407,
    2.909781,
    3.005679,
    3.077289,
    3.04223,
    3.751912,
    5.106617,
    7.007804,
    8.226671
  };
  //
  graph["EPPS16"] = TGraphAsymmErrors(10,hCrossSection_ppeguzkf_pdf_fx3003,hCrossSection_ppeguzkf_pdf_fy3003,hCrossSection_ppeguzkf_pdf_felx3003,hCrossSection_ppeguzkf_pdf_fehx3003,hCrossSection_ppeguzkf_pdf_fely3003,hCrossSection_ppeguzkf_pdf_fehy3003);
  graph.at("EPPS16").SetName("HIN-13007_EPPS16_XSec"); graph.at("EPPS16").SetTitle("HIN-13007_EPPS16_XSec");
  //
};

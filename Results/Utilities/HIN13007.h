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

void eta5TeVto8TeV(double& eta)
{
  const double eta5TeV = eta;
  const double shift = TMath::Log(8.16/5.02);
  double eta8TeV = 0.0;
  if (eta5TeV>=0.0) { eta8TeV = eta5TeV + shift; }
  if (eta5TeV <0.0) { eta8TeV = eta5TeV - shift; }
  eta = eta8TeV;
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
  for (uint i=0; i<10; i++) { eta5TeVto8TeV(Wasym_fx[i]); }
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
  graph["HIN-13007"] = TGraphAsymmErrors(10,Wasym_fx,Wasym_fy,Wasym_fex,Wasym_fex,Wasym_fey,Wasym_fey);
  graph.at("HIN-13007").SetName("HIN-13007"); graph.at("HIN-13007").SetTitle("HIN-13007");
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
  graph["HIN-13007"] = TGraphAsymmErrors(10,Wm_fx,Wm_fy,Wm_fex,Wm_fex,Wm_fey,Wm_fey);
  graph.at("HIN-13007").SetName("HIN-13007"); graph.at("HIN-13007").SetTitle("HIN-13007");
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
  graph["HIN-13007"] = TGraphAsymmErrors(10,Wp_fx,Wp_fy,Wp_fex,Wp_fex,Wp_fey,Wp_fey);
  graph.at("HIN-13007").SetName("HIN-13007"); graph.at("HIN-13007").SetTitle("HIN-13007");
  //
};

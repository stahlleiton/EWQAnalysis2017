// Auxiliary Headers
#include "Utilities/resultsUtils.h"
#include "../Utilities/EVENTUTILS.h"
// ROOT headers
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TEfficiency.h"
#include "TVectorD.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TPad.h"
#include "TLine.h"
// RooFit headers
#include "RooWorkspace.h"
#include "RooStringVar.h"
#include "RooRealVar.h"
// c++ headers
#include <iostream>
#include <map>
#include <array>
// CMS headers
#include "../Utilities/CMS/tdrstyle.C"
#include "../Utilities/CMS/CMS_lumi.C"


// ------------------ TYPE -------------------------------
typedef std::pair< std::string , uint              > StrIntPair;
typedef std::map < std::string , StrIntPair        > StrPairMap;
typedef std::map < std::string , TGraphAsymmErrors > GraphMap;
using Unc1DVec_t   =  std::map< std::string , TVectorD >;
//
using Unc1DMap2_t  =  std::map< std::string , std::map< std::string , std::map< std::string , Unc1DVec_t > > >;
using Unc1DMap_t   =  std::map< std::string , std::map< std::string , Unc1DMap2_t > >;
using EffVec_t     =  std::map< std::string , std::vector< TEfficiency > >;
using EffMap2_t    =  std::map< std::string , std::map< std::string , std::map< std::string , EffVec_t > > >;
using EffMap_t     =  std::map< std::string , std::map< std::string , EffMap2_t > >;
using BinMap_t     =  std::map< std::string , std::vector< double > >;
using BinMapMap_t  =  std::map< std::string , BinMap_t >;
using KeyPVec_t    =  std::vector< TKey* >;


const std::vector< std::string > CHG_  = { "Plus" , "Minus" };
/////////////////////
// OTHER FUNCTIONS //
/////////////////////

void drawCompareGraph(GraphMap& grMap, GraphMap& grCmpMap, const std::string& outDir, const std::string& chg, const StringMap& legLbl,
                      const std::string col="PA", const std::string whatDo="pull", const bool doCorrectYield=false);
void fillCompareGraph(TGraphAsymmErrors& graph, const BinMapMap& var, const std::string graphLbl);
bool extractInfo(BinQuadMap& inputVar, const std::string& workDirPath);
bool makePullGraph(GraphMap& grCmpMap, const GraphMap& grMap, const GraphMap& grMapRaw, const StrPairMap& cmpWS, const std::string whatDo="pull");
bool printCompareTable(GraphMap& grMap, const std::string& outDir, const std::string& chg, const std::string& refLbl, const std::string col, const bool doCorrectYield);
//
bool         getObjectsFromFile  ( EffMap_t& eff , Unc1DMap_t& unc , const std::string& filePath );
KeyPVec_t    getK                ( TList* list );

/////////////////////


void compareWS(const std::string mode = "Runs", const bool doCorrectYield=false)
{
  //
  StringMap refWS;
  StrPairMap cmpWS;
  if (mode=="Runs") {
    refWS = StringMap({
	{ "PA" , "/home/llr/cms/stahl/ElectroWeakAnalysis/NOMINAL/BOSONPTCORR/EWQAnalysis2017/Fitter/Output/NominalCM/METPF_RAW/DATA/W/PA/result/" }
      });
    cmpWS = StrPairMap({
	{ "pPb"  , { "/home/llr/cms/stahl/ElectroWeakAnalysis/NOMINAL/BOSONPTCORR/EWQAnalysis2017/Fitter/Output/NominalCM/METPF_RAW/DATA/W/pPb/result/" , 1 } },
        { "Pbp"  , { "/home/llr/cms/stahl/ElectroWeakAnalysis/NOMINAL/BOSONPTCORR/EWQAnalysis2017/Fitter/Output/NominalCM/METPF_RAW/DATA/W/Pbp/result/" , 1 } }
      });
  }
  else if (mode=="Runs2") {
    refWS = StringMap({
	{ "PA" , "/home/llr/cms/stahl/ElectroWeakAnalysis/EWQAnalysis2017/Fitter/Output/NominalCM/METPF_RAW/DATA/W/PA/result/" }
      });
    cmpWS = StrPairMap({
	{ "pPb"  , { "/home/llr/cms/stahl/ElectroWeakAnalysis/EWQAnalysis2017/Fitter/Output/NominalCM/METPF_RAW/DATA/W/pPb/result/" , 1 } },
        { "Pbp"  , { "/home/llr/cms/stahl/ElectroWeakAnalysis/EWQAnalysis2017/Fitter/Output/NominalCM/METPF_RAW/DATA/W/Pbp/result/" , 1 } }
      });
  }
  else if (mode=="PU") {
    refWS = StringMap({
	{ "WithPU" , "/home/llr/cms/stahl/ElectroWeakAnalysis/NOMINAL/BOSONPTCORR/EWQAnalysis2017/Fitter/Output/NominalCM/METPF_RAW/DATA/W/PA/result/" }
      });
    cmpWS = StrPairMap({
	{ "NoPU"  , { "/home/llr/cms/stahl/ElectroWeakAnalysis/TESTDS/FILTER/EWQAnalysis2017/Fitter/Output/NominalCM/METPF_RAW/DATA/W/PA/result/" , 1 } }
      });
  }
  else if (mode=="BKG") {
    refWS = StringMap({
	{ "REF" , "/home/llr/cms/stahl/ElectroWeakAnalysis/NOMINAL/BOSONPTCORR/EWQAnalysis2017/Fitter/Output/NominalCM/METPF_RAW/DATA/W/PA/result/" }
      });
    cmpWS = StrPairMap({
	//{ "VAR"  , { "/home/llr/cms/stahl/ElectroWeakAnalysis/EWQAnalysis2017/Fitter/Output/NominalCM/METPF_RAW/DATA/W/PA/result/" , 1 } }
	{ "VAR"  , { "/home/llr/cms/stahl/ElectroWeakAnalysis/NOMINAL/BOSONPTCORR/EWQAnalysis2017/Fitter/Output/NominalCM_TEST/METPF_RAW/DATA/W/PA/result/" , 1 } }
      });
  }
  /*
  else if (mode=="BKG") {
    refWS = StringMap({
	{ "OLD" , "/home/llr/cms/stahl/ElectroWeakAnalysis/EWQAnalysis2017/Fitter/Output/PreApproval_HomeWork/NominalCM/METPF_RAW/DATA/W/PA/result/" }
      });
    cmpWS = StrPairMap({
	{ "NEW"  , { "/home/llr/cms/stahl/ElectroWeakAnalysis/EWQAnalysis2017/Fitter/Output/NominalCM/METPF_RAW/DATA/W/PA/result/" , 1 } }
      });
  }
  */
  else { std::cout << "[ERROR] Mode (" << mode << ") is not valid!" << std::endl; return; }
  //
  StringMap legLabel;
  //
  // Extract the information
  BinPentaMap refInputVar , cmpInputVar;
  std::map< std::string , double > scaleMap;
  for (const auto& lblWS : refWS) { extractInfo(refInputVar[lblWS.first], lblWS.second); legLabel["Ref"] = Form("Ref (%s)", lblWS.first.c_str()); }
  for (const auto& lblWS : cmpWS) {
    auto& scale = scaleMap[lblWS.first];
    scale = 1.0;
    if (lblWS.first=="pPb")  { scale = ( (PA::LUMI::Data_pPb + PA::LUMI::Data_Pbp) / (PA::LUMI::Data_pPb) ); legLabel["pPb"] = Form("pPb (x%.2f)", scale); }
    if (lblWS.first=="Pbp")  { scale = ( (PA::LUMI::Data_pPb + PA::LUMI::Data_Pbp) / (PA::LUMI::Data_Pbp) ); legLabel["Pbp"] = Form("Pbp (x%.2f, inverse)", scale); }
    if (lblWS.first=="NoPU") { scale = 1.18; legLabel["NoPU"] = Form("NoPU (x%.2f)", scale); }
    if (lblWS.first=="Nom")  { legLabel["Nom"] = "Nominal"; }
    if (lblWS.first=="AllBkg") { legLabel["AllBkg"] = "Nominal + DY->Tau + VV"; }
    if (lblWS.first=="NEW") { legLabel["NEW"] = "NEW Results"; }
    if (lblWS.first=="OLD") { legLabel["OLD"] = "OLD Results"; }
    if (lblWS.first=="REF") { legLabel["REF"] = "Nominal"; }
    if (lblWS.first=="VAR") { legLabel["VAR"] = "With Boson PT Corr"; }
    extractInfo(cmpInputVar[lblWS.first], lblWS.second.first);
  }
  //
  // Extract the efficiency
  const std::string CWD = getcwd(NULL, 0);
  std::string preCWD = CWD; preCWD.erase(preCWD.find_last_of("/"), 10);
  const std::string effDir = Form("%s/Efficiency/Output/%s/", preCWD.c_str(), "NominalCM_WithBosonPT_WithHF");
  // Declare the efficiencies
  EffMap_t eff1D;
  Unc1DMap_t unc1D;
  // Extract information from the input file
  const std::string effFileName = effDir + "efficiencyTnP.root";
  if (!getObjectsFromFile(eff1D, unc1D, effFileName)) { return; }
  //
  // Set Style
  setStyle();
  //
  if (mode=="Runs2") {
    //
    BinPentaMap tmpVar;
    for (const auto& ch : cmpInputVar.at("pPb")) {
      //
      const auto& cmpVar_pPb = cmpInputVar.at("pPb").at(ch.first).at("N_WToMu");
      const auto& cmpVar_Pbp = cmpInputVar.at("Pbp").at(ch.first).at("N_WToMu");
      const auto& cmpVar_PA  = refInputVar.at("PA").at(ch.first).at("N_WToMu");
      //
      for (const auto& b : cmpVar_pPb.at("Val")) {
        double eff_pPb = 1.0 , eff_Pbp = 1.0 , eff_PA = 1.0;
        const double etaV = ( (b.first.etabin().high() + b.first.etabin().low()) / 2.0 ); // Mean value of eta bin
        if (doCorrectYield) {
          const std::string chg = (ch.first=="Pl" ? "Plus" : "Minus");
          const auto& effP_pPb = eff1D.at("EtaCM").at("MC_WToMuNu").at("pPb").at(chg).at("Total").at("TnP_Nominal")[0];
          const auto& effP_Pbp = eff1D.at("EtaCM").at("MC_WToMuNu").at("Pbp").at(chg).at("Total").at("TnP_Nominal")[0];
          const auto& effP_PA  = eff1D.at("EtaCM").at("MC_WToMuNu").at("PA" ).at(chg).at("Total").at("TnP_Nominal")[0];
          eff_pPb = effP_pPb.GetEfficiency(effP_pPb.FindFixBin(etaV));
          eff_Pbp = effP_Pbp.GetEfficiency(effP_Pbp.FindFixBin(-etaV));
          eff_PA  = effP_PA.GetEfficiency(effP_PA.FindFixBin(etaV));
        }
        if (ch.first=="Pl") {
          std::cout << "eta " << etaV << "  Yields: " << "pPb " << (cmpVar_pPb.at("Val").at(b.first)/eff_pPb) << "   " << "Pbp " << (cmpVar_Pbp.at("Val").at(b.first)/eff_Pbp) << "   PA " << (cmpVar_PA.at("Val").at(b.first)/eff_PA) << std::endl;
          std::cout << "eta " << etaV << "  Yields (scaled): " << "pPb " << (cmpVar_pPb.at("Val").at(b.first)*1.57/eff_pPb) << "   " << "Pbp " << (cmpVar_Pbp.at("Val").at(b.first)*2.77/eff_Pbp) << "   PA " << (cmpVar_PA.at("Val").at(b.first)/eff_PA) << std::endl;
        }
        tmpVar["SUM"][ch.first]["N_WToMu"]["Val"][b.first] = ( (cmpVar_pPb.at("Val").at(b.first)/eff_pPb) + (cmpVar_Pbp.at("Val").at(b.first)/eff_Pbp) )*eff_PA;
        tmpVar["SUM"][ch.first]["N_WToMu"]["Err_High"][b.first] = std::sqrt( std::pow(cmpVar_pPb.at("Err_High").at(b.first)/eff_pPb,2.0) + std::pow(cmpVar_Pbp.at("Err_High").at(b.first)/eff_Pbp,2.0) )*eff_PA;
        tmpVar["SUM"][ch.first]["N_WToMu"]["Err_Low" ][b.first] = std::sqrt( std::pow(cmpVar_pPb.at("Err_Low" ).at(b.first)/eff_pPb,2.0) + std::pow(cmpVar_Pbp.at("Err_Low" ).at(b.first)/eff_Pbp,2.0) )*eff_PA;
      }
    }
    //
    cmpInputVar.clear(); cmpInputVar = tmpVar;
    cmpWS.clear(); cmpWS["SUM"] = std::make_pair("", 2);
    scaleMap.clear(); scaleMap["SUM"] = 1.0;
    StringMap tmpLbl; tmpLbl["Ref"] = legLabel.at("Ref"); tmpLbl["SUM"] = "pPb+Pbp";
    legLabel.clear(); legLabel = tmpLbl;
  }
  //
  //
  // Create the output plots
  for (const auto& ref : refInputVar) {
    for (const auto& ch : ref.second) {
      // Create the Graphs
      GraphMap grMap, grMapRaw;
      //
      const auto& refVarRaw = ch.second.at("N_WToMu");
      fillCompareGraph(grMapRaw["Ref"], refVarRaw, Form("%s_%s", ref.first.c_str(), ch.first.c_str()));
      //
      auto refVar = ch.second.at("N_WToMu");
      if (doCorrectYield) {
        const std::string col = (ref.first=="pPb" ? "pPb" : (ref.first=="Pbp" ? "Pbp" : "PA"));
        const std::string chg = (ch.first=="Pl" ? "Plus" : "Minus");
        const auto& eff = eff1D.at("EtaCM").at("MC_WToMuNu").at(col).at(chg).at("Total").at("TnP_Nominal")[0];
        for (const auto& b : refVar.at("Val")) {
          const double etaV = (col=="Pbp" ? -1.0 : 1.0) * ( (b.first.etabin().high() + b.first.etabin().low()) / 2.0 ); // Mean value of eta bin
          const auto& effVal = eff.GetEfficiency(eff.FindFixBin(etaV));
          for (auto& v : refVar) { v.second.at(b.first) /= effVal; }
        }
      }
      for (const auto& b : refVar.at("Val")) {
        const double etaW = ( (b.first.etabin().high() - b.first.etabin().low()) / 2.0 )/0.1; // Width of eta bin
        for (auto& v : refVar) { v.second.at(b.first) /= etaW; }
      }
      if (scaleMap.count(ref.first)>0) {
        for (auto& v : refVar) { for (auto& b : v.second) { b.second *= scaleMap.at(ref.first); } }
      }
      fillCompareGraph(grMap["Ref"], refVar, Form("%s_%s", ref.first.c_str(), ch.first.c_str()));
      //
      for (const auto& cmp : cmpInputVar) {
        //
	const auto& cmpVarRaw = cmp.second.at(ch.first).at("N_WToMu");
	fillCompareGraph(grMapRaw[cmp.first], cmpVarRaw, Form("%s_%s", cmp.first.c_str(), ch.first.c_str()));
        //
	auto cmpVar = cmp.second.at(ch.first).at("N_WToMu");
        if (doCorrectYield) {
          const std::string col = (cmp.first=="pPb" ? "pPb" : (cmp.first=="Pbp" ? "Pbp" : "PA"));
          std::string chg = (ch.first=="Pl" ? "Plus" : "Minus");
          const auto& eff = eff1D.at("EtaCM").at("MC_WToMuNu").at(col).at(chg).at("Total").at("TnP_Nominal")[0];
          for (const auto& b : cmpVar.at("Val")) {
            const double etaV = (col=="Pbp" ? -1.0 : 1.0) * ( (b.first.etabin().high() + b.first.etabin().low()) / 2.0 ); // Mean value of eta bin
            const auto& effVal = eff.GetEfficiency(eff.FindFixBin(etaV));
            for (auto& v : cmpVar) { v.second.at(b.first) /= effVal; }
          }
        }
        for (const auto& b : cmpVar.at("Val")) {
          const double etaW = ( (b.first.etabin().high() - b.first.etabin().low()) / 2.0 )/0.1; // Width of eta bin
          for (auto& v : cmpVar) { v.second.at(b.first) /= etaW; }
        }
        if (scaleMap.count(cmp.first)>0) {
          for (auto& v : cmpVar) { for (auto& b : v.second) { b.second *= scaleMap.at(cmp.first); } }
        }
	fillCompareGraph(grMap[cmp.first], cmpVar, Form("%s_%s", cmp.first.c_str(), ch.first.c_str()));
      }
      // Set Output Directory
      const std::string CWD = getcwd(NULL, 0);
      const std::string outDir = CWD + "/Comparison/" + mode;
      // Draw the Graphs      
      GraphMap grCmpMap;
      std::string whatDo = "pull";
      if (mode=="BKG") { whatDo = "ratio"; }
      if (!makePullGraph(grCmpMap, grMap, grMapRaw, cmpWS, whatDo)) { return; }
      drawCompareGraph(grMap, grCmpMap, outDir, ch.first, legLabel, "PA", whatDo, doCorrectYield);
      printCompareTable(grMap, outDir, ch.first, ref.first, "PA", doCorrectYield);
    }
  }  
};




bool makePullGraph(GraphMap& grCmpMap, const GraphMap& grMap, const GraphMap& grMapRaw, const StrPairMap& cmpWS, const std::string whatDo)
{
  //
  for (const auto& gr : grMap) {
    if (gr.first=="Ref") continue;
    //
    auto& grCmp = grCmpMap[gr.first];
    // Set Bin Size
    const uint nBins = grMap.at("Ref").GetN();
    grCmp.Set(nBins);
    // Set Graph Name
    const std::string name = Form("gr_WToMu_Comparison_%s", gr.first.c_str());
    grCmp.SetName(name.c_str());
    // Type of Correlation
    const uint corrType = cmpWS.at(gr.first).second;
    //
    // Loop over each bin
    for (int iBin=0; iBin<grMap.at("Ref").GetN(); iBin++) {
      //
      // Extract the info of the reference graph
      double X_Val_Ref=0.0 , Y_Val_Ref=0.0;
      grMap.at("Ref").GetPoint(iBin, X_Val_Ref, Y_Val_Ref);
      const double X_Err_Ref      = grMap.at("Ref").GetErrorX(iBin);
      const double Y_Err_High_Ref = grMap.at("Ref").GetErrorYhigh(iBin);
      const double Y_Err_Low_Ref  = grMap.at("Ref").GetErrorYlow(iBin);
      //
      // Extract the info of the reference RAW graph
      double X_Val_RefRaw=0.0 , Y_Val_RefRaw=0.0;
      grMapRaw.at("Ref").GetPoint(iBin, X_Val_RefRaw, Y_Val_RefRaw);
      const double X_Err_RefRaw      = grMapRaw.at("Ref").GetErrorX(iBin);
      const double Y_Err_High_RefRaw = grMapRaw.at("Ref").GetErrorYhigh(iBin);
      const double Y_Err_Low_RefRaw  = grMapRaw.at("Ref").GetErrorYlow(iBin);
      if (X_Val_Ref!=X_Val_RefRaw) { std::cout << "[ERROR] X value of ref (" << X_Val_Ref << ") is different from raw (" << X_Val_RefRaw << ") !" << std::endl; return false; }
      if (X_Err_Ref!=X_Err_RefRaw) { std::cout << "[ERROR] X error of ref (" << X_Err_Ref << ") is different from raw (" << X_Err_RefRaw << ") !" << std::endl; return false; }
      //
      // Extract the info of the comparison graph
      double X_Val_Cmp=0.0 , Y_Val_Cmp=0.0;
      grMap.at(gr.first).GetPoint(iBin, X_Val_Cmp, Y_Val_Cmp);
      const double X_Err_Cmp      = gr.second.GetErrorX(iBin);
      const double Y_Err_High_Cmp = gr.second.GetErrorYhigh(iBin);
      const double Y_Err_Low_Cmp  = gr.second.GetErrorYlow(iBin);
      if (X_Val_Ref!=X_Val_Cmp) { std::cout << "[ERROR] X value of ref (" << X_Val_Ref << ") is different from cmp (" << X_Val_Cmp << ") !" << std::endl; return false; }
      if (X_Err_Ref!=X_Err_Cmp) { std::cout << "[ERROR] X error of ref (" << X_Err_Ref << ") is different from cmp (" << X_Err_Cmp << ") !" << std::endl; return false; }
      //
      // Extract the info of the comparison graph
      double X_Val_CmpRaw=0.0 , Y_Val_CmpRaw=0.0;
      grMapRaw.at(gr.first).GetPoint(iBin, X_Val_CmpRaw, Y_Val_CmpRaw);
      const double X_Err_CmpRaw      = grMapRaw.at(gr.first).GetErrorX(iBin);
      const double Y_Err_High_CmpRaw = grMapRaw.at(gr.first).GetErrorYhigh(iBin);
      const double Y_Err_Low_CmpRaw  = grMapRaw.at(gr.first).GetErrorYlow(iBin);
      if (X_Val_Ref!=X_Val_CmpRaw) { std::cout << "[ERROR] X value of ref (" << X_Val_Ref << ") is different from cmp raw (" << X_Val_CmpRaw << ") !" << std::endl; return false; }
      if (X_Err_Ref!=X_Err_CmpRaw) { std::cout << "[ERROR] X error of ref (" << X_Err_Ref << ") is different from cmp raw (" << X_Err_CmpRaw << ") !" << std::endl; return false; }
      //
      // Fill the information into the Pull Graph
      // X axis
      const double X_Val_Pull = X_Val_Ref;
      const double X_Err_Pull = X_Err_Ref;
      // Y axis
      const double  refF = (Y_Val_Ref/Y_Val_RefRaw);
      const double  cmpF = (Y_Val_Cmp/Y_Val_CmpRaw);
      //
      if (std::abs((refF*Y_Err_High_RefRaw)-Y_Err_High_Ref)>0.0001) {
        std::cout << "[ERROR] High Y error of ref (" << Y_Err_High_Ref << ") is different from scaled ref raw (" << refF*Y_Err_High_RefRaw << ") !" << std::endl; return false;
      }
      if (std::abs((refF*Y_Err_Low_RefRaw)-Y_Err_Low_Ref)>0.0001) {
        std::cout << "[ERROR] Low Y error of ref (" << Y_Err_Low_Ref << ") is different from scaled ref raw (" << refF*Y_Err_Low_RefRaw << ") !" << std::endl; return false;
      }
      if (std::abs((cmpF*Y_Err_High_CmpRaw)-Y_Err_High_Cmp)>0.0001) {
        std::cout << "[ERROR] High Y error of cmp (" << Y_Err_High_Cmp << ") is different from scaled cmp raw (" << cmpF*Y_Err_High_CmpRaw << ") !" << std::endl; return false;
      }
      if (std::abs((cmpF*Y_Err_Low_CmpRaw)-Y_Err_Low_Cmp)>0.0001) {
        std::cout << "[ERROR] Low Y error of cmp (" << Y_Err_Low_Cmp << ") is different from scaled cmp raw (" << cmpF*Y_Err_Low_CmpRaw << ") !" << std::endl; return false;
      }
      //
      const double diff = ((cmpF*Y_Val_CmpRaw) - (refF*Y_Val_RefRaw));
      double err_High = 0.0 , err_Low = 0.0;
      if (whatDo=="pull") {
        // Fully Uncorrelated
        if (corrType == 0) {
          err_High = std::sqrt(std::abs((Y_Err_High_Cmp*Y_Err_High_Cmp) + (Y_Err_High_Ref*Y_Err_High_Ref)));
          err_Low  = std::sqrt(std::abs((Y_Err_Low_Cmp *Y_Err_Low_Cmp ) + (Y_Err_Low_Ref *Y_Err_Low_Ref )));
        }
        // Fully Correlated
        else if (corrType == 1) {
          err_High = std::sqrt(std::abs((Y_Err_High_Cmp*Y_Err_High_Cmp) + (Y_Err_High_Ref*Y_Err_High_Ref) - 2.0*refF*cmpF*(Y_Err_High_RefRaw*Y_Err_High_RefRaw)));
          err_Low  = std::sqrt(std::abs((Y_Err_Low_Cmp *Y_Err_Low_Cmp ) + (Y_Err_Low_Ref *Y_Err_Low_Ref ) - 2.0*refF*cmpF*(Y_Err_Low_RefRaw *Y_Err_Low_RefRaw )));
        }
        // Partially Correlated
        else if (corrType == 2) {
          err_High = Y_Err_High_Cmp;
          err_Low  = Y_Err_Low_Cmp;
        }
        else { std::cout << "[ERROR] Correction type is invalid!" << std::endl; return false; }
      }
      if (whatDo=="ratio") {
        // Fully Uncorrelated
        if (corrType == 0) {
          double Y_Val_Pull = (Y_Val_Cmp/Y_Val_Ref);
          err_High = std::abs(Y_Val_Pull)*std::sqrt(std::pow(Y_Err_High_Cmp/Y_Val_Cmp, 2.0) + std::pow(Y_Err_High_Ref/Y_Val_Ref, 2.0));
          err_Low  = std::abs(Y_Val_Pull)*std::sqrt(std::pow(Y_Err_Low_Cmp/Y_Val_Cmp, 2.0) + std::pow(Y_Err_Low_Ref/Y_Val_Ref, 2.0));
        }
        // Fully Correlated
        else if (corrType == 1) {
          double Y_Val_Pull = (Y_Val_Cmp/Y_Val_Ref);
          err_High = std::abs(Y_Val_Pull)*std::sqrt(std::abs(std::pow(Y_Err_High_Cmp/Y_Val_Cmp, 2.0) + std::pow(Y_Err_High_Ref/Y_Val_Ref, 2.0) - (2.0/(Y_Val_Ref*Y_Val_Cmp))*std::pow(Y_Err_High_Cmp, 2.0)));
          err_Low  = std::abs(Y_Val_Pull)*std::sqrt(std::abs(std::pow(Y_Err_Low_Cmp/Y_Val_Cmp, 2.0) + std::pow(Y_Err_Low_Ref/Y_Val_Ref, 2.0) - (2.0/(Y_Val_Ref*Y_Val_Cmp))*std::pow(Y_Err_Low_Cmp, 2.0)));
        }
        // Partially Correlated
        else if (corrType == 2) {
          err_High = Y_Err_High_Cmp;
          err_Low  = Y_Err_Low_Cmp;
        }
        else { std::cout << "[ERROR] Correction type is invalid!" << std::endl; return false; }
      }
      double Y_Val_Pull = 0.0 , Y_Err_High_Pull = 0.0 , Y_Err_Low_Pull = 0.0;
      const double err = ((diff>0.0) ? err_High : err_Low);
      if (whatDo=="pull") {
        Y_Val_Pull      = ( (err>0.0) ? (diff/err)     : 0.0 );
        Y_Err_High_Pull = ( (err>0.0) ? (err_High/err) : 0.0 );
        Y_Err_Low_Pull  = ( (err>0.0) ? (err_Low /err) : 0.0 );
      }
      if (whatDo=="ratio") {
        Y_Val_Pull      = (Y_Val_Cmp/Y_Val_Ref) - 1.0;
        Y_Err_High_Pull = err_High;
        Y_Err_Low_Pull  = err_Low;
      }
      //
      // Fill the pull graph
      grCmp.SetPoint(iBin, X_Val_Pull, Y_Val_Pull);
      grCmp.SetPointError(iBin, X_Err_Pull, X_Err_Pull, Y_Err_Low_Pull, Y_Err_High_Pull);
    }
  }
  //
  return true;
};


void drawCompareGraph(GraphMap& grMap, GraphMap& grCmpMap, const std::string& outDir, const std::string& chg, const StringMap& legLbl,
                      const std::string col, const std::string whatDo, const bool doCorrectYield)
{
  //
  // Create Canvas
  TCanvas c("c", "c", 1000, 1000); c.cd();
  //
  // Create the Text Info
  TLatex tex; tex.SetNDC(); tex.SetTextSize(0.035); float dy = 0;
  std::vector< std::string > textToPrint;
  std::string sampleLabel = "W #rightarrow #mu + #nu_{#mu}";
  if (chg == "Pl") { sampleLabel = "W^{+} #rightarrow #mu^{+} + #nu_{#mu}"; }
  if (chg == "Mi") { sampleLabel = "W^{-} #rightarrow #mu^{-} + #nu_{#mu}"; }
  textToPrint.push_back(sampleLabel);
  textToPrint.push_back("p^{#mu}_{T} > 25 GeV/c");
  //
  // Initialize the Legend
  TLegend leg(0.2, 0.66, 0.4, 0.79);
  // Initialize the graph x range variables
  double xMin=0.0 , xMax=0.0 , yErrMin=9999999. , yErrMax=-1.;
  //
  // Format the Graphs
  //
  //  Main Frame
  for (auto& gr : grMap) {
    gr.second.SetTitle("");
    gr.second.SetMarkerColor(kBlue);
    gr.second.SetMarkerStyle(20);
    gr.second.GetXaxis()->SetTitle("");
    gr.second.GetXaxis()->SetTitleSize(0.05);
    gr.second.GetXaxis()->SetTitleFont(42);
    gr.second.GetXaxis()->SetTitleOffset(3);
    gr.second.GetXaxis()->SetLabelOffset(3);
    gr.second.GetXaxis()->SetLimits(-3.0, 2.0);
    if (doCorrectYield) { gr.second.GetYaxis()->SetTitle("Corrected Signal Yield"); }
    else                { gr.second.GetYaxis()->SetTitle("Raw Signal Yield");       }
    gr.second.GetYaxis()->SetLabelSize(0.044);
    gr.second.GetYaxis()->SetTitleSize(0.044);
    gr.second.GetYaxis()->SetTitleOffset(1.7);
    gr.second.GetYaxis()->SetTitleFont(42);
    gr.second.GetYaxis()->SetRangeUser(0.0, 10000.);
  }
  grMap.at("Ref").SetMarkerColor(kBlack);
  uint iGr = 0; std::vector< uint > COLOR = { kRed , (kGreen+2) , (kBlue+2) };
  for (auto& gr : grMap) { if (gr.first!="Ref") { gr.second.SetMarkerColor(COLOR[iGr]); iGr++; } }
  //
  //  Pull Frame
  for (auto& gr : grCmpMap) {
    gr.second.SetTitle("");
    gr.second.SetMarkerColor(kBlue);
    gr.second.SetMarkerStyle(20);
    gr.second.GetXaxis()->SetTitleOffset(1);
    gr.second.GetXaxis()->SetTitleSize(0.16);
    gr.second.GetXaxis()->SetLabelSize(0.14);
    gr.second.GetXaxis()->SetTitle("Muon #eta_{CM}");
    gr.second.GetXaxis()->SetLimits(-3.0, 2.0);
    if (whatDo=="pull" ) { gr.second.GetYaxis()->SetTitle("Pull"); }
    if (whatDo=="ratio") { gr.second.GetYaxis()->SetTitle("#frac{New-Old}{Old}"); }
    gr.second.GetYaxis()->CenterTitle(kTRUE);
    gr.second.GetYaxis()->SetTitleOffset(0.4);
    gr.second.GetYaxis()->SetTitleSize(0.16);
    gr.second.GetYaxis()->SetLabelSize(0.11);
    gr.second.GetYaxis()->SetNdivisions(404);
    TGaxis::SetMaxDigits(2); // to display powers of 10
    if (whatDo=="pull" ) { gr.second.GetYaxis()->SetRangeUser(-4.0, 4.0); }
    if (whatDo=="ratio") { gr.second.GetYaxis()->SetRangeUser(-0.02, 0.02); }
  }
  iGr = 0;
  for (auto& gr : grCmpMap) { if (gr.first!="Ref") { gr.second.SetMarkerColor(COLOR[iGr]); iGr++; } }
  //
  // Define the plotting pads
  TPad *pad1  = new TPad("pad1", "", 0, 0.23, 1, 1);  // Unique Pointer does produce Segmentation Fault, so don't use it
  TPad *pad2  = new TPad("pad2", "", 0, 0, 1, 0.228); // Unique Pointer does produce Segmentation Fault, so don't use it
  auto  pline = std::unique_ptr<TLine>(new TLine(-3.0, 0.0,  2.0, 0.0));
  // Format the pads
  c.cd();
  pad2->SetTopMargin(0.035);
  pad2->SetBottomMargin(0.4);
  pad2->SetFillStyle(4000); 
  pad2->SetFrameFillStyle(4000);
  pad1->SetBottomMargin(0.025);
  //
  // Create Legend
  for (auto& gr : grMap) { formatLegendEntry(*leg.AddEntry(&gr.second, legLbl.at(gr.first).c_str(), "pe")); }
  //
  // Draw the Graphs
  //
  // Main Frame
  pad1->Draw();
  pad1->cd();
  grMap.at("Ref").Draw("ap");
  for (auto& gr : grMap) { gr.second.Draw("samep"); }
  grMap.at("Ref").Draw("samep");
  //
  // Draw the Legend
  leg.Draw("same");
  //
  // Draw the text
  for (const auto& s: textToPrint) { tex.DrawLatex(0.22, 0.86-dy, s.c_str()); dy+=0.045; }
  //
  // Apply CMS style to pad
  int lumiId = 0;
  if (col=="pPb") { lumiId = 109; } else if (col=="Pbp") { lumiId = 110; } else if (col=="PA") { lumiId = 111; }
  CMS_lumi(pad1, lumiId, 33, "");
  gStyle->SetTitleFontSize(0.05);
  pad1->Update();
  //
  // Pull Frame
  c.cd();
  pad2->Draw();
  pad2->cd();
  grCmpMap.begin()->second.Draw("ap");
  for (auto& gr : grCmpMap) { gr.second.Draw("samep"); }
  //
  pline->Draw("same");
  pad2->Update();
  //
  // Save the pads
  //
  // Create Output Directory
  const std::string plotDir = outDir+"/Plots/" + col;
  makeDir(plotDir + "/png/");
  makeDir(plotDir + "/pdf/");
  makeDir(plotDir + "/root/");
  //
  // Save Canvas
  const std::string name = std::string(grMap.at("Ref").GetName()) + (doCorrectYield ? "_Corr" : "_Raw");
  c.SaveAs(( plotDir + "/png/"  + name + ".png"  ).c_str());
  c.SaveAs(( plotDir + "/pdf/"  + name + ".pdf"  ).c_str());
  c.SaveAs(( plotDir + "/root/" + name + ".root" ).c_str());
  //
  // Clean up memory
  c.Clear(); c.Close();
  //
};


void fillCompareGraph(TGraphAsymmErrors& graph, const BinMapMap& var, const std::string graphLbl)
{
  // Set Bin Size
  const uint nBins = var.at("Val").size();
  graph.Set(nBins);
  // Set Graph Name
  const std::string name = Form("gr_WToMu_%s", graphLbl.c_str());
  graph.SetName(name.c_str());
  // Determine index of bins
  const auto& binMap = var.at("Val");
  std::map< anabin<0> , uint > binIdx;
  uint iBin = 0; for (const auto& b : binMap) { binIdx[b.first] = iBin; iBin++; }
  // Fill the graph
  for (const auto& b : binMap) {
    const unsigned int iBin = binIdx.at(b.first);
    //
    // Extract the parameters needed for each axis
    //
    // X Value
    const double X = ( (b.first.etabin().high() + b.first.etabin().low()) / 2.0 ); // Mean value of eta bin
    // X Error
    const double Err_X = ( (b.first.etabin().high() - b.first.etabin().low()) / 2.0 ); // Width of eta bin
    double Err_X_High = Err_X;
    double Err_X_Low  = Err_X;
    // Y Value
    const double Y = var.at("Val").at(b.first);
    //
    // Compute total statistic error
    const double Err_Y_High = var.at("Err_High").at(b.first);
    const double Err_Y_Low  = var.at("Err_Low" ).at(b.first);
    //
    // Fill the nominal graph
    //
    graph.SetPoint(iBin, X, Y);
    graph.SetPointError(iBin, Err_X_Low, Err_X_High, Err_Y_Low, Err_Y_High);
  }
};


bool extractInfo(BinQuadMap& inputVar, const std::string& workDirPath)
{
  //
  // --------------------------------------------------------------------------------- //
  //
  // Get the list of input files
  //
  std::vector< std::string > inputFileNames;
  const std::string inputDirPath = workDirPath;
  if (!fileList(inputFileNames, inputDirPath)) { return false; };
  //
  // --------------------------------------------------------------------------------- //
  //
  // Loop over the input files
  //
  for (const auto& inputFileName : inputFileNames) {
    //
    std::cout << "Processing file: " << inputFileName << std::endl;
    //
    // Open input file
    const std::string inputFilePath = Form("%s/%s", inputDirPath.c_str(), inputFileName.c_str());
    TFile inputFile(inputFilePath.c_str(), "READ");
    //
    if (inputFile.IsOpen()==false || inputFile.IsZombie()==true) {
      std::cout << "[ERROR] The input file " << inputFilePath << " could not be open!" << std::endl; return false;
    }
    //
    // Extract the Workspace
    RooWorkspace* ws = (RooWorkspace*) inputFile.Get("workspace");
    if (ws == NULL) { std::cout << "[ERROR] File: " << inputFilePath << " does not have the workspace!" << std::endl; inputFile.Close(); return false; }
    //
    // Extract the information from the workspace
    const std::string COL   = (ws->obj("fitSystem")) ? ((RooStringVar*)ws->obj("fitSystem"))->getVal() : "";
    const std::string CHG   = (ws->obj("fitCharge")) ? ((RooStringVar*)ws->obj("fitCharge"))->getVal() : "";
    const std::string token = ( CHG + "_" + COL );
    //
    // Extract the Flag Information
    bool useEtaCM = false;
    if ( (ws->var("useEtaCM") != NULL) && (ws->var("useEtaCM")->getVal() == 1.0) ) { useEtaCM = true; }
    if (useEtaCM==false) { std::cout << "[ERROR] The eta variable is not CM!" << std::endl; inputFile.Close(); return false; }
    //
    // Extract the Dataset Variable Information
    RooRealVar*  etaVar = (RooRealVar*) ws->var("Muon_Eta");
    const double etaMin = ( etaVar  ? etaVar->getMin() : -1.0 );
    const double etaMax = ( etaVar  ? etaVar->getMax() : -1.0 );
    const double etaMin_LAB = (COL=="Pbp" ? -etaMax : etaMin);
    const double etaMax_LAB = (COL=="Pbp" ? -etaMin : etaMax);
    const double etaMin_CM  = PA::EtaLABtoCM(etaMin_LAB, true);
    const double etaMax_CM  = PA::EtaLABtoCM(etaMax_LAB, true);
    const auto bin = anabin<0>(etaMin_CM, etaMax_CM);
    if (std::abs(etaMax - etaMin) > 2.0) continue;
    //
    // Extract the Fitted Variable Information
    RooRealVar*  yieldVar   = (RooRealVar*) ws->var(Form("N_WToMu%s", token.c_str()));
    const double yieldVal   = ( yieldVar ? yieldVar->getVal()   : -1.0 );
    const double yieldErrLo = ( yieldVar ? yieldVar->getError() : -1.0 );
    const double yieldErrHi = ( yieldVar ? yieldVar->getError() : -1.0 );
    //
    // Fill the information
    inputVar[CHG]["N_WToMu"]["Val"][bin] = yieldVal;
    inputVar[CHG]["N_WToMu"]["Err_Low" ][bin] = yieldErrLo;
    inputVar[CHG]["N_WToMu"]["Err_High"][bin] = yieldErrHi;
    //
    // Clean up the memory
    delete ws;
    inputFile.Close();
  } // loop on the files
  //
  return true;
};


void createCompareTable(std::vector< std::string >& texTable, const std::vector< std::string >& colLbl, const std::vector< std::string >& colVar,
			   const std::vector< std::string >& colTitle, const GraphMap& grMap, const std::string& col)
{
  //
  const uint nCol = colVar.size();
  //
  texTable.push_back("  \\renewcommand{\\arraystretch}{1.5}");
  texTable.push_back(Form("  \\begin{tabular}{|c|*%dc|}", (nCol-1)));
  texTable.push_back("    \\hline");
  std::string tmp;
  tmp = ("    ");
  for (uint i = 0; i < nCol; i++) {
    tmp += colTitle[i];
    if (i<(nCol-1)) { tmp += " & "; }
    else { tmp += "\\\\"; }
  }
  texTable.push_back(tmp);
  texTable.push_back("    \\hline\\hline");
  //
  for (int iBin = 0; iBin < grMap.at("Ref").GetN(); iBin++) {
    tmp = ("    ");
    for (uint i = 0; i < nCol; i++) {
      //
      const auto&   v = colVar[i];
      const auto& lbl = colLbl[i];
      std::string val;
      //
      if (v=="Muon_Eta") {
	double x, y; grMap.at(lbl).GetPoint(iBin, x, y);
	const double errX = grMap.at(lbl).GetErrorX(iBin);
	const double min = x-errX;
	const double max = x+errX;
	val = Form("%s%.2f , %s%.2f", sgn(min), min, sgn(max) , max);
      }
      else if (v=="N_WToMu") {
	double x, y; grMap.at(lbl).GetPoint(iBin, x, y);
	const double errY = grMap.at(lbl).GetErrorY(iBin);
	val = Form("$%.2f \\pm %.2f$", y, errY);
      }
      else if (v=="Diff") {
	double x_Cmp, y_Cmp; grMap.at(lbl).GetPoint(iBin, x_Cmp, y_Cmp);
	double x_Ref, y_Ref; grMap.at("Ref").GetPoint(iBin, x_Ref, y_Ref);
	const double errY_Cmp = grMap.at(lbl).GetErrorY(iBin);
	const double errY_Ref = grMap.at("Ref").GetErrorY(iBin);
	val = Form("$%.2f \\pm %.2f$", (y_Ref - y_Cmp), errY_Cmp);
      }
      tmp += val;
      if (i<(nCol-1)) { tmp += " & "; } else { tmp += "\\\\"; }
    }
    texTable.push_back(tmp);
    texTable.push_back("    \\hline");
  }
  //
  texTable.push_back("  \\end{tabular}");
  //
};


void makeCompareTable(std::ofstream& file, GraphMap& grMap, const std::string& chg, const std::string& refLbl, const std::string col, const bool doCorrectYield)
{
  //
  // Initialize the latex table
  std::vector< std::string > texTable;
  //
  const std::string pTCUT = Form("$p_{T} > %.0f$~GeV/c", 25.0);
  //
  // Determine number of columns
  std::vector< std::string > colLbl , colVar , colTitle;
  //
  colLbl.push_back("Ref");
  colVar.push_back("Muon_Eta");
  colTitle.push_back("$\\eta_{LAB}$ Range");
  //
  colLbl.push_back("Ref");
  colVar.push_back("N_WToMu");
  colTitle.push_back(refLbl);
  //
  for (const auto& gr : grMap) {
    if (gr.first=="Ref") continue;
    colLbl.push_back(gr.first);
    colLbl.push_back(gr.first);
    colVar.push_back("N_WToMu");
    colVar.push_back("Diff");
    colTitle.push_back(gr.first);
    colTitle.push_back(Form("%s - Ref", gr.first.c_str()));
  }
  //
  texTable.push_back("\\begin{table}[h!]");
  texTable.push_back("  \\centering");
  if (colTitle.size()>5) { texTable.push_back("  \\resizebox{\\textwidth}{!}{"); }
  createCompareTable(texTable, colLbl, colVar, colTitle, grMap, col);
  texTable.push_back("  }");
  texTable.push_back(Form("  \\caption{%s}",
                          Form("Maximum systematic error of the measured observables determined for each systematic category in the %s collision system. All analysis cuts are applied%s.",
                               (col=="PA" ? "combined \\pPb and \\Pbp" : (col=="pPb" ? "\\pPb" : (col=="Pbp" ? "\\Pbp" : "????"))),
                               (pTCUT!="" ? Form(" including the muon %s cut", pTCUT.c_str()) : "")
                               )
                          )
                     );
  texTable.push_back("  \\label{tab:Systematics}");
  texTable.push_back("\\end{table}");
  //
  for (const auto& row : texTable) { file << row << std::endl; }
  file << std::endl; file << std::endl;
  //
};


bool printCompareTable(GraphMap& grMap, const std::string& outDir, const std::string& chg, const std::string& refLbl, const std::string col, const bool doCorrectYield)
{
  //
  std::cout << "[INFO] Filling the comparison table" << std::endl;
  //
  // Fill the tables
  //
  // Create Output Directory
  const std::string tableDir = outDir + "/Tables/Comparison/" + col;
  makeDir(tableDir);
  // Create Output Files for Full Systematics
  const std::string fileName_SYST = "systematic_" + col + "_" + chg + "_" + refLbl;
  std::ofstream file_SYST((tableDir + "/" + fileName_SYST + ".tex").c_str());
  if (file_SYST.is_open()==false) { std::cout << "[ERROR] File " << fileName_SYST << " was not created!" << std::endl; return false; }
  //
  makeCompareTable(file_SYST, grMap, chg, refLbl, col, doCorrectYield);
  //
  return true;
};


bool getObjectsFromFile(EffMap_t& eff, Unc1DMap_t& unc, const std::string& filePath)
{
  // Open the input file
  TFile file(filePath.c_str(), "READ");
  if (file.IsZombie()) { std::cout << "[ERROR] The input file " << filePath << " was not found!" << std::endl; return false; }
  file.cd();
  //
  //gain time, do not add the objects in the list in memory
  Bool_t status = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  //
  const std::string mainDirName = gDirectory->GetListOfKeys()->First()->GetName();
  TDirectory* mainDir = file.GetDirectory(mainDirName.c_str());
  if (mainDir==NULL) { std::cout << "[ERROR] Dir: " << mainDirName << " was not found!" << std::endl; return false; }
  mainDir->cd();
  // loop over all keys in this directory
  for (const auto& v : getK(mainDir->GetListOfKeys())) {
    if (v->ReadObj()->IsA()->InheritsFrom(TDirectory::Class()) == false) continue;
    TDirectory* varDir = mainDir->GetDirectory(v->GetName());
    varDir->cd();
    for (const auto& s : getK(varDir->GetListOfKeys())) {
      if (s->ReadObj()->IsA()->InheritsFrom(TDirectoryFile::Class()) == false) continue;
      TDirectory* sampleDir = varDir->GetDirectory(s->GetName());
      sampleDir->cd();
      for (const auto& c : getK(sampleDir->GetListOfKeys())) {
        if (c->ReadObj()->IsA()->InheritsFrom(TDirectoryFile::Class()) == false) continue;
        TDirectory* colDir = sampleDir->GetDirectory(c->GetName());
        colDir->cd();
        for (const auto& t : getK(colDir->GetListOfKeys())) {
          if (t->ReadObj()->IsA()->InheritsFrom(TDirectoryFile::Class()) == false) continue;
          TDirectory* typeDir = colDir->GetDirectory(t->GetName());
          typeDir->cd();
          for (const auto& co : getK(typeDir->GetListOfKeys())) {
            if (co->ReadObj()->IsA()->InheritsFrom(TDirectoryFile::Class())) {
              TDirectory* corDir = typeDir->GetDirectory(co->GetName());
              corDir->cd();
              for (auto& ch : CHG_) {
                for (const auto& key : getK(corDir->GetListOfKeys())) {
                  if (std::string(key->GetName()).find(ch)!=std::string::npos) {
                    if (key->ReadObj()->IsA()->InheritsFrom(TEfficiency::Class())) {
                      eff[v->GetName()][s->GetName()][c->GetName()][ch][t->GetName()][co->GetName()].push_back(*((TEfficiency*)key->ReadObj()));
                      const uint i = ( eff.at(v->GetName()).at(s->GetName()).at(c->GetName()).at(ch).at(t->GetName()).at(co->GetName()).size() - 1);
                      eff.at(v->GetName()).at(s->GetName()).at(c->GetName()).at(ch).at(t->GetName()).at(co->GetName())[i].SetName(key->GetName());
                    }
                    if (key->ReadObj()->IsA()->InheritsFrom(TVectorD::Class())) {
                      unc[v->GetName()][s->GetName()][c->GetName()][ch][t->GetName()][co->GetName()].ResizeTo( ((TVectorD*)key->ReadObj())->GetNrows() );
                      unc.at(v->GetName()).at(s->GetName()).at(c->GetName()).at(ch).at(t->GetName()).at(co->GetName()) = *((TVectorD*)key->ReadObj());
                    }
                  }
                }
                if (eff.at(v->GetName()).at(s->GetName()).at(c->GetName()).at(ch).at(t->GetName()).count("NoCorr")>0) {
                  eff[v->GetName()][s->GetName()][c->GetName()][ch][t->GetName()]["NoCorrOnly"] = eff.at(v->GetName()).at(s->GetName()).at(c->GetName()).at(ch).at(t->GetName()).at("NoCorr");
                  std::string name = eff.at(v->GetName()).at(s->GetName()).at(c->GetName()).at(ch).at(t->GetName()).at("NoCorr")[0].GetName();
                  name.replace(name.find("NoCorr"), name.find("NoCorr")+6, "NoCorrOnly");
                  eff.at(v->GetName()).at(s->GetName()).at(c->GetName()).at(ch).at(t->GetName()).at("NoCorrOnly")[0].SetName(name.c_str());
                }
                if (eff.at(v->GetName()).at(s->GetName()).at(c->GetName()).at(ch).at(t->GetName()).count("HFCorr")>0) {
                  eff[v->GetName()][s->GetName()][c->GetName()][ch][t->GetName()]["HFCorrOnly"] = eff.at(v->GetName()).at(s->GetName()).at(c->GetName()).at(ch).at(t->GetName()).at("HFCorr");
                  std::string name = eff.at(v->GetName()).at(s->GetName()).at(c->GetName()).at(ch).at(t->GetName()).at("HFCorr")[0].GetName();
                  name.replace(name.find("HFCorr"), name.find("HFCorr")+6, "HFCorrOnly");
                  eff.at(v->GetName()).at(s->GetName()).at(c->GetName()).at(ch).at(t->GetName()).at("HFCorrOnly")[0].SetName(name.c_str());
                }
              }
              typeDir->cd();
            }
            if (co->ReadObj()->IsA()->InheritsFrom(TVectorD::Class())) {
              for (auto& ch : CHG_) {
                if (std::string(co->GetName()).find(ch)!=std::string::npos) {
                  std::string corrName = "";
                  if (std::string(co->GetName()).find("TnP_Stat")!=std::string::npos) { corrName = "TnP_Stat"; }
                  if (std::string(co->GetName()).find("TnP_Syst")!=std::string::npos) { corrName = "TnP_Syst"; }
                  if (std::string(co->GetName()).find("TnP_Tot" )!=std::string::npos) { corrName = "TnP_Tot";  }
                  if (std::string(co->GetName()).find("MC_Syst" )!=std::string::npos) { corrName = "MC_Syst";  }
                  if (corrName!="") {
                    unc[v->GetName()][s->GetName()][c->GetName()][ch][t->GetName()][corrName].ResizeTo( ((TVectorD*)co->ReadObj())->GetNrows() );
                    unc.at(v->GetName()).at(s->GetName()).at(c->GetName()).at(ch).at(t->GetName()).at(corrName) = *((TVectorD*)co->ReadObj());
                  }
                }
              }
            }
          }
          colDir->cd();
        }
        sampleDir->cd();
      }
      varDir->cd();
    }
    mainDir->cd();
  }
  // Close File
  file.Close("R");
  // return
  return true;
};


KeyPVec_t getK(TList* list)
{
  TIter iter(list);
  KeyPVec_t out;
  while (TKey* key = (TKey*)iter()) { out.push_back(key); }
  return out;
};


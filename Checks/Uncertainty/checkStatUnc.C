#ifndef checkStatUnc_C
#define checkStatUnc_C

// Auxiliary Headers
#include "../../Fitter/Macros/Utilities/initClasses.h"
// ROOT headers
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TLatex.h"
// RooFit headers
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooAbsPdf.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooMCStudy.h"
// c++ headers
#include <dirent.h>
#include <iostream>
#include <string>


bool        fileList(std::vector< std::string >& fileNames, const std::string& dirPath);
void        printElectroWeakBinning(TCanvas& pad, const RooWorkspace& ws, const std::string& dsName);
void        setRange(RooPlot& frame);
void        formatFrame(RooPlot& frame);
RooRealVar  extractPar(RooPlot& frame);
RooRealVar  extractPar(RooPlot& frame, const RooRealVar& parIni, const RooDataSet& dataIni);

bool checkStatUnc(
                  const std::string workDirName = "NominalCM",
                  const std::string metTag      = "METPF_RAW",
                  const std::string dsTag       = "DATA",
                  const std::string colTag      = "PA"
                  )
{
  //
  // Set the CMS style
  setTDRStyle();
  //
  // --------------------------------------------------------------------------------- //
  //
  // Define the output file info
  const std::string CWD = getcwd(NULL, 0);
  const std::string outputDirPath  = Form("%s/Tree/%s/%s/%s/%s", CWD.c_str(), workDirName.c_str(), metTag.c_str(), dsTag.c_str(), colTag.c_str());
  const std::string outputFileName = "tree_allvars.root";
  const std::string outputFilePath = Form("%s/%s", outputDirPath.c_str(), outputFileName.c_str());
  //
  // --------------------------------------------------------------------------------- //
  //
  // Get the list of input files
  //
  std::vector< std::string > inputFileNames;
  std::string preCWD  = CWD; preCWD.erase(preCWD.find_last_of("/"), 100);
  std::string pre2CWD = preCWD; pre2CWD.erase(pre2CWD.find_last_of("/"), 100);
  const std::string inputDirPath = Form("%s/Fitter/Output/%s/%s/%s/%s/result", pre2CWD.c_str(), workDirName.c_str(), metTag.c_str(), dsTag.c_str(), colTag.c_str());
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
      std::cout << "[ERROR] The input file " << inputFilePath << " could not be created!" << std::endl; return false;
    }
    //
    // Extract the Workspace
    auto ws = (RooWorkspace*) inputFile.Get("workspace");
    if (ws == NULL) { std::cout << "[ERROR] File: " << inputFilePath << " does not have the workspace!" << std::endl; inputFile.Close(); return false; }
    //
    // Extract the information from the workspace
    const std::string DSTAG = (ws->obj("DSTAG"))     ? ((TObjString*)ws->obj("DSTAG"))->GetString().Data()     : "";
    const std::string CHA   = (ws->obj("channel"))   ? ((TObjString*)ws->obj("channel"))->GetString().Data()   : "";
    const std::string COL   = (ws->obj("fitSystem")) ? ((TObjString*)ws->obj("fitSystem"))->GetString().Data() : "";
    const std::string CHG   = (ws->obj("fitCharge")) ? ((TObjString*)ws->obj("fitCharge"))->GetString().Data() : "";
    const std::string OBJ   = (ws->obj("fitObject")) ? ((TObjString*)ws->obj("fitObject"))->GetString().Data() : "";
    // Check the information
    if (DSTAG.find(dsTag)==std::string::npos) { std::cout << "[ERROR] Workspace DSTAG " << DSTAG << " is not consistent with input dsTag " << dsTag << std::endl; inputFile.Close(); return false; }
    if (COL != colTag ) { std::cout << "[ERROR] Workspace COL " << COL << " is not consistent with input colTag " << colTag << std::endl; inputFile.Close(); return false; }
    if (OBJ != "W"    ) { std::cout << "[ERROR] Only W fits are currently supported in result macros!"      << std::endl; return false; }
    if (CHA != "ToMu" ) { std::cout << "[ERROR] Only muon channel is currently supported in result macros!" << std::endl; return false; }
    //
    const std::string token  = ( CHG + "_" + COL );
    const std::string tag    = ( OBJ + CHA + token );
    //
    bool useEtaCM = false;
    if ( (ws->var("useEtaCM") != NULL) && (ws->var("useEtaCM")->getVal() == 1.0) ) { useEtaCM = true; }
    //
    const std::string dsName = ( "d" + CHG + "_" + DSTAG );
    auto ds = ws->data(Form("%s", dsName.c_str()));
    if (ds == NULL) { std::cout << "[ERROR] The Dataset " << dsName << " was not found in the workspace" << std::endl; inputFile.Close(); return false; }
    //
    const std::string pdfName = ( "pdfMET_Tot" + tag );
    auto pdf = (RooAddPdf*) ws->pdf(Form("%s", pdfName.c_str()));
    if (pdf == NULL) { std::cout << "[ERROR] The PDF " << pdfName << " was not found in the workspace" << std::endl; inputFile.Close(); return false; }
    //
    auto var = (RooRealVar*) ws->var("MET");
    if (var == NULL) { std::cout << "[ERROR] The variable MET was not found in the workspace" << std::endl; inputFile.Close(); return false; }
    var->setBins(100);
    //
    auto par = (RooRealVar*) ws->var(Form("N_%s", tag.c_str()));
    if (par == NULL) { std::cout << "[ERROR] The parameter " << Form("N_%s", tag.c_str()) << " was not found in the workspace" << std::endl; inputFile.Close(); return false; }
    RooRealVar parOrig = *par;
    //
    RooMCStudy mcstudy(*pdf, *var, RooFit::Binned(kFALSE), RooFit::Silence(), RooFit::Extended(), RooFit::FitOptions(RooFit::Save(), RooFit::PrintLevel(-1)));
    //
    // ---------------------------------------------
    // Generate and fit 1000 samples of Poisson(nExpected) events
    mcstudy.generateAndFit(1000, ds->sumEntries());
    //
    // ------------------------------------------------
    // Make plots of the distributions of Signal Yield, the error on Signal Yield and the pull of Signal Yield
    auto framePull = mcstudy.plotPull(parOrig, RooFit::Bins(50), RooFit::Range(-5.0, 5.0), RooFit::FitGauss(kTRUE));
    //
    auto textBox = (TPaveText*)framePull->findObject("pullGauss_paramBox");
    if (textBox) {
      textBox->SetX1(0.18); textBox->SetX2(0.5); textBox->SetY1(0.63); textBox->SetY2(0.73); textBox->SetTextSize(0.025);
      textBox->SetBorderSize(0);
    }
    //
    // ------------------------------------------------
    // Draw all plots on a canvas
    auto cFig = std::unique_ptr<TCanvas>(new TCanvas("", "", 1000, 1000)); cFig->cd();
    framePull->Draw();
    //
    // Format the frame
    formatFrame(*framePull);
    //
    float xPos = 0.7, yPos = 0.89, dYPos = 0.045, dy = 0.025*8.;
    TLatex t = TLatex(); t.SetNDC(); t.SetTextSize(0.025);
    t.DrawLatex(xPos, yPos-dy, Form("N_{FIT} = %.0f #pm %.0f", parOrig.getVal(), parOrig.getError())); dy+=dYPos;
    auto framePar  = mcstudy.plotParam(parOrig);
    auto toyVal = extractPar(*framePar);
    t.DrawLatex(xPos, yPos-dy, Form("#mu_{TOY} = %.0f #pm %.0f", toyVal.getVal(), toyVal.getError())); dy+=dYPos;
    auto frameErr  = mcstudy.plotError(parOrig);
    auto toyErr = extractPar(*frameErr);
    t.DrawLatex(xPos, yPos-dy, Form("#sigma_{TOY} = %.0f #pm %.0f", toyErr.getVal(), toyErr.getError())); dy+=dYPos;
    //
    printElectroWeakBinning(*cFig, *ws, dsName);
    //
    int lumiId = 0;
    if (COL=="pPb") { lumiId = 115; } else if (COL=="Pbp") { lumiId = 116; } else if (COL=="PA") { lumiId = 117; }
    CMS_lumi(cFig.get(), lumiId, 33, "");
    //
    // Save the plot in different formats
    const std::string outputDir = CWD+"/MCStatStudy/" + workDirName+"/" + metTag+"/" + dsTag+"/" + colTag+"/";
    const std::string fileName = Form("PULL_%s_%s_%s_%s_%.0f_%.0f", "MET",
                                      tag.c_str(),
                                      dsTag.c_str(),
                                      ( ws->var("useEtaCM") ? "MuEtaCM" : "MuEta" ),
                                      (ws->var("Muon_Eta")->getMin()*100.0), (ws->var("Muon_Eta")->getMax()*100.0)
                                      );
    //
    gSystem->mkdir(Form("%splot/root/", outputDir.c_str()), kTRUE);
    cFig->SaveAs(Form("%splot/root/%s.root", outputDir.c_str(), fileName.c_str()));
    gSystem->mkdir(Form("%splot/png/", outputDir.c_str()), kTRUE);
    cFig->SaveAs(Form("%splot/png/%s.png", outputDir.c_str(), fileName.c_str()));
    gSystem->mkdir(Form("%splot/pdf/", outputDir.c_str()), kTRUE);
    cFig->SaveAs(Form("%splot/pdf/%s.pdf", outputDir.c_str(), fileName.c_str()));
    //
    cFig->Clear();
    cFig->Close();
    //
    // Clean up the memory
    delete ws;
    inputFile.Close();
    //
  }
  //
  return true;
};


bool fileList(std::vector< std::string >& fileNames, const std::string& dirPath)
{
  // Open the directory
  DIR * dpdf = opendir(dirPath.c_str());
  // Search for all the files inside the directory
  if (dpdf != NULL){
    struct dirent *epdf;
    while ((epdf = readdir(dpdf))){
      if (strcmp(epdf->d_name,".")!=0 && strcmp(epdf->d_name,"..")!=0 ) {
        //std::cout << "[INFO] Adding file: " << epdf->d_name << std::endl;
        fileNames.push_back(epdf->d_name);
      }
    }
  } else {
    std::cout << "[ERROR] Working directory ( " << dirPath << " ) was not found!" << endl; return false;
  }
  return true;
};


void formatFrame(RooPlot& frame)
{
  // Set the format of the frame
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  // General
  frame.SetTitle("");
  frame.SetMarkerColor(kBlue);
  frame.SetMarkerStyle(20);
  // X-axis
  frame.GetXaxis()->CenterTitle(kFALSE);
  frame.GetXaxis()->SetTitleOffset(1.1);
  frame.GetXaxis()->SetTitleSize(0.04);
  frame.GetXaxis()->SetLabelSize(0.035);
  frame.GetXaxis()->SetTitleFont(42);
  //frame.GetXaxis()->SetNdivisions(204);
  frame.GetXaxis()->SetLimits(-5.0, 5.0);
  // Y-axis
  frame.GetYaxis()->CenterTitle(kFALSE);
  frame.GetYaxis()->SetTitleOffset(1.6);
  frame.GetYaxis()->SetTitleSize(0.042);
  frame.GetYaxis()->SetLabelSize(0.035);
  frame.GetYaxis()->SetTitleFont(42);
  setRange(frame);
  //frame.GetYaxis()->SetNdivisions(204);
};


void setRange(RooPlot& frame)
{
  // Find maximum and minimum points of Plot to rescale Y axis
  const RooHist gr = *((RooHist*)frame.findObject(0, RooHist::Class()));
  Double_t YMax = -1e99 , YMin = 1e99;
  for (int i=0; i<=gr.GetN(); i++) { double x, y; gr.GetPoint(i, x, y); if (YMin>y) { YMin = y; }; if (YMax<y) { YMax = y; } }
  const double Yup = YMax + (YMax-YMin)*1.5;
  frame.GetYaxis()->SetRangeUser(YMin, Yup);
  return;
};


void printElectroWeakBinning(TCanvas& pad, const RooWorkspace& ws, const std::string& dsName)
{
  //
  char chgL = ' '; if (dsName.find("Pl")!=std::string::npos) { chgL = '+'; } else if (dsName.find("Mi")!=std::string::npos) { chgL = '-'; }
  const std::string text = Form("#font[62]{#scale[1.1]{W^{%c}#rightarrow#mu^{%c}+x}}", chgL, chgL);
  //
  pad.cd();
  float xPos = 0.2, yPos = 0.89, dYPos = 0.045, dy = 0.025;
  TLatex t = TLatex(); t.SetNDC(); t.SetTextSize(0.025);
  t.DrawLatex(xPos, yPos-dy, Form("%s", text.c_str())); dy+=dYPos;
  auto parIt = std::unique_ptr<TIterator>(((RooDataSet*)ws.data(dsName.c_str()))->get()->createIterator());
  for (RooRealVar* it = (RooRealVar*)parIt->Next(); it!=NULL; it = (RooRealVar*)parIt->Next() ) {
    if (std::string(it->GetName())=="MET" || std::string(it->GetName())=="Muon_Pt") continue;
    const std::string varName = it->GetName();
    double defaultMin = 0.0 , defaultMax = 100000.0;
    if (varName=="Muon_Eta") { defaultMin = -2.5; defaultMax = 2.5; }
    if (varName=="Muon_Iso") { defaultMin = -1.0; }
    if (ws.var(varName.c_str())) {
      double minVal = ws.var(varName.c_str())->getMin();
      double maxVal = ws.var(varName.c_str())->getMax();
      string fVarName = varEWQLabel.at(varName);
      const bool ispPb = ( dsName.find("_pPb")!=std::string::npos || dsName.find("_PA")!=std::string::npos );
      if (varName=="Muon_Eta" && ws.var("useEtaCM")!=NULL) {
        minVal = PA::EtaLABtoCM(ws.var("Muon_Eta")->getMin(), ispPb);
        maxVal = PA::EtaLABtoCM(ws.var("Muon_Eta")->getMax(), ispPb);
        fVarName = varEWQLabel.at("Muon_EtaCM");
      }
      //
      if (minVal!=defaultMin && maxVal==defaultMax) {
        t.DrawLatex(xPos, yPos-dy, Form("%.2f #leq %s", minVal, fVarName.c_str())); dy+=dYPos;
      }
      if (minVal==defaultMin && maxVal!=defaultMax) {
        t.DrawLatex(xPos, yPos-dy, Form("%s < %.2f", fVarName.c_str(), maxVal)); dy+=dYPos;
      }
      if (minVal!=defaultMin && maxVal!=defaultMax) {
        t.DrawLatex(xPos, yPos-dy, Form("%.2f #leq %s < %.2f", minVal, fVarName.c_str(), maxVal)); dy+=dYPos;
      }
    }
  }
  pad.Update();
  return;
};


RooRealVar extractPar(RooPlot& frame)
{
  const RooHist gr = *((RooHist*)frame.findObject(0, RooHist::Class()));
  double xMin, xMax, yDummy;  gr.GetPoint(0, xMin, yDummy);  gr.GetPoint(gr.GetN()-1, xMax, yDummy);
  TH1D hist("hist", "hist", gr.GetN(), xMin, xMax);
  for (int i = 1; i <= gr.GetN(); i++) {
    double x , y, yErr; gr.GetPoint(i-1, x, y); yErr = gr.GetErrorY(i-1);
    hist.SetBinContent(i, y);
    hist.SetBinError(i, yErr);
  }
  //
  RooRealVar var("var", "var", hist.GetMean(), "");
  var.setError(hist.GetRMS());
  return var;
};


RooRealVar extractPar(RooPlot& frame, const RooRealVar& parIni, const RooDataSet& dataIni)
{
  RooRealVar par = parIni;
  RooRealVar mean("mean", "Mean of pull", par.getVal(), -100000., 100000.);
  RooRealVar sigma("sigma", "Width of pull", par.getError(), 0.00001, 200.);
  RooGaussian gauss("gauss","Gaussian of pull", par, mean, sigma);
  auto data = std::unique_ptr<RooDataSet>(new RooDataSet(dataIni, "dataIni"));
  gauss.fitTo(*data, RooFit::Minos(0), RooFit::PrintLevel(-1));
  RooRealVar var("var", "var", mean.getVal(), "");
  var.setError(sigma.getVal());
  return var;
};




#endif // #ifndef checkStatUnc_C

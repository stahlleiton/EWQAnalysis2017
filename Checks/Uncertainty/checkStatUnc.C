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
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
// RooFit headers
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooAbsPdf.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooMCStudy.h"
#include "RooStringVar.h"
// c++ headers
#include <dirent.h>
#include <iostream>
#include <string>


typedef  std::map< std::string , std::map< std::string , std::vector< std::vector< double > > > >  ResultMap;


void        plotResult(const ResultMap& result, const std::string& colTag, const std::string& dsTag, const std::string& outputDir);
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
  const std::string outputDir = CWD+"/MCStatStudy/" + workDirName+"/" + metTag+"/" + dsTag+"/" + colTag+"/";
  //
  // --------------------------------------------------------------------------------- //
  //
  // Get the list of input files
  //
  std::vector< std::string > inputFileNames;
  std::string preCWD  = CWD; preCWD.erase(preCWD.find_last_of("/"), 100);
  std::string pre2CWD = preCWD; pre2CWD.erase(pre2CWD.find_last_of("/"), 100);
  const std::string inputDirPath = Form("%s/Fitter/Output/%s/%s/%s/W/%s/result", pre2CWD.c_str(), workDirName.c_str(), metTag.c_str(), dsTag.c_str(), colTag.c_str());
  if (!fileList(inputFileNames, inputDirPath)) { return false; };
  //
  // --------------------------------------------------------------------------------- //
  //
  // Initialize the Output Graphs
  ResultMap result;
  //
  // Loop over the input files
  //
  //
  for (const auto& inputFileName : inputFileNames) {
    //
    std::cout << "[INFO] Processing file: " << inputFileName << std::endl;
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
    const std::string DSTAG = (ws->obj("DSTAG")    ) ? ((RooStringVar*)ws->obj("DSTAG")    )->getVal() : "";
    const std::string CHA   = (ws->obj("channel")  ) ? ((RooStringVar*)ws->obj("channel")  )->getVal() : "";
    const std::string COL   = (ws->obj("fitSystem")) ? ((RooStringVar*)ws->obj("fitSystem"))->getVal() : "";
    const std::string CHG   = (ws->obj("fitCharge")) ? ((RooStringVar*)ws->obj("fitCharge"))->getVal() : "";
    const std::string OBJ   = (ws->obj("fitObject")) ? ((RooStringVar*)ws->obj("fitObject"))->getVal() : "";
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
    auto pdfP = (RooAddPdf*) ws->pdf(Form("%s", pdfName.c_str()));
    if (pdfP == NULL) { std::cout << "[ERROR] The PDF " << pdfName << " was not found in the workspace" << std::endl; inputFile.Close(); return false; }
    auto pdf = *pdfP;
    //
    auto varP = (RooRealVar*) ws->var("MET");
    if (varP == NULL) { std::cout << "[ERROR] The variable MET was not found in the workspace" << std::endl; inputFile.Close(); return false; }
    auto var = *varP;
    var.setBins(100);
    //
    auto parP = (RooRealVar*) ws->var(Form("N_%s", tag.c_str()));
    if (parP == NULL) { std::cout << "[ERROR] The parameter " << Form("N_%s", tag.c_str()) << " was not found in the workspace" << std::endl; inputFile.Close(); return false; }
    auto par = *parP;
    //
    bool ispPb = false; if (COL=="PA" || COL=="pPb") { ispPb = true; }
    double etaMin = ws->var("Muon_Eta")->getMin(); if (useEtaCM) { etaMin = PA::EtaLABtoCM(etaMin, ispPb); }
    double etaMax = ws->var("Muon_Eta")->getMax(); if (useEtaCM) { etaMax = PA::EtaLABtoCM(etaMax, ispPb); }
    //
    RooMCStudy mcstudy(pdf, var, RooFit::Binned(kTRUE), RooFit::Silence(), RooFit::Extended(), RooFit::FitOptions(RooFit::Save(), RooFit::PrintLevel(-1)));
    //
    // ---------------------------------------------
    // Generate and fit 1000 samples of Poisson(nExpected) events
    mcstudy.generateAndFit(10000, ds->sumEntries());
    //
    // ------------------------------------------------
    // Make plots of the distributions of Signal Yield, the error on Signal Yield and the pull of Signal Yield
    auto framePull = mcstudy.plotPull(par, RooFit::Bins(50), RooFit::Range(-5.0, 5.0), RooFit::FitGauss(kTRUE));
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
    t.DrawLatex(xPos, yPos-dy, Form("N_{FIT} = %.0f #pm %.0f", par.getVal(), par.getError())); dy+=dYPos;
    auto framePar  = mcstudy.plotParam(par);
    auto toyVal = extractPar(*framePar);
    t.DrawLatex(xPos, yPos-dy, Form("#mu_{TOY} = %.0f #pm %.0f", toyVal.getVal(), toyVal.getError())); dy+=dYPos;
    auto frameErr  = mcstudy.plotError(par);
    auto toyErr = extractPar(*frameErr);
    t.DrawLatex(xPos, yPos-dy, Form("#sigma_{TOY} = %.0f #pm %.0f", toyErr.getVal(), toyErr.getError())); dy+=dYPos;
    auto toyPull = extractPar(*framePull);
    //
    printElectroWeakBinning(*cFig, *ws, dsName);
    //
    int lumiId = 0;
    if (COL=="pPb") { lumiId = 115; } else if (COL=="Pbp") { lumiId = 116; } else if (COL=="PA") { lumiId = 117; }
    CMS_lumi(cFig.get(), lumiId, 33, "");
    //
    // Save the plot in different formats
    const std::string fileName = Form("PULL_%s_%s_%s_%s_%.0f_%.0f", "MET",
                                      tag.c_str(),
                                      dsTag.c_str(),
                                      ( ws->var("useEtaCM") ? "MuEtaCM" : "MuEta" ),
                                      (etaMin*100.0), (etaMax*100.0)
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
    // Add to the result container
    if (std::abs(etaMax - etaMin)<2.0) {
      result[CHG]["FIT"].push_back ( std::vector<double>({ ((etaMax + etaMin)/2.) , par.getVal()     , ((etaMax - etaMin)/2.) , par.getError()     }) );
      result[CHG]["TOY"].push_back ( std::vector<double>({ ((etaMax + etaMin)/2.) , toyVal.getVal()  , ((etaMax - etaMin)/2.) , toyVal.getError()  }) );
      result[CHG]["PULL"].push_back( std::vector<double>({ ((etaMax + etaMin)/2.) , toyPull.getVal() , ((etaMax - etaMin)/2.) , toyPull.getError() }) );
    }
    //
  }
  //
  plotResult(result, colTag, dsTag, outputDir);
  //
  // Return
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
  const std::string text = Form("#font[62]{#scale[1.1]{W^{%c}#rightarrow#mu^{%c}+#nu_{#mu}}}", chgL, chgL);
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


void updateYAxisRange(TGraphErrors& graph)
{
  // Find the max and min of graph
  double yMin = 9999999999. , yMax = -9999999999.;
  for (int i = 0; i < graph.GetN(); i++) {
    double x, y, yL, yH;
    graph.GetPoint(i, x, y);
    yL = y - graph.GetErrorY(i);
    yH = y + graph.GetErrorY(i);
    if (yH > yMax) { yMax = yH; }
    if (yL < yMin) { yMin = yL; }
  }
  // Set the up and down of y axis
  double fMax = 0.6 , fMin = 0.1;
  double yDown , yUp;
  yDown = (fMax*yMin - fMin*yMax)/(fMax - fMin);
  yUp   = yDown + (yMax-yMin)/(fMax-fMin);
  // Update the y range
  graph.GetYaxis()->SetRangeUser(yDown, yUp);
};


void plotResult(const ResultMap& result, const std::string& colTag, const std::string& dsTag, const std::string& outputDir)
{
  //
  // Fill the output graphs
  std::map< std::string , std::map< std::string , TGraphErrors > > graph;
  //
  for (const auto& ch : result) {
    for (const auto& t : ch.second) {
      const auto& b = t.second;
      auto& gr = graph[ch.first][t.first];
      gr.Set(b.size());
      for (uint i = 0; i < b.size(); i++) {
        gr.SetPoint(i, b[i][0], b[i][1]);
        gr.SetPointError(i, b[i][2], b[i][3]);
      }
      //
      gr.SetTitle("");
      gr.GetXaxis()->SetTitle("");
      gr.SetMarkerColor(kBlack);
      gr.SetMarkerStyle(kFullCircle);
      if (colTag=="PA" || colTag=="pPb") { gr.GetXaxis()->SetLimits(-3.0, 2.2); }
      if (colTag=="Pbp") { gr.GetXaxis()->SetLimits(-2.2, 3.0); }
      gr.GetXaxis()->SetTitleFont(42);
      gr.GetYaxis()->SetTitleFont(42);
      if (t.first=="PULL") {
        gr.GetXaxis()->SetTitle("Leading Muon #eta_{CM}");
        gr.GetXaxis()->SetTitleOffset(1);
        gr.GetXaxis()->SetTitleSize(0.16);
        gr.GetXaxis()->SetLabelSize(0.14);
        gr.GetYaxis()->SetTitle("Pull");
        gr.GetYaxis()->SetTitleOffset(0.4);
        gr.GetYaxis()->SetTitleSize(0.16);
        gr.GetYaxis()->SetLabelSize(0.1);
        gr.GetYaxis()->SetRangeUser(-3.0, 3.0);
      }
      else {
        gr.GetXaxis()->SetTitleSize(0.05);
        gr.GetXaxis()->SetTitleOffset(3);
        gr.GetXaxis()->SetLabelOffset(3);
        gr.GetYaxis()->SetTitle(Form("N_WToMu%s_%s", ch.first.c_str(), colTag.c_str()));
        gr.GetYaxis()->SetLabelSize(0.044);
        gr.GetYaxis()->SetTitleSize(0.044);
        gr.GetYaxis()->SetTitleOffset(1.7);
        updateYAxisRange(gr);
      }
      if (t.first!="TOY") { gr.SetMarkerColor(kBlue); }      
    }
    // Draw the pads
    TCanvas c("c", "c", 1000, 1000);
    // Define the plotting pads
    TPad* pad1 = new TPad("pad1", "", 0, 0.23, 1, 1);
    TPad* pad2 = new TPad("pad2", "", 0, 0, 1, 0.228);
    // Format the pads
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.4);
    pad2->SetFillStyle(4000); 
    pad2->SetFrameFillStyle(4000);
    pad1->SetBottomMargin(0.015);
    // Main Frame
    c.cd();
    pad1->Draw();
    pad1->cd();
    graph.at(ch.first).at("FIT").Draw("AP");
    graph.at(ch.first).at("TOY").Draw("P same");
    //
    TLatex t = TLatex(); t.SetNDC(); t.SetTextSize(0.045);
    char chgL = ' '; if (ch.first=="Pl") { chgL = '+'; } else if (ch.first=="Mi") { chgL = '-'; }
    t.DrawLatex(0.23, 0.80, Form("#font[62]{#scale[1.1]{W^{%c}#rightarrow#mu^{%c}+#nu_{#mu}}}", chgL, chgL));
    //
    TLegend leg(0.52, 0.74, 0.52+0.2, 0.74+0.10);
    leg.AddEntry(&graph.at(ch.first).at("FIT"), "Fitted", "p")->SetTextSize(0.040);
    leg.AddEntry(&graph.at(ch.first).at("TOY"), "Toy", "p")->SetTextSize(0.040);
    leg.Draw("same");
    //
    int lumiId = 0;
    if (colTag=="pPb") { lumiId = 115; } else if (colTag=="Pbp") { lumiId = 116; } else if (colTag=="PA") { lumiId = 117; }
    CMS_lumi(pad1, lumiId, 33, "");
    gStyle->SetTitleFontSize(0.05);
    pad1->Update();
    // Pull Frame
    c.cd();
    pad2->Draw();
    pad2->cd();
    graph.at(ch.first).at("PULL").Draw("AP");
    TLine line(graph.at(ch.first).at("PULL").GetXaxis()->GetXmin(), 0.0, graph.at(ch.first).at("PULL").GetXaxis()->GetXmax(), 0.0);
    line.SetLineColor(kBlue); line.SetLineWidth(3); line.SetLineStyle(2); line.Draw("same");
    pad2->Update();
    // Save the plot in different formats
    const std::string fileName = Form("RESULT_%s_%s_WToMu%s_%s", "MET", dsTag.c_str(), ch.first.c_str(), colTag.c_str());
    //
    c.SaveAs(Form("%splot/png/%s.png", outputDir.c_str(), fileName.c_str()));
    c.SaveAs(Form("%splot/pdf/%s.pdf", outputDir.c_str(), fileName.c_str()));
    c.Clear();
    c.Close();
  }
};


#endif // #ifndef checkStatUnc_C

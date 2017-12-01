//================================================================================================
//
// Perform fits to components of recoil from Z->mumu events
//
//  * Outputs a ROOT file of functions parametrizing the distribution of recoil components
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
// Auxiliary Headers
#include "infoRecoilPDF.h"
#include "../../Utilities/CMS/tdrstyle.C"
#include "../../Utilities/CMS/CMS_lumi.C"
// ROOT Headers
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>                    // file handle class
#include <TF1.h>                      // 1D function
#include <TH1.h>                      // 1D function
#include <TFitResult.h>               // class to handle fit results
#include <TGraphAsymmErrors.h>        // graph class
#include <TCanvas.h>                  // canvas
#include <TAxis.h>
#include <TPaveText.h>
#include <TKey.h>
#include <TClass.h>
#include <TVirtualFitter.h>
// C++ Headers
#include <iostream>                   // standard I/O
#include <fstream>                    // standard I/O
#include <sys/stat.h>

#endif

//=== FUNCTION DECLARATIONS ======================================================================================

// Axis and Text Utility functions
void updateYAxisRange ( std::unique_ptr<TGraphAsymmErrors>& graph );

std::string formatText(const std::string& text);


// generate web page
void makeHTML(const std::string outDir,
              const std::map< std::string , std::unique_ptr<TGraphAsymmErrors> >& u1Graph,
              const std::map< std::string , std::unique_ptr<TGraphAsymmErrors> >& u2Graph,
              const std::string uparName, const std::string uprpName);

inline bool fileExist (const std::string& name);

//--------------------------------------------------------------------------------------------------
// perform fit of recoil component
bool performFit(
                std::map< std::string , std::unique_ptr<TF1> >& uFit,
                std::map< std::string , TFitResultPtr >& uFitRes,
                TCanvas& c,
                std::map< std::string , std::unique_ptr<TGraphAsymmErrors> >& uGraph,
                const std::map< std::string , FcnInfo >& uFcn,
                const std::string uName,
                const std::string met,
                const std::string col,
                const std::string outputDir
                );


//=== MAIN MACRO ================================================================================================= 

void getRecoilPDF(
                  const bool isData = false,
                  uint pfumodel = 2, // u1 model (1 => single Gaussian, 2 => double Gaussian, 3 => triple Gaussian)
                  const bool applyHFCorr = true, // Only used for MC, in data it is ignored
                  const std::vector< std::string > metType = { "PF_RAW" /*, "PF_Type1" , "PF_NoHF_RAW" , "PF_NoHF_Type1" */},
                  const std::vector< std::string > COLL    = { "PA" /*, "Pbp" , "pPb" */}
                  )
{
  //
  const std::string uparName = "u1";
  const std::string uprpName = "u2";
  //
  std::string pfu12model = "singleGauss"; // Same model for u1 and u2 is used as set in fitRecoil.C
  if (pfumodel>3 || pfumodel<1){
    std::cout << "[ERROR] The supported models are: 1 => single Gaussian, 2 => double Gaussian, 3 => triple Gaussian" << std::endl;
    return;
  }
  if (pfumodel>1){
    if (pfumodel>2) pfu12model = "tripleGauss";
    else pfu12model = "doubleGauss";
  }
  //
  // Change the working directory
  const std::string CWD = getcwd(NULL, 0);
  const std::string mainDir = Form("%s/FitRecoil/", CWD.c_str());
  gSystem->mkdir(Form("%s/FitRecoil", CWD.c_str()), kTRUE);
  gSystem->ChangeDirectory(Form("%s/FitRecoil", CWD.c_str()));

  std::string dsLabel;
  if (isData) { dsLabel = "DATA"; }
  else { dsLabel = "MC_DYToMuMu_POWHEG"; }

  // Initialize Canvas
  setTDRStyle();
  auto c = std::unique_ptr<TCanvas>(new TCanvas("c", "c", 800, 800));

  // Step 1: Find all the input root files (Recoil Fit Graphs)
  for (const auto& met : metType) {
    for (const auto& col : COLL) {
      std::cout << "[INFO] Working with MET " << met << " and coll: " << col  << std::endl;
      //
      // Define directories and file names
      const std::string inputDir  = mainDir + dsLabel + "/" + ("MET_"+met) + "/" + col + "/" + (isData?"":(applyHFCorr?"HFCorr/":"noHFCorr/")) + pfu12model.c_str() + "/Fits/";
      const std::string outputDir = mainDir + dsLabel + "/" + ("MET_"+met) + "/" + col + "/" + (isData?"":(applyHFCorr?"HFCorr/":"noHFCorr/")) + pfu12model.c_str() + "/Results/";
      gSystem->mkdir(outputDir.c_str(), true);
      const std::string fileName  = "plots_RecoilPDF_" + met + "_" + col + ".root";
      //
      // Open the input file
      auto file = std::unique_ptr<TFile>(TFile::Open((inputDir + fileName).c_str(), "READ"));
      std::cout << "[INFO] Reading file: " << (inputDir + fileName) << std::endl;
      if (!file) { std::cout << "[ERROR] Input file " << (inputDir + fileName) << " was not found!" << std::endl; return; }
      //
      std::map< std::string , std::unique_ptr<TGraphAsymmErrors> > u1Graph;
      std::map< std::string , std::unique_ptr<TGraphAsymmErrors> > u2Graph;
      //
      // Loop over all TGraphs inside the file
      TIter next(file->GetListOfKeys());
      for (TKey* key = (TKey*)next(); key!=NULL; key = (TKey*)next() ) {
        TClass *cl = gROOT->GetClass(key->GetClassName());
        if (!cl->InheritsFrom("TGraphAsymmErrors")) continue;
        TGraphAsymmErrors *gr = (TGraphAsymmErrors*)key->ReadObj();
        std::string name = gr->GetName();
        std::cout << "[INFO] Importing Graph: " << name << std::endl;
        if (name.find("PFu1")!=std::string::npos) { u1Graph[name.substr(name.find("u1")+2)].reset(gr); }
        if (name.find("PFu2")!=std::string::npos) { u2Graph[name.substr(name.find("u2")+2)].reset(gr); }
      }
      //
      // Fit the Recoil graphs
      //
      // For u1 component
      std::map< std::string , std::unique_ptr<TF1> > u1Fit;
      std::map< std::string , TFitResultPtr > u1FitRes;
      if (!( performFit(u1Fit, u1FitRes, *c, u1Graph, u1FcnInfo.at(dsLabel).at(met).at("PA"), uparName, met, col, outputDir) )) { return; }
      // For u2 component
      std::map< std::string , std::unique_ptr<TF1> > u2Fit;
      std::map< std::string , TFitResultPtr > u2FitRes;
      if (!( performFit(u2Fit, u2FitRes, *c, u2Graph, u2FcnInfo.at(dsLabel).at(met).at("PA"), uprpName, met, col, outputDir) )) { return; }
      //
      // Save the fit results
      //
      const std::string outfname = (outputDir + Form("fits_RecoilPDF_%s_%s.root", met.c_str(), col.c_str()));
      auto outfile = std::unique_ptr<TFile>(new TFile(outfname.c_str(), "RECREATE"));
      for (const auto& graph  : u1Graph ) { if (graph.second ) graph.second->Write();  }
      for (const auto& graph  : u2Graph ) { if (graph.second ) graph.second->Write();  }
      for (const auto& fit    : u1Fit   ) { if (fit.second   ) fit.second->Write();    }
      for (const auto& fit    : u2Fit   ) { if (fit.second   ) fit.second->Write();    }
      for (const auto& fitres : u1FitRes) { if (fitres.second) fitres.second->Write(); }
      for (const auto& fitres : u2FitRes) { if (fitres.second) fitres.second->Write(); }
      outfile->Close();
      // Make Website with plots
      const std::string htmlDir = mainDir + dsLabel +"/"+ ("MET_"+met) +"/"+ col + "/" + (isData?"":(applyHFCorr?"HFCorr/":"noHFCorr/")) + pfu12model.c_str();
      makeHTML(htmlDir, u1Graph, u2Graph, uparName, uprpName);
      cout << "  <> Output saved in " << outputDir << endl;
      // Clean up
      file->Close();
    }
  }
  c->Close();
  return;
}
  
//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
bool performFit(
                std::map< std::string , std::unique_ptr<TF1> >& uFit,
                std::map< std::string , TFitResultPtr >& uFitRes,
                TCanvas& c,
                std::map< std::string , std::unique_ptr<TGraphAsymmErrors> >& uGraph,
                const std::map< std::string , FcnInfo >& uFcn,
                const std::string uName,
                const std::string met,
                const std::string col,
                const std::string outputDir
                )
{
  //
  // Create output directories
  //
  const std::string outDir = outputDir + uName + "/";
  gSystem->mkdir( (outDir + "png/").c_str(), true);
  gSystem->mkdir( (outDir + "pdf/").c_str(), true);
  //
  // Proceed to fit the recoil components
  //
  // For Error Band
  std::map< std::string , std::unique_ptr<TH1D> > uConfInt;
  //
  for (auto& graph : uGraph) {
    std::string varName = "";
    if (uFcn.count(graph.first)>0) { varName = graph.first; }
    else {
      // Default Fit Functions
      if (graph.first.find("mean")!=std::string::npos ) { varName = "mean1";  }
      if (graph.first.find("sigma")!=std::string::npos) { varName = "sigma1"; }
    }
    if (graph.first.find("dmean")!=std::string::npos || graph.first.find("rsigma")!=std::string::npos || graph.first.find("frac")!=std::string::npos) continue;
    if (uFcn.count(varName)==0) { std::cout << "[ERROR] Recoil parameter: " << graph.first << " was not found!" << std::endl; return false; }
    FcnInfo fcn = uFcn.at(varName);
    uFit[graph.first] = std::unique_ptr<TF1>(new TF1(Form("fcnPF%s%s", uName.c_str(), graph.first.c_str()), fcn.exp.c_str(), 0.0, 300.));
    if (uFit[graph.first]->GetNpar()!=int(fcn.par.size())) { std::cout << "[ERROR] Number of input parameters used for fitting: " << varName << " is inconsistent!" << std::endl; return false; }
    for (uint i=0; i<fcn.par.size(); i++) { uFit.at(graph.first)->SetParameter(i, fcn.par[i]); }
    for (uint i=0; i<fcn.par.size(); i++) { uFit.at(graph.first)->SetParName(i, FcnParName.at(fcn.exp)[i].c_str()); }
    for (uint i=0; i<fcn.min.size(); i++) { uFit.at(graph.first)->SetParLimits(i, fcn.min[i], fcn.max[i]); }
    graph.second->Fit(uFit.at(graph.first)->GetName(), "QMRN0W");
    uFitRes[graph.first] = graph.second->Fit(uFit.at(graph.first)->GetName(), "QMRN0SE");
    uFitRes.at(graph.first)->SetName(Form("fitresPF%s%s", uName.c_str(), graph.first.c_str()));
    /*Create a histogram to hold the confidence intervals*/
    uConfInt[graph.first].reset(new TH1D(Form("confintPF%s%s", uName.c_str(), graph.first.c_str()), "", graph.second->GetN(), graph.second->GetXaxis()->GetXmin(), graph.second->GetXaxis()->GetXmax()));
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(uConfInt.at(graph.first).get());
  }
  //
  // Compute the pull distribution
  //
  std::map< std::string , std::unique_ptr<TGraphAsymmErrors> > pull;
  for (auto& fit : uFit) {
    pull[fit.first].reset((TGraphAsymmErrors*)uGraph.at(fit.first)->Clone());
    pull.at(fit.first)->SetName(Form("pullPF%s%s", uName.c_str(), fit.first.c_str()));
    for (int i = 0; i < uGraph.at(fit.first)->GetN(); i++) {
      double xVal, yVal;
      uGraph.at(fit.first)->GetPoint(i, xVal, yVal);
      const double xErrLo = uGraph.at(fit.first)->GetErrorXlow(i);
      const double xErrHi = uGraph.at(fit.first)->GetErrorXhigh(i);
      const double yErrLo = uGraph.at(fit.first)->GetErrorYlow(i);
      const double yErrHi = uGraph.at(fit.first)->GetErrorYhigh(i);
      const double pull_xVal = xVal;
      const double pull_xErrLo = xErrLo;
      const double pull_xErrHi = xErrHi;
      const double yFitVal = fit.second->Eval(xVal);
      const double norm = ( (yVal>0) ? yErrLo : yErrHi );
      const double pull_yVal   = ( (norm>0.) ? ( (yVal - yFitVal) / norm ) : 0.0 );
      const double pull_yErrLo = ( (norm>0.) ? ( yErrLo / norm ) : 0.0 );
      const double pull_yErrHi = ( (norm>0.) ? ( yErrHi / norm ) : 0.0 );
      pull.at(fit.first)->SetPoint(i, pull_xVal, pull_yVal);
      pull.at(fit.first)->SetPointError(i, pull_xErrLo, pull_xErrHi, pull_yErrLo, pull_yErrHi);
    }
  }
  //
  // Plot the fit results
  //
  for (auto& graph : uGraph) {
    if (uFit.count(graph.first)==0) continue;
    // Create the text labels
    std::cout << "[INFO] Creating the text labels for the plots" << std::endl;
    const std::string xlabel = Form("p_{T}(ll) [GeV/c]");
    const std::string ylabel   = Form("%s(%s) [GeV/c]", formatText(graph.first).c_str(), formatText(uName).c_str());
    std::vector< std::string > parText;
    auto& fit = uFit.at(graph.first);
    for (int i = 0; i < fit->GetNpar(); i++) {
      double value = fit->GetParameter(i);
      double error = fit->GetParError(i);
      std::string name = fit->GetParName(i);
      parText.push_back(Form("%s = %.2f #pm %.2f" , name.c_str(), value, error));
      double min, max; fit->GetParLimits(i, min, max);
      if ( ( (std::abs(value - min)/error) < 3.0 ) || ( (std::abs(value - max)/error) < 3.0 ) ) { parText[i] += " (!)"; }
    }
    auto tb = std::unique_ptr<TPaveText>(new TPaveText(0.17, (0.88-0.04*(4+parText.size())), 0.45, 0.88, "NDC"));
    tb->SetTextColor(kBlack);
    tb->SetFillStyle(0);
    tb->SetBorderSize(0);
    tb->AddText(met.c_str());
    tb->AddText("   ");
    for (const auto& par : parText) { tb->AddText(par.c_str()); }
    tb->AddText("   ");
    tb->AddText(FcnExpText.at(fit->GetFormula()->GetTitle()).c_str());
    //
    // Format the graphs
    std::cout << "[INFO] Format the frames" << std::endl;
    //  Main Graph
    graph.second->SetTitle("");
    graph.second->GetXaxis()->SetTitle("");
    graph.second->SetMarkerColor(kBlack);
    graph.second->SetMarkerStyle(kFullCircle);
    graph.second->GetXaxis()->SetTitleSize(0.05);
    graph.second->GetXaxis()->SetTitleFont(42);
    graph.second->GetXaxis()->SetTitleOffset(3);
    graph.second->GetXaxis()->SetLabelOffset(3);
    graph.second->GetXaxis()->SetRangeUser(0.0, 200.);
    graph.second->GetYaxis()->SetTitle(ylabel.c_str());
    graph.second->GetYaxis()->SetLabelSize(0.044);
    graph.second->GetYaxis()->SetTitleSize(0.044);
    graph.second->GetYaxis()->SetTitleOffset(1.7);
    graph.second->GetYaxis()->SetTitleFont(42);
    updateYAxisRange(graph.second);
    //  Pull Graph
    pull.at(graph.first)->SetTitle("");
    pull.at(graph.first)->SetMarkerColor(kBlack);
    pull.at(graph.first)->SetMarkerStyle(kFullCircle);
    pull.at(graph.first)->GetXaxis()->SetTitleOffset(1);
    pull.at(graph.first)->GetXaxis()->SetTitleSize(0.16);
    pull.at(graph.first)->GetXaxis()->SetLabelSize(0.14);
    pull.at(graph.first)->GetXaxis()->SetTitle(xlabel.c_str());
    pull.at(graph.first)->GetXaxis()->SetRangeUser(0.0, 200.);
    pull.at(graph.first)->GetYaxis()->CenterTitle(kTRUE);
    pull.at(graph.first)->GetYaxis()->SetTitleOffset(0.4);
    pull.at(graph.first)->GetYaxis()->SetTitleSize(0.16);
    pull.at(graph.first)->GetYaxis()->SetLabelSize(0.1);
    pull.at(graph.first)->GetYaxis()->SetTitle("Pull");
    pull.at(graph.first)->GetYaxis()->SetRangeUser(-7.0, 7.0);
    // Confidence Interval
    uConfInt.at(graph.first)->SetStats(kFALSE);
    uConfInt.at(graph.first)->SetFillColor(kYellow);
    uConfInt.at(graph.first)->SetFillStyle(3225);
    //  Fit Function
    uFit.at(graph.first)->SetLineColor(kRed);
    uFit.at(graph.first)->SetLineWidth(2);
    uFit.at(graph.first)->SetLineStyle(1);
    // Define the plotting pads
    TPad *pad1  = new TPad("pad1", "", 0, 0.23, 1, 1);
    TPad *pad2  = new TPad("pad2", "", 0, 0, 1, 0.228);
    auto  pline = std::unique_ptr<TLine>(new TLine(0.0, 0.0,  200.0, 0.0));
    pline->SetLineColor(kBlue); pline->SetLineWidth(3); pline->SetLineStyle(2);
    // Format the pads
    c.cd();
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.4);
    pad2->SetFillStyle(4000); 
    pad2->SetFrameFillStyle(4000);
    pad1->SetBottomMargin(0.015);
    // Draw the frams
    std::cout << "[INFO] Drawing the frames" << std::endl;
    // Main Frame
    pad1->Draw();
    pad1->cd();
    graph.second->Draw("AP");
    uConfInt.at(graph.first)->Draw("e3 same");
    uFit.at(graph.first)->Draw("same");
    tb->Draw("same");
    graph.second->Draw("P same");
    // Apply CMS style to pad
    std::cout << "[INFO] Setting the CMS style on the plot" << std::endl;
    int lumiId = 0;
    const bool isData = (outputDir.find("DATA")!=std::string::npos);
    if (isData ) { if (col=="pPb") { lumiId = 109; } else if (col=="Pbp") { lumiId = 110; } else if (col=="PA") { lumiId = 111; } }
    if (!isData) { if (col=="pPb") { lumiId = 112; } else if (col=="Pbp") { lumiId = 113; } else if (col=="PA") { lumiId = 114; } }
    CMS_lumi(pad1, lumiId, 33, "");
    gStyle->SetTitleFontSize(0.05);
    pad1->SetLogy(false);
    pad1->SetLogx(false);
    pad1->Update();
    // Pull Frame
    c.cd();
    pad2->Draw();
    pad2->cd();
    pull[graph.first]->Draw("AP");
    pline->Draw("same");
    auto t = std::unique_ptr<TLatex>(new TLatex()); t->SetNDC(); t->SetTextSize(0.12);
    t->DrawLatex(0.7, 0.85, Form("#chi^{2}/ndof = %.0f / %d ", uFit.at(graph.first)->GetChisquare(), uFit.at(graph.first)->GetNDF()));
    pad2->Update();
    // Save the pads
    //
    c.SaveAs( (outDir + "png/" + Form("fitPF%s%s.png", uName.c_str(), graph.first.c_str())).c_str() );
    c.SaveAs( (outDir + "pdf/" + Form("fitPF%s%s.pdf", uName.c_str(), graph.first.c_str())).c_str() );
    c.Clear();
  }
  // Return
  return true;
}


//--------------------------------------------------------------------------------------------------
void makeHTML(
              const std::string outDir,
              const std::map< std::string , std::unique_ptr<TGraphAsymmErrors> >& u1Graph,
              const std::map< std::string , std::unique_ptr<TGraphAsymmErrors> >& u2Graph,
              const std::string uparName   = "u1",
              const std::string uprpName   = "u2"
              )
{

  ofstream htmlfile;
  std::string htmlfname;
  htmlfname = Form("%s/plots.html", outDir.c_str());
  std::cout << "[INFO] Proceed to create HLTM file: " << htmlfname << std::endl;
  htmlfile.open(htmlfname);
  if (!htmlfile.is_open()) { std::cout << "[ERROR] Was not able to create file: " << htmlfname << std::endl; return; }
  htmlfile << "<!DOCTYPE html" << endl;
  htmlfile << "    PUBLIC \"-//W3C//DTD HTML 3.2//EN\">" << endl;
  htmlfile << "<html>" << endl;
  htmlfile << "<head><title>Recoil Fits</title></head>" << endl;
  htmlfile << "<body bgcolor=\"EEEEEE\">" << endl;

  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl;
  uint counter = 0;
  for (const auto& gr : u1Graph) {
    if (counter%5 == 0) { htmlfile << "<tr>" << endl; }
    if (fileExist(outDir+"/Results/"+uparName+"/png/fitPF"+uparName+gr.first+".png")) {
      htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"Results/" << uparName << "/png/fitPF" << uparName << gr.first << ".png\"><img src=\"Results/" << uparName << "/png/fitPF" << uparName << gr.first << ".png\" alt=\"Results/" << uparName << "/png/fitPF" << uparName << gr.first << ".png\" width=\"100%\"></a></td>" << endl;
    }
    else {
      htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"Fits/" << uparName << "/png/pf" << uparName << gr.first << ".png\"><img src=\"Fits/" << uparName << "/png/pf" << uparName << gr.first << ".png\" alt=\"Fits/" << uparName << "/png/pf" << uparName << gr.first << ".png\" width=\"100%\"></a></td>" << endl;
    }
    if (counter%5 == 4) { htmlfile << "</tr>" << endl; }
    counter++;
  }
  if (u1Graph.count("sigma2")==0) {
    htmlfile << "<td width=\"20%\"></td>" << endl;
    htmlfile << "<td width=\"20%\"></td>" << endl;
    htmlfile << "<td width=\"20%\"></td>" << endl;
  }
  htmlfile << "</tr>" << endl;
  
  htmlfile << "</table>" << endl;
  htmlfile << "PF " << uparName << " fits:";
  htmlfile << " <a target=\"_blank\" href=\"pf" << uparName << "fits.html\">[linear scale]</a>";
  htmlfile << " <a target=\"_blank\" href=\"pf" << uparName << "fitslog.html\">[log scale]</a>" << endl;
  htmlfile << "<hr />" << endl;

  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl;
  htmlfile << "<tr>" << endl;
  counter = 0;
  for (const auto& gr : u2Graph) {
    if (counter%5 == 0) { htmlfile << "<tr>" << endl; }
    if (fileExist(outDir+"/Results/"+uprpName+"/png/fitPF"+uprpName+gr.first+".png")) {
      htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"Results/" << uprpName << "/png/fitPF" << uprpName << gr.first << ".png\"><img src=\"Results/" << uprpName << "/png/fitPF" << uprpName << gr.first << ".png\" alt=\"Results/" << uprpName << "/png/fitPF" << uprpName << gr.first << ".png\" width=\"100%\"></a></td>" << endl;
    }
    else {
      htmlfile << "<td width=\"20%\"><a target=\"_blank\" href=\"Fits/" << uprpName << "/png/pf" << uprpName << gr.first << ".png\"><img src=\"Fits/" << uprpName << "/png/pf" << uprpName << gr.first << ".png\" alt=\"Fits/" << uprpName << "/png/pf" << uprpName << gr.first << ".png\" width=\"100%\"></a></td>" << endl;
    }
    if (counter%5 == 4) { htmlfile << "</tr>" << endl; }
    counter++;
  }
  if (u2Graph.count("sigma2")==0) {
    htmlfile << "<td width=\"20%\"></td>" << endl;
    htmlfile << "<td width=\"20%\"></td>" << endl;
    htmlfile << "<td width=\"20%\"></td>" << endl;
  }
  htmlfile << "</tr>" << endl;
  htmlfile << "</table>" << endl;
  htmlfile << "PF " << uprpName << " fits:";
  htmlfile << " <a target=\"_blank\" href=\"pf" << uprpName << "fits.html\">[linear scale]</a>";
  htmlfile << " <a target=\"_blank\" href=\"pf" << uprpName << "fitslog.html\">[log scale]</a>" << endl;
  htmlfile << "<hr />" << endl;
  
  htmlfile << "</body>" << endl;
  htmlfile << "</html>" << endl;
  htmlfile.close(); 
};


void updateYAxisRange(std::unique_ptr<TGraphAsymmErrors>& graph)
{
  if (graph==NULL) return;
  // Find the max and min of graph
  double yMin = 9999999999. , yMax = -9999999999.;
  for (int i = 0; i < graph->GetN(); i++) {
    double x, y, yL, yH;
    graph->GetPoint(i, x, y);
    yL = y - graph->GetErrorY(i);
    yH = y + graph->GetErrorY(i);
    if (yH > yMax) { yMax = yH; }
    if (yL < yMin) { yMin = yL; }
  }
  // Set the up and down of y axis
  double fMax = 0.6 , fMin = 0.1;
  double yDown , yUp;
  yDown = (fMax*yMin - fMin*yMax)/(fMax - fMin);
  yUp   = yDown + (yMax-yMin)/(fMax-fMin);
  // Update the y range
  graph->GetYaxis()->SetRangeUser(yDown, yUp);
};


std::string formatText(const std::string& text)
{
  std::string out = text;
  if (out=="mean")  { out = "#minus#mu";   return out; }
  if (out=="sigma") { out = "#sigma"; return out; }
  if (out.find("mean")!=std::string::npos)  { out.replace(out.find("mean")  , std::string("mean").length()  , "#minus#mu_{"); out += "}"; }
  if (out.find("sigma")!=std::string::npos) { out.replace(out.find("sigma") , std::string("sigma").length() , "#sigma_{");    out += "}"; }
  if (out.find("frac")!=std::string::npos)  { out.replace(out.find("frac")  , std::string("frac").length()  , "f_{");         out += "}"; }
  if (out.find("u1")!=std::string::npos)    { out.replace(out.find("u1")    , std::string("u1").length()    , "u_{1}"); }
  if (out.find("u2")!=std::string::npos)    { out.replace(out.find("u2")    , std::string("u2").length()    , "u_{2}"); }
  return out;
};


inline bool fileExist (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
};

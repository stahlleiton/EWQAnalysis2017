#ifndef drawElectroWeakMETPlot_C
#define drawElectroWeakMETPlot_C

#include "../Utilities/initClasses.h"
#include "TGaxis.h"

void       printElectroWeakMETParameters ( TPad* pad, const RooWorkspace& ws, const std::string& pdfName, const uint& drawMode );
void       printElectroWeakBinning       ( TPad* pad, const RooWorkspace& ws, const std::string& dsName, const std::vector< std::string >& text, const uint& drawMode );
void       printElectroWeakLegend        ( TPad* pad, const RooPlot& frame, const StrMapMap& legInfo, const uint& drawMode );
RooHist*   makeRatioHist ( RooPlot* frame, const char* histname, const char* curvename, bool normalize, bool useAverage );
bool       getVar        ( std::vector<RooRealVar*>& varVec, const RooWorkspace& ws, const std::string& name, const std::string& pdfName );
void       parseVarName  ( const std::string& name, std::string& label );

bool drawElectroWeakMETPlot( RooWorkspace& ws,            // Local Workspace
                             // Select the type of datasets to fit
                             const std::string& fileName,
                             const std::string& outputDir,
                             const int& nBins
                             )
{
  // set the CMS style
  setTDRStyle();

  const std::string DSTAG = (ws.obj("DSTAG"))     ? ((TObjString*)ws.obj("DSTAG"))->GetString().Data()     : "";
  const std::string cha   = (ws.obj("channel"))   ? ((TObjString*)ws.obj("channel"))->GetString().Data()   : "";
  const std::string col   = (ws.obj("fitSystem")) ? ((TObjString*)ws.obj("fitSystem"))->GetString().Data() : "";
  const std::string chg   = (ws.obj("fitCharge")) ? ((TObjString*)ws.obj("fitCharge"))->GetString().Data() : "";
  const std::string obj   = (ws.obj("fitObject")) ? ((TObjString*)ws.obj("fitObject"))->GetString().Data() : "";

  std::string tag = ( obj + cha + chg + "_" + col );
  std::string dsName = ( "d" + chg + "_" + DSTAG );
  std::string pdfName = Form("pdfMET_Tot%s", tag.c_str());
  bool paperStyle = false;
  bool setLogScale = false;
  std::vector< double > range = { ws.var("MET")->getMin(), ws.var("MET")->getMax() };

  bool isMC = (DSTAG.find("MC")!=std::string::npos);
  bool isWeighted = ws.data(dsName.c_str())->isWeighted();
  int drawMode = 0;
  
  // Format Object name
  std::string process = "";
  char chgL = ' '; if (chg=="Pl") { chgL = '+'; } else if (chg=="Mi") { chgL = '-'; }
  if (obj=="WToTau") { process = Form("W^{%c}#rightarrow#tau^{%c}", chgL, chgL); } else if (obj=="W") { process = Form("W^{%c}", chgL);; } 
  else if (obj=="DYZ") { process = "Z/#gamma"; } else if (obj=="QCD") { process = "QCD"; }
  if (cha=="ToMu") { process += Form("#rightarrow#mu^{%c}+x", chgL); }
  process = Form("#font[62]{#scale[1.1]{%s}}", process.c_str());
  
  StrMapMap legInfo;

  std::map< std::string , RooPlot* > frame;
  std::map< std::string , TPad* > pad;

  // Create the main plot of the fit
  frame["MAIN"] = ws.var("MET")->frame( RooFit::Bins(nBins), RooFit::Range(range[0], range[1]) );
  if (ws.data(("CutAndCount_"+dsName).c_str())) {
    ws.data(("CutAndCount_"+dsName).c_str())->plotOn(frame["MAIN"], RooFit::Name(Form("plot_Tot%s", dsName.c_str())), RooFit::DataError(RooAbsData::SumW2), 
                                                     RooFit::XErrorSize(0), RooFit::MarkerColor(kBlack), RooFit::LineColor(kBlack), RooFit::MarkerSize(1.2));
  }
  else {
    ws.data(dsName.c_str())->plotOn(frame["MAIN"], RooFit::Name(Form("plot_Tot%s", dsName.c_str())), RooFit::DataError(RooAbsData::SumW2), 
                                    RooFit::XErrorSize(0), RooFit::MarkerColor(kBlack), RooFit::LineColor(kBlack), RooFit::MarkerSize(1.2));
  }
  legInfo["DATA"][Form("plot_Tot%s", dsName.c_str())] = "Data";

  if (ws.pdf(pdfName.c_str())) {
    RooArgList pdfList = ((RooAddPdf*)ws.pdf(pdfName.c_str()))->pdfList();
    if (pdfList.getSize()==1) {
      double norm = ws.data(dsName.c_str())->sumEntries();
      ws.pdf(pdfName.c_str())->plotOn(frame["MAIN"], RooFit::Name(Form("plot_%s", pdfName.c_str())),
                                      RooFit::Normalization(norm, RooAbsReal::NumEvent), RooFit::NumCPU(32),
                                      RooFit::LineColor(kBlack), RooFit::LineStyle(1), RooFit::Precision(1e-4)
                                      );
      legInfo["PDF"][Form("plot_%s", pdfName.c_str())] = "Total Fit";
      frame["PULL"] = (RooPlot*)frame["MAIN"]->emptyClone("PULL");
      RooHist *hpull = frame["MAIN"]->pullHist(0, 0, true);
      hpull->SetName("hpull");
      frame["PULL"]->addPlotable(hpull, "EP");
      drawMode = 1;
    }
    else {
      double norm = ws.data(dsName.c_str())->sumEntries();
      ws.pdf(pdfName.c_str())->plotOn(frame["MAIN"], RooFit::Name(Form("plot_%s", pdfName.c_str())),
                                      RooFit::Normalization(norm, RooAbsReal::NumEvent), RooFit::NumCPU(32),
                                      RooFit::LineColor(kBlack), RooFit::LineStyle(1), RooFit::Precision(1e-4)
                                      );
      frame["RATIO"] = (RooPlot*)frame["MAIN"]->emptyClone("RATIO");
      RooHist *hratio = makeRatioHist(frame["MAIN"], 0, 0, true, true);
      hratio->SetName("hratio");
      frame["RATIO"]->addPlotable(hratio, "EP");
      std::map< std::string , int > colorMap = { {"W" , kYellow} , {"DYZ" , kGreen+2} , {"WToTau" , kRed+1} , {"QCD" , kAzure-9} };
      TIterator* parIt = pdfList.createIterator();
      RooArgList* list = (RooArgList*)pdfList.Clone();
      std::map< double , RooAbsPdf* , std::greater< double > > pdfMap;
      for (RooAbsPdf* it = (RooAbsPdf*)parIt->Next(); it!=NULL; it = (RooAbsPdf*)parIt->Next() ) {
        std::string p = it->GetName() , obj;
        if (p.find("WToTau")!=std::string::npos) { obj = "WToTau"; } else if (p.find("W")!=std::string::npos) { obj = "W";} 
        else if (p.find("DYZ")!=std::string::npos) { obj = "DYZ"; } else if (p.find("QCD")!=std::string::npos) { obj = "QCD"; }
        double events = 0.;
        if (ws.var(("N_"+obj+cha+chg+"_"+col).c_str())) { events = ws.var(("N_"+obj+cha+chg+"_"+col).c_str())->getValV(); }
        else { events = ws.function(("N_"+obj+cha+chg+"_"+col).c_str())->getValV(); }
        pdfMap[events] = it;
      }
      for (auto& elem : pdfMap) {
        RooAbsPdf* it = elem.second;
        std::string p = it->GetName() , obj;
        if (p.find("WToTau")!=std::string::npos) { obj = "WToTau"; } else if (p.find("W")!=std::string::npos) { obj = "W";} 
        else if (p.find("DYZ")!=std::string::npos) { obj = "DYZ"; } else if (p.find("QCD")!=std::string::npos) { obj = "QCD"; }
        ws.pdf(pdfName.c_str())->plotOn(frame["MAIN"], RooFit::Name(Form("plot_%s", p.c_str())), RooFit::Components(*list),
                                        RooFit::Normalization(norm, RooAbsReal::NumEvent), RooFit::NumCPU(32),
                                        RooFit::FillStyle(1001), RooFit::FillColor(colorMap[obj]), RooFit::VLines(), RooFit::DrawOption("LCF"), RooFit::LineColor(kBlack), RooFit::LineStyle(1)
                                        );
        legInfo["TEMP"][Form("plot_%s", p.c_str())] = formatCut(obj);
        list->remove(*it);
        norm -= elem.first;
      }
      ws.data(dsName.c_str())->plotOn(frame["MAIN"], RooFit::Name(Form("plot_Tot%s", dsName.c_str())), RooFit::DataError(RooAbsData::SumW2), 
                                      RooFit::XErrorSize(0), RooFit::MarkerColor(kBlack), RooFit::LineColor(kBlack), RooFit::MarkerSize(1.2));
      drawMode = 2;
    }
  }

  TCanvas *cFig  = new TCanvas( Form("cMETFig_Tot%s", tag.c_str()), "cMETFig", 800, 800 );

  cFig->cd();

  if (drawMode==0) {
    TGaxis::SetMaxDigits(3);
    // Main Frame
    frame["MAIN"]->SetTitle("");
    frame["MAIN"]->GetYaxis()->SetTitleOffset(1.5);
    frame["MAIN"]->GetYaxis()->SetTitleSize(0.036);
    frame["MAIN"]->GetYaxis()->SetLabelSize(0.033);
    frame["MAIN"]->GetXaxis()->SetTitleOffset(1.1);
    frame["MAIN"]->GetXaxis()->SetTitleSize(0.036);
    frame["MAIN"]->GetXaxis()->SetLabelSize(0.033);
    frame["MAIN"]->GetXaxis()->SetTitle("|#slash{E}_{T}| (GeV/c)");
    pad["MAIN"] = new TPad( Form("padMAIN_Tot%s", tag.c_str()), "", 0, 0, 1, 1 );
  }
  else if (drawMode>0) {
    TGaxis::SetMaxDigits(3);
    // Main Frame
    frame["MAIN"]->SetTitle("");
    frame["MAIN"]->GetYaxis()->SetTitleOffset(1.5);
    frame["MAIN"]->GetYaxis()->SetTitleSize(0.036*(1./0.8));
    frame["MAIN"]->GetYaxis()->SetLabelSize(0.033*(1./0.8));
    frame["MAIN"]->GetXaxis()->SetTitleOffset(1.1);
    frame["MAIN"]->GetXaxis()->SetTitleSize(0.036);
    frame["MAIN"]->GetXaxis()->SetLabelSize(0.033);
    frame["MAIN"]->GetXaxis()->SetTitleOffset(3);
    frame["MAIN"]->GetXaxis()->SetLabelOffset(3);
    frame["MAIN"]->GetXaxis()->SetTitle("");
    pad["MAIN"] = new TPad( Form("padMAIN_Tot%s", tag.c_str()), "", 0, 0.2, 1, 1     );
    pad["MAIN"]->SetFixedAspectRatio(kTRUE);
    pad["MAIN"]->SetBottomMargin(0.015);
  }
  if (drawMode==1) {
    // Pull Frame
    frame["PULL"]->SetTitle("");
    frame["PULL"]->GetYaxis()->CenterTitle(kTRUE);
    frame["PULL"]->GetYaxis()->SetTitleOffset(0.4);
    frame["PULL"]->GetYaxis()->SetTitleSize(0.15);
    frame["PULL"]->GetYaxis()->SetLabelSize(0.13);
    frame["PULL"]->GetYaxis()->SetNdivisions(204);
    frame["PULL"]->GetYaxis()->SetTitle("Pull");
    frame["PULL"]->GetXaxis()->SetTitleOffset(1);
    frame["PULL"]->GetXaxis()->SetTitleSize(0.15);
    frame["PULL"]->GetXaxis()->SetLabelSize(0.15);
    frame["PULL"]->GetXaxis()->SetTitle("|#slash{E}_{T}| (GeV/c)");
    frame["PULL"]->GetYaxis()->SetRangeUser(-4.0, 4.0);
    pad["PULL"] = new TPad( Form("padPULL_Tot%s", tag.c_str()), "", 0, 0,    1, 0.20 );
    pad["PULL"]->SetFixedAspectRatio(kTRUE);
    pad["PULL"]->SetTopMargin(0.02);
    pad["PULL"]->SetBottomMargin(0.4);
    pad["PULL"]->SetFillStyle(4000);
    pad["PULL"]->SetFrameFillStyle(4000);
    pad["PULL"]->SetGridx(kTRUE);
    pad["PULL"]->SetGridy(kTRUE);
    // Draw the Pull
    pad["PULL"]->Draw();
    pad["PULL"]->cd();
    frame["PULL"]->Draw();
    printChi2(pad["PULL"], ws, *frame["MAIN"], "MET", dsName, pdfName, frame["MAIN"]->GetNbinsX());
    TLine *pline = new TLine(frame["PULL"]->GetXaxis()->GetXmin(), 0.0, frame["PULL"]->GetXaxis()->GetXmax(), 0.0);
    pline->Draw("same");
    pad["PULL"]->Update();
  }
  else if (drawMode==2) {
    // Pull Frame
    frame["RATIO"]->SetTitle("");
    frame["RATIO"]->GetYaxis()->CenterTitle(kTRUE);
    frame["RATIO"]->GetYaxis()->SetTitleOffset(0.4);
    frame["RATIO"]->GetYaxis()->SetTitleSize(0.15);
    frame["RATIO"]->GetYaxis()->SetLabelSize(0.13);
    frame["RATIO"]->GetYaxis()->SetNdivisions(204);
    frame["RATIO"]->GetYaxis()->SetTitle("#frac{DATA}{MC}");
    frame["RATIO"]->GetXaxis()->SetTitleOffset(1);
    frame["RATIO"]->GetXaxis()->SetTitleSize(0.15);
    frame["RATIO"]->GetXaxis()->SetLabelSize(0.15);
    frame["RATIO"]->GetXaxis()->SetTitle("|#slash{E}_{T}| (GeV/c)");
    frame["RATIO"]->GetYaxis()->SetRangeUser(0.0, 2.0);
    pad["RATIO"] = new TPad( Form("padRATIO_Tot%s", tag.c_str()), "", 0, 0,    1, 0.20 );
    pad["RATIO"]->SetFixedAspectRatio(kTRUE);
    pad["RATIO"]->SetTopMargin(0.02);
    pad["RATIO"]->SetBottomMargin(0.4);
    pad["RATIO"]->SetFillStyle(4000);
    pad["RATIO"]->SetFrameFillStyle(4000);
    pad["RATIO"]->SetGridx(kTRUE);
    pad["RATIO"]->SetGridy(kTRUE);
    // Draw the Ratio
    pad["RATIO"]->Draw();
    pad["RATIO"]->cd();
    frame["RATIO"]->Draw();
    TLine *pline = new TLine(frame["RATIO"]->GetXaxis()->GetXmin(), 1.0, frame["RATIO"]->GetXaxis()->GetXmax(), 1.0);
    pline->Draw("same");
    pad["RATIO"]->Update();
  }

  setRange(frame["MAIN"], ws, "MET", dsName, setLogScale);

  cFig->cd();
  pad["MAIN"]->Draw();
  pad["MAIN"]->cd();
  frame["MAIN"]->Draw();

  int lumiId = 0;
  if (col=="pPb") { lumiId = 109; } else if (col=="Pbp") { lumiId = 110; } else if (col=="PA") { lumiId = 111; }
  CMS_lumi(pad["MAIN"], lumiId, 33, "");

  printElectroWeakMETParameters(pad["MAIN"], ws, pdfName, drawMode);
  std::vector< std::string > text = { process };
  if (ws.obj(("CutAndCount_"+tag).c_str())) {
    text.push_back( formatCut( ((TObjString*)ws.obj(("CutAndCount_"+tag).c_str()))->GetString().Data(), varEWQLabel ) );
  }
  printElectroWeakBinning(pad["MAIN"], ws, dsName, text, drawMode);
  printElectroWeakLegend(pad["MAIN"], *frame["MAIN"], legInfo, drawMode);
  ws.import(*frame["MAIN"], Form("frame_Tot%s", tag.c_str()));

  pad["MAIN"]->SetLogy(setLogScale);
  pad["MAIN"]->Update();

  // Save the plot in different formats
  gSystem->mkdir(Form("%splot/root/", outputDir.c_str()), kTRUE);
  cFig->SaveAs(Form("%splot/root/%s.root", outputDir.c_str(), fileName.c_str()));
  gSystem->mkdir(Form("%splot/png/", outputDir.c_str()), kTRUE);
  cFig->SaveAs(Form("%splot/png/%s.png", outputDir.c_str(), fileName.c_str()));
  gSystem->mkdir(Form("%splot/pdf/", outputDir.c_str()), kTRUE);
  cFig->SaveAs(Form("%splot/pdf/%s.pdf", outputDir.c_str(), fileName.c_str()));

  cFig->Clear();
  cFig->Close();

  return true;
};


bool getVar(std::vector<RooRealVar*>& varVec, const RooWorkspace& ws, const std::string& name, const std::string& pdfName)
{
  varVec.clear();
  TIterator* parIt  = ws.allVars().selectByAttrib("Constant", kFALSE)->createIterator();
  for (RooRealVar* it = (RooRealVar*)parIt->Next(); it!=NULL; it = (RooRealVar*)parIt->Next() ) {
    std::string s(it->GetName());
    if((s.find("Pl")!=std::string::npos)!=(pdfName.find("Pl")!=std::string::npos)){ continue; }
    if (s.find(name)!=std::string::npos) { varVec.push_back(it);}
  }
  if (varVec.size()>0) return true;
  parIt = ws.allFunctions().selectByAttrib("Constant", kFALSE)->createIterator();
  for (RooRealVar* it = (RooRealVar*)parIt->Next(); it!=NULL; it = (RooRealVar*)parIt->Next() ) {
    std::string s(it->GetName());
    if((s.find("Pl")!=std::string::npos)!=(pdfName.find("Pl")!=std::string::npos)){ continue; }
    if (s.find(name)!=std::string::npos) { varVec.push_back(it); }
  }
  return (varVec.size()>0);
};


void parseVarName(const std::string& name, std::string& label)
{
  label = "";
  // Parse the parameter's labels
  stringstream ss(name); std::string s1, s2, s3;
  getline(ss, s1, '_'); getline(ss, s2, '_'); getline(ss, s3, '_');
  // Format QCD MET model parameters
  if (s1=="Alpha"){ s1="#alpha"; } else if (s1=="Beta"){ s1="#beta"; } 
  else if (s1=="XSection"){ s1="#sigma"; } else if (s1=="AccXEff"){ s1="#alphax#epsilon"; }
  // Format Object name
  std::string chg = ""; if (s2.find("Pl")!=std::string::npos) { chg = "+"; } else if (s2.find("Mi")!=std::string::npos) { chg = "-"; }
  if (s2.find("WToTau")!=std::string::npos) { s2 = "W#rightarrow#tau"; } else if (s2.find("W")!=std::string::npos) { s2 = "W"; } 
  else if (s2.find("DYZ")!=std::string::npos) { s2 = "Z/#gamma"; } else if (s2.find("QCD")!=std::string::npos) { s2 = "QCD"; }
  s2 = ( s2 + chg );
  if(s3!=""){ label = Form("%s_{%s}^{%s}", s1.c_str(), s2.c_str(), s3.c_str()); } else { label = Form("%s^{%s}", s1.c_str(), s2.c_str()); }
  return;
};


void printElectroWeakMETParameters(TPad* pad, const RooWorkspace& ws, const std::string& pdfName, const uint& drawMode)
{
  pad->cd();
  float xPos = 0.7, yPos = 0.74, dYPos = 0.045, dy = 0.025;
  TLatex t = TLatex(); t.SetNDC(); t.SetTextSize(0.023);
  if (drawMode>0) { dy = 0.065; dYPos *= (1./0.8); t.SetTextSize(0.023*(1./0.8)); }
  std::vector<RooRealVar*> vars; std::string label;
  if (getVar(vars, ws, "XSection_", pdfName)) {
    for (auto& v : vars) { parseVarName(v->GetName(),label); if(label!="") { t.DrawLatex(xPos, yPos-dy, Form("%s = %.3f#pm%.3f", label.c_str(), v->getValV(), v->getError())); dy+=dYPos; } }
  }
  if (vars.size()==0) {
    if (getVar(vars, ws, "N_", pdfName)) {
      for (auto& v : vars) { parseVarName(v->GetName(),label); if(label!="") { t.DrawLatex(xPos, yPos-dy, Form("%s = %.0f#pm%.0f", label.c_str(), v->getValV(), v->getError())); dy+=dYPos; } }
    }
  }
  TIterator* parIt;
  if (ws.pdf(pdfName.c_str())) { parIt = ws.pdf(pdfName.c_str())->getParameters(RooArgSet(*ws.var("MET")))->selectByAttrib("Constant", kFALSE)->createIterator(); }
  else { parIt = ws.allVars().selectByAttrib("Constant", kFALSE)->createIterator(); }
  for (RooRealVar* it = (RooRealVar*)parIt->Next(); it!=NULL; it = (RooRealVar*)parIt->Next() ) {
    // Parse the parameter's labels
    std::string label="", s(it->GetName());
    // Ignore dataset variables
    if(s=="MET" || s=="Muon_Pt" || s=="Muon_Eta" || s=="Muon_Iso" || s=="Muon_MT" || s=="Event_Type" || s=="Centrality"){ continue; }
    if((s.find("Pl")!=std::string::npos)!=(pdfName.find("Pl")!=std::string::npos)){ continue; }
    if(s.find("AccXEff")!=std::string::npos || s.find("XSection")!=std::string::npos) { continue; }
    parseVarName(it->GetName(), label); if (label=="") continue;
    // Print the parameter's results
    if(s.find("N")!=std::string::npos){ t.DrawLatex(xPos, yPos-dy, Form("%s = %.0f#pm%.0f", label.c_str(), it->getValV(), it->getError())); dy+=dYPos; }
    else { t.DrawLatex(xPos, yPos-dy, Form("%s = %.3f#pm%.3f", label.c_str(), it->getValV(), it->getError())); dy+=dYPos; }
  }
  delete parIt;
  pad->Update();
  return;
};


void printElectroWeakBinning(TPad* pad, const RooWorkspace& ws, const std::string& dsName, const std::vector< std::string >& text, const uint& drawMode)
{
  pad->cd();
  float xPos = 0.2, yPos = 0.89, dYPos = 0.045, dy = 0.025;
  TLatex t = TLatex(); t.SetNDC(); t.SetTextSize(0.025);
  if (drawMode>0) { dy *= (1./0.8); dYPos *= (1./0.8); t.SetTextSize(0.023*(1./0.8)); }
  t.DrawLatex(xPos, yPos-dy, Form("%s", text[0].c_str())); dy+=dYPos;
  TIterator* parIt = ((RooDataSet*)ws.data(dsName.c_str()))->get()->createIterator();
  for (RooRealVar* it = (RooRealVar*)parIt->Next(); it!=NULL; it = (RooRealVar*)parIt->Next() ) {
    if (std::string(it->GetName())=="MET" || std::string(it->GetName())=="Muon_Pt") continue;
    std::string varName = it->GetName();
    double defaultMin = 0.0 , defaultMax = 100000.0;
    if (varName=="Muon_Eta") { defaultMin = -2.5; defaultMax = 2.5; }
    if (ws.var(varName.c_str())) {
      if (ws.var(varName.c_str())->getMin()!=defaultMin && ws.var(varName.c_str())->getMax()==defaultMax) {
        t.DrawLatex(xPos, yPos-dy, Form("%g #leq %s", ws.var(varName.c_str())->getMin(), varEWQLabel[varName].c_str())); dy+=dYPos;
      }
      if (ws.var(varName.c_str())->getMin()==defaultMin && ws.var(varName.c_str())->getMax()!=defaultMax) {
        t.DrawLatex(xPos, yPos-dy, Form("%s < %g", varEWQLabel[varName].c_str(), ws.var(varName.c_str())->getMax())); dy+=dYPos;
      }
      if (ws.var(varName.c_str())->getMin()!=defaultMin && ws.var(varName.c_str())->getMax()!=defaultMax) {
        t.DrawLatex(xPos, yPos-dy, Form("%g #leq %s < %g", ws.var(varName.c_str())->getMin(), varEWQLabel[varName].c_str(), ws.var(varName.c_str())->getMax())); dy+=dYPos;
      }
    }
  }
  for (auto& txt : text) { if (text[0]!=txt) { t.DrawLatex(xPos, yPos-dy, Form("%s", txt.c_str())); dy+=dYPos; } }
  pad->Update();
  return;
};


void printElectroWeakLegend(TPad* pad, const RooPlot& frame, const StrMapMap& legInfo, const uint& drawMode)
{
  pad->cd();
  double ymax = 0.89, xmax = 0.65, dy = (0.89-0.64), dx = (0.65-0.48);
  TLegend* leg = new TLegend(xmax-dx, ymax-dy, xmax, ymax); leg->SetTextSize(0.03);
  if (drawMode>0) { dy *= (1./0.8); leg->SetTextSize(0.03*(1./0.8)); }
  std::map< std::string , std::string > drawOption = { { "DATA" , "pe" } , { "PDF" , "l" } , { "TEMP" , "fl" } };
  for (auto& map : legInfo) {
    for (auto& elem : map.second) {
      if (frame.findObject(elem.first.c_str())) { leg->AddEntry(frame.findObject(elem.first.c_str()), elem.second.c_str(), drawOption[map.first].c_str()); }
    }
  }
  leg->Draw("same");
  pad->Update();
  return;
};


RooHist* makeRatioHist(RooPlot* frame, const char* histname, const char* curvename, bool normalize, bool useAverage)
{
  // Find curve object
  RooCurve* curve = (RooCurve*) frame->findObject(curvename,RooCurve::Class()) ;
  if (!curve) {
    std::cout << "RooPlot::residHist(" << std::string(curvename) << ") cannot find curve" << std::endl;
    return 0;
  }
  // Find histogram object
  RooHist* hist = (RooHist*) frame->findObject(histname,RooHist::Class()) ;
  if (!hist) {
    std::cout << "RooPlot::residHist(" << std::string(histname) << ") cannot find histogram" << std::endl;
    return 0;
  }
  // Copy all non-content properties from hist1
  RooHist* histTMP = new RooHist(2);
  // Determine range of curve
  Double_t xstart,xstop,y ;
  curve->GetPoint(0,xstart,y) ;
  curve->GetPoint(curve->GetN()-1,xstop,y) ;
  // Add histograms, calculate Poisson confidence interval on sum value
  for(Int_t i=0 ; i<hist->GetN() ; i++) {
    Double_t x,point;
    hist->GetPoint(i,x,point) ;
    // Only calculate pull for bins inside curve range
    if (x<xstart || x>xstop) continue ;
    Double_t yy ;
    if (useAverage) {
      Double_t exl = hist->GetErrorXlow(i);
      Double_t exh = hist->GetErrorXhigh(i) ;
      if (exl<=0 ) exl = hist->GetErrorX(i);
      if (exh<=0 ) exh = hist->GetErrorX(i);
      if (exl<=0 ) exl = 0.5*hist->getNominalBinWidth();
      if (exh<=0 ) exh = 0.5*hist->getNominalBinWidth();
      yy = point / curve->average(x-exl,x+exh) ;
    } else {
      yy = point / curve->interpolate(x) ;
    }
    Double_t dyl = hist->GetErrorYlow(i)*(yy/point);
    Double_t dyh = hist->GetErrorYhigh(i)*(yy/point);
    histTMP->addBinWithError(x,yy,dyl,dyh);
  }
  return histTMP;
};



#endif // #ifndef drawElectroWeakMETPlot_C

#ifndef drawElectroWeakMETPlot_C
#define drawElectroWeakMETPlot_C

#include "../Utilities/initClasses.h"
#include "TGaxis.h"

void       printElectroWeakMETParameters ( TPad& pad , const RooWorkspace& ws , const std::string& pdfName , const uint& drawMode );
void       printElectroWeakBinning       ( TPad& pad , const RooWorkspace& ws , const std::string& dsName , const std::vector< std::string >& text , const uint& drawMode );
void       printElectroWeakLegend        ( TPad& pad , TLegend& leg , const RooPlot& frame , const StrMapMap_t& legInfo );
void       setRange                      ( RooPlot& frame , const RooWorkspace& ws , const std::string& varName , const std::string& dsName , const bool& setLogScale );
bool       getVar                        ( std::vector<RooRealVar>& varVec, const RooWorkspace& ws, const std::string& name, const std::string& pdfName );
void       parseVarName                  ( const std::string& name, std::string& label );
bool       printChi2                     ( TPad& pad , RooWorkspace& ws , const RooPlot& frame , const string& varLabel , const string& dataLabel , const string& pdfLabel );

bool drawElectroWeakMETPlot( RooWorkspace& ws,  // Local Workspace
                             // Select the type of datasets to fit
                             const std::string& fileName,
                             const std::string& outputDir,
                             const int& nBins,
                             const bool& yLogScale
                             )
{
  //
  // set the CMS style
  setTDRStyle();
  //
  const std::string DSTAG = (ws.obj("DSTAG"))     ? ((TObjString*)ws.obj("DSTAG"))->GetString().Data()     : "";
  const std::string cha   = (ws.obj("channel"))   ? ((TObjString*)ws.obj("channel"))->GetString().Data()   : "";
  const std::string col   = (ws.obj("fitSystem")) ? ((TObjString*)ws.obj("fitSystem"))->GetString().Data() : "";
  const std::string chg   = (ws.obj("fitCharge")) ? ((TObjString*)ws.obj("fitCharge"))->GetString().Data() : "";
  const std::string obj   = (ws.obj("fitObject")) ? ((TObjString*)ws.obj("fitObject"))->GetString().Data() : "";
  //
  const std::string tag = ( obj + cha + chg + "_" + col );
  const std::string dsName = ( "d" + chg + "_" + DSTAG );
  const std::string pdfName = Form("pdfMET_Tot%s", tag.c_str());
  const bool paperStyle = false;
  const bool setLogScale = yLogScale;
  const std::vector< double > range = { ws.var("MET")->getMin(), ws.var("MET")->getMax() };
  //
  if (ws.data(dsName.c_str())==NULL) { std::cout << "[ERROR] Dataset " << dsName << " was not found!" << std::endl; return false; }
  //
  const bool isMC = (DSTAG.find("MC")!=std::string::npos);
  const bool isWeighted = ws.data(dsName.c_str())->isWeighted();
  int drawMode = 0;
  //
  // Format Object name
  std::string process = "";
  char chgL = ' '; if (chg=="Pl") { chgL = '+'; } else if (chg=="Mi") { chgL = '-'; }
  if      (obj=="WToTau") { process = Form("W^{%c}#rightarrow#tau^{%c}", chgL, chgL); } else if (obj=="W") { process = Form("W^{%c}", chgL);; } 
  else if (obj=="DY"    ) { process = "Z/#gamma*"; }
  else if (obj=="QCD"   ) { process = "QCD";      }
  else if (obj=="TTbar" ) { process = "t#bar{t}"; }
  if (cha=="ToMu") { process += Form("#rightarrow#mu^{%c}+x", chgL); }
  process = Form("#font[62]{#scale[1.1]{%s}}", process.c_str());
  //
  StrMapMap_t legInfo;
  //
  std::map< std::string , std::unique_ptr<RooPlot> > frame;
  std::map< std::string , TPad* > pad; // Unique Pointer does produce Segmentation Fault, so don't use it
  //
  // Create the main plot of the fit
  frame["MAIN"] = std::unique_ptr<RooPlot>(ws.var("MET")->frame( RooFit::Bins(nBins), RooFit::Range(range[0], range[1]) ));
  if (ws.data(("CutAndCount_"+dsName).c_str())) {
    ws.data(("CutAndCount_"+dsName).c_str())->plotOn(frame.at("MAIN").get(), RooFit::Name(Form("plot_Tot%s", dsName.c_str())), RooFit::DataError(RooAbsData::SumW2), 
                                                     RooFit::XErrorSize(0), RooFit::MarkerColor(kBlack), RooFit::LineColor(kBlack), RooFit::MarkerSize(1.2));
  }
  else {
    ws.data(dsName.c_str())->plotOn(frame.at("MAIN").get(), RooFit::Name(Form("plot_Tot%s", dsName.c_str())), RooFit::DataError(RooAbsData::SumW2), 
                                    RooFit::XErrorSize(0), RooFit::MarkerColor(kBlack), RooFit::LineColor(kBlack), RooFit::MarkerSize(1.2));
  }
  legInfo["DATA"][Form("plot_Tot%s", dsName.c_str())] = "Data";
  //
  if (ws.pdf(pdfName.c_str())) {
    RooArgList pdfList = ((RooAddPdf*)ws.pdf(pdfName.c_str()))->pdfList();
    if (pdfList.getSize()==1) {
      double norm = ws.data(dsName.c_str())->sumEntries();
      ws.pdf(pdfName.c_str())->plotOn(frame.at("MAIN").get(), RooFit::Name(Form("plot_%s", pdfName.c_str())),
                                      RooFit::Normalization(norm, RooAbsReal::NumEvent), RooFit::NumCPU(32),
                                      RooFit::LineColor(kBlack), RooFit::LineStyle(1), RooFit::Precision(1e-4)
                                      );
      legInfo["PDF"][Form("plot_%s", pdfName.c_str())] = "Total Fit";
      frame["PULL"] = std::unique_ptr<RooPlot>((RooPlot*)frame.at("MAIN")->emptyClone("PULL"));
      RooHist* hPull = new RooHist(2.0); // !!!DONT USE UNIQUE POINTER!!!!, 2 represents the binWidth of MET
      if (!makePullHist(*hPull, *frame.at("MAIN").get(), "", "", true)) { return false; }
      hPull->SetName("hPull");
      frame.at("PULL")->addPlotable(hPull, "EP");
      drawMode = 1;
    }
    else {
      double norm = ws.data(dsName.c_str())->sumEntries();
      const std::map< std::string , int > colorMap = { {"W" , kYellow} , {"DY" , kGreen+2} , {"WToTau" , kRed+1} , {"QCD" , kAzure-9} , {"TTbar" , kOrange+2} };
      std::unique_ptr<TIterator> parIt = std::unique_ptr<TIterator>(pdfList.createIterator());
      std::unique_ptr<RooArgList> list = std::unique_ptr<RooArgList>((RooArgList*)pdfList.Clone());      
      if (list==NULL) { std::cout << "[ERROR] List of PDFs from " << pdfName << " is empty!" << std::endl; return false; }
      std::map< double , RooAbsPdf* , std::greater< double > > pdfMap;
      for (RooAbsPdf* it = (RooAbsPdf*)parIt->Next(); it!=NULL; it = (RooAbsPdf*)parIt->Next() ) {
        std::string obj = it->GetName(); obj = obj.substr(obj.find("_")+1); obj = obj.substr(0, obj.find(cha));
        double events = 0.;
        if (ws.var(("N_"+obj+cha+chg+"_"+col).c_str())) { events = ws.var(("N_"+obj+cha+chg+"_"+col).c_str())->getValV(); }
        else if (ws.function(("N_"+obj+cha+chg+"_"+col).c_str())) { events = ws.function(("N_"+obj+cha+chg+"_"+col).c_str())->getValV(); }
        else { std::cout << "[ERROR] The variable " << ("N_"+obj+cha+chg+"_"+col) << " was not found in the workspace" << std::endl; return false; }
        pdfMap[events] = it;
      }
      for (const auto& elem : pdfMap) {
        RooAbsPdf* it = elem.second;
        const std::string name = it->GetName();
        std::string obj = name; obj = obj.substr(obj.find("_")+1); obj = obj.substr(0, obj.find(cha));
        ws.pdf(pdfName.c_str())->plotOn(frame.at("MAIN").get(), RooFit::Name(Form("plot_%s", name.c_str())), RooFit::Components(*list),
                                        RooFit::Normalization(norm, RooAbsReal::NumEvent), RooFit::NumCPU(32),
                                        RooFit::FillStyle(1001), RooFit::FillColor(colorMap.at(obj)), RooFit::VLines(), RooFit::DrawOption("LF"), RooFit::LineColor(kBlack), RooFit::LineStyle(1)
                                        );
        legInfo["TEMP"][Form("plot_%s", name.c_str())] = formatCut(obj);
        list->remove(*it);
        norm -= elem.first;
      }
      ws.data(dsName.c_str())->plotOn(frame.at("MAIN").get(), RooFit::Name(Form("plot_Tot%s", dsName.c_str())), RooFit::DataError(RooAbsData::SumW2), 
                                      RooFit::XErrorSize(0), RooFit::MarkerColor(kBlack), RooFit::LineColor(kBlack), RooFit::MarkerSize(1.2));
      norm = ws.data(dsName.c_str())->sumEntries();
      ws.pdf(pdfName.c_str())->plotOn(frame.at("MAIN").get(), RooFit::Name(Form("plot_%s", pdfName.c_str())),
                                      RooFit::Normalization(norm, RooAbsReal::NumEvent), RooFit::NumCPU(32),
                                      RooFit::LineColor(kBlack), RooFit::LineStyle(1), RooFit::Precision(1e-4)
                                      );
      //
      frame["PULL"] = std::unique_ptr<RooPlot>((RooPlot*)frame.at("MAIN")->emptyClone("PULL"));
      RooHist* hPull = new RooHist(2.0); // !!!DONT USE UNIQUE POINTER!!!!, 2 represents the binWidth of MET
      if (!makePullHist(*hPull, *frame.at("MAIN"), "", "", true)) { return false; }
      hPull->SetName("hPull");
      frame.at("PULL")->addPlotable(hPull, "EP");
      drawMode = 1;
    }
  }
  //
  std::unique_ptr<TCanvas> cFig  = std::unique_ptr<TCanvas>(new TCanvas( Form("cMETFig_Tot%s", tag.c_str()), "cMETFig", 800, 800 ));
  cFig->cd();
  //
  std::unique_ptr<TLine> pLine;
  if (drawMode==0) {
    TGaxis::SetMaxDigits(3);
    // Main Frame
    frame.at("MAIN")->SetTitle("");
    frame.at("MAIN")->GetYaxis()->SetTitleOffset(1.5);
    frame.at("MAIN")->GetYaxis()->SetTitleSize(0.036);
    frame.at("MAIN")->GetYaxis()->SetLabelSize(0.033);
    frame.at("MAIN")->GetXaxis()->SetTitleOffset(1.1);
    frame.at("MAIN")->GetXaxis()->SetTitleSize(0.036);
    frame.at("MAIN")->GetXaxis()->SetLabelSize(0.033);
    frame.at("MAIN")->GetXaxis()->SetTitle("|#slash{E}_{T}| (GeV/c)");
    pad["MAIN"] = new TPad( Form("padMAIN_Tot%s", tag.c_str()), "", 0, 0, 1, 1 );
  }
  else if (drawMode>0) {
    TGaxis::SetMaxDigits(3);
    // Main Frame
    frame.at("MAIN")->SetTitle("");
    frame.at("MAIN")->GetYaxis()->SetTitleOffset(1.5);
    frame.at("MAIN")->GetYaxis()->SetTitleSize(0.036*(1./0.8));
    frame.at("MAIN")->GetYaxis()->SetLabelSize(0.033*(1./0.8));
    frame.at("MAIN")->GetXaxis()->SetTitleOffset(1.1);
    frame.at("MAIN")->GetXaxis()->SetTitleSize(0.036);
    frame.at("MAIN")->GetXaxis()->SetLabelSize(0.033);
    frame.at("MAIN")->GetXaxis()->SetTitleOffset(3);
    frame.at("MAIN")->GetXaxis()->SetLabelOffset(3);
    frame.at("MAIN")->GetXaxis()->SetTitle("");
    pad["MAIN"] = new TPad( Form("padMAIN_Tot%s", tag.c_str()), "", 0, 0.2, 1, 1 );
    pad.at("MAIN")->SetFixedAspectRatio(kTRUE);
    pad.at("MAIN")->SetBottomMargin(0.015);
  }
  if (drawMode==1) {
    // Pull Frame
    frame.at("PULL")->SetTitle("");
    frame.at("PULL")->GetYaxis()->CenterTitle(kTRUE);
    frame.at("PULL")->GetYaxis()->SetTitleOffset(0.4);
    frame.at("PULL")->GetYaxis()->SetTitleSize(0.15);
    frame.at("PULL")->GetYaxis()->SetLabelSize(0.10);
    frame.at("PULL")->GetYaxis()->SetNdivisions(206);
    frame.at("PULL")->GetYaxis()->SetTitle("Pull");
    frame.at("PULL")->GetXaxis()->SetTitleOffset(1);
    frame.at("PULL")->GetXaxis()->SetTitleSize(0.15);
    frame.at("PULL")->GetXaxis()->SetLabelSize(0.15);
    frame.at("PULL")->GetXaxis()->SetTitle("|#slash{E}_{T}| (GeV/c)");
    frame.at("PULL")->GetYaxis()->SetRangeUser(-6.0, 6.0);
    pad["PULL"] = new TPad( Form("padPULL_Tot%s", tag.c_str()), "", 0, 0, 1, 0.20 );
    pad.at("PULL")->SetFixedAspectRatio(kTRUE);
    pad.at("PULL")->SetTopMargin(0.02);
    pad.at("PULL")->SetBottomMargin(0.4);
    pad.at("PULL")->SetFillStyle(4000);
    pad.at("PULL")->SetFrameFillStyle(4000);
    pad.at("PULL")->SetGridx(kTRUE);
    pad.at("PULL")->SetGridy(kTRUE);
    // Draw the Pull
    pad.at("PULL")->Draw();
    pad.at("PULL")->cd();
    frame.at("PULL")->Draw();
    printChi2(*pad.at("PULL"), ws, *frame.at("MAIN"), "MET", dsName, pdfName);
    pLine = std::unique_ptr<TLine>(new TLine(frame.at("PULL")->GetXaxis()->GetXmin(), 0.0, frame.at("PULL")->GetXaxis()->GetXmax(), 0.0));
    pLine->Draw("same");
    pad.at("PULL")->Update();
  }
  //
  setRange(*frame.at("MAIN"), ws, "MET", dsName, setLogScale);
  //
  cFig->cd();
  pad.at("MAIN")->Draw();
  pad.at("MAIN")->cd();
  frame.at("MAIN")->Draw();
  //
  int lumiId = 0;
  if (col=="pPb") { lumiId = 115; } else if (col=="Pbp") { lumiId = 116; } else if (col=="PA") { lumiId = 117; }
  CMS_lumi(pad.at("MAIN"), lumiId, 33, "");
  //
  printElectroWeakMETParameters(*pad.at("MAIN"), ws, pdfName, drawMode);
  std::vector< std::string > text = { process };
  if (ws.obj(("CutAndCount_"+tag).c_str())) {
    text.push_back( formatCut( ((TObjString*)ws.obj(("CutAndCount_"+tag).c_str()))->GetString().Data(), varEWQLabel ) );
  }
  //
  printElectroWeakBinning(*pad.at("MAIN"), ws, dsName, text, drawMode);
  //
  double ymax = 0.89, xmax = 0.65, dy = (0.89-0.64), dx = (0.65-0.48);
  TLegend leg(xmax-dx, ymax-dy, xmax, ymax); leg.SetTextSize(0.03);
  if (drawMode>0) { dy *= (1./0.8); leg.SetTextSize(0.03*(1./0.8)); }
  printElectroWeakLegend(*pad.at("MAIN"), leg, *frame.at("MAIN"), legInfo);
  //
  ws.import(*frame.at("MAIN"), Form("frame_Tot%s", tag.c_str()));
  //
  pad.at("MAIN")->SetLogy(setLogScale);
  pad.at("MAIN")->Update();
  //
  // Save the plot in different formats
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
  return true;
};


bool getVar(std::vector<RooRealVar>& varVec, const RooWorkspace& ws, const std::string& name, const std::string& pdfName)
{
  varVec.clear();
  std::unique_ptr<TIterator> parIt = std::unique_ptr<TIterator>(ws.allVars().selectByAttrib("Constant", kFALSE)->createIterator());
  for (RooRealVar* it = (RooRealVar*)parIt->Next(); it!=NULL; it = (RooRealVar*)parIt->Next() ) {
    std::string s(it->GetName());
    if((s.find("Pl")!=std::string::npos)!=(pdfName.find("Pl")!=std::string::npos)){ continue; }
    if (s.find(name)!=std::string::npos) { varVec.push_back(*it);}
  }
  if (varVec.size()>0) return true;
  std::unique_ptr<TIterator> fncIt = std::unique_ptr<TIterator>(ws.allFunctions().selectByAttrib("Constant", kFALSE)->createIterator());
  for (RooRealVar* it = (RooRealVar*)fncIt->Next(); it!=NULL; it = (RooRealVar*)fncIt->Next() ) {
    std::string s(it->GetName());
    if((s.find("Pl")!=std::string::npos)!=(pdfName.find("Pl")!=std::string::npos)){ continue; }
    if (s.find(name)!=std::string::npos) { varVec.push_back(*it); }
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
  else if (s1=="Sigma0"){ s1="#sigma0"; } else if (s1=="Sigma1"){ s1="#sigma1"; } else if (s1=="Sigma2"){ s1="#sigma2"; }
  else if (s1=="XSection"){ s1="#sigma"; } else if (s1=="AccXEff"){ s1="#alphax#epsilon"; }
  // Format Object name
  std::string chg = ""; if (s2.find("Pl")!=std::string::npos) { chg = "+"; } else if (s2.find("Mi")!=std::string::npos) { chg = "-"; }
  if (s2.find("WToTau")!=std::string::npos) { s2 = "W#rightarrow#tau"; } else if (s2.find("W")!=std::string::npos) { s2 = "W"; } 
  else if (s2.find("DYZ")!=std::string::npos) { s2 = "Z/#gamma*"; } else if (s2.find("QCD")!=std::string::npos) { s2 = "QCD"; }
  else if (s2.find("TTbar")!=std::string::npos) { s2 = "t#bar{t}"; }
  s2 = ( s2 + chg );
  if(s3!=""){ label = Form("%s_{%s}^{%s}", s1.c_str(), s2.c_str(), s3.c_str()); } else { label = Form("%s^{%s}", s1.c_str(), s2.c_str()); }
  return;
};


void printElectroWeakMETParameters(TPad& pad, const RooWorkspace& ws, const std::string& pdfName, const uint& drawMode)
{
  pad.cd();
  float xPos = 0.7, yPos = 0.74, dYPos = 0.045, dy = 0.025;
  TLatex t = TLatex(); t.SetNDC(); t.SetTextSize(0.023);
  if (drawMode>0) { dy = 0.065; dYPos *= (1./0.8); t.SetTextSize(0.023*(1./0.8)); }
  std::vector<RooRealVar> vars; std::string label;
  if (vars.size()==0) {
    if (getVar(vars, ws, "N_", pdfName)) {
      for (const auto& v : vars) { parseVarName(v.GetName(), label); if(label!="") { t.DrawLatex(xPos, yPos-dy, Form("%s = %.0f#pm%.0f", label.c_str(), v.getValV(), v.getError())); dy+=dYPos; } }
    }
  }
  std::unique_ptr<TIterator> parIt;
  if (ws.pdf(pdfName.c_str())) {
    auto parList = std::unique_ptr<RooArgSet>(ws.pdf(pdfName.c_str())->getParameters(RooArgSet(*ws.var("MET"))));
    parIt = std::unique_ptr<TIterator>(parList->selectByAttrib("Constant", kFALSE)->createIterator());
  }
  else { parIt = std::unique_ptr<TIterator>(ws.allVars().selectByAttrib("Constant", kFALSE)->createIterator()); }
  for (RooRealVar* it = (RooRealVar*)parIt->Next(); it!=NULL; it = (RooRealVar*)parIt->Next() ) {
    // Parse the parameter's labels
    std::string label="", s(it->GetName());
    // Ignore dataset variables
    if(s=="MET" || s=="Muon_Pt" || s=="Muon_Eta" || s=="Muon_Iso" || s=="Muon_MT" || s=="Event_Type" || s=="Centrality") continue; 
    if((s.find("Pl")!=std::string::npos)!=(pdfName.find("Pl")!=std::string::npos)) continue;
    if(s.find("N_")!=std::string::npos) continue;
    parseVarName(it->GetName(), label); if (label=="") continue;
    // Print the parameter's results
    std::string txtLbl;
    if (s.find("Sigma2")!=std::string::npos) { txtLbl = Form("%s = %.3f#pm%.3f", label.c_str(), it->getValV()*1000., it->getError()*1000.); }
    else { txtLbl = Form("%s = %.3f#pm%.3f", label.c_str(), it->getValV(), it->getError()); }
    if (isAtLimit(*it)) { txtLbl += " (!)"; }
    t.DrawLatex(xPos, yPos-dy, txtLbl.c_str()); dy+=dYPos;
  }
  pad.Update();
  return;
};


void printElectroWeakBinning(TPad& pad, const RooWorkspace& ws, const std::string& dsName, const std::vector< std::string >& text, const uint& drawMode)
{
  pad.cd();
  float xPos = 0.2, yPos = 0.89, dYPos = 0.045, dy = 0.025;
  TLatex t = TLatex(); t.SetNDC(); t.SetTextSize(0.025);
  if (drawMode>0) { dy *= (1./0.8); dYPos *= (1./0.8); t.SetTextSize(0.023*(1./0.8)); }
  t.DrawLatex(xPos, yPos-dy, Form("%s", text[0].c_str())); dy+=dYPos;
  auto parIt = std::unique_ptr<TIterator>(((RooDataSet*)ws.data(dsName.c_str()))->get()->createIterator());
  for (RooRealVar* it = (RooRealVar*)parIt->Next(); it!=NULL; it = (RooRealVar*)parIt->Next() ) {
    if (std::string(it->GetName())=="MET" || std::string(it->GetName())=="Muon_Pt") continue;
    const std::string varName = it->GetName();
    double defaultMin = 0.0 , defaultMax = 100000.0;
    if (varName=="Muon_Eta") { defaultMin = -2.5; defaultMax = 2.5; }
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
  for (const auto& txt : text) { if (text[0]!=txt) { t.DrawLatex(xPos, yPos-dy, Form("%s", txt.c_str())); dy+=dYPos; } }
  pad.Update();
  return;
};


void printElectroWeakLegend(TPad& pad, TLegend& leg, const RooPlot& frame, const StrMapMap_t& legInfo)
{
  pad.cd();
  std::map< std::string , std::string > drawOption = { { "DATA" , "pe" } , { "PDF" , "l" } , { "TEMP" , "fl" } };
  for (const auto& map : legInfo) {
    for (const auto& elem : map.second) {
      if (frame.findObject(elem.first.c_str())) { leg.AddEntry(frame.findObject(elem.first.c_str()), elem.second.c_str(), drawOption[map.first].c_str()); }
    }
  }
  leg.Draw("same");
  pad.Update();
  return;
};


void setRange(RooPlot& frame, const RooWorkspace& ws, const std::string& varName, const std::string& dsName, const bool& setLogScale)
{
  // Find maximum and minimum points of Plot to rescale Y axis
  auto h = std::unique_ptr<TH1>(ws.data(dsName.c_str())->createHistogram("hist", *ws.var(varName.c_str()), RooFit::Binning(frame.GetNbinsX(), frame.GetXaxis()->GetXmin(), frame.GetXaxis()->GetXmax())));
  Double_t YMax = h->GetBinContent(h->GetMaximumBin());
  Double_t YMin = 1e99;
  for (int i=1; i<=h->GetNbinsX(); i++) if (h->GetBinContent(i)>0) YMin = min(YMin, h->GetBinContent(i));
  Double_t Yup(0.),Ydown(0.);
  if(setLogScale)
  {
    Yup = YMax*pow((YMax/0.1), 0.5);
    Ydown = 0.01;
  }
  else
  {
    Yup = YMax+(YMax-0.0)*0.7;
    Ydown = 0.0;
  }
  frame.GetYaxis()->SetRangeUser(Ydown,Yup);
  return;
};


bool printChi2(TPad& pad, RooWorkspace& ws, const RooPlot& frame, const string& varLabel, const string& dataLabel, const string& pdfLabel)
{
  //
  if (ws.data(dataLabel.c_str())==NULL) { std::cout << "[ERROR] Dataset " << dataLabel << " was not found!" << std::endl; return false; }
  if (ws.pdf (pdfLabel.c_str() )==NULL) { std::cout << "[ERROR] PDF "     << pdfLabel  << " was not found!" << std::endl; return false; }
  //
  pad.cd();
  //
  TH1D hData, hFit;
  if(!rooPlotToTH1(hData, hFit, frame)) { return false; }
  //
  auto parList = std::unique_ptr<RooArgSet>(ws.pdf(pdfLabel.c_str())->getParameters(*ws.data(dataLabel.c_str())));
  uint nFitPar = parList->selectByAttrib("Constant", kFALSE)->getSize();
  RooHist hPull(hData.GetBinWidth(1));
  if (!makePullHist(hPull, frame, "", "", true)) { return false; }
  double* ypulls = hPull.GetY();
  uint nFullBins = 0; double chi2=0;
  for (int i = 0; i < hData.GetNbinsX(); i++) {
    if ( (hData.GetBinError(i+1) != 0.0) && ((hData.GetBinContent(i+1) != 0.0) || (hFit.GetBinContent(i+1) != 0.0))) {
      chi2 += ypulls[i]*ypulls[i];
      nFullBins++;
    }
  }
  const int ndof = (nFullBins - nFitPar);
  std::cout << "[INFO] Using Standard method gives Chi2/NDoF " << chi2 << " / " << ndof << std::endl;
  //double chi2T = 0.; int ndofT = 0 , igood = -1;
  //hData.Chi2TestX(&hFit, chi2T, ndofT, igood, "UW");
  //std::cout << "[INFO] Using TH1 method gives Chi2/NDoF " << chi2T << " / " << ndofT << std::endl;
  //
  TLatex t = TLatex(); t.SetNDC(); t.SetTextSize(0.12);
  t.DrawLatex(0.76, 0.85, Form("#chi^{2}/ndof = %.0f / %d ", chi2, ndof));
  RooRealVar chi2Var((string("chi2_")+varLabel).c_str(),(string("chi2_")+varLabel).c_str(),chi2);
  RooRealVar ndofVar((string("ndof_")+varLabel).c_str(),(string("ndof_")+varLabel).c_str(),ndof);
  ws.import(chi2Var, kTRUE); ws.import(ndofVar, kTRUE);
  return true;
};



#endif // #ifndef drawElectroWeakMETPlot_C

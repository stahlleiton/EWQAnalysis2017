#ifndef drawElectroWeakMETPlot_C
#define drawElectroWeakMETPlot_C

#include "../Utilities/initClasses.h"
#include "TGaxis.h"
#include "TLegendEntry.h"

void       printElectroWeakMETParameters ( TPad& pad , const RooWorkspace& ws , const std::string& pdfName , const uint& drawMode );
void       printElectroWeakBinning       ( TPad& pad , const RooWorkspace& ws , const std::string& dsName , const std::vector< std::string >& text , const uint& drawMode , const int plotStyle );
void       printElectroWeakLegend        ( TPad& pad , TLegend& leg , const RooPlot& frame , const StringDiMap_t& legInfo , const double size=0.05 );
void       setRange                      ( RooPlot& frame , const RooWorkspace& ws , const std::string& varName , const std::string& dsName , const bool& setLogScale , const int& nBins );
bool       getVar                        ( std::vector<RooRealVar>& varVec, const RooWorkspace& ws, const std::string& name, const std::string& pdfName );
void       parseVarName                  ( const std::string& name, std::string& label );
bool       printGoF                      ( TPad& pad , RooWorkspace& ws , const RooPlot& frame , const string& varLabel , const string& dataLabel , const string& pdfLabel );
void       formatLegendEntry             ( TLegendEntry& e , double size=0.060 ) { e.SetTextSize(size); };

bool drawElectroWeakMETPlot( RooWorkspace& ws,  // Local Workspace
                             // Select the type of datasets to fit
                             const std::string& fileName,
                             const std::string& outputDir,
                             const bool& yLogScale,
                             const double maxRng = -1.0,
                             const bool doGoF = true,
                             const int  plotStyle = 0, // 3: Thesis , 2: Paper , 1: PAS , 0: AN
                                   bool redoFrame = false,
                             const bool doQCDHist = false
                             )
{
  //
  // set the CMS style
  setTDRStyle();
  //
  const std::string DSTAG = (ws.obj("DSTAG")    ) ? ((RooStringVar*)ws.obj("DSTAG")    )->getVal() : "";
  const std::string cha   = (ws.obj("channel")  ) ? ((RooStringVar*)ws.obj("channel")  )->getVal() : "";
  const std::string col   = (ws.obj("fitSystem")) ? ((RooStringVar*)ws.obj("fitSystem"))->getVal() : "";
  const std::string chg   = (ws.obj("fitCharge")) ? ((RooStringVar*)ws.obj("fitCharge"))->getVal() : "";
  const std::string obj   = (ws.obj("fitObject")) ? ((RooStringVar*)ws.obj("fitObject"))->getVal() : "";
  //
  const std::string tag = ( obj + cha + chg + "_" + col );
  const std::string dsName = ( "d" + chg + "_" + DSTAG );
  const std::string dsNameFit = ( (ws.data((dsName+"_FIT").c_str())!=NULL) ? (dsName+"_FIT") : dsName );
  const std::string pdfName = Form("pdfMET_Tot%s", tag.c_str());
  const bool setLogScale = yLogScale;
  //
  // Create the Range for Plotting
  const double binWidth = ws.var("MET")->getBinWidth(0);
  const double minRange = ws.var("MET")->getMin();
  const double maxRange = ( (maxRng>0.0) ? maxRng : ws.var("MET")->getMax() );
  const int    nBins    = int(std::round((maxRange - minRange)/binWidth));
  ws.var("MET")->setRange("METWindowPlot", minRange, maxRange);
  ws.var("MET")->setBins(nBins, "METWindowPlot");
  // BUG FIX
  const int oNBins = ws.var("MET")->getBins(); const double oMinRange = ws.var("MET")->getMin(); const double oMaxRange = ws.var("MET")->getMax();
  ws.var("MET")->setBins(nBins); ws.var("MET")->setRange(minRange, maxRange);
  //
  const bool useDS = (ws.data(dsName.c_str())!=NULL);
  const bool isMC = (DSTAG.find("MC")!=std::string::npos);
  const bool isWeighted = (useDS ? ws.data(dsName.c_str())->isWeighted() : false);
  int drawMode = 0;
  bool drawPull = false;  // false : Draw DATA/FIT , true : Draw the Pull
  //
  // Format Object name
  std::string process = "";
  std::string chgL = " "; if (chg=="Pl") { chgL = "+"; } else if (chg=="Mi") { chgL = "#font[122]{\55}"; }
  if      (obj=="WToTau" ) { process = Form("W^{%s}#rightarrow#tau^{%s}", chgL.c_str(), chgL.c_str()); }
  else if (obj=="DYToTau") { process = Form("Z/#gamma*#rightarrow#tau^{%s}", chgL.c_str()); }
  else if (obj=="W"      ) { process = Form("W^{%s}" , chgL.c_str()); }
  else if (obj=="WZ"     ) { process = Form("W^{%s}Z", chgL.c_str()); }
  else if (obj=="WW"     ) { process = "W^{+}W^{-}"; }
  else if (obj=="DY"     ) { process = "Z/#gamma*";  }
  else if (obj=="Z"      ) { process = "Z";          }
  else if (obj=="ZZ"     ) { process = "ZZ";         }
  else if (obj=="QCD"    ) { process = "QCD";        }
  else if (obj=="TTbar"  ) { process = "t#bar{t}";   }
  if (cha=="ToMu") {
    if (obj=="W") {
      if (plotStyle < 3) { process += Form("#kern[0.2]{#rightarrow}#kern[0.2]{#mu^{%s}}#kern[0.2]{%s}", chgL.c_str(), (chg=="Pl"?"#nu_{#mu}":"#bar{#nu}_{#mu}")); }
      else { process += Form(" #rightarrow #mu^{%s} + %s", chgL.c_str(), (chg=="Pl"?"#nu_{#mu}":"#bar{#nu}_{#mu}")); }
    }
    else { process += Form(" #rightarrow #mu^{%s} + x", chgL.c_str()); }
  }
  if (obj=="DY" ) { process = "Z/#gamma* #rightarrow #mu^{+} + #mu^{-}"; }
  if (obj=="Z"  ) { process = "Z #rightarrow #mu^{+} + #mu^{-}"; }
  process = Form("#font[62]{#scale[1.1]{%s}}", process.c_str());
  //
  StringDiMap_t legInfo;
  //
  std::map< std::string , std::unique_ptr<RooPlot> > frame;
  std::map< std::string , TPad* > pad; // Unique Pointer does produce Segmentation Fault, so don't use it
  //
  // Create the main plot of the fit
  const std::string frameName = Form("frame_Tot%s", tag.c_str());
  if (ws.obj(frameName.c_str())==NULL) { redoFrame = false; }
  //
  if (!redoFrame && ws.obj(frameName.c_str())!=NULL) { frame["MAIN"] = std::unique_ptr<RooPlot>((RooPlot*)ws.obj(frameName.c_str())); }
  else {
    if (useDS==false) { std::cout << "[ERROR] Dataset " << dsName << " was not found!" << std::endl; return false; }
    frame["MAIN"] = std::unique_ptr<RooPlot>(ws.var("MET")->frame( RooFit::Range("METWindowPlot") ));
    if (ws.data(("CutAndCount_"+dsName).c_str())) {
      ws.data(("CutAndCount_"+dsName).c_str())->plotOn(frame.at("MAIN").get(), RooFit::Name(Form("plot_Tot%s", dsName.c_str())), RooFit::Binning("METWindowPlot"),
                                                       RooFit::DataError(RooAbsData::SumW2), RooFit::XErrorSize(0),
                                                       RooFit::MarkerColor(kBlack), RooFit::LineColor(kBlack), RooFit::MarkerSize(1.55));
    }
    else {
      ws.data(dsName.c_str())->plotOn(frame.at("MAIN").get(), RooFit::Name(Form("plot_Tot%s", dsName.c_str())), RooFit::Binning("METWindowPlot"),
                                      RooFit::MarkerColor(kBlack), RooFit::LineColor(kBlack), RooFit::MarkerSize(1.55));
    }
    //
    if (ws.pdf(pdfName.c_str())) {
      RooArgList pdfList = ((RooAddPdf*)ws.pdf(pdfName.c_str()))->pdfList();
      if (pdfList.getSize()==1) {
        const double norm = ((RooAddPdf*)ws.pdf(pdfName.c_str()))->expectedEvents(RooArgSet(*ws.var("MET")));
        ws.pdf(pdfName.c_str())->plotOn(frame.at("MAIN").get(), RooFit::Name(Form("plot_%s", pdfName.c_str())), RooFit::Range("METWindowPlot"), RooFit::NormRange("METWindowPlot"),
                                        RooFit::Normalization(norm, RooAbsReal::NumEvent), RooFit::Precision(1e-7),
                                        RooFit::LineColor(kBlack), RooFit::LineStyle(1)
                                        );
      }
      else {
        double norm = ((RooAddPdf*)ws.pdf(pdfName.c_str()))->expectedEvents(RooArgSet(*ws.var("MET")));
        const std::map< std::string , int > colorMap = { {"W" , kYellow} , {"DY" , kGreen+2} , {"WToTau" , kRed+1} , {"QCD" , kAzure-9} , {"TTbar" , kOrange+1} ,
                                                         {"DYToTau" , kBlue+2} , {"WW" , kMagenta+1} };
        const std::map< std::string , double > pdfMapOrder = { {"W" , 9} , {"QCD" , 8} , {"DY" , 7} , {"WToTau" , 6} , {"DYToTau" , 5} , {"TTbar" , 4} , {"WW" , 3} };
        std::unique_ptr<TIterator> pdfIt = std::unique_ptr<TIterator>(pdfList.createIterator());
        std::unique_ptr<RooArgList> list = std::unique_ptr<RooArgList>((RooArgList*)pdfList.Clone());
        if (list==NULL) { std::cout << "[ERROR] List of PDFs from " << pdfName << " is empty!" << std::endl; return false; }
        std::map< double , RooAbsPdf* , std::greater< double > > pdfMap;
        std::map< std::string , double > pdfEvt;
        for (RooAbsPdf* it = (RooAbsPdf*)pdfIt->Next(); it!=NULL; it = (RooAbsPdf*)pdfIt->Next() ) {
          std::string obj = it->GetName(); obj = obj.substr(obj.find("_")+1); obj = obj.substr(0, obj.find(cha));
          double events = 0.;
          if (ws.var(("N_"+obj+cha+chg+"_"+col).c_str())) { events = ws.var(("N_"+obj+cha+chg+"_"+col).c_str())->getValV(); }
          else if (ws.function(("N_"+obj+cha+chg+"_"+col).c_str())) { events = ws.function(("N_"+obj+cha+chg+"_"+col).c_str())->getValV(); }
          else { std::cout << "[ERROR] The variable " << ("N_"+obj+cha+chg+"_"+col) << " was not found in the workspace" << std::endl; return false; }
          pdfMap[pdfMapOrder.at(obj)] = it;
          pdfEvt[obj] = events;
        }
        for (const auto& elem : pdfMap) {
          RooAbsPdf* ito = elem.second;
          const std::string nameo = ito->GetName();
          std::string objo = nameo; objo = objo.substr(objo.find("_")+1); objo = objo.substr(0, objo.find(cha));
          if (norm > 0.0) {
            RooArgList coef; RooArgList pdfs;
            std::unique_ptr<TIterator> tmpIt = std::unique_ptr<TIterator>(list->createIterator());
            for (RooAbsPdf* it = (RooAbsPdf*)tmpIt->Next(); it!=NULL; it = (RooAbsPdf*)tmpIt->Next() ) {
              std::string obj = it->GetName(); obj = obj.substr(obj.find("_")+1); obj = obj.substr(0, obj.find(cha));
              const std::string name = ("N_"+obj+cha+chg+"_"+col);
              const std::string pdfN = ("pdfMET_"+obj+cha+chg+"_"+col);
              if (obj=="QCD" && doQCDHist) {
                const auto pdfName  = ("pdfBIN_"+obj+cha+chg+"_"+col);
                const auto dataName = ("dsBIN_"+obj+cha+chg+"_"+col);
                if (ws.data(dataName.c_str())==NULL) {
                  std::unique_ptr<RooDataHist> data = std::unique_ptr<RooDataHist>(it->generateBinned(RooArgSet(*ws.var("MET")), 1000000, RooFit::Name(dataName.c_str()), RooFit::ExpectedData()));
                  ws.import(*data);
                }
                if (ws.pdf(pdfName.c_str())==NULL) {
                  std::unique_ptr<RooHistPdf> pdf = std::unique_ptr<RooHistPdf>(new RooHistPdf(pdfName.c_str(), pdfName.c_str(), *ws.var("MET"), *((RooDataHist*)ws.data(dataName.c_str()))));
                  ws.import(*pdf);
                }
                pdfs.add(*ws.pdf(pdfName.c_str()));
              }
              else { pdfs.add(*ws.pdf(pdfN.c_str())); }
              if (ws.var(name.c_str())) { coef.add(*ws.var(name.c_str())); }
              else if (ws.function(name.c_str())) { coef.add(*ws.function(name.c_str())); }
            }
            const std::string pdfPlotName = Form("pdfPlot%s_%s", (redoFrame?"RE":""), nameo.c_str());
            auto pf = RooAddPdf(pdfPlotName.c_str(), pdfPlotName.c_str(), pdfs, coef); if (ws.pdf(pdfPlotName.c_str())==NULL) { ws.import(pf); }
            ws.pdf(pdfPlotName.c_str())->plotOn(frame.at("MAIN").get(), RooFit::Name(Form("plot_%s", nameo.c_str())), RooFit::Range("METWindowPlot"), RooFit::NormRange("METWindowPlot"),
                                                RooFit::Normalization(norm, RooAbsReal::NumEvent), RooFit::Precision(1e-6),
                                                RooFit::FillStyle(1001), RooFit::FillColor(colorMap.at(objo)), RooFit::VLines(), RooFit::DrawOption("F")
                                                );
          }
          list->remove(*ito);
          norm -= pdfEvt.at(objo);
        }
        //
        ws.data(dsName.c_str())->plotOn(frame.at("MAIN").get(), RooFit::Name(Form("plot_Tot%s", dsName.c_str())), RooFit::Binning("METWindowPlot"),
                                        RooFit::MarkerColor(kBlack), RooFit::LineColor(kBlack), RooFit::MarkerSize(1.55));
        //
        norm = ws.data(dsNameFit.c_str())->sumEntries();
        //
        std::string pdfNameTot = Form("pdfPlot%s_pdfMETTot_%s%s%s_%s", (redoFrame?"RE":""), obj.c_str(), cha.c_str(), chg.c_str(), col.c_str());
        if (ws.pdf(pdfNameTot.c_str())==NULL) { pdfNameTot = pdfName; }
        ws.pdf(pdfNameTot.c_str())->plotOn(frame.at("MAIN").get(), RooFit::Name(Form("plot_%s", pdfName.c_str())), RooFit::Range("METWindowPlot"), RooFit::NormRange("METWindowPlot"),
                                           RooFit::Normalization(norm, RooAbsReal::NumEvent), RooFit::Precision(1e-7),
                                           RooFit::LineColor(kBlack), RooFit::LineStyle(1)
                                           );
      }
    }
    // Store the frame
    frame.at("MAIN")->SetTitle(frameName.c_str());
    if (ws.obj(frameName.c_str())==NULL) { ws.import(*frame.at("MAIN"), frame.at("MAIN")->GetTitle()); }
  }
  //
  legInfo["DATA"][Form("plot_Tot%s", dsName.c_str())] = ( isMC ? "Simulation" : "Data" );
  //
  if (ws.pdf(pdfName.c_str())) {
    RooArgList pdfList = ((RooAddPdf*)ws.pdf(pdfName.c_str()))->pdfList();
    if (pdfList.getSize()==1) {
      legInfo["PDF"][Form("plot_%s", pdfName.c_str())] = "Fit";
      frame["EXTRA"] = std::unique_ptr<RooPlot>((RooPlot*)frame.at("MAIN")->emptyClone("EXTRA"));
      RooHist* hPull = new RooHist(2.0); // !!!DONT USE UNIQUE POINTER!!!!, 2 represents the binWidth of MET
      if (!makePullHist (*hPull, *frame.at("MAIN").get(), "", "", true)) { return false; }; drawPull = true;
      hPull->SetName("hPull");
      frame.at("EXTRA")->addPlotable(hPull, "EP");
      drawMode = 1;
    }
    else {
      const std::vector< std::string > pdfOrder = { "W" , "QCD" , "DY" , "WToTau" , "DYToTau" , "TTbar" };
      for (const auto& pdfT : pdfOrder) {
        const std::string obj = pdfT;
        const std::string name = Form("pdfMETTot_%s%s%s_%s", obj.c_str(), cha.c_str(), chg.c_str(), col.c_str());
        legInfo["TEMP"][Form("plot_%s", name.c_str())] = formatCut(obj);
      }
      frame["EXTRA"] = std::unique_ptr<RooPlot>((RooPlot*)frame.at("MAIN")->emptyClone("EXTRA"));
      RooHist* hExtra = new RooHist(2.0); // !!!DONT USE UNIQUE POINTER!!!!, 2 represents the binWidth of MET
      if (drawPull) { if (!makePullHist (*hExtra, *frame.at("MAIN"), "", "", true)) { return false; } }
      else          { if (!makeRatioHist(*hExtra, *frame.at("MAIN"), "", "", true)) { return false; } }
      hExtra->SetName("hExtra");
      frame.at("EXTRA")->addPlotable(hExtra, "EP");
      drawMode = 1;
    }
  }
  //
  std::unique_ptr<TCanvas> cFig  = std::unique_ptr<TCanvas>(new TCanvas( Form("cMETFig_Tot%s", tag.c_str()), "cMETFig", 800, 800 ));
  cFig->cd();
  //
  std::unique_ptr<TLine> pLine;
  if (drawMode==0) {
    //TGaxis::SetMaxDigits(4);
    // Main Frame
    frame.at("MAIN")->SetTitle("");
    frame.at("MAIN")->SetMarkerSize(2.05);
    frame.at("MAIN")->GetYaxis()->CenterTitle(kTRUE);
    frame.at("MAIN")->GetXaxis()->CenterTitle(kTRUE);
    frame.at("MAIN")->GetYaxis()->SetTitleOffset(1.5);
    frame.at("MAIN")->GetYaxis()->SetTitleSize(0.036);
    frame.at("MAIN")->GetYaxis()->SetLabelSize(0.033);
    frame.at("MAIN")->GetXaxis()->SetTitleOffset(1.1);
    frame.at("MAIN")->GetXaxis()->SetTitleSize(0.036);
    frame.at("MAIN")->GetXaxis()->SetLabelSize(0.033);
    frame.at("MAIN")->GetXaxis()->SetTitle("|#slash{E}_{T}| (GeV/c)");
    if (plotStyle==1 || plotStyle==2 || plotStyle==3) { frame.at("MAIN")->GetXaxis()->SetTitle("p^{miss}_{T} (GeV/c)"); }
    //if (plotStyle==3) { frame.at("MAIN")->GetXaxis()->SetTitle("#slash{E}_{T} (GeV/c)"); }
    pad["MAIN"] = new TPad( Form("padMAIN_Tot%s", tag.c_str()), "", 0, 0, 1, 1 );
  }
  else if (drawMode>0) {
    //TGaxis::SetMaxDigits(4);
    // Main Frame
    frame.at("MAIN")->SetTitle("");
    frame.at("MAIN")->SetMarkerSize(2.05);
    frame.at("MAIN")->GetYaxis()->CenterTitle(kTRUE);
    frame.at("MAIN")->GetXaxis()->CenterTitle(kTRUE);
    frame.at("MAIN")->GetYaxis()->SetTitleOffset(0.9);
    frame.at("MAIN")->GetYaxis()->SetTitleSize(0.060*(1./0.8));
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
    if (drawPull==false && (plotStyle==1 || plotStyle==2 || plotStyle==3)) { pad.at("MAIN")->SetBottomMargin(0.016); }
    pad.at("MAIN")->SetFillStyle(0);
    pad.at("MAIN")->SetFrameFillStyle(0);
  }
  if (drawMode==1) {
    // Pull Frame
    frame.at("EXTRA")->SetTitle("");
    frame.at("EXTRA")->SetMarkerSize(2.25);
    frame.at("EXTRA")->GetYaxis()->CenterTitle(kTRUE);
    frame.at("EXTRA")->GetXaxis()->CenterTitle(kTRUE);
    frame.at("EXTRA")->GetYaxis()->SetTitleOffset(0.25);
    frame.at("EXTRA")->GetYaxis()->SetTitleSize(0.25);
    frame.at("EXTRA")->GetYaxis()->SetLabelSize(0.15);
    frame.at("EXTRA")->GetYaxis()->SetNdivisions(204);
    if (drawPull) { frame.at("EXTRA")->GetYaxis()->SetTitle("Pull"); }
    else          { frame.at("EXTRA")->GetYaxis()->SetTitle("#frac{Data}{Fit}"); }
    frame.at("EXTRA")->GetXaxis()->SetTitleOffset(0.7);
    frame.at("EXTRA")->GetXaxis()->SetTitleSize(0.25);
    frame.at("EXTRA")->GetXaxis()->SetLabelSize(0.15);
    frame.at("EXTRA")->GetXaxis()->SetTitle("|#slash{E}_{T}| (GeV/c)");
    if (plotStyle==1 || plotStyle==2 || plotStyle==3) { frame.at("EXTRA")->GetXaxis()->SetTitle("p^{miss}_{T} (GeV/c)"); }
    //if (plotStyle==3) { frame.at("EXTRA")->GetXaxis()->SetTitle("#slash{E}_{T} (GeV/c)"); }
    if (drawPull) { frame.at("EXTRA")->GetYaxis()->SetRangeUser(-6.0, 6.0); }
    else if (plotStyle==1 || plotStyle==2 || plotStyle==3) { frame.at("EXTRA")->GetYaxis()->SetRangeUser(0., 2.6); }
    else { frame.at("EXTRA")->GetYaxis()->SetRangeUser(-0.01, 2.1); }
    pad["EXTRA"] = new TPad( Form("padEXTRA_Tot%s", tag.c_str()), "", 0, 0, 1, 0.212 );
    pad.at("EXTRA")->SetFixedAspectRatio(kTRUE);
    pad.at("EXTRA")->SetTopMargin(0.02);
    if (drawPull==false && (plotStyle==1 || plotStyle==2 || plotStyle==3)) { pad.at("EXTRA")->SetTopMargin(0.0); }
    pad.at("EXTRA")->SetBottomMargin(0.5);
    //cFig->SetBottomMargin(0.5);
    pad.at("EXTRA")->SetFillStyle(4000);
    pad.at("EXTRA")->SetFrameFillStyle(4000);
    //pad.at("EXTRA")->SetGridx(kTRUE);
    //pad.at("EXTRA")->SetGridy(kTRUE);
    // Draw the Pull
    pad.at("EXTRA")->Draw();
    pad.at("EXTRA")->cd();
    frame.at("EXTRA")->Draw();
    if (doGoF) { printGoF(*pad.at("EXTRA"), ws, *frame.at("MAIN"), "MET", dsName, pdfName); }
    if (drawPull) { pLine = std::unique_ptr<TLine>(new TLine(frame.at("EXTRA")->GetXaxis()->GetXmin(), 0.0, frame.at("EXTRA")->GetXaxis()->GetXmax(), 0.0)); }
    else          { pLine = std::unique_ptr<TLine>(new TLine(frame.at("EXTRA")->GetXaxis()->GetXmin(), 1.0, frame.at("EXTRA")->GetXaxis()->GetXmax(), 1.0)); }
    pLine->Draw("same");
    pad.at("EXTRA")->Update();
  }
  //
  setRange(*frame.at("MAIN"), ws, "MET", dsName, setLogScale, ws.var("MET")->getBins("METWindowPlot"));
  //
  cFig->cd();
  pad.at("MAIN")->Draw();
  pad.at("MAIN")->cd();
  frame.at("MAIN")->Draw();
  //
  int lumiId = 0;
  if (isMC) { if (col=="pPb") { lumiId = 112; } else if (col=="Pbp") { lumiId = 113; } else if (col=="PA") { lumiId = 114; } }
  else { if (col=="pPb") { lumiId = 11820; } else if (col=="Pbp") { lumiId = 11821; } else if (col=="PA") { lumiId = 11822; } }
  CMS_lumi(pad.at("MAIN"), lumiId, 33, "", false, 0.8, false);
  //
  if (plotStyle==0) { printElectroWeakMETParameters(*pad.at("MAIN"), ws, pdfName, drawMode); }
  std::vector< std::string > text = { process };
  if (ws.obj(("CutAndCount_"+tag).c_str())) {
    text.push_back( formatCut( ((RooStringVar*)ws.obj(("CutAndCount_"+tag).c_str()))->getVal(), varEWQLabel ) );
  }
  //
  printElectroWeakBinning(*pad.at("MAIN"), ws, dsName, text, drawMode, plotStyle);
  //
  double xmin = 0.49 , xmax = 0.66 , ymin = 0.58 , ymax = 0.89;
  if (maxRng>0. && maxRng<120.) { ymax = 0.72; xmax = 0.40; ymin = 0.49; }
  double legSize = 0.047 , dy = (ymax-ymin);
  if (plotStyle==1 || plotStyle==2 || plotStyle==3) {
    if (legInfo.size()>2) { xmin = 0.74; ymin = 0.25; xmax = 0.90; ymax = 0.71; legSize = 0.06; }
    if (legInfo.size()<3) { xmin = 0.74; ymin = 0.35; xmax = 0.90; ymax = 0.71; legSize = 0.06; }
  }
  auto leg = std::unique_ptr<TLegend>(new TLegend(xmin, ymin, xmax, ymax));
  if (drawMode>0) { dy *= (1./0.8); }
  printElectroWeakLegend(*pad.at("MAIN"), *leg, *frame.at("MAIN"), legInfo, legSize);
  //
  pad.at("MAIN")->SetLogy(setLogScale);
  pad.at("MAIN")->Update();
  //
  // Save the plot in different formats
  gSystem->mkdir(Form("%splot/root/", outputDir.c_str()), kTRUE);
  cFig->SaveAs(Form("%splot/root/%s.root", outputDir.c_str(), fileName.c_str()));
  gSystem->mkdir(Form("%splot/pdf/", outputDir.c_str()), kTRUE);
  cFig->SaveAs(Form("%splot/pdf/%s.pdf", outputDir.c_str(), fileName.c_str()));
  gSystem->mkdir(Form("%splot/png/", outputDir.c_str()), kTRUE);
  cFig->SaveAs(Form("%splot/png/%s.png", outputDir.c_str(), fileName.c_str()));
  gSystem->mkdir(Form("%splot/C/", outputDir.c_str()), kTRUE);
  cFig->SaveAs(Form("%splot/C/%s.C", outputDir.c_str(), fileName.c_str()));
  //
  cFig->Clear();
  cFig->Close();
  //
  // Undo the changes
  ws.var("MET")->setBins(oNBins); ws.var("MET")->setRange(oMinRange, oMaxRange);
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
  if (s2.find("WToTau")!=std::string::npos) { s2 = "W#rightarrow#tau"; } else if (s2.find("WW")!=std::string::npos) { s2 = "WW"; }
  else if (s2.find("WZ")!=std::string::npos) { s2 = "WZ"; } else if (s2.find("W")!=std::string::npos) { s2 = "W#rightarrow#mu"; }
  else if (s2.find("DYToTau")!=std::string::npos) { s2 = "Z/#gamma*#rightarrow#tau"; }
  else if (s2.find("DY")!=std::string::npos) { s2 = "Z/#gamma*"; } else if (s2.find("Z")!=std::string::npos) { s2 = "Z"; }
  else if (s2.find("QCD")!=std::string::npos) { s2 = "QCD"; } else if (s2.find("TTbar")!=std::string::npos) { s2 = "t#bar{t}"; }
  s2 = ( s2 + chg );
  if(s3!=""){ label = Form("%s_{%s}^{%s}", s1.c_str(), s2.c_str(), s3.c_str()); } else { label = Form("%s^{%s}", s1.c_str(), s2.c_str()); }
  return;
};


void printElectroWeakMETParameters(TPad& pad, const RooWorkspace& ws, const std::string& pdfName, const uint& drawMode)
{
  pad.cd();
  float xPos = 0.69, yPos = 0.74, dYPos = 0.050, dy = 0.025;
  TLatex t = TLatex(); t.SetNDC(); t.SetTextSize(0.030);
  if (drawMode>0) { dy = 0.065; dYPos *= (1./0.8); t.SetTextSize(0.030*(1./0.8)); }
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
    txtLbl = Form("%s = %.3f#pm%.3f", label.c_str(), it->getValV(), it->getError());
    if (isParAtLimit(*it)) { txtLbl += " (!)"; }
    t.DrawLatex(xPos, yPos-dy, txtLbl.c_str()); dy+=dYPos;
  }
  pad.Update();
  return;
};


void printElectroWeakBinning(TPad& pad, const RooWorkspace& ws, const std::string& dsName, const std::vector< std::string >& text, const uint& drawMode, const int plotStyle)
{
  pad.cd();
  TLatex t = TLatex(); t.SetNDC(); t.SetTextSize(0.030);
  double xPos = 0.20, yPos = 0.89, dYPos = 0.050, dy = 0.035;
  if (plotStyle==1 || plotStyle==2 || plotStyle==3) { xPos = 0.22; yPos = 0.87; dYPos = 0.055; dy = 0.035; }
  t.SetTextSize(0.058*1.25); t.SetTextFont(61); t.DrawLatex(0.78, 0.82, "CMS"); t.SetTextFont(62);
  if (plotStyle!=2) { t.SetTextSize(0.044*1.25); t.SetTextFont(52); t.DrawLatex(0.69, 0.75, "Preliminary"); t.SetTextFont(62); }
  if (drawMode>0) { dy *= (1./0.8); dYPos *= (1./0.8); t.SetTextSize(0.040*(1./0.8)); }
  if (plotStyle==1 || plotStyle==2 || plotStyle==3) {
    t.SetTextSize(0.055*1.25); t.DrawLatex(xPos, 0.82, Form("%s", text[0].c_str())); dy+=dYPos;
    if (text[0].find("Z")==std::string::npos) {
      if (ws.var("Muon_Iso")->getMin()>0.3) { t.SetTextSize(0.035*1.25); t.DrawLatex(xPos, 0.75, Form("%.1f < I^{#mu} < %.1f", ws.var("Muon_Iso")->getMin(), ws.var("Muon_Iso")->getMax())); dy = (0.87-0.75) + dYPos; }
      else if (ws.var("Muon_Pt")->getMin()>18.) { t.SetTextSize(0.035*1.25); t.DrawLatex(xPos, 0.75, "p_{T}^{#mu} > 25 GeV/c"); dy = (0.87-0.75) + dYPos; }
    }
    else return;
  }
  else {
    t.SetTextSize(0.040); t.DrawLatex(xPos, yPos-dy, Form("%s", text[0].c_str())); dy+=dYPos;
  }
  std::vector<std::string> varNameList;
  if (ws.data(dsName.c_str())!=NULL) {
    auto parIt = std::unique_ptr<TIterator>(((RooDataSet*)ws.data(dsName.c_str()))->get()->createIterator());
    for (RooRealVar* it = (RooRealVar*)parIt->Next(); it!=NULL; it = (RooRealVar*)parIt->Next() ) {
      if (std::string(it->GetName())=="MET" || std::string(it->GetName())=="Muon_Pt") continue;
      varNameList.push_back(it->GetName());
    }
  }
  else { varNameList = std::vector<std::string>({ "Muon_Eta" }); }
  for (const auto& varName : varNameList) {
    if ((plotStyle==1 || plotStyle==2 || plotStyle==3) && varName!="Muon_Eta") continue;
    double defaultMin = 0.0 , defaultMax = 100000.0;
    if (varName=="Muon_Eta") { defaultMin = -2.5; defaultMax = 2.5; }
    if (varName=="Muon_Iso") { defaultMin = 0.0; }
    if (ws.var(varName.c_str())) {
      double minVal = ( (ws.var(varName.c_str())->getMin() <= 0.0) ? 0.0 : ws.var(varName.c_str())->getMin() ) ;
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
  //
  if (ws.data(dsName.c_str())!=NULL) {
    const double outTot = ws.data(dsName.c_str())->sumEntries();
    const double outCut = ( (ws.data((dsName+"_FIT").c_str())!=NULL) ? ws.data((dsName+"_FIT").c_str())->sumEntries() : outTot );
    if (outCut != outTot) { t.DrawLatex(xPos, yPos-dy, Form("Loss: (%.4f%%) %.0f evts", ((outTot-outCut)*100.0/outTot), (outTot-outCut))); }
  }
  //
  pad.Update();
  return;
};


void printElectroWeakLegend(TPad& pad, TLegend& leg, const RooPlot& frame, const StringDiMap_t& legInfo, const double size)
{
  pad.cd();
  std::map< std::string , std::string > drawOption = { { "DATA" , "pe" } , { "PDF" , "l" } , { "TEMP" , "f" } };
  const std::vector< std::string > pdfMapOrder = { "WToMu" , "QCD" , "DY" , "WToTau" , "DYToTau" , "TTbar" };
  for (const auto& map : legInfo) {
    if (map.first=="TEMP") {
      for (const auto& pdfM : pdfMapOrder) {
        std::pair< std::string , std::string > elem;
        for (const auto& el : map.second) { if (el.first.find(pdfM)!=std::string::npos) { elem = std::make_pair( el.first , el.second ); break; } }
        if (frame.findObject(elem.first.c_str())) { formatLegendEntry(*leg.AddEntry(frame.findObject(elem.first.c_str()), elem.second.c_str(), drawOption[map.first].c_str()), size); }
      }
    }
    else {
      for (const auto& elem : map.second) {
        if (frame.findObject(elem.first.c_str())) { formatLegendEntry(*leg.AddEntry(frame.findObject(elem.first.c_str()), elem.second.c_str(), drawOption[map.first].c_str()), size); }
      }
    }
  }
  leg.Draw("same");
  pad.Update();
  return;
};


void setRange(RooPlot& frame, const RooWorkspace& ws, const std::string& varName, const std::string& dsName, const bool& setLogScale, const int& nBins)
{
  // Find maximum and minimum points of Plot to rescale Y axis
  TH1D hData, hFit;
  if (ws.data(dsName.c_str())!=NULL) {
    auto h = std::unique_ptr<TH1>(ws.data(dsName.c_str())->createHistogram("hist", *ws.var(varName.c_str()), RooFit::Binning(nBins, frame.GetXaxis()->GetXmin(), frame.GetXaxis()->GetXmax())));
    hData = *((TH1D*)h.get());
  }
  else {
    if (!rooPlotToTH1(hData, hFit, frame)) { std::cout << "[ERROR] Could not find the RooHist from the frame!" << std::endl; return; }
  }
  Double_t YMax = hData.GetBinContent(hData.GetMaximumBin());
  Double_t YMin = 1e99;
  for (int i=1; i<=hData.GetNbinsX(); i++) if (hData.GetBinContent(i)>0) YMin = min(YMin, hData.GetBinContent(i));
  Double_t Yup(0.),Ydown(0.);
  if(setLogScale)
  {
    Yup = YMax*pow((YMax), 0.55);
    Ydown = 1.0;
  }
  else
  {
    Yup = YMax+(YMax-0.0)*0.55;
    Ydown = 0.0;
  }
  //Yup = 1000000.;
  frame.GetYaxis()->SetRangeUser(Ydown,Yup);
  //
  // Draw Lines for the MET range if cut
  if (ws.data((dsName+"_FIT").c_str())!=NULL) {
    double metMin = ws.var("MET")->getMin("METWindow");
    if (metMin > 0.0) {
      auto minline = new TLine(metMin, 0.0, metMin, (setLogScale?(Ydown*TMath::Power((Yup/Ydown),0.5)):(Ydown + (Yup-Ydown)*0.5)));
      minline->SetLineStyle(2); minline->SetLineColor(1); minline->SetLineWidth(3);
      frame.addObject(minline);
    }
    double metMax = ws.var("MET")->getMax("METWindow");
    auto maxline = new TLine(metMax, 0.0, metMax, (setLogScale?(Ydown*TMath::Power((Yup/Ydown),0.5)):(Ydown + (Yup-Ydown)*0.5)));
    maxline->SetLineStyle(2); maxline->SetLineColor(1); maxline->SetLineWidth(3);
    frame.addObject(maxline);
  }
  //
  return;
};


void addGoFToWS(RooWorkspace& ws, const std::string& parName, const double& parValue)
{
  if (ws.var(parName.c_str())) { ws.var(parName.c_str())->setVal(parValue);               }
  else                         { ws.factory(Form("%s[%.6f]", parName.c_str(), parValue)); }
};


bool printGoF(TPad& pad, RooWorkspace& ws, const RooPlot& frame, const string& varLabel, const string& dataLabel, const string& pdfLabel)
{
  //
  if (ws.data(dataLabel.c_str())==NULL && ws.var(Form("testStat_BCChi2_%s", varLabel.c_str()))!=NULL) {
    const int ndof_BCChi2 = int(ws.var(Form("ndofc_BCChi2_%s", varLabel.c_str()))->getVal());
    const double testStat_BCChi2 = ws.var(Form("testStat_BCChi2_%s", varLabel.c_str()))->getVal();
    TLatex t = TLatex(); t.SetNDC(); t.SetTextSize(0.15);
    t.DrawLatex(0.72, 0.83, Form("#chi^{2}/ndof = %.0f / %d ", testStat_BCChi2, ndof_BCChi2));
    return true;
  }
  //
  if (pdfLabel.find("QCDTo")!=std::string::npos && ws.var("ndof_MET")!=NULL && ws.var("chi2_MET")!=NULL) {
    const int ndof_Chi2 = int(ws.var("ndof_MET")->getVal());
    const double testStat_Chi2 = ws.var("chi2_MET")->getVal();
    TLatex t = TLatex(); t.SetNDC(); t.SetTextSize(0.15);
    t.DrawLatex(0.72, 0.83, Form("#chi^{2}/ndof = %.0f / %d ", testStat_Chi2, ndof_Chi2));
    return true;
  }
  //
  if (ws.data(dataLabel.c_str())==NULL) { std::cout << "[ERROR] Dataset " << dataLabel << " was not found!" << std::endl; return false; }
  if (ws.pdf (pdfLabel.c_str() )==NULL) { std::cout << "[ERROR] PDF "     << pdfLabel  << " was not found!" << std::endl; return false; }
  auto dataP = (RooDataSet*)ws.data(dataLabel.c_str());
  auto pdfP  = (RooAbsPdf* )ws.pdf (pdfLabel.c_str() );
  auto varP  = (RooRealVar*)ws.var (varLabel.c_str() );
  //
  // Find curve object
  //auto rFit = (RooCurve*) frame.findObject(0, RooCurve::Class());
  RooCurve* rFit = frame.getCurve(Form("plot_%s", pdfLabel.c_str()));
  if (!rFit) { std::cout << "[ERROR] The latest RooCurve was not found" << std::endl; return false; }
  // Find histogram object
  //auto rData = (RooHist*) frame.findObject(0, RooHist::Class());
  RooHist* rData = frame.getHist(Form("plot_Tot%s", dataLabel.c_str()));
  if (!rData) { std::cout << "[ERROR] The latest RooHist was not found" << std::endl; return false; }
  //
  pad.cd();
  //
  // Unbinned Goodness-of-Fit tests
  //
  RooFit::RooGoF GoF_Unbinned(dataP, pdfP, varP);
  GoF_Unbinned.setRange(varP->getMin(), varP->getMax());
  //GoF_Unbinned.setNtoys(100, true, RooFit::Extended(kTRUE), RooFit::Strategy(2), RooFit::NumCPU(32));
  //
  // Kolmogorov-Smirnov test
  double pvalue_KS = -1., testStat_KS = -1.;
  GoF_Unbinned.KSTest(pvalue_KS, testStat_KS);
  if (testStat_KS>=0.0) {
    std::cout << "[INFO] Using Kolmogorov-Smirnov test gives result " << testStat_KS << " ( " << pvalue_KS << " ) " << std::endl;
    addGoFToWS(ws, Form("testStat_KS_%s", varLabel.c_str()), testStat_KS);
    addGoFToWS(ws, Form("pvalue_KS_%s"  , varLabel.c_str()), pvalue_KS  );
  }
  //
  // Anderson-Darling test
  double pvalue_AD = -1., testStat_AD = -1.;
  GoF_Unbinned.ADTest(pvalue_AD, testStat_AD);
  if (testStat_AD>=0.0) {
    std::cout << "[INFO] Using Anderson-Darling test gives result " << testStat_AD << " ( " << pvalue_AD << " ) " << std::endl;
    addGoFToWS(ws, Form("testStat_AD_%s", varLabel.c_str()), testStat_AD);
    addGoFToWS(ws, Form("pvalue_AD_%s"  , varLabel.c_str()), pvalue_AD  );
  }
  //
  // Binned Goodness-of-Fit tests
  //
  RooFit::RooGoF GoF_Binned(rData, rFit);
  GoF_Binned.setRange(varP->getMin(), varP->getMax());
  GoF_Binned.setRebin(5, false); // We use 5 in approval plots
  //
  // Determine the number of free parameters
  auto parList = std::unique_ptr<RooArgSet>(ws.pdf(pdfLabel.c_str())->getParameters(*ws.data(dataLabel.c_str())));
  const int nFitPar = parList->selectByAttrib("Constant", kFALSE)->getSize();
  //
  // Baker-Cousins chi2 test
  int ndof_BCChi2 = -1;
  double pvalue_BCChi2 = -1., testStat_BCChi2 = -1.;
  GoF_Binned.BCChi2Test(pvalue_BCChi2, testStat_BCChi2, ndof_BCChi2, nFitPar);
  if (ndof_BCChi2>=0.0) {
    std::cout << "[INFO] Using Baker-Cousins chi2 test gives result " << testStat_BCChi2 << "/" << ndof_BCChi2 << " ( " << pvalue_BCChi2 << " ) " << std::endl;
    addGoFToWS(ws, Form("testStat_BCChi2_%s", varLabel.c_str()), testStat_BCChi2);
    addGoFToWS(ws, Form("pvalue_BCChi2_%s"  , varLabel.c_str()), pvalue_BCChi2  );
    addGoFToWS(ws, Form("ndofc_BCChi2_%s"   , varLabel.c_str()), ndof_BCChi2    );
  }
  //
  // Pearson chi2 test
  int ndof_PChi2 = -1;
  double pvalue_PChi2 = -1., testStat_PChi2 = -1.;
  GoF_Binned.PearsonChi2Test(pvalue_PChi2, testStat_PChi2, ndof_PChi2, nFitPar);
  if (ndof_PChi2>=0.0) {
    std::cout << "[INFO] Using Pearson chi2 test gives result " << testStat_PChi2 << "/" << ndof_PChi2 << " ( " << pvalue_PChi2 << " ) " << std::endl;
    addGoFToWS(ws, Form("testStat_PChi2_%s", varLabel.c_str()), testStat_PChi2);
    addGoFToWS(ws, Form("pvalue_PChi2_%s"  , varLabel.c_str()), pvalue_PChi2  );
    addGoFToWS(ws, Form("ndofc_PChi2_%s"   , varLabel.c_str()), ndof_PChi2    );
  }
  //
  // Default RooFit chi2 test (NOT RECOMMENDED)
  int ndof_RooFitChi2 = -1;
  double pvalue_RooFitChi2 = -1., testStat_RooFitChi2 = -1.;
  GoF_Binned.RooFitChi2Test(pvalue_RooFitChi2, testStat_RooFitChi2, ndof_RooFitChi2, nFitPar);
  if (ndof_RooFitChi2>=0.0) {
    std::cout << "[INFO] Using RooFit chi2 test gives result " << testStat_RooFitChi2 << "/" << ndof_RooFitChi2 << " ( " << pvalue_RooFitChi2 << " ) " << std::endl;
    addGoFToWS(ws, Form("testStat_RooFitChi2_%s", varLabel.c_str()), testStat_RooFitChi2);
    addGoFToWS(ws, Form("pvalue_RooFitChi2_%s"  , varLabel.c_str()), pvalue_RooFitChi2  );
    addGoFToWS(ws, Form("ndofc_RoofitChi2_%s"   , varLabel.c_str()), ndof_RooFitChi2    );
  }
  //
  TLatex t = TLatex(); t.SetNDC(); t.SetTextSize(0.15);
  if (ndof_BCChi2>0.) { t.DrawLatex(0.72, 0.83, Form("#chi^{2}/ndof = %.0f / %d ", testStat_BCChi2, ndof_BCChi2)); }
  //
  return true;
};



#endif // #ifndef drawElectroWeakMETPlot_C

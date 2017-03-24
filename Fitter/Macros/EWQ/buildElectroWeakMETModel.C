#ifndef buildElectroWeakMETModel_C
#define buildElectroWeakMETModel_C

#include "../Utilities/initClasses.h"

void setMETModelParameters ( GlobalInfo& info );
TH1* rebinhist             ( const TH1& hist, double xmin, double xmax );
bool histToPdf             ( RooWorkspace& ws, const string& pdfName, const string& dsName, const std::string& var, const std::vector< float >& range );
bool addMETModel           ( RooWorkspace& ws, const std::string& decay, const StrMapMap& models,  const GlobalInfo& info );


bool buildElectroWeakMETModel(RooWorkspace& ws, const StrMapMap& models, GlobalInfo&  info)
{
 // Initialize all the MET Model parameters needed for fitting
  setMETModelParameters(info);
  // Import the MET Models to the local workspace
  if (!addMETModel(ws, "WToMu", models, info)) { return false; }
  // Set Fixed parameters to constant (clean up)
  setFixedVarsToContantVars(ws);
  // save the initial values of the model we've just created
  ws.saveSnapshot("initialParameters",ws.allVars(),kTRUE);
  return true;
};

bool addMETModel(RooWorkspace& ws, const std::string& decay, const StrMapMap& models,  const GlobalInfo& info)
{
  std::string cha = info.Par.at("channel");
  for (auto& col : info.StrV.at("fitSystem")) {
    for (auto& chg : info.StrV.at("fitCharge")) {
      for (auto& obj : info.StrV.at("fitObject")) {
        std::string mainTag = obj + cha + chg;
        std::string mainLabel = mainTag + "_" + col;
        cout << Form("[INFO] Implementing %s MET Model for %s", mainTag.c_str(), col.c_str()) << endl;
        RooArgList pdfList;
        for (auto& model : models.at("Model_"+mainLabel)) {
          std::string tag = model.first + cha + chg;
          std::string label = tag + "_" + col;
          ws.import(RooStringVar(("Model_"+label).c_str(),"",model.second.c_str())); // Save the model name for bookkeeping
          switch(ModelDictionary["MET"][model.second])
            {
            case (int(METModel::CutAndCount)):
              {
                // check that all input parameters are defined
                if (!( 
                      info.Par.count("N_"+label) &&
                      info.Par.count("Cut_"+label)
                       )) {
                  cout << Form("[ERROR] Initial parameters where not found for %s CutAndCount Model in %s", tag.c_str(), col.c_str()) << endl; return false;
                }
                // create the variables for this model
                ws.factory( info.Par.at("N_"+label).c_str() );
                std::string cut = info.Par.at("Cut_"+label);
                ws.import(RooStringVar(("CutAndCount_"+label).c_str(),"",cut.c_str()));
                cout << Form("[INFO] %s in %s added for CutAndCount!", tag.c_str(), col.c_str()) << endl; break;
              }
            case (int(METModel::MultiJetBkg)):
              {
                // check that all input parameters are defined
                if (!( 
                      info.Par.count("N_"+label) &&
                      info.Par.count("Alpha_"+label) &&
                      info.Par.count("Beta_"+label) &&
                      info.Par.count("x0_"+label)
                       )) {
                  cout << Form("[ERROR] Initial parameters where not found for %s MultiJetBkg Model in %s", tag.c_str(), col.c_str()) << endl; return false;
                }
                // create the variables for this model
                /*
                std::string N = info.Par.at("N_"+label);
                if (model.second!=obj) { 
                  N = Form("RooFormulaVar::%s('@0*@1',{%s,%s})", ("N_"+label).c_str(), ("N_"+mainLabel).c_str(), Form("f_%s[0.2, 0.0, 10.0]", label.c_str()));
                  if (!ws.var(("N_"+mainLabel).c_str())) { ws.factory(info.Par.at("N_"+mainLabel).c_str()); }
                }
                ws.factory(N.c_str());
                */
                ws.factory( info.Par.at("N_"+label).c_str()     );
                ws.factory( info.Par.at("Alpha_"+label).c_str() );
                ws.factory( info.Par.at("Beta_"+label).c_str()  );
                ws.factory( info.Par.at("x0_"+label).c_str()    );
                // create the Gaussian PDFs for Constrain fits
                RooArgList pdfConstrains;
                if (info.Par.count("valAlpha_"+label) && info.Par.count("sigAlpha_"+label)) {
                  ws.factory(Form("Gaussian::Constr%s(%s,%s,%s)", Form("Alpha_%s", label.c_str()), Form("Alpha_%s", label.c_str()), info.Par.at("valAlpha_"+label).c_str(), info.Par.at("sigAlpha_"+label).c_str()));
                  pdfConstrains.add( *ws.pdf(Form("ConstrAlpha_%s", label.c_str())) );
                }
                if (info.Par.count("valBeta_"+label) && info.Par.count("sigBeta_"+label)) {
                  ws.factory(Form("Gaussian::Constr%s(%s,%s,%s)", Form("Beta_%s", label.c_str()), Form("Beta_%s", label.c_str()), info.Par.at("valBeta_"+label).c_str(), info.Par.at("sigBeta_"+label).c_str()));
                  pdfConstrains.add( *ws.pdf(Form("ConstrBeta_%s", label.c_str())) );
                }
                if (info.Par.count("valx0_"+label) && info.Par.count("sigx0_"+label)) {
                  ws.factory(Form("Gaussian::Constr%s(%s,%s,%s)", Form("x0_%s", label.c_str()), Form("x0_%s", label.c_str()), info.Par.at("valx0_"+label).c_str(), info.Par.at("sigx0_"+label).c_str()));
                  pdfConstrains.add( *ws.pdf(Form("Constrx0_%s", label.c_str())) );
                }
                if (pdfConstrains.getSize()>0) { ws.import(pdfConstrains, Form("pdfConstr%s", mainLabel.c_str())); }
                // create the PDF
                ws.factory(Form("RooGenericPdf::%s('exp(sqrt(@0+@1)*@2)*pow((@0+@1),@3)', {%s, %s, %s, %s})", Form("pdfMET_%s", label.c_str()), "MET", 
                                Form("x0_%s", label.c_str()), 
                                Form("Beta_%s", label.c_str()), 
                                Form("Alpha_%s", label.c_str())
                                ));
                ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("pdfMETTot_%s", label.c_str()),
                                Form("pdfMET_%s", label.c_str()),
                                Form("N_%s", label.c_str())
                                ));
                ws.pdf(Form("pdfMET_%s", label.c_str()))->setNormRange("METWindow");
                pdfList.add( *ws.pdf(Form("pdfMETTot_%s", label.c_str())) );
                cout << Form("[INFO] %s MultiJetBkg MET PDF in %s added!", tag.c_str(), col.c_str()) << endl; break;
              }
            case (int(METModel::Template)):
              {
                // check that all input parameters are defined
                if (!( 
                      info.Par.count("N_"+label)
                       )) {
                  cout << Form("[ERROR] Initial parameters where not found for %s Template Model in %s", obj.c_str(), col.c_str()) << endl; return false;
                }
                // create the variables for this model
                /*
                std::string N = info.Par.at("N_"+label);
                if (model.second!=obj) { 
                  N = Form("RooFormulaVar::%s('@0*@1',{%s,%s})", ("N_"+label).c_str(), ("N_"+mainLabel).c_str(), Form("f_%s[0.2, 0.0, 10.0]", label.c_str()));
                  if (!ws.var(("N_"+mainLabel).c_str())) { ws.factory(info.Par.at("N_"+mainLabel).c_str()); }
                }
                ws.factory(N.c_str());
                */
                ws.factory( info.Par.at("N_"+label).c_str() );
                // create the PDF
                std::string dsName = ( "d" + chg + "_" + "MC_" + model.first + "_" + info.Par.at("channelDS") + "_" + col );
                int nBins = min(int( round((info.Var.at("MET").at("Max") -  info.Var.at("MET").at("Min"))/info.Var.at("MET").at("binWidth"))*5. ), 1000);
                const std::vector< float > range = { float(nBins) , info.Var.at("MET").at("Min") , info.Var.at("MET").at("Max") };
                if (!histToPdf(ws, Form("pdfMET_%s", label.c_str()), dsName, "MET", range)) { return false; }
                ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("pdfMETTot_%s", label.c_str()),
                                Form("pdfMET_%s", label.c_str()),
                                Form("N_%s", label.c_str())
                                ));
                ws.pdf(Form("pdfMET_%s", label.c_str()))->setNormRange("METWindow");
                pdfList.add( *ws.pdf(Form("pdfMETTot_%s", label.c_str())) );
                cout << Form("[INFO] %s Template MET PDF in %s added!", tag.c_str(), col.c_str()) << endl; break;
              }
            default :
              {
                if (model.second=="") { cout << "[ERROR] MET Model for " << model.second << " was not defined (is empty)!" << endl; return false; }
                else { cout << "[ERROR] Selected MET Model: " << model.second << " has not been implemented" << endl; return false; }
              }
            }
        }
        if (pdfList.getSize()>0) {
          std::string pdfName = ( "pdfMET_Tot" + mainLabel );
          RooAbsPdf *themodel = new RooAddPdf(pdfName.c_str(), pdfName.c_str(), pdfList/*, yieldList*/ );
          ws.import(*themodel);
          delete themodel;
        }
      }
    }
  }
  return true;
};

bool histToPdf(RooWorkspace& ws, const string& pdfName, const string& dsName, const std::string& var, const std::vector< float >& range)
{
  if (ws.pdf(pdfName.c_str())) {
    cout << Form("[INFO] The %s Template has already been created!", pdfName.c_str()) << endl;
    return true; 
  }
  if (ws.data(dsName.c_str())==NULL) { std::cout << "[ERROR] DataSet " << dsName << " was not found!" << std::endl; return false; }
  if (ws.var(var.c_str())==NULL) { std::cout << "[ERROR] Variable " << var << " was not found!" << std::endl; return false; }
  // Create the histogram
  string histName = pdfName;
  histName.replace(histName.find("pdf"), string("pdf").length(), "h");
  TH1D* hist = (TH1D*)ws.data(dsName.c_str())->createHistogram(histName.c_str(), *ws.var(var.c_str()), RooFit::Binning(int(range[0]), range[1], range[2]));
  // Cleaning the input histogram
  // 1) Remove the Under and Overflow bins
  hist->ClearUnderflowAndOverflow();
  // 2) Set negative bin content to zero
  for (int i=0; i<=hist->GetNbinsX(); i++) { if (hist->GetBinContent(i)<0) { hist->SetBinContent(i, 0.0000000001); } }
  // 2) Reduce the range of histogram and rebin it
  TH1* hClean = rebinhist(*hist, range[1], range[2]);
  cout << Form("[INFO] Implementing %s Template", pdfName.c_str()) << endl;
  string dataName = pdfName;
  dataName.replace(dataName.find("pdf"), string("pdf").length(), "dh");
  RooDataHist* dataHist = new RooDataHist(dataName.c_str(), "", *ws.var(var.c_str()), hClean);
  if (dataHist==NULL) { std::cout << "[ERROR] DataHist used to create " << pdfName << " failed!" << std::endl; return false; } 
  if (dataHist->sumEntries()==0) { std::cout << "[ERROR] DataHist used to create " << pdfName << " is empty!" << std::endl; return false; } 
  if (fabs(dataHist->sumEntries()-hClean->GetSumOfWeights())>0.001) { std::cout << "[ERROR] DataHist used to create " << pdfName << "  " << " is invalid!  " << std::endl; return false; } 
  ws.import(*dataHist);
  RooHistPdf* pdf = new RooHistPdf(pdfName.c_str(), pdfName.c_str(), *ws.var(var.c_str()), *((RooDataHist*)ws.data(dataName.c_str())));
  //RooKeysPdf* pdf = new RooKeysPdf(pdfName.c_str(), pdfName.c_str(), *ws.var("ctauErr"), *((RooDataSet*)ws.data(dataName.c_str())),RooKeysPdf::NoMirror, isPbPb?0.4:0.4);
  if (pdf==NULL) { std::cout << "[ERROR] RooHistPDF " << pdfName << " is NULL!" << std::endl; return false; }
  ws.import(*pdf);
  delete pdf;
  delete dataHist;
  hist->Delete();
  return true;
};


TH1* rebinhist(const TH1& hist, double xmin, double xmax)
{
  TH1 *hcopy = (TH1*) hist.Clone("hcopy");
  // range of the new hist
  int imin = hcopy->FindBin(xmin);
  if (imin>=hcopy->GetNbinsX()) imin=1;
  int imax = hcopy->FindBin(0.999999*xmax);
  if (imax<=1) imax=hcopy->GetNbinsX();
  vector<double> newbins;
  newbins.push_back(hcopy->GetBinLowEdge(imin));
  for (int i=imin; i<=imax; i++) {    
    if (hcopy->GetBinContent(i)>0.000000000001) {
      newbins.push_back(hcopy->GetBinLowEdge(i)+hcopy->GetBinWidth(i));
    } else {
      int nrebin=2;
      for (i++; i<=imax; i++) {
        if (hcopy->GetBinContent(i)>0.000000000001) {
          newbins.push_back(hcopy->GetBinLowEdge(i)+hcopy->GetBinWidth(i));
          hcopy->SetBinContent(i,hcopy->GetBinContent(i)/nrebin);
          break;
        }
        nrebin++;
      }
    }
  }
  if (xmin < newbins[1]) newbins[0] = xmin;
  if (xmax > newbins[newbins.size()-2]) newbins[newbins.size()-1] = xmax;
  TH1 *ans = hcopy->Rebin(newbins.size()-1,"hnew",newbins.data());
  delete hcopy;
  return ans;
};

void setMETModelParameters(GlobalInfo& info)
{
  cout << "[INFO] Initializing MET Model parameters based on user input!" << endl;
  std::string cha = info.Par.at("channel");
  for (auto& col : info.StrV.at("fitSystem")) {
    for (auto& chg : info.StrV.at("fitCharge")) {
      double numEntries = info.Var.at("numEntries").at(chg);
      for (auto& obj : (info.Flag.at("incMCTemp") ? info.StrV.at("template") : info.StrV.at("fitObject"))) {
        std::string label = obj + cha + chg + "_" + col;
        std::string inputLabel = label;
        std::vector<std::string> tryChannel = { cha , "" };
        std::vector<std::string> trySystem  = ( info.Flag.at("doPA") ? std::vector<std::string>({col , "PA"}) : std::vector<std::string>({col}) );
        std::vector<std::string> tryCharge  = { chg , "" };
        if (obj=="W" || obj=="QCD"){
          for (auto& tryCha : tryChannel) {
            bool trySuccess = false;
            for (auto& tryCol : trySystem) {
              for (auto& tryChg : tryCharge) {
                if (info.Par.count("Model_"+inputLabel)==0) { inputLabel = (obj + tryCha + tryChg + "_" + tryCol); } else { trySuccess = true; break; }
                if (info.Par.count("Beta_"+inputLabel)==0) { inputLabel = (obj + tryCha + tryChg + "_" + tryCol); } else { trySuccess = true; break; }
              }
              if (trySuccess) break;
            }
            if (trySuccess) break;
          }
        }
        // NUMBER OF EVENTS
        if (info.Par.count("N_"+label)==0 || info.Par.at("N_"+label)=="") {
          info.Par["N_"+label] = Form("%s[%.4f,%.4f,%.4f]", ("N_"+label).c_str(), numEntries, /*-2.0*numEntries*/0.0, 2.0*numEntries);
        }
        // CUTS FOR CUT AND COUNT ALGO
        if (info.Par.count("Cut_"+label)==0 || info.Par.at("Cut_"+label)=="") {
          if (info.Par.count("Cut_"+inputLabel) && info.Par.at("Cut_"+inputLabel)!="") {
            std::string value = info.Par.at("Cut_"+inputLabel).substr( info.Par.at("Cut_"+inputLabel).find("["), info.Par.at("Cut_"+inputLabel).length() );
            info.Par["Cut_"+label] = Form("%s%s", ("Cut_"+label).c_str(), value.c_str());
          }
          else {
            info.Par["Cut_"+label] = "( 20 <= MET )&&( 40 <= Muon_MT )";
          }
        }
        // Multi Jet Model MODEL PARAMETERS
        if (info.Par.count("Beta_"+label)==0 || info.Par.at("Beta_"+label)=="") {
          if (info.Par.count("Beta_"+inputLabel) && info.Par.at("Beta_"+inputLabel)!="") {
            std::string value = info.Par.at("Beta_"+inputLabel).substr( info.Par.at("Beta_"+inputLabel).find("["), info.Par.at("Beta_"+inputLabel).length() );
            info.Par["Beta_"+label] = Form("%s%s", ("Beta_"+label).c_str(), value.c_str());
          }
          else {
            info.Par["Beta_"+label] = Form("%s[%.4f,%.4f,%.4f]", ("Beta_"+label).c_str(), -4.77, -10.0, -0.01);
          }
        }
        if (info.Par.count("Alpha_"+label)==0 || info.Par.at("Alpha_"+label)=="") {
          if (info.Par.count("Alpha_"+inputLabel) && info.Par.at("Alpha_"+inputLabel)!="") {
            std::string value = info.Par.at("Alpha_"+inputLabel).substr( info.Par.at("Alpha_"+inputLabel).find("["), info.Par.at("Alpha_"+inputLabel).length() );
            info.Par["Alpha_"+label] = Form("%s%s", ("Alpha_"+label).c_str(), value.c_str());
          }
          else {
            info.Par["Alpha_"+label] = Form("%s[%.4f,%.4f,%.4f]", ("Alpha_"+label).c_str(), 8.96, 0.0, 20.0);
          }
        }
        if (info.Par.count("x0_"+label)==0 || info.Par.at("x0_"+label)=="") {
          if (info.Par.count("x0_"+inputLabel) && info.Par.at("x0_"+inputLabel)!="") {
            std::string value = info.Par.at("x0_"+inputLabel).substr( info.Par.at("x0_"+inputLabel).find("["), info.Par.at("x0_"+inputLabel).length() );
            info.Par["x0_"+label] = Form("%s%s", ("x0_"+label).c_str(), value.c_str());
          }
          else {
            info.Par["x0_"+label] = Form("%s[%.4f,%.4f,%.4f]", ("x0_"+label).c_str(), 4.17, 0.0, 10.0);
          }
        }
        // Check Parameters for Constrain Fits
        std::vector<std::string > contrainLabel = { "val" , "sig" };
        for (auto& con : contrainLabel) {
          std::string name = Form("%sBeta_%s", con.c_str(), inputLabel.c_str());
          if (info.Par.count(name) && info.Par.at(name)!="") {
            std::string value = info.Par.at(name).substr( info.Par.at(name).find("["), info.Par.at(name).length() );
            info.Par[Form("%sBeta_%s", con.c_str(), label.c_str())] = Form("%sBeta_%s%s", con.c_str(), label.c_str(), value.c_str());
          }
          else break;
        }
        for (auto& con : contrainLabel) {
          std::string name = Form("%sAlpha_%s", con.c_str(), inputLabel.c_str());
          if (info.Par.count(name) && info.Par.at(name)!="") {
            std::string value = info.Par.at(name).substr( info.Par.at(name).find("["), info.Par.at(name).length() );
            info.Par[Form("%sAlpha_%s", con.c_str(), label.c_str())] = Form("%sAlpha_%s%s", con.c_str(), label.c_str(), value.c_str());
          }
          else break;
        }
        for (auto& con : contrainLabel) {
          std::string name = Form("%sx0_%s", con.c_str(), inputLabel.c_str());
          if (info.Par.count(name) && info.Par.at(name)!="") {
            std::string value = info.Par.at(name).substr( info.Par.at(name).find("["), info.Par.at(name).length() );
            info.Par[Form("%sx0_%s", con.c_str(), label.c_str())] = Form("%sx0_%s%s", con.c_str(), label.c_str(), value.c_str());
          }
          else break;
        }
      }
    }
  }
};


#endif // #ifndef buildElectroWeakMETModel_C

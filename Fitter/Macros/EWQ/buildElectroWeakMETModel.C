#ifndef buildElectroWeakMETModel_C
#define buildElectroWeakMETModel_C

#include "../Utilities/initClasses.h"

void setMETModelParameters ( GlobalInfo& info, const std::string& chg );
TH1* rebinhist             ( const TH1& hist, double xmin, double xmax, const std::string type="Old" );
bool histToPdf             ( RooWorkspace& ws, const string& pdfName, const string& dsName, const std::string& var, const std::vector< double >& range );
bool addMETModel           ( RooWorkspace& ws, const std::string& decay, const StringDiMap_t& models, const GlobalInfo& info, const std::string& chg );


bool buildElectroWeakMETModel(RooWorkspace& ws, const StringDiMap_t& models, GlobalInfo&  info, const std::string& chg)
{
 // Initialize all the MET Model parameters needed for fitting
  setMETModelParameters(info, chg);
  // Import the MET Models to the local workspace
  if (!addMETModel(ws, "WToMu", models, info, chg)) { return false; }
  // Set Fixed parameters to constant (clean up)
  setFixedVarsToContantVars(ws);  
  // save the initial values of the model we've just created
  ws.saveSnapshot("initialParameters",ws.allVars(),kTRUE);
  return true;
};

bool addMETModel(RooWorkspace& ws, const std::string& decay, const StringDiMap_t& models, const GlobalInfo& info, const std::string& chg)
{
  std::string cha = info.Par.at("channel");
  for (const auto& col : info.StrV.at("fitSystem")) {
    for (const auto& obj : info.StrV.at("fitObject")) {
      const std::string mainTag = obj + cha + chg;
      const std::string mainLabel = mainTag + "_" + col;
      std::cout << Form("[INFO] Implementing %s MET Model for %s", mainTag.c_str(), col.c_str()) << std::endl;
      RooArgList pdfList;
      for (const auto& model : models.at("Model_"+mainLabel)) {
	const std::string tag = model.first + cha + chg;
	const std::string label = tag + "_" + col;
	RooStringVar tmp; tmp.setVal(model.second.c_str()); tmp.SetTitle(("Model_"+label).c_str()); ws.import(*((TObject*)&tmp), tmp.GetTitle()); // Save the model name for bookkeeping
	// Create Models
	switch(ModelDictionary.at("MET").at(model.second))
	  {
	  case (int(METModel::CutAndCount)):
	    {
	      // check that all input parameters are defined
	      if (!( 
		    info.Par.count("N_"+label) &&
		    info.Par.count("Cut_"+label)
		     )) {
		std::cout << Form("[ERROR] Initial parameters where not found for %s CutAndCount Model in %s", tag.c_str(), col.c_str()) << std::endl; return false;
	      }
	      // create the variables for this model
	      ws.factory( info.Par.at("N_"+label).c_str() );
	      std::string cut = info.Par.at("Cut_"+label);
	      RooStringVar tmp; tmp.setVal(cut.c_str()); tmp.SetTitle(("CutAndCount_"+label).c_str()); ws.import(*((TObject*)&tmp), tmp.GetTitle());
	      std::cout << Form("[INFO] %s in %s added for CutAndCount!", tag.c_str(), col.c_str()) << std::endl; break;
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
		std::cout << Form("[ERROR] Initial parameters where not found for %s MultiJetBkg Model in %s", tag.c_str(), col.c_str()) << std::endl; return false;
	      }
	      // create the variables for this model
	      ws.factory( info.Par.at("N_"+label).c_str() );
	      RooArgList pdfConstrains;
	      std::vector< std::string > varNames = {"Alpha", "Beta", "x0"};
	      for (const auto v : varNames) {
		ws.factory( info.Par.at(v+"_"+label).c_str() );
		// create the Gaussian PDFs for Constrain fits
		if (info.Par.count("val"+v+"_"+label) && info.Par.count("sig"+v+"_"+label)) {
		  ws.factory(Form("Gaussian::Constr%s(%s,%s,%s)", (v+"_"+label).c_str(), (v+"_"+label).c_str(), info.Par.at("val"+v+"_"+label).c_str(), info.Par.at("sig"+v+"_"+label).c_str()));
		  pdfConstrains.add( *ws.pdf(("Constr"+v+"_"+label).c_str()) );
		}
	      }
	      if (pdfConstrains.getSize()>0) { ws.import(pdfConstrains, Form("pdfConstr%s", mainLabel.c_str())); }
	      //
	      // create the PDF
	      ws.factory(Form("RooGenericPdf::%s('TMath::Exp(@2*TMath::Sqrt(@0+@1))*(TMath::Power(TMath::Abs(@0+@1),@3))', {%s, %s, %s, %s})", Form("pdfMET_%s", label.c_str()), "MET", 
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
	      std::cout << Form("[INFO] %s MultiJetBkg MET PDF in %s added!", tag.c_str(), col.c_str()) << std::endl; break;
	    }
	  case (int(METModel::ModifiedRayleigh)):
	    {
	      // check that all input parameters are defined
	      if (!( 
		    info.Par.count("N_"+label) &&
		    info.Par.count("Sigma0_"+label) &&
		    info.Par.count("Sigma1_"+label) &&
		    info.Par.count("Sigma2_"+label)
		     )) {
		std::cout << Form("[ERROR] Initial parameters where not found for %s Modified Rayleigh Model in %s", tag.c_str(), col.c_str()) << std::endl; return false;
	      }
	      // create the variables for this model
	      ws.factory( info.Par.at("N_"+label).c_str() );
	      RooArgList pdfConstrains;
	      std::vector< std::string > varNames = {"Sigma0", "Sigma1", "Sigma2"};
	      for (const auto v : varNames) {
		ws.factory( info.Par.at(v+"_"+label).c_str() );
		// create the Gaussian PDFs for Constrain fits
		if (info.Par.count("val"+v+"_"+label) && info.Par.count("sig"+v+"_"+label)) {
		  ws.factory(Form("Gaussian::Constr%s(%s,%s,%s)", (v+"_"+label).c_str(), (v+"_"+label).c_str(), info.Par.at("val"+v+"_"+label).c_str(), info.Par.at("sig"+v+"_"+label).c_str()));
		  pdfConstrains.add( *ws.pdf(("Constr"+v+"_"+label).c_str()) );
		}
	      }
	      if (pdfConstrains.getSize()>0) { ws.import(pdfConstrains, Form("pdfConstr%s", mainLabel.c_str())); }
	      // create the PDF
	      const double metMax = 50.; // Optimized to reduce corralation (range set to 100 GeV/c)
	      ws.factory(Form("RooFormulaVar::NormMET('(@0 - %.6f)/%.6f',{MET})", metMax, metMax));
	      ws.factory(Form("RooGenericPdf::%s('@0*TMath::Exp(-1.0*((@0*@0)/(2.0*TMath::Power((@1 + (@2)*(@4) + (@3*@3)*((2.0*@4*@4)-1.0)), 2.0))))', {%s, %s, %s, %s, %s})", Form("pdfMET_%s", label.c_str()), "MET", 
			      Form("Sigma0_%s", label.c_str()),
			      Form("Sigma1_%s", label.c_str()),
			      Form("Sigma2_%s", label.c_str()),
			      "NormMET"
			      ));
	      ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("pdfMETTot_%s", label.c_str()),
			      Form("pdfMET_%s", label.c_str()),
			      Form("N_%s", label.c_str())
			      ));
	      ws.pdf(Form("pdfMET_%s", label.c_str()))->setNormRange("METWindow");
	      pdfList.add( *ws.pdf(Form("pdfMETTot_%s", label.c_str())) );
	      std::cout << Form("[INFO] %s Modified Rayleigh MET PDF in %s added!", tag.c_str(), col.c_str()) << std::endl; break;
	    }
	  case (int(METModel::Template)):
	    {
	      // check that all input parameters are defined
	      if (!( 
		    info.Par.count("N_"+label)
		     )) {
		std::cout << Form("[ERROR] %s was not found for %s Template Model in %s", ("N_"+label).c_str(), obj.c_str(), col.c_str()) << std::endl; return false;
	      }                
	      // create the Template
	      std::string dsName = ( "d" + chg + "_" + "MC_" + model.first + "_" + info.Par.at("channelDS") + "_" + col );
	      const std::vector< double > range = { double(ws.var("MET")->getBins("METWindow")) , ws.var("MET")->getMin("METWindow") , ws.var("MET")->getMax("METWindow") };
	      const bool proceed = histToPdf(ws, Form("pdfMET_%s", label.c_str()), dsName, "MET", range);
              if (proceed && (std::abs(info.Var.at("recoMCEntries").at(label)-ws.data(Form("dhMET_%s", label.c_str()))->sumEntries())>0.5)) {
                std::cout << "[WARNING] The number of events in " << Form("dhMET_%s", label.c_str()) << " changed from (" << info.Var.at("recoMCEntries").at(label) << ") to (" <<
                  ws.data(Form("dhMET_%s", label.c_str()))->sumEntries() << ")" << std::endl;
              }
	      //
	      if (proceed) {
		// create the variables for this model
		bool useMCRatio = false;
		const std::string refVarLabel = ( (obj=="W") ? ("NRef_W"+cha+chg+"_"+col) : ("N_W"+cha+chg+"_"+col) );
		if (info.Par.count(refVarLabel)>0) { RooWorkspace tmp; tmp.factory(info.Par.at(refVarLabel).c_str()); if (tmp.var(refVarLabel.c_str())->getVal()>0.0) { useMCRatio = true; } }
		//
		if ( (obj=="W" || obj=="QCD") && (useMCRatio) && (model.first!="QCD") && (model.first!="W") ) {
		  const std::string refLabel    = "W" + cha + chg + "_" + col;
		  const std::string refVarLabel = ( (obj=="W") ? ("NRef_"+refLabel) : ("N_"+refLabel) );
		  if (info.Par.count("N_"+refLabel)==0) { std::cout << "[ERROR] Parameter " << ("N_"+refLabel) << " was not found" << std::endl; return false; }
		  if (!ws.var(("rN_"+label).c_str())  ) { ws.factory(Form("%s[1.0]", ("rN_"+label).c_str())); }
		  if (!ws.var(refVarLabel.c_str())) { ws.factory(info.Par.at(refVarLabel).c_str()); }
		  //
		  RooWorkspace tmp; tmp.factory( info.Par.at("N_"+label).c_str() );
		  ws.var(("rN_"+label).c_str())->setVal( tmp.var(("N_"+label).c_str())->getVal() / ws.var(refVarLabel.c_str())->getVal() );
		  ws.var(("rN_"+label).c_str())->setConstant(kTRUE);
		  //
		  if (info.Par.count("rRN_"+label)>0 && info.Par.at("rRN_"+label)!="") { tmp.factory(info.Par.at("rRN_"+label).c_str()); }
		  if (tmp.var(("rRN_"+label).c_str())!=NULL && tmp.var(("rRN_"+label).c_str())->getVal()!=1.0) {
		    if (!ws.var(("rRN_"+label).c_str())) ws.factory(info.Par.at("rRN_"+label).c_str());
		    ws.var(("rRN_"+label).c_str())->setConstant(kTRUE);
		    ws.factory(Form("RooFormulaVar::%s('@0*@1*@2',{%s,%s,%s})", ("N_"+label).c_str(), ("rRN_"+label).c_str(), ("rN_"+label).c_str(), info.Par.at("N_"+refLabel).c_str()));
		  }
		  else {
		    ws.factory(Form("RooFormulaVar::%s('@0*@1',{%s,%s})", ("N_"+label).c_str(), ("rN_"+label).c_str(), info.Par.at("N_"+refLabel).c_str()));
		  }
		}
		else {
		  ws.factory( info.Par.at("N_"+label).c_str() );
		  if ( (obj=="Z" && model.first!="DY") ) { ws.var(("N_"+label).c_str())->setConstant(kTRUE); }
		}
		// create the PDF
		ws.factory(Form("RooExtendPdf::%s(%s,%s)", Form("pdfMETTot_%s", label.c_str()),
				Form("pdfMET_%s", label.c_str()),
				Form("N_%s", label.c_str())
				));
		ws.pdf(Form("pdfMET_%s", label.c_str()))->setNormRange("METWindow");
		pdfList.add( *ws.pdf(Form("pdfMETTot_%s", label.c_str())) );
		std::cout << Form("[INFO] %s Template MET PDF in %s added!", tag.c_str(), col.c_str()) << std::endl; break;
	      }
	      else {
		std::cout << Form("[INFO] %s Template MET PDF in %s was NOT created!", tag.c_str(), col.c_str()) << std::endl; break;
	      }
	    }
	  default :
	    {
	      if (model.second=="") { std::cout << "[ERROR] MET Model for " << model.second << " was not defined (is empty)!" << std::endl; return false; }
	      else { std::cout << "[ERROR] Selected MET Model: " << model.second << " has not been implemented" << std::endl; return false; }
	    }
	  }
      }
      if (pdfList.getSize()>0) {
	std::string pdfName = ( "pdfMET_Tot" + mainLabel );
        auto themodel = std::unique_ptr<RooAddPdf>(new RooAddPdf(pdfName.c_str(), pdfName.c_str(), pdfList));
        ws.import(*themodel);
      }
    }
  }
  return true;
};

bool histToPdf(RooWorkspace& ws, const string& pdfName, const string& dsName, const std::string& var, const std::vector< double >& range)
{
  //
  if (ws.pdf(pdfName.c_str())) { std::cout << Form("[INFO] The %s Template has already been created!", pdfName.c_str()) << std::endl; return true; }
  std::cout << Form("[INFO] Implementing %s Template", pdfName.c_str()) << std::endl;
  //
  if (ws.data(dsName.c_str())==NULL) { std::cout << "[WARNING] DataSet " << dsName << " was not found!" << std::endl; return false; }
  if (ws.data(dsName.c_str())->numEntries()<=2.0) { std::cout << "[WARNING] DataSet " << dsName << " has too few events!" << std::endl; return false; }
  if (ws.var(var.c_str())==NULL) { std::cout << "[WARNING] Variable " << var << " was not found!" << std::endl; return false; }
  // Create the histogram
  string histName = pdfName;
  histName.replace(histName.find("pdf"), string("pdf").length(), "h");
  std::unique_ptr<TH1D> hist = std::unique_ptr<TH1D>((TH1D*)ws.data(dsName.c_str())->createHistogram(histName.c_str(), *ws.var(var.c_str()), RooFit::Binning(int(range[0]), range[1], range[2])));
  if (hist==NULL) { std::cout << "[WARNING] Histogram " << histName << " is NULL!" << std::endl; return false; }
  // Cleaning the input histogram
  // 1) Remove the Under and Overflow bins
  hist->ClearUnderflowAndOverflow();
  // 2) Set negative bin content to zero
  for (int i=0; i<=hist->GetNbinsX(); i++) { if (hist->GetBinContent(i)<0.0) { hist->SetBinContent(i, 0.0); } }
  // 2) Reduce the range of histogram and rebin it
  hist.reset((TH1D*)rebinhist(*hist, range[1], range[2]));
  if (hist==NULL) { std::cout << "[WARNING] Cleaned Histogram of " << histName << " is NULL!" << std::endl; return false; }
  string dataName = pdfName;
  dataName.replace(dataName.find("pdf"), string("pdf").length(), "dh");
  std::unique_ptr<RooDataHist> dataHist = std::unique_ptr<RooDataHist>(new RooDataHist(dataName.c_str(), "", *ws.var(var.c_str()), hist.get()));
  if (dataHist==NULL) { std::cout << "[WARNING] DataHist used to create " << pdfName << " failed!" << std::endl; return false; } 
  if (dataHist->sumEntries()==0) { std::cout << "[WARNING] DataHist used to create " << pdfName << " is empty!" << std::endl; return false; } 
  if (std::abs(dataHist->sumEntries() - hist->GetSumOfWeights())>0.001) {
    std::cout << "[ERROR] DataHist (" << dataHist->sumEntries() << ")  used to create histogram (" << hist->GetSumOfWeights() << ") for PDF " << pdfName << "  " << " is invalid!  " << std::endl; return false;
  }
  ws.import(*dataHist);
  ws.var(var.c_str())->setBins(int(range[0])); // Bug Fix
  std::unique_ptr<RooHistPdf> pdf = std::unique_ptr<RooHistPdf>(new RooHistPdf(pdfName.c_str(), pdfName.c_str(), *ws.var(var.c_str()), *((RooDataHist*)ws.data(dataName.c_str()))));
  //std::unique_ptr<RooKeysPdf> pdf = std::unique_ptr<RooKeysPdf>(new RooKeysPdf(pdfName.c_str(), pdfName.c_str(), *ws.var(var.c_str()), *((RooDataSet*)ws.data(dsName.c_str())), RooKeysPdf::NoMirror, 0.4));
  if (pdf==NULL) { std::cout << "[WARNING] RooHistPDF " << pdfName << " is NULL!" << std::endl; return false; }
  ws.import(*pdf);
  return true;
};

TH1* rebinhist(const TH1& hist, double xmin, double xmax, const std::string type)
{
  std::unique_ptr<TH1> hcopy = std::unique_ptr<TH1>((TH1*) hist.Clone("hcopy"));
  // range of the new hist
  int imin = hcopy->FindBin(xmin);
  if (imin>=hcopy->GetNbinsX()) imin=1;
  int imax = hcopy->FindBin(0.999999*xmax);
  if (imax<=1) imax=hcopy->GetNbinsX();
  vector<double> newbins;
  newbins.push_back(hcopy->GetBinLowEdge(imin));
  for (int i=imin; i<=imax; i++) {
    if (hcopy->GetBinContent(i)>0.0) {
      newbins.push_back(hcopy->GetBinLowEdge(i)+hcopy->GetBinWidth(i));
    } else {
      int nrebin=2;
      for (i++; i<=imax; i++) {
        if (hcopy->GetBinContent(i)>0.0) {
          newbins.push_back(hcopy->GetBinLowEdge(i)+hcopy->GetBinWidth(i));
          const double newval = (hcopy->GetBinContent(i)/nrebin);
          hcopy->SetBinContent(i, newval);
          if (type=="New") { for (int j=1; j<nrebin; j++) { hcopy->SetBinContent(i-j, newval); } }
          break;
        }
        nrebin++;
      }
    }
  }
  if (type=="Old") {
    if (xmin < newbins[1]) newbins[0] = xmin;
    if (xmax > newbins[newbins.size()-1]) { newbins.push_back(xmax); }
    return hcopy->Rebin(newbins.size()-1, "hnew", newbins.data());
  }
  if (type=="New") {
    return (TH1*)hcopy->Clone("hnew");
  }
  return NULL;
};

void setMETModelParameters(GlobalInfo& info, const std::string& chg)
{
  std::cout << "[INFO] Initializing MET Model parameters based on user input!" << std::endl;
  std::string cha = info.Par.at("channel");
  for (const auto& col : info.StrV.at("fitSystem")) {
    const double numEntries = info.Var.at("numEntries").at(chg);
    for (const auto& mainObj : info.StrV.at("fitObject")) {
      for (const auto& obj : (info.Flag.at("incMCTemp_"+mainObj) ? info.StrV.at("template") : std::vector<std::string>({mainObj}))) {
	const std::string mainCha  = cha , mainChg  = chg , mainCol  = col;
	std::string foundCha = cha , foundChg = chg , foundCol = col;
	std::vector<std::string> tryChannel = { cha , ""  };
	std::vector<std::string> trySystem  = { col , "PA" , ""};
	std::vector<std::string> tryCharge  = { chg , ""  };
	for (const auto& tryCha : tryChannel) {
	  bool trySuccess = false;
	  for (const auto& tryCol : trySystem) {
	    for (const auto& tryChg : tryCharge) {
	      foundCha = tryCha; foundChg = tryChg; foundCol = tryCol;
	      const std::string tryLabel = obj + foundCha + foundChg + ( (tryCol!="") ? ("_" + foundCol) : "" );
	      if ( (info.Var.count(("rRN_"    + tryLabel).c_str())>0) || (info.Var.count(("Cut_"   + tryLabel).c_str())>0) ||
		   (info.Par.count(("Sigma0_" + tryLabel).c_str())>0) || (info.Par.count(("Alpha_" + tryLabel).c_str())>0) ) 
		{ trySuccess = true; break; }
	    }
	    if (trySuccess) break;
	  }
	  if (trySuccess) break;
	}
	//
	const std::string objLabel      = obj + mainCha  + mainChg  + ( (mainCol !="") ? ("_" + mainCol ) : "" );
	const std::string objFoundLabel = obj + foundCha + foundChg + ( (foundCol!="") ? ("_" + foundCol) : "" );
	//
	// NUMBER OF EVENTS
	if (info.Par.count("N_"+objLabel)==0 || info.Par.at("N_"+objLabel)=="") {
	  if (info.Par.count("N_"+objFoundLabel)==0 || info.Par.at("N_"+objFoundLabel)=="") {
	    double value = numEntries;
	    if ( info.Var.count("recoMCEntries")>0 && info.Var.at("recoMCEntries").count(objLabel)>0) { value = info.Var.at("recoMCEntries").at(objLabel); }
	    info.Par["N_"+objLabel] = Form("%s[%.10f,%.10f,%.10f]", ("N_"+objLabel).c_str(), value, 0.0, 10.0*numEntries);
	    if (obj==mainObj) {
	      info.Par["N_"+objLabel]    = Form("%s[%.10f,%.10f,%.10f]", ("N_"+objLabel).c_str(),    numEntries, 0.0, 10.0*numEntries);
	      info.Par["NRef_"+objLabel] = Form("%s[%.10f,%.10f,%.10f]", ("NRef_"+objLabel).c_str(), value,      0.0, 10.0*numEntries);
	    }
	  }
	  else {
	    std::string content = info.Par.at("N_"+objFoundLabel); content = content.substr( content.find("[") );
	    info.Par["N_"+objLabel] = Form("%s%s", ("N_"+objLabel).c_str(), content.c_str());;
	  }
	}
	// FACTOR OF RATIO OF RAW EVENTS (rRN)
	if (info.Var.count("rRN_"+objLabel)==0 || info.Var.at("rRN_"+objLabel).count("Val")==0 || info.Var.at("rRN_"+objLabel).at("Val")==0.0) {
	  if (info.Var.count("rRN_"+objFoundLabel)==0 || info.Var.at("rRN_"+objFoundLabel).count("Val")==0 || info.Var.at("rRN_"+objFoundLabel).at("Val")==0.0) {
	    info.Par["rRN_"+objLabel] = Form("%s[1.0]", ("rRN_"+objLabel).c_str());
	  }
	  else {
	    info.Par["rRN_"+objLabel] = Form("%s[%.10f]", ("rRN_"+objLabel).c_str(), info.Var.at("rRN_"+objFoundLabel).at("Val"));
	  }
	}
	else {
	  info.Par["rRN_"+objLabel] = Form("%s[%.10f]", ("rRN_"+objLabel).c_str(), info.Var.at("rRN_"+objLabel).at("Val"));
	}
	// CUTS FOR CUT AND COUNT ALGO
	if (info.Par.count("Cut_"+objLabel)==0 || info.Par.at("Cut_"+objLabel)=="") {
	  if (info.Par.count("Cut_"+objFoundLabel)==0 || info.Par.at("Cut_"+objFoundLabel)=="") {
	    info.Par["Cut_"+objLabel] = "( 20 <= MET )&&( 40 <= Muon_MT )";
	  }
	  else {
	    info.Par["Cut_"+objLabel] = info.Par.at("Cut_"+objFoundLabel);
	  }
	}
	// QCD MODEL PARAMETERS
	std::vector< std::string > varNames = {"Alpha", "Beta", "x0", "Sigma0", "Sigma1", "Sigma2"};
	for (const auto v : varNames) {
	  if (info.Par.count(v+"_"+objLabel)==0 || info.Par.at(v+"_"+objLabel)=="") {
	    if (info.Par.count(v+"_"+objFoundLabel)==0 || info.Par.at(v+"_"+objFoundLabel)=="") {
	      if (v=="Beta"   ) { info.Par[v+"_"+objLabel] = Form("%s[%.10f,%.10f,%.10f]", (v+"_"+objLabel).c_str(),  -3.100,   -20.000,    0.000); }
	      if (v=="Alpha"  ) { info.Par[v+"_"+objLabel] = Form("%s[%.10f,%.10f,%.10f]", (v+"_"+objLabel).c_str(),   6.000,   -10.000,   30.000); }
	      if (v=="x0"     ) { info.Par[v+"_"+objLabel] = Form("%s[%.10f,%.10f,%.10f]", (v+"_"+objLabel).c_str(),   3.000,     0.000,   30.000); }
	      if (v=="Sigma0" ) { info.Par[v+"_"+objLabel] = Form("%s[%.10f,%.10f,%.10f]", (v+"_"+objLabel).c_str(),  14.800,   -10.000,  100.000); }
	      if (v=="Sigma1" ) { info.Par[v+"_"+objLabel] = Form("%s[%.10f,%.10f,%.10f]", (v+"_"+objLabel).c_str(),   6.500,   -30.000,   50.000); }
	      if (v=="Sigma2" ) { info.Par[v+"_"+objLabel] = Form("%s[%.10f,%.10f,%.10f]", (v+"_"+objLabel).c_str(),   0.500,   -50.000,   50.000); }
	    }
	    else {
	      std::string content = info.Par.at(v+"_"+objFoundLabel); content = content.substr( content.find("[") );
	      info.Par[v+"_"+objLabel] = Form("%s%s", (v+"_"+objLabel).c_str(), content.c_str());
	    }
	  }
	  // Check Parameters for Constrain Fits
	  std::vector<std::string > constrainLabel = { "val" , "sig" };
	  for (const auto& con : constrainLabel) {
	    const std::string name = Form("%s%s_%s", con.c_str(), v.c_str(), objLabel.c_str());
	    if (info.Par.count(name) && info.Par.at(name)!="") {
	      std::string content = info.Par.at(name); content = content.substr( content.find("[") );
	      info.Par[Form("%s%s_%s", con.c_str(), v.c_str(), objLabel.c_str())] = Form("%s%s_%s%s", con.c_str(), v.c_str(), objLabel.c_str(), content.c_str());
	    }
	    else break;
	  }
	}
      }
    }
  }
};


#endif // #ifndef buildElectroWeakMETModel_C

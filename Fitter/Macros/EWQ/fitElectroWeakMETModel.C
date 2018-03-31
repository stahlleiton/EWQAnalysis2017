#ifndef fitElectroWeakMETModel_C
#define fitElectroWeakMETModel_C

#include "./drawElectroWeakMETPlot.C"
#include "./buildElectroWeakMETModel.C"
#include "../Utilities/initClasses.h"
#include "../../../Utilities/EVENTUTILS.h"

void   updateMETParameterRange    ( RooWorkspace& myws );
void   setMETGlobalParameterRange ( RooWorkspace& myws , const GlobalInfo& info );
void   setEWQCutParameters  ( GlobalInfo& info );
bool   setEWQModel          ( StringDiMap_t& model, GlobalInfo&  info );
int    importDataset        ( RooWorkspace& myws  , const std::map<string, RooWorkspace>& inputWS , const GlobalInfo& info, const std::string& chg );
void   setMETFileName       ( string& fileName, string& outputDir, const string& DSTAG, const string& plotLabel, const GlobalInfo& info );


bool fitElectroWeakMETModel( const RooWorkspaceMap_t& inputWorkspaces,    // Workspace with all the input RooDatasets
                             const GlobalInfo& inputInfo,     // Contains information on initial Parameters, cut values, flags, ...
                             const GlobalInfo& userInput,     // Contains information on initial Parameters, cut values, flags, ...
                             const std::string& outputDir,    // Path to output directory
                             // Select the type of datasets to fit
                             const std::string& DSTAG,        // Specifies the name of the dataset to fit
			     const bool saveAll=true
                             // Select the fitting options
                             //string inputFitDir = ""        // Location of the fit results
                             )
{
  // Set up the local workspace and the input information
  RooWorkspaceMap_t myws;
  GlobalInfo info(userInput);
  info.Copy(inputInfo, true); // Copy the user input information (avoid duplicating information in fitter)

  // Check the input settings
  // Figure out the collision system to fit
  info.Flag["fitMC"]     = (DSTAG.find("MC")!=std::string::npos);
  info.Flag.at("fitpPb") = (DSTAG.find("pPb")!=std::string::npos);
  info.Flag.at("fitPbp") = (DSTAG.find("Pbp")!=std::string::npos);
  info.Flag.at("fitPA")  = (DSTAG.find("PA")!=std::string::npos);
  info.StrV["fitSystem"].clear();
  for (const auto& col : info.StrV.at("system")  ) { if (info.Flag.at("fit"+col)) { info.StrV.at("fitSystem").push_back(col); } }
  // If we use Eta CM, then change back to Eta LAB because the variable used in all datasets is the Eta LAB
  if (info.Flag.at("useEtaCM")) {
    const bool ispPb = ( info.Flag.at("fitpPb") || info.Flag.at("fitPA") );
    std::cout << "[INFO] Using Muon Eta at Centre of Mass from " << (ispPb ? "p-Pb" : "Pb-p") << " LAB system"  << std::endl;
    std::cout << "CM: " << info.Var.at("Muon_Eta").at("Min") << "  " << info.Var.at("Muon_Eta").at("Max") << std::endl;
    info.Var.at("Muon_Eta").at("Min") = PA::EtaCMtoLAB(info.Var.at("Muon_Eta").at("Min"), ispPb);
    info.Var.at("Muon_Eta").at("Max") = PA::EtaCMtoLAB(info.Var.at("Muon_Eta").at("Max"), ispPb);
    std::cout << "LAB: " << info.Var.at("Muon_Eta").at("Min") << "  " << info.Var.at("Muon_Eta").at("Max") << std::endl;
  }

  // Set the range of all the parameters
  setEWQCutParameters(info);

  // Set models based on input files
  StringDiMap_t model;
  if (!setEWQModel(model, info)) { return false; }
  
  // Import the all the datasets needed for the fit
  bool doFit = false;
  info.Par["dsName_Pl"] = ("dPl_"+DSTAG);
  info.Par["dsName_Mi"] = ("dMi_"+DSTAG);
  // Add the main dataset to the list
  info.StrV["DSList"].push_back(DSTAG);
  // Check the datasets needed for the template
  for (const auto& tag : info.StrV.at("Tags")) {
    for (const auto& obj : info.StrV.at("fitObject")) {
      if (tag.find(obj)!=std::string::npos && (info.StrV.count("TEMPDS_"+tag)>0)) {
        if (!info.Flag.at("incMCTemp_"+obj)) { std::cout << "[ERROR] The input file for " << tag << " include templates but the input flag " << ("incMCTemp_"+obj) << " is false" << std::endl; return false; }
        break;
      }
    }
  }
  // Add the datasets needed for the template fits = true; } }
  for (const auto& tag : info.StrV.at("Tags")) {
    if (tag.find("Pl_")==std::string::npos) continue; // Only look at one charge since they are symmetric
    if (info.StrV.count("TEMPDS_"+tag)>0) {
      for (const auto& tempDS : info.StrV.at("TEMPDS_"+tag)) {
        if (std::find(info.StrV.at("DSList").begin(), info.StrV.at("DSList").end(), tempDS) == info.StrV.at("DSList").end()) { info.StrV.at("DSList").push_back(tempDS); }
      }
    }
  }

  // Proceed to import the list of datasets
  for (const auto& chg : info.StrV.at("fitCharge")) { 
    if ( !(myws[chg].data(info.Par.at(Form("dsName_%s", chg.c_str())).c_str())) ) {
      const int importID = importDataset(myws.at(chg), inputWorkspaces, info, chg);
      if (importID<0) { return false; }
      else if (importID==0) { doFit = false; }
    }
    info.Var["numEntries"][chg] = myws.at(chg).data(info.Par.at(Form("dsName_%s", chg.c_str())).c_str())->sumEntries(); if (info.Var.at("numEntries").at(chg)<=0) { doFit = false; }
  }

  // Store the number of MC ontries passing analysis cuts
  for (const auto& col : info.StrV.at("fitSystem")) {
    for (const auto& mainObj : info.StrV.at("fitObject")) {
      for (const auto& obj : (info.Flag.at("incMCTemp_"+mainObj) ? info.StrV.at("template") : std::vector<std::string>({mainObj}))) {
        std::string dsLabel = "MC_" + obj + "_" + info.Par.at("channelDS") + "_" + col;
        for (const auto& chg : info.StrV.at("fitCharge")) {
          if (myws.at(chg).data(Form("d%s_%s", chg.c_str(), dsLabel.c_str()))) {
            std::string label = obj + info.Par.at("channel") + chg + "_" + col;
            info.Var["recoMCEntries"][label] = myws.at(chg).data(Form("d%s_%s", chg.c_str(), dsLabel.c_str()))->sumEntries();
          }
        }
      }
    }
  }

  // Store the Luminosity if it is missing
  for (const auto& col : info.StrV.at("fitSystem")) {
    for (const auto& chg : info.StrV.at("fitCharge")) {
      if ( (col=="PA" || col=="pPb") && myws.at(chg).var("Luminosity_pPb")==NULL) { myws.at(chg).factory(Form("Luminosity_pPb[%.4f]", PA::LUMI::Data_pPb)); }
      if ( (col=="PA" || col=="Pbp") && myws.at(chg).var("Luminosity_Pbp")==NULL) { myws.at(chg).factory(Form("Luminosity_Pbp[%.4f]", PA::LUMI::Data_Pbp)); }
    }
  }

  // Set global parameters
  for (const auto& chg : info.StrV.at("fitCharge")) { setMETGlobalParameterRange(myws.at(chg), info); }

  // Update the MET Fit Range
  if (info.Var.at("MET").at("Max")>200.) { for (const auto& chg : info.StrV.at("fitCharge")) { updateMETParameterRange(myws.at(chg), info, chg, DSTAG, -1.); } }

  // Build the Fit Model
  for (const auto& chg : info.StrV.at("fitCharge")) { if (!buildElectroWeakMETModel(myws.at(chg), model, info, chg))  { return false; } }

  // Proceed to Fit and Save the results
  std::string cha = info.Par.at("channel");
  for (const auto& col : info.StrV.at("fitSystem")) {
    for (const auto& chg : info.StrV.at("fitCharge")) {
      for (const auto& obj : info.StrV.at("fitObject")) {
        // Save the info in the workspace
        RooStringVar tmp = RooStringVar();
        if (myws.at(chg).obj("DSTAG")) { ((RooStringVar*)myws.at(chg).obj("DSTAG"))->setVal(DSTAG.c_str()); }
        else { tmp.setVal(DSTAG.c_str()); tmp.SetTitle("DSTAG"); myws.at(chg).import(*((TObject*)&tmp), tmp.GetTitle()); }
        if (myws.at(chg).obj("channel")) { ((RooStringVar*)myws.at(chg).obj("channel"))->setVal(cha.c_str()); }
        else { tmp.setVal(cha.c_str()); tmp.SetTitle("channel"); myws.at(chg).import(*((TObject*)&tmp), tmp.GetTitle()); }
        if (myws.at(chg).obj("fitSystem")) { ((RooStringVar*)myws.at(chg).obj("fitSystem"))->setVal(col.c_str()); }
        else { tmp.setVal(col.c_str()); tmp.SetTitle("fitSystem"); myws.at(chg).import(*((TObject*)&tmp), tmp.GetTitle()); }
        if (myws.at(chg).obj("fitCharge")) { ((RooStringVar*)myws.at(chg).obj("fitCharge"))->setVal(chg.c_str()); }
        else { tmp.setVal(chg.c_str()); tmp.SetTitle("fitCharge"); myws.at(chg).import(*((TObject*)&tmp), tmp.GetTitle()); }
        if (myws.at(chg).obj("fitObject")) { ((RooStringVar*)myws.at(chg).obj("fitObject"))->setVal(obj.c_str()); }
        else { tmp.setVal(obj.c_str()); tmp.SetTitle("fitObject"); myws.at(chg).import(*((TObject*)&tmp), tmp.GetTitle()); }
        // Total PDF Name
        const std::string label = obj + cha + chg + "_" + col;
        const std::string pdfName = ( "pdfMET_Tot" + label );
        // Plot Name
        std::string modelN = info.Par.at("Model_"+label);
        modelN.erase(std::remove(modelN.begin(), modelN.end(), ' '), modelN.end());
        stringReplace( modelN, "[", "_" ); stringReplace( modelN, "]", "" ); stringReplace( modelN, "+", "_" ); stringReplace( modelN, ",", "" ); stringReplace( modelN, ";", "" );
        if (myws.at(chg).obj("modelName")) { ((RooStringVar*)myws.at(chg).obj("modelName"))->setVal(modelN.c_str()); }
        else { tmp.setVal(modelN.c_str()); tmp.SetTitle("modelName"); myws.at(chg).import(*((TObject*)&tmp), tmp.GetTitle()); }
        const std::string plotLabel = label + "_Model_" + modelN;
        // Output File Name
        std::string fileName = "";
        std::string outDir = outputDir;
        setMETFileName(fileName, outDir, DSTAG, plotLabel, info);
        // Dataset Name
        const std::string dsName = ( "d" + chg + "_" + DSTAG );
        // Check if the user wants to do binned fits and proceed to bin the data
	if (info.Flag.count("doBinnedFit")>0 && info.Flag.at("doBinnedFit")) { if (!createBinnedDataset(myws.at(chg))) { return false; } }
        // Set the name of the dataset to fit
        const std::string dsNameFit = ( (myws.at(chg).data((dsName+"_FIT").c_str())!=NULL) ? (dsName+"_FIT") : dsName );
        // check if we have already done this fit. If yes, do nothing and return true.
        bool found =  true; bool skipFit = false;
        std::unique_ptr<RooArgSet> newpars;
        if (myws.at(chg).pdf(pdfName.c_str())) { newpars = std::unique_ptr<RooArgSet>(myws.at(chg).pdf(pdfName.c_str())->getParameters(*myws.at(chg).data(dsName.c_str()))); }
        else { newpars = std::unique_ptr<RooArgSet>(new RooArgSet(myws.at(chg).allVars())); }
        found = found && isFitAlreadyFound(*newpars, Form("%sresult/%s.root", outDir.c_str(), ("FIT_"+fileName).c_str()), pdfName.c_str());
        if (found) {
          std::cout << "[INFO] This fit for " << label << " was already done, so I'll just go to the next one." << std::endl;
          continue;
        }
        // Fit the Datasets
        if (skipFit==false) {
          if (myws.at(chg).pdf(pdfName.c_str())) {
            bool isWeighted = myws.at(chg).data(dsNameFit.c_str())->isWeighted();
	    if (dsName.find("DATA")!=std::string::npos) { isWeighted = false; } // BUG FIX
            int numCores = info.Int.at("numCores");
            RooArgList* pdfConstrains = (RooArgList*)myws.at(chg).genobj(Form("pdfConstr%s", label.c_str()));
            std::unique_ptr<RooFitResult> fitResult;
            bool fitFailed = false;
            if (pdfConstrains!=NULL && pdfConstrains->getSize()>0) {
              std::cout << "[INFO] Fitting with constrain PDFs" << std::endl;
	      auto tmp = myws.at(chg).pdf(pdfName.c_str())->fitTo(*myws.at(chg).data(dsNameFit.c_str()), RooFit::Extended(kTRUE), RooFit::SumW2Error(isWeighted), RooFit::Strategy(2),
								  RooFit::Range("METWindow"), RooFit::ExternalConstraints(*pdfConstrains), RooFit::NumCPU(numCores), RooFit::Save());
	      fitResult.reset(tmp);
	      bool fitFailed = false; for (uint iSt = 0; iSt < fitResult->numStatusHistory(); iSt++) { if (fitResult->statusCodeHistory(iSt)!=0) { fitFailed = true; break; } }
              if (fitFailed) {
                std::cout << std::endl; std::cout << "[INFO] Using Strategy 2 with Minuit2 and Minimizer" << std::endl; std::cout << std::endl;
                tmp = myws.at(chg).pdf(pdfName.c_str())->fitTo(*myws.at(chg).data(dsNameFit.c_str()), RooFit::Extended(kTRUE), RooFit::SumW2Error(isWeighted),
							       RooFit::Minimizer("Minuit2","minimize"), RooFit::Strategy(2),
							       RooFit::Range("METWindow"), RooFit::ExternalConstraints(*pdfConstrains), RooFit::NumCPU(numCores), RooFit::Save());
                fitResult.reset(tmp);
                fitFailed = false; for (uint iSt = 0; iSt < fitResult->numStatusHistory(); iSt++) { if (fitResult->statusCodeHistory(iSt)!=0) { fitFailed = true; break; } }
              }
              if (fitFailed) {
                std::cout << std::endl; std::cout << "[INFO] Using Strategy 2 with Minuit2 and Scan" << std::endl; std::cout << std::endl;
                tmp = myws.at(chg).pdf(pdfName.c_str())->fitTo(*myws.at(chg).data(dsNameFit.c_str()), RooFit::Extended(kTRUE), RooFit::SumW2Error(isWeighted),
							       RooFit::Minimizer("Minuit2","scan"), RooFit::Strategy(2),
							       RooFit::Range("METWindow"), RooFit::ExternalConstraints(*pdfConstrains), RooFit::NumCPU(numCores), RooFit::Save());
                fitResult.reset(tmp);
                fitFailed = false; for (uint iSt = 0; iSt < fitResult->numStatusHistory(); iSt++) { if (fitResult->statusCodeHistory(iSt)!=0) { fitFailed = true; break; } }
              }
	    }
            else {
              auto tmp = myws.at(chg).pdf(pdfName.c_str())->fitTo(*myws.at(chg).data(dsNameFit.c_str()), RooFit::Extended(kTRUE), RooFit::SumW2Error(isWeighted), RooFit::Strategy(2),
								  RooFit::Range("METWindow"), RooFit::NumCPU(numCores), RooFit::Save());
              fitResult.reset(tmp);
              bool fitFailed = false; for (uint iSt = 0; iSt < fitResult->numStatusHistory(); iSt++) { if (fitResult->statusCodeHistory(iSt)!=0) { fitFailed = true; break; } }
              if (fitFailed) {
                std::cout << std::endl; std::cout << "[INFO] Using Strategy 2 with Minuit2 and Minimizer" << std::endl; std::cout << std::endl;
                tmp = myws.at(chg).pdf(pdfName.c_str())->fitTo(*myws.at(chg).data(dsNameFit.c_str()), RooFit::Extended(kTRUE), RooFit::SumW2Error(isWeighted),
							       RooFit::Minimizer("Minuit2","minimize"), RooFit::Strategy(2),
							       RooFit::Range("METWindow"), RooFit::NumCPU(numCores), RooFit::Save());
                fitResult.reset(tmp);
                fitFailed = false; for (uint iSt = 0; iSt < fitResult->numStatusHistory(); iSt++) { if (fitResult->statusCodeHistory(iSt)!=0) { fitFailed = true; break; } }
              }
              if (fitFailed) {
                std::cout << std::endl; std::cout << "[INFO] Using Strategy 2 with Minuit2 and Scan" << std::endl; std::cout << std::endl;
                tmp = myws.at(chg).pdf(pdfName.c_str())->fitTo(*myws.at(chg).data(dsNameFit.c_str()), RooFit::Extended(kTRUE), RooFit::SumW2Error(isWeighted),
							       RooFit::Minimizer("Minuit2","scan"), RooFit::Strategy(2),
							       RooFit::Range("METWindow"), RooFit::NumCPU(numCores), RooFit::Save());
                fitResult.reset(tmp);
                fitFailed = false; for (uint iSt = 0; iSt < fitResult->numStatusHistory(); iSt++) { if (fitResult->statusCodeHistory(iSt)!=0) { fitFailed = true; break; } }
              }
            }
            if (fitResult!=NULL) {
              fitResult->Print("v");
              fitResult->SetTitle(Form("fitResult_%s", pdfName.c_str()));
	      myws.at(chg).import(*fitResult, fitResult->GetTitle());
            }
            else { std::cout << "[ERROR] Fit Result returned by the PDF is NULL!" << std::endl; return false; }
            for (uint iSt = 0; iSt < fitResult->numStatusHistory(); iSt++) {
              if (fitResult->statusCodeHistory(iSt)!=0) {
                std::cout << "[ERROR] Fit failed in " << fitResult->statusLabelHistory(iSt) << " with status " << fitResult->statusCodeHistory(iSt) << " !" << std::endl; return false;
              }
            } 
          }
          else if ( myws.at(chg).obj(("CutAndCount_"+label).c_str()) ) {
            // cut and count
            cout << Form("[INFO] Using the CutAndCount Method with the following cut: %s", info.Par.at("Cut_"+label).c_str()) << std::endl;
            std::unique_ptr<RooDataSet> ds = std::unique_ptr<RooDataSet>((RooDataSet*)myws.at(chg).data(dsName.c_str())->reduce(RooFit::Cut(info.Par.at("Cut_"+label).c_str()),
																RooFit::Name(("CutAndCount_"+dsName).c_str())));
            myws.at(chg).var( ("N_"+label).c_str() )->setVal( ds->sumEntries() );
            myws.at(chg).var( ("N_"+label).c_str() )->setError( sqrt(ds->sumEntries()) );
            myws.at(chg).import(*ds);
            cout << Form("[INFO] Number of events that passed the cut: %.1f", myws.at(chg).var( ("N_"+label).c_str() )->getValV() ) << std::endl;
          }
          else {
            std::cout << "[ERROR] The PDF " << pdfName << " was not found!" << std::endl; return false;
          }
	  if (!drawElectroWeakMETPlot(myws.at(chg), ("PLOT_"+fileName), outDir, info.Flag.at("setLogScale"), 150., true, false)) { return false; }
          myws.at(chg).saveSnapshot("fittedParameters", myws.at(chg).allVars(), kTRUE);
          // Save the results
          if (!saveWorkSpace(myws.at(chg), Form("%sresult/", outDir.c_str()), Form("%s.root", ("FIT_"+fileName).c_str()), saveAll)) { return false; }
        }
      }
    }
  }
  return true;
};


void setEWQCutParameters(GlobalInfo& info)
{
  // Define the MET range
  if (info.Var.at("MET").at("Max")==100000.0) { info.Var.at("MET").at("Max") = 1000.0; }
  // Define the range for the Muon related parameters
  if (info.Flag.at("doMuon")) {
    // Define the Muon PT range
    if (info.Var.at("Muon_Pt").at("Min")==0.0) { info.Var.at("Muon_Pt").at("Min") = 25.0; }
    // Selecting W->MuNu Enhanced Events
    if (info.Flag.at("fitW")) {
      // Define the Muon Iso range
      if (info.Var.at("Muon_Iso").at("Max")==100000.0) { info.Var.at("Muon_Iso").at("Max") = 0.15; }
      if (info.Var.at("Muon_Iso").at("Min")>=info.Var.at("Muon_Iso").at("Max")) { info.Var.at("Muon_Iso").at("Min") = 0.0; }
      // Define the Event Type range
      if (info.Par.at("Event_Type")=="") { info.Par.at("Event_Type") = "Other"; }
    }
    // Selecting QCD Enhanced Events
    if (info.Flag.at("fitQCD")) {
      // Define the MET range (BUG FIX)
      info.Var.at("MET").at("Max") = 150.0;
      // Define the Muon Iso range
      if (info.Var.at("Muon_Iso").at("Min")<0.0) { info.Var.at("Muon_Iso").at("Min") = 0.15; }
      if (info.Var.at("Muon_Iso").at("Max")<=info.Var.at("Muon_Iso").at("Min")) { info.Var.at("Muon_Iso").at("Max") = 100000.0; }
      // Define the Event Type range
      if (info.Par.at("Event_Type")=="") { info.Par.at("Event_Type") = "Other"; }
    }
    // Selecting Z->MuMu Enhanced Events
    if (info.Flag.at("fitZ")) {
      // Define the MET range (BUG FIX)
      info.Var.at("MET").at("Max") = 150.0;
      // Define the Muon Iso range
      if (info.Var.at("Muon_Iso").at("Max")==100000.0) { info.Var.at("Muon_Iso").at("Max") = 0.15; }
      if (info.Var.at("Muon_Iso").at("Min")>=info.Var.at("Muon_Iso").at("Max")) { info.Var.at("Muon_Iso").at("Min") = 0.0; }
      // Define the Event Type range
      if (info.Par.at("Event_Type")=="") { info.Par.at("Event_Type") = "ZToMuMu"; }
    }
    // Print Informstion
    std::cout << "[INFO] Setting MET Magnitud range to min: " << info.Var.at("MET").at("Min") << " and max " << info.Var.at("MET").at("Max") << std::endl;
    std::cout << "[INFO] Setting Muon PT range to min: " << info.Var.at("Muon_Pt").at("Min") << " and max " << info.Var.at("Muon_Pt").at("Max") << std::endl;
    std::cout << "[INFO] Setting Muon Isolation range to min: " << info.Var.at("Muon_Iso").at("Min") << " and max " << info.Var.at("Muon_Iso").at("Max") << std::endl;
    std::cout << "[INFO] Setting Muon Trasnverse Mass range to min: " << info.Var.at("Muon_MT").at("Min") << " and max " << info.Var.at("Muon_MT").at("Max") << std::endl;
    std::cout << "[INFO] Setting Event Type range to : " << (info.Par.at("Event_Type")=="Other" ? "DrellYanVetoed" : info.Par.at("Event_Type")) << std::endl;
  }
  return;
};


bool setEWQModel(StringDiMap_t& model, GlobalInfo&  info)
{
  std::string cha = info.Par.at("channel");
  for (const auto& col : info.StrV.at("fitSystem")) {
    for (const auto& obj : info.StrV.at("fitObject")) {
      for (const auto& chg : info.StrV.at("fitCharge")) {
        std::string label = Form("Model_%s_%s", (obj+cha+chg).c_str(), col.c_str());
        std::string inputLabel = label;
        info.StrV["Tags"].push_back(obj+cha+chg+"_"+col);
        std::vector<std::string> tryChannel = { cha , "" };
        std::vector<std::string> trySystem  = ( info.Flag.at("doPA") ? std::vector<std::string>({col , "PA"}) : std::vector<std::string>({col}) );
        std::vector<std::string> tryCharge  = { chg , "" };
        for (const auto& tryCha : tryChannel) {
          bool trySuccess = false;
          for (const auto& tryCol : trySystem) {
            for (const auto& tryChg : tryCharge) {
              if (info.Par.count(inputLabel)==0) { inputLabel = ("Model_" + obj + tryCha + tryChg + "_" + tryCol); } else { trySuccess = true; break; }
            }
            if (trySuccess) break;
          }
          if (trySuccess) break;
        }
        if (info.Par.count(inputLabel)>0) {
          std::string value = info.Par.at(inputLabel);
          info.Par[label] = value;
          std::vector<std::string> k;
          if (value.find("+")!=std::string::npos) { if (!splitString(value, "+", k)) { return false; } }
          else { k.push_back(value); }
          for (auto& kk : k) {
            std::string modelName = kk;
            std::vector<std::string> p;
            if (kk.find("[")!=std::string::npos) { 
              modelName = kk.substr(0, kk.find("["));
              kk.erase(0, kk.find("[")+std::string("[").length());
              if (kk.find("]")==std::string::npos) { std::cout << "[ERROR] Missing ']' in model: " << value << std::endl; return false; }
              kk.erase(kk.find("]"),kk.length());
              if (kk.find(";")!=std::string::npos) { if (!splitString(kk, ";", p)) { return false; } }
              else if (kk.find(",")!=std::string::npos) { if (!splitString(kk, ",", p)) { return false; } }
              else { p.push_back(kk); }
            }
            else { p.push_back(obj); }
            for (const auto& ll : p) {
              std::string objectName  = Form("%s", (ll+cha+chg).c_str());;
              if (info.Flag.at("fitMC") && ll!=obj && modelName=="TEMP") continue;
              if (modelName=="TEMP") { modelName = "Template";    }
              if (modelName=="MJET") { modelName = "MultiJetBkg"; }
              if (ModelDictionary.at("MET").at(modelName)==0) {
                std::cout << "[ERROR] The MET " << obj << " model: " << modelName << " is invalid" << std::endl; return false;
              }
              model[label][ll] = modelName;
              if (modelName == "Template") {
                std::string dsTag = ( "MC_" + ll + "_" + info.Par.at("channelDS") + "_" + col );
                info.StrV[Form("TEMPDS_%s_%s", (obj+cha+chg).c_str(), col.c_str())].push_back(dsTag);
              }
            }
          }
        } else {
          std::cout << "[ERROR] " << (obj+cha+chg) << " " << "MET" << " model for " << col << " was not found in the initial parameters!" << std::endl; return false;
        }
      }
    }
  }
  return true;
};


int importDataset(RooWorkspace& myws  , const std::map<string, RooWorkspace>& inputWS , const GlobalInfo& info, const std::string& chg)
{
  // Define the selection string
  std::string cutDS = "";
  if (info.Par.at("Event_Type")!="") { cutDS += Form("(Event_Type==Event_Type::%s)&&", info.Par.at("Event_Type").c_str()); }
  for (const auto var : info.Var) {
    if (var.first!="MET" && var.first!="Muon_Pt" && var.first!="Muon_Eta" && var.first!="Muon_Iso" && var.first!="Muon_MT") continue; // Only cut on this variables
    if (var.second.at("Min")==var.second.at("Max")) { cutDS += Form("(%s == %g)", var.first.c_str(), var.second.at("Max")); }
    else { cutDS += Form("(%g <= %s && %s < %g)", var.second.at("Min"), var.first.c_str(), var.first.c_str(), var.second.at("Max")); }
    cutDS += "&&";
  }
  cutDS.erase(cutDS.size()-string("&&").length(), cutDS.size());
  TObjString tmp; tmp.SetString(cutDS.c_str()); myws.import(*((TObject*)&tmp), "Cut_DataSet"); // Save the cut expression for bookkeeping
  std::cout << "[INFO] Importing local RooDataSets with cuts: " << cutDS << std::endl;
  // Reduce and import the datasets
  for (const auto& label : info.StrV.at("DSList")) {
    // Extract the RooDatasets
    std::string dsType = "COR";
    if (inputWS.at(label).data(Form("dPl_%s_%s", dsType.c_str(), label.c_str()))==NULL) { dsType = "LUM"; }
    if (inputWS.at(label).data(Form("dPl_%s_%s", dsType.c_str(), label.c_str()))==NULL) { dsType = "SET"; }
    if (inputWS.at(label).data(Form("dPl_%s_%s", dsType.c_str(), label.c_str()))==NULL) { dsType = "RAW"; }
    if (inputWS.at(label).data(Form("dPl_%s_%s", dsType.c_str(), label.c_str()))==NULL) { std::cout << "[ERROR] Sample " << (dsType + "_" + label) << " was not found!" << std::endl; return -1; }
    const std::string extLabel = dsType + "_" + label;
    std::cout << "[INFO] Importing local RooDataSet " << extLabel << std::endl;
    //
    if (chg!="") {
      const std::string dsName = Form("d%s_%s", chg.c_str(), label.c_str());
      if (myws.data(dsName.c_str())!=NULL) continue;
      if ( inputWS.count(label)==0 || !(inputWS.at(label).data(Form("d%s_%s", chg.c_str(), extLabel.c_str())))){ 
        std::cout << "[ERROR] The dataset " <<  Form("d%s_%s", chg.c_str(), extLabel.c_str()) << " was not found!" << std::endl;
        return -1;
      }
      std::unique_ptr<RooDataSet> data = std::unique_ptr<RooDataSet>((RooDataSet*)inputWS.at(label).data(Form("d%s_%s", chg.c_str(), extLabel.c_str()))->reduce(cutDS.c_str()));
      if (data==NULL || data->sumEntries()==0){
        if (extLabel.find("MC_")!=std::string::npos) {
          std::cout << "[WARNING] No events from dataset " <<  Form("d%s_%s", chg.c_str(), extLabel.c_str()) << " passed the kinematic cuts!" << std::endl;
        }
        else { std::cout << "[ERROR] No events from dataset " <<  Form("d%s_%s", chg.c_str(), extLabel.c_str()) << " passed the kinematic cuts!" << std::endl; return -1; }
      }
      else {
        data->SetName(dsName.c_str());
        myws.import(*data);
      }
      std::cout << "[INFO] " << Form("%.0f", data->sumEntries()) << " weighted entries imported from local RooDataSet " << dsName << std::endl;
      //
      const std::vector< std::string > strLabels = { "METType" , "CorrectionApplied" , "RecoilMethod" };
      for (const auto& strL : strLabels) { if (inputWS.at(label).obj(strL.c_str()) && !myws.obj(strL.c_str())) { myws.import(*inputWS.at(label).obj(strL.c_str()), strL.c_str()); } }
      //
      // Set the range of each global parameter in the local roodataset
      if (myws.data(dsName.c_str())!=NULL) {
        const RooArgSet* row = myws.data(dsName.c_str())->get();
        for (const auto& var : info.Var) {
          if ( (var.first!="Event_Type") && row->find(Form("%s", var.first.c_str())) ) {
            ((RooRealVar*)row->find(Form("%s", var.first.c_str())))->setMin(var.second.at("Min"));
            ((RooRealVar*)row->find(Form("%s", var.first.c_str())))->setMax(var.second.at("Max"));
          }
        }
      }
    }
  }
  // Check if the user wants to use the Center of Mass Eta
  if (info.Flag.at("useEtaCM")) { myws.factory("useEtaCM[1.0]"); }
  // Set the range of each global parameter in the local workspace
  for (const auto& var : info.Var) {
    if ( myws.var(Form("%s", var.first.c_str())) && (var.second.count("Min")>0) ) {
      myws.var(Form("%s", var.first.c_str()))->setMin(var.second.at("Min"));
      myws.var(Form("%s", var.first.c_str()))->setMax(var.second.at("Max"));
    }
    else if ( (myws.var(Form("%s", var.first.c_str()))==NULL) && (var.second.count("Val")>0) ) {
      myws.factory(Form("%s[%.10f]", var.first.c_str(), var.second.at("Val")));
    }
  }
  if (info.Flag.at("useEtaCM")) {
    const bool ispPb = ( info.Flag.at("fitpPb") || info.Flag.at("fitPA") );
    std::cout << "[INFO] Analyzing bin: " << Form("%g < Muon_EtaCM < %g", PA::EtaLABtoCM(myws.var("Muon_Eta")->getMin(), ispPb), PA::EtaLABtoCM(myws.var("Muon_Eta")->getMax(), ispPb)) << std::endl;
  }
  else { std::cout << "[INFO] Analyzing bin: " << Form("%g < Muon_Eta < %g"  , myws.var("Muon_Eta")->getMin()  , myws.var("Muon_Eta")->getMax())   << std::endl; }
  return 1;
};


void setMETGlobalParameterRange(RooWorkspace& myws, const GlobalInfo& info)
{
  myws.var("MET")->setRange("METWindow", info.Var.at("MET").at("Min"), info.Var.at("MET").at("Max"));
  const int nBins = std::min(int( std::round((info.Var.at("MET").at("Max") - info.Var.at("MET").at("Min"))/info.Var.at("MET").at("binWidth")) ), 2000);
  myws.var("MET")->setBins(nBins, "METWindow");
  myws.var("MET")->setBins(nBins);
  return;
};


void setMETFileName(string& fileName, string& outputDir, const string& DSTAG, const string& plotLabel, const GlobalInfo& info)
{
  std::string dsTag  = DSTAG; dsTag.erase(dsTag.find("_MUON"), dsTag.length());
  std::string colTag = DSTAG; colTag = colTag.substr(colTag.find_last_of("_")+1);
  const std::string metTAG = "MET" + info.Par.at("METType");
  const std::string objTag = info.StrV.at("fitObject")[0];
  outputDir = Form("%s%s/%s/%s/%s/", outputDir.c_str(), metTAG.c_str(), dsTag.c_str(), objTag.c_str(), colTag.c_str());
  const bool ispPb = ( info.Flag.at("fitpPb") || info.Flag.at("fitPA") );
  const double etaMin = ( info.Flag.at("useEtaCM") ? PA::EtaLABtoCM(info.Var.at("Muon_Eta").at("Min"), ispPb) : info.Var.at("Muon_Eta").at("Min") );
  const double etaMax = ( info.Flag.at("useEtaCM") ? PA::EtaLABtoCM(info.Var.at("Muon_Eta").at("Max"), ispPb) : info.Var.at("Muon_Eta").at("Max") );
  const double isoMin = ( ( info.Var.at("Muon_Iso").at("Min") <= 0.0 ) ? 0.0 : info.Var.at("Muon_Iso").at("Min") );
  fileName = Form("%s_%s_%s_%s_%.0f_%.0f_MuIso_%.0f_%.0f", "MET", 
                  dsTag.c_str(),
                  plotLabel.c_str(),
                  ( info.Flag.at("useEtaCM") ? "MuEtaCM" : "MuEta" ),
                  (etaMin*100.0), (etaMax*100.0),
                  (isoMin*100.0), (info.Var.at("Muon_Iso").at("Max")*100.0)
                  );
  return;
};


#endif // #ifndef fitElectroWeakMETModel_C

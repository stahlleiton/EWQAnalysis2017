#ifndef fitElectroWeakMETModel_C
#define fitElectroWeakMETModel_C

#include "./drawElectroWeakMETPlot.C"
#include "./buildElectroWeakMETModel.C"
#include "../Utilities/initClasses.h"
#include "../Utilities/EVENTUTILS.h"

void   setMETGlobalParameterRange ( RooWorkspace& myws , GlobalInfo& info );
void   setEWQCutParameters  ( GlobalInfo& info );
bool   setEWQModel          ( StrMapMap& model, GlobalInfo&  info );
int    importDataset        ( RooWorkspace& myws  , const std::map<string, RooWorkspace>& inputWS , const GlobalInfo& info);
void   setMETFileName       ( string& fileName, string& outputDir, const string& DSTAG, const string& plotLabel, const GlobalInfo& info );
void   addLuminosity        ( GlobalInfo& info );


bool fitElectroWeakMETModel( const RooWorkspaceMap& inputWorkspaces,    // Workspace with all the input RooDatasets
                             const GlobalInfo& inputInfo,     // Contains information on initial Parameters, cut values, flags, ...
                             const GlobalInfo& userInput,     // Contains information on initial Parameters, cut values, flags, ...
                             const std::string& outputDir,    // Path to output directory
                             // Select the type of datasets to fit
                             const std::string& DSTAG//,                  // Specifies the name of the dataset to fit
                             // Select the fitting options
                             //string inputFitDir = ""        // Location of the fit results
                             )
{

  // Set up the local workspace and the input information
  RooWorkspace myws;
  GlobalInfo info(userInput);
  info.Copy(inputInfo, true); // Copy the user input information (avoid duplicating information in fitter)

  // Check the input settings
  // Figure out the collision system to fit
  info.Flag["fitMC"]  = (DSTAG.find("MC")!=std::string::npos);
  info.Flag["fitpPb"] = (DSTAG.find("pPb")!=std::string::npos);
  info.Flag["fitPbp"] = (DSTAG.find("Pbp")!=std::string::npos);
  info.Flag["fitPA"]  = (DSTAG.find("PA")!=std::string::npos);
  info.Flag["doPA"]   = (info.Flag.at("fitPbp") && info.Flag.at("fitpPb"));
  for (auto& col : info.StrV.at("system")  ) { if (info.Flag.at("fit"+col)) { info.StrV["fitSystem"].clear(); info.StrV.at("fitSystem").push_back(col); } }

  // Add Luminosity Variables
  addLuminosity(info);

  // Set the range of all the parameters
  setEWQCutParameters(info);

  // Set models based on input files
  StrMapMap model;
  if (!setEWQModel(model, info)) { return false; }
  
  // Import the all the datasets needed for the fit
  bool doFit = false;
  info.Par["dsName_Pl"] = ("dPl_"+DSTAG);
  info.Par["dsName_Mi"] = ("dMi_"+DSTAG);
  // Add the main dataset to the list
  info.StrV["DSList"].push_back(DSTAG);
  // Add the datasets needed for the template fits
  if (info.Flag.at("incMCTemp")) {
    for(auto& tag : info.StrV.at("Tags")) {
      if (tag.find("Pl_")==std::string::npos) continue; // Only look at one charge since they are symmetric
      if (info.StrV.count("TEMPDS_"+tag)>0) {
        for (auto& tempDS : info.StrV.at("TEMPDS_"+tag)) {
          if (std::find(info.StrV.at("DSList").begin(), info.StrV.at("DSList").end(), tempDS)== info.StrV.at("DSList").end()) { info.StrV.at("DSList").push_back(tempDS); }
        }
      }
    }
  }
  // Proceed to import the list of datasets
  if ( !(myws.data(info.Par.at("dsName_Pl").c_str())) || !(myws.data(info.Par.at("dsName_Mi").c_str())) ) {
    int importID = importDataset(myws, inputWorkspaces, info);
    if (importID<0) { return false; }
    else if (importID==0) { doFit = false; }
  }
  info.Var["numEntries"]["Pl"] = myws.data(info.Par.at("dsName_Pl").c_str())->sumEntries(); if (info.Var.at("numEntries").at("Pl")<=0) { doFit = false; }
  info.Var["numEntries"]["Mi"] = myws.data(info.Par.at("dsName_Mi").c_str())->sumEntries(); if (info.Var.at("numEntries").at("Mi")<=0) { doFit = false; }

  // Set global parameters
  setMETGlobalParameterRange(myws, info);

  // Build the Fit Model
  if (!buildElectroWeakMETModel(myws, model, info))  { return false; }

  // Proceed to Fit and Save the results
  std::string cha = info.Par.at("channel");
  for (auto& col : info.StrV.at("fitSystem")) {
    for (auto& chg : info.StrV.at("fitCharge")) {
      for (auto& obj : info.StrV.at("fitObject")) {
        // Save the info in the workspace
        TObjString tmp = TObjString();
        if (myws.obj("DSTAG")) { ((TObjString*)myws.obj("DSTAG"))->SetString(DSTAG.c_str()); }
        else { tmp.SetString(DSTAG.c_str()); myws.import(*((TObject*)&tmp), "DSTAG"); }
        if (myws.obj("channel")) { ((TObjString*)myws.obj("channel"))->SetString(cha.c_str()); }
        else { tmp.SetString(cha.c_str()); myws.import(*((TObject*)&tmp), "channel"); }
        if (myws.obj("fitSystem")) { ((TObjString*)myws.obj("fitSystem"))->SetString(col.c_str()); }
        else { tmp.SetString(col.c_str()); myws.import(*((TObject*)&tmp), "fitSystem"); }
        if (myws.obj("fitCharge")) { ((TObjString*)myws.obj("fitCharge"))->SetString(chg.c_str()); }
        else { tmp.SetString(chg.c_str()); myws.import(*((TObject*)&tmp), "fitCharge"); }
        if (myws.obj("fitObject")) { ((TObjString*)myws.obj("fitObject"))->SetString(obj.c_str()); }
        else { tmp.SetString(obj.c_str()); myws.import(*((TObject*)&tmp), "fitObject"); }
        // Total PDF Name
        std::string label = obj + cha + chg + "_" + col;
        std::string pdfName = ( "pdfMET_Tot" + label );
        // Plot Name
        std::string modelN = info.Par.at("Model_"+label);
        modelN.erase(std::remove(modelN.begin(), modelN.end(), ' '), modelN.end());
        stringReplace( modelN, "[", "_" ); stringReplace( modelN, "]", "" ); stringReplace( modelN, "+", "_" ); stringReplace( modelN, ",", "" ); stringReplace( modelN, ";", "" );
        if (myws.obj("modelName")) { ((TObjString*)myws.obj("modelName"))->SetString(modelN.c_str()); }
        else { tmp.SetString(modelN.c_str()); myws.import(*((TObject*)&tmp), "modelName"); }
        std::string plotLabel = label + "_Model_" + modelN;
        // Output File Name
        std::string fileName = "";
        std::string outDir = outputDir;
        setMETFileName(fileName, outDir, DSTAG, plotLabel, info);
        // Dataset Name
        std::string dsName = ( "d" + chg + "_" + DSTAG );
        // check if we have already done this fit. If yes, do nothing and return true.
        bool found =  true; bool skipFit = false;
        RooArgSet newpars;
        if (myws.pdf(pdfName.c_str())) { newpars = *(myws.pdf(pdfName.c_str())->getParameters(*myws.data(dsName.c_str()))); }
        else { newpars = myws.allVars(); }
        found = found && isFitAlreadyFound(newpars, Form("%sresult/%s.root", outDir.c_str(), ("FIT_"+fileName).c_str()), pdfName.c_str());
        if (found) {
          cout << "[INFO] This fit was already done, so I'll just go to the next one." << endl;
          return true;
        }
        // Fit the Datasets
        if (skipFit==false) {
          if (myws.pdf(pdfName.c_str())) {
            bool isWeighted = myws.data(dsName.c_str())->isWeighted();
            int numCores = info.Int.at("numCores");
            RooArgList* pdfConstrains = (RooArgList*)myws.genobj(Form("pdfConstr%s", label.c_str()));
            RooFitResult* fitResult;
            if (pdfConstrains!=NULL && pdfConstrains->getSize()>0) {
              fitResult = myws.pdf(pdfName.c_str())->fitTo(*myws.data(dsName.c_str()), RooFit::Extended(kTRUE), RooFit::SumW2Error(isWeighted), 
                                                           RooFit::Range("METWindow"), RooFit::ExternalConstraints(*pdfConstrains), RooFit::NumCPU(numCores), RooFit::Save());
            }
            else {
              fitResult = myws.pdf(pdfName.c_str())->fitTo(*myws.data(dsName.c_str()), RooFit::Extended(kTRUE), RooFit::SumW2Error(isWeighted), 
                                                           RooFit::Range("METWindow"), RooFit::NumCPU(numCores), RooFit::Save());
            }
            fitResult->Print("v");
            myws.import(*fitResult, Form("fitResult_%s", pdfName.c_str()));
          }
          else if ( myws.obj(("CutAndCount_"+label).c_str()) ) {
            // cut and count
            cout << Form("[INFO] Using the CutAndCount Method with the following cut: %s", info.Par.at("Cut_"+label).c_str()) << std::endl;
            RooDataSet* ds = (RooDataSet*)myws.data(dsName.c_str())->reduce(RooFit::Cut(info.Par.at("Cut_"+label).c_str()), RooFit::Name(("CutAndCount_"+dsName).c_str()));
            myws.var( ("N_"+label).c_str() )->setVal( ds->sumEntries() );
            myws.var( ("N_"+label).c_str() )->setError( sqrt(ds->sumEntries()) );
            myws.import(*ds);
            cout << Form("[INFO] Number of events that passed the cut: %.1f", myws.var( ("N_"+label).c_str() )->getValV() ) << std::endl;
          }
          else {
            std::cout << "[ERROR] The PDF " << pdfName << " was not found!" << std::endl; return false;
          }
          int nBins = min(int( round((info.Var.at("MET").at("Max") - info.Var.at("MET").at("Min"))/info.Var.at("MET").at("binWidth")) ), 1000);
          if (!drawElectroWeakMETPlot(myws, ("PLOT_"+fileName), outDir, nBins)) { return false; }
          myws.saveSnapshot("fittedParameters",myws.allVars(),kTRUE);
          // Save the results
          saveWorkSpace(myws, Form("%sresult/", outDir.c_str()), Form("%s.root", ("FIT_"+fileName).c_str()));
        }
      }
    }
  }
  return true;
};


void setEWQCutParameters(GlobalInfo& info)
{
  // Define the MET range
  if (info.Var.at("MET").at("Max")==100000.0) { info.Var.at("MET").at("Max") = 100.0; }
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
      // Define the Muon Iso range
      if (info.Var.at("Muon_Iso").at("Min")==0.0) { info.Var.at("Muon_Iso").at("Min") = 0.15; }
      if (info.Var.at("Muon_Iso").at("Max")<=info.Var.at("Muon_Iso").at("Min")) { info.Var.at("Muon_Iso").at("Max") = 100000.0; }
      // Define the Event Type range
      if (info.Par.at("Event_Type")=="") { info.Par.at("Event_Type") = "Other"; }
    }
    // Selecting DYZ->MuMu Enhanced Events
    if (info.Flag.at("fitDYZ")) {
      // Define the Muon Iso range
      if (info.Var.at("Muon_Iso").at("Max")==100000.0) { info.Var.at("Muon_Iso").at("Max") = 0.15; }
      if (info.Var.at("Muon_Iso").at("Min")>=info.Var.at("Muon_Iso").at("Max")) { info.Var.at("Muon_Iso").at("Min") = 0.0; }
      // Define the Event Type range
      if (info.Par.at("Event_Type")=="") { info.Par.at("Event_Type") = "DYZToMuMu"; }
    }
    // Print Informstion
    cout << "[INFO] Setting MET Magnitud range to min: " << info.Var.at("MET").at("Min") << " and max " << info.Var.at("MET").at("Max") << endl;
    cout << "[INFO] Setting Muon PT range to min: " << info.Var.at("Muon_Pt").at("Min") << " and max " << info.Var.at("Muon_Pt").at("Max") << endl;
    cout << "[INFO] Setting Muon Isolation range to min: " << info.Var.at("Muon_Iso").at("Min") << " and max " << info.Var.at("Muon_Iso").at("Max") << endl;
    cout << "[INFO] Setting Muon Trasnverse Mass range to min: " << info.Var.at("Muon_MT").at("Min") << " and max " << info.Var.at("Muon_MT").at("Max") << endl;
    cout << "[INFO] Setting Event Type range to : " << info.Par.at("Event_Type") << endl;
  }
  return;
};


bool setEWQModel(StrMapMap& model, GlobalInfo&  info)
{
  std::string cha = info.Par.at("channel");
  for (auto& col : info.StrV.at("fitSystem")) {
    for (auto& obj : info.StrV.at("fitObject")) {
      for (auto& chg : info.StrV.at("fitCharge")) {
        std::string label = Form("Model_%s_%s", (obj+cha+chg).c_str(), col.c_str());
        std::string inputLabel = label;
        info.StrV["Tags"].push_back(obj+cha+chg+"_"+col);
        std::vector<std::string> tryChannel = { cha , "" };
        std::vector<std::string> trySystem  = ( info.Flag.at("doPA") ? std::vector<std::string>({col , "PA"}) : std::vector<std::string>({col}) );
        std::vector<std::string> tryCharge  = { chg , "" };
        for (auto& tryCha : tryChannel) {
          bool trySuccess = false;
          for (auto& tryCol : trySystem) {
            for (auto& tryChg : tryCharge) {
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
            for (auto& ll : p) {
              std::string objectName  = Form("%s", (ll+cha+chg).c_str());;
              if (modelName=="TEMP") { modelName = "Template";    }
              if (modelName=="MJET") { modelName = "MultiJetBkg"; }
              if (ModelDictionary["MET"][modelName]==0) {
                cout << "[ERROR] The MET " << obj << " model: " << modelName << " is invalid" << endl; return false;
              }
              model[label][ll] = modelName;
              if (modelName == "Template") {
                std::string dsTag = ( "MC_" + ll + "_" + info.Par.at("channelDS") + "_" + col );
                info.StrV[Form("TEMPDS_%s_%s", (obj+cha+chg).c_str(), col.c_str())].push_back(dsTag);
              }
            }
          }
        } else {
          cout << "[ERROR] " << (obj+cha+chg) << " " << "MET" << " model for " << col << " was not found in the initial parameters!" << endl; return false;
        }
      }
    }
  }
  return true;
};


int importDataset(RooWorkspace& myws  , const std::map<string, RooWorkspace>& inputWS , const GlobalInfo& info)
{
  std::string cutDS = "";
  if (info.Par.at("Event_Type")!="") { cutDS += Form("(Event_Type==Event_Type::%s)&&", info.Par.at("Event_Type").c_str()); }
  for (auto var : info.Var) {
    if (var.first!="MET" && var.first!="Muon_Pt" && var.first!="Muon_Eta" && var.first!="Muon_Iso" && var.first!="Muon_MT") continue;
    if (var.second["Min"]==var.second["Max"]) { cutDS += Form("(%s == %.3f)", var.first.c_str(), var.second["Max"]); }
    else { cutDS += Form("(%.3f <= %s && %s < %.3f)", var.second["Min"], var.first.c_str(), var.first.c_str(), var.second["Max"]); }
    cutDS += "&&";
  }
  cutDS.erase(cutDS.size()-string("&&").length(), cutDS.size());
  TObjString tmp; tmp.SetString(cutDS.c_str()); myws.import(*((TObject*)&tmp), "Cut_DataSet"); // Save the cut expression for bookkeeping
  std::cout << "[INFO] Importing local RooDataSet with cuts: " << cutDS << std::endl;
  // Reduce and import the datasets
  for (auto& label : info.StrV.at("DSList")) {
    for (auto& chg : info.StrV.at("fitCharge")) {
      if ( inputWS.count(label)==0 || !(inputWS.at(label).data(Form("d%s_%s", chg.c_str(), label.c_str())))){ 
        std::cout << "[ERROR] The dataset " <<  Form("d%s_%s", chg.c_str(), label.c_str()) << " was not found!" << std::endl;
        return -1;
      }
      RooDataSet* data = (RooDataSet*)inputWS.at(label).data(Form("d%s_%s", chg.c_str(), label.c_str()))->reduce(cutDS.c_str());
      if (data->sumEntries()==0){ 
        std::cout << "[ERROR] No events from dataset " <<  Form("d%s_%s", chg.c_str(), label.c_str()) << " passed the kinematic cuts!" << std::endl; return -1;
      }
      else { myws.import(*data); }
      std::cout << "[INFO] " << data->numEntries() << " entries imported from local RooDataSet " << Form("d%s_%s", chg.c_str(), label.c_str()) << std::endl;
      delete data;
      myws.import(*((TObjString*)inputWS.at(label).obj("METType")), kTRUE);
      // Set the range of each global parameter in the local roodataset
      const RooArgSet* row = myws.data(Form("d%s_%s", chg.c_str(), label.c_str()))->get();
      for(auto& var : info.Var) {
        if ( (var.first!="Event_Type") && row->find(Form("%s", var.first.c_str())) ) {
          ((RooRealVar*)row->find(Form("%s", var.first.c_str())))->setMin(var.second.at("Min"));
          ((RooRealVar*)row->find(Form("%s", var.first.c_str())))->setMax(var.second.at("Max"));
        }
      }
    }
  }
  // Set the range of each global parameter in the local workspace
  for(auto& var : info.Var) {
    if ( myws.var(Form("%s", var.first.c_str())) ) {
      myws.var(Form("%s", var.first.c_str()))->setMin(var.second.at("Min"));
      myws.var(Form("%s", var.first.c_str()))->setMax(var.second.at("Max"));
    }
  }
  std::cout << "[INFO] Analyzing bin: " << Form("%.3f < Muon_Eta < %.3f", myws.var("Muon_Eta")->getMin(), myws.var("Muon_Eta")->getMax()) << std::endl;
  return 1;
};


void setMETGlobalParameterRange(RooWorkspace& myws, GlobalInfo& info)
{
  myws.var("MET")->setRange("METWindow", info.Var.at("MET").at("Min"), info.Var.at("MET").at("Max"));
  info.Par["METRange_Cut"] = Form("(MET>%.3f && MET<%.3f)", info.Var.at("MET").at("Min"), info.Var.at("MET").at("Max"));
  return;
};

void setMETFileName(string& fileName, string& outputDir, const string& DSTAG, const string& plotLabel, const GlobalInfo& info)
{
  std::string dsTag = DSTAG; dsTag.erase(dsTag.find("_MUON"), dsTag.length());
  std::string  metTAG = "MET" + info.Par.at("METType");
  outputDir = Form("%s%s/%s/", outputDir.c_str(), metTAG.c_str(), dsTag.c_str());
  fileName = Form("%s_%s_%s_MuEta_%.0f_%.0f_MuIso_%.0f_%.0f", "MET", 
                  dsTag.c_str(),
                  plotLabel.c_str(), 
                  (info.Var.at("Muon_Eta").at("Min")*10.0), (info.Var.at("Muon_Eta").at("Max")*10.0),
                  (info.Var.at("Muon_Iso").at("Min")*100.0), (info.Var.at("Muon_Iso").at("Max")*100.0)
                  );
  return;
};

void addLuminosity(GlobalInfo& info)
{
  if (info.Flag.at("fitMC")) {
    info.Par["Lumi_pPb"] = Form("%s[%.4f]", "Lumi_pPb", 1.0);
    info.Par["Lumi_Pbp"] = Form("%s[%.4f]", "Lumi_Pbp", 1.0);
    info.Par["Lumi_PA"] = Form("%s[%.4f]", "Lumi_PA", 1.0);
  }
  else {
    if (info.Par.count("Lumi_pPb")==0 || info.Par.at("Lumi_pPb")=="") { info.Par["Lumi_pPb"] = Form("%s[%.4f]", "Lumi_pPb", PA::Lumi_pPb); }
    if (info.Par.count("Lumi_Pbp")==0 || info.Par.at("Lumi_Pbp")=="") { info.Par["Lumi_Pbp"] = Form("%s[%.4f]", "Lumi_Pbp", PA::Lumi_Pbp); }
    if (info.Par.count("Lumi_PA")==0  || info.Par.at("Lumi_PA")=="" ) { info.Par["Lumi_PA"] = Form("%s[%.4f]", "Lumi_PA", (PA::Lumi_pPb+PA::Lumi_Pbp)); }
  }
}

#endif // #ifndef fitElectroWeakMETModel_C

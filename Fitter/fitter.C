#include "Macros/Utilities/initClasses.h"
#include "Macros/tree2DataSet.C"
#include "Macros/EWQ/fitElectroWeakMETModel.C"

bool checkSettings     ( const GlobalInfo& userInput);
bool parseFile         ( std::string FileName  , std::vector< StringMap >& data );
bool parseString       ( std::string input     , std::vector< double >& output );
bool iniWorkEnv        ( StringVectorMap& DIR  , const std::string& workDirName );
void findSubDir        ( std::vector< std::string >& dirlist, std::string dirname );
bool readFile          ( std::string FileName  , std::vector< std::vector< std::string > >& content, const int nCol=-1, int nRow=-1 );
bool getInputFileNames ( const std::string& InputTrees, std::map< std::string, std::vector< std::vector< std::string > > >& InputFileCollection );
bool setParameters     ( const StringMap& row , GlobalInfo& info , const std::string& Analysis );
bool addParameters     ( std::string InputFile , std::vector< GlobalInfo >& infoVector , const std::string& Analysis );

void fitter(
            const std::string workDirName = "QCDTemplate",// Working directory
            const std::bitset<1> useExt   = 0,          // Use external: (bit 0 (1)) Input DataSets
            // Select the type of datasets to fit
            const std::bitset<2> fitData  = 1,          // Fit Sample: (bit 0 (1)) Data , (bit 1 (2)) MC
            const std::bitset<3> fitColl  = 7,          // Fit System: (bit 0 (1)) pPb  , (bit 1 (2)) Pbp   , (bit 2 (4)) PA
            const std::bitset<3> fitChg   = 3,          // Fit Charge: (bit 0 (1)) Plus , (bit 1 (2)) Minus , (bit 2 (4)) Inclusive
            // Select the type of objects to fit
            const std::bitset<3> fitObj   = 2,          // Fit Objects: (bit 0 (1)) W , (bit 1 (2)) QCD
            // Select the fitting options
            const unsigned int   numCores = 32,         // Number of cores used for fitting
            const std::bitset<1> fitVar   = 1,          // Fit Variable: 1: MET
            const uint           varType  = 0,          // Type of MET to Fit: (0) PF Raw, (1) PF Type1, (2) NoHF Raw, (3) NoHF Type 1
            const std::string    Analysis = "WToMuNu",  // Type of Analysis
            // Select the drawing options
            const bool setLogScale  = true              // Draw plot with log scale
            )
{
  //  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  // -------------------------------------------------------------------------------
  // STEP 0: INITIALIZE THE FITTER WORK ENVIROMENT
  // The work enviroment is divided as follows:
  /*
    main |-> Macros: Contain all the macros
         |-> Input   |-> <WorkDir> : Contain Input File, Bin and Parameter List for a given work directory (e.g. 20160201)
	 |-> Output  |-> <WorkDir> : Contain Output Plots and Results for a given work directory (e.g. 20160201)
	 |-> DataSet : Contain all the datasets (MC and Data)
  */

  // Suppress Messages for RooFit
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Caching);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Plotting);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Integration);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::NumIntegration);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Minimization);
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;

  GlobalInfo userInput;
  // Store more information for fitting
  userInput.Par["extTreesFileDir"]    = "";
  userInput.Par["extDSDir_DATA"]      = "";
  userInput.Par["extDSDir_MC"]        = "";
  if (Analysis.find("Nu")!=std::string::npos) {
    userInput.Var["MET"]["type"]        = varType;
    userInput.Var["MET"]["binWidth"]    = 2.0;
    userInput.Par["extFitDir_MET"]      = "";
    userInput.Par["extInitFileDir_MET"] = "";
  }
  // Set all the Boolean Flags from the input settings
  userInput.StrV["setext"] = std::vector<std::string>({"ExtDS"});
  for (uint i=0; i<userInput.StrV.at("setext").size(); i++) { userInput.Flag["use"+userInput.StrV.at("setext")[i]] = useExt[i]; }
  userInput.StrV["sample"] = std::vector<std::string>({"Data", "MC"});
  for (uint i=0; i<userInput.StrV.at("sample").size(); i++) { userInput.Flag["fit"+userInput.StrV.at("sample")[i]] = fitData[i]; }
  userInput.StrV["system"] = std::vector<std::string>({"pPb", "Pbp", "PA"});
  for (uint i=0; i<userInput.StrV.at("system").size(); i++) { userInput.Flag["fit"+userInput.StrV.at("system")[i]] = fitColl[i]; }
  userInput.Flag["doPA"]  = (userInput.Flag.at("fitpPb") || userInput.Flag.at("fitPbp") || userInput.Flag.at("fitPA"));
  userInput.Flag["dopPb"] = (userInput.Flag.at("fitpPb") || userInput.Flag.at("fitPA"));
  userInput.Flag["doPbp"] = (userInput.Flag.at("fitPbp") || userInput.Flag.at("fitPA"));
  if (Analysis.find("Nu")!=std::string::npos) {
    userInput.StrV["variable"] = std::vector<std::string>({"MET"});
    userInput.StrV["METType"] = std::vector<std::string>({"PF_RAW","PF_Type1","PF_NoHF_RAW","PF_NoHF_Type1"});
    userInput.Par["METType"] = userInput.StrV.at("METType").at(varType);
  }
  if (Analysis.find("W")!=std::string::npos) {
    userInput.StrV["charge"] = std::vector<std::string>({"Plus","Minus","ChgInc"});
    userInput.StrV["object"] = std::vector<std::string>({"W", "QCD", "DY"});
    userInput.StrV["template"] = std::vector<std::string>({"W","QCD","DY","WToTau","TTbar"});
  }
  for (uint i=0; i<userInput.StrV.at("charge").size(); i++) { userInput.Flag["fit"+userInput.StrV.at("charge")[i]] = fitChg[i]; }
  for (const auto& chg : userInput.StrV.at("charge") ) {
    std::string label = ""; if (chg=="ChgInc") { label = ""; } else { label = chg.substr(0,2); }
    if (userInput.Flag.at("fit"+chg)) { userInput.StrV["fitCharge"].push_back(label); }
  }
  for (uint i=0; i<userInput.StrV.at("object").size(); i++) { userInput.Flag["fit"+userInput.StrV.at("object")[i]] = fitObj[i]; }
  for (const auto& obj : userInput.StrV.at("object") ) { if (userInput.Flag.at("fit"+obj)) { userInput.StrV["fitObject"].push_back(obj); } }
  for (uint i=0; i<userInput.StrV.at("variable").size(); i++) { userInput.Flag["fit"+userInput.StrV.at("variable")[i]] = fitVar[i]; }
  for (const auto& var : userInput.StrV.at("variable") ) { if (userInput.Flag.at("fit"+var)) { userInput.StrV["fitVariable"].push_back(var); } }
  for (const auto& obj : userInput.StrV.at("object") ) { userInput.Flag["incMCTemp_"+obj]   = false; } // Value set in STEP 1 (loading initial parameters)
  userInput.Flag["setLogScale"] = setLogScale;
  // Set all the Parameters from the input settings
  userInput.Par["Analysis"] = Analysis;
  userInput.Par["anaType"] = "", userInput.Par["chgLabel"] = "";
  if (userInput.Par.at("Analysis").find("W")!=std::string::npos) { userInput.Par.at("anaType") = "EWQ"; }
  userInput.Flag["doMuon"] = (Analysis.find("Mu")!=std::string::npos);
  if (userInput.Flag.at("doMuon")) { userInput.Par["channel"] = "ToMu"; userInput.Par["channelDS"] = "MUON"; }
  userInput.Flag["doElec"] = (Analysis.find("Ele")!=std::string::npos);
  if (userInput.Flag.at("doElec")) { userInput.Par["channel"] = "ToEl"; userInput.Par["channelDS"] = "ELEC"; }
  userInput.Int["numCores"] = numCores;
  bool fitTest = (workDirName=="Test");
  for (const auto& variab : userInput.StrV.at("variable")) { if (userInput.Flag.at("fit"+variab) || fitTest) { userInput.Par["extFitDir_"+variab] = ""; userInput.Par["extInitFileDir_"+variab] = "";} }
  // Check the User Input Settings
  if (!checkSettings(userInput)){ return; }
  
  // Set the Local Work Enviroment
  StringVectorMap DIR;
  if(!iniWorkEnv(DIR, workDirName)){ return; }
  /////////////////////
  std::map< std::string, std::string> inputFitDir;
  inputFitDir["MET"] = userInput.Par.at("extFitDir_MET");
  std::map< std::string, std::string> inputInitialFilesDir;
  inputInitialFilesDir["MET"] = userInput.Par.at("extInitFileDir_MET");
  inputInitialFilesDir["FILES"] = userInput.Par.at("extTreesFileDir");
  // Initiliaze all the input Fit and Initial File Directories
  std::vector< std::map< std::string, std::string> > inputFitDirs;
  inputFitDirs.push_back(inputFitDir);
  for(uint i=1; i<DIR["input"].size(); i++) {
    inputFitDirs.push_back(inputFitDir);
    for (const auto& iter : inputFitDirs[i]) {
      std::string key = iter.first;
      if (inputFitDirs[i][key]!="") {
        inputFitDirs[i][key] = DIR["input"][i];
        inputFitDirs[i][key].replace(inputFitDirs[i][key].find(DIR["input"][0]), DIR["input"][0].length(), inputFitDirs[0][key]);
      }
    }
  }
  std::vector< std::map< std::string, std::string> > inputInitialFilesDirs;
  inputInitialFilesDirs.push_back(inputInitialFilesDir);
  for(uint i=1; i<DIR["input"].size(); i++) {
    inputInitialFilesDirs.push_back(inputInitialFilesDir);
    for (const auto& iter : inputInitialFilesDirs[i]) {
      std::string key = iter.first;
      if (inputInitialFilesDirs[i][key]!="") {
        inputInitialFilesDirs[i][key] = DIR["input"][i];
        inputInitialFilesDirs[i][key].replace(inputInitialFilesDirs[i][key].find(DIR["input"][0]), DIR["input"][0].length(), inputInitialFilesDirs[0][key]);
      }
    }
  }
  if (userInput.Par.at("extTreesFileDir")=="") userInput.Par.at("extTreesFileDir") = DIR["input"][0];

  // -------------------------------------------------------------------------------
  // STEP 1: LOAD THE INITIAL PARAMETERS
  /*
    Input : List of initial parameters with format PT <tab> RAP <tab> CEN <tab> iniPar ... 
    Output: two vectors with one entry per kinematic bin filled with the cuts and initial parameters
  */

  std::vector< std::vector< GlobalInfo > > infoVectors; 
  std::map< std::string, std::map< std::string, bool > > VARMAP;
  for (const auto& var : userInput.StrV.at("variable")) {
    for (const auto& obj : userInput.StrV.at("object")) {
      VARMAP[var][obj] = (userInput.Flag.at("fit"+var) && userInput.Flag.at("fit"+obj));
    }
    if (VARMAP.at(var).at("W")) { VARMAP.at(var).at("QCD") = true; }
  }
  std::map< std::string , bool > COLMAP;
  for (const auto& col : userInput.StrV.at("system")) {
    COLMAP[col] = userInput.Flag.at("fit"+col);
  }
  for(uint j = 0; j < DIR["input"].size(); j++) {
    if (DIR["input"].size()>1 && j==0) continue; // First entry is always the main input directory
    std::vector< GlobalInfo > infoVector;
    for (const auto& VAR : VARMAP) {
      std::map< std::string , bool > PARMAP = VAR.second;
      for (const auto& PAR : PARMAP) {
        if (PAR.second) {
          for (const auto& COL : COLMAP) {
            if(COL.second) {
              std::string dir = DIR["input"][j];
              if (inputInitialFilesDirs[j][VAR.first]!="") { dir = inputInitialFilesDirs[j][VAR.first]; }
              std::string InputFile = "", name = (dir + "InitialParam_" + VAR.first + "_" + PAR.first);
              std::vector<std::string> tryChannel = { userInput.Par.at("channel") , "" };
              std::vector<std::string> trySystem  = { COL.first, "PA" };
              for (const auto& tryCha : tryChannel) {
                bool trySuccess = false;
                for (const auto& tryCol : trySystem) {
                  if (ifstream(InputFile).good()==false) { InputFile = (name + tryCha + "_" + tryCol + ".csv"); } else { trySuccess = true; break; }
                }
                if (trySuccess) break;
              }
              if (!addParameters(InputFile, infoVector, userInput.Par.at("Analysis"))) { return; }
              if (!userInput.Flag.at("incMCTemp_"+PAR.first)) { for (const auto& info : infoVector) { if(info.Flag.at("incMCTemp_"+PAR.first)) { userInput.Flag.at("incMCTemp_"+PAR.first) = true; break; } } }
            }
          }
        }
      }
    }
    infoVectors.push_back(infoVector);
  }

  // -------------------------------------------------------------------------------
  // STEP 2: CREATE/LOAD THE ROODATASETS
  /*
    Input : List of TTrees with format:  TAG <tab> FILE_NAME
    Output: Collection of RooDataSets splitted by tag name
  */

  const std::string InputTrees = userInput.Par.at("extTreesFileDir") + "InputTrees.txt";
  std::map< std::string, std::vector< std::vector< std::string > > > InputFileCollection;
  if(!getInputFileNames(InputTrees, InputFileCollection)){ return; }

  TObjArray* aDSTAG = new TObjArray(); // Array to store the different tags in the list of trees
  aDSTAG->SetOwner(true);
  
  std::map< std::string, RooWorkspace > Workspace;
  for (const auto& FileCollection : InputFileCollection) {
    // Get the file tag which has the following format: DSTAG_CHAN_COLL , i.e. DATA_MUON_Pbp
    string FILETAG = FileCollection.first;
    if (!FILETAG.size()) { std::cout << "[ERROR] FILETAG is empty!" << std::endl; return; }
    userInput.Par["localDSDir"] = DIR["dataset"][0];
    // Extract the filenames
    std::string dir = "";
    if ( (FILETAG.find("MUON")!=std::string::npos) &&  !userInput.Flag.at("doMuon") ) continue; // If we find Muon, check if the user wants Muon channel
    if ( (FILETAG.find("ELEC")!=std::string::npos) &&  !userInput.Flag.at("doElec") ) continue; // If we find Electron, check if the user wants Electron channel
    if ( (FILETAG.find("PA")!=std::string::npos)   &&  !userInput.Flag.at("doPA")   ) continue; // If we find PA, check if the user wants to do PA
    if ( (FILETAG.find("pPb")!=std::string::npos)  &&  !userInput.Flag.at("dopPb")  ) continue; // If we find pPb, check if the user wants pPb
    if ( (FILETAG.find("Pbp")!=std::string::npos)  &&  !userInput.Flag.at("doPbp")  ) continue; // If we find Pbp, check if the user wants Pbp
    bool fitDS = false;
    // If we have data, check if the user wants to fit data
    if ( (FILETAG.find("DATA")!=std::string::npos) ) {
      if (userInput.Flag.at("fitData")) {
        dir = userInput.Par.at("localDSDir");
        if (userInput.Flag.at("useExtDS") && userInput.Par.at("extDSDir_DATA")!="" && (existDir(userInput.Par.at("extDSDir_DATA"))==true)) { dir = userInput.Par.at("extDSDir_DATA"); }
        fitDS = true;
      }
    }
    // If we find MC, check if the user wants to fit MC
    if ( (FILETAG.find("MC")!=std::string::npos) ) {
      bool keep = false;
      if (userInput.Flag.at("fitMC")) { for (const auto& obj : userInput.StrV.at("object")) { if ( (FILETAG.find(obj)!=std::string::npos) && userInput.Flag.at("fit"+obj)) { keep = true; fitDS = true; break; } } }
      bool checkTemp = false; for(const auto& s : std::vector<std::string>({"W", "QCD", "DY"})) { if (userInput.Flag.at("incMCTemp_"+s)) { checkTemp = true; break; } }
      if (!keep && checkTemp) { for (const auto& tmp : userInput.StrV.at("template")) { if ( (FILETAG.find(tmp)!=std::string::npos) ) { keep = true; fitDS = false; break; } } }
      if (keep) {
        dir = userInput.Par.at("localDSDir");
        if (userInput.Flag.at("useExtDS") && userInput.Par.at("extDSDir_MC")!="" && (existDir(userInput.Par.at("extDSDir_MC"))==true)) { dir = userInput.Par.at("extDSDir_MC"); }
      }
    }
    if (dir!="") {
      StringVectorMap FileInfo;
      FileInfo["InputFileNames"].clear(); FileInfo["XSectionTags"].clear(); 
      for (const auto& row : FileCollection.second) {
        FileInfo.at("InputFileNames").push_back(row[0]);
        if (row.size()>1) { FileInfo.at("XSectionTags").push_back(row[1]); }
      }
      FileInfo["OutputFileDir"].push_back(dir);
      FileInfo["OutputFileDir"].push_back(DIR["dataset"][0]);
      if (FILETAG.find("PA")==std::string::npos) { FileInfo["DSNames"].push_back(FILETAG); }
      else {
        std::string NAMETAG = FILETAG; NAMETAG.erase(NAMETAG.find_last_of("_"), 5);
        if (userInput.Flag.at("dopPb")) FileInfo["DSNames"].push_back(NAMETAG+"_pPb");
        if (userInput.Flag.at("doPbp")) FileInfo["DSNames"].push_back(NAMETAG+"_Pbp");
      }
      // Produce the output datasets
      if(!tree2DataSet(Workspace, FileInfo, userInput.Par.at("METType"), userInput.Par.at("Analysis"))){ return; }
      if (fitDS) { for (const auto& DSTAG : FileInfo.at("DSNames")) { if (!aDSTAG->FindObject(DSTAG.c_str())) aDSTAG->Add(new TObjString(DSTAG.c_str())); } }
    }
  }
  if (Workspace.size()==0) {
    std::cout << "[ERROR] No tree files were found matching the user's input settings!" << std::endl; return;
  }

  // -------------------------------------------------------------------------------
  // STEP 3: COMBINE THE ROODATASETS
  /*
    Input : Collection of RooWorkspaces containing the pPb and Pbp RooDataSets
    Output: Collection of RooWorkspaces including the combined PA RooDataSets
  */

  if (userInput.Flag.at("fitPA")) {
    // Determine the sample tags (sample name without the collision tag)
    std::vector< std::string > sampleTags;
    for (const auto& w : Workspace) {
      std::string NAMETAG = w.first; NAMETAG.erase(NAMETAG.find_last_of("_"), 5);
      if (std::find(sampleTags.begin(), sampleTags.end(), NAMETAG)==sampleTags.end()) { sampleTags.push_back(NAMETAG); }
    }
    // Loop over each sample tag
    for (const auto& sample : sampleTags) {
      if (!createPADataset(Workspace, sample)) { return; }
      if (aDSTAG->FindObject((sample+"_pPb").c_str()) && !aDSTAG->FindObject((sample+"_PA").c_str())) aDSTAG->Add(new TObjString((sample+"_PA").c_str()));
    }
  }

  // -------------------------------------------------------------------------------  
  // STEP 4: FIT THE DATASETS
  /*
    Input : 
              -> The cuts and initial parameters per kinematic bin
	      -> The workspace with the full datasets included.
    Output: 
              -> Plots (png, pdf and root format) of each fit.
	      -> The local workspace used for each fit.
  */

  TIter nextDSTAG(aDSTAG);
  for(uint j = 0; j < infoVectors.size(); j++) {
    const int index = ( DIR["output"].size()>1 ? j+1 : j ); // First entry is always the main output directory
    const std::string outputDir = DIR["output"][index];
    //
    for (unsigned int i=0; i < infoVectors[j].size(); i++) {
      nextDSTAG.Reset();
      TObjString* soDSTAG(0x0);
      while ( (soDSTAG = static_cast<TObjString*>(nextDSTAG.Next())) )
        {
          const TString DSTAG = (TString)(soDSTAG->GetString());
          //
          if (Workspace.count(DSTAG.Data())>0) {
            //
            if (
                ( userInput.Flag.at("fitpPb") && DSTAG.Contains("pPb") ) ||
                ( userInput.Flag.at("fitPbp") && DSTAG.Contains("Pbp") ) ||
                ( userInput.Flag.at("fitPA" ) && DSTAG.Contains("PA" ) )
                )
              {
                //
                if (!fitElectroWeakMETModel( Workspace, infoVectors[j].at(i),
                                             userInput, 
                                             // Select the type of datasets to fit
                                             outputDir,
                                             DSTAG.Data()
                                             )
                    ) { return; }
              }
          } else {
            std::cout << "[ERROR] The workspace for " << DSTAG.Data() << " was not found!" << std::endl; return;
          }
        }
    }
  }
  aDSTAG->Delete();
  delete aDSTAG;
};


bool addParameters(std::string InputFile, std::vector< GlobalInfo >& infoVector, const std::string& Analysis)
{
  std::vector< StringMap >  data;
  if(!parseFile(InputFile, data)) { return false; }
  if (infoVector.size()==0) {
    for (const auto& row : data) {
      GlobalInfo info = GlobalInfo();
      if(!setParameters(row, info, Analysis)) { return false; }
      infoVector.push_back(info);
    }
  }
  else {
    if (data.size()!=infoVector.size()) { std::cout << "[ERROR] The initial parameters in file " << InputFile << " are not consistent with previous files!" << std::endl; return false; }
    for (unsigned int i=0; i<data.size(); i++) {
      GlobalInfo info = GlobalInfo();
      if (!setParameters(data.at(i), info, Analysis)) { return false; };
      if (info.Var != infoVector.at(i).Var) { std::cout << "[ERROR] The bins in file " << InputFile << " are not consistent with previous files!" << std::endl; return false; }
      infoVector.at(i).Copy(info, true);
    }
  }
  return true;
};


bool setParameters(const StringMap& row, GlobalInfo& info, const std::string& Analysis)
{
  // set initial values of variables
  if (Analysis.find("WToMuNu")!=std::string::npos) {
    info.Var["MET"]["Min"]        = 0.0;
    info.Var["MET"]["Max"]        = 100000.0;
    info.Var["Muon_Pt"]["Min"]    = 0.0;
    info.Var["Muon_Pt"]["Max"]    = 100000.0;
    info.Var["Muon_Eta"]["Min"]   = -2.5;
    info.Var["Muon_Eta"]["Max"]   = 2.5;
    info.Var["Muon_Iso"]["Min"]   = 0.0;
    info.Var["Muon_Iso"]["Max"]   = 100000.0;
    info.Var["Muon_MT"]["Min"]    = 0.0;
    info.Var["Muon_MT"]["Max"]    = 100000.0;
    info.Var["Centrality"]["Min"] = 0.0;
    info.Var["Centrality"]["Max"] = 100000.0;
    info.Par["Event_Type"] = "";
    info.Par["Muon_Chg"] = "";
    info.Par["Run_Type"] = "";
  }
  info.Par["Model"] = "";
  info.Par["Cut"] = "";
  for(const auto& s : std::vector<std::string>({"W", "QCD", "DY"})) { info.Flag["incMCTemp_"+s]   = false; }
  info.Flag["useEtaCM"] = (row.count("Muon_EtaCM") > 0);
  // set parameters from file
  for (const auto& col : row) {
    if (info.Flag.at("useEtaCM") && col.first=="Muon_Eta") continue;
    std::string colName = ( (col.first=="Muon_EtaCM") ? "Muon_Eta" : col.first );
    bool found = false;
    for (const auto& var : info.Var) {
      const std::string varName = var.first;
      if (colName==varName) {
        if (col.second=="") {
          std::cout << "[ERROR] Input column " << varName << " has invalid value: " << col.second << std::endl; return false;
        }
        std::vector<double> v;
        if (!parseString(col.second, v)) { return false; }
        if (v.size()!=2 && v.size()!=1) {
          std::cout << "[ERROR] Input column " << varName << " has incorrect number of values, it should have 1 or 2 values but has: " << v.size() << std::endl; return false;
        }
        info.Var.at(varName).at("Min") = v.at(0);
        info.Var.at(varName).at("Max") = v.at(v.size()-1);
        found = true; break;
      }
    }
    if (found==false) {
      for (const auto& par : info.Par) {
        std::string parName = par.first;
        if (colName.find(par.first)!=std::string::npos) {
          if (col.second=="") {
            std::cout << "[ERROR] Input column " << par.first << " has empty value" << std::endl; return false;
          }
          info.Par[colName] = col.second;
          for(const auto& s : std::vector<std::string>({"W", "QCD", "DY"})) {
            if ( (par.first=="Model") && (colName.find(s)!=std::string::npos)   && (col.second.find("TEMP")!=std::string::npos) ) { info.Flag.at("incMCTemp_"+s)   = true; }
          }
          found = true;
        }
      }
    }
    if (found==false) {
      if (col.second != "") {
        std::string value = col.second;
	// check that initial parameters format is correct: [ num, num, num ]
	if ((value.find("[")==std::string::npos)||(value.find("]")==std::string::npos)) {
	  // Special cases like parameter constrains could be set here but for now, let's keep it simple
          std::cout << "[ERROR] Either ']' or '[' are missing in the initial parameter values!" << std::endl; return false;
	} else {
	  value.erase(0, value.find("[")+std::string("[").length());
	  value.erase(value.find("]"), value.length());
	}
	std::vector<double> v;
	if (!parseString(value, v)) { return false; }
        if (v.size()>3 || v.size()<1) {
          std::cout << "[ERROR] Initial parameter " << col.first << " has incorrect number of values, it has: " << v.size() << std::endl; return false;
        }
	// everything seems alright, then proceed to save the values
	if (v.size()==1) {
	  // if only one value is given i.e. [ num ], consider it a constant value
	  info.Par[col.first] = Form("%s[%.6f]", col.first.c_str(), v.at(0));
	} else if (v.size()==2) {
          if (col.second.find("RNG")!=std::string::npos) {
            info.Par[col.first] = Form("%s[%.6f, %.6f]", col.first.c_str(), v.at(0), v.at(1));
          } 
          else { // For Constrained Fits
            info.Par[col.first] = Form("%s[%.6f, %.6f, %.6f]", col.first.c_str(), v.at(0), (v.at(0)-20.*v.at(1)), (v.at(0)+20.*v.at(1)));
            info.Par["val"+col.first] = Form("%s[%.6f]", ("val"+col.first).c_str(), v.at(0));
            info.Par["sig"+col.first] = Form("%s[%.6f]", ("sig"+col.first).c_str(), v.at(1));
          }
	} else if (v.size()==3) {
	  info.Par[col.first] = Form("%s[%.6f, %.6f, %.6f]", col.first.c_str(), v.at(0), v.at(1), v.at(2));
	}
      } else {
        info.Par[col.first] = "";
      }
    }
  }
  return true;
};


bool parseString(std::string input, std::vector< double >& output)
{
  // remove spaces from input string 
  input.erase(std::remove(input.begin(), input.end(), ' '), input.end());
  // proceed to parse input string
  char *end;
  while(input!="") {
    double d = strtod(input.c_str(), &end);
    if (end != input) {
      output.push_back(d);
    } else {
      cout << "[ERROR] The conversion from string to double failed!"; return false;
    }
    input = end;
    if (input!="") { input.erase(input.begin()); } // Delete the delimiter, should have length = 1
  }
  return true;
};


bool parseFile(std::string FileName, std::vector< StringMap >& data)
{
  std::vector< std::vector< std::string > > content, tmp; 
  if(!readFile(FileName, tmp, -1, 1)){ return false; }
  std::vector< std::string > header = tmp.at(0);
  if (header.size()==0) { std::cout << "[ERROR] The header is null!" << std::endl; return false; }
  if(!readFile(FileName, content, header.size())){ return false; }
  for (const auto& rHeader : header) {
    if (rHeader=="") { std::cout << "[ERROR] A column has no label!" << std::endl; return false; }
  }
  content.erase(content.begin()); // remove header
  for (const auto& row : content) {
    StringMap col;
    for (unsigned int i=0; i<header.size(); i++) {
      if (i<row.size()) {
	col[header.at(i)] = row.at(i);
      } else {
	col[header.at(i)] = "";
      }
    }
    data.push_back(col);
  }
  return true;
};					


bool getInputFileNames(const std::string& InputTrees, std::map< std::string, std::vector< std::vector< std::string > > >& InputFileCollection)
{
  std::vector< std::vector< std::string > > content; 
  if(!readFile(InputTrees, content, 3)){ return false; }
  for(auto& row : content) {
    for (auto& col : row) {
      col.erase(std::remove(col.begin(), col.end(), ' '), col.end());  // remove spaces
      col.erase(std::remove(col.begin(), col.end(), '\t'), col.end()); // remove tabs
    }
    if (row.at(0)!="" && row.at(1)=="") { std::cout << "[ERROR] There is an empty file name in your InputTrees.txt, please fix it" << std::endl; return false; }
    if (row.at(0)=="" && row.at(1)!="") { std::cout << "[ERROR] There is an empty file tag in your InputTrees.txt, please fix it" << std::endl; return false; }
    if (row.at(0)!="" && row.at(1)!="") {
      // store the filenames mapped by the tag
      std::vector< std::string > tmp;
      for (uint i = 1; i < row.size(); i++) {
        tmp.push_back(row.at(i));
      }
      InputFileCollection[row.at(0)].push_back(tmp);
    }
  }
  return true;
};


bool readFile(std::string FileName, std::vector< std::vector< std::string > >& content, const int nCol, int nRow)
{
  if (nCol==0 || nRow==0) { 
    std::cout << "[WARNING] Ignoring content of File: " << FileName << std::endl; return true; 
  }
  if (nRow!=1) { std::cout << "[INFO] Reading file: " << FileName << std::endl; }
  ifstream myfile(FileName.c_str());
  char delimiter = ' ';
  if (myfile.is_open()){ 
    std::string line, CHAR;
    while ( getline(myfile, line) ){
      std::stringstream row(line), tmp(line); tmp >> CHAR;
      if ((!tmp) || (CHAR.find('#')!=std::string::npos) || (CHAR.find("//")!=std::string::npos)) continue;
      if (delimiter == ' ' && line.find(',')!=std::string::npos) { delimiter = ','; }
      if (delimiter == ' ' && line.find(';')!=std::string::npos) { delimiter = ';'; }
      if (delimiter == ' ') { std::cout << "[ERROR] File: " << FileName << " has unknown delimiter!" << std::endl; return false; }
      if (nRow==0) break; else {nRow=nRow-1;}
      std::vector< std::string > cols; int i=0;
      while (true){
        std::string col; getline(row, col, delimiter);
	if ( (nCol>=0) ? (i>=nCol) : (col=="") ){ break; }
	cols.push_back(col);
	i++;
      }
      content.push_back(cols);
    }
  } else {
    std::cout << "[ERROR] File: " << FileName << " was not found!" << std::endl; return false;
  }
  return true;
};


bool iniWorkEnv(std::map< std::string, std::vector< std::string > >& DIR, const std::string& workDirName)
{
  std::cout << "[INFO] Initializing the work enviroment" << std::endl;
  DIR["main"].push_back(gSystem->ExpandPathName(gSystem->pwd()));
  DIR["macros"].push_back(DIR["main"][0] + "/Macros/");
  if (existDir(DIR["macros"][0].c_str())==false){ 
    std::cout << "[ERROR] Input directory: " << DIR["macros"][0] << " doesn't exist!" << std::endl;
    return false; 
  }
  DIR["input"].push_back(DIR["main"][0] + "/Input/" + workDirName + "/");
  if (existDir(DIR["input"][0])==false){ 
    std::cout << "[ERROR] Input directory: " << DIR["input"][0] << " doesn't exist!" << std::endl;
    return false; 
  } else {
    findSubDir(DIR["input"], DIR["input"][0]);
  }
  DIR["output"].push_back(DIR["main"][0] + "/Output/" + workDirName + "/");
  makeDir(DIR["output"][0]);
  for(uint j = 1; j < DIR["input"].size(); j++) {
    std::string subdir = DIR["input"][j];
    subdir.replace(subdir.find("/Input/"), std::string("/Input/").length(), "/Output/");
    makeDir(subdir);
    DIR["output"].push_back(subdir);
  } 
  DIR["dataset"].push_back(DIR["main"][0] + "/DataSet/");
  makeDir(DIR["dataset"][0]);
  return true;
};


void findSubDir(std::vector< std::string >& dirlist, std::string dirname)
{
  TSystemDirectory dir(dirname.c_str(), dirname.c_str());
  TList *subdirs = dir.GetListOfFiles();
  if (subdirs) {
    TSystemFile *subdir;
    TIter next(subdirs);
    while ((subdir=(TSystemFile*)next())) {
      if (subdir->IsDirectory() && std::string(subdir->GetName())!="." && std::string(subdir->GetName())!="..") {
        dirlist.push_back(dirname + subdir->GetName() + "/");
        std::cout << "[INFO] Input subdirectory: " << dirname + subdir->GetName() + "/" << " found!" << std::endl;
      }
    }
  }
  delete subdirs;
  return;
};


bool checkSettings(const GlobalInfo& userInput)
{ 
  std::cout << "[INFO] Checking user settings " << std::endl;

  if (userInput.Flag.at("fitW") && userInput.Flag.at("fitQCD")) { std::cout << "[ERROR] We can not fit QCD and W at the same time!" << std::endl; return false; }

  std::cout << "[INFO] All user setting are correct " << std::endl;
  return true;
};

#include "Macros/Utilities/initClasses.h"
#include "Macros/tree2DataSet.C"
#include "Macros/EWQ/fitElectroWeakMETModel.C"

bool checkSettings     ( const GlobalInfo& userInput);
bool parseFile         ( std::string FileName  , std::vector< StringMap >& data );
bool parseString       ( std::string input     , std::vector< double >& output );
bool iniWorkEnv        ( StringVectorMap& DIR  , const std::string& workDirName );
void findSubDir        ( std::vector< std::string >& dirlist, std::string dirname );
bool existDir          ( std::string dir );
bool readFile          ( std::string FileName  , std::vector< std::vector< std::string > >& content, const int nCol=-1, int nRow=-1 );
bool getInputFileNames ( const std::string& InputTrees, std::map< std::string, std::vector< std::string > >& InputFileCollection );
bool setParameters     ( const StringMap& row , GlobalInfo& info , const std::string& Analysis );
bool addParameters     ( std::string InputFile , std::vector< GlobalInfo >& infoVector , const std::string& Analysis );

void fitter(
            const std::string workDirName = "Test_TEMPMC",     // Working directory
            const std::bitset<1> useExt   = 0,          // Use external: (bit 0) Input DataSets
            // Select the type of datasets to fit
            const std::bitset<2> fitData  = 1,          // Fit Sample: (bit 0) Data , (bit 1) MC
            const std::bitset<3> fitColl  = 3,          // Fit System: (bit 0) pPb  , (bit 1) Pbp   , (bit 2) PA
            const std::bitset<3> fitChg   = 3,          // Fit Charge: (bit 0) Plus , (bit 1) Minus , (bit 2) Inclusive
            // Select the type of objects to fit
            const std::bitset<3> fitObj   = 1,          // Fit Objects: (bit 0) W , (bit 1) QCD , (bit 2) DYZ
            // Select the fitting options
            const unsigned int   numCores = 32,         // Number of cores used for fitting
            const std::bitset<1> fitVar   = 1,          // Fit Variable: 1: MET
            const uint           varType  = 1,          // Type of MET to Fit: (0) PF Raw, (1) PF Type1, (2) NoHF Raw, (3) NoHF Type 1
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
  if (userInput.Flag.at("fitPA")) { userInput.Flag.at("fitpPb") = true; userInput.Flag.at("fitPbp") = true; }
  userInput.Flag["doPA"] = (userInput.Flag.at("fitpPb") || userInput.Flag.at("fitPbp"));
  if (Analysis.find("Nu")!=std::string::npos) { 
    userInput.StrV["variable"] = std::vector<std::string>({"MET"});
    userInput.StrV["METType"] = std::vector<std::string>({"PFRaw","PFType1","NoHFRaw","NoHFType1"});
    userInput.Par["METType"] = userInput.StrV.at("METType").at(varType);
  }
  if (Analysis.find("W")!=std::string::npos) {
    userInput.StrV["charge"] = std::vector<std::string>({"Plus","Minus","ChgInc"});
    userInput.StrV["object"] = std::vector<std::string>({"W", "QCD", "DYZ"});
    userInput.StrV["template"] = std::vector<std::string>({"W","QCD","DYZ","WToTau"});
  }
  for (uint i=0; i<userInput.StrV.at("charge").size(); i++) { userInput.Flag["fit"+userInput.StrV.at("charge")[i]] = fitChg[i]; }
  for (auto& chg : userInput.StrV.at("charge") ) {
    std::string label = ""; if (chg=="ChgInc") { label = ""; } else { label = chg.substr(0,2); }
    if (userInput.Flag.at("fit"+chg)) { userInput.StrV["fitCharge"].push_back(label); }
  }
  for (uint i=0; i<userInput.StrV.at("object").size(); i++) { userInput.Flag["fit"+userInput.StrV.at("object")[i]] = fitObj[i]; }
  for (auto& obj : userInput.StrV.at("object") ) { if (userInput.Flag.at("fit"+obj)) { userInput.StrV["fitObject"].push_back(obj); } }
  for (uint i=0; i<userInput.StrV.at("variable").size(); i++) { userInput.Flag["fit"+userInput.StrV.at("variable")[i]] = fitVar[i]; }
  for (auto& var : userInput.StrV.at("variable") ) { if (userInput.Flag.at("fit"+var)) { userInput.StrV["fitVariable"].push_back(var); } }
  userInput.Flag["incMCTemp"]   = false; // Value set in STEP 1 (loading initial parameters)
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
  for (auto& variab : userInput.StrV.at("variable")) { if (userInput.Flag.at("fit"+variab) || fitTest) { userInput.Par["extFitDir_"+variab] = ""; userInput.Par["extInitFileDir_"+variab] = "";} }
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
    for(map<string, string>::iterator iter=inputFitDirs[i].begin(); iter!=inputFitDirs[i].end(); iter++) {
      string key = iter->first;
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
    for(map<string, string>::iterator iter=inputInitialFilesDirs[i].begin(); iter!=inputInitialFilesDirs[i].end(); iter++) {
      string key = iter->first;
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
  for (auto& var : userInput.StrV.at("variable")) {
    for (auto& obj : userInput.StrV.at("object")) {
      VARMAP[var][obj] = (userInput.Flag.at("fit"+var) && userInput.Flag.at("fit"+obj));
    }
  }
  std::map< std::string , bool > COLMAP;
  for (auto& col : userInput.StrV.at("system")) {
    COLMAP[col] = userInput.Flag.at("fit"+col); 
  }
  for(uint j = 0; j < DIR["input"].size(); j++) {
    if (DIR["input"].size()>1 && j==0) continue; // First entry is always the main input directory
    std::vector< GlobalInfo > infoVector;
    for(auto& VAR : VARMAP) {
      std::map< std::string , bool > PARMAP = VAR.second;
      for(auto& PAR : PARMAP) {
        if (PAR.second) {
          for(auto& COL : COLMAP) {
            if (userInput.Flag.at("fitPA") && COL.first!="PA") continue;
            if(COL.second) {
              std::string dir = DIR["input"][j];
              if (inputInitialFilesDirs[j][VAR.first]!="") { dir = inputInitialFilesDirs[j][VAR.first]; }
              std::string InputFile = "", name = (dir + "InitialParam_" + VAR.first + "_" + PAR.first);
              std::vector<std::string> tryChannel = { userInput.Par.at("channel") , "" };
              std::vector<std::string> trySystem  = ( userInput.Flag.at("doPA") ? std::vector<std::string>({COL.first , "PA"}) : std::vector<std::string>({COL.first}) );
              for (auto& tryCha : tryChannel) {
                bool trySuccess = false;
                for (auto& tryCol : trySystem) {
                  if (ifstream(InputFile).good()==false) { InputFile = (name + tryCha + "_" + tryCol + ".csv"); } else { trySuccess = true; break; }
                }
                if (trySuccess) break;
              }
              if (!addParameters(InputFile, infoVector, userInput.Par.at("Analysis"))) { return; }
              if (!userInput.Flag.at("incMCTemp")) { for (auto& info : infoVector) { if(info.Flag.at("incMCTemp")) { userInput.Flag.at("incMCTemp") = true; break; } } }
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

  const string InputTrees = userInput.Par.at("extTreesFileDir") + "InputTrees.txt";
  map<string, vector<string> > InputFileCollection;
  if(!getInputFileNames(InputTrees, InputFileCollection)){ return; }

  TObjArray* aDSTAG = new TObjArray(); // Array to store the different tags in the list of trees
  aDSTAG->SetOwner(true);
  std::map<string, RooWorkspace> Workspace;
  
  for(auto& FileCollection : InputFileCollection) {
    // Get the file tag which has the following format: DSTAG_CHAN_COLL , i.e. DATA_MUON_Pbp
    string FILETAG = FileCollection.first;
    if (!FILETAG.size()) { cout << "[ERROR] FILETAG is empty!" << endl; return; }
    userInput.Par["localDSDir"] = DIR["dataset"][0];
    // Extract the filenames
    StringVector InputFileNames = FileCollection.second;
    std::string dir = "";
    if ( (FILETAG.find("MUON")!=std::string::npos) &&  !userInput.Flag.at("doMuon") ) continue; // If we find Muon, check if the user wants Muon channel
    if ( (FILETAG.find("ELEC")!=std::string::npos) &&  !userInput.Flag.at("doElec") ) continue; // If we find Electron, check if the user wants Electron channel
    if ( (FILETAG.find("PA")!=std::string::npos)   &&  !userInput.Flag.at("doPA")   ) continue; // If we find PA, check if the user wants to do PA
    if ( (FILETAG.find("pPb")!=std::string::npos)  &&  !userInput.Flag.at("fitpPb") ) continue; // If we find pPb, check if the user wants pPb
    if ( (FILETAG.find("Pbp")!=std::string::npos)  &&  !userInput.Flag.at("fitPbp") ) continue; // If we find Pbp, check if the user wants Pbp
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
      if (userInput.Flag.at("fitMC")) { for (auto& obj : userInput.StrV.at("object")) { if ( (FILETAG.find(obj)!=std::string::npos) && userInput.Flag.at("fit"+obj)) { keep = true; fitDS = true; break; } } }
      if (!keep && userInput.Flag.at("incMCTemp")) { for (auto& tmp : userInput.StrV.at("template")) { if ( (FILETAG.find(tmp)!=std::string::npos) ) { keep = true; fitDS = false; break; } } }
      if (keep) {
        dir = userInput.Par.at("localDSDir");
        if (userInput.Flag.at("useExtDS") && userInput.Par.at("extDSDir_MC")!="" && (existDir(userInput.Par.at("extDSDir_MC"))==true)) { dir = userInput.Par.at("extDSDir_MC"); }
      }
    }
    if (dir!="") {
      StringVectorMap FileInfo;
      FileInfo["InputFileNames"] = InputFileNames;
      FileInfo["OutputFileDir"].push_back(dir);
      FileInfo["OutputFileDir"].push_back(DIR["dataset"][0]);
      FileInfo["VarType"].push_back(userInput.Par.at("METType"));
      if (FILETAG.find("PA")==std::string::npos) { FileInfo["DSNames"].push_back(FILETAG); }
      else {
        std::string NAMETAG = FILETAG; NAMETAG.erase(NAMETAG.find("PA"),string("PA").length());
        if (userInput.Flag.at("fitpPb")) FileInfo["DSNames"].push_back(NAMETAG+"pPb");
        if (userInput.Flag.at("fitPbp")) FileInfo["DSNames"].push_back(NAMETAG+"Pbp");
      }
      if(!tree2DataSet(Workspace, FileInfo, userInput.Par.at("Analysis"))){ return; }
      if (fitDS) { for (auto& DSTAG : FileInfo["DSNames"]) { if (!aDSTAG->FindObject(DSTAG.c_str())) aDSTAG->Add(new TObjString(DSTAG.c_str())); } }
    }
  }
  if (Workspace.size()==0) {
    cout << "[ERROR] No tree files were found matching the user's input settings!" << endl; return;
  }

  // -------------------------------------------------------------------------------  
  // STEP 3: FIT THE DATASETS
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
    int index = ( DIR["output"].size()>1 ? j+1 : j ); // First entry is always the main output directory
    std::string outputDir = DIR["output"][index];

    for (unsigned int i=0; i < infoVectors[j].size(); i++) {
      nextDSTAG.Reset();
      TObjString* soDSTAG(0x0);
      while ( (soDSTAG = static_cast<TObjString*>(nextDSTAG.Next())) )
        {
          TString DSTAG = (TString)(soDSTAG->GetString());
          TString wsName = DSTAG;
          
          if (Workspace.count(wsName.Data())>0) {
            bool ispPb = true, isPbp = true, proceed = false;
            if (userInput.Flag.at("fitpPb") && DSTAG.Contains("pPb")) { ispPb = true;  proceed = true; }
            if (userInput.Flag.at("fitPbp") && DSTAG.Contains("Pbp")) { isPbp = false; proceed = true; }

            if (proceed && !fitElectroWeakMETModel( Workspace, infoVectors[j].at(i),
                                                    userInput, 
                                                    // Select the type of datasets to fit
                                                    outputDir,
                                                    DSTAG.Data()
                                                    )
                ) { return; }
          } else {
            cout << "[ERROR] The workspace for " << wsName.Data() << " was not found!" << endl; return;
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
    for(auto& row : data) {
      GlobalInfo info = GlobalInfo();
      if(!setParameters(row, info, Analysis)) { return false; }
      infoVector.push_back(info);
    }
  }
  else {
    if (data.size()!=infoVector.size()) { cout << "[ERROR] The initial parameters in file " << InputFile << " are not consistent with previous files!" << endl; return false; }
    for (unsigned int i=0; i<data.size(); i++) {
      GlobalInfo info = GlobalInfo();
      if (!setParameters(data.at(i), info, Analysis)) { return false; };
      if (info.Var != infoVector.at(i).Var) { cout << "[ERROR] The bins in file " << InputFile << " are not consistent with previous files!" << endl; return false; }
      infoVector.at(i).Copy(info.Par, true);
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
  info.Flag["incMCTemp"] = false;
  // set parameters from file
  for (auto& col : row) {
    std::string colName = col.first;
    bool found = false;
    for (auto& var : info.Var) {
      std::string varName = var.first;
      if (colName==varName) {
        if (col.second=="") {
          cout << "[ERROR] Input column " << varName << " has invalid value: " << col.second << endl; return false;
        }
        std::vector<double> v;
        if (!parseString(col.second, v)) { return false; }
        if (v.size()!=2 && v.size()!=1) {
          cout << "[ERROR] Input column " << varName << " has incorrect number of values, it should have 1 or 2 values but has: " << v.size() << endl; return false;
        }
        info.Var[varName]["Min"] = v.at(0);
        info.Var[varName]["Max"] = v.at(v.size()-1);
        found = true; break;
      }
    }
    if (found==false) {
      for (auto& par : info.Par) {
        std::string parName = par.first;
        if (colName.find(par.first)!=std::string::npos) {
          if (col.second=="") {
            cout << "[ERROR] Input column " << par.first << " has empty value" << endl; return false;
          }
          info.Par[colName] = col.second;
          if ( (par.first=="Model") && (col.second.find("TEMP")!=std::string::npos) ) { info.Flag.at("incMCTemp") = true; }
          found = true;
        }
      }
    }
    if (found==false) {
      if (col.second != "") {
	string value = col.second;
	// check that initial parameters format is correct: [ num, num, num ]
	if ((value.find("[")==std::string::npos)||(value.find("]")==std::string::npos)) {
	  // Special cases like parameter constrains could be set here but for now, let's keep it simple
	  cout << "[ERROR] Either ']' or '[' are missing in the initial parameter values!" << endl; return false;
	} else {
	  value.erase(0, value.find("[")+string("[").length());
	  value.erase(value.find("]"), value.length());
	}
	std::vector<double> v;
	if (!parseString(value, v)) { return false; }
        if (v.size()>3 || v.size()<1) {
          cout << "[ERROR] Initial parameter " << col.first << " has incorrect number of values, it has: " << v.size() << endl; return false;
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
  if (header.size()==0) { cout << "[ERROR] The header is null!" << endl; return false; }
  if(!readFile(FileName, content, header.size())){ return false; }
  for(auto& rHeader : header) {
    if (rHeader=="") { cout << "[ERROR] A column has no label!" << endl; return false; }
  }
  content.erase(content.begin()); // remove header
  for(auto& row : content) {
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


bool getInputFileNames(const std::string& InputTrees, std::map< std::string, std::vector< std::string > >& InputFileCollection)
{
  std::vector< std::vector< std::string > > content; 
  if(!readFile(InputTrees, content, 2)){ return false; }
  for(auto& row : content) {
    // remove spaces
    row.at(0).erase(std::remove(row.at(0).begin(), row.at(0).end(), ' '), row.at(0).end());
    row.at(1).erase(std::remove(row.at(1).begin(), row.at(1).end(), ' '), row.at(1).end());
    // remove tabs
    row.at(0).erase(std::remove(row.at(0).begin(), row.at(0).end(), '\t'), row.at(0).end());
    row.at(1).erase(std::remove(row.at(1).begin(), row.at(1).end(), '\t'), row.at(1).end());
    if (row.at(0)!="" && row.at(1)=="") { cout << "[ERROR] There is an empty file name in your InputTrees.txt, please fix it" << endl; return false; }
    if (row.at(0)=="" && row.at(1)!="") { cout << "[ERROR] There is an empty file tag in your InputTrees.txt, please fix it" << endl; return false; }
    if (row.at(0)!="" && row.at(1)!="") {
      // store the filenames mapped by the tag
      InputFileCollection[row.at(0)].push_back(row.at(1));
    }
  }
  return true;
};


bool readFile(std::string FileName, std::vector< std::vector< std::string > >& content, const int nCol, int nRow)
{
  if (nCol==0 || nRow==0) { 
    cout << "[WARNING] Ignoring content of File: " << FileName << endl; return true; 
  }
  if (nRow!=1) { cout << "[INFO] Reading file: " << FileName << endl; }
  ifstream myfile(FileName.c_str());
  char delimiter = ' ';
  if (myfile.is_open()){ 
    std::string line, CHAR;
    while ( getline(myfile, line) ){
      std::stringstream row(line), tmp(line); tmp >> CHAR;
      if ((!tmp) || (CHAR.find('#')!=std::string::npos) || (CHAR.find("//")!=std::string::npos)) continue;
      if (delimiter == ' ' && line.find(',')!=std::string::npos) { delimiter = ','; }
      if (delimiter == ' ' && line.find(';')!=std::string::npos) { delimiter = ';'; }
      if (delimiter == ' ') { cout << "[ERROR] File: " << FileName << " has unknown delimiter!" << endl; return false; }
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
    cout << "[ERROR] File: " << FileName << " was not found!" << endl; return false;
  }
  return true;
};


bool iniWorkEnv(std::map< std::string, std::vector< std::string > >& DIR, const std::string& workDirName)
{
  cout << "[INFO] Initializing the work enviroment" << endl;
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
    string subdir = DIR["input"][j];
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
      if (subdir->IsDirectory() && string(subdir->GetName())!="." && string(subdir->GetName())!="..") {
        dirlist.push_back(dirname + subdir->GetName() + "/");
        cout << "[INFO] Input subdirectory: " << dirname + subdir->GetName() + "/" << " found!" << endl;
      }
    }
  }
  delete subdirs;
  return;
};


bool checkSettings(const GlobalInfo& userInput)
{ 
  cout << "[INFO] Checking user settings " << endl;

  if (userInput.Flag.at("fitW") && userInput.Flag.at("fitQCD")) { std::cout << "[ERROR] We can not fit QCD and W at the same time!" << std::endl; return false; }

  cout << "[INFO] All user setting are correct " << endl;
  return true;
};

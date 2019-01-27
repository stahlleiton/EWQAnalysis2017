// Auxiliary Headers
#include "../Utilities/EVENTUTILS.h"
#include "../Fitter/Macros/EWQ/drawElectroWeakMETPlot.C"
// ROOT headers
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TDirectory.h"
// RooFit headers
#include "RooWorkspace.h"
// c++ headers
#include <iostream>
#include <vector>
#include <map>
#include <string>
// CMS headers


// ------------------ TYPE -------------------------------


/////////////////////
// OTHER FUNCTIONS //
/////////////////////

bool fileList ( std::vector< std::string >& fileNames , const std::string& dirPath );

/////////////////////


void makeFitPlots(
                  const int plotSTYLE = 3, // 4: Thesis (LIN) , 3: Thesis (LOG) , 2: Paper , 1: PAS , 0: AN
                  const std::string workDirName="NominalCM",
                  const std::vector< std::string > collVec = { "PA" /*, "pPb" , "Pbp"  */}
                  )
{
  //
  // Define general info
  const std::string metTag = "METPF_RAW";
  const std::string dsTag  = "DATA";
  const std::string CWD = getcwd(NULL, 0);
  const std::string par = "W";
  std::string preCWD = CWD; preCWD.erase(preCWD.find_last_of("/"), 10);
  //
  // Loop over each collision system
  for (const auto& colTag : collVec) {
    //
    // Find the fit directory
    const std::string wsDirPath = Form("%s/Fitter/Output/%s/%s/%s/%s/%s/result", preCWD.c_str(), workDirName.c_str(), metTag.c_str(), dsTag.c_str(), par.c_str(), colTag.c_str());
    if (existDir(wsDirPath)==false) { std::cout << "[WARNING] Workspace directory " << wsDirPath << " was not found, will skip it!" << std::endl; return; }
    //
    // --------------------------------------------------------------------------------- //
    //
    // Get the list of input files
    //
    std::vector< std::string > inputFileNames;
    const std::string inputDirPath = wsDirPath;
    if (!fileList(inputFileNames, inputDirPath)) { return; };
    //
    // --------------------------------------------------------------------------------- //
    //
    // Loop over the input files
    //
    for (const auto& inputFileName : inputFileNames) {
      //
      std::cout << "Processing file: " << inputFileName << std::endl;
      //
      // Open input file
      const std::string inputFilePath = Form("%s/%s", inputDirPath.c_str(), inputFileName.c_str());
      TFile inputFile(inputFilePath.c_str(), "READ");
      //
      if (inputFile.IsOpen()==false || inputFile.IsZombie()==true) {
        std::cout << "[ERROR] The input file " << inputFilePath << " could not be open!" << std::endl; return;
      }
      //
      // Extract the Workspace
      RooWorkspace* ws = (RooWorkspace*) inputFile.Get("workspace");
      if (ws == NULL) { std::cout << "[ERROR] File: " << inputFilePath << " does not have the workspace!" << std::endl; inputFile.Close(); return; }
      //
      // Get the information needed
      //
      // Output File Name
      std::string fileName = inputFileName;
      fileName = fileName.substr(4, fileName.find(".root")-4);
      fileName = "PLOT_" + fileName;
      // Output Directory Name
      std::string outputDir = "";
      int plotStyle = plotSTYLE;
      if      (plotStyle==4) { outputDir =  Form("%s/Output/%s/%s/%s/FitPlots/ForThesis_LIN/%s/%s/",CWD.c_str(), workDirName.c_str(), metTag.c_str(), dsTag.c_str(), par.c_str(), colTag.c_str()); }
      else if (plotStyle==3) { outputDir =  Form("%s/Output/%s/%s/%s/FitPlots/ForThesis_LOG/%s/%s/",CWD.c_str(), workDirName.c_str(), metTag.c_str(), dsTag.c_str(), par.c_str(), colTag.c_str()); }
      else if (plotStyle==2) { outputDir =  Form("%s/Output/%s/%s/%s/FitPlots/ForPaper/%s/%s/", CWD.c_str(), workDirName.c_str(), metTag.c_str(), dsTag.c_str(), par.c_str(), colTag.c_str()); }
      else if (plotStyle==1) { outputDir =  Form("%s/Output/%s/%s/%s/FitPlots/ForPAS/%s/%s/", CWD.c_str(), workDirName.c_str(), metTag.c_str(), dsTag.c_str(), par.c_str(), colTag.c_str());   }
      else                   { outputDir =  Form("%s/Output/%s/%s/%s/FitPlots/ForAN/%s/%s/", CWD.c_str(), workDirName.c_str(), metTag.c_str(), dsTag.c_str(), par.c_str(), colTag.c_str());    }
      // Fit plotting settings
      const bool yLogScale = (plotStyle==0 || plotStyle==3); //
      double maxRng = ((plotStyle==1 || plotStyle==2 || plotStyle==4) ? 80. : 150.);
      const bool doGoF = true;
      const bool redoFrame = (plotStyle==1 || plotStyle==2 || plotStyle==4);
      const bool doQCDHist = (plotStyle==1 || plotStyle==2);
      if (plotSTYLE==4) { plotStyle = 3; }
      //
      // Draw the output plot
      if (!drawElectroWeakMETPlot(*ws, fileName, outputDir, yLogScale, maxRng, doGoF, plotStyle, redoFrame, doQCDHist)) { return; }
      //
    }
  }
};


bool fileList(std::vector< std::string >& fileNames, const std::string& dirPath)
{
  // Open the directory
  DIR * dpdf = opendir(dirPath.c_str());
  // Search for all the files inside the directory
  if (dpdf != NULL){
    struct dirent *epdf;
    while ((epdf = readdir(dpdf))){
      if (strcmp(epdf->d_name,".")!=0 && strcmp(epdf->d_name,"..")!=0 ) {
        //std::cout << "[INFO] Adding file: " << epdf->d_name << std::endl;
        fileNames.push_back(epdf->d_name);
      }
    }
  } else {
    std::cout << "[ERROR] Working directory ( " << dirPath << " ) was not found!" << endl; return false;
  }
  return true;
};

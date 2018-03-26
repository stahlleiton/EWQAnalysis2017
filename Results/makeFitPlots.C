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
                  const bool paperStyle = true,
                  const std::string workDirName="NominalCM",
                  const std::vector< std::string > collVec = { "PA" /*, "pPb" , "Pbp"  */}
                  )
{
  //
  // Define general info
  const std::string metTag = "METPF_RAW";
  const std::string dsTag  = "DATA";
  const std::string CWD = getcwd(NULL, 0);
  std::string preCWD = CWD; preCWD.erase(preCWD.find_last_of("/"), 10);
  //
  // Loop over each collision system
  for (const auto& colTag : collVec) {
    //
    // Find the fit directory
    const std::string wsDirPath = Form("%s/Fitter/Output/%s/%s/%s/W/%s/result", preCWD.c_str(), workDirName.c_str(), metTag.c_str(), dsTag.c_str(), colTag.c_str());
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
      if (paperStyle) { outputDir =  Form("%s/Output/%s/%s/%s/FitPlots/ForPaper/%s/", CWD.c_str(), workDirName.c_str(), metTag.c_str(), dsTag.c_str(), colTag.c_str()); }
      else            { outputDir =  Form("%s/Output/%s/%s/%s/FitPlots/ForAN/%s/", CWD.c_str(), workDirName.c_str(), metTag.c_str(), dsTag.c_str(), colTag.c_str());    }
      // Fit plotting settings
      const bool yLogScale = !paperStyle;
      const double maxRng = 150.;
      const bool doGoF = !paperStyle;
      //
      // Draw the output plot
      if (!drawElectroWeakMETPlot(*ws, fileName, outputDir, yLogScale, maxRng, doGoF, paperStyle)) { return; }
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

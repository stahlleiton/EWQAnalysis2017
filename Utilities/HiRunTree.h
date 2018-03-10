#ifndef HiRunTree_h
#define HiRunTree_h

// Header file for ROOT classes
#include <TROOT.h>
#include <TChain.h>
#include <TInterpreter.h>
#include <TDirectory.h>
#include <TSystem.h>
#include <TFile.h>

// Header file for c++ classes
#include <iostream>
#include <vector>
#include <tuple>
#include <map>
#include <string>

// Header file for the classes stored in the TChain


class HiRunTree {

public :

  HiRunTree();
  virtual ~HiRunTree();
  virtual Bool_t            GetTree         (const std::vector< std::string >&, const std::string& treeName="runAnalyzer");
  virtual Bool_t            GetTree         (const std::string&, const std::string& treeName="runAnalyzer");
  virtual Int_t             GetEntry        (Long64_t);
  virtual Long64_t          GetEntries      (void) { return fChain_->GetEntries(); }
  virtual Long64_t          GetTreeEntries  (void) { return fChain_->GetTree()->GetEntriesFast(); }
  virtual TChain*           Chain           (void) { return fChain_; }
  virtual void              Clear           (void);

  // VARIABLES GETTER
  Float_t                   xsec() { SetBranch("xsec"); return xsec_; }

 private:

  virtual Long64_t          LoadTree        (Long64_t);
  virtual void              SetBranch       (const std::string&);
  virtual void              InitTree        (void);
  virtual Int_t             LoadEntry       (void) { return fChain_->GetEntry(entry_); }

  TChain*                   fChain_;
  std::map<string, TChain*> fChainM_;
  Int_t                     fCurrent_ = -1;
  Long64_t                  entry_;

  // VARIABLES
  Float_t                   xsec_ = 0.;

  // VARIABLE BRANCHES
  TBranch                  *b_xsec;   //!
};

HiRunTree::HiRunTree() : fChain_(0)
{
}

HiRunTree::~HiRunTree()
{
  if (fChain_ && fChain_->GetCurrentFile()) delete fChain_->GetCurrentFile();
  for (auto& c : fChainM_) { if (c.second) { c.second->Reset(); delete c.second; } }
}

Bool_t HiRunTree::GetTree(const std::string& fileName, const std::string& treeName)
{
  std::vector<std::string> fileNames = {fileName};
  return GetTree(fileNames, treeName);
}

Bool_t HiRunTree::GetTree(const std::vector< std::string >& inFileName, const std::string& treeName)
{
  // Check the File Names
  std::vector< std::string > fileName = inFileName;
  for (auto& f : fileName) { if (f.find("/store/")!=std::string::npos && f.find("root://")==std::string::npos) { f = "root://cms-xrd-global.cern.ch/" + f.substr(f.find("/store/")); } }
  // Open the input files
  TFile *f = TFile::Open(fileName[0].c_str());
  if (!f || !f->IsOpen()) return false;
  // Extract the input TChains
  fChainM_.clear();
  TDirectory * dir;
  if (fileName[0].find("root://")!=std::string::npos) { dir = (TDirectory*)f->Get(treeName.c_str()); }
  else { dir = (TDirectory*)f->Get((fileName[0]+":/"+treeName).c_str()); }
  if (!dir) return false;
  if (dir->GetListOfKeys()->Contains("run")) { fChainM_["run"] = new TChain((treeName+"/run").c_str()  , "run"); }
  if (fChainM_.size()==0) return false;
  // Add the files in the TChain
  for (auto& c : fChainM_) {
    for (auto& f : fileName) { c.second->Add(Form("%s/%s/%s", f.c_str(), treeName.c_str(), c.first.c_str())); }; c.second->GetEntries();
  }
  for (auto& c : fChainM_) { if (!c.second) { std::cout << "[ERROR] fChain " << c.first << " was not created, some input files are missing" << std::endl; return false; } }
  // Initialize the input TChains (set their branches)
  InitTree();
  // Add Friend TChains
  fChain_ = (TChain*)fChainM_.begin()->second->Clone(Form("%s", treeName.c_str()));
  for (auto& c : fChainM_) {
    c.second->SetMakeClass(1); // For the proper setup.
    if (fChain_!=c.second) { fChain_->AddFriend(c.second, Form("%s", c.first.c_str()), kTRUE); } // Add the Friend TChain
  }
  if (fChain_ == 0) return false;
  // Set All Branches to Status 0
  fChain_->SetBranchStatus("*",0);
  //
  return true;
}

Int_t HiRunTree::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  entry_ = entry;
  if (LoadTree(entry_) < 0) return -1;
  Clear();
  return LoadEntry();
}

Long64_t HiRunTree::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain_) return -5;
  Long64_t centry = fChain_->LoadTree(entry);
  if (fChain_->GetTreeNumber() != fCurrent_) { fCurrent_ = fChain_->GetTreeNumber(); }
  return centry;
}

void HiRunTree::SetBranch(const std::string& n)
{
  if ( fChain_->GetBranch(n.c_str()) && (fChain_->GetBranchStatus(n.c_str()) == 0) ) {
    fChain_->SetBranchStatus(Form("*%s*", n.c_str()), 1);
    LoadEntry();
  }
}

void HiRunTree::InitTree(void)
{
  // INITIALIZE POINTERS
  xsec_ = 0.;

  if (fChainM_.size()==0) return;

  // SET VARIABLE BRANCHES
  if (fChainM_.count("run")>0) {
    if (fChainM_.at("run")->GetBranch("xsec")) fChainM_.at("run")->SetBranchAddress("xsec", &xsec_, &b_xsec);
  }
}

void HiRunTree::Clear(void)
{
  if (fChainM_.size()==0) return;

  // CLEAR VARIABLES
  xsec_ = 0.;
}

#endif

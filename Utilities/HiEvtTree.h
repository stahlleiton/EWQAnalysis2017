#ifndef HiEvtTree_h
#define HiEvtTree_h

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


class HiEvtTree {

public :

  HiEvtTree();
  virtual ~HiEvtTree();
  virtual Bool_t            GetTree         (const std::vector< std::string >&, const std::string& treeName="hiEvtAna");
  virtual Bool_t            GetTree         (const std::string&, const std::string& treeName="hiEvtAna");
  virtual Int_t             GetEntry        (Long64_t);
  virtual Long64_t          GetEntries      (void) { return fChain_->GetEntries(); }
  virtual Long64_t          GetTreeEntries  (void) { return fChain_->GetTree()->GetEntriesFast(); }
  virtual TChain*           Chain           (void) { return fChain_; }
  virtual void              Clear           (void);

  // VARIABLES GETTER
  UInt_t                    run()                 { SetBranch("run");                 return run_;                 }
  ULong64_t                 evt()                 { SetBranch("evt");                 return evt_;                 }
  UInt_t                    lumi()                { SetBranch("lumi");                return lumi_;                }
  Float_t                   vx()                  { SetBranch("vx");                  return vx_;                  }
  Float_t                   vy()                  { SetBranch("vy");                  return vy_;                  }
  Float_t                   vz()                  { SetBranch("vz");                  return vz_;                  }
  Int_t                     ProcessID()           { SetBranch("ProcessID");           return ProcessID_;           }
  Float_t                   pthat()               { SetBranch("pthat");               return pthat_;               }
  Float_t                   weight()              { SetBranch("weight");              return weight_;              }
  Float_t                   alphaQCD()            { SetBranch("alphaQCD");            return alphaQCD_;            }
  Float_t                   alphaQED()            { SetBranch("alphaQED");            return alphaQED_;            }
  Float_t                   qScale()              { SetBranch("qScale");              return qScale_;              }
  Int_t                     nMEPartons()          { SetBranch("nMEPartons");          return nMEPartons_;          }
  Int_t                     nMEPartonsFiltered()  { SetBranch("nMEPartonsFiltered");  return nMEPartonsFiltered_;  }
  std::pair<int, int>       pdfID()               { SetBranch("pdfID");               return GET(pdfID_);          }
  std::pair<float, float>   pdfX()                { SetBranch("pdfX");                return GET(pdfX_);           }
  std::pair<float, float>   pdfXpdf()             { SetBranch("pdfXpdf");             return GET(pdfXpdf_);        }
  std::vector<float>        ttbar_w()             { SetBranch("ttbar_w");             return GET(ttbar_w_);        }
  std::vector<int>          npus()                { SetBranch("npus");                return GET(npus_);           }
  std::vector<float>        tnpus()               { SetBranch("tnpus");               return GET(tnpus_);          }
  Int_t                     hiBin()               { SetBranch("hiBin");               return hiBin_;               }
  Float_t                   hiHF()                { SetBranch("hiHF");                return hiHF_;                }
  Float_t                   hiHFplus()            { SetBranch("hiHFplus");            return hiHFplus_;            }
  Float_t                   hiHFminus()           { SetBranch("hiHFminus");           return hiHFminus_;           }
  Float_t                   hiHFplusEta4()        { SetBranch("hiHFplusEta4");        return hiHFplusEta4_;        }
  Float_t                   hiHFminusEta4()       { SetBranch("hiHFminusEta4");       return hiHFminusEta4_;       }
  Float_t                   hiZDC()               { SetBranch("hiZDC");               return hiZDC_;               }
  Float_t                   hiZDCplus()           { SetBranch("hiZDCplus");           return hiZDCplus_;           }
  Float_t                   hiZDCminus()          { SetBranch("hiZDCminus");          return hiZDCminus_;          }
  Float_t                   hiHFhit()             { SetBranch("hiHFhit");             return hiHFhit_;             }
  Float_t                   hiHFhitPlus()         { SetBranch("hiHFhitPlus");         return hiHFhitPlus_;         }
  Float_t                   hiHFhitMinus()        { SetBranch("hiHFhitMinus");        return hiHFhitMinus_;        }
  Float_t                   hiET()                { SetBranch("hiET");                return hiET_;                }
  Float_t                   hiEE()                { SetBranch("hiEE");                return hiEE_;                }
  Float_t                   hiEB()                { SetBranch("hiEB");                return hiEB_;                }
  Float_t                   hiEEplus()            { SetBranch("hiEEplus");            return hiEEplus_;            }
  Float_t                   hiEEminus()           { SetBranch("hiEEminus");           return hiEEminus_;           }
  Int_t                     hiNpix()              { SetBranch("hiNpix");              return hiNpix_;              }
  Int_t                     hiNpixelTracks()      { SetBranch("hiNpixelTracks");      return hiNpixelTracks_;      }
  Int_t                     hiNtracks()           { SetBranch("hiNtracks");           return hiNtracks_;           }
  Int_t                     hiNtracksPtCut()      { SetBranch("hiNtracksPtCut");      return hiNtracksPtCut_;      }
  Int_t                     hiNtracksEtaCut()     { SetBranch("hiNtracksEtaCut");     return hiNtracksEtaCut_;     }
  Int_t                     hiNtracksEtaPtCut()   { SetBranch("hiNtracksEtaPtCut");   return hiNtracksEtaPtCut_;   }

 private:

  virtual Long64_t          LoadTree        (Long64_t);
  virtual void              SetBranch       (const std::string&);
  virtual void              InitTree        (void);
  virtual Int_t             LoadEntry       (void) { return fChain_->GetEntry(entry_); }
  virtual void              GenerateDictionaries (void);

  template <typename T> 
    T GET(T* x) { return ( (x) ? *x : T() ); }


  TChain*                   fChain_;
  std::map<string, TChain*> fChainM_;
  Int_t                     fCurrent_ = -1;
  Long64_t                  entry_;

  // VARIABLES
  UInt_t                    run_                = 0;
  ULong64_t                 evt_                = 0;
  UInt_t                    lumi_               = 0;
  Float_t                   vx_                 = -99.;
  Float_t                   vy_                 = -99.;
  Float_t                   vz_                 = -99.;
  Int_t                     ProcessID_          = -1;
  Float_t                   pthat_              = -1.;
  Float_t                   weight_             = 0.;
  Float_t                   alphaQCD_           = -1.;
  Float_t                   alphaQED_           = -1.;
  Float_t                   qScale_             = -1.;
  Int_t                     nMEPartons_         = -1;
  Int_t                     nMEPartonsFiltered_ = -1;
  std::pair<int, int>      *pdfID_;
  std::pair<float, float>  *pdfX_;
  std::pair<float, float>  *pdfXpdf_;
  std::vector<float>       *ttbar_w_; //weights for systematics
  std::vector<int>         *npus_;    //number of pileup interactions
  std::vector<float>       *tnpus_;   //true number of interactions
  Int_t                     hiBin_              = -1;
  Float_t                   hiHF_               = -1.;
  Float_t                   hiHFplus_           = -1.;
  Float_t                   hiHFminus_          = -1.;
  Float_t                   hiHFplusEta4_       = -1.;
  Float_t                   hiHFminusEta4_      = -1.;
  Float_t                   hiZDC_              = -1.;
  Float_t                   hiZDCplus_          = -1.;
  Float_t                   hiZDCminus_         = -1.;
  Float_t                   hiHFhit_            = -1.;
  Float_t                   hiHFhitPlus_        = -1.;
  Float_t                   hiHFhitMinus_       = -1.;
  Float_t                   hiET_               = -1.;
  Float_t                   hiEE_               = -1.;
  Float_t                   hiEB_               = -1.;
  Float_t                   hiEEplus_           = -1.;
  Float_t                   hiEEminus_          = -1.;
  Int_t                     hiNpix_             = -1;
  Int_t                     hiNpixelTracks_     = -1;
  Int_t                     hiNtracks_          = -1;
  Int_t                     hiNtracksPtCut_     = -1;
  Int_t                     hiNtracksEtaCut_    = -1;
  Int_t                     hiNtracksEtaPtCut_  = -1;

  // VARIABLE BRANCHES
  TBranch                  *b_run;   //!
  TBranch                  *b_evt;   //!
  TBranch                  *b_lumi;   //!
  TBranch                  *b_vx;   //!
  TBranch                  *b_vy;   //!
  TBranch                  *b_vz;   //!
  TBranch                  *b_ProcessID;   //!
  TBranch                  *b_pthat;   //!
  TBranch                  *b_weight;   //!
  TBranch                  *b_alphaQCD;   //!
  TBranch                  *b_alphaQED;   //!
  TBranch                  *b_qScale;   //!
  TBranch                  *b_nMEPartons;   //!
  TBranch                  *b_nMEPartonsFiltered;   //!
  TBranch                  *b_pdfID;   //!
  TBranch                  *b_pdfX;   //!
  TBranch                  *b_pdfXpdf;   //!
  TBranch                  *b_ttbar_w;   //!
  TBranch                  *b_npus;   //!
  TBranch                  *b_tnpus;   //!
  TBranch                  *b_hiBin;   //!
  TBranch                  *b_hiHF;   //!
  TBranch                  *b_hiHFplus;   //!
  TBranch                  *b_hiHFminus;   //!
  TBranch                  *b_hiHFplusEta4;   //!
  TBranch                  *b_hiHFminusEta4;   //!
  TBranch                  *b_hiZDC;   //!
  TBranch                  *b_hiZDCplus;   //!
  TBranch                  *b_hiZDCminus;   //!
  TBranch                  *b_hiHFhit;   //!
  TBranch                  *b_hiHFhitPlus;   //!
  TBranch                  *b_hiHFhitMinus;   //!
  TBranch                  *b_hiET;   //!
  TBranch                  *b_hiEE;   //!
  TBranch                  *b_hiEB;   //!
  TBranch                  *b_hiEEplus;   //!
  TBranch                  *b_hiEEminus;   //!
  TBranch                  *b_hiNpix;   //!
  TBranch                  *b_hiNpixelTracks;   //!
  TBranch                  *b_hiNtracks;   //!
  TBranch                  *b_hiNtracksPtCut;   //!
  TBranch                  *b_hiNtracksEtaCut;   //!
  TBranch                  *b_hiNtracksEtaPtCut;   //!
};

HiEvtTree::HiEvtTree() : fChain_(0)
{
}

HiEvtTree::~HiEvtTree()
{
  if (fChain_ && fChain_->GetCurrentFile()) delete fChain_->GetCurrentFile();
  for (auto& c : fChainM_) { if (c.second) { c.second->Reset(); delete c.second; } }
}

Bool_t HiEvtTree::GetTree(const std::string& fileName, const std::string& treeName)
{
  std::vector<std::string> fileNames = {fileName};
  return GetTree(fileNames, treeName);
}

Bool_t HiEvtTree::GetTree(const std::vector< std::string >& inFileName, const std::string& treeName)
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
  if (dir->GetListOfKeys()->Contains("HiTree")) { fChainM_["HiTree"] = new TChain((treeName+"/HiTree").c_str()  , "HiTree"); }
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

Int_t HiEvtTree::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  entry_ = entry;
  if (LoadTree(entry_) < 0) return -1;
  Clear();
  return LoadEntry();
}

Long64_t HiEvtTree::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain_) return -5;
  Long64_t centry = fChain_->LoadTree(entry);
  if (fChain_->GetTreeNumber() != fCurrent_) { fCurrent_ = fChain_->GetTreeNumber(); }
  return centry;
}

void HiEvtTree::SetBranch(const std::string& n)
{
  if ( fChain_->GetBranch(n.c_str()) && (fChain_->GetBranchStatus(n.c_str()) == 0) ) {
    fChain_->SetBranchStatus(Form("*%s*", n.c_str()), 1);
    LoadEntry();
  }
  if (n=="pdfX" || n=="pdfID" || n=="pdfXpdf") { SetBranch("first"); SetBranch("second"); }
}

void HiEvtTree::InitTree(void)
{
  // Generate the dictionary's needed
  GenerateDictionaries();

  // INITIALIZE POINTERS
  pdfID_              = 0;
  pdfX_               = 0;
  pdfXpdf_            = 0;
  ttbar_w_            = 0;
  npus_               = 0;
  tnpus_              = 0;
 
  // INITIALIZE VARIABLES 
  run_                = 0;
  evt_                = 0;
  lumi_               = 0;
  vx_                 = -99.;
  vy_                 = -99.;
  vz_                 = -99.;
  ProcessID_          = -1;
  pthat_              = -1.;
  weight_             = 0.;
  alphaQCD_           = -1.;
  alphaQED_           = -1.;
  qScale_             = -1.;
  nMEPartons_         = -1;
  nMEPartonsFiltered_ = -1;
  hiBin_              = -1;
  hiHF_               = -1.;
  hiHFplus_           = -1.;
  hiHFminus_          = -1.;
  hiHFplusEta4_       = -1.;
  hiHFminusEta4_      = -1.;
  hiZDC_              = -1.;
  hiZDCplus_          = -1.;
  hiZDCminus_         = -1.;
  hiHFhit_            = -1.;
  hiHFhitPlus_        = -1.;
  hiHFhitMinus_       = -1.;
  hiET_               = -1.;
  hiEE_               = -1.;
  hiEB_               = -1.;
  hiEEplus_           = -1.;
  hiEEminus_          = -1.;
  hiNpix_             = -1;
  hiNpixelTracks_     = -1;
  hiNtracks_          = -1;
  hiNtracksPtCut_     = -1;
  hiNtracksEtaCut_    = -1;
  hiNtracksEtaPtCut_  = -1;

  if (fChainM_.size()==0) return;

  // SET VARIABLE BRANCHES
  if (fChainM_.count("HiTree")>0) {
    if (fChainM_.at("HiTree")->GetBranch("run"))                fChainM_.at("HiTree")->SetBranchAddress("run",                &run_,                &b_run);
    if (fChainM_.at("HiTree")->GetBranch("evt"))                fChainM_.at("HiTree")->SetBranchAddress("evt",                &evt_,                &b_evt);
    if (fChainM_.at("HiTree")->GetBranch("lumi"))               fChainM_.at("HiTree")->SetBranchAddress("lumi",               &lumi_,               &b_lumi);
    if (fChainM_.at("HiTree")->GetBranch("vx"))                 fChainM_.at("HiTree")->SetBranchAddress("vx",                 &vx_,                 &b_vx);
    if (fChainM_.at("HiTree")->GetBranch("vy"))                 fChainM_.at("HiTree")->SetBranchAddress("vy",                 &vy_,                 &b_vy);
    if (fChainM_.at("HiTree")->GetBranch("vz"))                 fChainM_.at("HiTree")->SetBranchAddress("vz",                 &vz_,                 &b_vz);
    if (fChainM_.at("HiTree")->GetBranch("ProcessID"))          fChainM_.at("HiTree")->SetBranchAddress("ProcessID",          &ProcessID_,          &b_ProcessID);
    if (fChainM_.at("HiTree")->GetBranch("pthat"))              fChainM_.at("HiTree")->SetBranchAddress("pthat",              &pthat_,              &b_pthat);
    if (fChainM_.at("HiTree")->GetBranch("weight"))             fChainM_.at("HiTree")->SetBranchAddress("weight",             &weight_,             &b_weight);
    if (fChainM_.at("HiTree")->GetBranch("alphaQCD"))           fChainM_.at("HiTree")->SetBranchAddress("alphaQCD",           &alphaQCD_,           &b_alphaQCD);
    if (fChainM_.at("HiTree")->GetBranch("alphaQED"))           fChainM_.at("HiTree")->SetBranchAddress("alphaQED",           &alphaQED_,           &b_alphaQED);
    if (fChainM_.at("HiTree")->GetBranch("qScale"))             fChainM_.at("HiTree")->SetBranchAddress("qScale",             &qScale_,             &b_qScale);
    if (fChainM_.at("HiTree")->GetBranch("nMEPartons"))         fChainM_.at("HiTree")->SetBranchAddress("nMEPartons",         &nMEPartons_,         &b_nMEPartons);
    if (fChainM_.at("HiTree")->GetBranch("nMEPartonsFiltered")) fChainM_.at("HiTree")->SetBranchAddress("nMEPartonsFiltered", &nMEPartonsFiltered_, &b_nMEPartonsFiltered);
    if (fChainM_.at("HiTree")->GetBranch("pdfID"))              fChainM_.at("HiTree")->SetBranchAddress("pdfID",              &pdfID_,              &b_pdfID);
    if (fChainM_.at("HiTree")->GetBranch("pdfX"))               fChainM_.at("HiTree")->SetBranchAddress("pdfX",               &pdfX_,               &b_pdfX);
    if (fChainM_.at("HiTree")->GetBranch("pdfXpdf"))            fChainM_.at("HiTree")->SetBranchAddress("pdfXpdf",            &pdfXpdf_,            &b_pdfXpdf);
    if (fChainM_.at("HiTree")->GetBranch("ttbar_w"))            fChainM_.at("HiTree")->SetBranchAddress("ttbar_w",            &ttbar_w_,            &b_ttbar_w);
    if (fChainM_.at("HiTree")->GetBranch("npus"))               fChainM_.at("HiTree")->SetBranchAddress("npus",               &npus_,               &b_npus);
    if (fChainM_.at("HiTree")->GetBranch("tnpus"))              fChainM_.at("HiTree")->SetBranchAddress("tnpus",              &tnpus_,              &b_tnpus);
    if (fChainM_.at("HiTree")->GetBranch("hiBin"))              fChainM_.at("HiTree")->SetBranchAddress("hiBin",              &hiBin_,              &b_hiBin);
    if (fChainM_.at("HiTree")->GetBranch("hiHF"))               fChainM_.at("HiTree")->SetBranchAddress("hiHF",               &hiHF_,               &b_hiHF);
    if (fChainM_.at("HiTree")->GetBranch("hiHFplus"))           fChainM_.at("HiTree")->SetBranchAddress("hiHFplus",           &hiHFplus_,           &b_hiHFplus);
    if (fChainM_.at("HiTree")->GetBranch("hiHFminus"))          fChainM_.at("HiTree")->SetBranchAddress("hiHFminus",          &hiHFminus_,          &b_hiHFminus);
    if (fChainM_.at("HiTree")->GetBranch("hiHFplusEta4"))       fChainM_.at("HiTree")->SetBranchAddress("hiHFplusEta4",       &hiHFplusEta4_,       &b_hiHFplusEta4);
    if (fChainM_.at("HiTree")->GetBranch("hiHFminusEta4"))      fChainM_.at("HiTree")->SetBranchAddress("hiHFminusEta4",      &hiHFminusEta4_,      &b_hiHFminusEta4);
    if (fChainM_.at("HiTree")->GetBranch("hiZDC"))              fChainM_.at("HiTree")->SetBranchAddress("hiZDC",              &hiZDC_,              &b_hiZDC);
    if (fChainM_.at("HiTree")->GetBranch("hiZDCplus"))          fChainM_.at("HiTree")->SetBranchAddress("hiZDCplus",          &hiZDCplus_,          &b_hiZDCplus);
    if (fChainM_.at("HiTree")->GetBranch("hiZDCminus"))         fChainM_.at("HiTree")->SetBranchAddress("hiZDCminus",         &hiZDCminus_,         &b_hiZDCminus);
    if (fChainM_.at("HiTree")->GetBranch("hiHFhit"))            fChainM_.at("HiTree")->SetBranchAddress("hiHFhit",            &hiHFhit_,            &b_hiHFhit);
    if (fChainM_.at("HiTree")->GetBranch("hiHFhitPlus"))        fChainM_.at("HiTree")->SetBranchAddress("hiHFhitPlus",        &hiHFhitPlus_,        &b_hiHFhitPlus);
    if (fChainM_.at("HiTree")->GetBranch("hiHFhitMinus"))       fChainM_.at("HiTree")->SetBranchAddress("hiHFhitMinus",       &hiHFhitMinus_,       &b_hiHFhitMinus);
    if (fChainM_.at("HiTree")->GetBranch("hiET"))               fChainM_.at("HiTree")->SetBranchAddress("hiET",               &hiET_,               &b_hiET);
    if (fChainM_.at("HiTree")->GetBranch("hiEE"))               fChainM_.at("HiTree")->SetBranchAddress("hiEE",               &hiEE_,               &b_hiEE);
    if (fChainM_.at("HiTree")->GetBranch("hiEB"))               fChainM_.at("HiTree")->SetBranchAddress("hiEB",               &hiEB_,               &b_hiEB);
    if (fChainM_.at("HiTree")->GetBranch("hiEEplus"))           fChainM_.at("HiTree")->SetBranchAddress("hiEEplus",           &hiEEplus_,           &b_hiEEplus);
    if (fChainM_.at("HiTree")->GetBranch("hiEEminus"))          fChainM_.at("HiTree")->SetBranchAddress("hiEEminus",          &hiEEminus_,          &b_hiEEminus);
    if (fChainM_.at("HiTree")->GetBranch("hiNpix"))             fChainM_.at("HiTree")->SetBranchAddress("hiNpix",             &hiNpix_,             &b_hiNpix);
    if (fChainM_.at("HiTree")->GetBranch("hiNpixelTracks"))     fChainM_.at("HiTree")->SetBranchAddress("hiNpixelTracks",     &hiNpixelTracks_,     &b_hiNpixelTracks);
    if (fChainM_.at("HiTree")->GetBranch("hiNtracks"))          fChainM_.at("HiTree")->SetBranchAddress("hiNtracks",          &hiNtracks_,          &b_hiNtracks);
    if (fChainM_.at("HiTree")->GetBranch("hiNtracksPtCut"))     fChainM_.at("HiTree")->SetBranchAddress("hiNtracksPtCut",     &hiNtracksPtCut_,     &b_hiNtracksPtCut);
    if (fChainM_.at("HiTree")->GetBranch("hiNtracksEtaCut"))    fChainM_.at("HiTree")->SetBranchAddress("hiNtracksEtaCut",    &hiNtracksEtaCut_,    &b_hiNtracksEtaCut);
    if (fChainM_.at("HiTree")->GetBranch("hiNtracksEtaPtCut"))  fChainM_.at("HiTree")->SetBranchAddress("hiNtracksEtaPtCut",  &hiNtracksEtaPtCut_,  &b_hiNtracksEtaPtCut);
  }
}

void HiEvtTree::Clear(void)
{
  if (fChainM_.size()==0) return;

  // CLEAR VARIABLES
  run_                = 0;
  evt_                = 0;
  lumi_               = 0;
  vx_                 = -99.;
  vy_                 = -99.;
  vz_                 = -99.;
  ProcessID_          = -1;
  pthat_              = -1.;
  weight_             = 0.;
  alphaQCD_           = -1.;
  alphaQED_           = -1.;
  qScale_             = -1.;
  nMEPartons_         = -1;
  nMEPartonsFiltered_ = -1;
  hiBin_              = -1;
  hiHF_               = -1.;
  hiHFplus_           = -1.;
  hiHFminus_          = -1.;
  hiHFplusEta4_       = -1.;
  hiHFminusEta4_      = -1.;
  hiZDC_              = -1.;
  hiZDCplus_          = -1.;
  hiZDCminus_         = -1.;
  hiHFhit_            = -1.;
  hiHFhitPlus_        = -1.;
  hiHFhitMinus_       = -1.;
  hiET_               = -1.;
  hiEE_               = -1.;
  hiEB_               = -1.;
  hiEEplus_           = -1.;
  hiEEminus_          = -1.;
  hiNpix_             = -1;
  hiNpixelTracks_     = -1;
  hiNtracks_          = -1;
  hiNtracksPtCut_     = -1;
  hiNtracksEtaCut_    = -1;
  hiNtracksEtaPtCut_  = -1;

  // CLEAR PAIRS
  if (pdfID_)   { pdfID_->first   = -99;  pdfID_->second   = -99;  }
  if (pdfX_)    { pdfX_->first    = -99.; pdfX_->second    = -99.; }
  if (pdfXpdf_) { pdfXpdf_->first = -99.; pdfXpdf_->second = -99.; }

  // CLEAR VECTORS
  if (ttbar_w_) ttbar_w_->clear();
  if (npus_)    npus_->clear();
  if (tnpus_)   tnpus_->clear();
}

void HiEvtTree::GenerateDictionaries(void)
{
  std::vector< std::string > inVList = {
    "vector<int>",
    "vector<float>"
  };
  std::string CWD = getcwd(NULL, 0);
  gSystem->mkdir((CWD+"/cpp").c_str(), kTRUE);
  gSystem->ChangeDirectory((CWD+"/cpp").c_str());
  gInterpreter->AddIncludePath(Form("%s", (CWD+"/cpp").c_str())); // Needed to find the new dictionaries
  for (const auto& d : inVList) { gInterpreter->GenerateDictionary(d.c_str(), "vector"); }
  gSystem->ChangeDirectory(CWD.c_str());
}

#endif

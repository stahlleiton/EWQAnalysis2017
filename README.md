# EWQAnalysis2017
Temporal repository for the W analysis in pPb 8.16 TeV


To use the EWQForest trees in any root macro, you can do the following:

// Add the include file
#include "HiMuonTree.h"

// Add at the beginning
HiMuonTree muonTree = HiMuonTree();
if (!muonTree.GetTree("<PATH TO YOUR FILE LOCAL OR ROOTX>")) return;
Long64_t nentries = muonTree.GetEntries();

// Add in the event loop
if (muonTree.GetEntry(jentry)<0) break;
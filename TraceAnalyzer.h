//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon May  9 17:19:27 2016 by ROOT version 5.34/30
// from TTree tree/tree
// found on file: Good_events_not_normalized.root
//////////////////////////////////////////////////////////

#ifndef TraceAnalyzer_h
#define TraceAnalyzer_h

#include <iostream>

#include <TROOT.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH2.h>
#include <TSelector.h>
#include <TStyle.h>


// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class TraceAnalyzer : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   Float_t         de_l[16];
   Float_t         de_r[16];
   Float_t         cath;
   Float_t         stp0;
   Float_t         stp17;
   Int_t           evt;
   Int_t           seg[16];

   // List of branches
   TBranch        *b_de_l;   //!
   TBranch        *b_de_r;   //!
   TBranch        *b_cath;   //!
   TBranch        *b_stp0;   //!
   TBranch        *b_stp17;   //!
   TBranch        *b_evt;   //!
   TBranch        *b_seg;   //!

   TraceAnalyzer(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~TraceAnalyzer() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(TraceAnalyzer,0);
};

#endif

#ifdef TraceAnalyzer_cxx
void TraceAnalyzer::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("de_l", de_l, &b_de_l);
   fChain->SetBranchAddress("de_r", de_r, &b_de_r);
   fChain->SetBranchAddress("cath", &cath, &b_cath);
   fChain->SetBranchAddress("stp0", &stp0, &b_stp0);
   fChain->SetBranchAddress("stp17", &stp17, &b_stp17);
   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("seg", seg, &b_seg);
}

Bool_t TraceAnalyzer::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef TraceAnalyzer_cxx

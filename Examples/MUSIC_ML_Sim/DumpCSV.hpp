//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jul  2 12:29:48 2019 by ROOT version 6.14/06
// from TTree tree/MUSIC test
// found on file: /home/dasago/Data/MUSIC35_16C/16C13C_401Torr_filtered.root
//////////////////////////////////////////////////////////

#ifndef DumpCSV_h
#define DumpCSV_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TH2.h>
#include <TStyle.h>
#include <TStopwatch.h>
#include <fstream>

// Headers needed by this particular selector


class DumpCSV : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
   //   TTreeReaderValue<Float_t> side = {fReader, "side"};
   //   TTreeReaderValue<Float_t> sier = {fReader, "sier"};
   TTreeReaderArray<Float_t> de_l = {fReader, "de_l"};
   TTreeReaderArray<Float_t> de_r = {fReader, "de_r"};
   //   TTreeReaderArray<Int_t> seg = {fReader, "seg"};
   TTreeReaderValue<Float_t> stp0 = {fReader, "stp0"};
   TTreeReaderValue<Float_t> cath = {fReader, "cath"};
   //  TTreeReaderValue<Float_t> grid = {fReader, "grid"};
   TTreeReaderValue<Float_t> stp17 = {fReader, "stp17"};
   //   TTreeReaderValue<Float_t> tac = {fReader, "tac"};
   TTreeReaderValue<Int_t> reacStp = {fReader, "reacStp"};

   DumpCSV(TTree * /*tree*/ =0) { }
   virtual ~DumpCSV() { }
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
   
   TStopwatch StpWatch;
   int ProgressIndex;
   static const float ProgressFrac[6];
   long double TotalEntries;
   ofstream CSV;   

   ClassDef(DumpCSV,0);

};

#endif

#ifdef DumpCSV_cxx
void DumpCSV::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t DumpCSV::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.
  std::cout << "NOTIFY" << endl;
   return kTRUE;
}


#endif // #ifdef DumpCSV_cxx

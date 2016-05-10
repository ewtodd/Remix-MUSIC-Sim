#define TraceAnalyzer_cxx
// The class definition in TraceAnalyzer.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("TraceAnalyzer.C")
// Root > T->Process("TraceAnalyzer.C","some options")
// Root > T->Process("TraceAnalyzer.C+")
//

#include "TraceAnalyzer.h"

TCanvas* TV;
TFile* FSim;
TTree* SimT;
float sim_de_l[16];
float sim_de_r[16];
float sim_cath;
float sim_stp0;
float sim_stp17;

void TraceAnalyzer::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

   TV = new TCanvas("TV","Trace viewer");

   FSim = new TFile("TDB_4He_p_not_normalized.root");
   SimT = (TTree*)FSim->Get("simt");
   SimT->Print();
   SimT->SetBranchAddress("de_l", sim_de_l);
   SimT->SetBranchAddress("de_r", sim_de_r);
   SimT->SetBranchAddress("cath", &sim_cath);
   SimT->SetBranchAddress("stp0", &sim_stp0);
   SimT->SetBranchAddress("stp17", &sim_stp17);
}

void TraceAnalyzer::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t TraceAnalyzer::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either TraceAnalyzer::GetEntry() or TBranch::GetEntry()
   // to read either all or the required parts of the data. When processing
   // keyed objects with PROOF, the object is already loaded and is available
   // via the fObject pointer.
   //
   // This function should contain the "body" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.


  fChain->GetEntry(entry);
  
  // Fill experimental trace
  TGraph TraceExp;
  for (int stp=0; stp<16; stp++) {
    TraceExp.SetPoint(stp, stp+1, de_l[stp] + de_r[stp]);
  }

  // Fill simulated trace and compare with experimental trace
  TGraph TraceSim;  
  for (int sevt=0; sevt<SimT->GetEntries(); sevt++) {
    SimT->GetEntry(sevt);
    for (int stp=0; stp<16; stp++) {
      TraceSim.SetPoint(stp, stp+1, sim_de_l[stp] + sim_de_r[stp]); 
    }
    TV->cd();
    TraceExp.Draw("al*");
    TraceSim.Draw("l same");
    TV->Update();
    TV->WaitPrimitive();
  }



   return kTRUE;
}

void TraceAnalyzer::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void TraceAnalyzer::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

}

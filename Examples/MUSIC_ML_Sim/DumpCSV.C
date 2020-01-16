#define DumpCSV_cxx
// The class definition in DumpCSV.h has been generated automatically
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
// root> T->Process("DumpCSV.C")
// root> T->Process("DumpCSV.C","some options")
// root> T->Process("DumpCSV.C+")
//


#include "DumpCSV.h"

const float DumpCSV::ProgressFrac[6] = {0.01, 0.25, 0.5, 0.75, 0.9, 1.0};
Long64_t mainentry = 0;

void DumpCSV::Begin(TTree* tree)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();
  cout << option << endl;
  if (option=="1M")        mainentry = 1000000;
  else if (option=="2M")   mainentry = 2000000;
  else if (option=="3M")   mainentry = 3000000;
  else if (option=="4M")   mainentry = 4000000;
  else if (option=="5M")   mainentry = 5000000;
  else if (option=="6M")   mainentry = 6000000;
  else if (option=="7M")   mainentry = 7000000;
  else if (option=="8M")   mainentry = 8000000;
  else if (option=="9M")   mainentry = 9000000;
  else if (option=="10M")  mainentry = 10000000;
  else if (option=="11M")  mainentry = 11000000;
  else if (option=="12M")  mainentry = 12000000;
  else if (option=="13M")  mainentry = 13000000;
  else if (option=="14M")  mainentry = 14000000;
  else if (option=="15M")  mainentry = 15000000;
  else if (option=="16M")  mainentry = 16000000;
  else if (option=="17M")  mainentry = 17000000;
  else if (option=="18M")  mainentry = 18000000;
  else if (option=="19M")  mainentry = 19000000;
  else if (option=="20M")  mainentry = 20000000;

  TotalEntries = tree->GetEntries();
  cout << "Processing " << TotalEntries << " entries" << endl;
  StpWatch.Start();
  ProgressIndex = 0;
  CSV.open("17F_alpha_p_sim.csv", std::ofstream::out | std::ofstream::app);
}

void DumpCSV::SlaveBegin(TTree * /*tree*/)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();

}

Bool_t DumpCSV::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // When processing keyed objects with PROOF, the object is already loaded
  // and is available via the fObject pointer.
  //
  // This function should contain the \"body\" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.
  //
  // The processing can be stopped by calling Abort().
  //
  // Use fStatus to set the return value of TTree::Process().
  //
  // The return value is currently not used.
  
  // Simple progress monitor
  if (TotalEntries>1e4) {
    if ((long double)entry>=ProgressFrac[ProgressIndex]*TotalEntries) {
      cout << "\t" << ProgressFrac[ProgressIndex]*100 << "% processed (" 
	   << StpWatch.RealTime() << " s)" << endl;
      StpWatch.Start(kFALSE);
      if (ProgressIndex<6)
	ProgressIndex++;
    }
  }

  
  fReader.SetEntry(entry);
  //if (entry<10) 
    {
      CSV << mainentry << "," << *stp0 << ",";
      for (int i=0; i<16; i++) 
	CSV << de_l[i] << "," << de_r[i] << ",";
      CSV << *stp17 << "," << *cath << ",0," << *reacStp << endl;
    }
  
  mainentry++;

  return kTRUE;
}

void DumpCSV::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

}

void DumpCSV::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
  CSV.close();

}

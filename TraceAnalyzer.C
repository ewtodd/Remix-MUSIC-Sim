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

using namespace std;

// Reactions
// 0 - (alpha,alpha)
// 1 - (alpha,p)
// 2 - (alpha,n)
const int NumReac = 3;
string FileName[NumReac] = {"TDB_4He_4He_18000.root",
			    "TDB_4He_p_18000.root",
			    "TDB_4He_n_18000.root"};
string ReacName[NumReac] = {"(#alpha,#alpha)",
			    "(#alpha,p)",
			    "(#alpha,n)"};
int ReacEvts[NumReac] = {0,0,0};

TCanvas* TV;
TFile** FSim;
TTree** SimT;
TH2F* Bkgnd;
float sim_de_l[20];
float sim_de_r[20];
float sim_cath;
float sim_stp0;
float sim_stp17;


void TraceAnalyzer::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();

  Bkgnd = new TH2F("Bkgnd","",18,-0.5,17.5,1000,0,5);
  Bkgnd->GetXaxis()->SetTitle("Strip number");
  Bkgnd->GetXaxis()->CenterTitle();
  Bkgnd->GetYaxis()->SetTitle("#DeltaE [MeV]");
  Bkgnd->GetYaxis()->CenterTitle();
  Bkgnd->SetStats(0);

  TV = new TCanvas("TV","Trace viewer");
  TV->Divide(2,2);
  TV->cd(1)->SetGrid();
  TV->cd(2)->SetGrid();
  TV->cd(3)->SetGrid();
  Bkgnd->Draw();

  FSim = new TFile*[NumReac];
  SimT = new TTree*[NumReac];
  for (int r=0; r<NumReac; r++) {
    FSim[r] = new TFile(FileName[r].c_str());
    SimT[r] = (TTree*)FSim[r]->Get("simt");
    // NOTE: Different trees using the same address for their leaves
    // means we cannot get their entries in parallel. Access to the
    // leaves' information must be sequential.
    SimT[r]->SetBranchAddress("de_l", sim_de_l);
    SimT[r]->SetBranchAddress("de_r", sim_de_r);
    SimT[r]->SetBranchAddress("cath", &sim_cath);
    SimT[r]->SetBranchAddress("stp0", &sim_stp0);
    SimT[r]->SetBranchAddress("stp17", &sim_stp17);
    SimT[r]->Print();
  }
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
  cout << "Entry = " << entry << endl;
  // Fill experimental trace
  TGraph TraceExp;
  TGraph TraceLExp;
  TGraph TraceRExp;
  // Trace style
  TraceExp.SetLineColor(kBlack);
  TraceExp.SetLineWidth(2);
  TraceLExp.SetLineColor(kRed);
  TraceLExp.SetLineWidth(2);
  TraceRExp.SetLineColor(kBlue);
  TraceRExp.SetLineWidth(2);
  // Fill traces with experimental data
  TraceExp.SetPoint(0, 0, stp0);
  for (int stp=0; stp<16; stp++) {
    TraceExp.SetPoint(stp+1, stp+1, de_l[stp] + de_r[stp]);
    TraceLExp.SetPoint(stp+1, stp+1, de_l[stp]);
    TraceRExp.SetPoint(stp+1, stp+1, de_r[stp]);
  }
  TraceExp.SetPoint(17, 17, stp17);
  // Fill simulated trace and compare with experimental trace
  TGraph TraceSim;
  TGraph TraceRSim;
  TGraph TraceLSim;
  TraceSim.SetLineColor(kGray+2);
  TraceSim.SetLineWidth(2);
  TraceSim.SetLineStyle(2);
  TraceLSim.SetLineColor(kCyan);
  TraceLSim.SetLineWidth(2);
  TraceLSim.SetLineStyle(2);
  TraceRSim.SetLineColor(kMagenta);
  TraceRSim.SetLineWidth(2);
  TraceRSim.SetLineStyle(2);
  double ChiT = 0;
  double ChiL = 0;
  double ChiR = 0;
  double ChiLR = 0;
  double ChiTMin = 1000;
  double ChiLMin = 1000;
  double ChiRMin = 1000;
  double ChiLRMin = 1000;
  int NT = 0;
  int NL = 0;
  int NR = 0;
  int NLR = 0;
  int ReacBestFit = -1;
  for (int r=0; r<NumReac; r++) {
    cout << FileName[r] << endl;
    for (int sevt=0; sevt<SimT[r]->GetEntries(); sevt++) {
      SimT[r]->GetEntry(sevt);
      NT = 0;
      NL = 0;
      NR = 0;
      NLR = 0;
      TraceSim.SetPoint(0, 0, sim_stp0);
      for (int stp=0; stp<16; stp++) {
	TraceSim.SetPoint(stp+1, stp+1, sim_de_l[stp] + sim_de_r[stp]); 
	TraceRSim.SetPoint(stp+1, stp+1, sim_de_r[stp]); 
	TraceLSim.SetPoint(stp+1, stp+1, sim_de_l[stp]);
	// TraceLRSim.SetPoint(stp+1, stp+1, sim_de_l[stp]);
	// TraceLRSim.SetPoint(stp+1, stp+16, sim_de_l[stp]);
	// Chi square sums
	if (sim_de_l[stp]+sim_de_r[stp]>0) {
	  ChiT += pow(de_l[stp]+de_r[stp]-sim_de_l[stp]-sim_de_r[stp],2)/(sim_de_l[stp]+sim_de_r[stp]);
	  NT++;
	}
	if (sim_de_l[stp]>0) {
	  ChiL += pow(de_l[stp] - sim_de_l[stp],2)/sim_de_l[stp];
	  NL++;
	}
	if (sim_de_r[stp]>0) {
	  ChiR += pow(de_r[stp] - sim_de_r[stp],2)/sim_de_r[stp];
	  NR++;
	}
      }
      TraceSim.SetPoint(17, 17, sim_stp17);
      if (NT>1) 
	ChiT = ChiT/(NT-1);
      else 
	ChiT = 1000;
      if (NL>1)
	ChiL = ChiL/(NL-1);
      else 
	ChiL = 1000;
      if (NR>1)
	ChiR = ChiR/(NR-1);
      else
	ChiR = 1000;

      if (ChiT<ChiTMin && ChiL<ChiLMin && ChiR<ChiRMin) {
	cout << "ChiT = " << ChiT << "  " << NT;
	cout << "\t ChiL = " << ChiL << "  " << NL;
	cout << "\t ChiR = " << ChiR << "  " << NR << endl;
	if (ChiT<ChiTMin)
	  ChiTMin = ChiT;
	if (ChiL<ChiLMin)
	  ChiLMin = ChiL;
	if (ChiR<ChiRMin)
	  ChiRMin = ChiR;
	Bkgnd->SetTitle(Form("Exp entry = %d | %s | Sim entry = %d",entry,ReacName[r].c_str(),sevt));
	// TV->cd(1);
	// Bkgnd->Draw();
	// TraceExp.Draw("l* same");
	// TraceSim.Draw("l same");
	// TV->cd(2);
	// Bkgnd->Draw();
	// TraceLExp.Draw("l* same");
	// TraceLSim.Draw("l same");
	// TV->cd(3);
	// Bkgnd->Draw();
	// TraceRExp.Draw("l* same");
	// TraceRSim.Draw("l same");

	ReacBestFit = r;

	//	TV->Update();
	//	TV->WaitPrimitive();
      }
    }
  } // end for(r<NumReac)
  if (ReacBestFit>=0 && ReacBestFit<NumReac) {
    cout << " Best chi! " << ReacName[ReacBestFit] << "\n" << endl;
    ReacEvts[ReacBestFit] += 1;
  }
  else
    cout << " Bad fit!! " << ReacBestFit << endl;
  
  //  TV->Update();
  // TV->WaitPrimitive();
    
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

  cout << "Number of events for:" << endl;
  for (int r=0; r<NumReac; r++)
    cout << ReacName[r] << " = " << ReacEvts[r] << endl;
}

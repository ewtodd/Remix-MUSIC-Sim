{
  TChain* tree = new TChain("tree","tree");
  tree->AddFile("Good_events_an_ap_aa.root");
  tree->Process("TraceAnalyzer.C+");
}


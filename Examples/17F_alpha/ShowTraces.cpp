{

  int NTraces = 500;
  int stp = 4;
  int Kb = 43; // MeV
  int pressure = 470;
  
  TFile* TFap = new TFile("traces_gs.root");

  TGraph** Tap = new TGraph*[NTraces];
  for (int n=0; n<NTraces; n++) {
    Tap[n] = (TGraph*)TFap->Get(Form("traces/Trace%d",n));
    if (Tap[n]!=0) {
      Tap[n]->SetLineColor(kRed);
      Tap[n]->SetLineWidth(2);
      Tap[n]->SetMarkerColor(kRed);
    }
  }
  

  TCanvas* Can = new TCanvas("Can","Traces",0,0,1000,800);
  TH2F* HELoss = new TH2F("HELoss","^{88}Sr+^{4}He (600 Torr, E_{b} = 350 MeV)", 
			  18,-0.5,18-0.5, 400,0.0,30);
  HELoss->SetStats(0);
  Can->SetGrid();
  HELoss->GetXaxis()->SetTitle("Segment number");
  HELoss->GetXaxis()->CenterTitle();
  HELoss->GetYaxis()->SetTitle("#DeltaE [MeV]");
  HELoss->GetYaxis()->CenterTitle(); 
  HELoss->GetXaxis()->SetLimits(0,18);
  //  HELoss->GetYaxis()->SetLimits(0.0,10);
  HELoss->Draw();
  for (int n=0; n<NTraces; n++) 
    if (Tap[n]!=0)
      Tap[n]->Draw("l same");
  TLegend* Leg = new TLegend(0.174,0.156,0.346,0.374);
  //  Leg->AddEntry(Taa[0], "(#alpha,#alpha)", "l");
  Leg->AddEntry(Tap[0], "(#alpha,p)", "l");
  //Leg->AddEntry(Tan[0], "(#alpha,n)", "l");
  Leg->Draw();
}

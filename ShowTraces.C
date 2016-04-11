{
  int NTraces = 50;
  int stp = 4;
 
  TFile* TFaa = new TFile(Form("TracesR0_Stp%d_4He_4He.root",stp));
  TFile* TFap = new TFile(Form("TracesR0_Stp%d_4He_p.root",stp));
  TFile* TFan = new TFile(Form("TracesR0_Stp%d_4He_n.root",stp));

  TGraph** Taa = new TGraph*[NTraces];
  for (int n=0; n<NTraces; n++) {
    Taa[n] = (TGraph*)TFaa->Get(Form("Trace%d",n));
    if (Taa[n]!=0) {
      Taa[n]->SetLineColor(kGreen);
      Taa[n]->SetLineWidth(2);
      Taa[n]->SetMarkerColor(kGreen);
    }
  }
  TGraph** Tap = new TGraph*[NTraces];
  for (int n=0; n<NTraces; n++) {
    Tap[n] = (TGraph*)TFap->Get(Form("Trace%d",n));
    if (Tap[n]!=0) {
      Tap[n]->SetLineColor(kRed);
      Tap[n]->SetLineWidth(2);
      Tap[n]->SetMarkerColor(kRed);
    }
  }
  TGraph** Tan = new TGraph*[NTraces];
  for (int n=0; n<NTraces; n++) {
    Tan[n] = (TGraph*)TFan->Get(Form("Trace%d",n));
    if (Tan[n]!=0) {
      Tan[n]->SetLineColor(kBlue);
      Tan[n]->SetLineWidth(2);
      Tan[n]->SetMarkerColor(kBlue);
    }
  }


  TCanvas* Can = new TCanvas("Can","Traces",0,0,1000,800);
  TH2F* HELoss = new TH2F("HELoss","^{20}Ne+^{4}He (350 Torr, E_{b} = 55 MeV)", 18,-0.5, 18-0.5, 400,0.2,2);
  HELoss->SetStats(0);
  Can->SetGrid();
  HELoss->GetXaxis()->SetTitle("Segment number");
  HELoss->GetXaxis()->CenterTitle();
  HELoss->GetYaxis()->SetTitle("#DeltaE [MeV]");
  HELoss->GetYaxis()->CenterTitle(); 
  HELoss->GetXaxis()->SetLimits(0,18);
  HELoss->GetYaxis()->SetLimits(0.0,10);
  HELoss->Draw();
  for (int n=0; n<NTraces; n++) 
    Taa[n]->Draw("l same");
  for (int n=0; n<NTraces; n++) 
    Tap[n]->Draw("l same");
  for (int n=0; n<NTraces; n++) 
    Tan[n]->Draw("l same");

}

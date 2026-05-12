{

  int NTraces = 30;
  int stp = 3;
  int Kb = 350; // MeV
  int pressure = 600;
  
  // TFile* TFaa = new TFile(Form("traces_aa_st%d_500Torr_lise_ATIMA_HI.root",stp));
  //  TFile* TFap = new TFile(Form("Traces_Stp%d_ap_%dMeV_%dTorr.root",stp));
  // TFile* TFan = new TFile(Form("traces_an_st%d_500Torr_lise_ATIMA_HI.root",stp));
  //  TFile* TFab = new TFile(Form("Traces_Stp%d_aa_%dMeV_%dTorr.root",18));

  TFile* TFaa = new TFile(Form("traces_aa_st%d_500Torr_lise_ATIMA.root",stp));
  TFile* TFan = new TFile(Form("traces_an_st%d_500Torr_lise_ATIMA.root",stp));

  // TGraph** Tab = new TGraph*[NTraces];
  // for (int n=0; n<NTraces; n++) {
  //   Tab[n] = (TGraph*)TFab->Get(Form("traces/Trace%d",n));
  //   if (Tab[n]!=0) {
  //     Tab[n]->SetLineColor(kBlack);
  //     Tab[n]->SetLineWidth(2);
  //     Tab[n]->SetMarkerColor(kBlack);
  //   }
  // }

  TGraph** Taa = new TGraph*[NTraces];
  for (int n=0; n<NTraces; n++) {
    Taa[n] = (TGraph*)TFaa->Get(Form("traces/Trace_s%d_e%d",stp,n));
    if (Taa[n]!=0) {
      Taa[n]->SetLineColor(kGreen);
      Taa[n]->SetLineWidth(2);
      Taa[n]->SetMarkerColor(kGreen);
    }
  }
  
  // TGraph** Tap = new TGraph*[NTraces];
  // for (int n=0; n<NTraces; n++) {
  //   Tap[n] = (TGraph*)TFap->Get(Form("traces/Trace%d",n));
  //   if (Tap[n]!=0) {
  //     Tap[n]->SetLineColor(kRed);
  //     Tap[n]->SetLineWidth(2);
  //     Tap[n]->SetMarkerColor(kRed);
  //   }
  // }
  
  TGraph** Tan = new TGraph*[NTraces];
  for (int n=0; n<NTraces; n++) {
    Tan[n] = (TGraph*)TFan->Get(Form("traces/Trace_s%d_e%d",stp,n));
    if (Tan[n]!=0) {
      Tan[n]->SetLineColor(kBlue);
      Tan[n]->SetLineWidth(2);
      Tan[n]->SetMarkerColor(kBlue);
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
  // for (int n=0; n<NTraces; n++) 
  //   if (Tab[n]!=0)
  //     Tab[n]->Draw("l same");
  for (int n=0; n<NTraces; n++)
    if (Taa[n]!=0)
      Taa[n]->Draw("l same");
  // for (int n=0; n<NTraces; n++) 
  //   if (Tap[n]!=0)
  //     Tap[n]->Draw("l same");
  for (int n=0; n<NTraces; n++) 
    if (Tan[n]!=0)
      Tan[n]->Draw("l same");
  TLegend* Leg = new TLegend(0.174,0.156,0.346,0.374);
  Leg->AddEntry(Taa[0], "(#alpha,#alpha)", "l");
  //Leg->AddEntry(Tap[0], "(#alpha,p)", "l");
  Leg->AddEntry(Tan[0], "(#alpha,n)", "l");
  Leg->Draw();
}

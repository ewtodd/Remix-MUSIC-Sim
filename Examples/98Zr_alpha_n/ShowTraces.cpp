{

  int NTraces = 100;
  int stp = 4;
  int Kb = 400; // MeV
  int pressure = 500;
  
  // TFile* TFaa = new TFile(Form("Traces_Stp%d_aa_%dMeV_%dTorr.root",stp,Kb,pressure));
  // TFile* TFap = new TFile(Form("Traces_Stp%d_ap_%dMeV_%dTorr.root",stp,Kb,pressure));
  // TFile* TFan = new TFile(Form("Traces_Stp%d_an_%dMeV_%dTorr.root",stp,Kb,pressure));
  // TFile* TFab = new TFile(Form("Traces_Stp%d_aa_%dMeV_%dTorr.root",18,Kb,pressure));
  TFile* TFaa = new TFile("traces_aa_stp12.root");
  TFile* TFa2n = new TFile("traces_a2n_stp12.root");
  TFile* TFan = new TFile("traces_an_stp12.root");
  //  TFile* TFab = new TFile(Form("Traces_Stp%d_aa_%dMeV_%dTorr.root",18,Kb,pressure));

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
    Taa[n] = (TGraph*)TFaa->Get(Form("traces/Trace%d",n));
    if (Taa[n]!=0) {
      Taa[n]->SetLineColor(kGreen);
      Taa[n]->SetLineWidth(2);
      Taa[n]->SetMarkerColor(kGreen);
    }
  }
  
  TGraph** Ta2n = new TGraph*[NTraces];
  for (int n=0; n<NTraces; n++) {
    Ta2n[n] = (TGraph*)TFa2n->Get(Form("traces/Trace%d",n));
    if (Ta2n[n]!=0) {
      Ta2n[n]->SetLineColor(kRed);
      Ta2n[n]->SetLineWidth(2);
      Ta2n[n]->SetMarkerColor(kRed);
    }
  }
  
  TGraph** Tan = new TGraph*[NTraces];
  for (int n=0; n<NTraces; n++) {
    Tan[n] = (TGraph*)TFan->Get(Form("traces/Trace%d",n));
    if (Tan[n]!=0) {
      Tan[n]->SetLineColor(kBlue);
      Tan[n]->SetLineWidth(2);
      Tan[n]->SetMarkerColor(kBlue);
    }
  }


  TCanvas* Can = new TCanvas("Can","Traces",0,0,1000,800);
  TH2F* HELoss = new TH2F("HELoss","^{98}Zr+^{4}He (500 Torr, E_{b} = 400 MeV)", 
			  18,-0.5,18-0.5, 400,0.0,30);
  HELoss->SetStats(0);
  Can->SetGrid();
  HELoss->GetXaxis()->SetTitle("strip");
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
  for (int n=0; n<NTraces; n++) 
    if (Ta2n[n]!=0)
      Ta2n[n]->Draw("l same");
  for (int n=0; n<NTraces; n++) 
    if (Tan[n]!=0)
      Tan[n]->Draw("l same");
  TLegend* Leg = new TLegend(0.174,0.156,0.346,0.374);
  Leg->AddEntry(Taa[0], "(#alpha,#alpha)", "l");
  Leg->AddEntry(Tan[0], "(#alpha,n)", "l");
  Leg->AddEntry(Ta2n[0], "(#alpha,2n)", "l");
  Leg->Draw();
}

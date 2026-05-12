{

  int NTraces = 100;
  int stp = 6;
  int Kb = 79; // MeV
  int pressure = 220;
  
  TFile* TFap = new TFile("traces_ap_stp6_220Torr_lise_79MeV.root");
  TFile* TFaa = new TFile("traces_aa_stp6_220Torr_lise_79MeV.root");
  TFile* TFan = new TFile("traces_an_stp6_220Torr_lise_79MeV.root");
  TFile* TFub = new TFile("unreacted_37Cl_220Torr_lise_79MeV.root");

  TGraph** Tub = new TGraph*[NTraces];
  for (int n=0; n<NTraces; n++) {
    Tub[n] = (TGraph*)TFub->Get(Form("traces/Trace%d",n));
    if (Tub[n]!=0) {
      Tub[n]->SetLineColor(kBlack);
      Tub[n]->SetLineWidth(2);
      Tub[n]->SetMarkerColor(kBlack);
    }
  }

  TGraph** Taa = new TGraph*[NTraces];
  for (int n=0; n<NTraces; n++) {
    Taa[n] = (TGraph*)TFaa->Get(Form("traces/Trace%d",n));
    if (Taa[n]!=0) {
      Taa[n]->SetLineColor(kGray+1);
      Taa[n]->SetLineWidth(2);
      Taa[n]->SetMarkerColor(kGray+1);
    }
  }

  TGraph** Tap = new TGraph*[NTraces];
  for (int n=0; n<NTraces; n++) {
    Tap[n] = (TGraph*)TFap->Get(Form("traces/Trace%d",n));
    if (Tap[n]!=0) {
      Tap[n]->SetLineColor(kRed);
      Tap[n]->SetLineWidth(2);
      Tap[n]->SetMarkerColor(kRed);
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
  TH2F* HELoss = new TH2F("HELoss",Form("^{37}Cl+^{4}He (He at %d Torr, E_{beam} = %d MeV)",pressure, Kb),
			  18,-0.5,18-0.5, 400,1.0,5.0);
  HELoss->SetStats(0);
  Can->SetGrid();
  HELoss->GetXaxis()->SetTitle("Strip number (as in AnodeGeometry file)");
  HELoss->GetXaxis()->CenterTitle();
  HELoss->GetYaxis()->SetTitle("#DeltaE [MeV]");
  HELoss->GetYaxis()->CenterTitle(); 
  HELoss->GetXaxis()->SetLimits(1,18);
  //  HELoss->GetYaxis()->SetLimits(0.0,10);
  HELoss->Draw();
  for (int n=0; n<NTraces; n++) 
    if (Tub[n]!=0)
      Tub[n]->Draw("l same");
  for (int n=0; n<NTraces; n++)
    if (Taa[n]!=0)
      Taa[n]->Draw("l same");
  for (int n=0; n<NTraces; n++)
    if (Tap[n]!=0)
      Tap[n]->Draw("l same");
  for (int n=0; n<NTraces; n++) 
    if (Tan[n]!=0)
      Tan[n]->Draw("l same");
  TLegend* Leg = new TLegend(0.174,0.156,0.346,0.374);
  Leg->AddEntry(Tan[0], "(#alpha,n)", "l");
  Leg->AddEntry(Tap[0], "(#alpha,p)", "l");
  Leg->AddEntry(Taa[0], "(#alpha,#alpha)", "l");
  Leg->AddEntry(Tub[0], "beam", "l");
  Leg->Draw();
}

{

  int NTraces = 100;
  int stp = 4;
  int Kb = 700; // MeV
  int pressure = 480;
  
  TFile* TFpa = new TFile("traces_pa_stp4.root");
  TFile* TFpp = new TFile("traces_pp_stp4.root");
  TFile* TFub = new TFile("unreacted_70Zn.root");

  TGraph** Tub = new TGraph*[NTraces];
  for (int n=0; n<NTraces; n++) {
    Tub[n] = (TGraph*)TFub->Get(Form("traces/Trace%d",n));
    if (Tub[n]!=0) {
      Tub[n]->SetLineColor(kBlack);
      Tub[n]->SetLineWidth(2);
      Tub[n]->SetMarkerColor(kBlack);
    }
  }

  TGraph** Tpa = new TGraph*[NTraces];
  for (int n=0; n<NTraces; n++) {
    Tpa[n] = (TGraph*)TFpa->Get(Form("traces/Trace%d",n));
    if (Tpa[n]!=0) {
      Tpa[n]->SetLineColor(kRed);
      Tpa[n]->SetLineWidth(2);
      Tpa[n]->SetMarkerColor(kRed);
    }
  }
  
  TGraph** Tpp = new TGraph*[NTraces];
  for (int n=0; n<NTraces; n++) {
    Tpp[n] = (TGraph*)TFpp->Get(Form("traces/Trace%d",n));
    if (Tpp[n]!=0) {
      Tpp[n]->SetLineColor(kBlue);
      Tpp[n]->SetLineWidth(2);
      Tpp[n]->SetMarkerColor(kBlue);
    }
  }
  

  TCanvas* Can = new TCanvas("Can","Traces",0,0,1000,800);
  TH2F* HELoss = new TH2F("HELoss",Form("^{70}Zn+^{1}H (CH4 at %d Torr, E_{b} = %d MeV)",pressure, Kb),
			  18,-0.5,18-0.5, 400,0.0,50);
  HELoss->SetStats(0);
  Can->SetGrid();
  HELoss->GetXaxis()->SetTitle("Strip number (as in AnodeGeometry file)");
  HELoss->GetXaxis()->CenterTitle();
  HELoss->GetYaxis()->SetTitle("#DeltaE [MeV]");
  HELoss->GetYaxis()->CenterTitle(); 
  HELoss->GetXaxis()->SetLimits(0,18);
  //  HELoss->GetYaxis()->SetLimits(0.0,10);
  HELoss->Draw();
  for (int n=0; n<NTraces; n++) 
    if (Tub[n]!=0)
      Tub[n]->Draw("l same");
  for (int n=0; n<NTraces; n++)
    if (Tpa[n]!=0)
      Tpa[n]->Draw("l same");
  for (int n=0; n<NTraces; n++) 
    if (Tpp[n]!=0)
      Tpp[n]->Draw("l same");
  TLegend* Leg = new TLegend(0.174,0.156,0.346,0.374);
  Leg->AddEntry(Tpa[0], "(p,#alpha)", "l");
  Leg->AddEntry(Tpp[0], "(p,p)", "l");
  Leg->AddEntry(Tub[0], "beam", "l");
  Leg->Draw();
}

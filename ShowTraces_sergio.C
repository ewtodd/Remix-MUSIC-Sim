{

  int NTraces = 30;
  int stp = 5;
 

  TFile* TF1 = new TFile(Form("Traces_Stp%d_p_4He.root",stp));
  TFile* TF2 = new TFile(Form("Traces_Stp%d_p_4He.root",18));

  TGraph** Tr1 = new TGraph*[NTraces];
  for (int n=0; n<NTraces; n++) {
    Tr1[n] = (TGraph*)TF1->Get(Form("Trace%d",n));
    if (Tr1[n]!=0) {
      Tr1[n]->SetLineColor(kRed);
      Tr1[n]->SetLineWidth(2);
      Tr1[n]->SetMarkerColor(kRed);
    }
  }

  TGraph** Tr2 = new TGraph*[NTraces];
  for (int n=0; n<NTraces; n++) {
    Tr2[n] = (TGraph*)TF2->Get(Form("Trace%d",n));
    if (Tr2[n]!=0) {
      Tr2[n]->SetLineColor(kBlack);
      Tr2[n]->SetLineWidth(2);
      Tr2[n]->SetMarkerColor(kBlack);
    }
  }

  TCanvas* Can = new TCanvas("Can","Traces",0,0,1000,800);
  TH2F* HELoss = new TH2F("HELoss","^{31}P+p (CH4 at 300 Torr, E_{b} = 8 MeV/u)", 
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
    if (Tr2[n]!=0)
      Tr2[n]->Draw("l same");
  for (int n=0; n<NTraces; n++) 
    if (Tr1[n]!=0)
      Tr1[n]->Draw("l same");

  TLegend* Leg = new TLegend(0.174,0.156,0.346,0.374);
  Leg->AddEntry(Tr1[0], "(p,#alpha)", "l");
  Leg->AddEntry(Tr2[0], "beam", "l");
  Leg->Draw();
}

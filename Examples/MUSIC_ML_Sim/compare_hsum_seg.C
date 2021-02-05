{
  TCanvas* can = new TCanvas("can");
  can->Divide(3);
  
  // Experimental data
  TFile* f1 = new TFile("hsumseg_run21.root");
  TH2F* h1 = (TH2F*)f1->Get("hsum_seg");
  can->cd(1)->SetGrid();
  h1->Draw("colz");

  
  TH2F* hsim = (TH2F*)h1->Clone("hsim");
  hsim->Reset();
  hsim->SetTitle("uncalibrated simulation data");
  //  hsim->Draw();
 

  TFile* f2 = new TFile("fE37_5P406.root");
  TH2F* h2 = (TH2F*)f2->Get("hE37_5P406");
  can->cd(2)->SetGrid();
  h2->Draw("colz");

  // Fill hsim with the simulated data of h2 but "uncalibrate" it
  cout << h2->GetNbinsX() << endl;
  cout << h2->GetNbinsY() << endl;
  for (int bx=0; bx<h2->GetNbinsX(); bx++) { 
    int stp = h2->GetXaxis()->GetBinCenter(bx);
    for (int by=0; by<h2->GetNbinsY(); by++) {
      float de = h2->GetYaxis()->GetBinCenter(by);
      float counts = h2->GetBinContent(h2->GetBin(bx,by));
      float uncalde = 800*de - 100.0;
      hsim->Fill(stp, uncalde, counts);
    }
  }
  can->cd(3)->SetGrid();
  hsim->Draw("colz");
  
}

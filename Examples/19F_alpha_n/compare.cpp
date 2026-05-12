{
  gStyle->SetOptStat(0);
  TCanvas* can = new TCanvas("can","A=98 beams");
  can->Divide(2);
  
  TFile* fSr = new TFile("stp0_98Sr.root");
  TH1F* hSr = (TH1F*)fSr->Get("hSr");
  hSr->SetLineColor(3);

  TFile* fY = new TFile("stp0_98Y.root");
  TH1F* hY = (TH1F*)fY->Get("hY");
  hY->SetLineColor(2);

  TFile* fZr = new TFile("stp0_98Zr.root");
  TH1F* hZr = (TH1F*)fZr->Get("hZr");
  hZr->SetLineColor(1);
  
  TFile* fNb = new TFile("stp0_98Nb.root");
  TH1F* hNb = (TH1F*)fNb->Get("hNb");
  hNb->SetLineColor(4);

  can->cd(1);
  hSr->Draw();
  hY->Draw("same");
  hZr->Draw("same");
  hNb->Draw("same");

  auto leg = new TLegend(0.1,0.7,0.48,0.9);
  leg->AddEntry(hNb,"98Nb","l");
  leg->AddEntry(hZr,"98Zr","l");
  leg->AddEntry(hY,"98Y","l");
  leg->AddEntry(hSr,"98Sr","l");
  leg->Draw();
  
#if 0
  TFile* f2Sr = new TFile("stp0stp1_98Sr.root");
  TH1F* h2Sr = (TH1F*)f2Sr->Get("h2Sr");
  h2Sr->SetMarkerColor(3);

  TFile* f2Y = new TFile("stp0stp1_98Y.root");
  TH1F* h2Y = (TH1F*)f2Y->Get("h2Y");
  h2Y->SetMarkerColor(2);

  TFile* f2Zr = new TFile("stp0stp1_98Zr.root");
  TH1F* h2Zr = (TH1F*)f2Zr->Get("h2Zr");
  h2Zr->SetMarkerColor(1);
  
  TFile* f2Nb = new TFile("stp0stp1_98Nb.root");
  TH1F* h2Nb = (TH1F*)f2Nb->Get("h2Nb");
  h2Nb->SetMarkerColor(4);
#endif
  
  TFile* f3Sr = new TFile("stp0stp123_98Sr.root");
  TH1F* h3Sr = (TH1F*)f3Sr->Get("h3Sr");
  h3Sr->SetMarkerColor(3);
  
  TFile* f3Y = new TFile("stp0stp123_98Y.root");
  TH1F* h3Y = (TH1F*)f3Y->Get("h3Y");
  h3Y->SetMarkerColor(2);

  TFile* f3Zr = new TFile("stp0stp123_98Zr.root");
  TH1F* h3Zr = (TH1F*)f3Zr->Get("h3Zr");
  h3Zr->SetMarkerColor(1);
  
  TFile* f3Nb = new TFile("stp0stp123_98Nb.root");
  TH1F* h3Nb = (TH1F*)f3Nb->Get("h3Nb");
  h3Nb->SetMarkerColor(4);
  
  can->cd(2);
  h3Sr->Draw();
  h3Y->Draw("same");
  h3Zr->Draw("same");
  h3Nb->Draw("same");

  
}
